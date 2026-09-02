"""
MolSim Final — Backend Flask (v3.3)
• Input validation
• Improved error handling
• Type hints
"""

import json
import logging
import math
import os
import traceback
import base64
import binascii
import io
import struct
import tempfile
import time
from datetime import datetime
from zoneinfo import ZoneInfo
from typing import Tuple, Dict, Any, Optional, List
from urllib import request as urllib_request, parse as urllib_parse

from flask import Flask, render_template, request, jsonify, send_file, make_response
from flask_cors import CORS
from fpdf import FPDF
from werkzeug.exceptions import HTTPException, RequestEntityTooLarge

from analysis import compare, bulk_compare, build_chemical_space, get_mol, mol_to_png, validate_fingerprint_type, validate_metric
from chemo_suite.apps.mol_sim.pairwise import run_pairwise_compare
from chemo_suite.apps.mol_sim.batch import run_batch_compare
from chemo_suite.apps.nitro_ra.cpca import calculate_cpca
from chemo_suite.apps.nitro_ra.deep_pk import get_deep_pk_metabolism, submit_deep_pk_metabolism
from chemo_suite.apps.nitro_ra.metabolism import evaluate_metabolism
from chemo_suite.apps.nitro_ra.quantum import evaluate_quantum
from chemo_suite.apps.nitro_ra.nitrosamine_space import search_nitrosamine_space
from molsim_runtime import get_env_int, parse_cors_origins, make_error_payload

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)
REPORT_TIMEZONE = ZoneInfo("America/Sao_Paulo")


def get_report_generated_at() -> str:
    """Retorna o momento de geração no fuso horário de Brasília."""
    return datetime.now(REPORT_TIMEZONE).strftime("%d/%m/%Y %H:%M")


def ensure_report_generated_at(data: Dict[str, Any]) -> Dict[str, Any]:
    """Garante timestamp brasileiro sem sobrescrever um valor já informado."""
    normalized = dict(data)
    generated_at = normalized.get("generated_at")
    if not isinstance(generated_at, str) or not generated_at.strip():
        normalized["generated_at"] = get_report_generated_at()
    return normalized


class ExportReportPDF(FPDF):
    """PDF direto com rodapé uniforme em todas as páginas."""

    generated_at = "-"

    def footer(self):
        self.set_y(-10)
        self.set_font("Helvetica", "I", 7)
        self.set_text_color(100, 116, 139)
        self.cell(0, 5, f"PharmaSci · Gerado em {self.generated_at} · Página {self.page_no()}", align="C")


app = Flask(__name__)


MAX_CONTENT_LENGTH = get_env_int("MOLSIM_MAX_CONTENT_LENGTH", 2 * 1024 * 1024)
MAX_BATCH_ITEMS = get_env_int("MOLSIM_MAX_BATCH_ITEMS", 100)
MAX_SMILES_LENGTH = get_env_int("MOLSIM_MAX_SMILES_LENGTH", 4096)
MAX_NAME_LENGTH = get_env_int("MOLSIM_MAX_NAME_LENGTH", 256)
MAX_BATCH_REQUEST_BYTES = get_env_int("MOLSIM_MAX_BATCH_REQUEST_BYTES", 512 * 1024)
MAX_EXPORT_IMAGE_BYTES = get_env_int("MOLSIM_MAX_EXPORT_IMAGE_BYTES", 2 * 1024 * 1024)
MAX_EXPORT_IMAGE_WIDTH = get_env_int("MOLSIM_MAX_EXPORT_IMAGE_WIDTH", 4096)
MAX_EXPORT_IMAGE_HEIGHT = get_env_int("MOLSIM_MAX_EXPORT_IMAGE_HEIGHT", 4096)
PUBCHEM_LOOKUP_TIMEOUT = get_env_int("MOLSIM_PUBCHEM_TIMEOUT", 5)
PUBCHEM_CACHE_TTL_SECONDS = get_env_int("MOLSIM_PUBCHEM_CACHE_TTL_SECONDS", 3600)
PUBCHEM_CACHE = {}

app.config["MAX_CONTENT_LENGTH"] = MAX_CONTENT_LENGTH
app.config["PROPAGATE_EXCEPTIONS"] = False


CORS(app, resources={r"/*": {"origins": parse_cors_origins()}})


def _error_response(status_code: int, message: str, field: Optional[str] = None, code: str = "invalid_request"):
    return jsonify(make_error_payload(status_code, message, field=field, code=code)), status_code


def _is_api_route(path: str) -> bool:
    return path in API_JSON_ROUTES


API_JSON_ROUTES = {"/compare", "/bulk-compare", "/report-preview", "/export-pdf", "/lookup-name", "/nitro-ra/cpca", "/nitro-ra/analyze", "/nitro-ra/deep-pk"}


@app.before_request
def enforce_json_api_routes():
    if request.method not in {"POST", "PUT", "PATCH"}:
        return None
    if request.path not in API_JSON_ROUTES:
        return None
    content_type = request.content_type or ""
    if "application/json" not in content_type.lower():
        return _error_response(415, "Content-Type inválido. Use application/json.", "Content-Type")
    return None


def _guard_request_size():
    limit = app.config.get("MAX_CONTENT_LENGTH", MAX_CONTENT_LENGTH)
    try:
        content_length = request.content_length
        payload_size = len(request.get_data(cache=True, as_text=False))
    except RequestEntityTooLarge:
        return _error_response(413, "Tamanho da requisição excede o limite permitido.", "body")

    if content_length is not None and content_length > limit:
        return _error_response(413, "Tamanho da requisição excede o limite permitido.", "body")
    if payload_size > limit:
        return _error_response(413, "Tamanho da requisição excede o limite permitido.", "body")
    return None


def _parse_json_object_payload():
    guard = _guard_request_size()
    if guard:
        return None, guard
    if not request.is_json:
        return None, _error_response(415, "Content-Type deve ser application/json.", "Content-Type")
    try:
        payload = request.get_json(silent=False)
    except Exception:
        return None, _error_response(400, "JSON inválido.", "body")
    if payload is None:
        return None, _error_response(400, "Corpo JSON vazio.", "body")
    if not isinstance(payload, dict):
        return None, _error_response(400, "JSON deve ser um objeto.", "body")
    if len(payload) == 0:
        return None, _error_response(400, "Objeto JSON vazio.", "body")
    return payload, None


def _validate_required_string(data: Dict[str, Any], field: str, max_length: int) -> Tuple[Optional[str], Optional[str]]:
    value = data.get(field)
    if not isinstance(value, str):
        return None, f"{field} deve ser string."
    normalized = value.strip()
    if not normalized:
        return None, f"{field} não pode ser vazio."
    if len(normalized) > max_length:
        return None, f"{field} excede o limite de {max_length} caracteres."
    return normalized, None


def _validate_optional_string(data: Dict[str, Any], field: str, max_length: int) -> Optional[str]:
    if field not in data or data.get(field) is None:
        return None
    value = data.get(field)
    if not isinstance(value, str):
        return f"{field} deve ser string."
    if len(value.strip()) > max_length:
        return f"{field} excede o limite de {max_length} caracteres."
    return None


def _validate_optional_bool(data: Dict[str, Any], field: str) -> Optional[str]:
    if field not in data:
        return None
    value = data.get(field)
    if not isinstance(value, bool):
        return f"{field} deve ser booleano."
    return None


def _validate_optional_finite_number(data: Dict[str, Any], field: str) -> Optional[str]:
    if field not in data:
        return None
    value = data.get(field)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return f"{field} deve ser número."
    if not math.isfinite(float(value)):
        return f"{field} deve ser número finito."
    return None


def _decode_png_base64(value: Any, field: str) -> Tuple[Optional[bytes], Optional[Tuple[Any, int]]]:
    if value is None:
        return None, None
    if not isinstance(value, str):
        return None, _error_response(400, f"{field} deve ser string base64.", field)
    candidate = value.strip()
    if not candidate:
        return None, None
    if candidate.startswith("data:"):
        prefix = "data:image/png;base64,"
        if not candidate.lower().startswith(prefix):
            return None, _error_response(400, f"{field} deve ser data URI PNG base64.", field)
        candidate = candidate[len(prefix):].strip()
    try:
        decoded = base64.b64decode(candidate, validate=True)
    except (binascii.Error, ValueError):
        return None, _error_response(400, f"{field} contém base64 inválido.", field)
    if len(decoded) > MAX_EXPORT_IMAGE_BYTES:
        return None, _error_response(400, f"{field} excede o limite de bytes.", field)
    try:
        width, height = _extract_png_dimensions(decoded)
    except ValueError:
        return None, _error_response(400, f"{field} não é um PNG válido.", field)
    if width > MAX_EXPORT_IMAGE_WIDTH or height > MAX_EXPORT_IMAGE_HEIGHT:
        return None, _error_response(400, f"{field} excede o limite de dimensões permitido.", field)
    return decoded, None


def _get_pubchem_cache_key(smiles: Optional[str]) -> Optional[str]:
    if not isinstance(smiles, str):
        return None
    normalized = smiles.strip()
    return normalized or None


def _get_cached_pubchem_name(smiles: Optional[str]) -> Optional[str]:
    key = _get_pubchem_cache_key(smiles)
    if not key:
        return None
    cached = PUBCHEM_CACHE.get(key)
    if not cached:
        return None
    cached_at = cached.get('cached_at')
    if not isinstance(cached_at, (int, float)):
        return None
    if time.monotonic() - cached_at > PUBCHEM_CACHE_TTL_SECONDS:
        PUBCHEM_CACHE.pop(key, None)
        return None
    value = cached.get('value')
    return value if isinstance(value, str) and value.strip() else None


def _store_pubchem_name(smiles: Optional[str], name: Optional[str]) -> None:
    key = _get_pubchem_cache_key(smiles)
    if not key or not isinstance(name, str):
        return
    title = name.strip()
    if not title:
        return
    PUBCHEM_CACHE[key] = {'value': title, 'cached_at': time.monotonic()}


def get_pubchem_name_for_smiles(smiles: Optional[str]) -> Optional[str]:
    """Look up a molecule name in PubChem when a user provides a SMILES without a name."""
    normalized = _get_pubchem_cache_key(smiles)
    if not normalized:
        return None

    cached_name = _get_cached_pubchem_name(normalized)
    if cached_name:
        return cached_name

    try:
        encoded_smiles = urllib_parse.quote(normalized, safe='')
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/property/Title/JSON"
        req = urllib_request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib_request.urlopen(req, timeout=PUBCHEM_LOOKUP_TIMEOUT) as response:
            payload = json.loads(response.read().decode('utf-8', errors='ignore'))
            properties = payload.get('PropertyTable', {}).get('Properties', [])
            if not isinstance(properties, list) or not properties:
                return None
            first_item = properties[0]
            if not isinstance(first_item, dict):
                return None
            title = first_item.get('Title')
            if not isinstance(title, str):
                return None
            title = title.strip()
            if title:
                _store_pubchem_name(normalized, title)
            return title or None
    except Exception:
        logger.warning("PubChem lookup failed for SMILES %s", normalized, exc_info=True)
        return None


def apply_name_fallbacks(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Fill blank molecule names from PubChem only when the user left the field empty."""
    if not isinstance(payload, dict):
        return payload

    for key in ('name_ref', 'name_test', 'molecule_name', 'name', 'reference_name', 'target_name'):
        if key not in payload:
            continue
        current_value = payload.get(key)
        if isinstance(current_value, str) and current_value.strip():
            continue

        if key == 'name_ref':
            smiles_key = 'smiles_ref'
        elif key == 'name_test':
            smiles_key = 'smiles_test'
        else:
            smiles_key = 'smiles'

        if smiles_key in payload:
            resolved_name = get_pubchem_name_for_smiles(payload.get(smiles_key))
            if resolved_name:
                payload[key] = resolved_name

    for nested_key in ('molecule_1', 'molecule_2'):
        if nested_key in payload and isinstance(payload[nested_key], dict):
            payload[nested_key] = apply_name_fallbacks(payload[nested_key])

    return payload


DEFAULT_IMAGE_MARGIN = 2
MOLECULE_BOX_GAP = 8
MOLECULE_BOX_MIN_HEIGHT = 45
MOLECULE_BOX_MAX_HEIGHT = 75


def _extract_png_dimensions(img_bytes: bytes) -> Tuple[int, int]:
    if len(img_bytes) < 24 or img_bytes[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("Imagem PNG inválida.")
    width, height = struct.unpack(">II", img_bytes[16:24])
    if width <= 0 or height <= 0:
        raise ValueError("Dimensões PNG inválidas.")
    return width, height


def _fit_inside_box(content_w: float, content_h: float, box_w: float, box_h: float) -> Tuple[float, float]:
    scale = min(box_w / content_w, box_h / content_h)
    return content_w * scale, content_h * scale


def validate_compare_request(data: Dict[str, Any]) -> Tuple[bool, str, Optional[str]]:
    smiles_ref, err = _validate_required_string(data, "smiles_ref", MAX_SMILES_LENGTH)
    if err:
        return False, err, "smiles_ref"
    smiles_test, err = _validate_required_string(data, "smiles_test", MAX_SMILES_LENGTH)
    if err:
        return False, err, "smiles_test"

    metric = data.get("metric")
    if not isinstance(metric, str) or not metric.strip():
        return False, "metric deve ser string não vazia.", "metric"
    metric = metric.strip()
    if not validate_metric(metric):
        return False, f"Métrica inválida: {metric}", "metric"

    fp_type = data.get("fp_type", "Morgan2")
    if not isinstance(fp_type, str) or not fp_type.strip():
        return False, "fp_type deve ser string não vazia.", "fp_type"
    fp_type = fp_type.strip()
    if not validate_fingerprint_type(fp_type):
        return False, f"Tipo de fingerprint inválido: {fp_type}", "fp_type"

    name_ref_error = _validate_optional_string(data, "name_ref", MAX_NAME_LENGTH)
    if name_ref_error:
        return False, name_ref_error, "name_ref"
    name_test_error = _validate_optional_string(data, "name_test", MAX_NAME_LENGTH)
    if name_test_error:
        return False, name_test_error, "name_test"
    show_map_error = _validate_optional_bool(data, "show_similarity_map")
    if show_map_error:
        return False, show_map_error, "show_similarity_map"
    show_logd_error = _validate_optional_bool(data, "show_logd")
    if show_logd_error:
        return False, show_logd_error, "show_logd"

    data["smiles_ref"] = smiles_ref
    data["smiles_test"] = smiles_test
    data["metric"] = metric
    data["fp_type"] = fp_type
    return True, "", None


def validate_bulk_compare_request(data: Dict[str, Any]) -> Tuple[bool, str, Optional[str]]:
    content_length = request.content_length
    if content_length is not None and content_length > MAX_BATCH_REQUEST_BYTES:
        return False, "Payload excede o limite configurado para comparação em lote.", "body"

    ref_smiles, err = _validate_required_string(data, "ref_smiles", MAX_SMILES_LENGTH)
    if err:
        return False, err, "ref_smiles"

    smiles_list = data.get("smiles_list")
    if not isinstance(smiles_list, list):
        return False, "smiles_list deve ser uma lista.", "smiles_list"
    if len(smiles_list) == 0:
        return False, "smiles_list deve ser uma lista não vazia.", "smiles_list"
    if len(smiles_list) > MAX_BATCH_ITEMS:
        return False, f"smiles_list excede o limite de {MAX_BATCH_ITEMS} itens.", "smiles_list"

    normalized_smiles: List[str] = []
    for index, smiles in enumerate(smiles_list):
        if not isinstance(smiles, str):
            return False, f"Item {index} de smiles_list deve ser string.", "smiles_list"
        normalized = smiles.strip()
        if not normalized:
            return False, f"Item {index} de smiles_list não pode ser vazio.", "smiles_list"
        if len(normalized) > MAX_SMILES_LENGTH:
            return False, f"Item {index} de smiles_list excede o limite de {MAX_SMILES_LENGTH} caracteres.", "smiles_list"
        normalized_smiles.append(normalized)

    names_list = data.get("names_list")
    if names_list is not None:
        if not isinstance(names_list, list):
            return False, "names_list deve ser uma lista.", "names_list"
        if len(names_list) != len(smiles_list):
            return False, "names_list deve ter o mesmo tamanho de smiles_list.", "names_list"
        for index, name in enumerate(names_list):
            if not isinstance(name, str):
                return False, f"Item {index} de names_list deve ser string.", "names_list"
            if len(name.strip()) > MAX_NAME_LENGTH:
                return False, f"Item {index} de names_list excede o limite de {MAX_NAME_LENGTH} caracteres.", "names_list"

    metric = data.get("metric")
    if not isinstance(metric, str) or not metric.strip():
        return False, "metric deve ser string não vazia.", "metric"
    metric = metric.strip()
    if not validate_metric(metric):
        return False, f"Métrica inválida: {metric}", "metric"

    fp_type = data.get("fp_type", "Morgan2")
    if not isinstance(fp_type, str) or not fp_type.strip():
        return False, "fp_type deve ser string não vazia.", "fp_type"
    fp_type = fp_type.strip()
    if not validate_fingerprint_type(fp_type):
        return False, f"Tipo de fingerprint inválido: {fp_type}", "fp_type"

    show_chemical_space_error = _validate_optional_bool(data, "show_chemical_space")
    if show_chemical_space_error:
        return False, show_chemical_space_error, "show_chemical_space"

    data["ref_smiles"] = ref_smiles
    data["smiles_list"] = normalized_smiles
    data["metric"] = metric
    data["fp_type"] = fp_type
    return True, "", None


@app.route('/')
def index():
    response = make_response(render_template('index.html'))
    response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '0'
    return response


@app.route('/lookup-name', methods=['POST'])
def api_lookup_name():
    try:
        data, parse_error = _parse_json_object_payload()
        if parse_error:
            return parse_error

        if data is None:
            return jsonify({"name": None}), 200

        if 'smiles' not in data:
            return jsonify({"name": None}), 200

        raw_smiles = data.get('smiles')
        if raw_smiles is None:
            return jsonify({"name": None}), 200
        if not isinstance(raw_smiles, str):
            return _error_response(400, "smiles deve ser string.", "smiles")

        smiles = raw_smiles.strip()
        if not smiles:
            return jsonify({"name": None}), 200
        if len(smiles) > MAX_SMILES_LENGTH:
            return _error_response(400, f"smiles excede o limite de {MAX_SMILES_LENGTH} caracteres.", "smiles")

        name = get_pubchem_name_for_smiles(smiles)
        return jsonify({"name": name}), 200
    except Exception:
        logger.error("Erro ao consultar nome no PubChem: %s", traceback.format_exc())
        return _error_response(500, "Erro ao consultar nome no PubChem.", code="internal_error")


@app.route('/compare', methods=['POST'])
def api_compare():
    try:
        data, parse_error = _parse_json_object_payload()
        if parse_error:
            return parse_error

        valid, error_msg, field = validate_compare_request(data)
        if not valid:
            return _error_response(400, error_msg, field)

        data = apply_name_fallbacks(data)

        result, error = run_pairwise_compare(data, compare)

        if error:
            return _error_response(400, error, "smiles")

        return jsonify(result), 200
    except Exception:
        logger.error("Erro não tratado em /compare: %s", traceback.format_exc())
        return _error_response(500, "Erro interno do servidor.", code="internal_error")


@app.route('/bulk-compare', methods=['POST'])
def api_bulk_compare():
    try:
        data, parse_error = _parse_json_object_payload()
        if parse_error:
            return parse_error

        valid, error_msg, field = validate_bulk_compare_request(data)
        if not valid:
            return _error_response(400, error_msg, field)

        data = apply_name_fallbacks(data)
        # O nome do alvo é opcional: quando ausente, buscar no PubChem por SMILES.
        # Se a busca não retornar resultado, usar um rótulo neutro em vez de "Referência".
        raw_ref_name = data.get("ref_name")
        if isinstance(raw_ref_name, str) and raw_ref_name.strip():
            data["ref_name"] = raw_ref_name.strip()
        else:
            data["ref_name"] = get_pubchem_name_for_smiles(data.get("ref_smiles")) or "Alvo"

        results, error = run_batch_compare(data, bulk_compare, get_pubchem_name_for_smiles)

        if error:
            return _error_response(400, error, "smiles")

        payload = {"results": results, "ref_name": data.get("ref_name", "Alvo")}
        reference_mol, reference_error = get_mol(data["ref_smiles"])
        if reference_mol is not None and not reference_error:
            reference_png = mol_to_png(reference_mol, size=220)
            payload["reference_png"] = base64.b64encode(reference_png).decode("utf-8") if reference_png else None
        if data.get("show_chemical_space", False):
            payload["chemical_space"] = build_chemical_space(
                data["ref_smiles"],
                data["smiles_list"],
                data.get("names_list"),
                results or [],
                "MACCS",
                "Tanimoto",
                ref_name=data.get("ref_name", "Alvo"),
            )
        return jsonify(payload), 200
    except Exception:
        logger.error("Erro não tratado em /bulk-compare: %s", traceback.format_exc())
        return _error_response(500, "Erro interno do servidor.", code="internal_error")


@app.route('/nitro-ra/cpca', methods=['POST'])
def api_nitro_ra_cpca():
    try:
        data, parse_error = _parse_json_object_payload()
        if parse_error:
            return parse_error

        smiles, smiles_error = _validate_required_string(data, "smiles", MAX_SMILES_LENGTH)
        if smiles_error:
            return _error_response(400, smiles_error, "smiles")

        mdd_error = _validate_optional_finite_number(data, "mdd_mg")
        if mdd_error:
            return _error_response(400, mdd_error, "mdd_mg")
        mdd_mg = data.get("mdd_mg")
        if mdd_mg is not None and float(mdd_mg) <= 0:
            return _error_response(400, "mdd_mg deve ser maior que zero.", "mdd_mg")

        result = calculate_cpca(smiles, mdd_mg=mdd_mg)
        if result.get("status") == "invalid_smiles":
            return jsonify(result), 400
        return jsonify(result), 200
    except Exception:
        logger.error("Erro não tratado em /nitro-ra/cpca: %s", traceback.format_exc())
        return _error_response(500, "Erro interno ao calcular cPCA.", code="internal_error")


@app.route('/nitro-ra/analyze', methods=['POST'])
def api_nitro_ra_analyze():
    try:
        data, parse_error = _parse_json_object_payload()
        if parse_error:
            return parse_error

        smiles, smiles_error = _validate_required_string(data, "smiles", MAX_SMILES_LENGTH)
        if smiles_error:
            return _error_response(400, smiles_error, "smiles")

        modules = data.get("modules")
        if not isinstance(modules, list) or not modules:
            return _error_response(400, "Informe ao menos um módulo Nitro.RA.", "modules")
        allowed_modules = {"cpca", "quantum", "metabolism", "nitrosamine_space", "deep_pk"}
        selected_modules = []
        for module in modules:
            if not isinstance(module, str) or module not in allowed_modules:
                return _error_response(400, "Módulo Nitro.RA inválido.", "modules")
            if module not in selected_modules:
                selected_modules.append(module)
        data = apply_name_fallbacks(data)
        molecule_name = data.get("name") if isinstance(data.get("name"), str) and data.get("name").strip() else None
        mdd_error = _validate_optional_finite_number(data, "mdd_mg")

        if mdd_error:
            return _error_response(400, mdd_error, "mdd_mg")
        mdd_mg = data.get("mdd_mg")
        if mdd_mg is not None and float(mdd_mg) <= 0:
            return _error_response(400, "mdd_mg deve ser maior que zero.", "mdd_mg")
        if "deep_pk" in selected_modules and "metabolism" not in selected_modules:
            return _error_response(400, "O módulo Deep-PK requer que Metabolism também seja selecionado.", "modules")

        results = {}
        if "cpca" in selected_modules:
            results["cpca"] = calculate_cpca(smiles, mdd_mg=mdd_mg)
        if "quantum" in selected_modules:
            results["quantum"] = evaluate_quantum(smiles)
        if "metabolism" in selected_modules:
            results["metabolism"] = evaluate_metabolism(smiles)
        if "nitrosamine_space" in selected_modules:
            results["nitrosamine_space"] = search_nitrosamine_space(smiles)
        if "deep_pk" in selected_modules:
            results["deep_pk"] = submit_deep_pk_metabolism(smiles)

        response = {
            "module": "nitro_ra",
            "status": "ok",
            "smiles": smiles,
            "name": molecule_name,
            "modules": selected_modules,
            "results": results,
        }
        invalid_modules = {"invalid_smiles", "invalid_input"}
        if any(results.get(module, {}).get("status") in invalid_modules for module in selected_modules):
            response["status"] = "invalid_smiles"
            return jsonify(response), 400
        return jsonify(response), 200
    except Exception:
        logger.error("Erro não tratado em /nitro-ra/analyze: %s", traceback.format_exc())
        return _error_response(500, "Erro interno ao executar módulos Nitro.RA.", code="internal_error")


@app.route('/nitro-ra/deep-pk', methods=['POST'])
def api_deep_pk_submit():
    try:
        data, parse_error = _parse_json_object_payload()
        if parse_error:
            return parse_error

        smiles, smiles_error = _validate_required_string(data, "smiles", MAX_SMILES_LENGTH)
        if smiles_error:
            return _error_response(400, smiles_error, "smiles")
        result = submit_deep_pk_metabolism(smiles)
        if result.get("status") == "invalid_smiles":
            return jsonify(result), 400
        return jsonify(result), 200
    except Exception:
        logger.error("Erro não tratado em /nitro-ra/deep-pk: %s", traceback.format_exc())
        return _error_response(500, "Erro interno ao consultar o Deep-PK.", code="internal_error")


@app.route('/nitro-ra/deep-pk/<job_id>', methods=['GET'])
def api_deep_pk_status(job_id):
    try:
        result = get_deep_pk_metabolism(job_id)
        if result.get("status") == "invalid_job_id":
            return jsonify(result), 400
        return jsonify(result), 200
    except Exception:
        logger.error("Erro não tratado em /nitro-ra/deep-pk/<job_id>: %s", traceback.format_exc())
        return _error_response(500, "Erro interno ao consultar o status Deep-PK.", code="internal_error")


@app.route('/report-preview', methods=['POST'])
def report_preview():
    try:
        data, parse_error = _parse_json_object_payload()
        if parse_error:
            return parse_error
        data = ensure_report_generated_at(data)
        return render_template("report_preview.html", report=data), 200
    except Exception:
        logger.error("Erro ao gerar preview do relatório: %s", traceback.format_exc())
        return _error_response(500, "Erro interno do servidor.", code="internal_error")


@app.route('/export-pdf', methods=['POST'])
def export_pdf():
    try:
        data, parse_error = _parse_json_object_payload()
        if parse_error:
            return parse_error

        data = ensure_report_generated_at(data)
        similarity_error = _validate_optional_finite_number(data, "similarity")
        if similarity_error:
            return _error_response(400, similarity_error, "similarity")

        for field in ("name_ref", "name_test", "fp_type", "metric", "classification", "generated_at"):
            field_error = _validate_optional_string(data, field, MAX_NAME_LENGTH if field in ("name_ref", "name_test") else 128)
            if field_error:
                return _error_response(400, field_error, field)

        properties = data.get("properties", [])
        if not isinstance(properties, list):
            return _error_response(400, "properties deve ser uma lista.", "properties")
        for index, item in enumerate(properties):
            if not isinstance(item, dict):
                return _error_response(400, f"Item {index} de properties deve ser objeto.", "properties")

        validated_images: Dict[str, Optional[bytes]] = {}
        for image_field in ("png_ref", "png_test", "fingerprint_ref_png", "fingerprint_test_png"):
            image_bytes, image_error = _decode_png_base64(data.get(image_field), image_field)
            if image_error:
                return image_error
            validated_images[image_field] = image_bytes

        pdf = ExportReportPDF()
        pdf.generated_at = data["generated_at"]
        pdf.set_auto_page_break(auto=False)
        pdf.set_margins(10, 10, 10)

        def add_image_fixed_box(pdf_obj, img_bytes: Optional[bytes], x, y, box_w, box_h, margin=DEFAULT_IMAGE_MARGIN):
            if not img_bytes:
                return
            tmp_path = None
            try:
                png_w, png_h = _extract_png_dimensions(img_bytes)
                drawable_w = box_w - (2 * margin)
                drawable_h = box_h - (2 * margin)
                fitted_w, fitted_h = _fit_inside_box(float(png_w), float(png_h), drawable_w, drawable_h)
                offset_x = (drawable_w - fitted_w) / 2
                offset_y = (drawable_h - fitted_h) / 2

                with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                    tmp.write(img_bytes)
                    tmp.flush()
                    tmp_path = tmp.name

                pdf_obj.rect(x, y, box_w, box_h)
                pdf_obj.image(
                    tmp_path,
                    x=x + margin + offset_x,
                    y=y + margin + offset_y,
                    w=fitted_w,
                    h=fitted_h
                )
            except Exception as e:
                logger.warning(f"Erro ao adicionar imagem ao PDF: {str(e)}")
            finally:
                if tmp_path and os.path.exists(tmp_path):
                    os.remove(tmp_path)

        # ============================================================
        # PÁGINA 1 — RESUMO + ESTRUTURAS MOLECULARES
        # ============================================================
        pdf.add_page()

        pdf.set_font("Helvetica", "B", 15)
        pdf.cell(0, 8, "Análise de Similaridade Molecular", align='C')
        pdf.ln(2)

        pdf.set_font("Helvetica", "", 9)
        pdf.cell(0, 6, f"Método: {data.get('fp_type', '-')} | Métrica: {data.get('metric', '-')}", align='C')
        pdf.ln(3)

        pdf.set_font("Helvetica", "B", 11)
        pdf.cell(0, 7, f"{data.get('name_ref', 'Referência')} vs {data.get('name_test', 'Teste')}")
        pdf.ln()
        pdf.set_font("Helvetica", "B", 9)
        similarity = float(data.get("similarity", 0.0))
        pdf.cell(0, 6, f"Similaridade: {similarity:.4f} ({data.get('classification', '-')})")
        pdf.ln(4)

        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(0, 6, "Estruturas Moleculares")
        pdf.ln(1)

        box_y = pdf.get_y()
        page_bottom_limit = pdf.h - pdf.b_margin
        footer_space = 8

        usable_width = pdf.w - pdf.l_margin - pdf.r_margin
        box_w = (usable_width - MOLECULE_BOX_GAP) / 2

        left_x = pdf.l_margin
        right_x = left_x + box_w + MOLECULE_BOX_GAP

        available_height = page_bottom_limit - box_y - footer_space
        if available_height < MOLECULE_BOX_MIN_HEIGHT:
            box_h = max(20, available_height)
        else:
            box_h = min(MOLECULE_BOX_MAX_HEIGHT, available_height)

        pdf.set_font("Helvetica", "B", 8)
        pdf.set_xy(left_x, box_y - 4)
        pdf.cell(box_w, 3, data.get('name_ref', 'Referência'), align='C')

        pdf.set_xy(right_x, box_y - 4)
        pdf.cell(box_w, 3, data.get('name_test', 'Teste'), align='C')

        add_image_fixed_box(pdf, validated_images.get('png_ref'), left_x, box_y, box_w, box_h, margin=2)
        add_image_fixed_box(pdf, validated_images.get('png_test'), right_x, box_y, box_w, box_h, margin=2)

        pdf.set_y(box_y + box_h + 6)

        pdf.set_font("Helvetica", "I", 7)
        pdf.cell(0, 5, "Página 1/2", align='R')
        pdf.ln()

        # ============================================================
        # PÁGINA 2 — FINGERPRINTS + PROPRIEDADES
        # ============================================================
        pdf.add_page()

        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Fingerprints")
        pdf.ln(2)

        fp_box_w = 85
        fp_box_h = 36
        fp_box_y = pdf.get_y()

        pdf.set_font("Helvetica", "B", 8)
        pdf.set_xy(left_x, fp_box_y - 4)
        pdf.cell(fp_box_w, 3, f"Fingerprint: {data.get('name_ref', 'Referência')}", align='C')

        pdf.set_xy(right_x, fp_box_y - 4)
        pdf.cell(fp_box_w, 3, f"Fingerprint: {data.get('name_test', 'Teste')}", align='C')

        add_image_fixed_box(pdf, validated_images.get('fingerprint_ref_png'), left_x, fp_box_y, fp_box_w, fp_box_h, margin=2)
        add_image_fixed_box(pdf, validated_images.get('fingerprint_test_png'), right_x, fp_box_y, fp_box_w, fp_box_h, margin=2)

        pdf.set_y(fp_box_y + fp_box_h + 8)

        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Propriedades Físico-Químicas")
        pdf.ln(2)

        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(60, 7, "Propriedade", border=1, align='C')
        pdf.cell(40, 7, "Referência", border=1, align='C')
        pdf.cell(40, 7, "Teste", border=1, align='C')
        pdf.cell(40, 7, "Diferença", border=1, align='C')
        pdf.ln()

        pdf.set_font("Helvetica", "", 9)
        for prop in data.get('properties', []):
            pdf.cell(60, 7, str(prop.get('Propriedade', ''))[:20], border=1)
            pdf.cell(40, 7, str(prop.get('Referência', '')), border=1, align='C')
            pdf.cell(40, 7, str(prop.get('Teste', '')), border=1, align='C')
            pdf.cell(40, 7, str(prop.get('Diferença', '')), border=1, align='C')
            pdf.ln()

        pdf.ln(6)

        pdf_output = pdf.output()
        pdf_bytes = pdf_output.encode('latin-1') if isinstance(pdf_output, str) else pdf_output

        return send_file(
            io.BytesIO(pdf_bytes),
            mimetype='application/pdf',
            as_attachment=True,
            download_name='molsim_resultado.pdf'
        )
    except Exception:
        logger.error("Erro ao gerar PDF: %s", traceback.format_exc())
        return _error_response(500, "Erro interno do servidor.", code="internal_error")


@app.after_request
def apply_security_headers(response):
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["X-Frame-Options"] = "DENY"
    response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
    response.headers["Permissions-Policy"] = "geolocation=(), microphone=(), camera=()"
    return response


@app.route('/healthz', methods=['GET', 'HEAD'])
def healthz():
    return jsonify({"status": "ok"}), 200


@app.errorhandler(Exception)
def unhandled_exception(error):
    if isinstance(error, HTTPException):
        raise error
    if _is_api_route(request.path):
        logger.error("Erro não tratado em %s: %s", request.path, traceback.format_exc())
        return _error_response(500, "Erro interno do servidor.", code="internal_error")
    logger.error("Erro não tratado em %s: %s", request.path, traceback.format_exc())
    return _error_response(500, "Erro interno do servidor.", code="internal_error")


@app.errorhandler(400)
def bad_request(error):
    return _error_response(400, "Requisição inválida.", "body", code="invalid_request")


@app.errorhandler(404)
def not_found(error):
    return _error_response(404, "Rota não encontrada.", code="not_found")


@app.errorhandler(405)
def method_not_allowed(error):
    return _error_response(405, "Método HTTP não permitido para esta rota.", "method", code="method_not_allowed")


@app.errorhandler(415)
def unsupported_media_type(error):
    return _error_response(415, "Content-Type inválido. Use application/json.", "Content-Type")


@app.errorhandler(RequestEntityTooLarge)
def request_too_large(error):
    return _error_response(413, "Tamanho da requisição excede o limite permitido.", "body")


@app.errorhandler(413)
def request_too_large_legacy(error):
    return _error_response(413, "Tamanho da requisição excede o limite permitido.", "body")


@app.errorhandler(500)
def internal_error(error):
    logger.error("Erro interno do servidor: %s", traceback.format_exc())
    return _error_response(500, "Erro interno do servidor.", code="internal_error")


if __name__ == '__main__':
    raise SystemExit("Use 'python main.py' para iniciar a aplicação.")