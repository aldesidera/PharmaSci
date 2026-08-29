"""Optional Deep-PK CYP450 substrate/inhibitor integration.

The local Nitro.RA metabolism engine remains independent. This module only
adapts the public Deep-PK API into a small, explicit payload for CYP450
substrate and inhibitor endpoints; it does not generate metabolites.
"""

from __future__ import annotations

import json
import os
import re
import uuid
from typing import Any, Dict, List, Optional
from urllib import error as urllib_error
from urllib import parse as urllib_parse
from urllib import request as urllib_request

from rdkit import Chem


DEEP_PK_API_URL = os.getenv(
    "DEEP_PK_API_URL",
    "https://biosig.lab.uq.edu.au/deeppk/api/predict",
)
DEEP_PK_TIMEOUT_SECONDS = float(os.getenv("MOLSIM_DEEP_PK_TIMEOUT", "20"))
DEEP_PK_ISOFORMS = ("CYP1A2", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A4")
_DEEP_PK_JOB_ID_PATTERN = re.compile(r"^[A-Za-z0-9_.:-]{1,160}$")


class DeepPkError(RuntimeError):
    """Expected failure while communicating with or decoding Deep-PK."""

    def __init__(self, message: str, status_code: Optional[int] = None):
        super().__init__(message)
        self.status_code = status_code


def _canonicalize_smiles(smiles: str) -> Optional[str]:
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _decode_json_payload(raw: bytes) -> Dict[str, Any]:
    try:
        payload: Any = json.loads(raw.decode("utf-8"))
        if isinstance(payload, str):
            payload = json.loads(payload)
        if not isinstance(payload, dict):
            raise DeepPkError("A resposta do Deep-PK não possui um objeto JSON válido.")
        return payload
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise DeepPkError("A resposta do Deep-PK não pôde ser interpretada como JSON.") from exc


def _multipart_form(fields: Dict[str, str]) -> tuple[bytes, str]:
    boundary = f"----PharmaSciDeepPK{uuid.uuid4().hex}"
    chunks: List[bytes] = []
    for name, value in fields.items():
        chunks.extend(
            [
                f"--{boundary}\r\n".encode("utf-8"),
                f'Content-Disposition: form-data; name="{name}"\r\n\r\n'.encode("utf-8"),
                str(value).encode("utf-8"),
                b"\r\n",
            ]
        )
    chunks.append(f"--{boundary}--\r\n".encode("utf-8"))
    return b"".join(chunks), f"multipart/form-data; boundary={boundary}"


def _request_json(method: str, *, form: Optional[Dict[str, str]] = None, query: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    url = DEEP_PK_API_URL
    body = None
    headers = {"Accept": "application/json", "User-Agent": "PharmaSci-NitroRA/1.0"}
    if method.upper() == "GET":
        if query:
            body, content_type = _multipart_form(query)
            headers["Content-Type"] = content_type
    else:
        body, content_type = _multipart_form(form or {})
        headers["Content-Type"] = content_type

    req = urllib_request.Request(url, data=body, headers=headers, method=method.upper())
    try:
        with urllib_request.urlopen(req, timeout=DEEP_PK_TIMEOUT_SECONDS) as response:
            raw = response.read()
    except urllib_error.HTTPError as exc:
        raise DeepPkError(f"Deep-PK retornou HTTP {exc.code}.", status_code=exc.code) from exc
    except (urllib_error.URLError, TimeoutError, OSError) as exc:
        raise DeepPkError("Não foi possível conectar ao serviço Deep-PK.") from exc
    return _decode_json_payload(raw)


_DEEP_PK_PRESENTATIONS: Dict[str, Dict[str, str]] = {
    "not_selected": {
        "label": "Não selecionado",
        "tone": "neutral",
        "message": "O complemento externo não foi selecionado.",
        "action": "Selecione Deep-PK junto com Metabolism para iniciar a consulta.",
    },
    "running": {
        "label": "Processando",
        "tone": "info",
        "message": "A consulta foi enviada e aguarda o processamento dos endpoints CYP.",
        "action": "Aguarde a atualização do job ou consulte novamente.",
    },
    "ok": {
        "label": "Resultados disponíveis",
        "tone": "success",
        "message": "As previsões externas de substrato e inibição estão disponíveis.",
        "action": "Nenhuma ação adicional é necessária.",
    },
    "deep_pk_unavailable": {
        "label": "Indisponível",
        "tone": "error",
        "message": "Não foi possível acessar o serviço Deep-PK.",
        "action": "Verifique a conectividade e tente novamente mais tarde; o Metabolism local permanece válido.",
    },
    "deep_pk_error": {
        "label": "Resposta não processada",
        "tone": "error",
        "message": "O serviço retornou uma resposta que não pôde ser normalizada.",
        "action": "Repita a consulta; o resultado local permanece disponível.",
    },
    "deep_pk_timeout": {
        "label": "Tempo excedido",
        "tone": "warning",
        "message": "O Deep-PK excedeu o tempo de espera da consulta.",
        "action": "Nenhuma classificação externa foi apresentada; repita a consulta mais tarde.",
    },
    "invalid_smiles": {
        "label": "SMILES inválido",
        "tone": "error",
        "message": "O SMILES não pôde ser sanitizado localmente.",
        "action": "Corrija a estrutura antes de consultar o Deep-PK.",
    },
    "invalid_job_id": {
        "label": "Consulta inválida",
        "tone": "error",
        "message": "O identificador da consulta Deep-PK é inválido.",
        "action": "Inicie uma nova consulta a partir da análise local.",
    },
}


def normalize_deep_pk_message(result: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Add a stable, user-facing presentation contract to a Deep-PK payload."""

    payload = dict(result or {})
    status = str(payload.get("status") or "deep_pk_error")
    presentation = dict(_DEEP_PK_PRESENTATIONS.get(status, {
        "label": "Revisão necessária",
        "tone": "warning",
        "message": "O estado retornado pelo Deep-PK requer revisão.",
        "action": "Verifique a consulta e preserve o resultado local do Metabolism.",
    }))
    presentation["status"] = status
    payload["presentation"] = presentation
    return payload


def _base_result(smiles: str, status: str, message: str) -> Dict[str, Any]:
    result = {
        "module": "deep_pk_metabolism",
        "provider": "Deep-PK",
        "prediction_type": "metabolism",
        "status": status,
        "message": message,
        "smiles": smiles,
        "isoforms": [],
        "warnings": [
            "Os resultados Deep-PK são previsões externas de endpoints CYP e não substituem confirmação experimental, analítica ou regulatória."
        ],
    }
    return normalize_deep_pk_message(result)


def submit_deep_pk_metabolism(smiles: str) -> Dict[str, Any]:
    """Submit one canonical SMILES to Deep-PK's metabolism API."""

    canonical = _canonicalize_smiles(smiles)
    if canonical is None:
        result = _base_result(smiles if isinstance(smiles, str) else "", "invalid_smiles", "SMILES inválido ou não sanitizável.")
        result["warnings"] = ["O Deep-PK não foi consultado porque o SMILES não pôde ser sanitizado localmente."]
        return result

    try:
        payload = _request_json(
            "POST",
            form={"smiles": canonical, "pred_type": "metabolism"},
        )
    except DeepPkError as exc:
        result = _base_result(canonical, "deep_pk_unavailable", str(exc))
        result["http_status"] = exc.status_code
        return result

    job_id = payload.get("job_id")
    if not isinstance(job_id, str) or not _DEEP_PK_JOB_ID_PATTERN.fullmatch(job_id):
        result = _base_result(canonical, "deep_pk_error", "O Deep-PK não retornou um job_id válido.")
        return result

    return {
        **_base_result(canonical, "running", "Consulta Deep-PK enviada; aguardando os endpoints CYP."),
        "job_id": job_id,
    }


def _normalized_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def _find_endpoint_value(entry: Dict[str, Any], isoform: str, relation: str, field: str) -> Any:
    isoform_key = _normalized_key(isoform)
    relation_key = _normalized_key(relation)
    field_key = _normalized_key(field)
    for key, value in entry.items():
        normalized = _normalized_key(key)
        if isoform_key in normalized and relation_key in normalized and field_key in normalized:
            return value
    return None


_PREDICTION_LABELS = {
    "Substrate": "Substrato",
    "Non-Substrate": "Não substrato",
    "Inhibitor": "Inibidor",
    "Non-Inhibitor": "Não inibidor",
}

_INTERPRETATION_LABELS = {
    "High Non-Substrative Activity": "Alta atividade não substrativa",
    "High Non-Inhibition": "Alta ausência de inibição",
    "High Substrate Activity": "Alta atividade substrativa",
    "High Inhibition": "Alta atividade inibitória",
}


def _translate_prediction(value: Any) -> Optional[str]:
    if value is None:
        return None
    return _PREDICTION_LABELS.get(str(value), str(value))


def _translate_interpretation(value: Any) -> Optional[str]:
    if value is None:
        return None
    return _INTERPRETATION_LABELS.get(str(value), str(value))


def _probability_percent(value: Any) -> Optional[float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    if 0 <= value <= 1:
        return round(float(value) * 100, 1)
    return round(float(value), 1)


def _extract_isoforms(entry: Dict[str, Any]) -> List[Dict[str, Any]]:
    isoforms: List[Dict[str, Any]] = []
    for isoform in DEEP_PK_ISOFORMS:
        row: Dict[str, Any] = {"isoform": isoform}
        for relation in ("substrate", "inhibitor"):
            probability = _find_endpoint_value(entry, isoform, relation, "probability")
            row[relation] = {
                "prediction": _find_endpoint_value(entry, isoform, relation, "predictions"),
                "prediction_label": _translate_prediction(_find_endpoint_value(entry, isoform, relation, "predictions")),
                "probability": probability,
                "probability_percent": _probability_percent(probability),
                "interpretation": _find_endpoint_value(entry, isoform, relation, "interpretation"),
                "interpretation_label": _translate_interpretation(_find_endpoint_value(entry, isoform, relation, "interpretation")),
            }
        isoforms.append(row)
    return isoforms


def _first_molecule(payload: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    if "0" in payload and isinstance(payload["0"], dict):
        return payload["0"]
    for value in payload.values():
        if isinstance(value, dict) and any("prediction" in _normalized_key(key) for key in value):
            return value
    return None


def get_deep_pk_metabolism(job_id: str) -> Dict[str, Any]:
    """Retrieve and normalize one Deep-PK metabolism job."""

    if not isinstance(job_id, str) or not _DEEP_PK_JOB_ID_PATTERN.fullmatch(job_id):
        return _base_result("", "invalid_job_id", "job_id Deep-PK inválido.")

    try:
        payload = _request_json("GET", query={"job_id": job_id})
    except DeepPkError as exc:
        result = _base_result("", "deep_pk_unavailable", str(exc))
        result["job_id"] = job_id
        result["http_status"] = exc.status_code
        return result

    if str(payload.get("status", "")).lower() in {"running", "queued", "pending"}:
        return {
            **_base_result("", "running", "O Deep-PK ainda está processando o job."),
            "job_id": job_id,
        }

    entry = _first_molecule(payload)
    if entry is None:
        result = _base_result("", "deep_pk_error", "O Deep-PK não retornou endpoints de metabolismo reconhecíveis.")
        result["job_id"] = job_id
        return result

    isoforms = _extract_isoforms(entry)
    return {
        **_base_result(
            str(entry.get("SMILES", "")),
            "ok",
            "Endpoints Deep-PK de substrato e inibição CYP recuperados.",
        ),
        "job_id": job_id,
        "isoforms": isoforms,
    }
