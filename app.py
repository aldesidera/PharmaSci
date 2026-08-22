"""
MolSim Final — Backend Flask (v3.3)
• Input validation
• Improved error handling
• Type hints
"""

import logging
import json
from typing import Tuple, Dict, Any, Optional
from urllib import request as urllib_request, parse as urllib_parse

from flask import Flask, render_template, request, jsonify, send_file, make_response
from flask_cors import CORS
from analysis import compare, bulk_compare, validate_fingerprint_type
from fpdf import FPDF
import io
import traceback
import base64
import tempfile
import struct
import os

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)


def get_pubchem_name_for_smiles(smiles: Optional[str]) -> Optional[str]:
    """Look up a molecule name in PubChem when a user provides a SMILES without a name."""
    if not smiles:
        return None

    try:
        encoded_smiles = urllib_parse.quote(smiles, safe='')
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/property/Title/JSON"
        req = urllib_request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib_request.urlopen(req, timeout=10) as response:
            payload = json.loads(response.read().decode('utf-8', errors='ignore'))
            properties = payload.get('PropertyTable', {}).get('Properties', [])
            if not properties:
                return None
            title = properties[0].get('Title')
            if not isinstance(title, str):
                return None
            title = title.strip()
            return title or None
    except Exception:
        logger.warning("PubChem lookup failed for SMILES %s", smiles, exc_info=True)
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


def validate_compare_request(data: Dict[str, Any]) -> Tuple[bool, str]:
    required_fields = ['smiles_ref', 'smiles_test', 'metric']
    for field in required_fields:
        if field not in data or not data[field]:
            return False, f"Campo obrigatório ausente: {field}"

    fp_type = data.get('fp_type', 'Morgan2')
    if not validate_fingerprint_type(fp_type):
        return False, f"Tipo de fingerprint inválido: {fp_type}"

    metric = data.get('metric')
    if metric not in ['Tanimoto', 'Dice']:
        return False, f"Métrica inválida: {metric}"

    return True, ""


def validate_bulk_compare_request(data: Dict[str, Any]) -> Tuple[bool, str]:
    required_fields = ['ref_smiles', 'smiles_list', 'metric']
    for field in required_fields:
        if field not in data or not data[field]:
            return False, f"Campo obrigatório ausente: {field}"

    if not isinstance(data['smiles_list'], list) or len(data['smiles_list']) == 0:
        return False, "smiles_list deve ser uma lista não vazia"

    fp_type = data.get('fp_type', 'Morgan2')
    if not validate_fingerprint_type(fp_type):
        return False, f"Tipo de fingerprint inválido: {fp_type}"

    metric = data.get('metric')
    if metric not in ['Tanimoto', 'Dice']:
        return False, f"Métrica inválida: {metric}"

    return True, ""


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
        data = request.get_json(silent=True) or {}
        smiles = str(data.get('smiles', '') or '').strip()
        if not smiles:
            return jsonify({"name": None}), 200
        name = get_pubchem_name_for_smiles(smiles)
        return jsonify({"name": name}), 200
    except Exception as e:
        logger.error(f"Erro ao consultar nome no PubChem: {traceback.format_exc()}")
        return jsonify({"name": None, "error": str(e)[:100]}), 500


@app.route('/compare', methods=['POST'])
def api_compare():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"error": "Requisição JSON vazia"}), 400

        data = apply_name_fallbacks(data)

        valid, error_msg = validate_compare_request(data)
        if not valid:
            return jsonify({"error": error_msg}), 400

        fp_type = data.get('fp_type', 'Morgan2')
        result, error = compare(
            data['smiles_ref'],
            data['smiles_test'],
            data.get('name_ref') or 'Molécula Referência',
            data.get('name_test') or 'Molécula Teste',
            fp_type,
            data['metric'],
            show_map=data.get('show_logd', True)
        )

        if error:
            return jsonify({"error": error}), 400

        return jsonify(result), 200
    except Exception as e:
        error_trace = traceback.format_exc()
        logger.error(f"Erro não tratado em /compare: {error_trace}")
        return jsonify({"error": f"Erro interno do servidor: {str(e)[:100]}"}), 500


@app.route('/bulk-compare', methods=['POST'])
def api_bulk_compare():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"error": "Requisição JSON vazia"}), 400

        data = apply_name_fallbacks(data)

        valid, error_msg = validate_bulk_compare_request(data)
        if not valid:
            return jsonify({"error": error_msg}), 400

        names_list = data.get('names_list')
        if isinstance(names_list, list):
            resolved_names = []
            for index, name in enumerate(names_list):
                candidate = name if isinstance(name, str) and name.strip() else None
                if candidate:
                    resolved_names.append(candidate)
                    continue
                smiles = data['smiles_list'][index] if index < len(data['smiles_list']) else None
                pubchem_name = get_pubchem_name_for_smiles(smiles)
                resolved_names.append(pubchem_name or f"Mol_{index + 1}")
            data['names_list'] = resolved_names

        results, error = bulk_compare(
            data['ref_smiles'],
            data['smiles_list'],
            data.get('names_list'),
            data.get('fp_type', 'Morgan2'),
            data['metric']
        )

        if error:
            return jsonify({"error": error}), 400

        return jsonify({"results": results}), 200
    except Exception as e:
        error_trace = traceback.format_exc()
        logger.error(f"Erro não tratado em /bulk-compare: {error_trace}")
        return jsonify({"error": f"Erro interno do servidor: {str(e)[:100]}"}), 500


@app.route('/report-preview', methods=['POST'])
def report_preview():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"error": "Requisição JSON vazia"}), 400
        return render_template("report_preview.html", report=data), 200
    except Exception as e:
        error_trace = traceback.format_exc()
        logger.error(f"Erro ao gerar preview do relatório: {error_trace}")
        return jsonify({"error": f"Erro ao gerar preview: {str(e)[:100]}"}), 500


@app.route('/export-pdf', methods=['POST'])
def export_pdf():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"error": "Requisição JSON vazia"}), 400

        pdf = FPDF()
        pdf.set_auto_page_break(auto=False)
        pdf.set_margins(10, 10, 10)

        def add_image_fixed_box(pdf_obj, img_b64, x, y, box_w, box_h, margin=DEFAULT_IMAGE_MARGIN):
            if not img_b64:
                return
            tmp_path = None
            try:
                img_bytes = base64.b64decode(img_b64)
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
        pdf.cell(0, 8, "Análise de Similaridade Molecular", ln=True, align='C')
        pdf.ln(2)

        pdf.set_font("Helvetica", "", 9)
        pdf.cell(0, 6, f"Método: {data.get('fp_type', '-')} | Métrica: {data.get('metric', '-')}", ln=True, align='C')
        pdf.ln(3)

        pdf.set_font("Helvetica", "B", 11)
        pdf.cell(0, 7, f"{data.get('name_ref', 'Referência')} vs {data.get('name_test', 'Teste')}", ln=True)
        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(0, 6, f"Similaridade: {data.get('similarity', 0.0):.4f} ({data.get('classification', '-')})", ln=True)
        pdf.ln(4)

        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(0, 6, "Estruturas Moleculares", ln=True)
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

        add_image_fixed_box(pdf, data.get('png_ref'), left_x, box_y, box_w, box_h, margin=2)
        add_image_fixed_box(pdf, data.get('png_test'), right_x, box_y, box_w, box_h, margin=2)

        pdf.set_y(box_y + box_h + 6)

        pdf.set_font("Helvetica", "I", 7)
        pdf.cell(0, 5, "Página 1/2", ln=True, align='R')

        # ============================================================
        # PÁGINA 2 — FINGERPRINTS + PROPRIEDADES
        # ============================================================
        pdf.add_page()

        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Fingerprints", ln=True)
        pdf.ln(2)

        fp_box_w = 85
        fp_box_h = 36
        fp_box_y = pdf.get_y()

        pdf.set_font("Helvetica", "B", 8)
        pdf.set_xy(left_x, fp_box_y - 4)
        pdf.cell(fp_box_w, 3, f"Fingerprint: {data.get('name_ref', 'Referência')}", align='C')

        pdf.set_xy(right_x, fp_box_y - 4)
        pdf.cell(fp_box_w, 3, f"Fingerprint: {data.get('name_test', 'Teste')}", align='C')

        add_image_fixed_box(pdf, data.get('fingerprint_ref_png'), left_x, fp_box_y, fp_box_w, fp_box_h, margin=2)
        add_image_fixed_box(pdf, data.get('fingerprint_test_png'), right_x, fp_box_y, fp_box_w, fp_box_h, margin=2)

        pdf.set_y(fp_box_y + fp_box_h + 8)

        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Propriedades Físico-Químicas", ln=True)
        pdf.ln(2)

        pdf.set_font("Helvetica", "B", 9)
        pdf.cell(60, 7, "Propriedade", border=1, align='C')
        pdf.cell(40, 7, "Referência", border=1, align='C')
        pdf.cell(40, 7, "Teste", border=1, align='C')
        pdf.cell(40, 7, "Diferença", border=1, ln=True, align='C')

        pdf.set_font("Helvetica", "", 9)
        for prop in data.get('properties', []):
            pdf.cell(60, 7, str(prop.get('Propriedade', ''))[:20], border=1)
            pdf.cell(40, 7, str(prop.get('Referência', '')), border=1, align='C')
            pdf.cell(40, 7, str(prop.get('Teste', '')), border=1, align='C')
            pdf.cell(40, 7, str(prop.get('Diferença', '')), border=1, ln=True, align='C')

        pdf.ln(6)
        pdf.set_font("Helvetica", "I", 7)
        pdf.cell(0, 5, f"Gerado em {data.get('generated_at', 'data não informada')}", ln=True, align='C')
        pdf.cell(0, 5, "Página 2/2", ln=True, align='R')

        pdf_output = pdf.output(dest='S')
        pdf_bytes = pdf_output.encode('latin-1') if isinstance(pdf_output, str) else pdf_output

        return send_file(
            io.BytesIO(pdf_bytes),
            mimetype='application/pdf',
            as_attachment=True,
            download_name='molsim_resultado.pdf'
        )
    except Exception as e:
        error_trace = traceback.format_exc()
        logger.error(f"Erro ao gerar PDF: {error_trace}")
        return jsonify({"error": f"Erro ao gerar PDF: {str(e)[:100]}"}), 500


@app.errorhandler(404)
def not_found(error):
    return jsonify({"error": "Rota não encontrada"}), 404


@app.errorhandler(500)
def internal_error(error):
    return jsonify({"error": "Erro interno do servidor"}), 500


if __name__ == '__main__':
    logger.info("🚀 MolSim v3.3 iniciando em http://localhost:5000")
    app.run(debug=False, host='0.0.0.0', port=5000)