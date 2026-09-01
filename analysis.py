"""
MolSim Final v3.3 — Análise Molecular Avançada
• Visualização 3D das moléculas
• Layout mais compacto e profissional
• Similarity Map aprimorado
• Suporte a PDF com PNGs
• Type hints, constantes e logging (v3.3)
"""

import logging
from functools import lru_cache
from typing import Tuple, Dict, Optional, List, Any
import base64
import io
import math
import re
import numpy as np
from PIL import Image, ImageChops, ImageDraw, ImageFont

from chemo_suite.core.chemical_space import (
    DESCRIPTOR_KEYS,
    calculate_multimodal_space,
    classical_mds,
    descriptor_vector,
    normalized_stress,
    select_nearest_indices,
)
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors, Crippen, MACCSkeys
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.Chem import Mol
from rdkit.Geometry import Point2D
from rdkit.Geometry import Point2D

# ============================================================================
# CONSTANTES
# ============================================================================

FINGERPRINT_BITS = 2048
MORGAN_RADIUS = 2
SIMILARITY_MAP_FINGERPRINT = "Morgan2"
SIMILARITY_MAP_METRIC = "Tanimoto"
MACC_KEYS = 166

LOGD_PH_COEFFICIENT = 0.12
NEUTRAL_PH = 7
PH_RANGE_START = 0
PH_RANGE_END = 15

SVG_SIZE_MOLECULE = 400
SVG_SIZE_SIMILARITY_MAP = 800
PNG_SIZE_MOLECULE = 400
PNG_SIZE_SIMILARITY_MAP = 800
FINGERPRINT_IMAGE_SIZE = 300

SIMILARITY_THRESHOLDS: Dict[str, Dict[str, float]] = {
    "Tanimoto": {
        "Muito Alta": 0.85,
        "Alta": 0.65,
        "Moderada": 0.40,
        "Baixa": 0.20,
        "Muito Baixa": 0.0,
    },
    "Dice": {
        "Muito Alta": 0.90,
        "Alta": 0.70,
        "Moderada": 0.50,
        "Baixa": 0.30,
        "Muito Baixa": 0.0,
    },
}

VALID_FINGERPRINT_TYPES = {"Morgan2", "RDKit", "MACCS"}
VALID_METRICS = {"Tanimoto", "Dice"}

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def validate_fingerprint_type(fp_type: str) -> bool:
    return fp_type in VALID_FINGERPRINT_TYPES


def validate_metric(metric: str) -> bool:
    return metric in VALID_METRICS


def _is_valid_molecule(mol: Optional[Mol]) -> bool:
    return isinstance(mol, Mol) and mol.GetNumAtoms() > 0


def _normalize_smiles(smiles: Optional[str]) -> str:
    if not isinstance(smiles, str):
        return ""
    return smiles.strip()


def get_mol(smiles: str) -> Tuple[Optional[Mol], Optional[str]]:
    normalized = _normalize_smiles(smiles)
    if not normalized:
        return None, "SMILES inválido – valor vazio."

    try:
        mol = Chem.MolFromSmiles(normalized)
        if mol is None:
            return None, "SMILES inválido – não foi possível parsear."

        Chem.SanitizeMol(mol)
        AllChem.Compute2DCoords(mol)
        logger.info(f"Molécula carregada com sucesso: {normalized[:50]}")
        return mol, None
    except Exception as e:
        err_str = str(e).lower()
        if "valence" in err_str:
            error_msg = "Erro de Valência: átomo com número inválido de ligações."
        elif any(x in err_str for x in ["ring", "kekul", "aromatic"]):
            error_msg = "Erro de anel ou Kekulização: estrutura aromática inválida."
        else:
            error_msg = f"Erro ao processar SMILES: {str(e)[:100]}"
        logger.error(f"Erro ao carregar molécula: {error_msg}")
        return None, error_msg


def estimate_pka(mol: Optional[Mol]) -> Optional[Dict[str, Optional[float]]]:
    if not mol:
        return None

    try:
        pka_acid = None
        pka_basic = None

        sulfonic_acid = Chem.MolFromSmarts("S(=O)(=O)([OH])")
        carboxylic_acid = Chem.MolFromSmarts("[#6](=O)[OH]")
        aromatic_phenol = Chem.MolFromSmarts("[c][OH]")
        aliphatic_alcohol = Chem.MolFromSmarts("[CX4][OH]")
        aromatic_amine = Chem.MolFromSmarts("[c][N;H2,H1,H0]")
        primary_amine = Chem.MolFromSmarts("[N;H2;!$(N~[C]=[O]);!$(N~[c])]")
        secondary_amine = Chem.MolFromSmarts("[N;H1;!$(N~[C]=[O]);!$(N~[c])]")
        tertiary_amine = Chem.MolFromSmarts("[N;H0;!$(N~[C]=[O]);!$(N~[c])]")
        amide = Chem.MolFromSmarts("[NX3](=O)")
        imidazole = Chem.MolFromSmarts("[nH]1cccc1")
        pyridine = Chem.MolFromSmarts("[nH0]1ccccc1")
        amidine = Chem.MolFromSmarts("[C](=[NH2])[NH2]")
        guanidine = Chem.MolFromSmarts("[N]=C(N)N")

        if mol.HasSubstructMatch(sulfonic_acid):
            pka_acid = 1.5
        elif mol.HasSubstructMatch(carboxylic_acid):
            pka_acid = 4.5
        elif mol.HasSubstructMatch(aromatic_phenol):
            pka_acid = 10.0
        elif mol.HasSubstructMatch(aliphatic_alcohol):
            pka_acid = 15.0

        if mol.HasSubstructMatch(guanidine):
            pka_basic = 12.5
        elif mol.HasSubstructMatch(amidine):
            pka_basic = 11.5
        elif mol.HasSubstructMatch(imidazole):
            pka_basic = 6.8
        elif mol.HasSubstructMatch(pyridine):
            pka_basic = 5.2
        elif mol.HasSubstructMatch(aromatic_amine):
            pka_basic = 5.0
        elif mol.HasSubstructMatch(secondary_amine):
            pka_basic = 10.5
        elif mol.HasSubstructMatch(primary_amine):
            pka_basic = 9.5
        elif mol.HasSubstructMatch(tertiary_amine):
            pka_basic = 9.0
        elif mol.HasSubstructMatch(amide):
            pka_basic = 15.0

        if pka_acid is None and pka_basic is None:
            return None

        return {"pKa ácido": round(pka_acid, 2) if pka_acid is not None else None,
                "pKa básico": round(pka_basic, 2) if pka_basic is not None else None}
    except Exception as e:
        logger.warning(f"Não foi possível estimar pKa: {str(e)}")
        return None


def _estimate_water_solubility_ph7(mol: Mol, mw: float, logp: float, rotatable_bonds: int, pka: Optional[Dict[str, Any]]) -> float:
    """Estima solubilidade em água a pH 7 em mg/L; não substitui medida experimental."""
    heavy_atoms = max(mol.GetNumHeavyAtoms(), 1)
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    aromatic_proportion = aromatic_atoms / heavy_atoms
    log_s_neutral = (
        0.16
        - 1.5 * logp
        - 0.01 * (mw - 40.0)
        + 0.066 * rotatable_bonds
        + 0.066 * aromatic_proportion
    )
    ionization_factor = 1.0
    if pka:
        pka_acid = pka.get("pKa ácido")
        pka_basic = pka.get("pKa básico")
        if isinstance(pka_acid, (int, float)):
            ionization_factor = 1.0 + 10 ** (7.0 - float(pka_acid))
        elif isinstance(pka_basic, (int, float)):
            ionization_factor = 1.0 + 10 ** (float(pka_basic) - 7.0)
    log_s_ph7 = log_s_neutral + math.log10(max(ionization_factor, 1.0))
    solubility_mg_l = (10 ** min(log_s_ph7, 12.0)) * mw * 1000.0
    return round(max(solubility_mg_l, 0.0), 2)


@lru_cache(maxsize=512)
def _compute_property_payload(smiles_key: str) -> Optional[Dict[str, Any]]:
    try:
        mol = Chem.MolFromSmiles(smiles_key)
        if mol is None:
            return None
        mol_work = Chem.RemoveHs(mol)
        mw = Descriptors.MolWt(mol_work)
        logp = Crippen.MolLogP(mol_work)
        tpsa = Descriptors.TPSA(mol_work)
        hbd = Descriptors.NumHDonors(mol_work)
        hba = Descriptors.NumHAcceptors(mol_work)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol_work)
        pka = estimate_pka(mol_work)
        water_solubility_ph7 = _estimate_water_solubility_ph7(mol_work, mw, logp, rotatable_bonds, pka)

        return {
            "Massa Molecular (g/mol)": round(mw, 2),
            "Coeficiente de Partição (LogP)": round(logp, 3),
            "Área de Superfície Polar (Å²)": round(tpsa, 2),
            "pKa ácido": round(pka["pKa ácido"], 2) if pka and pka.get("pKa ácido") is not None else "N/A",
            "pKa básico": round(pka["pKa básico"], 2) if pka and pka.get("pKa básico") is not None else "N/A",
            "Doadores de H (HBD)": int(hbd),
            "Receptores de H (HBA)": int(hba),
            "Ligações rotacionáveis (RotB)": int(rotatable_bonds),
            "Solubilidade em água (pH 7, estimada) (mg/L)": water_solubility_ph7,
        }
    except Exception:
        return None


@lru_cache(maxsize=512)
def _compute_fingerprint(smiles_key: str, fp_type: str) -> Optional[Any]:
    try:
        mol = Chem.MolFromSmiles(smiles_key)
        if mol is None:
            return None
        if not isinstance(mol, Mol):
            return None
        mol_work = Chem.RemoveHs(mol)
        if mol_work is None:
            return None

        if fp_type == "Morgan2":
            generator = getattr(rdFingerprintGenerator, "GetMorganGenerator", None)
            if callable(generator):
                morgan = generator(radius=MORGAN_RADIUS, fpSize=FINGERPRINT_BITS, includeChirality=True)
                return morgan.GetFingerprint(mol_work)
            return rdMolDescriptors.GetMorganFingerprintAsBitVect(
                mol_work, radius=MORGAN_RADIUS, nBits=FINGERPRINT_BITS
            )
        if fp_type == "RDKit":
            return Chem.RDKFingerprint(mol_work)
        if fp_type == "MACCS":
            return MACCSkeys.GenMACCSKeys(mol_work)
        return None
    except Exception:
        return None


def get_properties(mol: Optional[Mol]) -> Optional[Dict[str, Any]]:
    if not mol:
        logger.warning("Tentativa de calcular propriedades com mol=None")
        return None

    try:
        smiles_key = Chem.MolToSmiles(mol, canonical=True)
        cached = _compute_property_payload(smiles_key)
        if cached is not None:
            return cached
        logger.warning("Propriedades não cacheadas para SMILES: %s", smiles_key)
        return None
    except Exception as e:
        logger.error(f"Erro ao calcular propriedades: {str(e)}")
        return None


def get_fingerprint(mol: Optional[Mol], fp_type: str) -> Optional[Any]:
    if not mol:
        logger.warning("Tentativa de gerar fingerprint com mol=None")
        return None
    if not validate_fingerprint_type(fp_type):
        logger.error(f"Tipo de fingerprint inválido: {fp_type}")
        return None

    try:
        smiles_key = Chem.MolToSmiles(mol, canonical=True)
        cached = _compute_fingerprint(smiles_key, fp_type)
        if cached is not None:
            return cached
        logger.warning("Fingerprint não cacheado para %s / %s", smiles_key, fp_type)
        return None
    except Exception as e:
        logger.error(f"Erro ao gerar fingerprint {fp_type}: {str(e)}")
        return None


def fingerprint_to_png(mol: Optional[Mol], fp_type: str, size: int = FINGERPRINT_IMAGE_SIZE) -> Optional[str]:
    """Gera uma imagem textual simples do fingerprint em PNG/base64."""
    try:
        fp = get_fingerprint(mol, fp_type)
        if fp is None:
            return None

        on_bits = list(fp.GetOnBits()) if hasattr(fp, "GetOnBits") else []
        bit_count = len(on_bits)

        canvas = Image.new("RGB", (size, 120), "white")
        draw = ImageDraw.Draw(canvas)
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 18)
            small_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14)
        except OSError:
            font = ImageFont.load_default()
            small_font = font
        draw.text((12, 12), f"Fingerprint: {fp_type}", fill="#172033", font=font)
        draw.text((12, 48), f"Bits ativos: {bit_count}", fill="#334155", font=small_font)
        draw.text((12, 76), f"Tamanho: {size}px", fill="#64748b", font=small_font)
        output = io.BytesIO()
        canvas.save(output, format="PNG", optimize=True)
        return base64.b64encode(output.getvalue()).decode("utf-8")
    except Exception as e:
        logger.error(f"Erro ao gerar imagem do fingerprint: {str(e)}")
        return None


def calc_similarity(fp1: Any, fp2: Any, metric: str) -> float:
    if not validate_metric(metric):
        logger.error(f"Métrica inválida: {metric}")
        return 0.0

    try:
        if metric == "Dice":
            return round(DataStructs.DiceSimilarity(fp1, fp2), 4)
        return round(DataStructs.TanimotoSimilarity(fp1, fp2), 4)
    except Exception as e:
        logger.error(f"Erro ao calcular similaridade: {str(e)}")
        return 0.0


def classify_similarity(similarity: float, metric: str) -> str:
    thresholds = SIMILARITY_THRESHOLDS.get(metric, SIMILARITY_THRESHOLDS["Tanimoto"])
    for classification, threshold in sorted(thresholds.items(), key=lambda x: x[1], reverse=True):
        if similarity >= threshold:
            return classification
    return "Muito Baixa"


def _mol_cache_key(mol: Optional[Mol]) -> Optional[str]:
    if not mol:
        return None
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


@lru_cache(maxsize=256)
def _mol_to_svg_cached(smiles_key: str, size: int) -> str:
    try:
        mol = Chem.MolFromSmiles(smiles_key)
        if mol is None or not _is_valid_molecule(mol):
            return ""
        mol_no_h = Chem.RemoveHs(mol)
        if not _is_valid_molecule(mol_no_h):
            return ""
        drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
        options = drawer.drawOptions()
        options.addStereoAnnotation = True
        options.prepareMolsBeforeDrawing = True
        options.bondLineWidth = 2.5
        options.minFontSize = 14
        options.annotationFontScale = 0.85
        drawer.DrawMolecule(mol_no_h)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace(f'width="{size}px"', 'width="100%"').replace(f'height="{size}px"', 'height="100%"')
    except Exception as e:
        logger.error(f"Erro ao gerar SVG cacheado: {str(e)}")
        return ""


@lru_cache(maxsize=256)
def _trim_png_whitespace(encoded: bytes, padding: int = 12) -> bytes:
    """Trim empty canvas around a molecule while preserving a small readable margin."""
    try:
        image = Image.open(io.BytesIO(encoded)).convert("RGBA")
        white = Image.new("RGBA", image.size, (255, 255, 255, 255))
        composited = Image.alpha_composite(white, image).convert("RGB")
        background = Image.new("RGB", composited.size, (255, 255, 255))
        diff = ImageChops.difference(composited, background)
        bbox = diff.getbbox()
        if not bbox:
            return encoded
        left = max(0, bbox[0] - padding)
        top = max(0, bbox[1] - padding)
        right = min(image.width, bbox[2] + padding)
        bottom = min(image.height, bbox[3] + padding)
        cropped = composited.crop((left, top, right, bottom))
        output = io.BytesIO()
        cropped.save(output, format="PNG", optimize=True)
        return output.getvalue()
    except Exception as exc:
        logger.warning(f"Não foi possível recortar margens do PNG molecular: {exc}")
        return encoded


@lru_cache(maxsize=256)
def _mol_to_png_cached(smiles_key: str, size: int) -> bytes:
    try:
        mol = Chem.MolFromSmiles(smiles_key)
        if mol is None or not _is_valid_molecule(mol):
            return b""
        mol_no_h = Chem.RemoveHs(mol)
        if not _is_valid_molecule(mol_no_h):
            return b""
        canvas_size = max(size, 500)
        drawer = rdMolDraw2D.MolDraw2DCairo(canvas_size, canvas_size)
        options = drawer.drawOptions()
        options.addStereoAnnotation = True
        options.prepareMolsBeforeDrawing = True
        options.bondLineWidth = 2.2
        options.minFontSize = 18
        options.annotationFontScale = 0.9
        options.fixedBondLength = 30
        drawer.DrawMolecule(mol_no_h)
        drawer.FinishDrawing()
        return _trim_png_whitespace(drawer.GetDrawingText())
    except Exception as e:
        logger.error(f"Erro ao gerar PNG cacheado: {str(e)}")
        return b""


def mol_to_svg(mol: Optional[Mol], size: int = SVG_SIZE_MOLECULE) -> str:
    try:
        if not _is_valid_molecule(mol):
            return ""
        smiles_key = _mol_cache_key(mol)
        if not smiles_key:
            return ""
        return _mol_to_svg_cached(smiles_key, size)
    except Exception as e:
        logger.error(f"Erro ao gerar SVG: {str(e)}")
        return ""


def _get_similarity_map_drawer(size: int):
    canvas_height = max(480, int(size * 0.625))
    drawer = rdMolDraw2D.MolDraw2DSVG(size, canvas_height)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 1.6
    opts.minFontSize = 14
    opts.padding = 0.22
    opts.additionalAtomLabelPadding = 0.06
    return drawer


def _similarity_map_fingerprint(mol, idx):
    return SimilarityMaps.GetMorganFingerprint(
        mol,
        idx,
        radius=MORGAN_RADIUS,
        useFeatures=True,
        useChirality=True,
    )


def _crop_similarity_map_svg(svg: str, canvas_width: int, canvas_height: int) -> str:
    """Ajusta o viewBox ao conteúdo do heatmap sem remover a margem de segurança."""
    try:
        bounds = [float(canvas_width), float(canvas_height), 0.0, 0.0]
        found = False

        def include(x: float, y: float, radius: float = 0.0) -> None:
            nonlocal found
            bounds[0] = min(bounds[0], x - radius)
            bounds[1] = min(bounds[1], y - radius)
            bounds[2] = max(bounds[2], x + radius)
            bounds[3] = max(bounds[3], y + radius)
            found = True

        number_pattern = r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?"
        for path_data in re.findall(r"<path\b[^>]*\bd=['\"]([^'\"]+)['\"]", svg):
            values = [float(value) for value in re.findall(number_pattern, path_data)]
            for index in range(0, len(values) - 1, 2):
                include(values[index], values[index + 1])

        for cx, cy, radius in re.findall(
            rf"<circle\b[^>]*\bcx=['\"]({number_pattern})['\"][^>]*\bcy=['\"]({number_pattern})['\"][^>]*\br=['\"]({number_pattern})['\"]",
            svg,
        ):
            include(float(cx), float(cy), float(radius))

        for x1, y1, x2, y2 in re.findall(
            rf"<line\b[^>]*\bx1=['\"]({number_pattern})['\"][^>]*\by1=['\"]({number_pattern})['\"][^>]*\bx2=['\"]({number_pattern})['\"][^>]*\by2=['\"]({number_pattern})['\"]",
            svg,
        ):
            include(float(x1), float(y1))
            include(float(x2), float(y2))

        for x, y in re.findall(rf"<text\b[^>]*\bx=['\"]({number_pattern})['\"][^>]*\by=['\"]({number_pattern})['\"]", svg):
            include(float(x), float(y), 8.0)

        if not found:
            return svg

        min_x, min_y, max_x, max_y = bounds
        content_width = max(max_x - min_x, 1.0)
        content_height = max(max_y - min_y, 1.0)
        margin = max(12.0, min(28.0, 0.07 * max(content_width, content_height)))
        min_x = max(0.0, min_x - margin)
        min_y = max(0.0, min_y - margin)
        max_x = min(float(canvas_width), max_x + margin)
        max_y = min(float(canvas_height), max_y + margin)
        cropped_width = max(max_x - min_x, 1.0)
        cropped_height = max(max_y - min_y, 1.0)
        return re.sub(
            r"viewBox=['\"][^'\"]+['\"]",
            f'viewBox="{min_x:.2f} {min_y:.2f} {cropped_width:.2f} {cropped_height:.2f}"',
            svg,
            count=1,
        )
    except Exception as exc:
        logger.warning(f"Não foi possível recortar o viewBox do heatmap: {exc}")
        return svg


@lru_cache(maxsize=128)
def _similarity_map_cached(ref_key: str, test_key: str, size: int) -> Optional[str]:
    try:
        mol_ref = Chem.MolFromSmiles(ref_key)
        mol_test = Chem.MolFromSmiles(test_key)
        if mol_ref is None or mol_test is None or not _is_valid_molecule(mol_ref) or not _is_valid_molecule(mol_test):
            return None

        mol_ref_no_h = Chem.RemoveHs(mol_ref)
        mol_test_no_h = Chem.RemoveHs(mol_test)
        if not _is_valid_molecule(mol_ref_no_h) or not _is_valid_molecule(mol_test_no_h):
            return None

        drawer = _get_similarity_map_drawer(size)

        try:
            SimilarityMaps.GetSimilarityMapForFingerprint(
                mol_ref_no_h,
                mol_test_no_h,
                _similarity_map_fingerprint,
                metric=DataStructs.TanimotoSimilarity,
                draw2d=drawer,
            )
        except Exception as inner_e:
            logger.warning(f"SimilarityMaps falhou: {str(inner_e)}")
            return None

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        if '<svg ' in svg and 'preserveAspectRatio=' not in svg:
            svg = svg.replace('<svg ', '<svg preserveAspectRatio="xMidYMid meet" ', 1)
        svg = re.sub(r"width=['\"][^'\"]+['\"]", 'width="100%"', svg, count=1)
        svg = re.sub(r"height=['\"][^'\"]+['\"]", 'height="100%"', svg, count=1)
        # WeasyPrint pode tratar width/height percentuais do SVG como dimensões
        # intrínsecas quando ele está dentro de um flex container. O estilo
        # inline força o mapa a ocupar o quadro sem ampliar o viewBox de forma
        # arbitrária, preservando o contorno e todas as informações desenhadas.
        svg = svg.replace('<svg ', '<svg style="display:block;width:100%;height:100%;" ', 1)
        return _crop_similarity_map_svg(svg, size, max(480, int(size * 0.625)))
    except Exception as e:
        logger.error(f"Erro ao gerar SimilarityMap cacheado: {str(e)}")
        return None


def generate_similarity_map(mol_ref: Optional[Mol], mol_test: Optional[Mol], size: int = SVG_SIZE_SIMILARITY_MAP) -> Optional[str]:
    try:
        if not _is_valid_molecule(mol_ref) or not _is_valid_molecule(mol_test):
            return None
        ref_key = _mol_cache_key(mol_ref)
        test_key = _mol_cache_key(mol_test)
        if not ref_key or not test_key:
            return None
        return _similarity_map_cached(ref_key, test_key, size)
    except Exception as e:
        logger.error(f"Erro ao gerar SimilarityMap: {str(e)}")
        return None


def mol_to_png(mol: Optional[Mol], size: int = PNG_SIZE_MOLECULE) -> Optional[bytes]:
    try:
        if not _is_valid_molecule(mol):
            return None
        smiles_key = _mol_cache_key(mol)
        if not smiles_key:
            return None
        encoded = _mol_to_png_cached(smiles_key, size)
        return None if not encoded else encoded
    except Exception as e:
        logger.error(f"Erro ao gerar PNG: {str(e)}")
        return None


@lru_cache(maxsize=128)
def _similarity_map_png_cached(ref_key: str, test_key: str, size: int) -> Optional[bytes]:
    try:
        mol_ref = Chem.MolFromSmiles(ref_key)
        mol_test = Chem.MolFromSmiles(test_key)
        if mol_ref is None or mol_test is None or not _is_valid_molecule(mol_ref) or not _is_valid_molecule(mol_test):
            return None

        mol_ref_no_h = Chem.RemoveHs(mol_ref)
        mol_test_no_h = Chem.RemoveHs(mol_test)
        if not _is_valid_molecule(mol_ref_no_h) or not _is_valid_molecule(mol_test_no_h):
            return None

        canvas_size = max(size, 800)
        canvas_height = max(480, int(canvas_size * 0.625))
        drawer = rdMolDraw2D.MolDraw2DCairo(canvas_size, canvas_height)
        opts = drawer.drawOptions()
        opts.bondLineWidth = 1.6
        opts.minFontSize = 16
        opts.annotationFontScale = 0.9
        opts.padding = 0.22
        opts.additionalAtomLabelPadding = 0.06

        try:
            SimilarityMaps.GetSimilarityMapForFingerprint(
                mol_ref_no_h,
                mol_test_no_h,
                _similarity_map_fingerprint,
                metric=DataStructs.TanimotoSimilarity,
                draw2d=drawer,
            )
        except Exception as inner_e:
            logger.warning(f"SimilarityMaps PNG falhou: {str(inner_e)}")
            return None

        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        logger.error(f"Erro ao gerar SimilarityMap PNG cacheado: {str(e)}")
        return None


def generate_similarity_map_png(mol_ref: Optional[Mol], mol_test: Optional[Mol], size: int = PNG_SIZE_SIMILARITY_MAP) -> Optional[bytes]:
    try:
        if not _is_valid_molecule(mol_ref) or not _is_valid_molecule(mol_test):
            return None
        ref_key = _mol_cache_key(mol_ref)
        test_key = _mol_cache_key(mol_test)
        if not ref_key or not test_key:
            return None
        return _similarity_map_png_cached(ref_key, test_key, size)
    except Exception as e:
        logger.error(f"Erro ao gerar SimilarityMap PNG: {str(e)}")
        return None


def get_3d_molblock(mol: Optional[Mol]) -> Optional[str]:
    try:
        if not _is_valid_molecule(mol):
            return None

        mol3d = Chem.AddHs(mol)
        if not _is_valid_molecule(mol3d):
            logger.warning("Molécula 3D inválida após AddHs().")
            return None

        embed_status = AllChem.EmbedMolecule(mol3d, randomSeed=42)

        if embed_status != 0:
            logger.warning("Falha ao embutir molécula em 3D")
            return None

        try:
            AllChem.MMFFOptimizeMolecule(mol3d, maxIters=500)
        except Exception as optimize_error:
            logger.warning(f"Falha na otimização MMFF: {str(optimize_error)}")

        return Chem.MolToMolBlock(mol3d)
    except Exception as e:
        logger.error(f"Erro ao gerar 3D: {str(e)}")
        return None


def calc_logd_vs_ph(mol: Optional[Mol], ph_range: Optional[List[int]] = None) -> Optional[List[Dict[str, Any]]]:
    if not mol:
        return None
    if ph_range is None:
        ph_range = list(range(PH_RANGE_START, PH_RANGE_END))
    try:
        logp = Crippen.MolLogP(mol)
        return [
            {"pH": ph, "LogD": round(logp + (NEUTRAL_PH - ph) * LOGD_PH_COEFFICIENT, 3)}
            if ph < NEUTRAL_PH else
            {"pH": ph, "LogD": round(logp - (ph - NEUTRAL_PH) * LOGD_PH_COEFFICIENT, 3)}
            for ph in ph_range
        ]
    except Exception as e:
        logger.error(f"Erro ao calcular LogD vs pH: {str(e)}")
        return None


def _safe_property_delta(reference: Any, candidate: Any) -> Any:
    if isinstance(reference, (int, float)) and isinstance(candidate, (int, float)):
        try:
            return round(abs(float(reference) - float(candidate)), 2)
        except (TypeError, ValueError):
            return "-"
    return "-"


def _append_warning(warnings: List[Dict[str, str]], code: str, message: str) -> None:
    if any(item.get("code") == code for item in warnings):
        return
    warnings.append({"code": code, "message": message})


def compare(
    smiles_ref: str,
    smiles_test: str,
    name_ref: str,
    name_test: str,
    fp_type: str,
    metric: str,
    show_map: bool = True
) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    mol_ref, err_ref = get_mol(smiles_ref)
    mol_test, err_test = get_mol(smiles_test)

    if err_ref or err_test:
        return None, err_ref or err_test
    if not mol_ref or not mol_test:
        return None, "SMILES inválido"

    props_ref = get_properties(mol_ref)
    props_test = get_properties(mol_test)
    if not props_ref or not props_test:
        return None, "Não foi possível calcular as propriedades das moléculas."

    fp_ref = get_fingerprint(mol_ref, fp_type)
    fp_test = get_fingerprint(mol_test, fp_type)

    if not fp_ref or not fp_test:
        return None, f"Erro ao gerar fingerprint {fp_type}"

    similarity = calc_similarity(fp_ref, fp_test, metric)
    classification = classify_similarity(similarity, metric)
    warnings: List[Dict[str, str]] = []

    if any("pKa" in key for key in props_ref.keys()):
        _append_warning(
            warnings,
            "pka_heuristic_estimate",
            "pKa é uma estimativa heurística aproximada e não uma previsão validada."
        )

    props_data = []
    for key in props_ref.keys():
        ref_val = props_ref.get(key)
        test_val = props_test.get(key)
        diff = _safe_property_delta(ref_val, test_val)
        props_data.append({"Propriedade": key, "Referência": ref_val, "Teste": test_val, "Diferença": diff})

    similarity_map = generate_similarity_map(mol_ref, mol_test) if show_map else None
    if show_map and similarity_map is None:
        _append_warning(
            warnings,
            "similarity_map_unavailable",
            "Mapa de similaridade indisponível; o score e as propriedades foram calculados normalmente."
        )

    png_ref = mol_to_png(mol_ref)
    png_test = mol_to_png(mol_test)
    similarity_map_png = generate_similarity_map_png(mol_ref, mol_test) if show_map else None
    molblock_ref = get_3d_molblock(mol_ref)
    molblock_test = get_3d_molblock(mol_test)
    if molblock_ref is None:
        _append_warning(
            warnings,
            "molblock_ref_unavailable",
            "Molblock 3D da referência indisponível; usando visualização 2D como fallback."
        )
    if molblock_test is None:
        _append_warning(
            warnings,
            "molblock_test_unavailable",
            "Molblock 3D do teste indisponível; usando visualização 2D como fallback."
        )

    logd_ref = calc_logd_vs_ph(mol_ref)
    logd_test = calc_logd_vs_ph(mol_test)

    result = {
        "name_ref": name_ref,
        "name_test": name_test,
        "smiles_ref": smiles_ref,
        "smiles_test": smiles_test,
        "similarity": similarity,
        "classification": classification,
        "properties": props_data,
        "svg_ref": mol_to_svg(mol_ref),
        "svg_test": mol_to_svg(mol_test),
        "png_ref": base64.b64encode(png_ref).decode("utf-8") if png_ref else None,
        "png_test": base64.b64encode(png_test).decode("utf-8") if png_test else None,
        "fingerprint_ref_png": fingerprint_to_png(mol_ref, fp_type),
        "fingerprint_test_png": fingerprint_to_png(mol_test, fp_type),
        "similarity_map": similarity_map,
        "similarity_map_png": base64.b64encode(similarity_map_png).decode("utf-8") if similarity_map_png else None,
        "molblock_ref": molblock_ref,
        "molblock_test": molblock_test,
        "fp_type": fp_type,
        "metric": metric,
        "similarity_map_fingerprint": SIMILARITY_MAP_FINGERPRINT,
        "similarity_map_metric": SIMILARITY_MAP_METRIC,
        "logd_ref": logd_ref,
        "logd_test": logd_test,
        "warnings": warnings,
    }
    return result, None


def _chemical_space_descriptor_vector(properties: Dict[str, Any]) -> List[Optional[float]]:
    vector = descriptor_vector(properties)
    return vector.tolist() if vector is not None else [None] * len(DESCRIPTOR_KEYS)


def build_chemical_space(
    ref_smiles: str,
    smiles_list: List[str],
    names_list: Optional[List[str]],
    results: List[Dict[str, Any]],
    fp_type: str = "Morgan2",
    metric: str = "Tanimoto",
    display_limit: int = 10,
) -> Dict[str, Any]:
    """Build a deterministic multimodal map and keep the nearest test molecules."""
    entries = [{"name": "Referência", "smiles": ref_smiles, "role": "reference"}]
    for index, item in enumerate(results):
        if item and not item.get("error") and item.get("smiles"):
            entries.append({
                "name": item.get("name") or f"Mol_{index + 1}",
                "smiles": item["smiles"],
                "role": "test",
                "similarity_to_reference": item.get("similarity"),
                "properties": item.get("properties") or {},
            })

    fingerprints = []
    property_vectors = []
    valid_entries = []
    for entry in entries:
        mol, error = get_mol(entry["smiles"])
        fp = get_fingerprint(mol, fp_type) if mol else None
        properties = entry.get("properties") or (get_properties(mol) if mol else None) or {}
        if error or mol is None or fp is None:
            continue
        fingerprints.append(fp)
        property_vectors.append(_chemical_space_descriptor_vector(properties))
        valid_entries.append({**entry, "properties": properties})

    count = len(valid_entries)
    limit = max(0, int(display_limit))
    if count == 0:
        return {
            "points": [],
            "method": "MDS clássico",
            "weights": {"structural": 0.6, "physicochemical": 0.4},
            "display_limit": limit,
            "descriptors": list(DESCRIPTOR_KEYS),
        }

    metrics = calculate_multimodal_space(fingerprints, property_vectors, metric)
    selected_indices = select_nearest_indices(
        metrics["reference_global_distances"],
        display_limit=limit,
        reference_index=0,
    )
    selected_distances = metrics["global_distances"][np.ix_(selected_indices, selected_indices)]
    coordinates = classical_mds(selected_distances)
    stress = normalized_stress(selected_distances, coordinates)

    points = []
    for local_index, original_index in enumerate(selected_indices):
        entry = valid_entries[original_index]
        point = {
            "name": entry["name"],
            "smiles": entry["smiles"],
            "role": entry["role"],
            "x": round(float(coordinates[local_index, 0]), 6),
            "y": round(float(coordinates[local_index, 1]), 6),
            "similarity_to_reference": round(float(1.0 - metrics["reference_structural_distances"][original_index]), 6),
            "structural_distance": round(float(metrics["reference_structural_distances"][original_index]), 6),
            "physicochemical_distance": round(float(metrics["normalized_physicochemical_distances"][0, original_index]), 6),
            "global_distance": round(float(metrics["reference_global_distances"][original_index]), 6),
            "rotatable_bonds": entry["properties"].get("Ligações rotacionáveis (RotB)"),
        }
        points.append(point)

    return {
        "points": points,
        "method": "MDS clássico sobre distância multimodal quadrática",
        "distance_formula": "sqrt(0,6 × (1 − similaridade MACCS/fingerprint)² + 0,4 × Dist.FQ_normalizada²)",
        "weights": {"structural": 0.6, "physicochemical": 0.4},
        "descriptors": ["Massa Molecular", "LogP", "TPSA", "HBD", "HBA", "RotB"],
        "normalization": "Z-score populacional no conjunto referência + lote; Dist.FQ = norma Euclidiana dos seis descritores; divisor = maior Dist.FQ em relação à referência",
        "reference_included": True,
        "display_limit": limit,
        "displayed_candidates": max(0, len(points) - 1),
        "total_valid_points": count,
        "fq_normalizer": round(float(metrics["fq_normalizer"]), 6),
        "mds_stress": round(float(stress), 6),
        "fingerprint": fp_type,
        "metric": metric,
    }

def bulk_compare(
    ref_smiles: str,
    smiles_list: List[str],
    names_list: Optional[List[str]] = None,
    fp_type: str = "Morgan2",
    metric: str = "Tanimoto"
) -> Tuple[Optional[List[Dict[str, Any]]], Optional[str]]:
    mol_ref, err = get_mol(ref_smiles)
    if err or not mol_ref:
        return None, err or "SMILES de referência inválido"

    fp_ref = get_fingerprint(mol_ref, fp_type)
    if not fp_ref:
        return None, f"Erro ao gerar fingerprint {fp_type}"

    fps = []
    results = []
    for i, smiles in enumerate(smiles_list):
        name = names_list[i] if names_list and i < len(names_list) else f"Mol_{i+1}"
        mol, err_mol = get_mol(smiles)
        if err_mol or not mol:
            results.append({"name": name, "similarity": None, "classification": None, "error": err_mol})
            continue
        fp = get_fingerprint(mol, fp_type)
        if not fp:
            results.append({"name": name, "similarity": None, "classification": None, "error": "Erro fingerprint"})
            continue
        fps.append(fp)
        properties = get_properties(mol) or {}
        mol_svg = mol_to_svg(mol, size=220)
        mol_png = mol_to_png(mol, size=220)
        results.append({
            "name": name,
            "smiles": smiles,
            "error": None,
            "svg": mol_svg,
            "png": base64.b64encode(mol_png).decode("utf-8") if mol_png else None,
            "properties": properties,
        })

    if not fps:
        return results, None

    similarities = DataStructs.BulkDiceSimilarity(fp_ref, fps) if metric == "Dice" else DataStructs.BulkTanimotoSimilarity(fp_ref, fps)

    sim_idx = 0
    for res in results:
        if res.get("error") is None and res.get("similarity") is None:
            sim = round(similarities[sim_idx], 4)
            res["similarity"] = sim
            res["classification"] = classify_similarity(sim, metric)
            sim_idx += 1

    return results, None