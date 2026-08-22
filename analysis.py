"""
MolSim Final v3.3 — Análise Molecular Avançada
• Visualização 3D das moléculas
• Layout mais compacto e profissional
• Similarity Map aprimorado
• Suporte a PDF com PNGs
• Type hints, constantes e logging (v3.3)
"""

import logging
from typing import Tuple, Dict, Optional, List, Any
import base64

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors, Crippen, MACCSkeys
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps
from rdkit.Chem import Mol

# ============================================================================
# CONSTANTES
# ============================================================================

FINGERPRINT_BITS = 2048
MORGAN_RADIUS = 2
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


def get_mol(smiles: str) -> Tuple[Optional[Mol], Optional[str]]:
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None, "SMILES inválido – não foi possível parsear."

        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)
        logger.info(f"Molécula carregada com sucesso: {smiles[:50]}")
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


def get_properties(mol: Optional[Mol]) -> Optional[Dict[str, Any]]:
    if not mol:
        logger.warning("Tentativa de calcular propriedades com mol=None")
        return None

    try:
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        pka = estimate_pka(mol)

        return {
            "Massa Molecular (g/mol)": round(mw, 2),
            "Coeficiente de Partição (LogP)": round(logp, 3),
            "Área de Superfície Polar (Å²)": round(tpsa, 2),
            "pKa ácido": round(pka["pKa ácido"], 2) if pka and pka.get("pKa ácido") is not None else "N/A",
            "pKa básico": round(pka["pKa básico"], 2) if pka and pka.get("pKa básico") is not None else "N/A",
            "Doadores de H (HBD)": int(hbd),
            "Receptores de H (HBA)": int(hba),
            "Ligações Rotacionais": int(rot_bonds),
        }
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
        if fp_type == "Morgan2":
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                mol, radius=MORGAN_RADIUS, nBits=FINGERPRINT_BITS
            )
        elif fp_type == "RDKit":
            fp = Chem.RDKFingerprint(mol)
        elif fp_type == "MACCS":
            fp = MACCSkeys.GenMACCSKeys(mol)
        else:
            return None
        return fp
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

        drawer = rdMolDraw2D.MolDraw2DCairo(size, 120)
        drawer.DrawString(f"Fingerprint: {fp_type}")
        drawer.DrawString(f"Bits ativos: {bit_count}")
        drawer.DrawString(f"Tamanho: {size}px")
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode("utf-8")
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


def mol_to_svg(mol: Optional[Mol], size: int = SVG_SIZE_MOLECULE) -> str:
    try:
        if not mol:
            return ""
        mol_no_h = Chem.RemoveHs(mol)
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
        logger.error(f"Erro ao gerar SVG: {str(e)}")
        return ""


def _get_similarity_map_drawer(size: int):
    drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 3
    opts.minFontSize = 14
    return drawer


def _similarity_map_fingerprint(mol, idx):
    return SimilarityMaps.GetMorganFingerprint(
        mol,
        idx,
        radius=MORGAN_RADIUS,
        useFeatures=True,
        useChirality=True,
    )


def generate_similarity_map(mol_ref: Optional[Mol], mol_test: Optional[Mol], size: int = SVG_SIZE_SIMILARITY_MAP) -> Optional[str]:
    try:
        if not mol_ref or not mol_test:
            return None

        mol_ref_no_h = Chem.RemoveHs(mol_ref)
        mol_test_no_h = Chem.RemoveHs(mol_test)

        drawer = _get_similarity_map_drawer(size)

        try:
            SimilarityMaps.GetSimilarityMapForFingerprint(
                mol_ref_no_h,
                mol_test_no_h,
                _similarity_map_fingerprint,
                draw2d=drawer,
            )
        except Exception as inner_e:
            logger.warning(f"SimilarityMaps falhou: {str(inner_e)}")
            return None

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace(f'width="{size}px"', 'width="100%"').replace(f'height="{size}px"', 'height="100%"')
    except Exception as e:
        logger.error(f"Erro ao gerar SimilarityMap: {str(e)}")
        return None


def mol_to_png(mol: Optional[Mol], size: int = PNG_SIZE_MOLECULE) -> Optional[bytes]:
    try:
        if not mol:
            return None

        mol_no_h = Chem.RemoveHs(mol)

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
        return drawer.GetDrawingText()
    except Exception as e:
        logger.error(f"Erro ao gerar PNG: {str(e)}")
        return None


def generate_similarity_map_png(mol_ref: Optional[Mol], mol_test: Optional[Mol], size: int = PNG_SIZE_SIMILARITY_MAP) -> Optional[bytes]:
    try:
        if not mol_ref or not mol_test:
            return None

        mol_ref_no_h = Chem.RemoveHs(mol_ref)
        mol_test_no_h = Chem.RemoveHs(mol_test)

        canvas_size = max(size, 800)
        drawer = rdMolDraw2D.MolDraw2DCairo(canvas_size, canvas_size)
        opts = drawer.drawOptions()
        opts.bondLineWidth = 3
        opts.minFontSize = 16
        opts.annotationFontScale = 0.9

        try:
            SimilarityMaps.GetSimilarityMapForFingerprint(
                mol_ref_no_h,
                mol_test_no_h,
                _similarity_map_fingerprint,
                draw2d=drawer,
            )
        except Exception as inner_e:
            logger.warning(f"SimilarityMaps PNG falhou: {str(inner_e)}")
            return None

        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        logger.error(f"Erro ao gerar SimilarityMap PNG: {str(e)}")
        return None


def get_3d_molblock(mol: Optional[Mol]) -> Optional[str]:
    try:
        if not mol:
            return None

        mol3d = Chem.AddHs(mol)
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
    fp_ref = get_fingerprint(mol_ref, fp_type)
    fp_test = get_fingerprint(mol_test, fp_type)

    if not fp_ref or not fp_test:
        return None, f"Erro ao gerar fingerprint {fp_type}"

    similarity = calc_similarity(fp_ref, fp_test, metric)
    classification = classify_similarity(similarity, metric)

    props_data = []
    for key in props_ref.keys():
        ref_val = props_ref[key]
        test_val = props_test[key]
        diff = round(abs(ref_val - test_val), 2) if isinstance(ref_val, (int, float)) else "-"
        props_data.append({"Propriedade": key, "Referência": ref_val, "Teste": test_val, "Diferença": diff})

    similarity_map = generate_similarity_map(mol_ref, mol_test) if show_map else None
    png_ref = mol_to_png(mol_ref)
    png_test = mol_to_png(mol_test)
    similarity_map_png = generate_similarity_map_png(mol_ref, mol_test) if show_map else None
    molblock_ref = get_3d_molblock(mol_ref)
    molblock_test = get_3d_molblock(mol_test)
    logd_ref = calc_logd_vs_ph(mol_ref)
    logd_test = calc_logd_vs_ph(mol_test)

    return {
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
        "logd_ref": logd_ref,
        "logd_test": logd_test,
    }, None


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
        properties = get_properties(mol)
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