"""Renderizador SVG do fluxograma CPCA do Nitro.RA.

O desenho é deliberadamente determinístico: a lógica permanece no motor cPCA e
este módulo somente converte o decision_trace em uma visualização auditável.
"""
from __future__ import annotations

from html import escape
from typing import Any, Dict, Iterable, Mapping

W, H = 1200, 1832
CX = 420
DW, DH = 420, 122
RX, RW, RH = 920, 270, 76
LINE = "#8f625b"
ACCENT = "#f97316"
ACCENT_STRONG = "#ea580c"
ACCENT_TEXT = "#9a3412"
CORAL = "#ff6f61"
REVIEW = "#b45309"
NAVY = "#3d2020"


def _trace_map(result: Mapping[str, Any]) -> Dict[str, Mapping[str, Any]]:
    return {str(item.get("id")): item for item in (result.get("decision_trace") or []) if isinstance(item, Mapping)}


def _cls(base: str, active: bool = False, terminal: bool = False, review: bool = False) -> str:
    values = [base]
    if active:
        values.append("active")
    if terminal:
        values.append("terminal")
    if review:
        values.append("review")
    return " ".join(values)


def _diamond(cx: int, cy: int, lines: Iterable[str], node_id: str, active: bool = False) -> str:
    points = f"{cx},{cy- DH/2} {cx+DW/2},{cy} {cx},{cy+DH/2} {cx-DW/2},{cy}"
    lines = list(lines)
    start = cy - (len(lines)-1)*15
    text = "".join(f'<text x="{cx}" y="{start+i*30}" class="node-text">{escape(line)}</text>' for i, line in enumerate(lines))
    return f'<g id="{node_id}" class="{_cls("decision", active)}"><polygon points="{points}"/>{text}</g>'


def _result(y: int, category: str, ai: str, node_id: str, active: bool = False, review: bool = False) -> str:
    cls = _cls("result", active, terminal=active, review=review)
    return f'<g id="{node_id}" class="{cls}"><rect x="{RX-RW/2}" y="{y-RH/2}" width="{RW}" height="{RH}" rx="12"/><text x="{RX}" y="{y-4}" class="result-title">{escape(category)}</text><text x="{RX}" y="{y+23}" class="result-subtitle">{escape(ai)}</text></g>'


def _path(points: str, path_id: str, active: bool = False, review: bool = False) -> str:
    cls = _cls("connector", active, review=review)
    marker = "arrow-review" if review else ("arrow-active" if active else "arrow-neutral")
    return f'<path id="{path_id}" class="{cls}" d="{points}" marker-end="url(#{marker})"/>'


def render_cpca_flowchart(result: Mapping[str, Any]) -> str:
    trace = _trace_map(result)
    category = result.get("potency_category")
    score = result.get("potency_score")
    status = str(result.get("status") or "manual_review")

    def gate_active(gate_id: str) -> bool:
        return gate_id in trace

    def branch_active(gate_id: str, branch: str) -> bool:
        item = trace.get(gate_id) or {}
        state = item.get("status")
        if gate_id == "alpha_hydrogen_presence":
            return (branch == "yes" and state == "passed") or (branch == "no" and state == "triggered")
        if gate_id == "alpha_hydrogen_activation":
            return (branch == "yes" and state == "passed") or (branch == "no" and state == "triggered")
        if gate_id == "tertiary_alpha_carbon":
            return (branch == "yes" and state == "triggered") or (branch == "no" and state == "passed")
        return False

    score_path = lambda wanted: score is not None and ((wanted == "ge4" and score >= 4) or (wanted == "eq3" and score == 3) or (wanted == "eq2" and score == 2) or (wanted == "le1" and score <= 1))
    svg = [f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-labelledby="cpca-flow-title cpca-flow-desc">
<style>
svg {{ background:#f8fafc; font-family:DejaVu Sans,Arial,sans-serif; }}
.header {{ fill:{NAVY}; font-size:24px; font-weight:600; letter-spacing:.15px; }}
.subheader {{ fill:#8f625b; font-size:14px; }}
.start {{ fill:#fff1eb; stroke:{ACCENT_STRONG}; stroke-width:2.2; }}
.decision {{ fill:#fffaf8; stroke:{LINE}; stroke-width:2; }}
.decision.active {{ fill:#fff1eb; stroke:{ACCENT_STRONG}; stroke-width:3; }}
.node-text {{ fill:{NAVY}; font-size:14px; font-weight:500; letter-spacing:.05px; text-anchor:middle; dominant-baseline:middle; }}
.result {{ fill:#fff; stroke:#d8b8ad; stroke-width:1.6; }}
.result.active {{ fill:#fff1eb; stroke:{ACCENT_STRONG}; stroke-width:3; }}
.result.review {{ fill:#fff7ed; stroke:{REVIEW}; stroke-width:2; }}
.result-title {{ fill:{NAVY}; font-size:16px; font-weight:600; text-anchor:middle; }}
.result-subtitle {{ fill:#6f514b; font-size:13px; font-weight:500; text-anchor:middle; }}
.connector {{ fill:none; stroke:{LINE}; stroke-width:2.5; stroke-linecap:round; stroke-linejoin:round; }}
.connector.active {{ stroke:{ACCENT_STRONG}; stroke-width:3.5; }}
.connector.review {{ stroke:{REVIEW}; }}
.branch-label {{ fill:#6f514b; font-size:12px; font-weight:500; text-anchor:middle; }}
.branch-label.active {{ fill:{ACCENT_TEXT}; }}
.branch-label.review {{ fill:{REVIEW}; }}
.footer {{ fill:#8f625b; font-size:11px; text-anchor:middle; }}
</style>
<title id="cpca-flow-title">Fluxograma decisório CPCA</title>
<desc id="cpca-flow-desc">Fluxograma CPCA traduzido para português e destacado conforme o caminho calculado pelo motor.</desc>
<defs><marker id="arrow-neutral" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0,0 L7,3.5 L0,7 Q2,3.5 0,0" fill="{LINE}"/></marker><marker id="arrow-active" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0,0 L7,3.5 L0,7 Q2,3.5 0,0" fill="{ACCENT_STRONG}"/></marker><marker id="arrow-review" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0,0 L7,3.5 L0,7 Q2,3.5 0,0" fill="{REVIEW}"/></marker></defs>
<text x="600" y="34" class="header" text-anchor="middle">Fluxograma decisório CPCA</text>
<text x="600" y="56" class="subheader" text-anchor="middle">Caminho calculado pelo motor cPCA do Nitro.RA</text>
<rect x="260" y="84" width="320" height="48" rx="24" class="start"/><text x="420" y="114" class="node-text">Nitrosamina analisada</text>
''']
    ys = [205, 370, 535, 700, 865, 1030, 1195, 1360]
    qs = [
        (['A molécula possui', 'hidrogênios nos', 'carbonos alfa?'], 'alpha_hydrogen_presence'),
        (['Há mais de um hidrogênio', 'em um ou ambos os lados', 'do grupo N-nitroso?'], 'alpha_hydrogen_activation'),
        (['Existe carbono alfa', 'terciário?'], 'tertiary_alpha_carbon'),
        (['Calcular a pontuação de potência', 'conforme o Anexo A'], 'score_calc'),
        (['A pontuação de potência', 'é ≥ 4?'], 'score_ge4'),
        (['A pontuação de potência', 'é igual a 3?'], 'score_eq3'),
        (['A pontuação de potência', 'é igual a 2?'], 'score_eq2'),
        (['A pontuação de potência', 'é ≤ 1?'], 'score_le1'),
    ]
    for y, (lines, node_id) in zip(ys, qs):
        active = gate_active(node_id) or (node_id == 'score_calc' and score is not None) or (node_id == 'score_ge4' and score is not None) or (node_id == 'score_eq3' and score is not None) or (node_id == 'score_eq2' and score is not None) or (node_id == 'score_le1' and score is not None)
        svg.append(_diamond(CX, y, lines, node_id, active))

    # Vertical main flow.
    svg.append(_path('M420,132 L420,144 L420,144', 'start-to-q1'))
    for i in range(len(ys)-1):
        active = i < 3 and branch_active(['alpha_hydrogen_presence','alpha_hydrogen_activation','tertiary_alpha_carbon'][i], 'yes') or (i >= 3 and score is not None)
        svg.append(_path(f'M{CX},{ys[i]+DH/2} L{CX},{ys[i+1]-DH/2}', f'main-{i}', bool(active)))

    # Direct branches to results, all independent and aligned.
    branch_data = [
        (ys[0], 'cat5-alpha', 'Categoria 5', 'AI = 1500 ng/dia', branch_active('alpha_hydrogen_presence','no'), 'Não'),
        (ys[1], 'cat5-activation', 'Categoria 5', 'AI = 1500 ng/dia', branch_active('alpha_hydrogen_activation','no'), 'Não'),
        (ys[2], 'cat5-tertiary', 'Categoria 5', 'AI = 1500 ng/dia', branch_active('tertiary_alpha_carbon','yes'), 'Sim'),
        (ys[4], 'cat4', 'Categoria 4', 'AI = 1500 ng/dia', score_path('ge4'), 'Sim'),
        (ys[5], 'cat3', 'Categoria 3', 'AI = 400 ng/dia', score_path('eq3'), 'Sim'),
        (ys[6], 'cat2', 'Categoria 2', 'AI = 100 ng/dia', score_path('eq2'), 'Sim'),
        (ys[7], 'cat1', 'Categoria 1', 'AI = 18 ng/dia', score_path('le1'), 'Sim'),
    ]
    for y, node_id, cat, ai, active, label in branch_data:
        svg.append(_path(f'M{CX+DW/2},{y} L{RX-RW/2},{y}', node_id+'-path', active))
        svg.append(_result(y, cat, ai, node_id, active))
        svg.append(f'<text x="{(CX+DW/2+RX-RW/2)/2}" y="{y-10}" class="branch-label{" active" if active else ""}">{label}</text>')

    # Score <= 1: No -> manual review from the lower vertex, independent from Category 1.
    svg.append(_path(f'M{CX},{ys[7]+DH/2} L{CX},{1535} L{RX-RW/2},{1535}', 'manual-review-path', status == 'manual_review', review=True))
    svg.append(_result(1535, 'Revisão manual necessária', 'Fora do caminho padrão', 'manual-review', status == 'manual_review', review=True))
    svg.append(f'<text x="{CX+20}" y="{ys[7]+DH/2+24}" class="branch-label review">Não</text>')

    # Main-path labels
    labels = ['Sim','Sim','Não','Não','Não','Não','Não']
    for i, label in enumerate(labels):
        y = (ys[i]+ys[i+1])//2
        active = i < 2 and branch_active(['alpha_hydrogen_presence','alpha_hydrogen_activation'][i], 'yes') or i == 2 and branch_active('tertiary_alpha_carbon','no') or i >= 3 and score is not None
        svg.append(f'<text x="{CX+18}" y="{y}" class="branch-label{" active" if active else ""}">{label}</text>')
    svg.append('<line x1="90" y1="1735" x2="1110" y2="1735" stroke="#d5e0e8" stroke-width="1.5"/>')
    svg.append('<text x="600" y="1760" class="footer">Referência visual: EMA — Apêndice 2, Figura 2. Visualização baseada no decision_trace do motor CPCA.</text></svg>')
    return ''.join(svg)
