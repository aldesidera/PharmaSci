from pathlib import Path

W, H = 1900, 2900
bg = '#f8fafc'
navy = '#12304a'
blue = '#1d5f8a'
line = '#718494'
green = '#17845b'
green_fill = '#e8f6ef'
slate_fill = '#eef3f7'
orange = '#b45309'
orange_fill = '#fff4e5'
muted = '#526574'
font = 'DejaVu Sans, sans-serif'
parts = []

def esc(t): return t.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')
def text(x,y,s,size=26,fill=navy,weight='400',anchor='middle'):
    parts.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-family="{font}" font-size="{size}px" font-weight="{weight}" fill="{fill}">{esc(s)}</text>')
def rounded(x,y,w,h,fill,stroke=navy,r=16,sw=3):
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
def diamond(cx,cy,w,h,fill=slate_fill,stroke=line,sw=3):
    pts=f'{cx},{cy-h/2} {cx+w/2},{cy} {cx},{cy+h/2} {cx-w/2},{cy}'
    parts.append(f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
def connector(points,color=line,sw=3):
    d='M '+' L '.join(f'{x},{y}' for x,y in points)
    parts.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{sw}" stroke-linecap="round" stroke-linejoin="round" marker-end="url(#arrow)"/>')
def branch_label(x,y,s,color=muted): text(x,y,s,21,color,'700')

parts.append(f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">
<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L8,4 L0,8 Q2.5,4 0,0 z" fill="{line}"/></marker></defs>
<rect width="100%" height="100%" fill="{bg}"/>
''')
# Header
text(950,76,'Fluxograma decisório CPCA',42,navy,'700')
text(950,118,'Abordagem de categorização de potência carcinogênica para N-nitrosaminas',24,muted,'400')
text(950,158,'Prévia em português — caminho demonstrativo destacado em verde',20,green,'700')
parts.append('<line x1="140" y1="198" x2="1760" y2="198" stroke="#d5e0e8" stroke-width="2"/>')

cx=720; dw=620; dh=190; rx=1430; rw=430; rh=112
start_y=255
ys=[475,755,1035,1315,1595,1875,2155,2435]
questions=[
    ['A molécula possui','hidrogênios nos','carbonos alfa?'],
    ['Há mais de um hidrogênio','em um ou ambos os lados','do grupo N-nitroso?'],
    ['Existe carbono alfa','terciário?'],
    ['Calcular o Potency Score','conforme o Anexo A'],
    ['O Potency Score','é ≥ 4?'],
    ['O Potency Score','é igual a 3?'],
    ['O Potency Score','é igual a 2?'],
    ['O Potency Score','é ≤ 1?'],
]
# start node
rounded(cx-205,start_y,410,70,'#e7f0f8',blue,35,3); text(cx,start_y+45,'Nitrosamina analisada',27,navy,'700')
# decision nodes
for i,(y,lines) in enumerate(zip(ys,questions)):
    is_score_calc=i==3
    diamond(cx,y,dw,dh,green_fill if is_score_calc else slate_fill,green if is_score_calc else line,4 if is_score_calc else 3)
    base=y-(len(lines)-1)*17
    for j,s in enumerate(lines): text(cx,base+j*34,s,23,navy,'700')
# result cards
results={
    475: ('Categoria 5','AI = 1500 ng/dia',False),
    755: ('Categoria 5','AI = 1500 ng/dia',False),
    1035: ('Categoria 5','AI = 1500 ng/dia',False),
    1595: ('Categoria 4','AI = 1500 ng/dia',False),
    1875: ('Categoria 3','AI = 400 ng/dia',True),
    2155: ('Categoria 2','AI = 100 ng/dia',False),
    2435: ('Categoria 1','AI = 18 ng/dia',False),
}
for y,(cat,ai,selected) in results.items():
    rounded(rx-rw/2,y-rh/2,rw,rh,green_fill if selected else '#ffffff',green if selected else '#9aabb8',14,4 if selected else 2)
    text(rx,y-7,cat,27,green if selected else navy,'700')
    text(rx,y+30,ai,23,green if selected else muted,'600')
# manual review
rounded(rx-rw/2,2635,rw,96,orange_fill,orange,14,3)
text(rx,2674,'Revisão manual necessária',23,orange,'700')
text(rx,2704,'fora do caminho padrão',19,orange,'400')

# Exact vertical flow, vertex-to-vertex
connector([(cx,start_y+70),(cx,ys[0]-dh/2)],line,3)
for a,b in zip(ys[:-1],ys[1:]):
    color=green if a in (475,755,1035,1315,1595) else line
    connector([(cx,a+dh/2),(cx,b-dh/2)],color,3)

# Side branches: every outcome is independent and goes right from the correct node.
# Q1: Não -> Cat5; Q2: Não -> Cat5; Q3: Sim -> Cat5
branch_specs=[
    (ys[0],'Não',results[ys[0]][2],0),
    (ys[1],'Não',False,1),
    (ys[2],'Sim',False,2),
    (ys[4],'Sim',False,4),
    (ys[5],'Sim',True,5),
    (ys[6],'Sim',False,6),
    (ys[7],'Sim',False,7),
]
for y,label,selected,_ in branch_specs:
    color=green if selected else line
    connector([(cx+dw/2,y),(rx-rw/2,y)],color,3)
    branch_label((cx+dw/2+rx-rw/2)/2,y-15,label,green if selected else muted)
# Q8 Não -> review: saída inferior independente, sem competir com o ramo Sim
connector([(cx,ys[7]+dh/2),(cx,2683),(rx-rw/2,2683)],orange,2.5)
branch_label(cx+28,ys[7]+dh/2+42,'Não',orange)
# Main-path labels, centered between nodes and never over a shape
main_labels=[(ys[0]+ys[1])//2,'Sim',(ys[1]+ys[2])//2,'Sim',(ys[2]+ys[3])//2,'Não',(ys[3]+ys[4])//2,'Não',(ys[4]+ys[5])//2,'Não',(ys[5]+ys[6])//2,'Não',(ys[6]+ys[7])//2,'Não']
for i in range(0,len(main_labels),2):
    y=main_labels[i]; s=main_labels[i+1]
    branch_label(cx+32,y,s,green if y in (615,895,1175) else muted)
# Footer
parts.append('<line x1="140" y1="2790" x2="1760" y2="2790" stroke="#d5e0e8" stroke-width="2"/>')
text(950,2822,'Referência visual: EMA — Appendix 2, Figura 2. A decisão deve ser calculada pelo motor CPCA.',17,muted,'400')
parts.append('</svg>')
Path('/home/ubuntu/PharmaSci_git/artifacts/cpca_fluxograma_profissional.svg').write_text(''.join(parts),encoding='utf-8')
