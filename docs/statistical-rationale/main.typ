#import "report-theme.typ": report-accent, report-theme

#show: report-theme.with(
  title: "Racional estatístico dos espaços químicos",
  author: "Manus AI",
  rhythm: "report",
  running-header: true,
)

#set par(first-line-indent: 0em)
#show link: set text(fill: report-accent)

#page(margin: (top: 30%, x: 2.2cm), numbering: none, header: none)[
  #align(center)[
    #text(size: 28pt, weight: "bold", fill: report-accent)[Racional estatístico]
    #v(0.25em)
    #text(size: 19pt, weight: "bold")[Dos espaços químicos]
    #v(0.8em)
    #text(size: 13pt, fill: luma(80))[Mol.Sim Batch · Mol.Sim Par a Par · Nitro.RA]
    #v(2em)
    #line(length: 42%, stroke: 0.7pt + report-accent)
    #v(2em)
    #text(size: 11pt)[Documento técnico de abordagem e interpretação]
    #v(0.4em)
    #text(size: 10pt, fill: luma(85))[PharmaSci · versão revisada · 04 de setembro de 2026]
  ]
]

#page(numbering: none, header: none)[
  #outline(title: [Sumário], indent: 1.5em)
]
#counter(page).update(1)

= Objetivo e escopo

Este documento descreve o racional estatístico atualmente adotado para as três formas de exploração de espaço químico disponíveis no PharmaSci: o mapa multimodal do *Mol.Sim Batch*, a comparação estrutural *Par a Par* e os mapas relativos do *Nitro.RA*. O objetivo é registrar, de forma auditável, quais dados entram em cada abordagem, como as distâncias são construídas, qual algoritmo de projeção é utilizado e quais interpretações não devem ser feitas.

A expressão “espaço químico” não representa uma única análise no aplicativo. O *Batch* e os espaços PubChem/EMA do *Nitro.RA* geram mapas 2D por MDS clássico sobre uma matriz de distância global. O *Par a Par* não gera um MDS: ele apresenta uma matriz/heatmap de similaridade estrutural entre duas moléculas, usando Morgan2 e Tanimoto. Essa distinção é central para evitar que uma métrica de um módulo seja atribuída a outro.

#block(fill: rgb("eef5ff"), stroke: (left: 3pt + report-accent), inset: 10pt, radius: 5pt)[
  *Regra de separação científica.* O Batch representa relações dentro do lote fornecido pelo usuário; o Par a Par representa uma comparação estrutural direta; o Nitro.RA cria dois espaços externos independentes, PubChem e EMA, com bibliotecas, perfis e interpretações próprios. Nenhum desses mapas transfere automaticamente atividade, toxicidade, AI ou equivalência regulatória.
]

== Protocolo de avaliação estatística

A avaliação foi conduzida como uma análise determinística de distância, e não como um modelo preditivo supervisionado. O procedimento foi definido antes da interpretação visual: primeiro são validados os SMILES e os descritores; depois são calculadas as dissimilaridades estruturais e físico-químicas; em seguida as componentes são combinadas; por fim, a matriz selecionada é projetada em duas dimensões. Essa ordem evita escolher vizinhos com base apenas na aparência do gráfico.

#enum(
  [*Validação e elegibilidade.* Cada molécula é convertida pelo RDKit, recebe SMILES canônico e só entra no cálculo se possuir fingerprint e os seis descritores finitos. O alvo é mantido como referência explícita.],
  [*Padronização.* Os descritores são transformados em z-scores para impedir que a massa molecular, o LogP ou qualquer outra variável domine a norma apenas por sua unidade ou amplitude.],
  [*Duas fontes de evidência.* A distância estrutural representa a informação do fingerprint; a distância físico-química representa a posição relativa nos seis descritores. As duas componentes permanecem observáveis separadamente.],
  [*Combinação pré-especificada.* A distância global usa média quadrática ponderada, com 0,60 para estrutura e 0,40 para descritores. O peso não é ajustado depois de observar o resultado.],
  [*Seleção e projeção.* Os candidatos são ordenados pela distância global em relação ao alvo; somente depois os até 10 mais próximos são projetados por MDS clássico. O stress é calculado sobre as distâncias do conjunto exibido.],
)

O racional estatístico é, portanto, de *triagem multimodal*: ele organiza relações químicas já observadas nos dados de entrada, mas não estima causalidade, potência carcinogênica, AI ou atividade farmacológica. A diferença entre o perfil local do Batch/PubChem e o perfil fixo da EMA é uma decisão de comparabilidade: o primeiro descreve a consulta atual; o segundo permite comparar consultas dentro de uma biblioteca regulatória versionada.

= Visão comparativa das três abordagens

A tabela resume as diferenças de entrada, métrica e saída. Ela é uma visão de arquitetura estatística, não uma tentativa de fundir os módulos em um único modelo.

#table(
  columns: (1.35fr, 1.55fr, 1.7fr, 1.55fr),
  inset: 7pt,
  align: (left, left, left, left),
  fill: (x, y) => if y == 0 { rgb("dbeafe") } else { none },
  [*Abordagem*], [*Biblioteca / conjunto*], [*Métrica de entrada*], [*Saída principal*],
  [*Mol.Sim Batch*], [Lote fornecido pelo usuário + referência], [MACCS/Tanimoto no mapa; distância multimodal quadrática], [MDS 2D, até 10 vizinhos e tabela de distâncias],
  [*Mol.Sim Par a Par*], [Duas moléculas fornecidas pelo usuário], [Morgan2/Tanimoto], [Heatmap de similaridade estrutural; sem MDS],
  [*Nitro.RA PubChem*], [Até 40 CIDs recuperados e filtrados], [MACCS/Tanimoto + seis descritores], [MDS relativo, até 10 candidatos e ranking],
  [*Nitro.RA EMA*], [Perfil fixo do Apêndice I da EMA], [MACCS/Tanimoto + perfil z-score versionado], [MDS relativo, até 10 referências e ranking],
)

== Vocabulário comum

A *similaridade estrutural* $S_(i,j)$ varia entre 0 e 1 quando calculada por Tanimoto em fingerprints binários. Para construir uma matriz de distâncias, o aplicativo usa a dissimilaridade estrutural:

$ d_"estrutural"(i,j) = 1 - S_(i,j) $

Assim, similaridade 1 corresponde a distância 0, enquanto similaridade 0 corresponde a distância 1. A transformação não muda a ordenação dos pares; apenas orienta a métrica para que valores maiores representem maior afastamento.

Os seis descritores físico-químicos compartilhados pela análise multimodal são massa molecular (MW), LogP, TPSA, HBD, HBA e ligações rotacionáveis (RotB). O Par a Par é uma exceção deliberada: seu heatmap é estrutural e não incorpora esses seis descritores.

= Mol.Sim Batch — espaço do lote do usuário

== Escopo e entrada

O Batch começa com uma molécula de referência e uma lista de moléculas fornecida pelo usuário. O espaço químico só é calculado quando a opção correspondente é ativada e é independente de qualquer biblioteca PubChem ou EMA. Cada entrada válida é sanitizada pelo RDKit, recebe o fingerprint e a métrica enviados pelo endpoint — na operação atual do aplicativo, MACCS/Tanimoto — e tem seus seis descritores calculados. A função interna permanece parametrizável para testes e extensões, mas o fluxo de produção do Batch fixa MACCS/Tanimoto para manter a comparação multimodal reproduzível.

A referência é sempre incluída como primeiro ponto. Entradas inválidas ou com erro são descartadas do espaço calculado, mas continuam identificáveis no resultado geral da comparação quando aplicável.

== Padronização dos descritores

Como os descritores possuem unidades e amplitudes distintas, o Batch ajusta um perfil populacional sobre a referência e o lote. A transformação é:

$ z_(i,k) = (x_(i,k) - μ_k) / σ_k $

em que $μ_k$ e $σ_k$ são a média e o desvio padrão populacional do descritor $k$. Valores ausentes são preenchidos pela mediana da coluna; quando a coluna é constante, o desvio é substituído por 1 para evitar divisão por zero. Nesse último caso, a coluna não cria separação entre as moléculas.

== Distância multimodal e seleção

A distância físico-química é a norma Euclidiana dos seis deltas padronizados. Ela é normalizada pelo maior valor observado em relação à referência:

$ d^* _"FQ"(i,j) = d_"FQ"(i,j) / max_l d_"FQ"(r,l) $

A distância global é uma combinação quadrática, não uma soma linear:

$ d_"global"(i,j) = sqrt(0.60 d_"estrutural"(i,j)^2 + 0.40 d^* _"FQ"(i,j)^2) $

A componente estrutural recebe peso 0,60 e a físico-química, 0,40. A seleção mantém a referência e até 10 moléculas com menor distância global em relação a ela. O ranking científico é, portanto, equivalente a ordenar por maior similaridade global.

== Projeção MDS e interpretação

O MDS clássico recebe a matriz completa de distâncias globais entre as estruturas selecionadas: referência–teste, teste–teste e diagonal nula. A matriz é duplamente centralizada:

$ B = -1/2 J D^2 J $

Os dois maiores autovalores não negativos e seus autovetores geram as coordenadas 2D. O MDS é uma projeção; seus eixos não são propriedades químicas e podem ser rotacionados, refletidos ou transladados sem modificar as distâncias.

O *stress MDS* quantifica a discrepância entre as distâncias originais e as distâncias Euclidianas no plano. Valores menores indicam melhor preservação relativa das relações; o indicador não valida toxicidade, equivalência farmacológica ou relevância regulatória.

#figure(
  image("chemical-space-example.png", width: 100%),
  caption: [Exemplo de projeção MDS sobre uma matriz de distância global. A proximidade visual deve ser conferida junto às métricas da tabela.],
) <fig:batch-mds>

= Mol.Sim Par a Par — comparação estrutural direta

== Escopo independente

O Par a Par compara duas moléculas fornecidas pelo usuário. Ele não utiliza a biblioteca PubChem, a biblioteca EMA, a seleção de vizinhos do Batch ou a distância multimodal de seis descritores. Seu objetivo é mostrar, de forma direta, quais regiões estruturais são responsáveis pela semelhança calculada.

== Morgan2 e Tanimoto

O fingerprint é Morgan2, isto é, um fingerprint circular de raio 2. A similaridade é calculada por Tanimoto entre os fingerprints binários. Para bits de presença, a interpretação é:

$ T(A,B) = c / (a + b - c) $

em que $a$ e $b$ são os números de bits ativos em cada molécula e $c$ é o número de bits ativos compartilhados. A escala varia de 0, nenhuma interseção observada, a 1, fingerprints idênticos.

== Heatmap e mapa de contribuição

O heatmap mostra a similaridade estrutural da comparação, enquanto o mapa molecular destaca regiões que contribuem para a comparação do fingerprint. A implementação utiliza Tanimoto explicitamente nas rotinas de similaridade do mapa. O enquadramento da estrutura é adaptativo para reduzir cortes do contorno em moléculas pequenas ou grandes.

Não há eixo MDS, stress, distância físico-química ou ranking de vizinhos no Par a Par. O valor estrutural pode ser acompanhado pelas propriedades físico-químicas exibidas no relatório, mas essas propriedades não alteram o heatmap.

= Nitro.RA — dois espaços químicos independentes

O Nitro.RA possui dois espaços químicos que não devem ser combinados: um espaço relativo baseado em candidatos PubChem e um espaço de referência baseado no perfil fixo do Apêndice I da EMA. Ambos usam MACCS/Tanimoto e a mesma combinação multimodal de seis descritores, mas diferem na origem da biblioteca e no tratamento estatístico do perfil.

#block(fill: rgb("fff7ed"), stroke: (left: 3pt + rgb("d97706")), inset: 10pt, radius: 5pt)[
  *Regra do alvo.* Quando o alvo já está presente em uma biblioteca, sua estrutura é removida dos candidatos por comparação de SMILES canônico. O alvo continua como referência do mapa, mas não é contado como vizinho nem incluído duas vezes no ranking ou nos cálculos comparativos.
]

= Nitro.RA PubChem — lote relativo recuperado

== Consulta e filtros

A busca consulta o endpoint de similaridade 2D do PubChem com até 40 CIDs. Os candidatos recuperados são convertidos pelo RDKit e filtrados pela presença do grupo N-nitroso usando o padrão SMARTS `[N;X3][N;X2]=O`. O próprio alvo é excluído quando seu SMILES canônico coincide com o de um candidato PubChem.

Após o filtro estrutural, são mantidas as moléculas com os seis descritores e fingerprints válidos. O candidato é ranqueado pela distância global; até 10 candidatos são exibidos. A tabela é apresentada por maior similaridade global de cima para baixo, enquanto a coluna de distância global cresce no sentido equivalente de menor para maior distância.

== Perfil estatístico relativo

O perfil z-score é ajustado sobre o conjunto usado naquela análise, formado pela referência e pelos candidatos PubChem válidos. A distância físico-química é normalizada em relação ao maior valor observado contra a referência nesse conjunto. Portanto, a escala é relativa à consulta e pode mudar quando o lote recuperado ou os candidatos elegíveis mudam.

A matriz selecionada é projetada por MDS clássico após a seleção dos 10 menores valores de distância global. O alvo é transladado para coordenadas `(0, 0)` somente para facilitar a leitura relativa do mapa; essa translação não altera nenhuma distância entre pontos.

== Limitações

PubChem é uma fonte externa sujeita a disponibilidade, cobertura e resposta da consulta. A ausência de uma estrutura no lote não prova sua ausência no universo químico. Similaridade 2D, descritores e MDS são ferramentas de triagem e não transferem automaticamente AI, cPCA ou conclusão regulatória.

= Nitro.RA EMA — perfil fixo do Apêndice I

== Construção da biblioteca

O espaço EMA usa um snapshot local do Apêndice I da EMA. O perfil contém apenas a planilha de nitrosaminas, deduplica as estruturas por SMILES canônico e calcula os seis descritores válidos. O alvo é removido da biblioteca quando já estiver listado, mas permanece como referência do mapa. O tamanho da biblioteca de origem, a versão, a data e a URL são mantidos no resultado para auditoria.

== Perfil z-score versionado

Ao contrário do PubChem, o perfil EMA é pré-ajustado e fixo. A média e o desvio padrão dos seis descritores são calculados na biblioteca de referência; o divisor físico-químico é fixado pelo maior valor de distância observado no perfil. Uma nova consulta não recalibra o perfil com os candidatos exibidos.

Essa escolha melhora a comparabilidade entre consultas dentro da mesma versão da biblioteca, mas não transforma o espaço EMA em uma representação completa do universo químico. Outliers não são removidos automaticamente e a atualização do snapshot muda o perfil de forma versionada.

== Parâmetros do perfil EMA atualmente instalado

O perfil versionado utilizado pelo aplicativo contém 243 estruturas únicas da planilha `N-nitrosamines`. A tabela registra os centros e escalas populacionais usados na transformação; ela permite reproduzir a padronização sem recalcular o perfil a partir dos candidatos de uma consulta.

#table(
  columns: (1.25fr, 1.15fr, 1.15fr),
  inset: 6pt,
  align: (left, right, right),
  fill: (x, y) => if y == 0 { rgb("dbeafe") } else { none },
  [*Descritor*], [*Centro μ*], [*Escala σ*],
  [MW], [348,519], [196,712],
  [LogP], [2,417], [1,941],
  [TPSA], [91,105], [68,751],
  [HBD], [1,309], [2,329],
  [HBA], [5,342], [3,520],
  [RotB], [6,103], [3,733],
)

O divisor global físico-químico fixado no snapshot é $18.009001$. A versão de referência é `EMA/42261/2025 Rev.13`, com atualização registrada em 24/06/2026. Esses valores são parâmetros de normalização, não limites toxicológicos e não devem ser interpretados como médias populacionais gerais fora da biblioteca EMA.

== Projeção e interpretação

A referência e os candidatos EMA válidos entram na distância multimodal, o alvo é transladado para `(0, 0)` no mapa e até 10 menores distâncias são exibidas. A translação é feita após o MDS e altera apenas a origem visual, não as relações entre pontos.
 PubChem e EMA aparecem em gráficos e tabelas separados. A proximidade com uma nitrosamina do Apêndice I não atribui automaticamente ao alvo o AI, a categoria cPCA ou qualquer decisão regulatória da referência.

= Comparação final e limitações

#table(
  columns: (1.5fr, 1.35fr, 1.35fr, 1.8fr),
  inset: 7pt,
  align: (left, left, left, left),
  fill: (x, y) => if y == 0 { rgb("dbeafe") } else { none },
  [*Dimensão*], [*Batch*], [*Par a Par*], [*Nitro.RA PubChem / EMA*],
  [Projeção], [MDS clássico], [Nenhuma; heatmap], [MDS clássico relativo],
  [Fingerprint], [MACCS/Tanimoto no mapa], [Morgan2/Tanimoto], [MACCS/Tanimoto],
  [Descritores], [MW, LogP, TPSA, HBD, HBA, RotB], [Não entram no heatmap], [Os mesmos seis descritores],
  [Biblioteca], [Lote do usuário], [Duas moléculas], [PubChem relativo ou EMA fixo],
  [Exibição], [Referência + até 10 vizinhos], [Duas estruturas], [Referência + até 10 candidatos],
)

A maior limitação comum é que uma projeção 2D comprime relações multidimensionais. O mapa deve ser lido junto com as distâncias quantitativas e com a qualidade dos dados de entrada. Nenhuma das abordagens substitui avaliação toxicológica, confirmação analítica, modelagem mecanística ou decisão regulatória.

A distância global é uma decisão de modelo com pesos explícitos, não uma verdade universal. Os valores 0,60 e 0,40 devem ser reavaliados caso o objetivo científico, o conjunto de referência ou a política de validação mude. O mesmo vale para a escolha de fingerprint, para o limite de 10 estruturas e para a versão do perfil EMA.

= Referências

[1] Gower, J. C. (1966). *Some Distance Properties of Latent Root and Vector Methods Used in Multivariate Analysis*. Biometrika, 53(3–4), 325–338. #link("https://doi.org/10.2307/2333639")[doi:10.2307/2333639]

[2] R Core Team. *Classical (Metric) Multidimensional Scaling*. R Documentation. #link("https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cmdscale.html")[R Documentation]

[3] Bajusz, D.; Rácz, A.; Héberger, K. (2015). *Why is Tanimoto index an appropriate choice for fingerprint-based similarity calculations?* Journal of Cheminformatics, 7, 20. #link("https://doi.org/10.1186/s13321-015-0069-3")[doi:10.1186/s13321-015-0069-3]

[4] RDKit. *The RDKit Book — Molecular Fingerprints and Similarity*. #link("https://www.rdkit.org/docs/RDKit_Book.html")[RDKit documentation]

[5] PubChem. *PUG REST — Programmatic Access*. #link("https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest")[PubChem PUG REST]

[6] EMA. *Questions and answers on nitrosamine impurities in human medicinal products — Appendix I*. #link("https://www.ema.europa.eu/en/human-regulatory-overview/post-authorisation/referral-assessment-reports-and-article-5-procedures/nitrosamine-impurities")[EMA nitrosamine guidance]

= Resumo operacional

#enum(
  [Identificar o módulo: Batch, Par a Par, Nitro.RA PubChem ou Nitro.RA EMA.],
  [No Batch e no Nitro.RA, calcular MW, LogP, TPSA, HBD, HBA e RotB; no Par a Par, manter o heatmap Morgan2/Tanimoto separado.],
  [Padronizar os descritores com o perfil apropriado: ajustado ao conjunto no Batch/PubChem ou fixo e versionado no EMA.],
  [Calcular a dissimilaridade estrutural por $1-S$ e a distância físico-química Euclidiana dos seis descritores quando aplicável.],
  [Combinar as componentes por distância quadrática com pesos 0,60 e 0,40.],
  [Selecionar até 10 menores distâncias globais e aplicar MDS clássico nos espaços Batch/PubChem/EMA.],
  [Interpretar o resultado como triagem quantitativa, sem transferir automaticamente atividade, toxicidade, AI ou decisão regulatória.],
)
