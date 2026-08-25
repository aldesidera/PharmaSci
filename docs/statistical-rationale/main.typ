#import "report-theme.typ": report-accent, report-theme

#show: report-theme.with(
  title: "Racional estatístico do espaço químico 2D",
  author: "Manus AI",
  rhythm: "report",
  running-header: true,
)

#page(margin: (top: 30%, x: 2.2cm), numbering: none, header: none)[
  #set par(first-line-indent: 0em)
  #align(center)[
    #text(size: 26pt, weight: "bold", fill: report-accent)[Racional estatístico do espaço químico 2D]
    #v(0.5em)
    #text(size: 14pt, fill: luma(80))[Mol.Sim — comparação molecular em lote]
    #v(2em)
    #line(length: 40%, stroke: 0.5pt + luma(160))
    #v(2em)
    #text(size: 12pt)[Manus AI  •  25 de agosto de 2026]
  ]
]

#page(numbering: none, header: none)[
  #outline(title: [Sumário], indent: 1.5em)
]
#counter(page).update(1)

= Objetivo e escopo

Este documento registra o racional estatístico utilizado pelo módulo *Mol.Sim* para construir um espaço químico bidimensional a partir de uma comparação molecular em lote. O objetivo não é produzir uma classificação farmacológica automática, mas oferecer uma representação auditável da proximidade estrutural e físico-química entre a molécula de referência e as moléculas testadas.

A análise combina duas fontes de informação. A primeira é a similaridade estrutural calculada sobre fingerprints moleculares. A segunda é a distância entre descritores físico-químicos padronizados. A combinação atualmente adotada no aplicativo é de 60% para a componente estrutural e 40% para a componente físico-química.

= Padronização dos descritores por Z-score

Os descritores apresentam escalas e unidades diferentes. Massa molecular, LogP, TPSA, HBD, HBA e RotB não podem ser combinados diretamente sem padronização, pois uma variável com maior amplitude numérica dominaria a distância. Por isso, cada descritor é convertido em Z-score usando a média e o desvio padrão populacional calculados sobre o conjunto formado pela referência e pelo lote.

A transformação aplicada a um valor $x_i$ é:

$ z_i = (x_i - μ) / σ $

em que $μ$ é a média do descritor no conjunto analisado e $σ$ é o desvio padrão populacional. O uso do desvio populacional corresponde a considerar o lote fornecido para aquela análise como o universo de comparação, e não como uma amostra destinada a estimar uma população externa.

Quando um descritor não está disponível, o valor ausente é substituído pela mediana observada para aquele descritor no conjunto. Quando os valores são constantes, o desvio padrão é substituído por 1 apenas para evitar divisão por zero; nessa situação, a coluna não gera diferença entre as moléculas.

#table(
  columns: (2.3fr, 1fr, 2.7fr),
  inset: 7pt,
  fill: luma(242),
  [*Descritor*], [*Unidade/forma*], [*Uso no modelo*],
  [Massa molecular], [g/mol], [Z-score],
  [LogP], [adimensional], [Z-score],
  [TPSA], [Å²], [Z-score],
  [HBD e HBA], [contagem], [Z-score],
  [RotB], [contagem], [Z-score; inclui ligações rotacionáveis],
)

= Distância físico-química

Para cada par de moléculas $i$ e $j$, calcula-se a diferença entre os Z-scores correspondentes. Como a planilha utiliza variações absolutas, a implementação usa diferenças absolutas antes da agregação. A distância físico-química é a norma Euclidiana desses seis deltas:

$ d_"FQ"(i,j) = sqrt(sum_(k=1)^6 (z_(i, k) - z_(j, k))^2) $

Os seis descritores são MM, LogP, TPSA, HBD, HBA e RotB. Portanto, $d_"FQ"$ é sempre maior ou igual a zero. Em seguida, a distância é normalizada pelo maior valor da distância físico-química em relação à referência no lote:

$ d^* _"FQ"(i,j) = d_"FQ"(i,j) / max_l d_"FQ"(r,l), quad 0 ≤ d^* _"FQ" ≤ 1 $

A normalização coloca a componente físico-química em uma escala comparável à distância estrutural, que também varia entre zero e um para fingerprints binários.

= Similaridade, dissimilaridade e distância global

O fingerprint produz uma *similaridade*: quanto maior o valor, mais características moleculares são compartilhadas. Para construir uma matriz de distâncias, o aplicativo converte essa medida em dissimilaridade por complemento:

$ d_"estrutural"(i,j) = 1 - S_(i,j) $

Assim, similaridade 1 corresponde a distância 0, enquanto similaridade 0 corresponde a distância 1. Essa transformação não muda a ordenação dos pares; apenas inverte a orientação para que valores maiores signifiquem maior afastamento.

A distância global utilizada pelo espaço químico é:

$ d_"global"(i,j) = 0.60 d_"estrutural"(i,j) + 0.40 d^* _"FQ"(i,j) $

A escolha da dissimilaridade, em vez da similaridade, é necessária porque o MDS recebe uma matriz de dissimilaridades/distâncias e procura pontos cujas distâncias geométricas reproduzam essa matriz. Se a similaridade fosse fornecida diretamente, moléculas muito parecidas seriam tratadas como geometricamente afastadas, invertendo o significado do mapa. A similaridade original continua preservada nos resultados, tooltips e tabelas; somente o insumo da projeção é convertido para distância.

#table(
  columns: (2.2fr, 1.2fr, 3fr),
  inset: 7pt,
  fill: luma(242),
  [*Métrica*], [*Faixa típica*], [*Interpretação*],
  [Similaridade], [0–1], [1 = mais semelhante],
  [Dissimilaridade estrutural], [0–1], [0 = mais semelhante; 1 = mais distante],
  [Dist.FQ normalizada], [0–1], [0 = mesmas propriedades no modelo],
  [Distância global], [0–1], [0 = maior proximidade composta],
)

= Como funciona o MDS

*MDS* significa *Multidimensional Scaling*, ou escala multidimensional. É uma técnica de ordenação que recebe uma matriz de dissimilaridades e procura representar seus objetos em um espaço com menos dimensões. Neste caso, a matriz é convertida em duas coordenadas para formar o gráfico X/Y.

A entrada do MDS não é a tabela de similaridades isoladas em relação à referência. É a matriz completa de distâncias globais entre os pares do conjunto: referência–teste, teste–teste e referência–referência. Essa característica permite que o mapa represente também a organização interna do lote.

No MDS clássico, a matriz de distâncias $D$ é elevada ao quadrado e duplamente centralizada para produzir uma matriz de produtos escalares:

$ B = -1/2 J D^2 J $

A decomposição espectral de $B$ fornece autovalores e autovetores. Os dois maiores autovalores não negativos são usados para construir as coordenadas:

$ X_(2D) = V_(2D) sqrt(Λ_(2D)) $

O resultado é um conjunto de pontos cuja distância Euclidiana tenta reproduzir as dissimilaridades originais. Como a representação é bidimensional, pode haver alguma distorção quando muitas relações não cabem exatamente em um plano. O gráfico deve, portanto, ser interpretado junto com as métricas quantitativas.

Os eixos são coordenadas matemáticas, não propriedades químicas. Podem assumir valores positivos ou negativos, e a orientação pode ser refletida ou rotacionada sem alterar as distâncias. O sinal de X ou Y não significa similaridade positiva/negativa, nem atividade biológica maior/menor.

#figure(
  image("chemical-space-example.png", width: 100%),
  caption: [Exemplo do espaço químico 2D. O gráfico é uma projeção MDS sobre a matriz de distância global; os valores `G` correspondem à distância global em relação à referência.],
) <fig:chemical-space>

#quote(block: true)[
  *Como ler o gráfico:* pontos próximos representam moléculas próximas segundo a combinação definida de estrutura e propriedades. A referência é destacada em laranja; as moléculas comparadoras aparecem em azul. As coordenadas MDS podem ser negativas, mas as métricas de similaridade e distância mostradas nos rótulos/tabelas permanecem nas suas escalas próprias.
]

= Racional estatístico e limitações

O procedimento de Z-score reduz o efeito das unidades e das diferenças de escala entre os descritores, mas a média e o desvio padrão dependem do conjunto analisado. Por isso, a posição de uma molécula pode mudar quando o lote é alterado. A inclusão de RotB é deliberada: o número de ligações rotacionáveis captura uma dimensão de flexibilidade conformacional simples, complementar à informação de fingerprint e aos descritores de polaridade, tamanho e lipofilicidade.

A combinação 60:40 é uma decisão explícita do modelo. Ela dá maior influência à estrutura, mas permite que diferenças físico-químicas alterem a distância final. A escolha não deve ser interpretada como uma validação universal desses pesos; idealmente, eles devem ser avaliados em um conjunto de referência adequado ao objetivo regulatório ou toxicológico específico.

O MDS é uma ferramenta exploratória e de triagem. Proximidade no plano não constitui prova de equivalência farmacológica, toxicológica ou regulatória. A interpretação deve considerar a estrutura molecular, a qualidade dos fingerprints, as propriedades calculadas, as incertezas dos descritores e as evidências experimentais disponíveis.

= Referências

[1] Gower, J. C. (1966). *Some Distance Properties of Latent Root and Vector Methods Used in Multivariate Analysis*. Biometrika, 53(3–4), 325–338. #link("https://doi.org/10.2307/2333639")[doi:10.2307/2333639]

[2] R Core Team. *Classical (Metric) Multidimensional Scaling*. R Documentation, função `cmdscale`. A documentação descreve MDS clássico como uma representação de dissimilaridades por pontos cujas distâncias são aproximadamente iguais às dissimilaridades de entrada e destaca a indeterminação por translação, rotação e reflexão. #link("https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cmdscale.html")[R Documentation]

[3] Bajusz, D.; Rácz, A.; Héberger, K. (2015). *Why is Tanimoto index an appropriate choice for fingerprint-based similarity calculations?* Journal of Cheminformatics, 7, 20. #link("https://doi.org/10.1186/s13321-015-0069-3")[doi:10.1186/s13321-015-0069-3]

[4] Andrade, C. (2021). *Z-Scores, Standard Scores, and Composite Test Scores*. Indian Journal of Psychological Medicine, 43(5), 451–453. #link("https://pmc.ncbi.nlm.nih.gov/articles/PMC8826187/")[PMC8826187]

= Apêndice: resumo operacional

#enum(
  [Calcular os descritores MM, LogP, TPSA, HBD, HBA e RotB para a referência e o lote.],
  [Substituir ausências pela mediana do descritor e calcular Z-scores com média e desvio padrão populacional.],
  [Calcular a similaridade de fingerprint e convertê-la em dissimilaridade por $1-S$.],
  [Calcular a Dist.FQ Euclidiana e normalizá-la em relação ao maior valor observado contra a referência.],
  [Combinar as componentes com pesos 0,60 e 0,40.],
  [Aplicar MDS clássico à matriz completa de distâncias globais e representar os dois primeiros eixos.],
)
