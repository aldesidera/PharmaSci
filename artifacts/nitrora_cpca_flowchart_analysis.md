# Análise do fluxograma CPCA do EMA para o Nitro.RA

## Fonte oficial

EMA, Appendix 2 to Questions and Answers for marketing authorisation holders/applicants on nitrosamine impurities. Figura 2, página 4: `Flowchart to Predict the Potency Category of an N-nitrosamine`.

Fonte: https://www.ema.europa.eu/en/documents/other/appendix-2-carcinogenic-potency-categorisation-approach-n-nitrosamines_en.pdf

## Estrutura decisória observada

1. A N-nitrosamina possui hidrogênios em seus carbonos alfa? Não -> Categoria de potência 5, AI 1500 ng/dia.
2. Possui mais de um hidrogênio em um ou ambos os lados do grupo N-nitroso? Não -> Categoria 5, AI 1500 ng/dia.
3. Possui carbono alfa terciário? Sim -> Categoria 5, AI 1500 ng/dia.
4. Calcular o Potency Score. Score >= 4? Sim -> Categoria 4, AI 1500 ng/dia.
5. Score = 3? Sim -> Categoria 3, AI 400 ng/dia.
6. Score = 2? Sim -> Categoria 2, AI 100 ng/dia.
7. Score <= 1? Sim -> Categoria 1, AI 18 ng/dia.

O fluxograma possui nós de decisão, ramos Sim/Não, caixas de resultado e notas de rodapé para a definição do carbono alfa terciário e referência ao Anexo A para cálculo do Potency Score.

## Recomendação

Reproduzir como diagrama vetorial declarativo em Mermaid para o painel web e gerar uma versão SVG/PNG para o PDF. A lógica decisória deve continuar no motor Python do cPCA; o diagrama deve ser apenas uma visualização do `decision_trace` retornado pelo motor. Isso evita divergência entre o gráfico e o cálculo.

D2 é uma alternativa para diagrama mais editorial, mas Mermaid é mais simples de integrar ao frontend existente e ao renderer local. O Nitro.RA deve apresentar a versão em português, preservar a referência oficial do EMA e destacar o caminho efetivamente percorrido pelo alvo.

Não usar uma imagem estática como única fonte da lógica. O diagrama deve receber estado: nós avaliados, respostas Sim/Não, score, categoria, AI, exclusões e necessidade de revisão manual.

## Limitações importantes

A implementação deve reproduzir a lógica do documento EMA, mas não deve ser apresentada como decisão regulatória automática. Casos fora do escopo, múltiplos grupos N-nitroso, estruturas não suportadas ou divergência entre regra e fonte devem resultar em `manual_review`.
