# Recomendações para Mol.SM e Nitro.RA

## Conclusão central

Mol.SM e Nitro.RA não devem receber a mesma camada tecnológica. Mol.SM é um módulo de comparação estrutural e físico-química; Nitro.RA é um sistema de avaliação toxicológica, regulatória e mecanística. A geometria 3D genérica ETKDGv3/UFF não deve ser exibida em nenhum dos dois como resultado principal.

## Mol.SM

### Melhorias de maior valor

1. Padronização molecular antes da comparação: sanitização, sais, tautômeros, estereoquímica, carga, InChIKey e hash da forma analisada.
2. Fingerprints explicitamente versionados: MACCS/Tanimoto para comparação principal e Morgan2/Tanimoto somente quando a análise exigir maior sensibilidade estrutural.
3. Distinção entre similaridade estrutural, distância físico-química e distância global. O gráfico deve mostrar os três componentes separadamente quando solicitado.
4. Diagnóstico de validade: estruturas inválidas, propriedades ausentes, duplicatas, outliers e moléculas fora do domínio do lote.
5. Benchmark estatístico: bootstrap ou intervalo de confiança da similaridade, análise de sensibilidade dos pesos 0,6/0,4 e comparação com pesos alternativos.
6. Curadoria de lote com datamol/RDKit e clustering/diversidade para evitar que o resultado seja dominado por duplicatas.
7. Heatmap explicável: manter MACCS/Tanimoto como métrica oficial do heatmap e mostrar bits/subestruturas responsáveis pela similaridade sem adicionar QM.

### Ferramentas indicadas

Datamol para preparação, normalização, clustering e manipulação de lotes; RDKit para fingerprints e descritores; scikit-learn para bootstrap, MDS, métricas e validação. ROBERT e Chemprop não são prioridades no Mol.SM, pois o módulo não é originalmente um preditor.

### O que não adicionar

Não adicionar DFT, HOMO/LUMO, xTB ou graphpancake ao fluxo padrão do Mol.SM. Essas features podem introduzir complexidade e alterar a noção de similaridade sem melhorar a tarefa principal. Um módulo separado de análise eletrônica só seria justificável se houver uma pergunta química explícita.

## Nitro.RA

### Melhorias de maior valor

1. Padronização molecular e forma química: registrar sais, tautômeros, estereoquímica, carga e forma efetivamente usada no CPCA e na busca EMA/PubChem.
2. Workflow CPCA visual e auditável: mostrar grupos estruturais identificados, regra aplicada, categoria, justificativa e referência EMA, sem transformar a classificação em uma caixa preta.
3. Biblioteca EMA versionada: conservar nome, CAS, SMILES, AI, categoria CPCA, versão da fonte, data de atualização e motivo de exclusão do alvo.
4. Espaço químico com diagnóstico: separar similaridade MACCS/Tanimoto, distância físico-química, distância global, MDS e stress; informar se o alvo foi excluído da biblioteca.
5. Incerteza e limites: deixar claro que proximidade química não transfere automaticamente Acceptable Intake, CPCA ou toxicidade.
6. Metabolism em camadas: separar regra estrutural, evidência de acessibilidade/reatividade, Deep-PK externo e conclusão integrada.
7. QM seletivo e mecanístico: utilizar QM somente para priorizar reatividade relativa ou sítios candidatos, nunca para substituir CPCA, EMA ou evidência toxicológica.
8. Proveniência de fontes externas: registrar endpoint, timestamp, payload resumido, timeout, cache, status e versão do conector.

### Ferramentas indicadas

Datamol/RDKit para padronização e busca; ROBERT para eventual modelo de priorização se houver dataset toxicológico curado; AQME + xTB/CREST para análises QM em lote; cclib para leitura de ORCA/Gaussian; MORFEUS somente para descritores 3D/QM com hipótese mecanística clara.

### Graphpancake no Nitro.RA

O conceito mais útil é a representação em três níveis: atributos de átomos, ligações e molécula inteira, associada a proveniência. Não é necessário copiar a infraestrutura completa. Um grafo eletrônico somente deve ser criado quando houver resultados QM armazenados e uma hipótese sobre reatividade, por exemplo cargas, dipolo, gap ou propriedades locais no entorno da função N-nitroso.

### O que não adicionar

Não usar HOMO/LUMO isoladamente como marcador de carcinogenicidade, potência ou AI. Não substituir CPCA/EMA por distância global, PubChem, Deep-PK ou QM. Não executar DFT, NBO ou xTB durante uma requisição interativa sem fila, cache e timeout.

## Arquitetura comum recomendada

| Camada | Mol.SM | Nitro.RA |
|---|---|---|
| Identidade | SMILES canônico, InChIKey, estereoquímica | Igual, mais forma regulatória e sal |
| Estrutura 2D | MACCS, Morgan2, subestruturas | SMARTS N-nitroso, grupos CPCA |
| Descritores | MW, LogP, TPSA, HBD, HBA, RotB | Mesmo conjunto, com interpretação de triagem |
| 3D/QM | Não usar no padrão | Opcional, somente para hipótese mecanística |
| Estatística | Similaridade, MDS, bootstrap, clustering | Espaço PubChem/EMA com diagnóstico e exclusão do alvo |
| Modelo ML | Não prioritário | Somente após dataset toxicológico curado |
| Relatório | Comparabilidade e incerteza estatística | Evidência, fonte, regra, status externo e limitação |

## Ordem recomendada

### Mol.SM

Primeiro padronizar moléculas e duplicatas; depois validar a distância global e seus pesos; em seguida adicionar intervalos de confiança, clustering/diversidade e explicação do heatmap. Somente depois avaliar novos fingerprints.

### Nitro.RA

Primeiro fortalecer o workflow CPCA/EMA e a proveniência; depois revisar o espaço PubChem/EMA com diagnóstico; em seguida estabilizar estados Deep-PK e Metabolism; somente então testar QM seletivo para reatividade.

## Referências

[1] Awesome Python Chemistry. https://github.com/lmmentel/awesome-python-chemistry
[2] AQME. https://github.com/jvalegre/aqme
[3] Chemprop. https://github.com/chemprop/chemprop
[4] ROBERT. https://github.com/jvalegre/robert
[5] MORFEUS. https://github.com/digital-chemistry-laboratory/morfeus
[6] Sil S, Maskeri MA, Scheidt KA. graphpancake. https://doi.org/10.1186/s13321-026-01182-w
