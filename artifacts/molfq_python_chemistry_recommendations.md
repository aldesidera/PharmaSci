# Recomendações críticas para o Mol.FQ

## Tese central

A geometria 3D simples ETKDGv3/UFF não deve permanecer como resultado principal do Mol.FQ. Ela só é útil como etapa intermediária para gerar descritores com significado físico ou eletrônico. Coordenadas e distâncias genéricas, sem relação demonstrada com o alvo do modelo, acrescentam complexidade e pouco valor ao usuário.

## Ferramentas prioritárias

AQME é a melhor referência para uma camada QM real: geração/seleção de conformeros, minimização, workflows com xTB/CREST e extração de descritores. Deve ser executado em modo batch/assíncrono, não durante cada consulta web.

Chemprop é a principal opção para um preditor calibrado de propriedades, mas só depois de existir dataset curado. O modelo deve ser comparado contra ESOL/RDKit e avaliado por divisão por scaffold, conjunto externo, incerteza e domínio de aplicabilidade.

ROBERT deve ser usado como pipeline de benchmarking e seleção de modelos, não como dependência de produção. Ele pode comparar descritores RDKit, Mordred, fingerprints e features QM, gerando métricas e relatórios reprodutíveis.

cclib deve normalizar a leitura de saídas de ORCA, Gaussian e outros motores QM, evitando que o Pharma.Sci implemente parsers próprios para cada programa.

MORFEUS só deve ser usado se forem escolhidos descritores com interpretação clara, como SASA, Sterimol, buried volume, dispersion descriptors ou xTB electronic descriptors. Não deve ser adicionado apenas para exibir coordenadas 3D.

datamol pode melhorar a preparação de dados, padronização, manipulação de moléculas, clustering, conformeros e curadoria, mas é uma camada de conveniência, não um modelo preditivo.

## Melhorias recomendadas

1. Padronizar entrada molecular: sanitização, dessalinização opcional, forma neutra/ionizada, tautomeria, estereoquímica, carga formal, InChIKey e hash da forma molecular usada.
2. Trabalhar com um ensemble de microespécies dependente de pH para solubilidade e LogD, em vez de uma única estrutura canônica.
3. Implementar um pipeline de features por níveis: RDKit/Mordred como baseline; descritores específicos de solubilidade; QM somente quando houver necessidade e infraestrutura.
4. Treinar Mol.FQ primeiro com ESOL/GSE/RDKit contra um modelo Chemprop ou ROBERT validado. Só aceitar o modelo avançado se superar o baseline em teste externo.
5. Usar AQME + xTB/CREST para calcular descriptors QM úteis: cargas, dipolo, energia HOMO/LUMO, gap, polarizabilidade, energia relativa de conformeros e descritores de solvatação quando disponíveis.
6. Usar cclib para converter os resultados QM em um schema interno comum.
7. Apresentar no relatório: valor baseline, valor calibrado, intervalo de incerteza, domínio de aplicabilidade, versão do modelo e proveniência. Não apresentar coordenadas cruas como resultado científico.

## O que descartar ou manter oculto

A tabela de coordenadas ETKDG/UFF, a distância média entre átomos e o rótulo genérico “geometria 3D enriquecida” devem ser removidos do fluxo padrão do usuário. Podem permanecer como artefatos internos de diagnóstico ou como entrada para um descritor específico.

Não se deve introduzir NBO/DLPNO-MP2, GNN ou DFT completo antes de existir um dataset que permita demonstrar ganho preditivo. Também não se deve usar HOMO/LUMO isoladamente como evidência de toxicidade, carcinogenicidade ou decisão regulatória.

## Ordem de implementação

Fase A: normalização de moléculas, microespécies, dataset e proveniência.

Fase B: baseline robusto de propriedades e validação por scaffold.

Fase C: Chemprop/ROBERT para o primeiro modelo calibrado de solubilidade.

Fase D: AQME + xTB/CREST em lote, com cclib e cache, gerando apenas descritores QM que provarem ganho.

Fase E: integração seletiva de QM ao modelo, com ablação baseline versus baseline+QM.

## Referências

[1] Awesome Python Chemistry. https://github.com/lmmentel/awesome-python-chemistry
[2] AQME. https://github.com/jvalegre/aqme
[3] Chemprop. https://github.com/chemprop/chemprop
[4] ROBERT. https://github.com/jvalegre/robert
[5] MORFEUS. https://github.com/digital-chemistry-laboratory/morfeus
[6] Sil S, Maskeri MA, Scheidt KA. graphpancake. https://doi.org/10.1186/s13321-026-01182-w
[7] graphpancake repository. https://github.com/sneha-sil/graphpancake
