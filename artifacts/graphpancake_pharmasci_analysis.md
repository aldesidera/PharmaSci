# Análise do artigo graphpancake para o Pharma.Sci

## Evidências principais

O artigo apresenta o graphpancake como uma biblioteca Python que transforma resultados de cálculos de estrutura eletrônica em grafos moleculares hierárquicos. As fontes incluem coordenadas 3D, termodinâmica, NPA/JANPA e NBO, com atributos em nós, ligações e nível molecular. O trabalho define quatro níveis: DFT, NPA, NBO e QM.

O pipeline descrito combina geração de conformeros com GFN2-xTB/ORCA GOAT, otimização/frequências em r2SCAN-3c, termodinâmica com Shermo, energia de ponto único DLPNO-MP2, NBO 7.0 e NPA/JANPA. Isso é uma referência arquitetural, não uma dependência que deva ser executada no fluxo web interativo do Pharma.Sci.

A implementação usa CLI, processamento em lote paralelo, YAML para configuração de recursos, cache/lazy loading, controle de memória e SQLite para armazenar grafos e permitir consultas e exportações. O artigo relata que os recursos quânticos melhoraram R² em aproximadamente 0,3–0,5 em relação a recursos derivados de SMILES em tarefas de regressão, mas também informa que os ganhos não foram consistentes em classificação e que DFT isolado não foi suficiente para superar recursos de SMILES em MPNN.

Os autores observam que falhas nos cálculos removeram moléculas dos datasets, introduzindo possível viés para moléculas mais rígidas. Também destacam custo computacional, dependência de dados experimentais e limitações de tamanho dos datasets.

## Mapeamento para o Pharma.Sci

O motor compartilhado `chemo_suite/core/chemical_space.py` já oferece um ponto de extensão claro: hoje usa seis descritores (`MW`, `LogP`, `TPSA`, `HBD`, `HBA`, `RotB`), perfil z-score, distância estrutural por fingerprint, distância físico-química e distância global ponderada. A extensão recomendada é adicionar um perfil opcional de descritores eletrônicos/3D sem alterar o modo clássico por padrão.

`analysis.py` ainda depende majoritariamente de descritores RDKit 2D, fingerprints e heurísticas de pKa/solubilidade. O módulo `mol_fq/solubility.py` já possui contrato para artefato calibrado versionado, manifest, métricas e fallback explícito para ESOL/GSE. Esse contrato pode ser reutilizado para modelos de propriedades quânticas e de solubilidade.

`core/conformer.py` atualmente gera conformero com ETKDGv3 e UFF, mas não persiste uma geometria/ensemble rico nem extrai descritores eletrônicos. `apps/nitro_ra/quantum.py` ainda é placeholder, portanto o artigo sugere uma direção para esse módulo, mas não deve ser integrado diretamente sem uma camada assíncrona, cache e validação.

## Recomendação de produto

Prioridade 1: criar uma camada de proveniência e qualidade por molécula, com versão do RDKit, método de conformação, força/campo ou nível QM, status, custo, falha, timestamp, hash do SMILES normalizado e fonte dos dados.

Prioridade 2: separar claramente os modos `rápido`, `enriquecido` e `QM`. O modo rápido mantém o app atual; o enriquecido usa conformero ETKDG/UFF, geometria 3D e cargas/descritores semiempíricos quando disponíveis; o modo QM completo deve ser batch/assíncrono e nunca bloquear a geração interativa do relatório.

Prioridade 3: introduzir uma representação molecular comum com níveis de complexidade. O nível básico usa RDKit; o intermediário acrescenta coordenadas, distâncias, ângulos, propriedades atômicas e geometria; o avançado acrescenta cargas NPA, ordens de ligação e orbitais quando houver software e licença disponíveis.

Prioridade 4: expandir o espaço químico com uma camada opcional de distância eletrônica. A distância global atual deve permanecer como baseline comparável. A nova distância deve ser reportada separadamente e somente combinada após validação em dataset, evitando alterar retroativamente os resultados publicados.

Prioridade 5: usar os descritores QM como candidatos para modelos calibrados de Mol.FQ, metabolismo e Nitro.RA, sempre comparando baseline RDKit/SMILES contra modelo híbrido e quantificando ganho, incerteza, domínio de aplicabilidade e custo.

Prioridade 6: adotar o padrão de reprodutibilidade do artigo: CLI, configuração YAML, processamento em lote, SQLite/Parquet, cache, lazy loading, logs estruturados e testes de integridade das geometrias.

## Limitações e cautelas

O artigo não justifica colocar DFT/NBO/DLPNO-MP2 em todas as consultas web. O pipeline é caro, depende de ORCA, JANPA, Shermo e, em parte, NBO 7.0; falhas podem enviesar o conjunto. Além disso, o ganho reportado é tarefa- e dataset-dependente. Para Nitro.RA, HOMO/LUMO isolados não devem ser tratados como evidência regulatória ou mecanismo toxicológico suficiente.

## Referências

[1] Sil S, Maskeri MA, Scheidt KA. graphpancake: A Python package for representing organic molecules as molecular graphs utilizing electronic structure theory. Journal of Cheminformatics, 2026. https://doi.org/10.1186/s13321-026-01182-w

[2] graphpancake repository. https://github.com/sneha-sil/graphpancake

[3] graphpancake documentation. https://graphpancake.readthedocs.io/
