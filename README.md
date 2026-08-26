# 🧪 MolSim v10 — build estabilizado

Versão ativa de trabalho do MolSim/MolSim_ver10, com foco em estabilização de API, validação de payloads, segurança operacional, consistência do fluxo pair/batch/report/export e evolução incremental do Nitro.RA.

Esta documentação reflete o estado atual do projeto em desenvolvimento: funcional, mais resiliente e com contratos de erro padronizados, sem reescrever a lógica científica principal.

---

## ✅ Estado atual

- Backend Flask com validação centralizada de JSON
- Contrato de erro padronizado: `{"error":{"code":"invalid_request","message":"...","field":"..."}}`
- Rotas críticas protegidas contra payload malformado
- `/healthz` disponível
- CORS restrito por variável de ambiente
- limites de conteúdo configuráveis por ambiente
- export PDF validado com base64/PNG/dimensões/tamanho
- suporte temporário para `show_logd` e compatibilidade com `show_similarity_map`
- cliente mais consistente em pair/batch/preview/export

---

## 📋 Funcionalidades principais

- ✅ Comparação individual entre duas moléculas
- ✅ Comparação em lote (batch)
- ✅ Validação de fingerprint e métrica
- ✅ Preview de relatório visual
- ✅ Export para PDF via preview
- ✅ Fallback de nomes via PubChem quando o campo está vazio
- ✅ Health check /healthz
- ✅ Interface responsiva e melhor consistência do estado do cliente
- ✅ Nitro.RA: seleção de cPCA, Quantum, Metabolism e busca opcional de nitrosaminas para o mesmo SMILES
- ✅ Nitro.RA: quatro abas de resultados independentes; cPCA e Metabolism funcionais, Quantum preparado para evolução
- ✅ Nitro.RA cPCA: estrutura química, limite FDA, motivos do resultado e lookup do AI EMA Appendix 1
- ✅ Nitro.RA espaço químico: busca PubChem, filtro N-nitroso por RDKit, ranking MACCS/Tanimoto, descritores MW/LogP/TPSA/HBD/HBA/RotB, distância físico-química e mapa PCA
- ✅ Nitro.RA Metabolism: predição estrutural CYP450 por α-hidroxilação, sítios alfa vulneráveis, metabólitos de Fase I e intermediários diazônio hipotéticos

---

## 🚀 Instalação e execução

> Entrypoint oficial único: `main.py`  
> Não execute `app.py` diretamente.

### Windows (PowerShell)

```powershell
cd C:\caminho\para\MolSim_ver10
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python main.py
```

Abra no navegador:

```text
http://127.0.0.1:5000
```

### Linux/macOS

```bash
cd /caminho/para/MolSim_ver10
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python main.py
```

---

## 🔧 Variáveis de ambiente úteis

```bash
MOLSIM_CORS_ORIGINS=http://127.0.0.1:3000,http://localhost:3000
MOLSIM_MAX_CONTENT_LENGTH=2097152
MOLSIM_MAX_BATCH_ITEMS=100
MOLSIM_MAX_SMILES_LENGTH=4096
MOLSIM_MAX_NAME_LENGTH=256
MOLSIM_MAX_BATCH_REQUEST_BYTES=524288
MOLSIM_MAX_EXPORT_IMAGE_BYTES=2097152
MOLSIM_MAX_EXPORT_IMAGE_WIDTH=4096
MOLSIM_MAX_EXPORT_IMAGE_HEIGHT=4096

# Busca opcional de espaço químico no PubChem
MOLSIM_PUBCHEM_SPACE_THRESHOLD=50
MOLSIM_PUBCHEM_SPACE_MAX_RECORDS=100
MOLSIM_PUBCHEM_SPACE_MAX_CANDIDATES=10
MOLSIM_PUBCHEM_SPACE_TIMEOUT=12
MOLSIM_PUBCHEM_SPACE_CACHE_TTL_SECONDS=3600
```

Observações:
- Se `MOLSIM_CORS_ORIGINS` não for informado, o app restringe o acesso aos domínios locais padrão (`127.0.0.1` e `localhost`).
- O app não abre em `0.0.0.0` por padrão em execução local.

---

## 📁 Estrutura relevante

```text
MolSim_ver10/
├── main.py
├── app.py
├── analysis.py
├── chemo_suite/
│   ├── main.py
│   ├── core/
│   │   ├── parser.py
│   │   └── conformer.py
│   └── apps/
│       ├── mol_sim/
│       │   ├── pairwise.py
│       │   └── batch.py
│       └── nitro_ra/
│           ├── cpca.py
│           ├── quantum.py
│           ├── metabolism.py
│           ├── deep_pk.py
│           └── nitrosamine_space.py
├── requirements.txt
├── static/
├── templates/
│   ├── index.html
│   └── report_preview.html
├── tests/
│   ├── test_phase1_app.py
│   ├── test_modular_structure.py
│   ├── test_cpca.py
│   ├── test_metabolism.py
│   ├── test_deep_pk.py
│   └── test_main_cli.py
├── chemo_suite/apps/nitro_ra/data/ema_appendix1.json

├── README.md
└── ...
```

---

## 🧪 Como usar o app

1. Insira os SMILES da referência e de comparação.
2. Opcionalmente preencha nomes para as moléculas.
3. Escolha fingerprint e métrica.
4. Clique em comparar.
5. Visualize o score, classificação, mapa e propriedades.
6. Use o preview de relatório para revisar o conteúdo e exportar para PDF.
7. No modo batch, informe a referência e a lista de moléculas.
8. Para o Nitro.RA, ative o módulo, informe um SMILES, marque um ou mais checkboxes e execute as análises.
9. Navegue pelas abas cPCA, Espaço Químico, Quantum e Metabolism para consultar os resultados separadamente; módulos futuros ficam explicitamente identificados como em desenvolvimento.
10. Marque `Espaço Químico` quando quiser consultar, sob demanda, compostos estruturalmente semelhantes no PubChem; a busca interna mantém o filtro específico para candidatos N-nitroso.
11. Informe a dose diária máxima em mg/dia quando desejar a conversão do AI para ppm.
12. Marque `Deep-PK` junto com `Metabolism` para consultar, depois da análise local, somente os endpoints externos de substrato e inibição CYP.
13. Leia a tabela Deep-PK como uma previsão probabilística complementar; ela não gera metabólitos nem intermediários e não substitui confirmação experimental.

---

## Fase 6 — Operação e manutenção

Objetivos desta etapa:
- reduzir risco operacional em rotas críticas e no ciclo de erro do cliente/servidor;
- manter o app seguro em produção local com limites e headers consistentes;
- preservar compatibilidade com a API já estabilizada;
- seguir sem reescrever a lógica científica principal.

Itens reforçados nesta etapa:
- respostas internas de erro sem vazar stack trace ao cliente;
- tratamento explícito de rotas JSON e requisições problemáticas;
- health check com suporte a `HEAD` para diagnósticos e monitoramento leve;
- documentação de critérios operacionais para continuidade do roadmap.

---

## Nitro.RA — cPCA

O primeiro módulo científico funcional do Nitro.RA implementa uma triagem estrutural baseada na [orientação RAIL da FDA](https://www.fda.gov/media/170794/download), usando o fluxo da Figura 1 e as pontuações do Appendix A. O motor registra a versão da regra, a fonte regulatória, os centros N-nitroso detectados, os carbonos alfa, os hidrogênios alfa, as features estruturais e a justificativa do resultado.

A implementação é deliberadamente conservadora. Estruturas sem centro N-nitroso retornam `not_nitrosamine`; estruturas inválidas retornam `invalid_smiles`; estruturas fora do escopo cPCA ou com padrão não mapeado retornam `manual_review` ou `not_applicable`. A [EMA mantém orientação própria e documentos atualizados para nitrosaminas](https://www.ema.europa.eu/en/human-regulatory-overview/post-authorisation/pharmacovigilance-post-authorisation/referral-procedures-human-medicines/nitrosamine-impurities/nitrosamine-impurities-guidance-marketing-authorisation-holders), portanto o resultado não deve ser tratado como limite universal entre jurisdições.

O cPCA é uma estimativa de triagem para apoio técnico e não substitui avaliação toxicológica, dados composto-específicos, read-across ou decisão regulatória final.

---

## Nitro.RA — Metabolism CYP450

O módulo **Metabolism** implementa uma primeira camada local e determinística de predição estrutural de Fase I. Ele reconhece padrões N-nitroso alifáticos e aromáticos e, quando encontra carbonos alifáticos sp3 com pelo menos um hidrogênio diretamente adjacentes ao nitrogênio N-nitroso, gera hipóteses de **α-hidroxilação CYP450**. Uma N-nitrosamina aromática pode ser reconhecida sem possuir carbono α sp3 elegível; nesse caso o estado é `no_alpha_sites`, sem geração artificial de produtos. Essa escolha é coerente com a literatura que trata a α-hidroxilação como uma etapa crítica da ativação metabólica de N-nitrosaminas, mas a regra não substitui um modelo treinado, dados enzimáticos ou confirmação experimental [1].

Para cada sítio elegível, o endpoint retorna o índice do átomo, hidrogênios α, indicação de anel, contexto enzimático de **CYP2E1/CYP3A4**, regra aplicada, SMILES do produto α-hidroxilado, SVG estrutural e propriedades físico-químicas quando disponíveis. O resultado também inclui uma representação de fragmento diazônio como **intermediário mecanístico hipotético**; essa estrutura é uma hipótese de triagem e não deve ser interpretada como espécie isolada, produto experimental confirmado ou decisão de risco.

A implementação atual usa `prediction_mode: rule_based` e o identificador versionado `CYP450_ALPHA_HYDROXYLATION_N_NITROSO`. O campo `reaction_smarts` mantém a transformação auditável; a aplicação do OH é localizada pelo índice do carbono α para não confundir múltiplas correspondências em substratos assimétricos. O código expõe `predict_cyp450_metabolism(smiles)` e mantém `evaluate_metabolism(smiles)` para compatibilidade com a rota existente. BioTransformer 3.0 e SyGMa são referências externas importantes para comparação e evolução futura: BioTransformer oferece uma opção CYP450 baseada em regras e/ou aprendizado de máquina [2] [3], enquanto SyGMa usa regras de Fase I/Fase II e pontuação empírica para ordenar produtos potenciais [4]. Eles não são chamados automaticamente pelo fluxo local atual, evitando dependência de serviço externo e mantendo o resultado reproduzível.

Estados relevantes do módulo são `ok`, `invalid_smiles`, `not_nitrosamine` e `no_alpha_sites`. Em todos os casos, o painel informa os limites da inferência e preserva a distinção entre uma hipótese computacional e uma conclusão toxicológica.

## Nitro.RA — Deep-PK complementar

O complemento **Deep-PK** é opcional e só é consultado quando o usuário marca `Deep-PK` junto com `Metabolism`. A interface primeiro mostra a estrutura molecular, os sítios α, as hipóteses de α-hidroxilação e os intermediários locais; abaixo desse conteúdo, exibe uma tabela externa com as previsões de **Substrato?**, **Probabilidade**, **Inibidor?** e **Probabilidade** para as isoformas CYP documentadas pelo serviço. A confiança textual retornada pelo Deep-PK aparece sob cada classificação.

O Deep-PK não substitui o motor local e não é usado para gerar os SMILES dos metabólitos ou do diazônio-surrogate. A resposta é marcada como proveniente de serviço externo, pode permanecer em processamento por meio de `job_id` e pode retornar `deep_pk_unavailable`, `deep_pk_error` ou `deep_pk_timeout`. Falhas externas não invalidam o resultado local de Metabolism. O fluxo envia o SMILES canônico ao serviço público somente após ação explícita do usuário.

### Endpoints Deep-PK

#### `/nitro-ra/deep-pk`

- Recebe `POST` em JSON com `smiles`.
- Sanitiza e canonicaliza o SMILES localmente antes do envio.
- Submete `pred_type=metabolism` em `multipart/form-data` ao Deep-PK.
- Retorna `status: running` e um `job_id`, ou um estado explícito de indisponibilidade/erro.

#### `/nitro-ra/deep-pk/<job_id>`

- Recebe `GET` para consultar o `job_id` assíncrono.
- Envia o identificador ao Deep-PK no formato multipart documentado.
- Normaliza os endpoints de substrato e inibição para CYP1A2, CYP2C19, CYP2C9, CYP2D6 e CYP3A4.
- Não retorna produtos metabólicos nem interpreta a classificação como prova de biotransformação.

A URL padrão é `https://biosig.lab.uq.edu.au/deeppk/api/predict` e pode ser sobrescrita com `DEEP_PK_API_URL`. O timeout individual pode ser ajustado por `MOLSIM_DEEP_PK_TIMEOUT`; a consulta continua opcional e não é executada para o usuário que não marcar o checkbox. A documentação oficial do serviço deve ser consultada antes de uso operacional [5].

### Referências científicas

[1] [Chakravarti et al. — Computational Prediction of Metabolic α-Carbon Hydroxylation Potential of N-Nitrosamines](https://pmc.ncbi.nlm.nih.gov/articles/PMC10283024/).
[2] [BioTransformer — Overview and CYP450 Help](https://biotransformer.ca/help).
[3] [Wishart et al. — BioTransformer 3.0, Nucleic Acids Research](https://pmc.ncbi.nlm.nih.gov/articles/PMC9252798/).
[4] [Ridder & Wagener — SyGMa, PubMed PMID 18311745](https://pubmed.ncbi.nlm.nih.gov/18311745/).
[5] [Deep-PK — API Documentation](https://biosig.lab.uq.edu.au/deeppk/api_docs).

---

## Fase 7 — Critérios finais e validação

Critérios finais de pronto para a etapa de encerramento:
- API com payloads malformados rejeitados sem 500;
- contratos de erro estáveis em rotas críticas;
- `healthz` funcional e sem dependência de PubChem;
- CORS e host default em modo local seguro;
- export PDF e preview protegidos contra entrada inválida;
- validação automatizada de regressões em `tests/` e sem quebra de fluxos estabilizados.

Atenção:
- itens de refinamento de UX visual e RDKit mais aprofundado continuam como evolução futura fora do escopo da estabilização atual;
- a base atual permanece priorizando compatibilidade e segurança operacional.

---

## 📊 Fluxos de API principais

### /compare
- Requer JSON válido em `application/json`
- Valida `smiles_ref`, `smiles_test`, `fp_type`, `metric`
- Rejeita payloads malformados com 400/415

### /bulk-compare
- Valida `ref_smiles`, `smiles_list`, `names_list`, métricas e limites
- Aplica limites configuráveis para batch e payload

### /nitro-ra/cpca
- Requer JSON válido em `application/json`
- Recebe `smiles` e aceita `mdd_mg` opcional
- Retorna evidências estruturais, contagem de centros N-nitroso, Potency Score, categoria e AI em ng/dia
- Retorna `manual_review` ou `not_applicable` para estruturas fora do escopo suportado, sem inventar uma categoria
- Usa a conversão `ppm = AI (ng/dia) / dose diária máxima (mg)` quando `mdd_mg` é informado

### /nitro-ra/analyze
- Requer JSON válido em `application/json`
- Recebe `smiles`, uma lista `modules` com `cpca`, `quantum`, `metabolism` e/ou `nitrosamine_space`, e `mdd_mg` opcional
- Retorna um objeto `results` separado por módulo, preservando o resultado de cada análise para as abas da interface
- Quantum continua retornando `status: not_implemented`, sem produzir valores fictícios
- Metabolism retorna hipóteses locais `rule_based` de α-hidroxilação CYP450; não representa confirmação experimental nem substitui BioTransformer, SyGMa ou avaliação toxicológica
- Quando cPCA é selecionado, o resultado inclui `structure_svg`, `canonical_smiles` e o objeto `ema`
- O objeto `ema` consulta o snapshot `EMA/42261/2025 Rev.13`, atualizado em 24/06/2026, por SMILES canônico
- Uma estrutura ausente do Appendix 1 retorna `ema.status: not_listed`; o sistema não infere AI EMA a partir da categoria FDA
- `nitrosamine_space` consulta, sob demanda, a similaridade 2D do [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest), filtra `[N;X3][N;X2]=O`, seleciona até 10 candidatos e retorna descritores, Tanimoto, distância físico-química e pontos PCA
- O PubChem fornece o lote de CIDs; o ranking de Tanimoto é recalculado localmente com MACCS, usando o pipeline de similaridade do Mol.Sim, e a seleção final também considera a distância físico-química min–max calculada com MW, LogP, TPSA, HBD, HBA e RotB
- Falhas de rede/timeout retornam `status: pubchem_unavailable`; nenhum resultado retorna `status: no_nitrosamines`; esses estados não bloqueiam o cPCA nem inventam candidatos
- Metabolism pode retornar `not_nitrosamine`, `no_alpha_sites` ou `invalid_smiles`; N-nitrosos aromáticos reconhecidos sem carbono α sp3 retornam `no_alpha_sites`; quando `ok`, informa `alpha_sites`, `metabolites`, `reactive_intermediates`, `rule_id`, contexto enzimático e avisos de triagem
- O complemento Deep-PK não é incluído em `modules` do endpoint local; quando selecionado na UI, o frontend chama `/nitro-ra/deep-pk` depois da resposta local e consulta `/nitro-ra/deep-pk/<job_id>` até concluir ou atingir o limite de espera

### /nitro-ra/deep-pk
- Requer JSON válido em `application/json` e recebe `smiles`
- Canonicaliza o SMILES localmente e submete `pred_type=metabolism` ao Deep-PK
- Retorna `running` com `job_id`, ou estados `deep_pk_unavailable`, `deep_pk_error` e `invalid_smiles`

### /nitro-ra/deep-pk/<job_id>
- Consulta um job assíncrono do Deep-PK via `GET`
- Retorna cinco isoformas com objetos separados de `substrate` e `inhibitor`, cada um com `prediction`, `probability` e `interpretation`
- Não altera nem mistura as hipóteses estruturais locais de Metabolism

### /report-preview
- Gera o preview do relatório visual
- Não deve falhar em 500 para payload inválido

### /export-pdf
- Valida `similarity`, tipos, valores finitos, base64, PNG e limites
- Rejeita erros de entrada antes da geração do PDF

### /healthz
- Retorna `200` e `{"status":"ok"}` sem depender de PubChem

---

## 🐛 Resolução de problemas comuns

### Content-Type inválido
- Use `Content-Type: application/json` nas requisições.

### JSON inválido
- Verifique a sintaxe do corpo JSON.

### 500 em payload malformado
- Isso foi tratado na versão estabilizada; use o contrato padronizado de erro.

### PubChem indisponível
- A busca de espaço químico é opcional e só é executada quando o checkbox correspondente está marcado. HTTP 503/504, timeout ou falha de conexão aparecem na aba como `PubChem indisponível`, sem lançar erro no frontend; o cPCA e os demais módulos selecionados continuam preservados.
- A busca usa cache temporário e limita o lote inicial. A ausência de um composto no lote retornado não prova ausência no universo químico; o resultado é uma triagem estrutural e não substitui avaliação toxicológica, confirmação analítica, read-across ou decisão regulatória.
- O lookup de nome usado em outros fluxos é tolerante a falha e não bloqueia o fluxo principal.

### Preview/PDF não abre
- Verifique popup blockers e use o botão de export diretamente no app.

---

## 🔗 Recursos úteis

- RDKit: https://www.rdkit.org/
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
- PubChem PUG REST: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- SMILES: https://en.wikipedia.org/wiki/Simplified_molecular_input_line_entry_system

---

## 📝 Status e histórico

**Versão ativa:** v10
**Status:** estabilizado e validado em fluxo principal
**Escopo atual:** hardening de API, cliente e preview/export, sem reescrita da lógica científica principal

### Evolução relevante na v10
- validação centralizada de JSON
- `/healthz`
- CORS restringido
- `MAX_CONTENT_LENGTH` por ambiente
- export PDF validado
- tratamento consistente de erros
- UX/cliente mais resiliente

---

**MolSim v10 — estabilizado para uso local e validação contínua.**


## Espaço químico 2D no Mol.Sim

No modo **Batch**, o checkbox **Construir espaço químico 2D** habilita a geração de um gráfico com a referência e, no máximo, os **10 vizinhos mais próximos** do lote. O posicionamento usa MDS clássico sobre a distância multimodal quadrática `sqrt(0,6 × D_struct² + 0,4 × D_FQ²)`, em que `D_struct = 1 − Tanimoto` e `D_FQ` é a distância físico-química normalizada. O espaço químico usa MACCS/Tanimoto, independentemente do fingerprint escolhido para a comparação estrutural principal.

A componente físico-química utiliza massa molecular, LogP, TPSA, HBD, HBA e **ligações rotacionáveis (RotB)**. Cada descritor é convertido em Z-score populacional no conjunto referência + lote; a Dist.FQ é a norma Euclidiana das diferenças dos seis Z-scores e é normalizada pelo maior valor observado em relação à referência. Os pKa continuam disponíveis na tabela físico-química do Mol.Sim, mas não entram na Dist.FQ. O endpoint `/bulk-compare` preserva o contrato anterior quando `show_chemical_space` não é enviado ou é `false`; quando `true`, retorna `chemical_space.points` com a referência, até 10 vizinhos, coordenadas MDS, similaridade MACCS, Dist.FQ normalizada, distância global, limite de exibição, quantidade total válida e stress MDS.

No Nitro.RA, o módulo `nitrosamine_space` produz dois espaços independentes. O primeiro é **PubChem — vizinhança relativa**: consulta até 40 CIDs, filtra N-nitroso, calcula a distância para todos os candidatos válidos e exibe as 10 menores distâncias. O segundo é **EMA Appendix 1 — referência**: usa a folha `N-nitrosamines` do snapshot local, deduplicada por SMILES canônico, e aplica um perfil z-score fixo versionado em `chemo_suite/apps/nitro_ra/data/ema_chemical_space_profile.json`. Os dois gráficos usam MACCS/Tanimoto, seis descritores, distância multimodal quadrática e MDS clássico, mas não se deve comparar diretamente os valores absolutos entre as duas referências.

Os gráficos são ferramentas de **triagem visual**: proximidade no plano representa a distância definida para aquela referência e não constitui prova de equivalência farmacológica, toxicológica ou regulatória. A proximidade de uma estrutura EMA não transfere automaticamente AI, categoria cPCA ou conclusão regulatória ao alvo. A interpretação deve considerar fonte, versão do snapshot, cobertura, scores, estruturas e propriedades físico-químicas.
