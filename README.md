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
- ✅ Nitro.RA: quatro abas de resultados independentes; cPCA funcional e Quantum/Metabolism preparados para evolução
- ✅ Nitro.RA cPCA: estrutura química, limite FDA, motivos do resultado e lookup do AI EMA Appendix 1
- ✅ Nitro.RA espaço químico: busca PubChem, filtro N-nitroso por RDKit, ranking Morgan2/Tanimoto, distância global de descritores e mapa PCA

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
9. Navegue pelas abas cPCA, Quantum, Metabolism e Nitrosaminas e Espaço Químico para consultar os resultados separadamente; módulos futuros ficam explicitamente identificados como em desenvolvimento.
10. Marque `Incluir busca e comparação de Nitrosaminas` quando quiser consultar, sob demanda, compostos semelhantes no PubChem e visualizar o espaço químico.
11. Informe a dose diária máxima em mg/dia quando desejar a conversão do AI para ppm.

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
- Módulos ainda não implementados retornam `status: not_implemented`, sem produzir valores fictícios
- Quando cPCA é selecionado, o resultado inclui `structure_svg`, `canonical_smiles` e o objeto `ema`
- O objeto `ema` consulta o snapshot `EMA/42261/2025 Rev.13`, atualizado em 24/06/2026, por SMILES canônico
- Uma estrutura ausente do Appendix 1 retorna `ema.status: not_listed`; o sistema não infere AI EMA a partir da categoria FDA
- `nitrosamine_space` consulta, sob demanda, a similaridade 2D do [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest), filtra `[N;X3][N;X2]=O`, seleciona até 10 candidatos e retorna descritores, Tanimoto, distância global e pontos PCA
- O PubChem fornece o lote de CIDs; o ranking final de Tanimoto é recalculado localmente com Morgan2, usando o pipeline de descritores do Mol.Sim, e a seleção final também considera a distância global normalizada
- Falhas de rede/timeout retornam `status: pubchem_unavailable`; nenhum resultado retorna `status: no_nitrosamines`; esses estados não bloqueiam o cPCA nem inventam candidatos

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

No modo **Batch**, o checkbox **Construir espaço químico 2D** habilita a geração de um gráfico com a referência e as moléculas válidas do lote. O posicionamento usa MDS clássico sobre uma distância composta por **60% de distância estrutural** — o complemento de Tanimoto ou Dice entre fingerprints — e **40% de distância físico-química normalizada**.

A componente físico-química reproduz a aba **Dist.FQ** do Excel e utiliza massa molecular, LogP, TPSA, HBD, HBA e **ligações rotacionáveis (RotB)**. Cada descritor é convertido em Z-score populacional no conjunto referência + lote; a Dist.FQ é a norma Euclidiana das diferenças absolutas desses seis Z-scores e é normalizada pelo maior valor observado em relação à referência. Os pKa continuam disponíveis na tabela físico-química do Mol.Sim, mas não entram na Dist.FQ desta réplica do modelo. O endpoint `/bulk-compare` preserva o contrato anterior quando o campo `show_chemical_space` não é enviado ou é `false`; quando `true`, retorna também `chemical_space.points`, com `x`, `y`, nome, SMILES, papel da molécula, similaridade em relação à referência, Dist.FQ normalizada e distância global.

O gráfico é uma ferramenta de **triagem visual**: proximidade no plano representa a distância composta definida acima e não constitui prova de equivalência farmacológica, toxicológica ou regulatória. A projeção deve ser interpretada junto com os scores, as estruturas e a tabela de propriedades físico-químicas.
