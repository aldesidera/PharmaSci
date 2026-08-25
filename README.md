# 🧪 MolSim v10 — build estabilizado

Versão ativa de trabalho do MolSim/MolSim_ver10, com foco em estabilização de API, validação de payloads, segurança operacional e consistência do fluxo pair/batch/report/export.

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
│           └── metabolism.py
├── requirements.txt
├── static/
├── templates/
│   ├── index.html
│   └── report_preview.html
├── tests/
│   ├── test_phase1_app.py
│   └── test_modular_structure.py
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
- O lookup de nome é tolerante a falha e não bloqueia o fluxo principal.

### Preview/PDF não abre
- Verifique popup blockers e use o botão de export diretamente no app.

---

## 🔗 Recursos úteis

- RDKit: https://www.rdkit.org/
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
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
