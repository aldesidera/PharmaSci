# 🧪 MolSim Final v2.2

Versão única simplificada para análise de similaridade molecular com Flask backend e HTML/CSS/JS frontend profissional.

---

## 📋 Características

- ✅ **Fingerprints:** Morgan2 ECFP, Morgan2 FCFP, RDKit, MACCS
- ✅ **Métricas:** Tanimoto, Dice
- ✅ **Propriedades:** 7 propriedades físico-químicas principais
- ✅ **Exportação:** Preview HTML/CSS + impressão para PDF com ajuste visual
- ✅ **Interface:** Dark mode profissional com acentos vibrantes
- ✅ **Backend:** Flask puro
- ✅ **SMILES:** Campos vazios (sem exemplos pré-carregados)
- ✅ **Design:** Inspirado em softwares científicos profissionais

---

## 🚀 Instalação Rápida

### Linux/Mac

```bash
# 1. Extrair arquivo
tar -xzf molsim_final.tar.gz
cd molsim_final

# 2. Instalar dependências
pip install -r requirements.txt

# 3. Executar
python app.py

# 4. Abrir no navegador
http://localhost:5000
```

### Windows

```bash
# 1. Extrair arquivo
# Use WinRAR, 7-Zip ou similar para extrair molsim_final.tar.gz

# 2. Abrir PowerShell e navegar para a pasta
cd molsim_final

# 3. Instalar dependências
pip install -r requirements.txt

# 4. Executar
python app.py

# 5. Abrir navegador
http://localhost:5000
```

---

## 📁 Estrutura de Arquivos

```
molsim_final/
├── app.py                 # Backend Flask (endpoints REST)
├── analysis.py            # Lógica de análise molecular
├── requirements.txt       # Dependências Python
├── templates/
│   ├── index.html        # Frontend (dark mode profissional)
│   └── report_preview.html # Template de relatório imprimível
├── run.sh                # Script de inicialização (Linux/Mac)
└── README.md             # Este arquivo
```

---

## 🔧 Como Usar

1. **Insira SMILES** das duas moléculas (campos vazios)
2. **Digite nomes** (opcional) para identificar as moléculas
3. **Escolha Fingerprint:**
   - **Morgan2 ECFP:** Extended Connectivity (sem features)
   - **Morgan2 FCFP:** Feature-based (com features)
   - **RDKit:** Topológico padrão
   - **MACCS:** Keys estruturais (166 bits)
4. **Escolha Métrica:** Tanimoto ou Dice
5. **Clique "Analisar"**
6. **Veja os resultados:**
   - Similaridade e classificação
   - Estruturas moleculares (sem H)
   - Tabela de propriedades
   - Propriedades físico-químicas e pKa aproximado
7. **Abra o preview de relatório**, ajuste layout com sliders e use **Imprimir / Salvar PDF**

---

## 📊 Propriedades Físico-Químicas Calculadas

| Propriedade | Descrição | Unidade |
|------------|-----------|---------|
| Peso Molecular | Massa molar da molécula | g/mol |
| LogP | Coeficiente de partição (lipofilicidade) | - |
| Classe de Lipofilicidade | Classificação: Hidrofílica, Baixa, Moderada, Alta | - |
| LogS | Solubilidade (modelo ESOL) | log(mol/L) |
| Classe de Solubilidade | Classificação: Muito Solúvel, Solúvel, Pouco Solúvel, Insolúvel | - |
| TPSA | Área de Superfície Polar | Ų |
| Lipinski Ro5 | Passa no Rule of Five? | Sim/Não |

---

## 🧬 Fingerprints Disponíveis

### Morgan2 ECFP (Extended Connectivity)
- **Raio:** 2
- **Bits:** 2048
- **Features:** Não (apenas conectividade)
- **Uso:** Análise geral de similaridade estrutural
- **Melhor para:** Comparações rápidas e genéricas

### Morgan2 FCFP (Feature-based)
- **Raio:** 2
- **Bits:** 2048
- **Features:** Sim (aromáticos, H-donors, H-acceptors, etc)
- **Uso:** Análise de propriedades farmacofóricas
- **Melhor para:** Descoberta de fármacos e QSAR

### RDKit Topológico
- **Bits:** 2048
- **Tipo:** Topológico
- **Uso:** Referência padrão RDKit
- **Melhor para:** Validação e comparação com literatura

### MACCS Keys
- **Bits:** 166
- **Tipo:** Chaves estruturais predefinidas
- **Uso:** Análise de padrões estruturais
- **Melhor para:** Classificação de compostos

---

## 📊 Propriedades atuais do relatório

O relatório reúne propriedades relevantes para comparação estrutural e físico-química, incluindo:

- Massa Molecular (g/mol)
- pKa aproximado
- LogP
- TPSA
- HBD/HBA
- Ligações rotacionais

---

## 🎨 Design Profissional

Inspirado em softwares científicos profissionais como:
- ChemComp MOE
- DECTRIS Cloud
- NanoLabo
- PyMOL

**Características:**
- Dark mode com gradientes azul/preto
- Acentos em ciano e verde fluorescente
- Padrões geométricos de fundo
- Tipografia grande e impactante
- Animações suaves
- Responsivo para desktop, tablet e mobile

---

## 🧪 Exemplos de SMILES

```
Aspirina:        CC(=O)OC1=CC=CC=C1C(=O)O
Paracetamol:     CC(=O)NC1=CC=C(O)C=C1
Ibuprofeno:      CC(C)Cc1ccc(cc1)C(C)C(=O)O
Naproxeno:       COc1ccc2cc(ccc2c1)C(C)C(=O)O
Diclofenaco:     OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl
```

---

## 🐛 Troubleshooting

### "Port already in use"
```bash
# Usar porta diferente
python app.py --port 5001
```

### "ModuleNotFoundError: rdkit"
```bash
# Reinstalar RDKit
pip install rdkit-pypi numpy<2
```

### "Canvas is already in error state"
- Isso foi corrigido na v2.1
- Se persistir, limpe o cache do navegador (Ctrl+Shift+Delete)

### "SMILES inválido"
- Verifique a sintaxe SMILES em: https://www.chemspider.com/
- Ou use: https://pubchem.ncbi.nlm.nih.gov/

---

## 📥 Exportação PDF (via preview visual)

O preview de relatório permite ajustar via mouse (sliders) antes de salvar em PDF:
- ✅ Tamanho de título e texto
- ✅ Altura e espaçamento das caixas de moléculas
- ✅ Salvamento local do layout para próximos relatórios

Depois, use o botão **Imprimir / Salvar PDF** no próprio preview.

O conteúdo exportado contém:
- ✅ Método (Fingerprint + Métrica)
- ✅ Similaridade e classificação
- ✅ Estruturas moleculares (sem H)
- ✅ Tabela de propriedades
- ✅ Gráfico LogD vs pH

---

## 🔗 Recursos Úteis

- **RDKit:** https://www.rdkit.org/
- **ChemSpider:** https://www.chemspider.com/
- **PubChem:** https://pubchem.ncbi.nlm.nih.gov/
- **SMILES Tutorial:** https://en.wikipedia.org/wiki/Simplified_molecular_input_line_entry_system

---

## 📝 Versão e Histórico

**Versão:** 2.2 | **Status:** ✅ Pronto para Produção | **Data:** Março 2026

### Mudanças v2.2
- ✅ Corrigido imports (rdMolDescriptors)
- ✅ Adicionado Morgan2 ECFP e FCFP
- ✅ Corrigido erro do Canvas
- ✅ Botão PDF em posição correta
- ✅ Gráfico LogD vs pH com Chart.js

### Mudanças v2.1
- ✅ Dark mode profissional
- ✅ Acentos vibrantes (ciano + verde)
- ✅ Padrões geométricos
- ✅ Gráfico LogD vs pH interativo

### Mudanças v2.0
- ✅ Redesign completo
- ✅ Interface profissional
- ✅ Propriedades físico-químicas expandidas

---

## 📧 Suporte

Se encontrar problemas:
1. Verifique o arquivo `requirements.txt`
2. Reinstale as dependências: `pip install -r requirements.txt --force-reinstall`
3. Limpe o cache: `pip cache purge`
4. Tente novamente: `python app.py`

---

**Desenvolvido com ❤️ para análise molecular profissional**
