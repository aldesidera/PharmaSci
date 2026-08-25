#!/bin/bash

echo "🧪 MolSim - Iniciando..."
echo ""

# Verificar se Python está instalado
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 não encontrado"
    exit 1
fi

# Verificar se as dependências estão instaladas
echo "📦 Verificando dependências..."
python3 -c "import flask; import rdkit; import fpdf" 2>/dev/null

if [ $? -ne 0 ]; then
    echo "📥 Instalando dependências..."
    pip install -r requirements.txt
fi

echo ""
echo "✅ Iniciando servidor (entrypoint único: main.py)..."
echo "🌐 Acesse: http://127.0.0.1:5000"
echo ""

python3 main.py
