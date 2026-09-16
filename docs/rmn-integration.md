# Integração RMN baseada em RDKit e NMRShiftDB2

## Estado da implementação

O PharmaSci agora possui o endpoint `POST /rmn/predict`. Ele recebe um objeto JSON com `smiles` e, opcionalmente, `nucleus` — por exemplo, `13C` ou `1H`. O SMILES é validado e canonicalizado pelo RDKit antes da consulta ao NMRShiftDB2.

```json
{
  "smiles": "CCO",
  "nucleus": "13C"
}
```

A resposta inclui o SMILES canônico, o núcleo, a URL de origem, os identificadores de espectro/molécula, solvente, temperatura, frequência do equipamento quando informados e uma lista de picos com ppm, intensidade, multiplicidade e referências atômicas.

## Proveniência

O campo `provenance` pode ser `measured`, `hose-predicted` ou `unavailable`. O NMRShiftDB2 retorna um espectro medido quando há registro compatível; caso contrário, o endpoint público retorna uma previsão baseada em HOSE. A interface deve exibir essa distinção explicitamente e nunca rotular uma previsão HOSE como dado experimental.

## Exemplo de resposta resumida

```json
{
  "smiles_canonical": "CCO",
  "nucleus": "13C",
  "provenance": "measured",
  "spectrum_id": "nmrshiftdb10014832",
  "solvent": "Unreported",
  "peaks": [
    {"ppm": 18.6, "intensity": 0.78, "multiplicity": "Q", "atom_refs": ["a2"]},
    {"ppm": 57.4, "intensity": 1.0, "multiplicity": "T", "atom_refs": ["a1"]}
  ]
}
```

## Limites da primeira fase

Esta integração é um provedor de dados/predição, não um modelo próprio treinado pelo PharmaSci. O resultado depende da cobertura do NMRShiftDB2 e das condições experimentais disponíveis. Solvente, temperatura, frequência, estado físico, pH e método devem ser mantidos como metadados quando retornados. O conjunto SDF público é útil para índice estrutural e IDs, mas não deve ser tratado sozinho como tabela completa de deslocamentos: as listas de espectros são obtidas pelos endpoints de espectro/CML.

O próximo passo é criar a camada de visualização RMN no frontend e uma camada de cache local com TTL, preservando a URL de origem e a proveniência. Em seguida, deve ser criado um conjunto de validação por núcleo e classe química, com MAE/RMSE, cobertura de sinais e separação rigorosa entre moléculas de treino e teste antes de introduzir um modelo ML próprio.

## Referências

[1]: https://sourceforge.net/p/nmrshiftdb2/wiki/AutomationInterfaces/ "NMRShiftDB2 — Automation Interfaces"
[2]: https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdbhtml/using.html "NMRShiftDB2 — Using the database"
[3]: https://www.rdkit.org/docs/GettingStartedInPython.html "RDKit — Getting Started in Python"
