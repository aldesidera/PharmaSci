# Política de cache e fallback offline de RMN

O Mol.FQ utiliza um cache local persistente para resultados do NMRShiftDB2. A chave é formada pelo SMILES canonicalizado pelo RDKit e pelo núcleo analisado, por exemplo `CCO|13C`.

## Padrões

- TTL do cache fresco: 900 segundos (15 minutos).
- Janela máxima de fallback stale-if-error: 604800 segundos (7 dias).
- Arquivo padrão: `.cache/nmrshiftdb.json`.
- Caminho alternativo: variável `PHARMASCI_NMR_CACHE_PATH`.
- TTL alternativo: `PHARMASCI_NMR_CACHE_TTL_SECONDS`.
- Janela stale alternativa: `PHARMASCI_NMR_STALE_MAX_AGE_SECONDS`.

## Estados

`live` significa que a resposta veio de uma consulta atual ao NMRShiftDB2. `fresh-cache` significa que a consulta externa não foi repetida porque o resultado ainda estava dentro do TTL. `stale-cache` significa que o provedor estava indisponível e o resultado anterior foi reutilizado dentro da janela de fallback.

O fallback preserva a proveniência original do espectro: `measured` continua sendo dado medido e `hose-predicted` continua sendo previsão HOSE. O fato de o resultado ter vindo do cache é informado separadamente por `cache_status`.

## Segurança operacional

O cache é gravado de forma atômica por arquivo temporário e `os.replace`, evitando deixar um JSON parcialmente escrito. O conteúdo deve ser considerado dado derivado de uma fonte externa e não substitui a confirmação da origem, do solvente, da temperatura, da frequência e do método experimental.

Quando não existe entrada válida no cache e o NMRShiftDB2 está indisponível, o Mol.FQ retorna o núcleo como `unavailable` e não fabrica sinais. Quando o resultado está fora da janela stale, a consulta falha de forma explícita e o usuário é informado.

A configuração pode desabilitar a persistência definindo `PHARMASCI_NMR_CACHE_TTL_SECONDS=0`, mantendo a consulta sem armazenamento local.
