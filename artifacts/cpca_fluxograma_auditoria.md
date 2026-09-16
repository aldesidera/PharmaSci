# Auditoria estrutural do fluxograma CPCA

## Problemas identificados na versão anterior

1. O primeiro ramo estava rotulado como `Sim` para Categoria 5, mas no fluxo oficial a saída para Categoria 5 ocorre quando a resposta à primeira pergunta é `Não`.
2. Os ramos de Categoria 5 das três primeiras perguntas foram ligados por uma linha lateral compartilhada que atravessa o desenho e cria ambiguidade visual. A Figura 2 do EMA usa saídas laterais independentes, todas para a direita.
3. Os conectores verticais usavam posições antigas de losangos menores; depois do aumento das formas, começavam e terminavam dentro dos nós.
4. A seta inicial também terminava dentro do primeiro losango.
5. A saída de Categoria 1 e a revisão manual não estavam semanticamente ligadas ao nó `Score ≤ 1`.
6. Os rótulos `Sim/Não` não estavam sempre próximos do ramo correspondente.
7. O fluxo precisa reproduzir a ordem oficial: decisões no centro, resultado lateral à direita; não é necessário criar saídas à esquerda.

## Correção planejada

- usar coluna central de decisões;
- usar uma caixa lateral independente para cada resultado de cada decisão;
- ligar cada ramo lateral diretamente ao quadro correspondente;
- usar conectores verticais entre as decisões com início no vértice inferior e término no vértice superior;
- usar rótulos de ramo junto ao segmento horizontal;
- criar um nó explícito para `Score ≤ 1` e ligar `Sim` à Categoria 1;
- ligar `Não` desse nó a `Revisão manual necessária`;
- destacar apenas um caminho demonstrativo, sem alterar a semântica dos ramos.
