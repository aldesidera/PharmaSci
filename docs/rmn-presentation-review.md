# Revisão visual da apresentação RMN

A validação no navegador do Mol.FQ confirmou que a nova seção RMN apresenta um resumo executivo com sinais retornados, núcleos disponíveis e origem temporal, seguido de painéis separados para ¹H e ¹³C.

A hierarquia visual ficou adequada ao tema verde do Mol.FQ, com badges distintos para previsão HOSE e consulta atual. A primeira renderização ainda mostrava valores de ppm com precisão excessiva, causando rolagem horizontal na tabela de ¹³C. Também havia repetição do aviso de ausência de registro medido, uma vez por núcleo.

A correção aplicada limita os valores numéricos do frontend a três casas decimais e deduplica os avisos agregados pelo motor. A validação automatizada posterior passou com 116 testes.

A versão do PDF com espectros foi renderizada em quatro páginas. A página de ¹H mostra o resumo, badges, gráfico de bastões com eixo invertido de 8,53 a 1,26 ppm e tabela com quatro casas decimais. A página anterior permanece legível, sem colisões de layout. A página de ¹³C será verificada na mesma rodada para confirmar duas casas decimais e paginação.
