# Validação visual do módulo RMN

A interface foi aberta a partir da versão atualizada do PharmaSci e o chip **RMN** foi exibido corretamente ao lado de Mol.Sim e Nitro.RA.

Com o SMILES `CCO` e núcleo `13C`, o frontend exibiu a badge verde **Experimental**, os metadados do NMRShiftDB2, dois picos e o gráfico vertical de deslocamentos químicos. O retorno veio com `provenance=measured` e espectro `nmrshiftdb10014832`.

Com o SMILES `C1=CN(C=N1)N=O` e núcleo `13C`, o frontend exibiu a badge laranja **Predição HOSE**, três picos, o SMILES canônico `O=Nn1ccnc1`, ausência explícita de solvente/temperatura e o aviso de que não foi encontrado registro medido. O gráfico foi renderizado normalmente.

A distinção visual entre experimental e HOSE está funcional e não houve deslocamento do resultado para outro módulo durante a alternância.
