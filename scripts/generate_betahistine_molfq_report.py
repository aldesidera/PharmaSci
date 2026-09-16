from pathlib import Path
from app import app
from weasyprint import HTML

out_dir = Path('/home/ubuntu/PharmaSci_git/artifacts/molfq_betahistine')
out_dir.mkdir(parents=True, exist_ok=True)
smiles = 'CNCCC1=CC=CC=N1'
analyses = ['summary', 'structural', 'ionization', 'lipophilicity', 'solubility', 'nmr', 'geometry_3d']
client = app.test_client()
response = client.post('/mol-fq/report-preview', json={'smiles': smiles, 'name': 'Betahistina', 'analyses': analyses})
if response.status_code != 200:
    raise RuntimeError(f'Falha ao gerar preview: {response.status_code} {response.get_data(as_text=True)}')
html_path = out_dir / 'relatorio_molfq_betahistina.html'
pdf_path = out_dir / 'relatorio_molfq_betahistina.pdf'
html_path.write_bytes(response.data)
HTML(string=response.get_data(as_text=True), base_url='/home/ubuntu/PharmaSci_git').write_pdf(str(pdf_path))
print(html_path)
print(pdf_path)
