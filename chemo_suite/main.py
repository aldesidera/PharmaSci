import os
from typing import Optional


def run_molsim_web(host: Optional[str] = None, port: Optional[int] = None) -> None:
    from app import app, logger

    project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    effective_host = host or os.environ.get("MOLSIM_HOST", "127.0.0.1")
    effective_port = port or int(os.environ.get("MOLSIM_PORT", "5000"))
    logger.info("📦 Runtime ativo: MolSim_ver10 (%s)", project_dir)
    logger.info("🚀 MolSim v3.3 iniciando em http://%s:%s", effective_host, effective_port)
    app.run(debug=False, host=effective_host, port=effective_port)
