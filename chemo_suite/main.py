import os
from pathlib import Path
from typing import Optional

from .port_management import configured_port, prepare_web_port


def run_molsim_web(host: Optional[str] = None, port: Optional[int] = None) -> None:
    """Start the shared Flask application for every web module."""
    from app import app, logger

    project_dir = Path(__file__).resolve().parents[1]
    effective_host = host or os.environ.get("MOLSIM_HOST", "127.0.0.1")
    effective_port = port or configured_port()
    stopped = prepare_web_port(effective_port, project_dir)
    if stopped:
        logger.info(
            "Instância anterior encerrada na porta %s: %s",
            effective_port,
            ", ".join(map(str, stopped)),
        )

    logger.info("📦 Runtime ativo: MolSim_ver10 (%s)", project_dir)
    logger.info("🚀 MolSim v3.3 iniciando em http://%s:%s", effective_host, effective_port)
    app.run(debug=False, host=effective_host, port=effective_port)
