import logging
import os
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


def get_env_int(env_name: str, default: int, minimum: int = 1) -> int:
    raw_value = os.environ.get(env_name)
    if raw_value is None:
        return default
    try:
        parsed = int(raw_value)
    except ValueError:
        logger.warning("Valor inválido para %s: %s. Usando %s.", env_name, raw_value, default)
        return default
    if parsed < minimum:
        logger.warning("Valor de %s abaixo do mínimo (%s). Usando %s.", env_name, minimum, default)
        return default
    return parsed


def parse_cors_origins() -> List[str]:
    raw_origins = os.environ.get("MOLSIM_CORS_ORIGINS", "").strip()
    if not raw_origins:
        return [r"^https?://(127\.0\.0\.1|localhost)(:\d+)?$"]
    origins = [origin.strip() for origin in raw_origins.split(",") if origin.strip()]
    return origins or [r"^https?://(127\.0\.0\.1|localhost)(:\d+)?$"]


def make_error_payload(status_code: int, message: str, field: Optional[str] = None, code: str = "invalid_request") -> Dict[str, Any]:
    payload: Dict[str, Any] = {"error": {"code": code, "message": message}}
    if field:
        payload["error"]["field"] = field
    return payload
