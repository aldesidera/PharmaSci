"""Safe port management for PharmaSci web launchers."""

from __future__ import annotations

import os
import re
import signal
import subprocess
import time
from pathlib import Path
from typing import Optional

_PORT_PID_PATTERN = re.compile(r"pid=(\d+)")
_FALSE_VALUES = {"0", "false", "no", "off"}


def auto_restart_enabled() -> bool:
    """Return whether automatic replacement of an old web instance is enabled."""
    return os.environ.get("MOLSIM_AUTO_RESTART_PORT", "1").strip().lower() not in _FALSE_VALUES


def configured_port(default: int = 5000) -> int:
    """Read and validate the shared web port configuration."""
    value = os.environ.get("MOLSIM_PORT", str(default))
    port = int(value)
    if not 1 <= port <= 65535:
        raise ValueError("MOLSIM_PORT deve estar entre 1 e 65535")
    return port


def _port_listener_pids(port: int) -> list[int]:
    try:
        output = subprocess.check_output(
            ["ss", "-ltnp", f"sport = :{int(port)}"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return []

    pids: list[int] = []
    for match in _PORT_PID_PATTERN.finditer(output):
        pid = int(match.group(1))
        if pid not in pids:
            pids.append(pid)
    return pids


def _process_cmdline(pid: int) -> str:
    try:
        raw = Path(f"/proc/{pid}/cmdline").read_bytes()
    except (OSError, ValueError):
        return ""
    return raw.replace(b"\x00", b" ").decode(errors="replace").strip()


def _process_cwd(pid: int) -> Optional[Path]:
    try:
        return Path(os.readlink(f"/proc/{pid}/cwd")).resolve()
    except (OSError, ValueError):
        return None


def _is_pharmasci_server(pid: int, project_dir: Path) -> bool:
    if pid == os.getpid():
        return False
    return "main.py" in _process_cmdline(pid) and _process_cwd(pid) == project_dir.resolve()


def stop_previous_pharmasci_server(port: int, project_dir: Optional[Path] = None) -> list[int]:
    """Stop only older PharmaSci ``main.py`` instances listening on ``port``."""
    root = (project_dir or Path(__file__).resolve().parents[1]).resolve()
    stopped: list[int] = []
    for pid in _port_listener_pids(port):
        if not _is_pharmasci_server(pid, root):
            continue
        try:
            os.kill(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            continue
        stopped.append(pid)

    deadline = time.monotonic() + 2.0
    while stopped and time.monotonic() < deadline:
        if not any(Path(f"/proc/{pid}").exists() for pid in stopped):
            break
        time.sleep(0.05)

    for pid in stopped:
        if Path(f"/proc/{pid}").exists():
            try:
                os.kill(pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
    return stopped


def prepare_web_port(port: int, project_dir: Optional[Path] = None) -> list[int]:
    """Replace an older PharmaSci web process before any module is served."""
    if not auto_restart_enabled():
        return []
    return stop_previous_pharmasci_server(port, project_dir)
