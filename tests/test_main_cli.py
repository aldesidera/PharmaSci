import json

import main as main_module
from chemo_suite import port_management
from main import main


def test_main_cli_cpca_isolated(capsys):
    exit_code = main([
        "--app",
        "nitro_ra",
        "--module",
        "cpca",
        "--smiles",
        "O=NN1CCCCC1",
        "--mdd-mg",
        "10",
    ])

    assert exit_code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "ok"
    assert payload["potency_category"] == 3
    assert payload["ai_ng_day"] == 400.0
    assert payload["ppm_limit"] == 40.0


def test_port_listener_pids_parses_ss_output(monkeypatch):
    monkeypatch.setattr(
        port_management.subprocess,
        "check_output",
        lambda *args, **kwargs: (
            'LISTEN 0 128 127.0.0.1:5000 0.0.0.0:* users:(("python",pid=321,fd=3))\n'
            'LISTEN 0 128 127.0.0.1:5000 0.0.0.0:* users:(("python",pid=321,fd=4))\n'
        ),
    )

    assert port_management._port_listener_pids(5000) == [321]


def test_stop_previous_server_only_kills_matching_project_process(monkeypatch, tmp_path):
    monkeypatch.setattr(port_management, "_port_listener_pids", lambda port: [101, 202])
    monkeypatch.setattr(port_management, "_is_pharmasci_server", lambda pid, root: pid == 101)
    monkeypatch.setattr(port_management.os, "kill", lambda pid, sig: None)

    stopped = port_management.stop_previous_pharmasci_server(5000, tmp_path)

    assert stopped == [101]


def test_web_entrypoint_delegates_to_shared_launcher(monkeypatch):
    calls = []
    monkeypatch.setattr(main_module, "run_molsim_web", lambda: calls.append(("run",)))

    assert main_module.main([]) == 0
    assert calls == [("run",)]


def test_shared_web_launcher_prepares_port_before_flask(monkeypatch):
    from types import SimpleNamespace
    import sys
    from chemo_suite import main as launcher

    calls = []
    fake_app = SimpleNamespace(run=lambda **kwargs: calls.append(("run", kwargs)))
    fake_logger = SimpleNamespace(info=lambda *args, **kwargs: None)
    monkeypatch.setitem(sys.modules, "app", SimpleNamespace(app=fake_app, logger=fake_logger))
    monkeypatch.setenv("MOLSIM_PORT", "5123")
    monkeypatch.setattr(launcher, "prepare_web_port", lambda port, project_dir: calls.append(("stop", port)) or [777])

    launcher.run_molsim_web()

    assert calls[0] == ("stop", 5123)
    assert calls[1][0] == "run"


def test_isolated_cpca_does_not_restart_port(monkeypatch):
    called = []
    monkeypatch.setattr(port_management, "prepare_web_port", lambda port, project_dir: called.append(port))

    assert main(["--app", "nitro_ra", "--module", "cpca", "--smiles", "O=NN1CCCCC1"]) == 0
    assert called == []
