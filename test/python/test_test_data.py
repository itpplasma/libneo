"""test_data.download retries transient server errors and rejects altered files."""
import hashlib
import http.server
import threading

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import test_data  # noqa: E402
from test_data import download  # noqa: E402

PAYLOAD = b"pinned test data\n"
SHA = hashlib.sha256(PAYLOAD).hexdigest()


@pytest.fixture
def server(monkeypatch):
    monkeypatch.setattr(test_data.time, "sleep", lambda seconds: None)
    state = {"failures": 0, "body": PAYLOAD, "requests": 0}

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            state["requests"] += 1
            if state["requests"] <= state["failures"]:
                self.send_error(504)
                return
            self.send_response(200)
            self.end_headers()
            self.wfile.write(state["body"])

        def log_message(self, *args):
            pass

    httpd = http.server.HTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}/file", state
    httpd.shutdown()


def test_retries_gateway_timeouts(server, tmp_path):
    url, state = server
    state["failures"] = 2
    out = download(url, tmp_path / "data.nc", SHA)
    assert out.read_bytes() == PAYLOAD
    assert state["requests"] == 3


def test_rejects_altered_file(server, tmp_path):
    url, state = server
    state["body"] = b"truncated"
    with pytest.raises(RuntimeError, match="SHA-256"):
        download(url, tmp_path / "data.nc", SHA, attempts=2)
    assert not (tmp_path / "data.nc").exists()
    assert not list(tmp_path.glob("*.part"))


def test_replaces_stale_copy(server, tmp_path):
    url, state = server
    (tmp_path / "data.nc").write_bytes(b"stale")
    download(url, tmp_path / "data.nc", SHA)
    assert (tmp_path / "data.nc").read_bytes() == PAYLOAD
