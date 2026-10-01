#!/usr/bin/env python3
"""Exercise the production gate against an isolated concurrent-dispatch API."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import textwrap
import unittest


FAKE_GH = r'''#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys

p = Path(os.environ["GATE_FIXTURE"])
data = json.loads(p.read_text())
args = sys.argv[1:]
data.setdefault("calls", []).append(args)

def finish(output="", status=0):
    p.write_text(json.dumps(data))
    if output:
        print(output)
    sys.exit(status)

if args[0] == "api" or args[:2] == ["workflow", "run"]:
    fields = {}
    for index, arg in enumerate(args[:-1]):
        if arg in ("-f", "-F"):
            key, value = args[index + 1].split("=", 1)
            fields[key] = True if arg == "-F" and value == "true" else value
    if args[0] == "api":
        if fields.get("ref") != "main":
            finish("wrong downstream ref", 1)
        if fields.get("return_run_details") is not True:
            finish()
        if fields.get("inputs[libneo_ref]") != data["ref"]:
            finish("wrong candidate input", 1)
        if data.get("full") and fields.get("inputs[full]") != "true":
            finish("missing full-tier input", 1)
    elif fields.get("libneo_ref") != data["ref"]:
        finish("wrong candidate input", 1)
    failures = data.get("dispatch_failures", 0)
    if failures:
        data["dispatch_failures"] = failures - 1
        finish("transient dispatch error", 1)
    index = data.get("dispatched", 0)
    data["dispatched"] = index + 1
    run_id = data.get("returned_ids", [101 + index])[0]
    data.setdefault("runs", {})[str(run_id)] = data["verdicts"][index]
    if args[0] == "api":
        query = args[args.index("--jq") + 1]
        finish(str(run_id) if "workflow_run_id" in query and run_id else "")
    finish(f"https://github.com/itpplasma/Fixture/actions/runs/{run_id}")

if args[:2] == ["run", "list"]:
    # A different caller dispatches run999 while our own run101 is in flight.
    finish("999" if data.get("dispatched") else "100")

if args[:2] == ["run", "view"]:
    run_id = args[2]
    failures = data.get("status_failures", 0)
    if failures:
        data["status_failures"] = failures - 1
        finish("transient status error", 1)
    verdict = data.get("runs", {}).get(run_id, data["latest_verdict"])
    finish("completed " + verdict)

finish("unexpected gh operation", 1)
'''


class DownstreamGateTests(unittest.TestCase):
    def exercise(self, verdicts, latest="success", **options):
        workflow = (Path(__file__).resolve().parents[1] /
                    ".github/workflows/downstream-gate.yml")
        dispatch_step = workflow.read_text().split(
            "      - name: Dispatch downstream CI", 1
        )[1]
        production_script = textwrap.dedent(
            dispatch_step.split("        run: |\n", 1)[1]
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "gh"
            executable.write_text(FAKE_GH)
            executable.chmod(0o755)
            sleeper = root / "sleep"
            sleeper.write_text("#!/bin/sh\nexit 0\n")
            sleeper.chmod(0o755)
            (root / "ci").mkdir()
            blocking = options.pop("blocking", "yes")
            full = options.get("full", False)
            (root / "ci/downstreams").write_text(
                f"Fixture main.yml {blocking} {'yes' if full else 'no'}\n"
            )
            fixture = root / "fixture.json"
            fixture.write_text(json.dumps({
                "ref": "a" * 40, "verdicts": verdicts,
                "latest_verdict": latest, **options,
            }))
            script = root / "gate.sh"
            script.write_text(production_script)
            environment = dict(os.environ, PATH=f"{root}:{os.environ['PATH']}",
                               REF="a" * 40, FULL=str(full).lower(),
                               GATE_FIXTURE=str(fixture))
            result = subprocess.run(["bash", "-e", str(script)], cwd=root,
                                    env=environment, capture_output=True, text=True)
            return result, json.loads(fixture.read_text())

    def test_failed_created_run_cannot_borrow_concurrent_success(self):
        result, _ = self.exercise(["failure"])
        self.assertNotEqual(result.returncode, 0, result.stdout)

    def test_successful_created_run_ignores_concurrent_failure(self):
        result, _ = self.exercise(["success"], latest="failure")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_missing_or_invalid_created_id_fails_closed(self):
        for run_id in (None, 0, -1, "unknown"):
            with self.subTest(run_id=run_id):
                result, _ = self.exercise(["success"], returned_ids=[run_id])
                self.assertNotEqual(result.returncode, 0, result.stdout)

    def test_cancelled_created_run_uses_its_replacement(self):
        result, data = self.exercise(["cancelled", "success"])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(data["dispatched"], 2)

    def test_status_read_failure_cannot_become_success(self):
        result, _ = self.exercise(["failure"], status_failures=2)
        self.assertNotEqual(result.returncode, 0, result.stdout)

    def test_dispatch_retry_and_full_input_preserve_contract(self):
        result, data = self.exercise(["success"], dispatch_failures=2, full=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(data["dispatched"], 1)

    def test_report_only_failure_stays_nonblocking(self):
        result, _ = self.exercise(["failure"], blocking="no")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
