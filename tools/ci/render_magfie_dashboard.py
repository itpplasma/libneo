#!/usr/bin/env python3
"""Generate a static dashboard for magfie PNG artifacts.

The script keeps previously published branch dashboards when the caller
pre-populates ``--output`` with the current site contents (e.g. via wget).
It then refreshes the entry for the requested branch and rebuilds the
landing pages.
"""

from __future__ import annotations

import argparse
import html
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


def _copy_pngs(src_root: Path, dest_root: Path) -> list[Path]:
    copied: list[Path] = []
    if not src_root.exists():
        return copied

    for path in sorted(src_root.rglob("*.png")):
        if not path.is_file():
            continue
        try:
            rel_path = path.relative_to(src_root)
        except ValueError:
            rel_path = Path(path.name)
        destination = dest_root / rel_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        copied.append(rel_path if isinstance(rel_path, Path) else Path(rel_path))
    return copied


def _render_gallery(rel_paths: Iterable[Path]) -> str:
    items = []
    for rel_path in sorted(rel_paths):
        rel_str = Path("images") / rel_path
        rel_str = rel_str.as_posix()
        items.append(
            "<figure>"
            f"<img src='{html.escape(rel_str)}' alt='{html.escape(rel_str)}'"
            " loading='lazy'>"
            f"<figcaption>{html.escape(rel_str)}</figcaption>"
            "</figure>"
        )
    if not items:
        return "<p>No PNG outputs were produced for this run.</p>"
    return "<div class='gallery'>" + "".join(items) + "</div>"


def _build_branch_html(branch: str, commit: str, run_id: str, repo: str,
                       timestamp: str, gallery_html: str) -> str:
    return f"""<!DOCTYPE html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <title>libneo magfie PNG dashboard</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; }}
    a {{ color: #0366d6; }}
    .meta {{ margin-bottom: 1.5rem; }}
    .gallery {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 1.5rem; }}
    figure {{ margin: 0; }}
    figcaption {{ margin-top: 0.5rem; font-size: 0.9rem; word-break: break-word; }}
    img {{ width: 100%; height: auto; border: 1px solid #ccd; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
  </style>
</head>
<body>
  <h1>libneo magfie PNG dashboard</h1>
  <div class='meta'>
    <p><strong>Branch:</strong> {html.escape(branch)}</p>
    <p><strong>Commit:</strong> <code>{html.escape(commit)}</code></p>
    <p><strong>Workflow run:</strong> <a href='https://github.com/{html.escape(repo)}/actions/runs/{html.escape(run_id)}'>{html.escape(run_id)}</a></p>
    <p><strong>Generated:</strong> {html.escape(timestamp)}</p>
  </div>
  {gallery_html}
</body>
</html>
"""


def _build_test_overview(metadata: dict[str, dict[str, str]], timestamp: str) -> str:
    rows = []
    for branch, info in sorted(metadata.items()):
        link = f"{info['path'].rstrip('/')}/"
        updated = info.get("updated", "-")
        commit = info.get("commit", "-")[:12]
        run_id = info.get("run_id", "-")
        repo = info.get("repo", "")
        if repo and run_id and run_id != "-":
            workflow_link = f"https://github.com/{html.escape(repo)}/actions/runs/{html.escape(run_id)}"
            workflow_text = f"run {html.escape(run_id)}"
            workflow_cell = f"<a href='{workflow_link}'>{workflow_text}</a>"
        else:
            workflow_cell = "-"
        status = "available" if info.get("has_pngs") else "missing"
        rows.append(
            "<tr>"
            f"<td><a href='{html.escape(link)}'>{html.escape(branch)}</a></td>"
            f"<td>{html.escape(updated)}</td>"
            f"<td><code>{html.escape(commit)}</code></td>"
            f"<td>{workflow_cell}</td>"
            f"<td>{html.escape(status)}</td>"
            "</tr>"
        )

    if not rows:
        rows.append("<tr><td colspan='5'>No branch dashboards published yet.</td></tr>")

    table_rows = "\n".join(rows)
    return f"""<!DOCTYPE html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <title>libneo Branch Dashboards</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 1rem; }}
    th, td {{ border: 1px solid #ccc; padding: 0.5rem; text-align: left; }}
    th {{ background-color: #f5f5f5; }}
    code {{ font-family: Consolas, monospace; }}
  </style>
</head>
<body>
  <h1>Branch Dashboards</h1>
  <p>Generated: {html.escape(timestamp)}</p>
  <table>
    <thead>
      <tr><th>Branch</th><th>Updated (UTC)</th><th>Commit</th><th>Workflow</th><th>Status</th></tr>
    </thead>
    <tbody>
{table_rows}
    </tbody>
  </table>
</body>
</html>
"""


def _build_root_stub(timestamp: str) -> str:
    return f"""<!DOCTYPE html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <title>libneo GitHub Pages</title>
  <meta http-equiv='refresh' content='0; url=./test/' />
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; }}
  </style>
</head>
<body>
  <h1>libneo GitHub Pages</h1>
  <p>This site hosts automatically generated test dashboards.</p>
  <p>You will be redirected to <a href='./test/'>/test/</a> momentarily.</p>
  <p>Generated: {html.escape(timestamp)}</p>
</body>
</html>
"""


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png-root", default="png-artifacts",
                        help="Directory that holds the downloaded PNG artifacts")
    parser.add_argument("--output", default="dashboard",
                        help="Output directory for the generated dashboard")
    parser.add_argument("--branch", required=True, help="Branch name for metadata")
    parser.add_argument("--commit", required=True, help="Commit SHA for metadata")
    parser.add_argument("--run-id", required=True, help="Workflow run identifier")
    parser.add_argument("--repo", required=True, help="<owner>/<repo> for links")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()

    png_root = Path(args.png_root)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    test_root = output_dir / "test"
    branch_path = test_root / Path(args.branch)
    branch_images = branch_path / "images"

    branch_exists = branch_path.exists()
    if branch_exists:
        shutil.rmtree(branch_path)
    branch_images.mkdir(parents=True, exist_ok=True)

    copied = _copy_pngs(png_root, branch_images)

    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    gallery_html = _render_gallery(copied)
    branch_html = _build_branch_html(
        branch=args.branch,
        commit=args.commit,
        run_id=args.run_id,
        repo=args.repo,
        timestamp=timestamp,
        gallery_html=gallery_html,
    )
    branch_path.joinpath("index.html").write_text(branch_html, encoding="utf-8")

    metadata_path = test_root / "branches.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    else:
        metadata = {}

    metadata[args.branch] = {
        "path": args.branch,
        "updated": timestamp,
        "commit": args.commit,
        "has_pngs": bool(copied),
        "run_id": args.run_id,
        "repo": args.repo,
    }

    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")

    test_overview_html = _build_test_overview(metadata, timestamp)
    test_root.joinpath("index.html").write_text(test_overview_html, encoding="utf-8")

    root_index = output_dir / "index.html"
    if root_index.exists():
        existing = root_index.read_text(encoding="utf-8")
        if "libneo Test Dashboards" in existing:
            root_index.unlink()
    if not root_index.exists():
        root_index.write_text(_build_root_stub(timestamp), encoding="utf-8")


if __name__ == "__main__":
    main()
