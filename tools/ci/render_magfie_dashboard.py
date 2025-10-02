#!/usr/bin/env python3
"""Generate a static dashboard for test PNG artifacts.

The script keeps previously published branch dashboards when the caller
pre-populates ``--output`` with the current site contents (e.g. via wget).
It then refreshes the entry for the requested branch and rebuilds the
landing pages.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

# Common CSS for consistent styling across pages
_BASE_STYLE = """
body { font-family: Arial, sans-serif; margin: 2rem; }
a { color: #0366d6; text-decoration: none; }
a:hover { text-decoration: underline; }
code { font-family: Consolas, Monaco, monospace; background: #f6f8fa; padding: 2px 4px; border-radius: 3px; }
"""


def _get_pr_info(branch: str, repo: str) -> dict | None:
    """Fetch PR information for a branch using GitHub CLI.

    Returns dict with keys: number, url, title, draft, or None if no PR exists.
    """
    if not os.environ.get("GH_TOKEN"):
        return None

    try:
        result = subprocess.run(
            ["gh", "pr", "list", "--repo", repo, "--head", branch,
             "--json", "number,url,title,isDraft", "--jq", ".[0]"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None

        data = json.loads(result.stdout)
        if not data:
            return None

        return {
            "number": data.get("number"),
            "url": data.get("url"),
            "title": data.get("title"),
            "draft": data.get("isDraft", False),
        }
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError,
            json.JSONDecodeError, FileNotFoundError):
        return None


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
            f"<a href='{html.escape(rel_str)}' target='_blank'>"
            f"<img src='{html.escape(rel_str)}' alt='{html.escape(rel_str)}'"
            " loading='lazy'>"
            "</a>"
            f"<figcaption>{html.escape(rel_str)}</figcaption>"
            "</figure>"
        )
    if not items:
        return "<p>No PNG outputs were produced for this run.</p>"
    return "<div class='gallery'>" + "".join(items) + "</div>"


def _build_branch_html(branch: str, commit: str, run_id: str, repo: str,
                       timestamp: str, gallery_html: str, pr_info: dict | None) -> str:
    pr_section = ""
    if pr_info:
        pr_url = pr_info.get("url", "")
        pr_number = pr_info.get("number", "")
        pr_title = pr_info.get("title", "")
        pr_draft = pr_info.get("draft", False)
        draft_badge = " <span style='background:#6a737d;color:white;padding:2px 6px;border-radius:3px;font-size:0.85em'>DRAFT</span>" if pr_draft else ""
        pr_section = f"<p><strong>Pull Request:</strong> <a href='{html.escape(pr_url)}'>#{html.escape(str(pr_number))} {html.escape(pr_title)}</a>{draft_badge}</p>"

    return f"""<!DOCTYPE html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <title>libneo Test Dashboard – {html.escape(branch)}</title>
  <style>
    {_BASE_STYLE}
    .meta {{ margin-bottom: 1.5rem; }}
    .gallery {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 1.5rem; }}
    figure {{ margin: 0; }}
    figure a {{ display: block; cursor: pointer; }}
    figcaption {{ margin-top: 0.5rem; font-size: 0.9rem; word-break: break-word; }}
    img {{ width: 100%; height: auto; border: 1px solid #ccd; box-shadow: 0 2px 4px rgba(0,0,0,0.1); transition: box-shadow 0.2s, transform 0.2s; }}
    img:hover {{ box-shadow: 0 4px 8px rgba(0,0,0,0.2); transform: scale(1.02); }}
    .back-link {{ display: inline-block; margin-bottom: 1rem; }}
  </style>
</head>
<body>
  <a href='../index.html' class='back-link'>← Back to all branches</a>
  <h1>Test Dashboard – {html.escape(branch)}</h1>
  <div class='meta'>
    <p><strong>Branch:</strong> {html.escape(branch)}</p>
    {pr_section}
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

        status = "available" if info.get("has_pngs") else "no artifacts"

        pr_info = info.get("pr_info")
        branch_display = html.escape(branch)
        if pr_info:
            pr_number = pr_info.get("number", "")
            pr_url = pr_info.get("url", "")
            pr_draft = pr_info.get("draft", False)
            draft_badge = " <span style='background:#6a737d;color:white;padding:2px 4px;border-radius:3px;font-size:0.75em;font-weight:bold'>DRAFT</span>" if pr_draft else ""
            branch_display = f"{html.escape(branch)} (<a href='{html.escape(pr_url)}'>#{html.escape(str(pr_number))}</a>){draft_badge}"

        rows.append(
            "<tr>"
            f"<td><a href='{html.escape(link)}'>{branch_display}</a></td>"
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
  <title>libneo Test Dashboards</title>
  <style>
    {_BASE_STYLE}
    h1 {{ margin-bottom: 0.5rem; }}
    .subtitle {{ color: #586069; margin-bottom: 1rem; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 1rem; }}
    th, td {{ border: 1px solid #ccc; padding: 0.5rem; text-align: left; }}
    th {{ background-color: #f5f5f5; font-weight: 600; }}
  </style>
</head>
<body>
  <h1>Test Dashboards</h1>
  <p class='subtitle'>Automated test artifacts for all branches and pull requests</p>
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
    {_BASE_STYLE}
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
    pr_info = _get_pr_info(args.branch, args.repo)

    gallery_html = _render_gallery(copied)
    branch_html = _build_branch_html(
        branch=args.branch,
        commit=args.commit,
        run_id=args.run_id,
        repo=args.repo,
        timestamp=timestamp,
        gallery_html=gallery_html,
        pr_info=pr_info,
    )
    branch_path.joinpath("index.html").write_text(branch_html, encoding="utf-8")

    metadata_path = test_root / "branches.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    else:
        metadata = {}

    branch_metadata = {
        "path": args.branch,
        "updated": timestamp,
        "commit": args.commit,
        "has_pngs": bool(copied),
        "run_id": args.run_id,
        "repo": args.repo,
    }
    if pr_info:
        branch_metadata["pr_info"] = pr_info

    metadata[args.branch] = branch_metadata

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
