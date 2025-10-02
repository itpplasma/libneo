#!/usr/bin/env python3
"""Generate a static dashboard for magfie PNG artifacts."""

from __future__ import annotations

import argparse
import html
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
            rel_path = path.name
        destination = dest_root / rel_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        copied.append(rel_path if isinstance(rel_path, Path) else Path(rel_path))
    return copied


def _render_gallery(rel_paths: Iterable[Path]) -> str:
    items = []
    for rel_path in sorted(rel_paths):
        rel_posix = Path("images") / rel_path
        rel_str = rel_posix.as_posix()
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


def _build_html(branch: str, commit: str, run_id: str, repo: str,
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
    images_dir = output_dir / "images"

    if output_dir.exists():
        shutil.rmtree(output_dir)
    images_dir.mkdir(parents=True, exist_ok=True)

    copied = _copy_pngs(png_root, images_dir)

    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    gallery_html = _render_gallery(copied)
    html_text = _build_html(
        branch=args.branch,
        commit=args.commit[:12],
        run_id=args.run_id,
        repo=args.repo,
        timestamp=timestamp,
        gallery_html=gallery_html,
    )

    output_dir.joinpath("index.html").write_text(html_text, encoding="utf-8")


if __name__ == "__main__":
    main()
