#!/usr/bin/env python3
"""Find the admitted, owned native artifacts for a queue-bot main push."""
import json
import os
from urllib.request import Request, urlopen


def publication_run(runs, tree, jobs, target_tree):
    # API order is newest first. Never borrow an older success after a newer
    # matching failure, nor artifacts from a different source tree.
    for run in runs:
        if tree(run['head_sha']) != target_tree:
            continue
        if run['conclusion'] != 'success':
            raise RuntimeError('The matching queue validation did not succeed')
        native = {'Build and test (Release)', 'Build and test (Debug)', 'Build Sphinx docs'}
        found = {job['name']: job for job in jobs(run['id']) if job['name'] in native}
        if set(found) != native:
            raise RuntimeError('Matching queue run is missing native/doc jobs')
        for job in found.values():
            if job['conclusion'] != 'success' or 'self-hosted' not in job['labels']:
                raise RuntimeError('Publication requires successful owned native/doc jobs')
        return run['id']
    raise RuntimeError('No completed queue validation matches the merged source tree')


def main():
    repo = os.environ['GITHUB_REPOSITORY']
    def api(path):
        req = Request('https://api.github.com/repos/' + repo + path,
                      headers={'Authorization': 'Bearer ' + os.environ['GH_TOKEN'],
                               'Accept': 'application/vnd.github+json'})
        with urlopen(req, timeout=30) as response:
            return json.load(response)
    def tree(sha):
        return api('/git/commits/' + sha)['tree']['sha']
    runs = api('/actions/workflows/main.yml/runs?event=merge_group&status=completed&per_page=100')['workflow_runs']
    run_id = publication_run(runs, tree,
                             lambda rid: api(f'/actions/runs/{rid}/jobs?per_page=100')['jobs'],
                             tree(os.environ['GITHUB_SHA']))
    with open(os.environ['GITHUB_OUTPUT'], 'a') as output:
        print(f'run-id={run_id}', file=output)
    print(f'Publishing exact-tree owned queue artifacts from run {run_id}')


if __name__ == '__main__':
    main()
