#!/usr/bin/env python3
import unittest
from publication_source import publication_run


class PublicationSource(unittest.TestCase):
    def setUp(self):
        self.runs = [{'id': 2, 'head_sha': 'other', 'conclusion': 'success'},
                     {'id': 1, 'head_sha': 'matching', 'conclusion': 'success'}]
        self.jobs = [{'name': name, 'conclusion': 'success', 'labels': ['self-hosted', 'itp-build']}
                     for name in ['Build and test (Release)', 'Build and test (Debug)', 'Build Sphinx docs']]

    def select(self):
        return publication_run(self.runs, lambda sha: sha, lambda _: self.jobs, 'matching')

    def test_exact_tree_selects_owned_run(self):
        self.assertEqual(self.select(), 1)

    def test_other_tree_cannot_supply_artifacts(self):
        self.runs.pop()
        with self.assertRaises(RuntimeError):
            self.select()

    def test_hosted_native_artifacts_rejected(self):
        self.jobs[0]['labels'] = ['ubuntu-latest']
        with self.assertRaises(RuntimeError):
            self.select()

    def test_failure_cannot_borrow_older_success(self):
        self.runs.insert(0, {'id': 3, 'head_sha': 'matching', 'conclusion': 'failure'})
        with self.assertRaises(RuntimeError):
            self.select()


if __name__ == '__main__':
    unittest.main()
