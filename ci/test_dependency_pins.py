#!/usr/bin/env python3
"""Exercise CMake dependency selection across reconfiguration and overrides."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class DependencyPinsTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        (self.root / "cmake").mkdir()
        helper = Path(__file__).resolve().parents[1] / "cmake/DependencyPins.cmake"
        helper = Path(os.environ.get("DEPENDENCY_PINS_HELPER", helper))
        shutil.copyfile(helper, self.root / "cmake/DependencyPins.cmake")
        (self.root / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.22)
project(pin_selection NONE)
include(cmake/DependencyPins.cmake)
file(WRITE "${CMAKE_BINARY_DIR}/selected" "${FORTIO_REF}\n${FORTNUM_REF}\n")
''')
        self.set_pins("1", "2")

    def set_pins(self, fortio, fortnum):
        (self.root / "fpm.toml").write_text(
            f'fortio = {{ rev = "{fortio * 40}" }}\n'
            f'fortnum = {{ rev = "{fortnum * 40}" }}\n'
        )

    def run_cmake(self, *arguments):
        result = subprocess.run(["cmake", *arguments], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return (self.root / "build/selected").read_text().splitlines()

    def test_manifest_update_reconfigures_existing_build(self):
        selected = self.run_cmake("-S", str(self.root), "-B", str(self.root / "build"))
        self.assertEqual(selected, ["1" * 40, "2" * 40])
        self.set_pins("3", "4")
        selected = self.run_cmake("--build", str(self.root / "build"))
        self.assertEqual(selected, ["3" * 40, "4" * 40])

    def test_explicit_override_survives_manifest_update(self):
        selected = self.run_cmake(
            "-S", str(self.root), "-B", str(self.root / "build"),
            "-DFORTIO_REF=review-candidate",
        )
        self.assertEqual(selected, ["review-candidate", "2" * 40])
        self.set_pins("3", "4")
        selected = self.run_cmake("--build", str(self.root / "build"))
        self.assertEqual(selected, ["review-candidate", "4" * 40])


if __name__ == "__main__":
    unittest.main()
