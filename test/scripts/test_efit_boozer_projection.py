"""Exercise the actual converter projection with an analytic noncircular map."""
import math
import pathlib
import subprocess
import sys
import tempfile

binary = pathlib.Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory(dir=pathlib.Path.cwd()) as temporary:
    root = pathlib.Path(temporary)
    (root / "efit_to_boozer.inp").write_text("1\n4\n16384\n4\n3\n32\n1\n")
    subprocess.run([str(binary)], cwd=root, check=True, capture_output=True)
    lines = iter((root / "fromefit_neo_lhs.bc").read_text().splitlines())
    for _ in range(6):
        next(lines)
    for surface in range(3):
        next(lines); next(lines)
        profile = list(map(float, next(lines).split()))
        s = profile[0]
        next(lines)
        for mode in range(33):
            values = list(map(float, next(lines).split()))[2:]
            expected = [0.0] * 8
            if mode == 0:
                expected[0] = 6.2
                expected[6] = 1e-4
            elif mode == 1:
                expected[0] = .9 * math.sqrt(s)
                expected[3] = 1.4 * math.sqrt(s)
                expected[6] = .15e-4
            elif mode == 2:
                expected[0] = .07 * s
                expected[2] = .09 * s
                expected[7] = .05e-4
            elif mode == 3:
                expected[1] = .04 * s
                expected[3] = -.025 * s
            # The displacement itself is defined in the original angle;
            # only independently known R/Z/B coefficients are checked.
            for index in [0, 1, 2, 3, 6, 7]:
                tolerance = 5e-8 if index < 4 else 5e-13
                assert abs(values[index] - expected[index]) < tolerance, (mode, index, values[index], expected[index])
        # Keep the independently defined physical covariant-current integral.
        average = 0.0
        for j in range(1, 16385):
            theta = 2 * math.pi * j / 16384
            angle = theta + .18 * math.sin(theta) + .035 * math.sin(2 * theta)
            jac = 1 + .18 * math.cos(theta) + .07 * math.cos(2 * theta)
            weight = jac * (1 + 1e-5 * (.4 + math.cos(17 * theta)))
            b = 1 + .15 * math.cos(angle) + .05 * math.sin(2 * angle)
            average += weight / b**2
        expected_volume = -(average / 16384) * 1e-6 * (2 * math.pi)**2
        rounding = 5e-9 * abs(expected_volume) + 2e-16
        assert abs(profile[5] - expected_volume) < rounding, (profile[5], expected_volume)
print("Noncircular phase-Jacobian projection and physical-current integral: PASS")
