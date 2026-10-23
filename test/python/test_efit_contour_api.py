"""Public contour API checked against an analytic Solov'ev equilibrium.

The polynomial solves Delta*psi = -mu0 R^2 p' - FF' exactly.  Its contour
quadrature below uses neither the native splines nor the Boozer construction.
All files are generated synthetic data, with the same polynomial as the native
test_efit_sfl test, rather than machine or experiment fixtures.
"""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest


from python.tests.efit_solovev import (
    R0, B0, C, PSI_EDGE, _psi, _grad, _axis, _contour_oracle, _write_solovev,
)


@pytest.fixture
def boozer_api(monkeypatch):
    extension_path = os.environ.get("LIBNEO_EFIT_EXTENSION")
    if extension_path:
        spec = importlib.util.spec_from_file_location("_efit_to_boozer", extension_path)
        extension = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(extension)
    else:
        extension = pytest.importorskip("_efit_to_boozer")
    return _load_api(monkeypatch, extension)


def _load_api(monkeypatch, extension):
    monkeypatch.setitem(sys.modules, "_efit_to_boozer", extension)
    source = Path(__file__).resolve().parents[2]/"python/efit_to_boozer/boozer.py"
    spec = importlib.util.spec_from_file_location("contour_api_under_test", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("psi_sign,f_sign", [(1, 1), (-1, 1), (1, -1), (-1, -1)])
def test_contour_signed_field_oracle(tmp_path, boozer_api, psi_sign, f_sign):
    from scipy.integrate import quad

    gfile = tmp_path/"solovev.g"
    _write_solovev(gfile, psi_sign, f_sign)
    transform, metadata = boozer_api.get_boozer_transform(
        np.array([.15, .65, 1.0]), 32, gfile=gfile,
        n_theta=128, n_q=32, return_metadata=True)
    assert len(transform) == 7
    np.testing.assert_allclose(metadata["axis"], _axis(), atol=2e-8, rtol=0)
    radius, dtheta, _, theta_b, theta_geom, _, field = transform
    assert theta_b.shape == theta_geom.shape == (2, 65)
    assert np.all(np.diff(theta_b, axis=1) > 0)
    np.testing.assert_allclose(theta_b[:, -1]-theta_b[:, 0], 2*np.pi, atol=1e-14)
    psi_axis = _psi(*_axis())
    toroidal_flux = quad(lambda p: _contour_oracle(p)[0], psi_axis, PSI_EDGE)[0]
    assert metadata["psitor_max"] == pytest.approx(
        psi_sign*toroidal_flux*1e8, rel=3e-6)
    for k in range(2):
        angle = theta_b[k, :-1]
        geom = angle+dtheta[k](angle)
        r = metadata["axis"][0]+radius[k](angle)*np.cos(geom)
        z = metadata["axis"][1]+radius[k](angle)*np.sin(geom)
        psi = _psi(r, z)
        # Check the returned geometry independently before using its flux value.
        assert np.ptp(psi)/(PSI_EDGE-_psi(*_axis())) < 2e-7
        q, current_i, current_g = _contour_oracle(np.mean(psi), psi_sign, f_sign)
        partial_flux = quad(lambda p: _contour_oracle(p)[0], psi_axis, np.mean(psi))[0]
        assert partial_flux/toroidal_flux == pytest.approx([.15, .65][k], rel=3e-6)
        assert metadata["iota"][k] == pytest.approx(1/q, rel=3e-7)
        assert metadata["curr_i"][k] == pytest.approx(current_i*1e6, rel=3e-7)
        assert metadata["curr_g"][k] == pytest.approx(current_g*1e6, rel=3e-7)
        pr, pz = _grad(r, z)
        expected_b = np.sqrt(pr*pr+pz*pz+current_g**2)/r
        np.testing.assert_allclose(field[k](angle)*1e-4, expected_b, rtol=3e-7)


def test_contour_quality_error_is_recoverable(tmp_path, boozer_api):
    gfile = tmp_path/"inconsistent_pressure.g"
    _write_solovev(gfile, pprime_scale=.7)
    args = (np.array([.3, 1.0]), 16)
    kwargs = dict(gfile=gfile, n_theta=128, n_q=32)
    with pytest.raises(ValueError, match="Grad-Shafranov"):
        boozer_api.get_boozer_transform(*args, **kwargs)
    transform, metadata = boozer_api.get_boozer_transform(
        *args, **kwargs, allow_poor_equilibrium=True,
        override_reason="synthetic inconsistent profile regression", return_metadata=True)
    assert len(transform[0]) == 1
    assert np.isfinite(metadata["iota"]).all()
    assert metadata["overridden"]
    assert metadata["status_code"] != 0
    # A rejected initialization must not prevent the next valid conversion.
    _write_solovev(gfile)
    valid = boozer_api.get_boozer_transform(*args, **kwargs)
    assert len(valid[0]) == 1


@pytest.mark.parametrize("selection", ["method", "positional", "deck"])
def test_explicit_legacy_dispatch(tmp_path, monkeypatch, selection):
    monkeypatch.chdir(tmp_path)
    calls = []
    extension = SimpleNamespace(efit_to_boozer=SimpleNamespace(
        init=lambda **kwargs: calls.append(kwargs)))
    api = _load_api(monkeypatch, extension)

    def legacy_angles(s, num_theta):
        calls.append(float(s))
        angle = np.linspace(0, 2*np.pi, 2*num_theta+1)
        return (np.ones_like(angle), angle, angle, angle[:-1],
                np.zeros_like(angle), np.full_like(angle, 2e4))

    monkeypatch.setattr(api, "_get_legacy_angles_G_and_B0", legacy_angles)
    args, kwargs = [np.array([.3, 1.0]), 8], {}
    if selection == "method":
        kwargs["method"] = "legacy"
    elif selection == "positional":
        args.append(2)
    else:
        (tmp_path/"efit_to_boozer.inp").write_text("100\n32\n64\n64\n32\n12\n1e10\n2\n")
        args.append(0)
    result = api.get_boozer_transform(*args, **kwargs)
    assert len(result) == 7
    assert calls == [{"method": "legacy"}, .3]
    assert result[0][0](1.0) == pytest.approx(1.0)


def test_zero_gs_sources_require_explicit_flux_units(tmp_path, boozer_api):
    from libneo.eqdsk_base import read_eqdsk, write_eqdsk

    gfile = tmp_path/"zero_gs_sources.g"
    _write_solovev(gfile)
    equilibrium = read_eqdsk(gfile)
    for profile in ("ptotprof", "dpressdpsiprof", "fdfdpsiprof"):
        equilibrium[profile].fill(0.0)
    write_eqdsk(gfile, equilibrium)
    options = dict(n_theta=128, n_q=32, allow_poor_equilibrium=True,
                   override_reason="synthetic zero GS source regression")
    with pytest.raises(ValueError, match="cannot infer poloidal-flux units"):
        boozer_api.get_contour_transform(gfile, [.3], 16, **options)
    metadata = boozer_api.get_contour_transform(gfile, [.3], 16, exp_bpol=0, **options)
    assert metadata["overridden"]
    assert np.all(metadata["iota"] > 0)
    assert np.isfinite(metadata["spectra"]).all()


@pytest.mark.parametrize("record,expected", [
    ("1e10", None), ("1e300", None), ("-1e10", -1e10), ("2e10", 2e10),
    ("-2e10", -2e10), ("-1e300", -1e300),
    ("NaN", "error"), ("Infinity", "error"), ("-Infinity", "error"),
])
def test_deck_flux_bound_sentinels(tmp_path, monkeypatch, record, expected):
    monkeypatch.chdir(tmp_path)
    api = _load_api(monkeypatch, SimpleNamespace())
    (tmp_path/"field_divB0.inp").write_text("0\n1\n1.0\n72\n0.99\n4\n'field.g'\n")
    (tmp_path/"efit_to_boozer.inp").write_text(f"100\n8\n32\n16\n8\n6\n{record}\n")
    if expected == "error":
        with pytest.raises(ValueError):
            api._contour_input(None, None)
    else:
        gfile, bound = api._contour_input(None, None)
        assert gfile == str(tmp_path/"field.g")
        assert bound == expected


@pytest.mark.parametrize("helper", ["get_boozer_harmonics_divide_f_by_B0",
                                    "get_boozer_harmonics_divide_f_by_B0_1D"])
def test_harmonic_division_preserves_equilibrium_options(monkeypatch, helper):
    api = _load_api(monkeypatch, SimpleNamespace())
    stor, field, calls = np.array([.2, .6, 1.0]), np.array([2.0, 4.0]), []
    options = dict(method="legacy", allow_poor_equilibrium=True,
                   override_reason="synthetic harmonic regression")

    def fake_field(s, num_theta, **kwargs):
        np.testing.assert_array_equal(s, stor)
        assert num_theta == 8
        calls.append(kwargs)
        return [lambda angle: field[0], lambda angle: field[1]]

    monkeypatch.setattr(api, "get_B0_of_s_theta_boozer", fake_field)
    zero = [lambda angle: 0.0]*2
    f = lambda s, theta, phi: (s+1)*(3+2*np.cos(theta))
    result = getattr(api, helper)(f, stor, 8, 1, 1, 0, zero, zero, **options)
    assert calls == [options]
    expected = ((stor[:-1]+1)/field)[:, None]*np.array([1.0, 3.0, 1.0])
    np.testing.assert_allclose(result, expected, atol=1e-14)


@pytest.mark.parametrize("psi_sign,f_sign", [(1, 1), (-1, 1), (1, -1), (-1, -1)])
def test_chartmap_signed_field_oracle(tmp_path, monkeypatch, boozer_api, psi_sign, f_sign):
    netcdf = pytest.importorskip("netCDF4")
    from scipy.integrate import quad

    monkeypatch.setitem(sys.modules, "efit_to_boozer.boozer", boozer_api)
    from libneo.eqdsk_to_boozer_chartmap import convert_eqdsk_to_chartmap

    gfile, output = tmp_path/"solovev.g", tmp_path/"chartmap.nc"
    _write_solovev(gfile, psi_sign, f_sign)
    torflux = convert_eqdsk_to_chartmap(
        gfile, output, nrho=12, ntheta=16, nzeta=8, ntheta_int=128, n_q=32)
    with netcdf.Dataset(output) as dataset:
        x, y, z, bmod = [np.asarray(dataset.variables[name][:])
                          for name in ("x", "y", "z", "Bmod")]
        r, z = np.hypot(x, y)*.01, z*.01
        psi = _psi(r, z)
        assert np.max(np.ptp(psi, axis=(0, 1)))/(PSI_EDGE-_psi(*_axis())) < 2e-7
        oracle = np.array([_contour_oracle(float(np.mean(psi[..., k])), psi_sign, f_sign)
                           for k in range(12)])
        np.testing.assert_allclose(dataset.variables["B_theta"][:], oracle[:, 1]*1e6,
                                   rtol=3e-6)
        np.testing.assert_allclose(dataset.variables["B_phi"][:], oracle[:, 2]*1e6,
                                   rtol=3e-7)
        pr, pz = _grad(r, z)
        f_squared = (R0*B0)**2-4*C*(psi-PSI_EDGE)
        np.testing.assert_allclose(bmod*1e-4, np.sqrt(pr*pr+pz*pz+f_squared)/r,
                                   rtol=3e-7)
        expected_flux = quad(lambda p: _contour_oracle(p)[0], _psi(*_axis()), PSI_EDGE)[0]
        assert torflux == pytest.approx(-f_sign*expected_flux*1e8, rel=3e-6)
        assert float(dataset.torflux) == torflux


@pytest.mark.parametrize("mode", ["uninitialized", "missing", "malformed", "quality", "valid"])
def test_native_init_and_magdata_are_recoverable(tmp_path, boozer_api, mode):
    """A native status must become a Python exception, never terminate Python."""
    if mode != "uninitialized":
        if mode == "malformed":
            (tmp_path/"field.g").write_text(
                f"{'truncated synthetic':48s}{0:4d}{65:4d}{65:4d}\n1.0\n")
        elif mode != "missing":
            _write_solovev(tmp_path/"field.g", pprime_scale=.7 if mode == "quality" else 1)
        if mode == "quality":
            _write_solovev(tmp_path/"valid.g")
        (tmp_path/"efit_to_boozer.inp").write_text("100\n32\n64\n64\n32\n12\n1e10\n")
        (tmp_path/"field_divB0.inp").write_text(
            "0\n1\n1.0\n72\n0.99\n4\n'field.g'\n'unused'\n'unused'\n'unused'\n0\n0\n1\n")
    script = """
import importlib.util, json, numpy as np, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location('_efit_to_boozer', sys.argv[1])
ext = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ext)
api, mode = ext.efit_to_boozer, sys.argv[2]
if mode == 'quality':
    deck = Path('field_divB0.inp')
    text = deck.read_text()
    deck.write_text(text.replace("'field.g'", "'valid.g'"))
    assert api.init() is None
    api.magdata(1, np.array(.35), np.array(0.), np.array(.7))
    deck.write_text(text)
try:
    if mode != 'uninitialized':
        assert api.init() is None
    values = api.magdata(1, np.array(.35), np.array(0.), np.array(.7))
except ValueError as error:
    assert mode != 'valid', str(error)
    if mode == 'quality':
        assert 'Grad-Shafranov' in str(error), str(error)
    try:
        api.magdata(1, np.array(.35), np.array(0.), np.array(.7))
    except ValueError:
        pass
    else:
        raise AssertionError('Rejected init left stale magdata state accessible')
    print('RESULT:'+json.dumps({'error': str(error)}))
else:
    assert mode == 'valid', 'Native failure was silently accepted'
    print('RESULT:'+json.dumps({'values': [float(x) for x in values]}))
"""
    result = subprocess.run(
        [sys.executable, "-c", script, boozer_api._efit_to_boozer.__file__, mode],
        cwd=tmp_path, capture_output=True, text=True, timeout=60,
        env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    assert result.returncode == 0, result.stdout+result.stderr
    records = [line[7:] for line in result.stdout.splitlines() if line.startswith("RESULT:")]
    assert len(records) == 1, result.stdout+result.stderr
    payload = json.loads(records[0])
    if mode == "valid":
        values = payload["values"]
        q, _, current, _, _, bmod, _, r, _, _, z, *_ = values
        r, z = r*.01, z*.01
        q_expected, current_i, current_g = _contour_oracle(_psi(r, z))
        assert q == pytest.approx(q_expected, rel=3e-6)
        assert current == pytest.approx((current_g+current_i/q_expected)*1e6, rel=3e-6)
        pr, pz = _grad(r, z)
        assert bmod*1e-4 == pytest.approx(np.sqrt(pr*pr+pz*pz+current_g**2)/r, rel=3e-6)


@pytest.mark.parametrize("mode", ["magdata_s", "magdata_theta", "contour_s",
                                    "contour_edge", "deck_psimax",
                                    "init_exp_bpol", "contour_exp_bpol"])
def test_raw_native_input_rejection_survives_fast_math(tmp_path, boozer_api, mode):
    """Exercise compiled guards directly, without the Python validation layer."""
    _write_solovev(tmp_path/"field.g")
    (tmp_path/"efit_to_boozer.inp").write_text("100\n32\n64\n64\n32\n12\n1e10\n")
    (tmp_path/"field_divB0.inp").write_text(
        "0\n1\n1.0\n72\n0.99\n4\n'field.g'\n'unused'\n'unused'\n'unused'\n0\n0\n1\n")
    script = """
import importlib.util, numpy as np, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location('_efit_to_boozer', sys.argv[1])
ext = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ext)
api, mode = ext.efit_to_boozer, sys.argv[2]
deck = Path('efit_to_boozer.inp')
valid_deck = deck.read_text()
assert api.init() is None
bad_values = (-2, 2) if 'exp_bpol' in mode else (np.nan, np.inf, -np.inf)
for bad in bad_values:
    try:
        if mode == 'magdata_s':
            api.magdata(1, np.array(bad), np.array(0.), np.array(.7))
        elif mode == 'magdata_theta':
            api.magdata(1, np.array(.35), np.array(0.), np.array(bad))
        elif mode == 'deck_psimax':
            records = valid_deck.splitlines()
            records[6] = str(bad)
            deck.write_text('\\n'.join(records)+'\\n')
            api.init()
        elif mode == 'init_exp_bpol':
            api.init(exp_bpol=bad)
        else:
            s = np.array([bad if mode == 'contour_s' else .35])
            use_edge = mode == 'contour_edge'
            edge = bad if use_edge else 0.
            kwargs = {'exp_bpol': bad} if mode == 'contour_exp_bpol' else {}
            api.contour_transform('field.g', s, 8, 64, 16, 1, use_edge, edge,
                                  False, '', **kwargs)
    except ValueError:
        pass
    else:
        raise AssertionError(f'{mode} accepted {bad!r}')
    if mode in ('deck_psimax', 'init_exp_bpol'):
        deck.write_text(valid_deck)
        assert api.init() is None
    values = api.magdata(1, np.array(.35), np.array(0.), np.array(.7))
    assert len(values) == 16 and np.isfinite(values).all()
    if mode.startswith('contour_'):
        result = api.contour_transform('field.g', np.array([.35]), 8, 64, 16,
                                        1, False, 0., False, '')
        assert result[3] == 0 and np.isfinite(result[0]).all()
print('RAW_INPUT_REJECTION_OK')
"""
    result = subprocess.run(
        [sys.executable, "-c", script, boozer_api._efit_to_boozer.__file__, mode],
        cwd=tmp_path, capture_output=True, text=True, timeout=60,
        env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})
    assert result.returncode == 0, result.stdout+result.stderr
    assert "RAW_INPUT_REJECTION_OK" in result.stdout


@pytest.mark.parametrize("backend", ["native", "cli"])
def test_compiled_flux_bound_sentinels(tmp_path, boozer_api, backend):
    _write_solovev(tmp_path/"field.g")
    (tmp_path/"field_divB0.inp").write_text(
        "0\n1\n1.0\n72\n0.99\n4\n'field.g'\n'unused'\n'unused'\n'unused'\n0\n0\n1\n")
    records = ("1e10", "1e300", "2e10", "NaN", "Infinity", "-Infinity")
    environment = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    if backend == "cli":
        binary = os.environ.get("EFIT_TO_BOOZER_BINARY")
        if binary is None:
            binary = Path(boozer_api._efit_to_boozer.__file__).resolve().parent/"efit_to_boozer.x"
        binary = Path(binary).resolve()
        if not binary.is_file():
            pytest.skip("set EFIT_TO_BOOZER_BINARY to exercise the shared CLI loader")
        for record in records:
            (tmp_path/"efit_to_boozer.inp").write_text(f"100\n8\n32\n16\n8\n6\n{record}\n")
            result = subprocess.run([str(binary)], cwd=tmp_path, capture_output=True,
                                    text=True, timeout=60, env=environment)
            output = result.stdout+result.stderr
            if record in ("1e10", "1e300"):
                assert result.returncode == 0, output
            else:
                assert result.returncode > 0, output
                assert ("psi_edge" if record == "2e10" else "finite") in output, output
        return
    script = """
import importlib.util, numpy as np, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location('_efit_to_boozer', sys.argv[1])
ext = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ext)
api, deck = ext.efit_to_boozer, Path('efit_to_boozer.inp')
for record in ('1e10', '1e300', '2e10', 'NaN', 'Infinity', '-Infinity'):
    deck.write_text('100\\n8\\n32\\n16\\n8\\n6\\n'+record+'\\n')
    try:
        assert api.init() is None
    except ValueError:
        assert record not in ('1e10', '1e300'), record
        try:
            api.magdata(1, np.array(.35), np.array(0.), np.array(.7))
        except ValueError:
            pass
        else:
            raise AssertionError('Rejected bound left stale initialized state')
    else:
        assert record in ('1e10', '1e300'), record
        assert np.isfinite(api.magdata(1, np.array(.35), np.array(0.), np.array(.7))).all()
deck.write_text('100\\n8\\n32\\n16\\n8\\n6\\n1e10\\n')
assert api.init() is None
assert np.isfinite(api.magdata(1, np.array(.35), np.array(0.), np.array(.7))).all()
print('BOUND_SENTINELS_OK')
"""
    result = subprocess.run(
        [sys.executable, "-c", script, boozer_api._efit_to_boozer.__file__],
        cwd=tmp_path, capture_output=True, text=True, timeout=60, env=environment)
    assert result.returncode == 0, result.stdout+result.stderr
    assert "BOUND_SENTINELS_OK" in result.stdout


def test_native_cli_historical_scan_sentinel_uses_file_edge(tmp_path):
    from libneo.boozer import BoozerFile
    from libneo.eqdsk_to_boozer_chartmap import _write_field_divB0_inp
    from scipy.integrate import quad

    binary = Path(os.environ.get("EFIT_TO_BOOZER_BINARY",
                                 Path(__file__).resolve().parents[2]/"build/efit_to_boozer.x"))
    if not binary.exists():
        pytest.skip("efit_to_boozer.x not built")
    gfile = tmp_path/"solovev.g"
    _write_solovev(gfile)
    # Exactly the historical seven-record deck, without optional selectors.
    (tmp_path/"efit_to_boozer.inp").write_text("100\n32\n64\n64\n8\n12\n1e10\n")
    _write_field_divB0_inp(tmp_path/"field_divB0.inp", str(gfile))
    result = subprocess.run([str(binary.resolve())], cwd=tmp_path,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout+result.stderr
    output = BoozerFile(str(tmp_path/"fromefit.bc"))
    expected_flux = 2*np.pi*quad(lambda p: _contour_oracle(p)[0],
                                _psi(*_axis()), PSI_EDGE)[0]
    assert output.flux == pytest.approx(expected_flux, rel=3e-6)


def test_explicit_units_reject_nonfinite_gfile_data(tmp_path, boozer_api):
    from libneo.eqdsk_base import read_eqdsk, write_eqdsk

    valid = tmp_path/"valid.g"
    _write_solovev(valid)
    for field, bad in (("PsiVs", np.nan), ("dpressdpsiprof", np.inf)):
        data = read_eqdsk(valid)
        data[field].flat[0] = bad
        write_eqdsk(tmp_path/f"{field}.g", data)
    script = """
import importlib.util, numpy as np, sys
spec = importlib.util.spec_from_file_location('_efit_to_boozer', sys.argv[1])
ext = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ext)
api = ext.efit_to_boozer
for path in ('PsiVs.g', 'dpressdpsiprof.g'):
    for units in (0, 1):
        try:
            api.contour_transform(path, np.array([.3]), 8, 64, 16, 1, False,
                0., True, 'synthetic finite-input regression', exp_bpol=units)
        except ValueError as error:
            assert 'nonfinite G-EQDSK' in str(error), str(error)
        else:
            raise AssertionError('Explicit units accepted nonfinite equilibrium data')
result = api.contour_transform('valid.g', np.array([.3]), 8, 64, 16,
                                1, False, 0., False, '', exp_bpol=0)
assert result[3] == 0 and np.isfinite(result[0]).all()
print('NONFINITE_GFILE_REJECTED')
"""
    result = subprocess.run(
        [sys.executable, "-c", script, boozer_api._efit_to_boozer.__file__],
        cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout+result.stderr
    assert "NONFINITE_GFILE_REJECTED" in result.stdout
