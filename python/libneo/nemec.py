from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


TWOPI = 2.0 * np.pi
SOURCE_NOTES = (
    "NEMEC manual and readnemec.f90 in the Strumberger output bundle: "
    "R/Z and all field components use phase m*theta - n*zeta*nfp, with "
    "left-handed flux coordinates and hphip=(1/2*pi)*dPhi/ds. "
    "VMEC/STELLOPT wout stores xn=n*nfp and uses cos/sin(m*theta - xn*zeta); "
    "DESC's VMEC reader/writer uses the same phase and signgs=-1. "
    "SFINCS read_NEMEC_file.f90 flips n and several signs for its own "
    "downstream convention; those flips are not applied for VMEC wout output."
)


@dataclass(frozen=True)
class NemecData:
    nfp: int
    ns: int
    mpol: int
    ntor: int
    iasym: int
    phiedge: float
    xm: np.ndarray
    xn: np.ndarray
    rmnc: np.ndarray
    zmns: np.ndarray
    rmns: np.ndarray
    zmnc: np.ndarray
    lmns: np.ndarray
    lmnc: np.ndarray
    bsupumnc: np.ndarray
    bsupvmnc: np.ndarray
    bsupumns: np.ndarray
    bsupvmns: np.ndarray
    bsubumnc: np.ndarray
    bsubvmnc: np.ndarray
    bsubsmns: np.ndarray
    bsubumns: np.ndarray
    bsubvmns: np.ndarray
    bsubsmnc: np.ndarray
    iotaf: np.ndarray
    iotas: np.ndarray
    mass: np.ndarray
    pres: np.ndarray
    phipf: np.ndarray
    phips: np.ndarray
    buco: np.ndarray
    bvco: np.ndarray
    phi: np.ndarray
    vp: np.ndarray
    over_r: np.ndarray
    jcuru: np.ndarray
    jcurv: np.ndarray
    specw: np.ndarray


def read_nemec(path: str | Path) -> NemecData:
    values = np.fromstring(Path(path).read_text(), sep=" ")
    if values.size < 8:
        raise ValueError("NEMEC file is too short")

    _, nfp_f, nrho_f, mpol_f, ntor_f, mpnt_f, iasym_f, phiedge = values[:8]
    nfp = int(round(nfp_f))
    nrho = int(round(nrho_f))
    mpol = int(round(mpol_f))
    ntor = int(round(ntor_f))
    mpnt = int(round(mpnt_f))
    iasym = int(round(iasym_f))
    modes = _mode_numbers(nfp, mpol, ntor)
    if mpnt != modes.shape[0]:
        raise ValueError(f"NEMEC mpnt={mpnt} does not match mpol/ntor modes={modes.shape[0]}")

    fourier_size = nrho * mpnt * 16
    profile_size = (nrho - 1) * 12
    expected = 8 + fourier_size + profile_size
    if values.size < expected:
        raise ValueError(f"NEMEC file has {values.size} values, expected at least {expected}")

    fourier = values[8:8 + fourier_size].reshape(nrho, mpnt, 16)
    profile = values[8 + fourier_size:expected].reshape(nrho - 1, 12)

    iota = _with_axis(profile[:, 0])
    hphip = _with_axis(profile[:, 3], first_value=profile[0, 3])
    return NemecData(
        nfp=nfp,
        ns=nrho,
        mpol=mpol,
        ntor=ntor,
        iasym=iasym,
        phiedge=float(phiedge),
        xm=modes[:, 0].astype(float),
        xn=modes[:, 1].astype(float),
        rmnc=fourier[:, :, 0],
        zmns=fourier[:, :, 1],
        rmns=fourier[:, :, 2],
        zmnc=fourier[:, :, 3],
        bsupumnc=fourier[:, :, 4],
        bsupvmnc=fourier[:, :, 5],
        bsupumns=fourier[:, :, 6],
        bsupvmns=fourier[:, :, 7],
        lmns=fourier[:, :, 8],
        lmnc=fourier[:, :, 9],
        bsubumnc=fourier[:, :, 10],
        bsubvmnc=fourier[:, :, 11],
        bsubsmns=fourier[:, :, 12],
        bsubumns=fourier[:, :, 13],
        bsubvmns=fourier[:, :, 14],
        bsubsmnc=fourier[:, :, 15],
        iotaf=iota,
        iotas=iota,
        mass=_with_axis(profile[:, 1]),
        pres=_with_axis(profile[:, 2]),
        phipf=-TWOPI * hphip,
        phips=-TWOPI * hphip,
        buco=_with_axis(profile[:, 4]),
        bvco=_with_axis(profile[:, 5]),
        phi=_with_axis(profile[:, 6]),
        vp=_with_axis(profile[:, 7]),
        over_r=_with_axis(profile[:, 8]),
        jcuru=_with_axis(profile[:, 9]),
        jcurv=_with_axis(profile[:, 10]),
        specw=_with_axis(profile[:, 11]),
    )


def write_vmec_wout(data: NemecData, path: str | Path) -> None:
    out = Path(path)
    with Dataset(out, "w") as ds:
        ds.source_conventions = SOURCE_NOTES
        ds.createDimension("dim_00100", 100)
        ds.createDimension("dim_00200", 200)
        ds.createDimension("dim_00020", 20)
        ds.createDimension("dim_00001", 1)
        ds.createDimension("mn_mode", data.xm.size)
        ds.createDimension("mn_mode_nyq", data.xm.size)
        ds.createDimension("n_tor", data.ntor + 1)
        ds.createDimension("radius", data.ns)
        ds.createDimension("time", 100)

        _scalar(ds, "version_", 9.0)
        _char(ds, "input_extension", "nemec2vmec", "dim_00100")
        _char(ds, "mgrid_file", "", "dim_00200")
        _char(ds, "pcurr_type", "", "dim_00020")
        _char(ds, "pmass_type", "", "dim_00020")
        _char(ds, "piota_type", "", "dim_00020")
        _char(ds, "mgrid_mode", "", "dim_00001")

        for name, value in _scalar_values(data).items():
            _scalar(ds, name, value)

        _array(ds, "xm", ("mn_mode",), data.xm)
        _array(ds, "xn", ("mn_mode",), data.xn)
        _array(ds, "xm_nyq", ("mn_mode_nyq",), data.xm)
        _array(ds, "xn_nyq", ("mn_mode_nyq",), data.xn)

        axis_modes = data.ntor + 1
        _array(ds, "raxis_cc", ("n_tor",), data.rmnc[0, :axis_modes])
        _array(ds, "zaxis_cs", ("n_tor",), data.zmns[0, :axis_modes])
        _array(ds, "raxis_cs", ("n_tor",), data.rmns[0, :axis_modes])
        _array(ds, "zaxis_cc", ("n_tor",), data.zmnc[0, :axis_modes])

        for name, values in _profile_arrays(data).items():
            _array(ds, name, ("radius",), values)
        _array(ds, "fsqt", ("time",), np.zeros(100, dtype=float))
        _array(ds, "wdot", ("time",), np.zeros(100, dtype=float))

        _array(ds, "rmnc", ("radius", "mn_mode"), data.rmnc)
        _array(ds, "zmns", ("radius", "mn_mode"), data.zmns)
        _array(ds, "lmns", ("radius", "mn_mode"), data.lmns)
        _array(ds, "rmns", ("radius", "mn_mode"), data.rmns)
        _array(ds, "zmnc", ("radius", "mn_mode"), data.zmnc)
        _array(ds, "lmnc", ("radius", "mn_mode"), data.lmnc)

        bmnc, bmns, gmnc, gmns = _magnetic_spectra(data)
        _array(ds, "gmnc", ("radius", "mn_mode_nyq"), gmnc)
        _array(ds, "gmns", ("radius", "mn_mode_nyq"), gmns)
        _array(ds, "bmnc", ("radius", "mn_mode_nyq"), bmnc)
        _array(ds, "bmns", ("radius", "mn_mode_nyq"), bmns)
        _array(ds, "bsubumnc", ("radius", "mn_mode_nyq"), data.bsubumnc)
        _array(ds, "bsubvmnc", ("radius", "mn_mode_nyq"), data.bsubvmnc)
        _array(ds, "bsubsmns", ("radius", "mn_mode_nyq"), data.bsubsmns)
        _array(ds, "bsubumns", ("radius", "mn_mode_nyq"), data.bsubumns)
        _array(ds, "bsubvmns", ("radius", "mn_mode_nyq"), data.bsubvmns)
        _array(ds, "bsubsmnc", ("radius", "mn_mode_nyq"), data.bsubsmnc)
        _array(ds, "bsupumnc", ("radius", "mn_mode_nyq"), data.bsupumnc)
        _array(ds, "bsupvmnc", ("radius", "mn_mode_nyq"), data.bsupvmnc)
        _array(ds, "bsupumns", ("radius", "mn_mode_nyq"), data.bsupumns)
        _array(ds, "bsupvmns", ("radius", "mn_mode_nyq"), data.bsupvmns)


def _mode_numbers(nfp: int, mpol: int, ntor: int) -> np.ndarray:
    modes = []
    for m in range(mpol):
        n_min = 0 if m == 0 else -ntor
        for n in range(n_min, ntor + 1):
            modes.append((m, n * nfp))
    return np.asarray(modes, dtype=int)


def _magnetic_spectra(data: NemecData) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cos_phase, sin_phase = _projection_grid(data)
    bmnc = np.zeros_like(data.rmnc)
    bmns = np.zeros_like(data.rmnc)
    gmnc = np.zeros_like(data.rmnc)
    gmns = np.zeros_like(data.rmnc)

    for js in range(data.ns):
        bsupu = _series(data.bsupumnc[js], data.bsupumns[js], cos_phase, sin_phase)
        bsupv = _series(data.bsupvmnc[js], data.bsupvmns[js], cos_phase, sin_phase)
        bsubu = _series(data.bsubumnc[js], data.bsubumns[js], cos_phase, sin_phase)
        bsubv = _series(data.bsubvmnc[js], data.bsubvmns[js], cos_phase, sin_phase)
        bdotb = np.maximum(bsupu * bsubu + bsupv * bsubv, 0.0)
        bmnc[js], bmns[js] = _project(np.sqrt(bdotb), cos_phase, sin_phase, data.xm, data.xn)

        lamu = _lambda_theta(data, js, cos_phase, sin_phase)
        numerator = (1.0 + lamu) * data.phipf[js]
        sqrtg = np.zeros_like(bsupv)
        np.divide(numerator, bsupv, out=sqrtg, where=np.abs(bsupv) > 1.0e-14)
        gmnc[js], gmns[js] = _project(sqrtg, cos_phase, sin_phase, data.xm, data.xn)

    return bmnc, bmns, gmnc, gmns


def _projection_grid(data: NemecData) -> tuple[np.ndarray, np.ndarray]:
    ntheta = max(64, 4 * data.mpol + 1)
    nzeta = max(64, 4 * data.ntor + 1)
    theta = np.linspace(0.0, TWOPI, ntheta, endpoint=False)
    zeta = np.linspace(0.0, TWOPI / max(data.nfp, 1), nzeta, endpoint=False)
    angle = (
        data.xm[:, None, None] * theta[None, :, None]
        - data.xn[:, None, None] * zeta[None, None, :]
    )
    return np.cos(angle), np.sin(angle)


def _series(cos_coeff: np.ndarray, sin_coeff: np.ndarray, cos_phase: np.ndarray, sin_phase: np.ndarray) -> np.ndarray:
    return np.tensordot(cos_coeff, cos_phase, axes=(0, 0)) + np.tensordot(sin_coeff, sin_phase, axes=(0, 0))


def _lambda_theta(data: NemecData, js: int, cos_phase: np.ndarray, sin_phase: np.ndarray) -> np.ndarray:
    return (
        np.tensordot(data.xm * data.lmns[js], cos_phase, axes=(0, 0))
        - np.tensordot(data.xm * data.lmnc[js], sin_phase, axes=(0, 0))
    )


def _project(
    values: np.ndarray,
    cos_phase: np.ndarray,
    sin_phase: np.ndarray,
    xm: np.ndarray,
    xn: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    cos_coeff = 2.0 * np.mean(values[None, :, :] * cos_phase, axis=(1, 2))
    sin_coeff = 2.0 * np.mean(values[None, :, :] * sin_phase, axis=(1, 2))
    axis = (xm == 0.0) & (xn == 0.0)
    cos_coeff[axis] *= 0.5
    sin_coeff[axis] = 0.0
    return cos_coeff, sin_coeff


def _with_axis(values: np.ndarray, first_value: float = 0.0) -> np.ndarray:
    return np.concatenate(([float(first_value)], np.asarray(values, dtype=float)))


def _scalar(ds: Dataset, name: str, value: float | int) -> None:
    dtype = "i4" if isinstance(value, int) else "f8"
    var = ds.createVariable(name, dtype)
    var[...] = value


def _char(ds: Dataset, name: str, value: str, dim: str) -> None:
    var = ds.createVariable(name, "S1", (dim,))
    encoded = value.encode("ascii")[: ds.dimensions[dim].size]
    arr = np.full(ds.dimensions[dim].size, b" ", dtype="S1")
    arr[:len(encoded)] = np.frombuffer(encoded, dtype="S1")
    var[:] = arr


def _array(ds: Dataset, name: str, dims: tuple[str, ...], values: np.ndarray) -> None:
    var = ds.createVariable(name, "f8", dims)
    var[:] = np.asarray(values, dtype=float)


def _scalar_values(data: NemecData) -> dict[str, float | int]:
    r_boundary, z_boundary = _boundary_extrema_samples(data)
    return {
        "wb": 0.0,
        "wp": 0.0,
        "gamma": 0.0,
        "rmax_surf": float(np.max(r_boundary)),
        "rmin_surf": float(np.min(r_boundary)),
        "zmax_surf": float(np.max(z_boundary)),
        "nfp": data.nfp,
        "ns": data.ns,
        "mpol": data.mpol,
        "ntor": data.ntor,
        "mnmax": data.xm.size,
        "mnmax_nyq": data.xm.size,
        "niter": 0,
        "itfsq": 0,
        "lasym__logical__": 1 if data.iasym else 0,
        "lrecon__logical__": 0,
        "lfreeb__logical__": 1,
        "lrfp__logical__": 0,
        "ier_flag": 0,
        "aspect": 0.0,
        "betatotal": 0.0,
        "betapol": 0.0,
        "betator": 0.0,
        "betaxis": 0.0,
        "b0": 0.0,
        "rbtor0": 0.0,
        "rbtor": 0.0,
        "signgs": -1,
        "IonLarmor": 0.0,
        "volavgB": 0.0,
        "ctor": 0.0,
        "Aminor_p": 0.0,
        "Rmajor_p": 0.0,
        "volume_p": 0.0,
        "ftolv": 0.0,
        "fsql": 0.0,
        "fsqr": 0.0,
        "fsqz": 0.0,
        "nextcur": 0,
    }


def _boundary_extrema_samples(data: NemecData) -> tuple[np.ndarray, np.ndarray]:
    theta = np.linspace(0.0, TWOPI, 257, endpoint=False)
    zeta = np.linspace(0.0, TWOPI / max(data.nfp, 1), 257, endpoint=False)
    angle = (
        data.xm[:, None, None] * theta[None, :, None]
        - data.xn[:, None, None] * zeta[None, None, :]
    )
    r = np.tensordot(data.rmnc[-1], np.cos(angle), axes=(0, 0))
    z = np.tensordot(data.zmns[-1], np.sin(angle), axes=(0, 0))
    if data.iasym:
        r += np.tensordot(data.rmns[-1], np.sin(angle), axes=(0, 0))
        z += np.tensordot(data.zmnc[-1], np.cos(angle), axes=(0, 0))
    return r, z


def _profile_arrays(data: NemecData) -> dict[str, np.ndarray]:
    zeros = np.zeros(data.ns, dtype=float)
    q_factor = np.zeros(data.ns, dtype=float)
    np.divide(1.0, data.iotaf, out=q_factor, where=data.iotaf != 0.0)
    return {
        "iotaf": data.iotaf,
        "q_factor": q_factor,
        "presf": data.pres,
        "phi": data.phi,
        "phipf": data.phipf,
        "chi": zeros,
        "chipf": data.iotaf * data.phipf,
        "jcuru": data.jcuru,
        "jcurv": data.jcurv,
        "iotas": data.iotas,
        "mass": data.mass,
        "pres": data.pres,
        "beta_vol": zeros,
        "buco": data.buco,
        "bvco": data.bvco,
        "vp": data.vp,
        "specw": data.specw,
        "phips": data.phips,
        "over_r": data.over_r,
        "jdotb": zeros,
        "bdotb": zeros,
        "bdotgradv": zeros,
        "DMerc": zeros,
        "DShear": zeros,
        "DWell": zeros,
        "DCurr": zeros,
        "DGeod": zeros,
        "equif": zeros,
    }
