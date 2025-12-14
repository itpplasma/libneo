from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from libneo.vmec import VMECGeometry


@dataclass(frozen=True)
class PortSpec:
    phi: float
    z: float
    radius: float
    length: float
    inboard: bool = False


def wall_mesh_from_wout(
    wout_path: Path,
    *,
    s_index: int = -1,
    n_theta: int = 128,
    n_zeta: int = 96,
    use_asym: bool = True,
) -> "trimesh.Trimesh":
    import trimesh

    geom = VMECGeometry.from_file(str(wout_path))
    ns = geom.rmnc.shape[1]
    if s_index < 0:
        s_index = ns + s_index
    if s_index < 0 or s_index >= ns:
        raise ValueError("s_index out of range")

    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, n_zeta, endpoint=False)

    xyz = np.zeros((n_zeta, n_theta, 3), dtype=float)
    for i, zet in enumerate(zeta):
        r, z, _ = geom.coords(s_index, theta, float(zet), use_asym=use_asym)
        xyz[i, :, 0] = r * np.cos(zet)
        xyz[i, :, 1] = r * np.sin(zet)
        xyz[i, :, 2] = z

    vertices = xyz.reshape((-1, 3))

    faces: list[tuple[int, int, int]] = []
    for i in range(n_zeta):
        inext = (i + 1) % n_zeta
        for j in range(n_theta):
            jnext = (j + 1) % n_theta
            v00 = i * n_theta + j
            v10 = inext * n_theta + j
            v01 = i * n_theta + jnext
            v11 = inext * n_theta + jnext
            faces.append((v00, v10, v11))
            faces.append((v00, v11, v01))

    return trimesh.Trimesh(vertices=vertices, faces=np.asarray(faces), process=False)


def _wrap_angle_pi(angle: np.ndarray) -> np.ndarray:
    return (angle + np.pi) % (2.0 * np.pi) - np.pi


def _find_surface_point(
    mesh: "trimesh.Trimesh",
    phi_target: float,
    z_target: float,
    *,
    phi_tol: float = 0.15,
    z_tol: float = 0.3,
    inboard: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    centers = mesh.triangles_center
    phi = np.arctan2(centers[:, 1], centers[:, 0])
    z = centers[:, 2]
    r = np.sqrt(centers[:, 0] ** 2 + centers[:, 1] ** 2)

    dphi = np.abs(_wrap_angle_pi(phi - float(phi_target)))
    dz = np.abs(z - float(z_target))
    mask = (dphi <= phi_tol) & (dz <= z_tol)
    if not np.any(mask):
        raise RuntimeError("no candidate faces for requested (phi,z)")

    if inboard:
        idx = int(np.argmin(np.where(mask, r, np.inf)))
        sign = -1.0
    else:
        idx = int(np.argmax(np.where(mask, r, -np.inf)))
        sign = 1.0

    center = centers[idx].copy()
    direction = sign * np.array([np.cos(phi_target), np.sin(phi_target), 0.0], dtype=float)
    nrm = float(np.linalg.norm(direction))
    if nrm == 0.0:
        raise RuntimeError("degenerate port direction")
    direction /= nrm
    return center, direction


def cut_cylindrical_ports(mesh: "trimesh.Trimesh", ports: list[PortSpec]) -> "trimesh.Trimesh":
    import trimesh

    centers = mesh.triangles_center
    keep = np.ones((mesh.faces.shape[0],), dtype=bool)

    for p in ports:
        center, direction = _find_surface_point(
            mesh, p.phi, p.z, inboard=p.inboard
        )
        v = centers - center[None, :]
        t = v @ direction
        radial = np.linalg.norm(v - t[:, None] * direction[None, :], axis=1)
        half_len = 0.5 * float(p.length)
        in_cyl = (np.abs(t) <= half_len) & (radial <= float(p.radius))
        keep &= ~in_cyl

    faces = mesh.faces[keep]
    return trimesh.Trimesh(vertices=mesh.vertices.copy(), faces=faces, process=False)
