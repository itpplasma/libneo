# %% Standard imports
import os
from pathlib import Path

import pytest
import numpy as np
import matplotlib.pyplot as plt

# Import modules to test
from libneo import StorGeom2MarsCoords
from libneo import order_monotonically

mars_dir = "/proj/plasma/DATA/DEMO/MARS/MARSQ_KNTV10_NEO2profs_KEYTORQ_1/"
pytestmark = pytest.mark.skipif(not os.path.exists(mars_dir), reason="MARS data not available on this system")

sqrtspol = np.linspace(0, 1, 5)
chi = np.linspace(-np.pi, np.pi, 10)
chi[-1] = chi[-1] - 0.1
coordsystem = {'radius': sqrtspol.tolist(), 'angle': np.meshgrid(chi, sqrtspol)[0].tolist()}
new_coordsystem = {'radius': (0.5 * sqrtspol).tolist(), 'angle': np.mod(np.meshgrid(chi, sqrtspol)[0], 2 * np.pi).tolist()}


def _figure_path(name: str) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    target = repo_root / "build" / "test" / "python"
    target.mkdir(parents=True, exist_ok=True)
    return target / name


def _save_current_figure(name: str) -> None:
    plt.savefig(_figure_path(name), dpi=150, bbox_inches="tight")
    plt.close()


def test_order_monotonically():
    input_coords = {'radius': [1, 3, 2], 'angle': [[3, 1, 2], [7, 9, 8], [4, 5, 6]]}
    new_coords = {'radius': [1, 2, 3], 'angle': [[1, 2, 3], [4, 5, 6], [7, 8, 9]]}
    coordsystem_sorted, new_coordsystem_sorted = order_monotonically(input_coords, new_coords)
    assert is_monotonically_increasing(coordsystem_sorted['radius'])
    assert not is_monotonically_increasing(new_coordsystem_sorted['radius'])
    assert is_monotonically_increasing(coordsystem_sorted['angle'])
    assert not is_monotonically_increasing(new_coordsystem_sorted['angle'])
    assert np.concatenate(coordsystem_sorted['angle']).tolist() == [1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert np.concatenate(new_coordsystem_sorted['angle']).tolist() == [2, 3, 1, 7, 8, 9, 4, 6, 5]


def is_monotonically_increasing(arr):
    return np.all(np.diff(arr) > 0)


def test_MarsCoords2StorThetageom_init():
    mars = StorGeom2MarsCoords(mars_dir)
    assert mars is not None


def test_MarsCoords2StorThetageom_load_points_from():
    mars_coords, stor_geom_coords = StorGeom2MarsCoords.load_points_from(mars_dir)
    assert len(mars_coords['angle']) == len(mars_coords['radius'])
    assert len(mars_coords['radius']) == len(stor_geom_coords['radius'])
    for i in range(len(mars_coords['radius'])):
        assert isinstance(mars_coords['radius'][i], float)
        assert len(mars_coords['angle'][i]) == len(stor_geom_coords['angle'][i])


def test_MarsCoord2StorThetageom_stor2sqrtspol():
    converter = StorGeom2MarsCoords(mars_dir)
    stor2sqrtspol = converter._stor2sqrtspol
    assert stor2sqrtspol(0) == 0
    assert stor2sqrtspol(1) == 1
    assert stor2sqrtspol(0.5) != 0.5
    from neo2_mars import mars_sqrtspol2stor
    sqrtspol_values = np.linspace(0, 1, 100)
    stor = mars_sqrtspol2stor(mars_dir, sqrtspol_values)
    assert np.allclose(sqrtspol_values, stor2sqrtspol(stor), atol=1e-2)


def test_MarsCoords2StorThetageom_visual_check():
    converter = StorGeom2MarsCoords(mars_dir)
    test_stor = np.linspace(0.0, 1, 5)
    test_theta_geom = np.linspace(-np.pi, np.pi, 800)
    plt.figure()
    for stor in test_stor:
        sqrtspol_val, chi_val = converter(stor, test_theta_geom)
        plt.plot(test_theta_geom, chi_val, '.', label=f"stor={stor}")
    plt.xlabel('theta_geom [1]')
    plt.ylabel('chi [1]')
    plt.legend()
    _save_current_figure('coordinate_converter_chi_vs_theta_geom.png')


def test_MarsCoords2StorThetageom_coords_domain_visual_check():
    converter = StorGeom2MarsCoords(mars_dir)
    mars_coords = converter.mars_coords
    stor_geom_coords = converter.stor_geom_coords
    fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    colors = ['b', 'y', 'g', 'r']

    stor_geom_sorted, mars_sorted = order_monotonically(stor_geom_coords, mars_coords)
    for idx, i in enumerate(range(0, len(mars_sorted['radius']), 100)):
        ax[0].plot(stor_geom_sorted['angle'][i], mars_sorted['angle'][i], '.' + colors[idx % len(colors)],
                   label=f"sqrtspol={mars_sorted['radius'][i]:.2f}")
        discontinuity = np.where(np.abs(np.diff(mars_sorted['angle'][i])) > np.pi)[0]
        if len(discontinuity) > 0:
            ax[0].axvline(x=stor_geom_sorted['angle'][i][discontinuity[0]], color=colors[idx % len(colors)], linestyle='--',
                          label='discontinuity by data')
            ax[0].axvline(x=stor_geom_sorted['angle'][i][discontinuity[0] + 1], color=colors[idx % len(colors)], linestyle='--')
    ax[0].set_xlabel('theta_geom [1]')
    ax[0].set_ylabel('chi [1]')
    ax[0].legend()

    mars_sorted, stor_geom_sorted = order_monotonically(mars_coords, stor_geom_coords)
    for idx, i in enumerate(range(0, len(mars_sorted['radius']), 100)):
        ax[1].plot(mars_sorted['angle'][i], stor_geom_sorted['angle'][i], '.' + colors[idx % len(colors)],
                   label=f"stor={stor_geom_sorted['radius'][i]:.2f}")
        discontinuity = np.where(np.abs(np.diff(stor_geom_sorted['angle'][i])) > np.pi)[0]
        if len(discontinuity) > 0:
            ax[1].axvline(x=mars_sorted['angle'][i][discontinuity[0]], color=colors[idx % len(colors)], linestyle='--',
                          label='discontinuity by data')
            ax[1].axvline(x=mars_sorted['angle'][i][discontinuity[0] + 1], color=colors[idx % len(colors)], linestyle='--')
    ax[1].set_xlabel('chi [1]')
    ax[1].set_ylabel('theta_geom [1]')
    ax[1].legend()
    _save_current_figure('coordinate_converter_coordinate_domains.png')


def test_MarsCoord2StorThetageom_performance_profile():
    import cProfile
    import pstats
    cProfile.run('test_MarsCoords2StorThetageom_performance()', filename='mars_coords2stor_thetageom.profile')
    p = pstats.Stats('mars_coords2stor_thetageom.profile')
    p.strip_dirs().sort_stats('cumulative').print_stats(20)


def test_MarsCoords2StorThetageom_performance():
    converter = StorGeom2MarsCoords(mars_dir)
    test_stor = np.linspace(0.0, 1, 100)
    test_theta_geom = np.linspace(-np.pi, np.pi, 100)
    for _ in range(20):
        for stor in test_stor:
            converter(stor, test_theta_geom)


if __name__ == "__main__":
    test_MarsCoords2StorThetageom_init()
    test_MarsCoords2StorThetageom_load_points_from()
    test_MarsCoord2StorThetageom_stor2sqrtspol()
    print("All tests passed!")
    test_MarsCoords2StorThetageom_visual_check()
    test_MarsCoords2StorThetageom_coords_domain_visual_check()
