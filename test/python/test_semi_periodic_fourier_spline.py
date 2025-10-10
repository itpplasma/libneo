# %% Standard imports
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pytest

# modules to test
from libneo import SemiPeriodicFourierSpline
from libneo import SemiPeriodicFourierSpline_ExceedingNyquistLimit


def _figure_path(name: str) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    target = repo_root / "build" / "test" / "python"
    target.mkdir(parents=True, exist_ok=True)
    return target / name


def _save_current_figure(name: str) -> None:
    plt.savefig(_figure_path(name), dpi=150, bbox_inches="tight")
    plt.close()


s = np.linspace(0, 1, 5)
angle = np.linspace(-np.pi, np.pi, 10)


def test_init():
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F)


def test_SemiPeriodicFourierSpline_ExceedingNyquistLimit():
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    with pytest.raises(SemiPeriodicFourierSpline_ExceedingNyquistLimit):
        SemiPeriodicFourierSpline(s, angle, F, max_modenumber=len(angle))


def test_fourier_coefs_shape():
    max_m = 4
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F, max_modenumber=max_m)
    coefs = spline._fourier_coefs
    assert coefs.shape == (len(s), max_m + 1)


def test_get_fourier_coefs_splines():
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F)
    assert np.allclose(spline._fourier_coefs_splines(s), spline._fourier_coefs)
    assert spline._fourier_coefs_splines(s).shape == (len(s), len(spline._modenumbers))


def test_call():
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F)
    assert np.allclose(spline(s, angle), F)
    assert np.allclose(spline(s, angle + 2 * np.pi), spline(s, angle))
    assert np.allclose(spline(s, angle - 2 * np.pi), spline(s, angle))
    r_new = np.random.uniform(0, 1, 5)
    angle_new = np.random.uniform(-np.pi, np.pi, 5)
    assert np.allclose(spline(r_new, angle_new, grid=False), trial_trigonometric_func(r_new, angle_new))


def test_get_fourier_coefs_visual_check():
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F)
    plt.figure()
    r_plot = np.linspace(0, 1, 30)
    coefs_plot = spline._fourier_coefs_splines(r_plot)
    for i, m in enumerate(spline._modenumbers):
        plt.plot(r_plot, coefs_plot[:, i].real, '-', label=f'real m={m}')
        plt.plot(r_plot, coefs_plot[:, i].imag, '.--', label=f'imag m={m}')
    plt.xlabel('s [1]')
    plt.ylabel('fm [1]')
    plt.title('Fourier coefficients of trigonometric trial function depending on s')
    plt.legend()
    _save_current_figure('semi_periodic_fourier_coefs.png')


def test_trial_theta_geom_func_visual_check():
    local_angle = np.linspace(-np.pi, np.pi, 100)
    Angle, S = np.meshgrid(local_angle, s)
    Theta_geom = theta_geom_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, local_angle, Theta_geom)
    plt.figure()
    angle_plot = np.linspace(-np.pi, np.pi, 19)
    s_plot = np.linspace(0, 1, 9)
    Angle_plot, S_plot = np.meshgrid(angle_plot, s_plot)
    R, Z = polar2cart(S_plot, theta_geom_func(S_plot, Angle_plot))
    plt.plot(R.flatten(), Z.flatten(), 'o', label='via base func')
    R, Z = polar2cart(S_plot, spline(s_plot, angle_plot))
    plt.plot(R.flatten(), Z.flatten(), 'x', label='via spline')
    plt.xlabel('R [1]')
    plt.ylabel('Z [1]')
    plt.legend()
    plt.title('Theta_geom rays (uniform grid)')
    _save_current_figure('semi_periodic_theta_geom_uniform.png')


def test_trial_theta_geom_func_different_angle_per_surface_visual_check():
    local_angle = np.linspace(-np.pi, np.pi, 100)
    Angle, S = np.meshgrid(local_angle, s)
    Angle[:, 1:-1] = Angle[:, 1:-1] * (S[:, 1:-1] + 1) / 2
    Theta_geom = theta_geom_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, Angle, Theta_geom)
    plt.figure()
    angle_plot = np.linspace(-np.pi, np.pi, 19)
    s_plot = np.linspace(0, 1, 9)
    Angle_plot, S_plot = np.meshgrid(angle_plot, s_plot)
    R, Z = polar2cart(S_plot, theta_geom_func(S_plot, Angle_plot))
    plt.plot(R.flatten(), Z.flatten(), 'o', label='via base func')
    R, Z = polar2cart(S_plot, spline(s_plot, angle_plot))
    plt.plot(R.flatten(), Z.flatten(), 'x', label='via spline')
    plt.xlabel('R [1]')
    plt.ylabel('Z [1]')
    plt.legend()
    plt.title('Theta_geom rays (uneven grid)')
    _save_current_figure('semi_periodic_theta_geom_uneven.png')


def trial_trigonometric_func(s_val, angle_val):
    return s_val * np.cos(angle_val) * np.sin(angle_val)


def theta_geom_func(s_val, chi):
    return chi


def polar2cart(rho, theta):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y


if __name__ == '__main__':
    test_init()
    test_SemiPeriodicFourierSpline_ExceedingNyquistLimit()
    test_fourier_coefs_shape()
    test_get_fourier_coefs_splines()
    test_call()
    test_get_fourier_coefs_visual_check()
    test_trial_theta_geom_func_visual_check()
    test_trial_theta_geom_func_different_angle_per_surface_visual_check()
