# %% Standard imports
import numpy as np
import matplotlib.pyplot as plt
import pytest

# modules to test
from libneo import SemiPeriodicFourierSpline
from libneo import SemiPeriodicFourierSpline_ExceedingNyquistLimit

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
        spline = SemiPeriodicFourierSpline(s, angle, F, max_modenumber=len(angle))

def test_fourier_coefs_shape():
    max_m = 4
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F, max_modenumber=max_m)
    coefs = spline._fourier_coefs
    assert coefs.shape == (len(s), max_m+1)

def test_get_fourier_coefs():
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F)
    assert np.allclose(spline._get_fourier_coefs(s), spline._fourier_coefs)
    assert spline._get_fourier_coefs(s).shape == (len(s), len(spline.modenumbers))

def test_call():
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F)
    assert np.allclose(spline(s, angle), F)
    assert np.allclose(spline(s, angle+2*np.pi), spline(s, angle))
    assert np.allclose(spline(s, angle-2*np.pi), spline(s, angle))
    r_new = np.random.uniform(0, 1, 5)
    angle_new = np.random.uniform(-np.pi, np.pi, 5)
    assert np.allclose(spline(r_new, angle_new,grid=False), trial_trigonometric_func(r_new, angle_new))

def test_get_fourier_coefs_visual_check():
    Angle, S = np.meshgrid(angle, s)
    F = trial_trigonometric_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, F)
    plt.figure()
    r_plot = np.linspace(0, 1, 30)
    coefs_plot = spline._get_fourier_coefs(r_plot)
    for i,m in enumerate(spline.modenumbers):
        plt.plot(r_plot, coefs_plot[:,i].real, '-', label=f'real m={m}')
        plt.plot(r_plot, coefs_plot[:,i].imag, '.--', label=f'imag m={m}')
    plt.xlabel('s [1]')
    plt.ylabel('fm [1]')
    plt.title('Fourier coefficients of trigonometric trial function depending on s')
    plt.legend()
    plt.show()

def test_trial_theta_geom_func_visual_check():
    angle = np.linspace(-np.pi, np.pi, 100)
    Angle, S = np.meshgrid(angle, s)
    Theta_geom = theta_geom_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, angle, Theta_geom)
    plt.figure()
    angle_plot = np.linspace(-np.pi, np.pi,19)
    s_plot = np.linspace(0, 1, 9)
    Angle_plot, S_plot = np.meshgrid(angle_plot, s_plot)
    R, Z = polar2cart(S_plot, theta_geom_func(S_plot, Angle_plot))
    plt.plot(R.flatten(), Z.flatten(), 'o', label='via base func')
    R, Z = polar2cart(S_plot, spline(s_plot, angle_plot))
    plt.plot(R.flatten(), Z.flatten(), 'x', label='via spline')
    plt.xlabel('R [1]')
    plt.ylabel('Z [1]')
    plt.legend()
    plt.title('Theta_geom rays as function of equidistant input angle \n splined form uniform datapoints N='+str(len(angle))+'x'+str(len(s))+' [angle x s]')
    plt.show()

def test_trial_theta_geom_func_different_angle_per_surface_visual_check():
    angle = np.linspace(-np.pi, np.pi, 100)
    Angle, S = np.meshgrid(angle, s)
    Angle[:,1:-1] = Angle[:,1:-1]*(S[:,1:-1] + 1)/2 # uneven distributed angles, but same endpoints
    Theta_geom = theta_geom_func(S, Angle)
    spline = SemiPeriodicFourierSpline(s, Angle, Theta_geom)
    plt.figure()
    angle_plot = np.linspace(-np.pi, np.pi,19)
    s_plot = np.linspace(0, 1, 9)
    Angle_plot, S_plot = np.meshgrid(angle_plot, s_plot)
    R, Z = polar2cart(S_plot, theta_geom_func(S_plot, Angle_plot))
    plt.plot(R.flatten(), Z.flatten(), 'o', label='via base func')
    R, Z = polar2cart(S_plot, spline(s_plot, angle_plot))
    plt.plot(R.flatten(), Z.flatten(), 'x', label='via spline')
    plt.xlabel('R [1]')
    plt.ylabel('Z [1]')
    plt.legend()
    plt.title('Theta_geom rays as function of equidistant input angle \n splined from uneven distributed datapoints N='+str(len(angle))+'x'+str(len(s))+' [angle x s]')
    plt.show()

def trial_trigonometric_func(s,angle):
    return s*np.cos(angle)*np.sin(angle)

def theta_geom_func(s,chi):
    theta_geom = chi
    return theta_geom

def polar2cart(rho,theta):
    x = rho*np.cos(theta)
    y = rho*np.sin(theta)
    return x,y

if __name__ == '__main__':
    test_init()
    test_SemiPeriodicFourierSpline_ExceedingNyquistLimit()
    test_fourier_coefs_shape()
    test_get_fourier_coefs()
    test_call()
    print("All tests passed.")
    test_trial_theta_geom_func_visual_check()
    test_trial_theta_geom_func_different_angle_per_surface_visual_check()
    test_get_fourier_coefs_visual_check()