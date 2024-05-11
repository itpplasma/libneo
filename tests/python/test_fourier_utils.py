# %% Standard imports
import numpy as np
import matplotlib.pyplot as plt

# modules to test
from libneo import fourier_coefs_half, fourier_coefs_full, get_half_fft

x = np.linspace(0, 2*np.pi, 100)

def test_fourier_coefs_half():
    m = np.arange(0, 5)
    fm = fourier_coefs_half(trial_trigonometric_func(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert not is_imag_part_vanisching(f_eval)
    assert np.allclose(np.real(f_eval), trial_trigonometric_func(new_x))

def test_fourier_coefs_full():
    m = np.arange(-4, 5)
    fm = fourier_coefs_full(trial_trigonometric_func(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert is_imag_part_vanisching(f_eval)
    assert np.allclose(np.real(f_eval), trial_trigonometric_func(new_x))

def test_fourier_coefs_half_on_complex_function():
    m = np.arange(0, 5)
    fm = fourier_coefs_half(complex_trial_func(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert not is_imag_part_vanisching(f_eval)
    assert not np.allclose(np.real(f_eval), np.real(complex_trial_func(new_x)))
    assert not np.allclose(np.imag(f_eval), np.imag(complex_trial_func(new_x)))

def test_fourier_coefs_full_on_complex_function():
    m = np.arange(-4, 5)
    fm = fourier_coefs_full(complex_trial_func(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    fft_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(fft_eval)
    assert not is_imag_part_vanisching(fft_eval)
    assert np.allclose(np.real(fft_eval), np.real(complex_trial_func(new_x)))
    assert np.allclose(np.imag(fft_eval), np.imag(complex_trial_func(new_x)))

def test_get_half_fft_coefficients_on_trigonometric_func():
    grid = np.linspace(0, 2*np.pi, 20)
    fm, m = get_half_fft(trial_trigonometric_func(grid), grid)
    # The coefficients from the trial function sin(x)*sin(x)*cos(x) + 1
    assert np.allclose(fm.real, np.array([1,0.25,0,-0.25,0,0,0,0,0,0]))
    assert np.allclose(fm.imag, np.array([0,0,0,0,0,0,0,0,0,0]))

def test_get_half_fft_even_odd_grid():
    even_grid = np.linspace(-np.pi, np.pi, 20)
    odd_grid = np.linspace(-np.pi, np.pi, 21)
    f_even = get_fft_series_lambda_function(trial_func, even_grid)
    f_odd = get_fft_series_lambda_function(trial_func, odd_grid)
    assert np.allclose(f_even(even_grid), trial_func(even_grid))
    assert not np.allclose(f_even(odd_grid), trial_func(odd_grid))
    assert np.allclose(f_odd(odd_grid), trial_func(odd_grid))
    assert not np.allclose(f_odd(even_grid), trial_func(even_grid))

def test_get_half_fft_on_trigonometric_func():
    grid = np.linspace(0, 2*np.pi, 20)
    fm, m = get_half_fft(trial_trigonometric_func(grid), grid)
    grid_eval = np.linspace(0, 2*np.pi, 100)
    fft_eval = fm.dot(np.exp(1.0j * np.outer(m, grid_eval)))
    assert is_complex(fft_eval)
    assert not is_imag_part_vanisching(fft_eval)
    assert np.allclose(np.real(fft_eval), trial_trigonometric_func(grid_eval))

def test_get_half_fft_on_func():
    grid = np.linspace(-np.pi, np.pi, 20)
    between_grid = grid[:-1] + 0.5*np.diff(grid)
    fft_series = get_fft_series_lambda_function(trial_func, grid)
    assert np.allclose(fft_series(grid), trial_func(grid))
    assert not np.allclose(fft_series(between_grid), trial_func(between_grid))

def test_get_half_fft_phaseshift():
    grid = np.linspace(0, 2*np.pi, 10)
    f_minus_pi = get_fft_series_lambda_function(trial_func, grid - np.pi)
    f_plus_pi = get_fft_series_lambda_function(trial_func, grid + np.pi)
    grid_eval = add_intermediate_steps(grid)
    assert np.allclose(f_minus_pi(grid_eval), f_plus_pi(grid_eval))
    grid_plus_pi = grid + np.pi
    assert np.allclose(f_plus_pi(grid_plus_pi), trial_func(grid_plus_pi))
    grid_plus_pi_between = get_intermediate_steps(grid_plus_pi)
    assert not np.allclose(f_plus_pi(grid_plus_pi_between), trial_func(grid_plus_pi_between))   

def test_get_half_fft_phaseshift_visual_check():
    grid = np.linspace(0, 2*np.pi, 10)
    f = get_fft_series_lambda_function(trial_func, grid)
    f_minus_pi = get_fft_series_lambda_function(trial_func, grid - np.pi)
    f_plus_pi = get_fft_series_lambda_function(trial_func, grid + np.pi)
    grid_eval = add_intermediate_steps(grid)
    plt.figure()
    plt.plot(grid_eval, trial_func(grid_eval), 'ro-', label='base function')
    plt.plot(grid_eval, f(grid_eval), 'D-', label='no shift')
    plt.plot(grid_eval, f_plus_pi(grid_eval), 'x-', label='+pi shifted')
    plt.plot(grid_eval, f_minus_pi(grid_eval),'+--', label='-pi shifted')
    plt.legend()
    plt.xlabel('x [1]')
    plt.ylabel('f(x) [1]')
    plt.title('Reconstruction of trial function from ffts on shifted intervals')
    plt.show()

def test_get_half_fft_even_odd_grid_visual_check():
    even_grid = np.linspace(-np.pi, np.pi, 10)
    odd_grid = np.linspace(-np.pi, np.pi, 11)
    f_even = get_fft_series_lambda_function(trial_func, even_grid)
    f_odd = get_fft_series_lambda_function(trial_func, odd_grid)
    grid_eval = np.sort(np.concatenate([even_grid, odd_grid]))
    plt.figure()
    plt.plot(grid_eval, trial_func(grid_eval), 'D-', label='base function')
    plt.plot(grid_eval, f_even(grid_eval), 'x-', label='fft on even grid')
    plt.plot(grid_eval, f_odd(grid_eval), '+-', label='fft on odd grid')
    plt.legend()
    plt.xlabel('x [1]')
    plt.ylabel('f(x) [1]')
    plt.title('Reconstruction of trial function from ffts on even and odd grid')
    plt.show()

def get_fft_series_lambda_function(func, grid):
    fm, m = get_half_fft(func(grid), grid)
    return lambda x: np.real(fm.dot(np.exp(1.0j * np.outer(m, x))))

def add_intermediate_steps(grid):
    return np.sort(np.concatenate([grid,get_intermediate_steps(grid)]))

def get_intermediate_steps(grid):
    return grid[:-1] + 0.5*np.diff(grid)

def test_fourier_coefs_visual_check():
    m = np.arange(-4, 5)
    fm = fourier_coefs_full(trial_trigonometric_func(x)[:-1], x, m)
    m_half = np.arange(0, 5)
    fm_half = fourier_coefs_half(trial_trigonometric_func(x)[:-1], x, m_half)
    x_rand = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = np.real(fm.dot(np.exp(1.0j * np.outer(m, x_rand))))
    f_eval_half = np.real(fm_half.dot(np.exp(1.0j * np.outer(m_half, x_rand))))
    plt.figure()
    plt.plot(x, trial_trigonometric_func(x), label='trial function')
    plt.plot(x_rand, f_eval, 'o', label='full fourier eval')
    plt.plot(x_rand, f_eval_half, 'x', label='half fourier eval')
    plt.legend()
    plt.xlabel('x [1]')
    plt.ylabel('f(x) [1]')
    plt.title('Reconstruction of trial function from fourier coef computation')
    plt.show()

def test_get_half_fft_visual_check():
    make_visual_check_get_half_fft_on(trial_trigonometric_func)
    make_visual_check_get_half_fft_on(trial_func)

def make_visual_check_get_half_fft_on(func):
    grid = np.linspace(0, 2*np.pi, 20)
    fm, m = get_half_fft(func(grid), grid)
    grid_eval = np.sort(np.concatenate([grid,grid[:-1] + 0.5*np.diff(grid)]))
    fourier_eval = fm.dot(np.exp(1.0j * np.outer(m, grid_eval)))
    plt.figure()
    plt.plot(grid, func(grid), 'D', label='trial data')
    plt.plot(grid_eval, np.real(fourier_eval), '-x', label='half fft eval')
    plt.legend()
    plt.xlabel('x [1]')
    plt.ylabel('f(x) [1]')
    plt.title('Reconstruction of trial function from fft results')
    plt.show()

def test_get_half_fft_coefficients_visual_check():
    grid = np.linspace(0, 2*np.pi, 20)
    fm, m = get_half_fft(trial_trigonometric_func(grid), grid)
    plt.figure()
    plt.plot(m, fm.real, '-D', label='real-part fft coefs')
    plt.plot(m, fm.imag, '-x', label='imag-part fft coefs')
    plt.legend()
    plt.xlabel('modenumber m [1]')
    plt.ylabel('fm [1]')
    plt.title('Fourier coefficients of trigonometric function from FFT')
    plt.show()

def trial_trigonometric_func(x):
    return np.sin(x)*np.sin(x)*np.cos(x) + 1

def trial_func(x):
    return np.mod(x, 2*np.pi)

def complex_trial_func(x):
    return np.sin(x) + np.sin(x)*np.cos(x)*1.0j

def is_imag_part_vanisching(x):
    return np.allclose(np.imag(x), 0)

def is_complex(x):
    return np.iscomplexobj(x)

if __name__ == '__main__':
    test_fourier_coefs_half()
    test_fourier_coefs_full()
    test_fourier_coefs_half_on_complex_function()
    test_fourier_coefs_full_on_complex_function()
    test_get_half_fft_coefficients_on_trigonometric_func()
    test_get_half_fft_even_odd_grid()
    test_get_half_fft_on_trigonometric_func()
    test_get_half_fft_on_func()
    test_get_half_fft_phaseshift()
    print("All tests passed!")
    test_get_half_fft_phaseshift_visual_check()
    test_get_half_fft_even_odd_grid_visual_check()
    test_fourier_coefs_visual_check()
    test_get_half_fft_visual_check()
    test_get_half_fft_coefficients_visual_check()