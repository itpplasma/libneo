# %% Standard imports
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# modules to test
from libneo import fourier_coefs_half, fourier_coefs_full, get_half_fft

def save_plot(fig, filename):
    """Helper to save plots to build/test/python directory"""
    output = Path(f"build/test/python/{filename}")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150, bbox_inches='tight')
    return output

def test_fourier_coefs_half():
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(0, 5)
    fm = fourier_coefs_half(trial_trigonometric_func(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert not is_imag_part_vanisching(f_eval)
    assert np.allclose(np.real(f_eval), trial_trigonometric_func(new_x))

def test_fourier_coefs_full():
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(-4, 5)
    fm = fourier_coefs_full(trial_trigonometric_func(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert is_imag_part_vanisching(f_eval)
    assert np.allclose(np.real(f_eval), trial_trigonometric_func(new_x))

def test_fourier_coefs_half_on_complex_function():
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(0, 5)
    fm = fourier_coefs_half(complex_trial_func(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert not is_imag_part_vanisching(f_eval)
    assert not np.allclose(np.real(f_eval), np.real(complex_trial_func(new_x)))
    assert not np.allclose(np.imag(f_eval), np.imag(complex_trial_func(new_x)))

def test_fourier_coefs_full_on_complex_function():
    x = np.linspace(0, 2*np.pi, 100)
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
    fig = plt.figure()
    plt.plot(grid_eval, trial_func(grid_eval), 'ro-', label='base function')
    plt.plot(grid_eval, f(grid_eval), 'D-', label='no shift')
    plt.plot(grid_eval, f_plus_pi(grid_eval), 'x-', label='+pi shifted')
    plt.plot(grid_eval, f_minus_pi(grid_eval),'+--', label='-pi shifted')
    plt.legend()
    plt.xlabel('x [1]')
    plt.ylabel('f(x) [1]')
    plt.title('Reconstruction of trial function from ffts on shifted intervals')
    save_plot(fig, "fourier_phaseshift_visual.png")
    plt.show()
    plt.close(fig)

def test_get_half_fft_even_odd_grid_visual_check():
    even_grid = np.linspace(-np.pi, np.pi, 10)
    odd_grid = np.linspace(-np.pi, np.pi, 11)
    f_even = get_fft_series_lambda_function(trial_func, even_grid)
    f_odd = get_fft_series_lambda_function(trial_func, odd_grid)
    grid_eval = np.sort(np.concatenate([even_grid, odd_grid]))
    fig = plt.figure()
    plt.plot(grid_eval, trial_func(grid_eval), 'D-', label='base function')
    plt.plot(grid_eval, f_even(grid_eval), 'x-', label='fft on even grid')
    plt.plot(grid_eval, f_odd(grid_eval), '+-', label='fft on odd grid')
    plt.legend()
    plt.xlabel('x [1]')
    plt.ylabel('f(x) [1]')
    plt.title('Reconstruction of trial function from ffts on even and odd grid')
    save_plot(fig, "fourier_even_odd_grid_visual.png")
    plt.show()
    plt.close(fig)

def get_fft_series_lambda_function(func, grid):
    fm, m = get_half_fft(func(grid), grid)
    return lambda x: np.real(fm.dot(np.exp(1.0j * np.outer(m, x))))

def add_intermediate_steps(grid):
    return np.sort(np.concatenate([grid,get_intermediate_steps(grid)]))

def get_intermediate_steps(grid):
    return grid[:-1] + 0.5*np.diff(grid)

def test_fourier_coefs_visual_check():
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(-4, 5)
    fm = fourier_coefs_full(trial_trigonometric_func(x)[:-1], x, m)
    m_half = np.arange(0, 5)
    fm_half = fourier_coefs_half(trial_trigonometric_func(x)[:-1], x, m_half)
    x_rand = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = np.real(fm.dot(np.exp(1.0j * np.outer(m, x_rand))))
    f_eval_half = np.real(fm_half.dot(np.exp(1.0j * np.outer(m_half, x_rand))))
    fig = plt.figure()
    plt.plot(x, trial_trigonometric_func(x), label='trial function')
    plt.plot(x_rand, f_eval, 'o', label='full fourier eval')
    plt.plot(x_rand, f_eval_half, 'x', label='half fourier eval')
    plt.legend()
    plt.xlabel('x [1]')
    plt.ylabel('f(x) [1]')
    plt.title('Reconstruction of trial function from fourier coef computation')
    save_plot(fig, "fourier_coefs_visual.png")
    plt.show()
    plt.close(fig)

def test_get_half_fft_visual_check():
    make_visual_check_get_half_fft_on(trial_trigonometric_func)
    make_visual_check_get_half_fft_on(trial_func)

def make_visual_check_get_half_fft_on(func):
    grid = np.linspace(0, 2*np.pi, 20)
    fm, m = get_half_fft(func(grid), grid)
    grid_eval = np.sort(np.concatenate([grid,grid[:-1] + 0.5*np.diff(grid)]))
    fourier_eval = fm.dot(np.exp(1.0j * np.outer(m, grid_eval)))
    fig = plt.figure()
    plt.plot(grid, func(grid), 'D', label='trial data')
    plt.plot(grid_eval, np.real(fourier_eval), '-x', label='half fft eval')
    plt.legend()
    plt.xlabel('x [1]')
    plt.ylabel('f(x) [1]')
    plt.title('Reconstruction of trial function from fft results')
    func_name = func.__name__
    save_plot(fig, f"fourier_half_fft_{func_name}_visual.png")
    plt.show()
    plt.close(fig)

def test_get_half_fft_coefficients_visual_check():
    grid = np.linspace(0, 2*np.pi, 20)
    fm, m = get_half_fft(trial_trigonometric_func(grid), grid)
    fig = plt.figure()
    plt.plot(m, fm.real, '-D', label='real-part fft coefs')
    plt.plot(m, fm.imag, '-x', label='imag-part fft coefs')
    plt.legend()
    plt.xlabel('modenumber m [1]')
    plt.ylabel('fm [1]')
    plt.title('Fourier coefficients of trigonometric function from FFT')
    save_plot(fig, "fourier_coefficients_visual.png")
    plt.show()
    plt.close(fig)

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

def test_angle_conversion_example_visual_check():
    chi = np.linspace(-np.pi, np.pi, 50)
    boundary = -np.pi + 2
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    geom = trial_geom_func(chi, boundary)
    ax[0].plot(chi, geom, 'b-', label='geom angle function')
    ax[0].set_xlabel('chi [1]')
    ax[0].set_ylabel('geom [1]')

    geom_ordered, chi_ordered = sort_by_first_list(geom, chi)
    chi_ordered_unwrapped = np.unwrap(chi_ordered)
    geom_equi = np.linspace(-np.pi, np.pi, 50)
    chi_ordered_unwrapped_equi = np.interp(geom_equi, geom_ordered, chi_ordered_unwrapped)
    ax[1].plot(geom_ordered, chi_ordered_unwrapped, 'ko', label=' unwrapped chi angle function')
    ax[1].plot(geom_equi, chi_ordered_unwrapped_equi, 'g.', label='interpolated unwrapped chi function')
    ax[1].plot(geom_ordered, chi_ordered, 'b-', label='chi angle function')

    delta_ordered_unwrapped_equi = chi_ordered_unwrapped_equi - geom_equi
    ax[1].plot(geom_equi, delta_ordered_unwrapped_equi, 'm.', label='chi_unwrapped - geom')
    ax[1].axhline(y=-np.pi, color='r', linestyle='--', label='period')
    ax[1].axhline(y=np.pi, color='r', linestyle='--')
    ax[1].legend()
    ax[1].set_xlabel('geom [1]')
    ax[1].set_ylabel('chi [1]')
    ax[1].axis('equal')

    delta_m, m = get_half_fft(delta_ordered_unwrapped_equi, geom_equi)
    delta_fourier = lambda geom: np.real(delta_m.dot(np.exp(1.0j * np.outer(m, geom))))
    geom_eval = np.linspace(-np.pi, np.pi, 100)
    chi_eval = np.mod(delta_fourier(geom_eval) + geom_eval + np.pi, 2*np.pi) - np.pi
    ax[0].plot(chi_eval, geom_eval, 'g.', label='fourier reconstruction (geom_equi)')
    ax[0].axvline(x=-np.pi, color='r', linestyle='--', label='period')
    ax[0].axvline(x=np.pi, color='r', linestyle='--')
    ax[0].legend()
    ax[0].axis('equal')
    plt.suptitle('Conversion between geom and chi angle functions')
    save_plot(fig, "fourier_angle_conversion_visual.png")
    plt.show()
    plt.close(fig)

def trial_geom_func(x, boundary: float=-np.pi):
    angle = (x - np.pi - boundary) + np.sin((x - np.pi - boundary)) + 0.3*np.sin(2*(x - np.pi - boundary))
    return np.mod(angle + np.pi, 2*np.pi) - np.pi

def sort_by_first_list(list1, list2):
    if isinstance(list1, np.ndarray):
        list1 = list1.tolist()
    if isinstance(list2, np.ndarray):
        list2 = list2.tolist()
    total_list = list(zip(list1, list2))
    total_list.sort(key=lambda x: x[0])
    tupel1, tupel2 = zip(*total_list)
    return list(tupel1), list(tupel2)


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
    test_angle_conversion_example_visual_check()