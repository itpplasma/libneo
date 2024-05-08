# %% Standard imports
import numpy as np

# modules to test
from libneo import fourier_coefs_half, fourier_coefs_full

def trial_function(x):
    return np.sin(x)*np.sin(x)*np.cos(x) + 1

def complex_trial_function(x):
    return np.sin(x) + np.sin(x)*np.cos(x)*1.0j

def test_fourier_coefs_half():
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(0, 5)
    fm = fourier_coefs_half(trial_function(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert not is_imag_part_vanisching(f_eval)
    assert np.allclose(np.real(f_eval), trial_function(new_x))

def test_fourier_coefs_full():
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(-4, 5)
    fm = fourier_coefs_full(trial_function(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert is_imag_part_vanisching(f_eval)
    assert np.allclose(np.real(f_eval), trial_function(new_x))

def test_fourier_coefs_half_on_complex_function():
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(0, 5)
    fm = fourier_coefs_half(complex_trial_function(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert not is_imag_part_vanisching(f_eval)
    assert not np.allclose(np.real(f_eval), np.real(complex_trial_function(new_x)))
    assert not np.allclose(np.imag(f_eval), np.imag(complex_trial_function(new_x)))

def test_fourier_coefs_full_on_complex_function():
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(-4, 5)
    fm = fourier_coefs_full(complex_trial_function(x)[:-1], x, m)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = fm.dot(np.exp(1.0j * np.outer(m, new_x)))
    assert is_complex(f_eval)
    assert not is_imag_part_vanisching(f_eval)
    assert np.allclose(np.real(f_eval), np.real(complex_trial_function(new_x)))
    assert np.allclose(np.imag(f_eval), np.imag(complex_trial_function(new_x)))

def test_fourier_coefs_visual_check():
    import matplotlib.pyplot as plt
    x = np.linspace(0, 2*np.pi, 100)
    m = np.arange(-4, 5)
    fm = fourier_coefs_full(trial_function(x)[:-1], x, m)
    m_half = np.arange(0, 5)
    fm_half = fourier_coefs_half(trial_function(x)[:-1], x, m_half)
    new_x = np.sort(np.random.uniform(0, 2*np.pi, 20))
    f_eval = np.real(fm.dot(np.exp(1.0j * np.outer(m, new_x))))
    f_eval_half = np.real(fm_half.dot(np.exp(1.0j * np.outer(m_half, new_x))))
    plt.figure()
    plt.plot(new_x, trial_function(new_x), label='trial function')
    plt.plot(new_x, f_eval, 'o', label='full fourier')
    plt.plot(new_x, f_eval_half, 'x', label='half fourier')
    plt.legend()
    plt.show()

def is_imag_part_vanisching(x):
    return np.allclose(np.imag(x), 0)

def is_complex(x):
    return np.iscomplexobj(x)

if __name__ == '__main__':
    test_fourier_coefs_half()
    test_fourier_coefs_full()
    test_fourier_coefs_half_on_complex_function()
    test_fourier_coefs_full_on_complex_function()
    print("All tests passed!")
    #test_fourier_coefs_visual_check()