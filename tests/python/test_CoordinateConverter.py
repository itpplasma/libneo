# %% Standard imports
import pytest
import numpy as np
import matplotlib.pyplot as plt

# Import modules to test
from libneo import PolarCoordConverter, MarsCoords2StorThetageom

mars_dir = "/proj/plasma/DATA/DEMO/MARS/MARSQ_KNTV10_NEO2profs_KEYTORQ_1/"
sqrtspol = np.linspace(0,1,5)
chi = np.linspace(-np.pi,np.pi,10)
chi[-1] = chi[-1] - 2*np.finfo(chi[-1]).eps
coordsystem = {'radius': sqrtspol, 'angle': np.meshgrid(chi, sqrtspol)[0]}
new_coordsystem = {'radius': 1/2*sqrtspol, 'angle': np.mod(np.meshgrid(chi + np.sin(chi), sqrtspol)[0], 2*np.pi)}

def test_PolarCoordConverter_convert_unity():
    radius, angle = 0.5, -1/2*np.pi
    new_radius, new_angle = PolarCoordConverter.convert(coordsystem, coordsystem, radius, angle)
    assert new_radius == radius
    assert new_angle == angle

def test_PolarCoordConverter_convert():
    radius, angle = 0.5, -1/2*np.pi
    new_radius, new_angle = PolarCoordConverter.convert(coordsystem, new_coordsystem, radius, angle)
    assert new_radius == radius * 1/2
    assert new_angle == np.mod(angle, 2*np.pi)

def test_PolarCoordConverter_convert_output_shape():
    radius, angle = 0.5*np.ones_like(coordsystem['radius']), chi
    new_radius, new_angle = PolarCoordConverter.convert(coordsystem, new_coordsystem, radius, angle)
    assert new_radius.shape == radius.shape
    assert new_angle.shape == (len(radius), len(angle))

def test_PolarCoordConverter_convert_visual_check():
    radius, angle = 0.9*coordsystem['radius'], np.linspace(-np.pi, np.pi, 100)
    new_radius, new_angle = PolarCoordConverter.convert(coordsystem, new_coordsystem, radius, angle)
    plt.figure()
    for i in range(len(radius)):
        plt.plot(angle, new_angle[i], '.', label=f"radius={radius[i]}")
    plt.xlabel('angle [1]')
    plt.ylabel('new_angle [1]')
    plt.legend()

def test_PolarCoordConverter_convert_invers_visual_check():
    radius, angle = 0.9*new_coordsystem['radius'], np.linspace(0, 2*np.pi, 100)
    new_radius, new_angle = PolarCoordConverter.convert(new_coordsystem, coordsystem, radius, angle)
    plt.figure()
    for i in range(len(radius)):
        plt.plot(angle, new_angle[i], '.', label=f"radius={radius[i]}")
    plt.xlabel('angle [1]')
    plt.ylabel('new_angle [1]')
    plt.legend()

def test_PolarCoordConverter_order_monotonically():
    coordsystem = {'radius': np.array([1, 3, 2]), 'angle': np.array([[3, 1, 2], [7, 9, 8], [4, 5 ,6]])}
    new_coordsystem = {'radius': np.array([1, 2, 3]), 'angle': np.array([[1, 2, 3], [4, 5, 6], [7, 8 ,9]])}
    coordsystem, new_coordsystem = PolarCoordConverter.order_monotonically(coordsystem, new_coordsystem)
    assert is_monotonically_increasing(coordsystem['radius'])
    assert not is_monotonically_increasing(new_coordsystem['radius'])
    assert is_monotonically_increasing(coordsystem['angle'])
    assert not is_monotonically_increasing(new_coordsystem['angle'])
    assert np.all(coordsystem['angle'].flatten() == np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]))
    assert np.all(new_coordsystem['angle'].flatten() == np.array([2, 3, 1, 7, 8, 9, 4, 6, 5]))

def is_monotonically_increasing(arr):
    return np.all(np.diff(arr) > 0)


def test_init():
    mars = MarsCoords2StorThetageom(mars_dir)
    assert mars is not None

def test_coords_shape():
    mars = MarsCoords2StorThetageom(mars_dir)
    mars_coords = mars.coords
    new_coords = mars.new_coords
    assert mars_coords['angle'].shape[0] == len(mars_coords['radius'])
    assert mars_coords['angle'].shape == new_coords['angle'].shape
    assert mars_coords['radius'].shape == new_coords['radius'].shape

if __name__ == "__main__":
    #test_PolarCoordConverter_convert_unity()
    #test_PolarCoordConverter_convert()
    #test_PolarCoordConverter_convert_output_shape()
    test_PolarCoordConverter_convert_visual_check()
    test_PolarCoordConverter_convert_invers_visual_check()
    test_PolarCoordConverter_order_monotonically()
    test_init()
    test_coords_shape()
# %%