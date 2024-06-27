# %% Standard imports
import pytest
import numpy as np
import matplotlib.pyplot as plt

# Import modules to test
from libneo import StorGeom2MarsCoords
from libneo import convert, order_monotonically, shift_angle_to_period, shift_angle_away_from_discontinuity

mars_dir = "/proj/plasma/DATA/DEMO/MARS/MARSQ_KNTV10_NEO2profs_KEYTORQ_1/"
sqrtspol = np.linspace(0,1,5)
chi = np.linspace(-np.pi,np.pi,10)
chi[-1] = chi[-1] - 0.1
coordsystem = {'radius': sqrtspol.tolist(), 'angle': np.meshgrid(chi, sqrtspol)[0].tolist()}
new_coordsystem = {'radius': (1/2*sqrtspol).tolist(), 'angle': np.mod(np.meshgrid(chi, sqrtspol)[0], 2*np.pi).tolist()}

def test_PolarCoordConverter_convert_unity():
    radius, angle = 0.5, -1/2*np.pi
    new_radius, new_angle = convert(coordsystem, coordsystem, radius, angle)
    assert new_radius == radius
    assert new_angle == angle

def test_PolarCoordConverter_convert():
    radius, angle = 0.5, -1/2*np.pi
    new_radius, new_angle = convert(coordsystem, new_coordsystem, radius, angle)
    assert new_radius == radius * 1/2
    assert new_angle == np.mod(angle, 2*np.pi)

def test_PolarCoordConverter_convert_different_number_of_coordinate_points_per_radius():
    import copy
    radius, angle = 0.5, -1/2*np.pi
    modified_coordsystem = copy.deepcopy(coordsystem)
    modified_coordsystem['angle'][2] = modified_coordsystem['angle'][2][:-1]
    modified_new_coordsystem = copy.deepcopy(new_coordsystem)
    modified_new_coordsystem['angle'][2] = modified_new_coordsystem['angle'][2][:-1]
    new_radius, new_angle = convert(modified_coordsystem, modified_new_coordsystem, radius, angle)
    assert new_radius == radius * 1/2
    assert new_angle == np.mod(angle, 2*np.pi)
    for i in range(len(modified_coordsystem['angle'])):
        if i == 2:
            continue
        assert modified_coordsystem['angle'][i] == coordsystem['angle'][i]
        assert modified_new_coordsystem['angle'][i] == new_coordsystem['angle'][i]
        assert modified_coordsystem['angle'][2] != coordsystem['angle'][i]
        assert modified_new_coordsystem['angle'][2] != new_coordsystem['angle'][i]
        assert modified_coordsystem['angle'][2] != modified_coordsystem['angle'][i]
        assert modified_new_coordsystem['angle'][2] != modified_new_coordsystem['angle'][i]

def test_PolarCoordConverter_convert_output_shape():
    radius, angle = 0.5*np.ones_like(coordsystem['radius']), chi
    assert pytest.raises(ValueError, convert, coordsystem, new_coordsystem, radius, angle)
    radius, angle = 0.5*np.ones_like(new_coordsystem['radius']), chi[:5]
    new_radius, new_angle = convert(coordsystem, new_coordsystem, radius, angle)
    assert new_radius.shape == radius.shape
    assert new_angle.shape == angle.shape
    angle = np.meshgrid(chi, radius)[0]
    assert pytest.raises(ValueError, convert, coordsystem, new_coordsystem, radius, angle[:-1])
    new_radius, new_angle = convert(coordsystem, new_coordsystem, radius, angle)
    assert new_radius.shape == radius.shape
    assert new_angle.shape == (len(radius), len(angle[0]))

def test_PolarCoordConverter_order_monotonically():
    coordsystem = {'radius': [1, 3, 2], 'angle': [[3, 1, 2], [7, 9, 8], [4, 5 ,6]]}
    new_coordsystem = {'radius': [1, 2, 3], 'angle': [[1, 2, 3], [4, 5, 6], [7, 8 ,9]]}
    coordsystem, new_coordsystem = order_monotonically(coordsystem, new_coordsystem)
    assert is_monotonically_increasing(coordsystem['radius'])
    assert not is_monotonically_increasing(new_coordsystem['radius'])
    assert is_monotonically_increasing(coordsystem['angle'])
    assert not is_monotonically_increasing(new_coordsystem['angle'])
    assert np.concatenate(coordsystem['angle']).tolist() == [1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert np.concatenate(new_coordsystem['angle']).tolist() == [2, 3, 1, 7, 8, 9, 4, 6, 5]

def is_monotonically_increasing(arr):
    return np.all(np.diff(arr) > 0)

def test_PolarCoordConverter_shift_angle_to_period():
    angle = np.linspace(-np.pi, np.pi, 10) + 0.1
    coordsystem_angle = np.linspace(-np.pi, np.pi, 30)
    angle_in_period = shift_angle_to_period(angle, coordsystem_angle)
    assert np.allclose(angle_in_period[angle < np.pi], angle[angle < np.pi])
    assert np.allclose(angle_in_period[angle > np.pi], angle[angle > np.pi] - 2*np.pi)

def test_PolarCoordConverter_shift_angle_away_from_discontinuity(): 
    angle = np.linspace(0, 2*np.pi, 11)
    coord_angle = np.linspace(0, 2*np.pi, 10)
    new_coord_angle = np.concatenate([np.linspace(0, np.pi, 5), np.linspace(-np.pi, -1e-10, 5)])
    original_angle = angle.copy()
    angle = shift_angle_away_from_discontinuity(angle, coord_angle, new_coord_angle)
    assert np.allclose(original_angle[0:5], angle[0:5])
    assert ~np.allclose(original_angle[5], angle[5])
    assert np.allclose(angle[5], coord_angle[4])
    assert np.allclose(original_angle[6:], angle[6:])

def test_MarsCoords2StorThetageom_init():
    mars = StorGeom2MarsCoords(mars_dir)
    assert mars is not None

def test_MarsCoords2StorThetageom_coords_shape():
    converter = StorGeom2MarsCoords(mars_dir)
    mars_coords = converter.mars_coords
    stor_geom_coords = converter.stor_geom_coords
    assert len(mars_coords['angle']) == len(mars_coords['radius'])
    assert len(mars_coords['radius']) == len(stor_geom_coords['radius'])
    for i in range(len(mars_coords['radius'])):
        assert isinstance(mars_coords['radius'][i],float)
        assert len(mars_coords['angle'][i]) == len(stor_geom_coords['angle'][i])

def test_MarsCoord2StorThetageom_stor2sqrtspol():
    converter = StorGeom2MarsCoords(mars_dir)
    stor2sqrtspol = converter._stor2sqrtspol
    assert stor2sqrtspol(0) == 0
    assert stor2sqrtspol(1) == 1
    assert not stor2sqrtspol(0.5) == 0.5
    from neo2_mars import mars_sqrtspol2stor
    sqrtspol = np.linspace(0, 1, 100)
    stor = mars_sqrtspol2stor(mars_dir, sqrtspol)
    assert np.allclose(sqrtspol, stor2sqrtspol(stor),atol=1e-2)

def test_PolarCoordConverter_convert_visual_check():
    radius, angle = 0.9*np.array(new_coordsystem['radius']), np.linspace(-np.pi, np.pi, 100)
    angle = np.meshgrid(angle, radius)[0]
    new_radius, new_angle = convert(coordsystem, new_coordsystem, radius, angle)
    plt.figure()
    for i in range(len(radius)):
        plt.plot(angle[i], new_angle[i], '.', label=f"radius={radius[i]}")
    plt.xlabel('angle [1]')
    plt.ylabel('new_angle [1]')
    plt.axvline(x=np.min(np.array(coordsystem['angle'])), color='k', linestyle='--', label='period by provided data')
    plt.axvline(x=np.max(np.array(coordsystem['angle'])), color='k', linestyle='--')
    c, new_c = order_monotonically(coordsystem, new_coordsystem)
    discontinuity = np.where(np.abs(np.diff(new_c['angle'][0])) > np.pi)[0]
    plt.axvline(x=c['angle'][0][discontinuity[0]], color='r', linestyle='--', label='discontinuity by provided data')
    plt.axvline(x=c['angle'][0][discontinuity[0]+1], color='r', linestyle='--')
    plt.legend()
    plt.show()

def test_PolarCoordConverter_convert_invers_visual_check():
    radius, angle = 0.9*np.array(new_coordsystem['radius']), np.linspace(0, 2*np.pi, 100)
    angle = np.meshgrid(angle, radius)[0]
    new_radius, new_angle = convert(new_coordsystem, coordsystem, radius, angle)
    plt.figure()
    for i in range(len(radius)):
        plt.plot(angle[i], new_angle[i], '.', label=f"radius={radius[i]}")
    plt.xlabel('new_angle [1]')
    plt.ylabel('angle [1]')
    plt.axvline(x=np.min(np.array(new_coordsystem['angle'])), color='k', linestyle='--', label='period by provided data')
    plt.axvline(x=np.max(np.array(new_coordsystem['angle'])), color='k', linestyle='--')
    new_c, c = order_monotonically(new_coordsystem, coordsystem)
    discontinuity = np.where(np.abs(np.diff(c['angle'][0])) > np.pi)[0]
    plt.axvline(x=new_c['angle'][0][discontinuity[0]], color='r', linestyle='--', label='discontinuity by provided data')
    plt.axvline(x=new_c['angle'][0][discontinuity[0]+1], color='r', linestyle='--')
    plt.legend()
    plt.show()

def test_MarsCoords2StorThetageom_visual_check():
    converter = StorGeom2MarsCoords(mars_dir)
    test_stor = np.linspace(0.0, 1, 5)
    test_theta_geom = np.linspace(-np.pi, np.pi, 800)
    plt.figure()
    for stor in test_stor:
        sqrtspol, chi = converter(stor, test_theta_geom)
        plt.plot(test_theta_geom, chi, '.', label=f"stor={stor}")
    plt.xlabel('theta_geom [1]')
    plt.ylabel('chi [1]')
    plt.legend()
    plt.show()

def test_MarsCoords2StorThetageom_coords_domain_visual_check():
    converter = StorGeom2MarsCoords(mars_dir)
    mars_coords = converter.mars_coords
    stor_geom_coords = converter.stor_geom_coords
    fig, ax = plt.subplots(1,2, figsize=(10,5))
    color = ['b', 'y', 'g', 'r']

    stor_geom_coords, mars_coords = order_monotonically(stor_geom_coords, mars_coords)
    color_count = 0
    for i in range(0,len(mars_coords['radius']),100):
        ax[0].plot(stor_geom_coords['angle'][i], mars_coords['angle'][i], '.' + color[color_count], label=f"sqrtspol={mars_coords['radius'][i]}")
        discontinuity = np.where(np.abs(np.diff(mars_coords['angle'][i])) > np.pi)[0]
        if len(discontinuity) > 0:
            ax[0].axvline(x=stor_geom_coords['angle'][i][discontinuity[0]], color=color[color_count], linestyle='--', label='discontinuity by data')
            ax[0].axvline(x=stor_geom_coords['angle'][i][discontinuity[0]+1], color=color[color_count], linestyle='--')
        color_count += 1
    ax[0].set_xlabel('theta_geom [1]')
    ax[0].set_ylabel('chi [1]')
    ax[0].legend()

    mars_coords, stor_geom_coords = order_monotonically(mars_coords, stor_geom_coords)
    color_count = 0
    for i in range(0,len(mars_coords['radius']),100):
        ax[1].plot(mars_coords['angle'][i], stor_geom_coords['angle'][i], '.' + color[color_count], label=f"stor={stor_geom_coords['radius'][i]}")
        discontinuity = np.where(np.abs(np.diff(stor_geom_coords['angle'][i])) > np.pi)[0]
        if len(discontinuity) > 0:
            ax[1].axvline(x=mars_coords['angle'][i][discontinuity[0]], color=color[color_count], linestyle='--', label='discontinuity by data')
            ax[1].axvline(x=mars_coords['angle'][i][discontinuity[0]+1], color=color[color_count], linestyle='--')
        color_count += 1
    ax[1].set_xlabel('chi [1]')
    ax[1].set_ylabel('theta_geom [1]')
    ax[1].legend()
    plt.show()

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
            sqrtspol, chi = converter(stor, test_theta_geom)


if __name__ == "__main__":
    test_PolarCoordConverter_convert_unity()
    test_PolarCoordConverter_convert()
    test_PolarCoordConverter_convert_different_number_of_coordinate_points_per_radius()
    test_PolarCoordConverter_convert_output_shape()
    test_PolarCoordConverter_order_monotonically()
    test_PolarCoordConverter_shift_angle_to_period()
    test_PolarCoordConverter_shift_angle_away_from_discontinuity()
    test_MarsCoords2StorThetageom_init()
    test_MarsCoords2StorThetageom_coords_shape()
    test_MarsCoord2StorThetageom_stor2sqrtspol()
    print("All tests passed!")

    #test_MarsCoord2StorThetageom_performance_profile()

    test_PolarCoordConverter_convert_visual_check()
    test_PolarCoordConverter_convert_invers_visual_check()
    test_MarsCoords2StorThetageom_visual_check()
    test_MarsCoords2StorThetageom_coords_domain_visual_check()