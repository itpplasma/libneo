# %%
import numpy as np

class StorGeom2MarsCoords:

    def __init__(self, mars_dir: str, max_poloidal_mode: int=None):
        self.max_poloidal_mode = max_poloidal_mode
        self.stor_geom_coords, self.mars_coords = self.load_points_from(
            mars_dir)
        self._set_conversion_func()

    @staticmethod
    def load_points_from(mars_dir: str):
        from omfit_classes.omfit_mars import OMFITmars
        from neo2_mars import mars_sqrtspol2stor
        mars = OMFITmars(mars_dir)
        sqrtspol = mars['sim0']['s'].values
        R, Z = mars.get_RZ()
        R = R[sqrtspol<=1.0,:]
        Z = Z[sqrtspol<=1.0,:]
        sqrtspol = sqrtspol[sqrtspol<=1.0]
        chi = mars['sim0']['R']['chi'].values
        chi_end = chi[-1] - 2*np.finfo(chi[-1]).eps
        interp = lambda fp, xp, x: np.interp(x, xp, fp)
        R[:,-1] = np.apply_along_axis(interp, axis=1, arr=R, xp=chi, x=chi_end)
        Z[:,-1] = np.apply_along_axis(interp, axis=1, arr=Z, xp=chi, x=chi_end)
        chi[-1] = chi_end
        chi = np.meshgrid(chi, sqrtspol)[0]
        stor = mars_sqrtspol2stor(mars_dir,sqrtspol)
        geom = np.arctan2(Z-Z[0,0],R-R[0,0])
        sqrtspol = sqrtspol.tolist()
        chi = chi.tolist()
        geom = geom.tolist()
        stor = stor.tolist()

        lower_bound = -np.pi
        upper_bound = +np.pi
        n = len(geom[0]) + 2
        geom[0] = np.linspace(lower_bound, upper_bound, n).tolist()
        chi[0] = np.linspace(chi[0][0], chi[0][-1], n).tolist() # at the axis the two coordinate systems agree up to a constant shift
        for i in range(1,len(sqrtspol)):
            temp_chi = chi[i].copy()
            geom_unwrapped = np.unwrap(geom[i])
            lower_bound_in_unwrapped_system = np.mod(-np.pi - geom_unwrapped[0], 2*np.pi) + geom_unwrapped[0]
            upper_bound_in_unwrapped_system = np.mod(np.pi - geom_unwrapped[0], 2*np.pi) + geom_unwrapped[0]
            geom[i].append(lower_bound)
            chi[i].append(np.interp(lower_bound_in_unwrapped_system, geom_unwrapped, temp_chi))
            geom[i].append(upper_bound)
            chi[i].append(np.interp(upper_bound_in_unwrapped_system, geom_unwrapped, temp_chi))

        mars_coords = {'radius': sqrtspol, 'angle': chi}
        stor_geom_coords = {'radius': stor, 'angle': geom}
        return order_monotonically(stor_geom_coords, mars_coords)

    def _set_conversion_func(self):
        from scipy.interpolate import CubicSpline
        self._stor2sqrtspol = CubicSpline(self.stor_geom_coords['radius'],
            self.mars_coords['radius'])
        self._stor_geom2chi = self._get_chi_func_over_stor_geom()

    def _get_chi_func_over_stor_geom(self):
        from scipy.interpolate import CubicSpline
        from .SemiPeriodicFourierSpline import SemiPeriodicFourierSpline
        stor = np.array(self.stor_geom_coords['radius'])
        chi_unwrapped = np.unwrap(self.mars_coords['angle'])
        geom = np.array(self.stor_geom_coords['angle'])
        angle_delta = chi_unwrapped - geom
        angle_delta_spline = SemiPeriodicFourierSpline(
            stor, geom, angle_delta, self.max_poloidal_mode)
        chi_spline_over_stor_geom = lambda stor, geom: np.mod(geom + angle_delta_spline(stor, geom, grid=False) + np.pi, 2*np.pi) - np.pi
        return chi_spline_over_stor_geom

    def __call__(self, stor, geom):
        sqrtspol = self._stor2sqrtspol(stor)
        chi = self._stor_geom2chi(stor, geom)
        return sqrtspol, chi


def order_monotonically(coordsystem: dict, new_coordsystem: dict)->dict:
    radius, angle = coordsystem['radius'].copy(), coordsystem['angle'].copy()
    new_radius, new_angle = new_coordsystem['radius'].copy(), new_coordsystem['angle'].copy()
    ordered_by_radius = sorted(list(zip(radius, angle, new_radius, new_angle)), key=lambda x: x[0])
    radius, angle, new_radius, new_angle = zip(*ordered_by_radius)
    radius, angle, new_radius, new_angle = list(radius), list(angle), list(new_radius), list(new_angle)
    for i in range(len(radius)):
        angle_i, new_angle_i = angle[i], new_angle[i]
        ordered_by_angle = sorted(list(zip(angle_i, new_angle_i)),
            key=lambda x: x[0])
        angle_i, new_angle_i = zip(*ordered_by_angle)
        angle[i], new_angle[i] = list(angle_i), list(new_angle_i)
    ordered_coord = {'radius': radius, 'angle': angle}
    ordered_new_coord = {'radius': new_radius, 'angle': new_angle}
    return ordered_coord, ordered_new_coord
