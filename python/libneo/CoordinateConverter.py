# %%
from abc import ABC, abstractmethod
import numpy as np

class PolarCoordConverter(ABC):
    @abstractmethod
    def _load_points_from(src: str):
        pass

    @staticmethod
    def convert(coordsystem: dict, new_coordsystem: dict, radius, angle)->np.ndarray:
        if not isinstance(radius, np.ndarray):
            radius = np.array([radius])
        if not isinstance(angle, np.ndarray):
            angle = np.array([angle])
        if ((angle.ndim == 1) & ((len(radius) != len(angle)) and (len(radius) > 1))):
            raise ValueError("The length of radius and angle should be the same.")
        if ((angle.ndim == 2) & (len(radius) != angle.shape[0])):
            raise ValueError("The angle needs to be of shape (len(radius),n_angle).")
        #coordsystem, new_coordsystem = PolarCoordConverter.order_monotonically(coordsystem, new_coordsystem)
        new_radius = np.interp(radius, coordsystem['radius'], new_coordsystem['radius'])
        new_angle = []
        for k in range(len(radius)):
            new_angle_over_radius = PolarCoordConverter.get_new_angle_on_coordsystem_radius_levels(angle[k], coordsystem, new_coordsystem)
            interp = lambda fp, xp, x: np.interp(x, xp, fp)
            new_angle.append(np.apply_along_axis(interp, axis=0, arr=new_angle_over_radius, xp=coordsystem['radius'], x=radius[k]))
        return new_radius, np.array(new_angle)

    @staticmethod
    def order_monotonically(coordsystem: dict, new_coordsystem: dict)->dict:
        radius, angle = coordsystem['radius'].copy(), coordsystem['angle'].copy()
        new_radius, new_angle = new_coordsystem['radius'].copy(), new_coordsystem['angle'].copy()
        ordered_by_radius = sorted(list(zip(radius, angle, new_radius, new_angle)), key=lambda x: x[0])
        radius, angle, new_radius, new_angle = zip(*ordered_by_radius)
        radius, angle, new_radius, new_angle = list(radius), list(angle), list(new_radius), list(new_angle)
        for i in range(len(radius)):
            angle_i, new_angle_i = angle[i], new_angle[i]
            ordered_by_angle = sorted(list(zip(angle_i, new_angle_i)), key=lambda x: x[0])
            angle_i, new_angle_i = zip(*ordered_by_angle)
            angle[i], new_angle[i] = list(angle_i), list(new_angle_i)
        ordered_coord = {'radius': radius, 'angle': angle}
        ordered_new_coord = {'radius': new_radius, 'angle': new_angle}
        return ordered_coord, ordered_new_coord

    @staticmethod
    def get_new_angle_on_coordsystem_radius_levels(angle, coordsystem, new_coordsystem):
        new_angle_over_radius = []
        for i in range(len(coordsystem['radius'])):
            #angle = PolarCoordConverter.shift_angle_to_period(angle, coordsystem['angle'][i])
            #angle = PolarCoordConverter.shift_angle_away_from_discontinuity(angle, coordsystem['angle'][i], new_coordsystem['angle'][i])
            new_angle_over_radius.append(np.interp(angle, coordsystem['angle'][i], new_coordsystem['angle'][i]))
        return np.unwrap(np.array(new_angle_over_radius), axis=0)

    @staticmethod
    def shift_angle_to_period(angle, coordsystem_angle):
        period_start = min(coordsystem_angle)
        period_end = max(np.unwrap(coordsystem_angle))
        period = period_end - period_start
        angle_in_period = np.mod(angle - period_start, period) + period_start
        return angle_in_period

    @staticmethod
    def shift_angle_away_from_discontinuity(angle, coordsystem_angle, new_coordsystem_angle):
        dicontinuity = np.where(np.abs(np.diff(new_coordsystem_angle)) > np.pi)[0]
        if len(dicontinuity) == 0:
            return angle
        angle_before_discontinuity = coordsystem_angle[dicontinuity[0]]
        angle_after_discontinuity = coordsystem_angle[dicontinuity[0]+1]
        I = is_in_discontinuity(angle, angle_before_discontinuity, angle_after_discontinuity)
        angle_shifted = angle.copy()
        if isinstance(angle, np.ndarray):
            angle_shifted[I] = angle_before_discontinuity
        elif I:
            angle_shifted = angle_before_discontinuity
        return angle_shifted
    
def is_in_discontinuity(angle, angle_before_discontinuity, angle_after_discontinuity):
    angles_over_dicontinuity = np.unwrap(np.array([angle_before_discontinuity, 
                                                   angle_after_discontinuity]))
    return np.logical_and(angles_over_dicontinuity[0] < angle, 
                          angle < angles_over_dicontinuity[1])


class MarsCoords2StorThetageom(PolarCoordConverter):
    def __init__(self, mars_dir:str):
        self._load_points_from(mars_dir)

    def _load_points_from(self, mars_dir: str):
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
        geom[sqrtspol==0,:] = chi[0].copy()
        sqrtspol = sqrtspol.tolist()
        chi = chi.tolist()
        geom = geom.tolist()
        stor = stor.tolist()

        geom_lower_bound = -np.pi
        geom_upper_bound = np.pi - 1e-10
        geom_wrap = np.unwrap(geom, axis=1)
        for i, s in enumerate(sqrtspol):
            if s == 0:
                continue
            temp_chi = chi[i].copy()
            chi[i].append(np.interp(geom_lower_bound, geom_wrap[i,:]-2*np.pi, temp_chi))
            chi[i].append(np.interp(geom_upper_bound, geom_wrap[i,:], temp_chi))
            geom[i].append(geom_lower_bound)
            geom[i].append(geom_upper_bound)
        mars_coords = {'radius': sqrtspol, 'angle': chi}
        stor_geom_coords = {'radius': stor, 'angle': geom}
        self.stor_geom_coords, self.mars_coords = PolarCoordConverter.order_monotonically(stor_geom_coords, mars_coords)

    def _set_splines(self):
        from scipy.interpolate import CubicSpline
        self.stor2sqrtspol = CubicSpline(self.stor_geom_coords['radius'], self.mars_coords['radius'])
        geom2chi = []
        for i in range(len(self.stor_geom_coords['radius'])):
            spline = CubicSplines(self.stor_geom_coords['angle'][i], self.mars_coords['angle'][i],'periodic')
            

    def stor_thetageom2mars(self, stor, theta_geom):
        sqrtspol, chi = PolarCoordConverter.convert(self.stor_geom_coords, self.mars_coords, stor, theta_geom)
        return sqrtspol, chi

if __name__=="__main__":
    mars_dir = "/proj/plasma/DATA/DEMO/MARS/MARSQ_KNTV10_NEO2profs_KEYTORQ_1/"

    