# %%
from abc import ABC, abstractmethod
import numpy as np

class PolarCoordConverter(ABC):
    @abstractmethod
    def _load_points_from(src: str):
        pass

    @staticmethod
    def convert(coordsystem: dict, new_coordsystem: dict, radius: float, angle: float)->np.ndarray:
        coordsystem, new_coordsystem = PolarCoordConverter.order_monotonically(coordsystem, new_coordsystem)
        new_radius = np.interp(radius, coordsystem['radius'], new_coordsystem['radius'])
        new_angle_over_old_radius = []
        for i in range(len(coordsystem['radius'])):
            new_angle_over_old_radius.append(np.interp(angle, coordsystem['angle'][i], new_coordsystem['angle'][i]))
        interp = lambda fp, xp, x: np.interp(x, xp, fp)
        new_angle = np.apply_along_axis(interp, axis=0, arr=new_angle_over_old_radius, xp=coordsystem['radius'], x=radius)
        return new_radius, new_angle

    @staticmethod
    def order_monotonically(coord: dict, new_coord: dict)->dict:
        ordered_coord, ordered_new_coord = coord.copy(), new_coord.copy()
        id_radius = np.argsort(ordered_coord['radius'])
        ordered_coord['radius'] = ordered_coord['radius'][id_radius]
        ordered_coord['angle'] = ordered_coord['angle'][id_radius]
        ordered_new_coord['radius'] = ordered_new_coord['radius'][id_radius]
        ordered_new_coord['angle'] = ordered_new_coord['angle'][id_radius]
        id_angle = np.argsort(ordered_coord['angle'], axis=1)
        ordered_coord['angle'] = np.take_along_axis(ordered_coord['angle'], id_angle, axis=1)
        ordered_new_coord['angle'] = np.take_along_axis(ordered_new_coord['angle'], id_angle, axis=1)
        return ordered_coord, ordered_new_coord

def len_including_scalar(arr):
    try:
        return len(arr)
    except TypeError:
        return None

            
def interpolate2D(xy:np.ndarray, z:np.ndarray, x0y0:np.ndarray)->float:
    from numpy import interp
    z_at_x0_along_y = interp(x0[0], xy[0], z)
    z_at_x0_y0 = interp(x0y0[1], xy[1], z_at_x0_along_y)
    return z_at_x0_y0

def convert2Dcoords(old_coords:np.ndarray, new_coords:np.ndarray, point_oldcoords:np.ndarray)->np.ndarray: 
    point_new_coords = np.array([np.nan, np.nan])
    point_new_coords[0] = interpolate2D(xy=old_coords,z=new_coords[0], x0y0=point_oldcoords)
    point_new_coords[1] = interpolate2D(xy=old_coords,z=new_coords[1], x0y0=point_oldcoords)
    return point_new_coords

class MarsCoords2StorThetageom(PolarCoordConverter):
    def __init__(self, src:str):
        self._load_points_from(src)

    def _load_points_from(self, mars_dir: str):
        from omfit_classes.omfit_mars import OMFITmars
        from neo2_mars import mars_sqrtspol2stor
        mars = OMFITmars(mars_dir)
        sqrtspol = mars['sim0']['s']
        R, Z = mars.get_RZ()
        R = R[sqrtspol<=1.0,:]
        Z = Z[sqrtspol<=1.0,:]
        sqrtspol = sqrtspol[sqrtspol<=1.0]
        chi = mars['sim0']['R']['chi'].values
        chi = np.meshgrid(chi,sqrtspol)[0]
        stor = mars_sqrtspol2stor(mars_dir,sqrtspol)
        theta_geom = np.arctan2(Z-Z[0],R-R[0])
        theta_geom = np.mod(theta_geom,2*np.pi)
        idx = np.argsort(theta_geom, axis=1)
        theta_geom = np.take_along_axis(theta_geom, idx, axis=1)
        chi = np.take_along_axis(chi, idx, axis=1)
        self.idx = idx
        self.coords = {'radius': sqrtspol, 'angle': chi}
        self.new_coords = {'radius': stor, 'angle': theta_geom}

if __name__=="__main__":
    mars_dir = "/proj/plasma/DATA/DEMO/MARS/MARSQ_KNTV10_NEO2profs_KEYTORQ_1/"

    