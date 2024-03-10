from __future__ import print_function, absolute_import, division
import _interpolate
import f90wrap.runtime
import logging
import numpy

class Interpolate(f90wrap.runtime.FortranModule):
    """
    Module interpolate
    
    
    Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 5-405
    
    """
    @f90wrap.runtime.register_class("interpolate.SplineData1D")
    class SplineData1D(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=splinedata1d)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 8-14
        
        """
        def __init__(self, handle=None):
            """
            self = Splinedata1D()
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 8-14
            
            
            Returns
            -------
            this : Splinedata1D
            	Object to be constructed
            
            
            Automatically generated constructor for splinedata1d
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _interpolate.f90wrap_splinedata1d_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Splinedata1D
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 8-14
            
            Parameters
            ----------
            this : Splinedata1D
            	Object to be destructed
            
            
            Automatically generated destructor for splinedata1d
            """
            if self._alloc:
                _interpolate.f90wrap_splinedata1d_finalise(this=self._handle)
        
        @property
        def order(self):
            """
            Element order ftype=integer  pytype=int
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 9
            
            """
            return _interpolate.f90wrap_splinedata1d__get__order(self._handle)
        
        @order.setter
        def order(self, order):
            _interpolate.f90wrap_splinedata1d__set__order(self._handle, order)
        
        @property
        def num_points(self):
            """
            Element num_points ftype=integer  pytype=int
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 10
            
            """
            return _interpolate.f90wrap_splinedata1d__get__num_points(self._handle)
        
        @num_points.setter
        def num_points(self, num_points):
            _interpolate.f90wrap_splinedata1d__set__num_points(self._handle, num_points)
        
        @property
        def periodic(self):
            """
            Element periodic ftype=logical pytype=bool
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 11
            
            """
            return _interpolate.f90wrap_splinedata1d__get__periodic(self._handle)
        
        @periodic.setter
        def periodic(self, periodic):
            _interpolate.f90wrap_splinedata1d__set__periodic(self._handle, periodic)
        
        @property
        def x_min(self):
            """
            Element x_min ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 12
            
            """
            return _interpolate.f90wrap_splinedata1d__get__x_min(self._handle)
        
        @x_min.setter
        def x_min(self, x_min):
            _interpolate.f90wrap_splinedata1d__set__x_min(self._handle, x_min)
        
        @property
        def h_step(self):
            """
            Element h_step ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 13
            
            """
            return _interpolate.f90wrap_splinedata1d__get__h_step(self._handle)
        
        @h_step.setter
        def h_step(self, h_step):
            _interpolate.f90wrap_splinedata1d__set__h_step(self._handle, h_step)
        
        @property
        def coeff(self):
            """
            Element coeff ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 14
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata1d__array__coeff(self._handle)
            if array_handle in self._arrays:
                coeff = self._arrays[array_handle]
            else:
                coeff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata1d__array__coeff)
                self._arrays[array_handle] = coeff
            return coeff
        
        @coeff.setter
        def coeff(self, coeff):
            self.coeff[...] = coeff
        
        def __str__(self):
            ret = ['<splinedata1d>{\n']
            ret.append('    order : ')
            ret.append(repr(self.order))
            ret.append(',\n    num_points : ')
            ret.append(repr(self.num_points))
            ret.append(',\n    periodic : ')
            ret.append(repr(self.periodic))
            ret.append(',\n    x_min : ')
            ret.append(repr(self.x_min))
            ret.append(',\n    h_step : ')
            ret.append(repr(self.h_step))
            ret.append(',\n    coeff : ')
            ret.append(repr(self.coeff))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("interpolate.SplineData2D")
    class SplineData2D(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=splinedata2d)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 16-22
        
        """
        def __init__(self, handle=None):
            """
            self = Splinedata2D()
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 16-22
            
            
            Returns
            -------
            this : Splinedata2D
            	Object to be constructed
            
            
            Automatically generated constructor for splinedata2d
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _interpolate.f90wrap_splinedata2d_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Splinedata2D
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 16-22
            
            Parameters
            ----------
            this : Splinedata2D
            	Object to be destructed
            
            
            Automatically generated destructor for splinedata2d
            """
            if self._alloc:
                _interpolate.f90wrap_splinedata2d_finalise(this=self._handle)
        
        @property
        def order(self):
            """
            Element order ftype=integer  pytype=int
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata2d__array__order(self._handle)
            if array_handle in self._arrays:
                order = self._arrays[array_handle]
            else:
                order = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata2d__array__order)
                self._arrays[array_handle] = order
            return order
        
        @order.setter
        def order(self, order):
            self.order[...] = order
        
        @property
        def num_points(self):
            """
            Element num_points ftype=integer  pytype=int
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 18
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata2d__array__num_points(self._handle)
            if array_handle in self._arrays:
                num_points = self._arrays[array_handle]
            else:
                num_points = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata2d__array__num_points)
                self._arrays[array_handle] = num_points
            return num_points
        
        @num_points.setter
        def num_points(self, num_points):
            self.num_points[...] = num_points
        
        @property
        def periodic(self):
            """
            Element periodic ftype=logical pytype=bool
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 19
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata2d__array__periodic(self._handle)
            if array_handle in self._arrays:
                periodic = self._arrays[array_handle]
            else:
                periodic = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata2d__array__periodic)
                self._arrays[array_handle] = periodic
            return periodic
        
        @periodic.setter
        def periodic(self, periodic):
            self.periodic[...] = periodic
        
        @property
        def h_step(self):
            """
            Element h_step ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 20
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata2d__array__h_step(self._handle)
            if array_handle in self._arrays:
                h_step = self._arrays[array_handle]
            else:
                h_step = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata2d__array__h_step)
                self._arrays[array_handle] = h_step
            return h_step
        
        @h_step.setter
        def h_step(self, h_step):
            self.h_step[...] = h_step
        
        @property
        def x_min(self):
            """
            Element x_min ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 21
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata2d__array__x_min(self._handle)
            if array_handle in self._arrays:
                x_min = self._arrays[array_handle]
            else:
                x_min = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata2d__array__x_min)
                self._arrays[array_handle] = x_min
            return x_min
        
        @x_min.setter
        def x_min(self, x_min):
            self.x_min[...] = x_min
        
        @property
        def coeff(self):
            """
            Element coeff ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 22
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata2d__array__coeff(self._handle)
            if array_handle in self._arrays:
                coeff = self._arrays[array_handle]
            else:
                coeff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata2d__array__coeff)
                self._arrays[array_handle] = coeff
            return coeff
        
        @coeff.setter
        def coeff(self, coeff):
            self.coeff[...] = coeff
        
        def __str__(self):
            ret = ['<splinedata2d>{\n']
            ret.append('    order : ')
            ret.append(repr(self.order))
            ret.append(',\n    num_points : ')
            ret.append(repr(self.num_points))
            ret.append(',\n    periodic : ')
            ret.append(repr(self.periodic))
            ret.append(',\n    h_step : ')
            ret.append(repr(self.h_step))
            ret.append(',\n    x_min : ')
            ret.append(repr(self.x_min))
            ret.append(',\n    coeff : ')
            ret.append(repr(self.coeff))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("interpolate.SplineData3D")
    class SplineData3D(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=splinedata3d)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 24-30
        
        """
        def __init__(self, handle=None):
            """
            self = Splinedata3D()
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 24-30
            
            
            Returns
            -------
            this : Splinedata3D
            	Object to be constructed
            
            
            Automatically generated constructor for splinedata3d
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _interpolate.f90wrap_splinedata3d_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Splinedata3D
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 24-30
            
            Parameters
            ----------
            this : Splinedata3D
            	Object to be destructed
            
            
            Automatically generated destructor for splinedata3d
            """
            if self._alloc:
                _interpolate.f90wrap_splinedata3d_finalise(this=self._handle)
        
        @property
        def order(self):
            """
            Element order ftype=integer  pytype=int
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata3d__array__order(self._handle)
            if array_handle in self._arrays:
                order = self._arrays[array_handle]
            else:
                order = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata3d__array__order)
                self._arrays[array_handle] = order
            return order
        
        @order.setter
        def order(self, order):
            self.order[...] = order
        
        @property
        def num_points(self):
            """
            Element num_points ftype=integer  pytype=int
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 26
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata3d__array__num_points(self._handle)
            if array_handle in self._arrays:
                num_points = self._arrays[array_handle]
            else:
                num_points = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata3d__array__num_points)
                self._arrays[array_handle] = num_points
            return num_points
        
        @num_points.setter
        def num_points(self, num_points):
            self.num_points[...] = num_points
        
        @property
        def periodic(self):
            """
            Element periodic ftype=logical pytype=bool
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 27
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata3d__array__periodic(self._handle)
            if array_handle in self._arrays:
                periodic = self._arrays[array_handle]
            else:
                periodic = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata3d__array__periodic)
                self._arrays[array_handle] = periodic
            return periodic
        
        @periodic.setter
        def periodic(self, periodic):
            self.periodic[...] = periodic
        
        @property
        def h_step(self):
            """
            Element h_step ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 28
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata3d__array__h_step(self._handle)
            if array_handle in self._arrays:
                h_step = self._arrays[array_handle]
            else:
                h_step = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata3d__array__h_step)
                self._arrays[array_handle] = h_step
            return h_step
        
        @h_step.setter
        def h_step(self, h_step):
            self.h_step[...] = h_step
        
        @property
        def x_min(self):
            """
            Element x_min ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 29
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata3d__array__x_min(self._handle)
            if array_handle in self._arrays:
                x_min = self._arrays[array_handle]
            else:
                x_min = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata3d__array__x_min)
                self._arrays[array_handle] = x_min
            return x_min
        
        @x_min.setter
        def x_min(self, x_min):
            self.x_min[...] = x_min
        
        @property
        def coeff(self):
            """
            Element coeff ftype=real(dp) pytype=float
            
            
            Defined at /home/ert/code/libneo/build/interpolate.f90.i line 30
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _interpolate.f90wrap_splinedata3d__array__coeff(self._handle)
            if array_handle in self._arrays:
                coeff = self._arrays[array_handle]
            else:
                coeff = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _interpolate.f90wrap_splinedata3d__array__coeff)
                self._arrays[array_handle] = coeff
            return coeff
        
        @coeff.setter
        def coeff(self, coeff):
            self.coeff[...] = coeff
        
        def __str__(self):
            ret = ['<splinedata3d>{\n']
            ret.append('    order : ')
            ret.append(repr(self.order))
            ret.append(',\n    num_points : ')
            ret.append(repr(self.num_points))
            ret.append(',\n    periodic : ')
            ret.append(repr(self.periodic))
            ret.append(',\n    h_step : ')
            ret.append(repr(self.h_step))
            ret.append(',\n    x_min : ')
            ret.append(repr(self.x_min))
            ret.append(',\n    coeff : ')
            ret.append(repr(self.coeff))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def construct_splines_1d(x_min, x_max, y, order, periodic):
        """
        spl = construct_splines_1d(x_min, x_max, y, order, periodic)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 36-53
        
        Parameters
        ----------
        x_min : float
        x_max : float
        y : float array
        order : int
        periodic : bool
        
        Returns
        -------
        spl : Splinedata1D
        
        """
        spl = _interpolate.f90wrap_construct_splines_1d(x_min=x_min, x_max=x_max, y=y, \
            order=order, periodic=periodic)
        spl = f90wrap.runtime.lookup_class("interpolate.SplineData1D").from_handle(spl, \
            alloc=True)
        return spl
    
    @staticmethod
    def destroy_splines_1d(self):
        """
        destroy_splines_1d(self)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 55-57
        
        Parameters
        ----------
        spl : Splinedata1D
        
        """
        _interpolate.f90wrap_destroy_splines_1d(spl=self._handle)
    
    @staticmethod
    def evaluate_splines_1d(self, x):
        """
        y = evaluate_splines_1d(self, x)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 59-73
        
        Parameters
        ----------
        spl : Splinedata1D
        x : float
        
        Returns
        -------
        y : float
        
        """
        y = _interpolate.f90wrap_evaluate_splines_1d(spl=self._handle, x=x)
        return y
    
    @staticmethod
    def evaluate_splines_1d_der(self, x):
        """
        y, dy = evaluate_splines_1d_der(self, x)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 75-92
        
        Parameters
        ----------
        spl : Splinedata1D
        x : float
        
        Returns
        -------
        y : float
        dy : float
        
        """
        y, dy = _interpolate.f90wrap_evaluate_splines_1d_der(spl=self._handle, x=x)
        return y, dy
    
    @staticmethod
    def evaluate_splines_1d_der2(self, x):
        """
        y, dy, d2y = evaluate_splines_1d_der2(self, x)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 94-115
        
        Parameters
        ----------
        spl : Splinedata1D
        x : float
        
        Returns
        -------
        y : float
        dy : float
        d2y : float
        
        """
        y, dy, d2y = _interpolate.f90wrap_evaluate_splines_1d_der2(spl=self._handle, \
            x=x)
        return y, dy, d2y
    
    @staticmethod
    def construct_splines_2d(x_min, x_max, y, order, periodic):
        """
        spl = construct_splines_2d(x_min, x_max, y, order, periodic)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 117-158
        
        Parameters
        ----------
        x_min : float array
        x_max : float array
        y : float array
        order : int array
        periodic : bool array
        
        Returns
        -------
        spl : Splinedata2D
        
        """
        spl = _interpolate.f90wrap_construct_splines_2d(x_min=x_min, x_max=x_max, y=y, \
            order=order, periodic=periodic)
        spl = f90wrap.runtime.lookup_class("interpolate.SplineData2D").from_handle(spl, \
            alloc=True)
        return spl
    
    @staticmethod
    def evaluate_splines_2d(self, x):
        """
        y = evaluate_splines_2d(self, x)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 160-182
        
        Parameters
        ----------
        spl : Splinedata2D
        x : float array
        
        Returns
        -------
        y : float
        
        """
        y = _interpolate.f90wrap_evaluate_splines_2d(spl=self._handle, x=x)
        return y
    
    @staticmethod
    def destroy_splines_2d(self):
        """
        destroy_splines_2d(self)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 184-186
        
        Parameters
        ----------
        spl : Splinedata2D
        
        """
        _interpolate.f90wrap_destroy_splines_2d(spl=self._handle)
    
    @staticmethod
    def construct_splines_3d(x_min, x_max, y, order, periodic):
        """
        spl = construct_splines_3d(x_min, x_max, y, order, periodic)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 188-255
        
        Parameters
        ----------
        x_min : float array
        x_max : float array
        y : float array
        order : int array
        periodic : bool array
        
        Returns
        -------
        spl : Splinedata3D
        
        """
        spl = _interpolate.f90wrap_construct_splines_3d(x_min=x_min, x_max=x_max, y=y, \
            order=order, periodic=periodic)
        spl = f90wrap.runtime.lookup_class("interpolate.SplineData3D").from_handle(spl, \
            alloc=True)
        return spl
    
    @staticmethod
    def evaluate_splines_3d(self, x):
        """
        y = evaluate_splines_3d(self, x)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 257-287
        
        Parameters
        ----------
        spl : Splinedata3D
        x : float array
        
        Returns
        -------
        y : float
        
        """
        y = _interpolate.f90wrap_evaluate_splines_3d(spl=self._handle, x=x)
        return y
    
    @staticmethod
    def evaluate_splines_3d_der2(self, x, dy, d2y):
        """
        y = evaluate_splines_3d_der2(self, x, dy, d2y)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 289-391
        
        Parameters
        ----------
        spl : Splinedata3D
        x : float array
        dy : float array
        d2y : float array
        
        Returns
        -------
        y : float
        
        """
        y = _interpolate.f90wrap_evaluate_splines_3d_der2(spl=self._handle, x=x, dy=dy, \
            d2y=d2y)
        return y
    
    @staticmethod
    def destroy_splines_3d(self):
        """
        destroy_splines_3d(self)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 393-395
        
        Parameters
        ----------
        spl : Splinedata3D
        
        """
        _interpolate.f90wrap_destroy_splines_3d(spl=self._handle)
    
    @staticmethod
    def disp_3d(self):
        """
        disp_3d(self)
        
        
        Defined at /home/ert/code/libneo/build/interpolate.f90.i lines 397-405
        
        Parameters
        ----------
        spl : Splinedata3D
        
        """
        _interpolate.f90wrap_disp_3d(spl=self._handle)
    
    _dt_array_initialisers = []
    

interpolate = Interpolate()

