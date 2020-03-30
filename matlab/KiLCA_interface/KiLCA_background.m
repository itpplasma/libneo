classdef KiLCA_background < handle & blueprint & hdf5_output
%classdef KiLCA_background < handle & blueprint & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% class containing the information of the background in KiLCA.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) Rtor, rpl, Btor
% *) profpath, flag_recalc, flag_back, splinedeg, vgalsys, vscale, mi, ce, ci
% *) flag_deb
% READONLY:
% *) INDICES, BLUEPRINT, READY
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_background(R, r)
% *) function c = plain(obj)
%##########################################################################
% comments:
%--------------------------------------------------------------------------
% if plasma profiles are steep at r=0 (their derivatives are not zero) then dPhi0 is not zero at r=0
% what leads to a singular oscillaltions in the f0 moments.
% 
% For the low frequency (< 1kHz) computations it is advantageous to use a moving frame.
% The wave frequency in a moving frame is om' = om - kz*V_gal_sys, kz=n/rtor, n - toroidal number.
% Usually the frequency is fixed in the lab frame (flab = om/2/pi),
% so we recalculate the frequency om' in the moving frame in the code for each n harmonic.
%##########################################################################

%author:   Philipp Ulbl
%created:  08.08.2019
%modified: 26.03.2020

    properties (Transient, SetAccess = 'protected')
        INDICES = [2:4, 7:15, 18];      %indices of parameters in blueprint files
        BLUEPRINT = 'background.in';    %name of blueprint file
        SEP = '#'
    end
    
    properties
        Rtor            %big torus radius (cm) of the machine: default=none
        rpl             %plasma radius (cm)
        Btor = 2e4     	%toroidal magnetic field (G) at the center: default=2.0e4
        profpath = './profiles/'	%path to input background profiles: default=./profiles/
        flag_recalc = 1             %1 - if background must be recalculated (7 input profiles are needed), 0 - otherwise (11 input profiles are needed), -1 - from interface
        flag_back = 'f' %flag for background ('f'-full, 'w'-wkb, 'h'-hom): default=f
        splinedeg = 9   %splines degree: >= NC + 2N+1, where N - order of flr expansion, NC - spl degree for C matrices, must be odd!!!
        vgalsys = -1e9%V_gal_sys is a velocity (cm/c) of a moving frame: default=-0.5e7
        vscale = 1      %V_scale: scale factor for the Vz velocity profile: Vz = V_scale*Vz - V_gal_sys
        mi = 2          %m_i: ions mass in units of proton mass: default=2
        ce = 1          %collisions coefficient for electrons: default=1.0
        ci = 1          %collisions coefficient for ions: default=1.0
        flag_deb = 0    %flag for debugging mode (additional checks are performed in the code): default=1
    end
    
    methods
        function obj = KiLCA_background(R, r)
            %##############################################################
            %function obj = KiLCA_background(R, r)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the background class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % R     ... torus radius (large)
            % r     ... plasma radius (small)
            %##############################################################
            obj.Rtor = R;
            obj.rpl = r;
            
            obj.READY = true;
        end
        
        function export2HDF5(obj, fname, loc)
            %##############################################################
            %function export2HDF5(obj, fname, loc)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % exports most important content of this class to hdf5file.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % fname  ... name of hdf5 file with path
            % loc    ... location of this sub-hierarchy in hdf5tree
            %##############################################################
            
            obj.writeHDF5(fname, loc, 'Rtor', 'big torus radius', 'cm');
            obj.writeHDF5(fname, loc, 'rpl', 'plasma radius (LCFS)', 'cm');
            obj.writeHDF5(fname, loc, 'Btor', 'toroidal magnetic field at the center', 'G');
            obj.writeHDF5(fname, loc, 'vgalsys', 'velocity of the moving frame', 'cm s^{-1}');
        end
    end
end