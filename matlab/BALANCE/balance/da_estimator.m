classdef da_estimator < handle & hdf5_output
%##########################################################################
%classdef da_estimator < handle & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is used to make an estimation for the anomalous diffusion
% coefficient. It uses equilibrium data and already preprocessed profiles.
% It loads data depending on the input of the load method: if astra files
% are found in the path, astra output is loaded for Da. This can also be
% specified to be used as P in the estimation for Da which is called 
% "combined" here. If no astra dat is found, it looks for power data in
% form of 2 files for the radiated power and the input power. If nothing is
% found, the "universal constant" is used.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% CONSTANT:
% *) DA_UNIV
% *) e
% READONLY:
% *) profpath, fluxdatapath
% *) type
% *) r, psi_pol_norm, V, S
% *) ne, Te, dTe
% *) Pinp, Prad, Peff
% *) chi_e, Q_e
% *) Da
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = da_estimator(profpath, fluxdatapath)
% *) function loadEstimation(obj, shot, time, path, combine)
% *) function export2HDF5(obj, fname, loc)
%########################################################################## 

%author:   Philipp Ulbl
%created:  30.03.2020

    properties(Constant)
        DA_UNIV = 1e4   %universal constant
        e = 1.602e-19   %electron charge in SI (to calc W to eV/s)
    end
    
    properties(SetAccess = private)
        
        profpath        %path to profile output
        fluxdatapath    %path to fluxdata
        
        type            %type of Da: const/est/combined/astra
        
        r               %small radius
        psi_pol_norm    %normalized poloidal flux
        V               %flux surface volume
        S               %flux surface area
        
        ne              %density profile
        Te              %electron temperature profile
        dTe             %radial derivative of Te
        
        Pinp = []       %input power
        Prad = []       %radiated power in plasma
        Peff = []       %effective power on electrons
        
        chi_e = []      %astra chi of electrons
        Q_e = []        %astra Q of electrons
        
        Da              %anomalous diffusion coefficient
    end
    
    methods
        function obj = da_estimator(profpath, fluxdatapath)
            %##############################################################
            %function obj = da_estimator(profpath, fluxdatapath)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the class. loads profiles and equilibrium
            % data.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % profpath     ... path of profiles
            % fluxdatapath ... path of fluxdata with equil_r_q_psi.dat
            %##############################################################
            
            if(exist(profpath, 'dir') ~= 7), error(['profilepath not found in: ', profpath]); end
            if(exist(fluxdatapath, 'dir') ~= 7), error(['fluxdatapath not found in: ', fluxdatapath]); end
            
            %set paths
            obj.profpath = profpath;
            obj.fluxdatapath = fluxdatapath;
            
            %load profiles and equi
            equilrqpsi = dlmread([obj.fluxdatapath, 'equil_r_q_psi.dat'], '', 3, 0);
            obj.r = equilrqpsi(:, 1); %equivalent radius
            obj.psi_pol_norm = equilrqpsi(:, 3)./equilrqpsi(end, 3); %normalized poloidal flux
            obj.V = equilrqpsi(:, 7); %flux surface volume
            
            %compute derivative with smoothing
            smo = smooth2level(obj.V, obj.r, 2, 5e-2);
            obj.V = smo{1};
            obj.S = smo{2};
            
            %load raw profiles
            ne_raw = load([profpath, 'n.dat']);
            Te_raw = load([profpath, 'Te.dat']);
            
            %interp on r
            obj.ne = interp1(ne_raw(:, 1), ne_raw(:, 2), obj.r, 'pchip');
            obj.Te = interp1(Te_raw(:, 1), Te_raw(:, 2), obj.r, 'pchip');
            obj.dTe = gradient(obj.Te, obj.r);
        end
        
        function loadEstimation(obj, shot, time, path, combine)
            %##############################################################
            %function loadEstimation(obj, shot, time, path, combine)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % loads input for the estimation and calculates Da. needs shot
            % number and time to search for files in path: 1) if astra
            % file is found (ASTRA_PB_AUG_#<shot>_t<time in s>) it will be 
            % priorized and the flag combined is used. If not, files 
            % <shot>_BPD_Prad.dat and <shot>_TOT_P_TOT.dat are searched
            % for. If found they are loaded and the time slices are
            % interpolated to time. If no useful file is found in path, the
            % universal constant is used (const. in this class).
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % shot    ... shot number
            % time    ... shot time in ms
            % path    ... path of data
            % combine ... flag if astra data is found to use pure astra or
            %             combined Da
            %##############################################################
            
            astraname = [path, 'ASTRA_PB_AUG_#', num2str(shot),'_t', num2str(time/1000)];
            radpowername = [path, num2str(shot),'_BPD_Prad.dat'];
            inppowername = [path, num2str(shot),'_TOT_P_TOT.dat'];
            
            %check if astraname exists otherwise go for powernames
            if(exist(astraname, 'file'))
                
                %load data
                raw = importdata(astraname);
                %transpose
                dat = raw.data';
                %delete nan
                dat = dat(~isnan(dat));
                %reshape file which has 19columns but written in rows
                dat = reshape(dat, 19, numel(dat)/19)';

                %extract poloidal rho
                rho_pol = dat(:, 3);
                %interpolate with psi=rho^2. chi from SI to cgs and factor 3/2
                obj.chi_e = interp1(rho_pol.^2, dat(:, 4), obj.psi_pol_norm, 'pchip') .* 1e4;
                %P from MW to W and to eV/s
                obj.Q_e = interp1(rho_pol.^2, dat(:, 6), obj.psi_pol_norm, 'pchip') .* 1e6 ./ obj.e;
                
                %use full astra or combined (estimation from astra Q)
                if(~combine)
                    obj.Da = obj.chi_e .* 1.5;
                    obj.type = 'astra';
                else
                    obj.Da = -obj.Q_e ./ obj.ne ./ obj.dTe ./ obj.S;
                    obj.type = 'combined';
                end
                
                %other quantities not needed
                %T_e = interp1(rho_pol.^2, dat(:, 8), obj.psi_pol_norm, 'pchip') .* 1e3;
                %n_e = interp1(rho_pol.^2, dat(:, 10), obj.psi_pol_norm, 'pchip') .* 1e13;
                %V = interp1(rho_pol.^2, dat(:, 13), obj.psi_pol_norm, 'pchip') .* 1e6;
                
            elseif(exist(radpowername, 'file') && exist(inppowername, 'file'))
                
                %load radiated power time trace
                P_rad_raw = load(radpowername);
                %interp on time and convert to eV/s
                obj.Prad = interp1(P_rad_raw(:, 1), P_rad_raw(:, 2), time/1000) / obj.e;
                
                %load input power time trace
                P_inp_raw = load(inppowername);
                %interp on time and convert to eV/s
                obj.Pinp = interp1(P_inp_raw(:, 1), P_inp_raw(:, 2), time/1000) / obj.e;
                
                %distribute evenly to ions and electrons
                obj.Peff = 0.5 .* (obj.P_inp - obj.P_rad);
                
                obj.Da = -obj.Peff ./ obj.ne ./ obj.dTe ./ obj.S;
                obj.type = 'est';
                
            %if no files can be found use universal constant
            else
                obj.Da = obj.DA_UNIV .* ones(size(obj.r));
                obj.type = 'const';
            end
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
            
            obj.writeHDF5(fname, loc, 'type', obj.type, 'string');
            obj.writeHDF5(fname, loc, 'r', 'small radius', 'cm');
            obj.writeHDF5(fname, loc, 'V', 'flux surface volume', 'cm^3');
            obj.writeHDF5(fname, loc, 'S', 'flux surface area', 'cm^2');
            
            obj.writeHDF5(fname, loc, 'ne', 'density profile', 'cm^{-3}');
            obj.writeHDF5(fname, loc, 'Te', 'electron temperature profile', 'eV');
            obj.writeHDF5(fname, loc, 'dTe', 'radial derivative of electron temperature profile', 'eV cm^{-1}');
            obj.writeHDF5(fname, loc, 'Da', 'anomalous diffusion coefficient', 'cm^2 s^{-1}');
            
            if(~isempty(obj.Pinp))
                obj.writeHDF5(fname, loc, 'Pinp', 'input power', 'eV s^{-1}');
            end
            if(~isempty(obj.Prad))
                obj.writeHDF5(fname, loc, 'Prad', 'radiated power in plasma', 'eV s^{-1}');
            end
            if(~isempty(obj.Peff))
                obj.writeHDF5(fname, loc, 'Peff', 'effective power on electrons', 'eV s^{-1}');
            end
            
            if(~isempty(obj.chi_e))
                obj.writeHDF5(fname, loc, 'chi_e', 'astra chi of electrons', 'cm^2 s^{-1}');
            end
            if(~isempty(obj.Q_e))
                obj.writeHDF5(fname, loc, 'Q_e', 'astra Q of electrons', 'eV s^{-1}');
            end
        end
    end
end