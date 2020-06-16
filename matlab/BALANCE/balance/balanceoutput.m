classdef balanceoutput < handle & hdf5_output
%##########################################################################
%classdef balanceoutput < handle & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is a container with all output of the balance code.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) OUTFILES
% *) m, n
% *) fort1000, fort5000
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = balanceoutput(m, n)
% *) loadOutput(obj, path)
%########################################################################## 

%author:   Philipp Ulbl
%created:  14.01.2020

    properties (Constant=true)
        OUTFILES = {'fort*','Brvac.dat','cond_e.dat',...
            'equipotentials.dat','equisource.dat',...
            'par_current_e.dat','par_current_i.dat',...
            'init_params.dat', 'timstep_evol.dat', 'amn_theta.dat',...
            'br_abs_res.dat'}; %cell array of all output files of balance code
    end
    
    properties
        path            %path of output
        
        m               %modenumber m
        n               %modenumber n
        r_res           %location of this resonant surface
        
        fort1000        %fort.1000 file -> profiles
        fort5000        %fort.5000 file -> D-coefficients, etc..
        
        time_step       %time steps in time evoluton. only there if time evolution used
        time            %total time in time evoluton. only there if time evolution used
        Br_abs_ant      %absolute radial magnetic field at antenna. only there if time evolution used
        antenna_factor  %scaling factor Cmn^2 in time evolution. only there if time evolution used
        Br_abs_res      %absolute radial magnetic field at resonance. only there if time evolution used
    end
    
    properties(SetAccess = private)
        Brvac_file
        cond_e_file
        equipotentials_file
        equisource_file
        par_current_e_file
        par_current_i_file
    end
    
    properties(Dependent)
        r               %radius
        Brvac           %vacuum radial magnetic field
        
        I10             %I factor 10 (Heyn14 eq. 50)
        I11             %I factor 11 (Heyn14 eq. 50)
        I21             %I factor 21 (Heyn14 eq. 50)
        I31             %I factor 22 (Heyn14 eq. 50)
        
        psi0            %
        phi0            %
        psi1            %
        phi1            %
        Br_Abs          %absolute radial magnetic field
        Br_minus_ckpEs_over_wexB_Abs %
        equisource      %
        Jmpe            %parallel current electrons (Heyn14 eq. 60)
        Jmpi            %parallel current ions (Heyn14 eq. 60)
        
        Jpe             %parallel current electrons (full)
        Jpi             %parallel current ions (full)
    end
    %getter
    methods
       function q = get.r(obj)
           q = obj.saveGetProp('Brvac', 1);
       end
       function q = get.Brvac(obj)
           q = obj.saveGetProp('Brvac', 2);
       end
       function q = get.I10(obj)
           q1 = obj.saveGetProp('cond_e', 2);
           q2 = obj.saveGetProp('cond_e', 3);
           q = q1 + 1i.* q2;
       end
       function q = get.I11(obj)
           q1 = obj.saveGetProp('cond_e', 4);
           q2 = obj.saveGetProp('cond_e', 5);
           q = q1 + 1i.* q2;
       end
       function q = get.I21(obj)
           q1 = obj.saveGetProp('cond_e', 6);
           q2 = obj.saveGetProp('cond_e', 7);
           q = q1 + 1i.* q2;
       end
       function q = get.I31(obj)
           q1 = obj.saveGetProp('cond_e', 8);
           q2 = obj.saveGetProp('cond_e', 9);
           q = q1 + 1i.* q2;
       end
       function q = get.psi0(obj)
           q = obj.saveGetProp('equipotentials', 2);
       end
       function q = get.phi0(obj)
           q = obj.saveGetProp('equipotentials', 3);
       end
       function q = get.psi1(obj)
           q1 = obj.saveGetProp('equipotentials', 4);
           q2 = obj.saveGetProp('equipotentials', 5);
           q = q1 + 1i.* q2;
       end
       function q = get.phi1(obj)
           q1 = obj.saveGetProp('equipotentials', 6);
           q2 = obj.saveGetProp('equipotentials', 7);
           q = q1 + 1i.* q2;
       end
       function q = get.Br_Abs(obj)
           q = obj.saveGetProp('equipotentials', 8);
       end
       function q = get.Br_minus_ckpEs_over_wexB_Abs(obj)
           q = obj.saveGetProp('equipotentials', 9);
       end
       function q = get.equisource(obj)
           q = obj.saveGetProp('equisource', 1:4);
       end
       function q = get.Jmpe(obj)
           q1 = obj.saveGetProp('par_current_e', 2);
           q2 = obj.saveGetProp('par_current_e', 3);
           q = q1 + 1i.* q2;
       end
       function q = get.Jpe(obj)
           q1 = obj.saveGetProp('par_current_e', 4);
           q2 = obj.saveGetProp('par_current_e', 5);
           q = q1 + 1i.* q2;
       end
       function q = get.Jmpi(obj)
           q1 = obj.saveGetProp('par_current_i', 2);
           q2 = obj.saveGetProp('par_current_i', 3);
           q = q1 + 1i.* q2;
       end
       function q = get.Jpi(obj)
           q1 = obj.saveGetProp('par_current_i', 4);
           q2 = obj.saveGetProp('par_current_i', 5);
           q = q1 + 1i.* q2;
       end
    end
    
    methods
        function obj = balanceoutput(m, n, r_res)
            %##############################################################
            %function obj = balanceoutput(m, n, r_res)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the balanceoutput class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % m     ... modenumber m
            % n     ... modenumber n
            % r_res ... location of this resonant surface
            %##############################################################    
            
            obj.m = m;
            obj.n = n;
            obj.r_res = r_res;
        end
        
        function loadOutput(obj, path)
            %##############################################################
            %function loadOutput(obj, path)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % loads all output of the balance code using specific classes.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path  ... mpath of the output files
            %##############################################################    
            
            obj.path = path;
            
            if(exist([path, 'timstep_evol.dat'], 'file'))
                            
                %read time evol file
                time_evol =  load([path, 'timstep_evol.dat']);
                obj.time_step = time_evol(:, 2);
                obj.time = [0; time_evol(:, end)]; %add t=0 at start (to get same length as Br_abs_ant)
                %not needed columns:
                %time_evol(:, 1) = N
                %time_evol(:, 3) = timscale_dqle
                %time_evol(:, 4) = timscal(1)
                %time_evol(:, 5) = rate_dql
                
                %read content of br_abs_res.dat file
                brabsres_raw =  load([path, 'br_abs_res.dat']);
                obj.antenna_factor = brabsres_raw(:, 3);
                obj.Br_abs_res = brabsres_raw(:, 4);
            else
                obj.time = [];
            end
            
            %if time evolution
            if(~isempty(obj.time))
                %save as cell array of objects
                obj.fort1000 = cell(1, numel(obj.time));
                obj.fort5000 = cell(1, numel(obj.time));
                
                %initialize br antenna
                obj.Br_abs_ant = nan(numel(obj.time), 1);
                %get object for each time and read br antenna
                for k=1:numel(obj.time)
                   obj.fort1000{k} = f1000([path, 'fort.', sprintf('%04d', k+999)]);
                   obj.fort5000{k} = f5000([path, 'fort.', sprintf('%04d', k+4999)]);
                   obj.Br_abs_ant(k) = obj.fort5000{k}.getLastBrAbs();
                end

            else %load only single file
                %fort 1000 + 5000
                obj.fort1000 = f1000([path, 'fort.1000']);
                obj.fort5000 = f5000([path, 'fort.5000']);
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
            
            obj.writeHDF5(fname, loc, 'time', 'total time', 's');
            obj.writeHDF5(fname, loc, 'time_step', 'time step', 's');
            obj.writeHDF5(fname, loc, 'Br_abs_ant', 'absolute radial magnetic field at antenna (unscaled)', 'G');
            obj.writeHDF5(fname, loc, 'antenna_factor', 'scaling factor Cmn^2 in time evolution.', '1');
            obj.writeHDF5(fname, loc, 'Br_abs_res', 'absolute radial magnetic field at resonance (scaled)', 'G');
            
        end
    end
    
    methods(Access=private)
        
        function loadPropFromFile(obj, file)
            %loads full file into properties
            obj.([file, '_file']) = load([obj.path, file, '.dat']);
            
        end 
        function q = saveGetProp(obj, prop, col)
            %checks if property has been loaded and loads it if not
            
            if(isempty(obj.([prop, '_file'])))
                obj.loadPropFromFile(prop);
            end
            
            q = obj.([prop, '_file'])(:, col);
        end
    end
end
