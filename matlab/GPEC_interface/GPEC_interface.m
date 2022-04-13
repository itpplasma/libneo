classdef GPEC_interface < handle
%##########################################################################
%classdef GPEC_interface < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is used as a (by now) minimalist interface to run the code
% GPEC.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) GPECHOME
% READONLY:
% *) LIB_GPEC, path_run
% *) shot, time, name
% *) coil, equil
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = GPEC_interface(runpath, shot, time, name)
% *) function write(obj)
% *) function run(obj)
%########################################################################## 

%author:   Philipp Ulbl
%created:  31.03.2020

    properties
        GPECHOME = '/proj/plasma/CODE/GPEC/';   %path to GPEC src files
    end

    properties(SetAccess='private')
        LIB_GPEC        %path to matlab lib with this file
        
        path_run        %path of run
        
        shot            %shot number
        time            %time of shot in ms
        name            %name of the run
        
        coil  = []      %inputfile for coil
        equil = []      %inputfile for equi
        dcon = []       %inputfile for dcon
        gpec = []       %inputfile for gpec
        vac = []        %inputfile for vac
    end
    
    methods(Access=public)
        function obj = GPEC_interface(runpath, shot, time, name)
            %##############################################################
            %function obj = GPEC_interface(runpath, shot, time, name)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the GPEC_interface class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % runpath   ... path of run of the balance code
            % shot      ... shot number 
            % time      ... time in ms
            % name      ... name of the run (optional)
            %##############################################################    

            %get location of this file
            s = fileparts(which(mfilename));      %filename
            obj.LIB_GPEC = [s, '/'];
            
            %add path to InputFile class
            addpath(genpath(obj.LIB_GPEC));
            
            obj.path_run = runpath;
            obj.shot = shot;
            obj.time = time;
            
            if(nargin < 4 || isempty(name))
                obj.name = 'yet another unnamed run';
            else
                obj.name = name;
            end
        end
        
        function setEqui(obj, fname)
            %##############################################################
            %function setEqui(obj, fname)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % sets equilibrium data.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % fname   ... name + path of equilibrium g-file
            %##############################################################
            
            %check if file exists
            if(~exist(fname, 'file'))
                error('equilibrium file not found.')
            end
            
            %load inputfile and write equil path
            obj.equil = InputFile([obj.LIB_GPEC, 'template/equil.in']);
            obj.equil.read();
            obj.equil.EQUIL_CONTROL.eq_filename = fname;
        end
        
        function setCoil(obj, fname)
            %##############################################################
            %function setCoil(obj, fname)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % sets coil data.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % fname   ... name + path of coilfile
            %##############################################################
            
            %check if file exists
            if(~exist(fname, 'file'))
                error('coilfile not found.')
            end
            
            %load coil file using AUG_coil class
            rawcoil = AUG_coil(fname);
            rawcoil.read();
            %convert to cell array
            current = num2cell([rawcoil.Iu; rawcoil.Il]);
            
            %load inputfile and write coil data
            obj.coil = InputFile([obj.LIB_GPEC, '/template/coil.in']);
            obj.coil.read();
            obj.coil.COIL_CONTROL.coil_cur = current;
        end
        
        function loadDGV(obj)
            %##############################################################
            %function loadDGV(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % loads dcon, gpec and vac namelists from templates.
            %##############################################################
            
            obj.dcon = InputFile([obj.LIB_GPEC, 'template/dcon.in']);
            obj.dcon.read();
            
            obj.gpec = InputFile([obj.LIB_GPEC, 'template/gpec.in']);
            obj.gpec.read();
            
            obj.vac = InputFile([obj.LIB_GPEC, 'template/vac.in']);
            obj.vac.read();
        end
        
        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % creates directory structure and copies all input files
            %##############################################################

            %load files if not done
            if(isempty(obj.gpec) || isempty(obj.dcon) || isempty(obj.vac))
                obj.loadDGV();
            end
            
            %create path
            system(['mkdir -p ', obj.path_run]);
            %delete old outputs
            system(['rm ', obj.path_run, '*.nc 2>/dev/null']);
            system(['rm ', obj.path_run, '*.bin 2>/dev/null']);
            system(['rm ', obj.path_run, '*.out 2>/dev/null']);
            system(['rm ', obj.path_run, '*.log 2>/dev/null']);
            
            %write namelists
            if(~isempty(obj.coil))
                obj.coil.write([obj.path_run, 'coil.in']);
            end
            if(~isempty(obj.equil))
                obj.equil.write([obj.path_run, 'equil.in']);
            end
            obj.dcon.write([obj.path_run, 'dcon.in']);
            obj.gpec.write([obj.path_run, 'gpec.in']);
            obj.vac.write([obj.path_run, 'vac.in']);
        end
        
        function run(obj)
            %##############################################################
            %function run(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % runs dcon and gpec
            %##############################################################

            %run dcon
            obj.run_code('dcon');
            %run gpec
            obj.run_code('gpec');
            
            %future: obj.run_code('pentrc');
        end
    end
    
    methods(Access=private)
        function run_code(obj, code)
            %function to run code by name. calcs time and outputs logfile
                        
            mpath = pwd();
            cd(obj.path_run);
            
            %create softlink to executable
            sourcepath = [obj.GPECHOME, code, '/', code];
            system(['ln -sfT ', sourcepath, ' ./', code]);
            
            %get time and create log file
            start_time = datetime;
            logfile = [code, '_', strrep(datestr(start_time), ' ', '_'), '.log'];
            fid = fopen(logfile, 'w');

            disp(['Start of ', upper(code), ' at ', datestr(start_time)])
            disp(['Shot: ', num2str(obj.shot), ', Time: ', num2str(obj.time), 'ms, Name: ', obj.name])
            %run dcon
            [stat, res] = system([obj.getConf(), '; ./', code]);
            %write to log file
            fprintf(fid, '%s\n', res);
            if(stat ~= 0)
                error(['Error in ', upper(code), '. Result = ', res, '. See log file in output directory.'])
            end
            disp(['Finished ', upper(code), ' at ', datestr(datetime)])
            disp(['Total runtime was ', string(datetime-start_time)])
            
            cd(mpath);
        end
        
        function conf = getConf(obj)
            %returns config lines
            
            conf = ['source /afs/itp.tugraz.at/opt/intel/2018.1/bin/compilervars.sh intel64 ; ',...
                    'export FC=ifort ; ',...
                    'export F77=ifort ; ',...
                    'export F90=ifort ; ',...
                    'export CC=icc ; ',...
                    'export CXX=icc ; ',...
                    'export PROJLIBS=/proj/plasma/Libs/intel ; ',...
                    'export NETCDFHOME=$PROJLIBS/NetCDF ; ',...
                    'export CFLAGS=-I$NETCDFHOME/include/ ; ',...
                    'export FFLAGS=-I$NETCDFHOME/include/ ; ',...
                    'export CXXFLAGS=-I$NETCDFHOME/include/ ; ',...
                    'export LDFLAGS=-L$NETCDFHOME/lib/ ; ',...
                    'export LD_LIBRARY_PATH=$NETCDFHOME/lib/:$LD_LIBRARY_PATH ; ',...
                    'export F90HOME=/afs/itp.tugraz.at/opt/intel/2018.1 ; ',...
                    'export GPECHOME=', obj.GPECHOME ' ; ',...
                    'ulimit -s unlimited'];
        end
    end
end