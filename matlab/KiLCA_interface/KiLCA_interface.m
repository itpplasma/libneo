classdef KiLCA_interface < handle
%classdef KiLCA_interface
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is intended to be used to use the code KiLCA within MATLAB. 
% The only pre-requisite is a compiled version of the code as an executable
% file. This class can then be used to create a setup, generate the file
% structure the code needs, run the code and get the output of the KiLCA
% including automatic plots of the most important quantities and radial 
% profiles of those to be used in post processing.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) EXEC_PATH, BLUE_PATH, PROF_PATH
% *) antenna, background, eigmode, modes, output, zones
% READONLY:
% *) path, pathofrun, pathofmfile, pathofprofiles
% *) run_type, run_stat, run_res
% *) rdy_torun, rdy_fromrun
% *) backgrounddata, lineardata
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_interface(path, rtype)
% *) function set_antenna(obj, ra, nmod)
% *) function set_background(obj, Rtor, rpl)
% *) function set_zones(obj, r, b, m)
% *) function set_ASDEX(obj, nmodes)
% *) function write(obj)
% *) function res = run(obj)
%##########################################################################

    %author:   Philipp Ulbl
    %created:  05.02.2019
    %modified: 22.08.2019

    properties (Access = public)
        
        %These 3 paths are !important!
        
        %path to executable
        EXEC_PATH = '~/KiLCA/KiLCA-2.4.2/exe/KiLCA_Normal_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit';            
        %EXEC_PATH = '~/KiLCA/KiLCA-2.4.2/exe/KiLCA_Normal_V_2.4.2_MDNO_NC_POLYNOMIAL_Release_64bit';
        %EXEC_PATH = '~/KiLCA/KiLCA-2.4.2/exe/KiLCA_Normal_V_2.4.2_MDYES_FPGEN_POLYNOMIAL_Release_64bit';
        
        %EXEC_PATH = '~/KiLCA/KiLCA-2.4.2/exe/KiLCA_EigParam_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit';            
        %EXEC_PATH = '~/KiLCA/KiLCA-2.4.2/exe/KiLCA_EigParam_V_2.4.2_MDNO_NC_POLYNOMIAL_Release_64bit';
        %EXEC_PATH = '~/KiLCA/KiLCA-2.4.2/exe/KiLCA_EigParam_V_2.4.2_MDYES_FPGEN_POLYNOMIAL_Release_64bit';
        
        BLUE_PATH = 'blueprints/'   %path to folder that contains blueprints
        PROF_PATH = 'profiles/'     %path to folder that contains profiles
        
        %------------------------------------------------------------------
        
        antenna                     %antenna settings represented by KiLCA_antenna class
        
        background                  %background settings represented by KiLCA_background class
        
        eigmode = KiLCA_eigmode();  %eigmode settings represented by KiLCA_eigmode class
        
        modes = KiLCA_modes();      %modes settings represented by KiLCA_modes class
        
        output = KiLCA_output();    %output settings represented by KiLCA_output class
        
        zones = {};                 %zones settings represented by cell array of KiLCA_zone class

    end
    
    properties (SetAccess = private)
   
        path = [];              %path of local KiLCA data (where data is stored)
        pathofmfile = [];       %path of the script that created this object
        pathofprofiles = [];    %path of profiles (local copy from PROF_PATH)
        pathofrun = [];         %path of run (where the generated run file is located for this instance)
        
        %------------------------------------------------------------------
        
        run_type = [];          %type of run: vacuum, imhd or flre
        run_stat = [];          %Command exit status, returned as either 0 or a nonzero integer. When the command is successful, status is 0. Otherwise, status is a nonzero integer.
        run_res = [];           %result of KiLCA run (console output)
        
        %------------------------------------------------------------------
        
        rdy_torun = false;      %flag: this interface is ready to run KiLCA
        rdy_fromrun = false;    %flag: run finished without errors
        
        %------------------------------------------------------------------
        
        backgrounddata = [];    %result: backgrounddata
        lineardata = {};        %result: if run not vacuum: 1 element for each mode
        dispersiondata = {};    %result: if run not vacuum + output disp=2
        
        postprocessors = {};    %processes results for given mode
    end
    
    methods (Access = 'public')
        
        function obj = KiLCA_interface(path, rtype)
            %##############################################################
            %function obj = KiLCA_interface(path, rtype)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the KiLCA_interface class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path  ... path to save and config KiLCA runs
            % rtype ... type of run: vacuum, imhd or flre
            %##############################################################
            
            %set all relevant paths and specify run type
            obj.path = path;
            obj.pathofprofiles = [path, obj.PROF_PATH];
            obj.pathofmfile = [pwd(), '/'];
            
            if strcmp(rtype, 'vacuum') || strcmp(rtype, 'flre') || strcmp(rtype, 'imhd')
                obj.pathofrun = [path, rtype, '/'];
                obj.run_type = rtype;
            end
        end
        
        function set_antenna(obj, ra, nmod)
            %##############################################################
            %function set_antenna(obj, ra, nmod)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % initializes the antenna with the minimum amount of needed
            % parameters to run.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % ra    ... position of antenna in cm (small radius)
            % nmod  ... number of modes to calculate (< # modes in .modes)
            %##############################################################
            
            obj.antenna = KiLCA_antenna(ra, nmod);
        end
        
        function set_background(obj, Rtor, rpl)
            %##############################################################
            %function set_background(obj, Rtor, rpl)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % initializes the background with the minimum amount of needed
            % parameters to run.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % Rtor  ... big torus radius in cm
            % rpl   ... plasma radius in cm (small radius)
            %##############################################################
            
            obj.background = KiLCA_background(Rtor, rpl);
        end
        
        function set_zones(obj, r, b, m)
            %##############################################################
            %function set_zones(obj, r, b, m)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % initializes the zones of the run.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % r     ... position of zone boundaries
            % b     ... type of zone boundaries (equally long to r)
            %           (center, infinity, interface, antenna)
            % m     ... media between zone boundaries (1 element less than r)
            %           (vacuum, medium, imhd, rmhd, flre)
            %##############################################################
            
            %check if r is a vector (even no scalar!)
            if sum(size(r)>1) ~= 1
                error('r must be a vector.');
            end
            
            %check for dimensions
            if size(r) ~= size(b)
                error('size of r and b does not match');
            end
            if numel(r) ~= (numel(m)+1) %4 boundaries, 3 mediums
                error('size of r and m does not match');
            end
            
            %check if all r values are ascending
            if ~all(r(2:end) - r(1:(end-1))>0)
                error('r values are not ascending')
            end
            
            %create zones
            for k = 1:numel(m)
                obj.zones{k} = KiLCA_zone(r(k), b{k}, m{k}, r(k+1), b{k+1}); 
            end
        end
        
        function set_ASDEX(obj, nmodes)
            %##############################################################
            %function set_ASDEX(obj, nmodes)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % initializes the class for a standard run on ASDEX parameters.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % nmodes ... number of modes to calculate (default = 0)
            %##############################################################
            
            %default # modes = 0
            if nargin < 2, nmodes = 0; end

            %set antenna
            obj.set_antenna(70, nmodes);
            %set background
            obj.set_background(170.05, 67);
            %set toroidal magnetic field
            obj.background.Btor = -17563.3704;
            
            %set zones with ASDEX setup
            r = [3.0, 67.0, 70.0, 80.0];
            b = {'center', 'interface', 'antenna', 'idealwall'};
            m = {obj.run_type, 'vacuum', 'vacuum'};
            obj.set_zones(r, b, m);
        end
        
        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % creates directory structure, copies all needed files. Needs
            % to be called before .run()
            %##############################################################
            
            %create paths if not already there
            system(['mkdir -p ', obj.path]);
            system(['mkdir -p ', obj.pathofprofiles]);
            system(['mkdir -p ', obj.pathofrun]);

            %make local copy of profiles
            system(['cp ', obj.PROF_PATH, '* ', obj.pathofprofiles]);
            %delete all existing input files
            system(['rm -f ', obj.path, '*.in']);
            system(['rm -f ', obj.pathofrun, '*.in']);
            
            %create symbolic link to exe
            system(['ln -sf ', obj.EXEC_PATH, ' ', obj.pathofrun, 'run_local']);
            %create symbolic link to profiles directory
            system(['ln -sfT ', obj.pathofprofiles, ' ', obj.pathofrun, 'profiles']);
            
            %create antenna input file
            obj.antenna.write(obj.BLUE_PATH, obj.pathofrun);
            
            %create background input file and symbolic link
            obj.background.write(obj.BLUE_PATH, obj.path);
            system(['ln -sf ', [obj.path, obj.background.BLUEPRINT], ' ', obj.pathofrun, obj.background.BLUEPRINT]);
            
            %create eigmode input file and symbolic link
            obj.eigmode.write(obj.BLUE_PATH, obj.path);
            system(['ln -sf ', [obj.path, obj.eigmode.BLUEPRINT], ' ', obj.pathofrun, obj.eigmode.BLUEPRINT]);
            
            %create output input file and symbolic link
            obj.output.write(obj.BLUE_PATH, obj.path);
            system(['ln -sf ', [obj.path, obj.output.BLUEPRINT], ' ', obj.pathofrun, obj.output.BLUEPRINT]);
            
            %create modes input file and symbolic link
            obj.modes.write(obj.BLUE_PATH, obj.path);
            system(['ln -sf ', [obj.path, obj.modes.BLUEPRINT], ' ', obj.pathofrun, obj.modes.BLUEPRINT]);
            
            %create input file for each zone object
            for k = 1:numel(obj.zones)
                obj.zones{k}.write(obj.BLUE_PATH, obj.pathofrun, k);
            end

            %set flag
            obj.rdy_torun = true;
        end
        
        function stat = run(obj)
            %##############################################################
            %function res = run(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % runs KiLCA if ready to run. deletes old output files and sets
            % status variables in this class after run.
            %##############################################################
            % output:
            %--------------------------------------------------------------
            % stat   ... output status: 0 if no errors
            %##############################################################
            
            %check if ready
            if ~obj.rdy_torun
                error('interface not ready.');
            end
            
            %check if exe can be found
            if exist(obj.EXEC_PATH, 'file')~=2
                error('KiLCA executable not found.');
            end
            
            %remove old outputs if there
            system(['rm -rf ', obj.pathofrun, 'background-data/']);
            system(['rm -rf ', obj.pathofrun, 'linear-data/']);
            system(['rm -rf ', obj.pathofrun, 'dispersion-data/']);
            
            %change directory to path of run, run KiLCA exe and change back
            cd(obj.pathofrun);
            disp('%######################################################')
            disp(['KiLCA start with type = ', obj.run_type])
            [obj.run_stat, obj.run_res] = system('./run_local');
            disp(['KiLCA end with status = ', num2str(obj.run_stat)])
            disp('%######################################################')
            cd(obj.pathofmfile);
            
            %if run successful
            if obj.run_stat == 0
                obj.rdy_fromrun = true; 
                %add result: background data to class
                obj.backgrounddata = KiLCA_data_background([obj.pathofrun, '/background-data/']);
                %add result: linear data to class (1 cell element for each mode)
                for k = 1:obj.antenna.nmod
                    %directory name is complicated..
                    dirname = ['m_', num2str(obj.modes.m(k)), ...
                               '_n_', num2str(obj.modes.n(k)), ...
                               '_flab_[', num2str(obj.antenna.flab(1)), ...
                               ',', num2str(obj.antenna.flab(2)), ']/'];
                    %initialize lineardata class
                    obj.lineardata{k} = KiLCA_data_linear([obj.pathofrun, 'linear-data/', dirname]);
                end
                %add result: dispersion data to class (1 cell element for each mode)
                for k = 1:obj.antenna.nmod
                    %directory name is complicated..
                    dirname = ['m_', num2str(obj.modes.m(k)), ...
                               '_n_', num2str(obj.modes.n(k)), ...
                               '_flab_[', num2str(obj.antenna.flab(1)), ...
                               ',', num2str(obj.antenna.flab(2)), ']/'];
                    %initialize lineardata class
                    obj.dispersiondata{k} = KiLCA_data_dispersion([obj.pathofrun, 'dispersion-data/', dirname]);
                end
            end
            
            stat = obj.run_stat;
        end
        
        function post(obj, imode)
            %##############################################################
            %function res = run(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % runs KiLCA if ready to run. deletes old output files and sets
            % status variables in this class after run.
            %##############################################################
            % output:
            %--------------------------------------------------------------
            % stat   ... output status: 0 if no errors
            %##############################################################
            
            if(numel(imode) > numel(obj.lineardata))
                error('too many modes to process specified in imode.');
            end
            if(any(imode) < 1 || any(imode) > numel(obj.lineardata))
                error('imode must contain valid indices of modes.');
            end
            if(obj.rdy_fromrun == false)
                error('KiLCA must be run first to be able to postprocess.'); 
            end
            
            for k = 1:numel(imode)
                obj.postprocessors{k} = KiLCA_postprocessor(obj, imode(k));
            end
        end
    end
end
