classdef balanceoutput < handle
%##########################################################################
%classdef balanceoutput < handle
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
        OUTFILES = {'fort.1000','fort.5000','Brvac.dat','cond_e.dat',...
            'equipotentials.dat','equisource.dat',...
            'par_current_e.dat','par_current_i.dat',...
            'init_params.dat'}; %cell array of all output files of balance code
    end
    
    properties
        m               %modenumber m
        n               %modenumber n
        
        fort1000        %fort.1000 file -> profiles
        fort5000        %fort.5000 file -> D-coefficients, etc..
        
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
    
    methods
        function obj = balanceoutput(m, n)
            %##############################################################
            %function obj = balanceoutput(m, n)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the balanceoutput class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % m     ... modenumber m
            % n     ... modenumber n
            %##############################################################    
            
            obj.m = m;
            obj.n = n;
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
            
            %fort 1000 + 5000
            obj.fort1000 = f1000([path, obj.OUTFILES{1}]);
            obj.fort5000 = f5000([path, obj.OUTFILES{2}]);
            %obj.fort5000 = f5000([path, 'fort.5000.0.000']);
            
            %Brvac
            raw = load([path, obj.OUTFILES{3}]);
            obj.r = raw(:, 1);
            obj.Brvac = raw(:, 2);
            
            %cond_e
            raw = load([path, obj.OUTFILES{4}]);
            obj.I10 = raw(:, 2) + 1i .* raw(:, 3);
            obj.I11 = raw(:, 4) + 1i .* raw(:, 5);
            obj.I21 = raw(:, 6) + 1i .* raw(:, 7);
            obj.I31 = raw(:, 8) + 1i .* raw(:, 9);
            
            %equipotentials
            raw = load([path, obj.OUTFILES{5}]);
            obj.psi0 = raw(:, 1);
            obj.phi0 = raw(:, 2);
            obj.psi1 = raw(:, 3) + 1i .* raw(:, 4);
            obj.phi1 = raw(:, 5) + 1i .* raw(:, 6);
            obj.Br_Abs = raw(:, 7);
            obj.Br_minus_ckpEs_over_wexB_Abs = raw(:, 8);
            
            %equisource
            raw = load([path, obj.OUTFILES{6}]);
            obj.equisource = raw(:, 2:end);
            
            %parallel current electrons
            raw = load([path, obj.OUTFILES{7}]);
            obj.Jmpe = raw(:, 2) + 1i .* raw(:, 3);
            obj.Jpe  = raw(:, 4) + 1i .* raw(:, 5);
            
            %parallel current ions
            raw = load([path, obj.OUTFILES{8}]);
            obj.Jmpi = raw(:, 2) + 1i .* raw(:, 3);
            obj.Jpi  = raw(:, 4) + 1i .* raw(:, 5);
        end
    end
end
