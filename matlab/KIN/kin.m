classdef kin < handle
%classdef kin < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is intended to be used to create .kin files for the use in
% e.g. GPEC. This file-type represents density/temperature/we profiles as
% also used in KiLCA (there each profile has a single file, Er instead of we).
%
% Planned:
% +) add reading of existent kin files
% +) add functionalities to create the .kin directly out of
%    single profile files and vice versa.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% s, ni, ne, ti, te, we
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = kin(fname)
% *) function set(obj, s, ni, ne, ti, te, we)
% *) function write(obj)
% *) function plotprofile(obj, var)
% *) function plotsummary(obj)
%##########################################################################

%author: Philipp Ulbl
%created: 14.05.2019

    properties (Access = public)
        fname = {}; %name of file with path
        
        s       %flux surface label s -> [0, 1]
        ni      %ion density in m^-3
        ne      %electon density in m^-3
        ti      %ion temperature in eV
        te      %electron density in eV
        we      %electric rotational frequency in rad/s
    end
    
    methods (Access = public)
        
        function obj = kin(fname)
            %##############################################################
            %function obj = kin(fname)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % creates an instance of the class. 
            %##############################################################
            % fname ... path of the file with filename without extension
            %##############################################################
            
            %splits the extenson from the filename and path if there
            [p, f, ~] = fileparts(fname);
            
            %reconstruct path+name, add extension .kin
            obj.fname = [p, '/', f, '.kin'];
        end
        
        function set(obj, s, ni, ne, ti, te, we)
            %##############################################################
            %function set(obj, s, ni, ne, ti, te, we)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % fills the class with all important properties.
            %##############################################################
            % s     ... flux surface label
            % ni    ... ion density in m^-3
            % ne    ... electron density in m^-3
            % ti    ... ion temperature in eV
            % te    ... electron temperature in eV
            % we    ... electric rotational frequency in rad/s
            %##############################################################
            
            obj.s   = s;
            obj.ni  = ni;
            obj.ne  = ne;
            obj.ti  = ti;
            obj.te  = te;
            obj.we  = we;
        end
        
        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes the data in the class into a .kin file. checks for
            % completeness of data beforehand.
            %##############################################################
                 
            %check format of vectors
            obj.check();
            
            %check if file exists and delete
            if isfile(obj.fname)
                delete(obj.fname);
                warning(['existing file <', obj.fname ,'> overwritten.']);
            end

            %construct dataset out of column or row vectors
            if (size(obj.s, 1) > size(obj.s, 2))
                topr = [obj.s; obj.ni; obj.ne; obj.ti; obj.te; obj.we];
            else
                topr = [obj.s; obj.ni; obj.ne; obj.ti; obj.te; obj.we]';
            end
            
            %open file with option to append data
            fid = fopen(obj.fname, 'a');
                        
            %print header line
            head = {'psi', 'ni[m^-3]', 'ne[m^-3]', 'ti[eV]', 'te[eV]', 'we[rad/s]'};
            cellfun(@(s) fprintf(fid, '%16s ' , s), head);
            
            %print each line
            for k = 1:numel(obj.s)
                fmt = '%16.10e ';
                fprintf(fid, '\n');
                fprintf(fid, fmt, topr(k, :));
            end
            
            %close file
            fclose(fid);
        end
        
        function plotprofile(obj, var)
            %##############################################################
            %function plotprofile(obj, var)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots a single profile over the flux surface label.
            %##############################################################
            % var   ... property to plot
            %           can be: ni, ne, ti, te, we
            %##############################################################
            
            switch var
                case 'ni'
                    y = obj.ni; l = 'm^{-3}';
                case 'ne'
                    y = obj.ne; l = 'm^{-3}';
                case 'ti'
                    y = obj.ti; l = 'eV';
                case 'te'
                    y = obj.te; l = 'eV';
                case 'we'
                    y = obj.we; l = 'rad s^{-1}';
                otherwise
                    error(['property ', var, ' not found.'])
            end
            
            plot(obj.s, y, 'DisplayName', var);
            xlabel('s')
            ylabel([var, ' / ', l])
            legend
        end
        
        function plotsummary(obj)
            %##############################################################
            %function plotsummary(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % makes a 3-part summary plot: densities, temperatures and we
            % in 3-subplots.
            %##############################################################
            
            %plot both densities in 1 axis
            subplot(3, 1, 1)
            obj.plotprofile('ni')
            hold on
            obj.plotprofile('ne')
            hold off
            ylabel('n / m^{-3}')
            
            %plot both temperatures in 1 axis
            subplot(3, 1, 2)
            obj.plotprofile('ti')
            hold on
            obj.plotprofile('te')
            hold off
            ylabel('T / eV')
            
            %plot we
            subplot(3, 1, 3)
            obj.plotprofile('we')
        end
        
%         function fromsinglefiles(obj, p)
%             %here comes a function which reads single profiles and puts
%             %them together into a kin file
%         end
    end
    
    methods(Access = private)
        
        function check(obj)
            
            %construct cell-array with vectors that must be checked
            tochk = {obj.s; obj.ni; obj.ne; obj.ti; obj.te; obj.we};
            %get size of each vector
            siz = cell2mat(cellfun(@size, tochk, 'UniformOutput', false));
            
            %throw error if they are not the same in size
            if (~all(siz==siz(1, :)))
                error('Size of profiles does not match.');
            end
        end
        
    end
end