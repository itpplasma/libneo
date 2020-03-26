classdef AUG_coil < handle
%##########################################################################
%classdef AUG_coil < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% Class to manage coil files from AUG experiments. In these files the
% vector [Iu, Il] is writen in a row (upper and lower coil currents).
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) WINDINGS
% *) path, Il, Iu
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = AUG_coil(file)
% *) read(obj)
% *) write(obj)
% *) export2Kisslinger(obj, kfile, flip)
%########################################################################## 

%author:   Philipp Ulbl
%created:  08.01.2020

    properties
        WINDINGS = 5        %windings in coil (default=5, AUG)
            
        path                %path to coil file
        Il                  %lower coil currents
        Iu                  %upper coil currents
    end
    
    methods
        function obj = AUG_coil(file)
            %##############################################################
            %function obj = AUG_coil(file)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the AUG_coil class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % file  ... path to coil file (experiment)
            %############################################################## 
            
            obj.path = file;
        end
        
        function read(obj, flip)
            %##############################################################
            %function read(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % reads the coil file. If needed, the sequence of Iu and Il can 
            % be flipped: instead of [Iu, Il] -> [Il, Iu] is read in.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % flip  ... bool to indicate if sequence of lower and upper
            %           coils should be flipped. default = false. 
            %############################################################## 
            
            %default for flip = false
            if(nargin < 2 || isempty(flip))
                flip = false;
            end
            
            %read file
            raw = load(obj.path);
            
            %extract currents
            if(flip == false)
                obj.Iu = raw(1:end/2);          %upper
                obj.Il = raw((end/2+1):end);    %lower
            else
                obj.Il = raw(1:end/2);          %lower
                obj.Iu = raw((end/2+1):end);    %upper
            end
                
        end
        
        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes the coil file.
            %##############################################################
            
            obj.writefile(obj.path, [obj.Il, obj.Iu]);
        end
        
        function export2Kisslinger(obj, kfile, flip)
            %##############################################################
            %function export2Kisslinger(obj, kfile, flip)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % exports the coil file to Kisslinger format. For this, coil
            % currents are multiplied by WINDINGS and divided by a factor
            % of 10 which comes from the cgs conversion (without c). If
            % needed, the sequence of Iu and Il can be flipped: instead of
            % [Iu, Il] -> [Il, Iu] is written in the Kisslinger file.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % kfile ... path to new kisslinger coil file
            % flip  ... bool to indicate if sequence of lower and upper
            %           coils should be flipped. default = false. 
            %############################################################## 
            
            %default for flip = false
            if(nargin < 3 || isempty(flip))
                flip = false;
            end
            
            %flip lower/upper coils in file
            if(flip == false)
                data = [obj.Iu, obj.Il];
            else
                data = [obj.Il, obj.Iu];
            end
            
            %multiply coil currents by windings
            %divide by cgs factor 10 (without c)
            data = data .* obj.WINDINGS / 10;
            
            %write to kisslinger file
            obj.writefile(kfile, data);
        end
    end
    
    methods(Access = 'private')
        
        function writefile(obj, file, data)
            %private function to write file (interface for both write and
            %export2Kisslinger)
            
            %open file to write
            fid = fopen(file, 'w');
            %write file
            fprintf(fid, '%f ', data);
            %close file
            fclose(fid);
        end
    end
end