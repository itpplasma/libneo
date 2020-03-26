classdef balanceoptions < handle & blueprint
%##########################################################################
%classdef balanceoptions < handle & blueprint
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is an interface for the balance.in input file for the Balance
% code.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) INDICES, BLUEPRINT
% *) path_flre, path_vac
% *) Btor, Rtor, rmin, rmax, rsep
% *) npointmin, gg_fac, gg_width, gg_rres, Nstorage, tmax_fac
% *) antenna_fac, iboutype, iwrite, epsilon dperp, icoll, Z_i, am
% *) rb_cut_in, re_cut_in, rb_cut_out, re_cut_out, write_formfactors
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = balanceoptions(pf, pv)
% *) write(obj, path_from, path_to)
% *) obj2 = copy(obj)
%########################################################################## 

%author:   Philipp Ulbl
%created:  07.01.2020

    properties (Transient, SetAccess = 'protected')
        INDICES = 1:26;             %indices of parameters in blueprint files
        BLUEPRINT = 'balance.in';   %name of blueprint file
        SEP = '  '
    end

    properties
        path_flre=''                            %path to flre project
        path_vac=''                             %path to vacuum project
        Btor=-17410.9079                        %toroidal field at center
        Rtor=170.61                             %big torus radius
        rmin=3                                  %minimum small radius
        rmax=70                                 %maximum small radius
        rsep=67                                 %separatrix position
        npointmin=3000                          %
        gg_fac=[50.0, 200.0]                    %factor for gengrid
        gg_width=2.0                            %width for gengrid
        gg_rres=95.34                           %r_res for gengrid
        Nstorage=300                            %
        tmax_fac=5e1                            %
        antenna_fac=[1, 0.125, 0.0625, 0.25, 1] %
        iboutype=1                              %
        iwrite=[0, 1]                           %
        epsilon=1e-6                            %
        dperp=1e4                               %
        icoll=0                                 %
        Z_i=1.0                                 %
        am=2.0                                  %
        rb_cut_in=20.0                          %
        re_cut_in=25.0                          %
        rb_cut_out=67.5                         %
        re_cut_out=68.0                         %
        write_formfactors=[false, true]         %
    end
    
    methods
        function obj = balanceoptions(pf, pv)
            %##############################################################
            %function obj = balanceoptions(pf, pv)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the balanceoptions class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % pf    ... path to flre
            % pv    ... path to vacuum
            %##############################################################    
            
            %set filepaths
            obj.path_flre = pf;
            obj.path_vac = pv;
            
            %superclass property
            obj.READY = true;
        end
        
        function write(obj, path_from, path_to)
            %##############################################################
            %function write(obj, path_from, path_to)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % write content of class to file field_divB0.inp in path_to. A 
            % "Blueprint" is needed (existing, working file) which is
            % located in path_from.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path_from ... path of blueprint
            % path_to   ... path to write new file in
            %##############################################################    
            
            %create temporary 2nd instance, change props there and write it
            %(to not overwrite this instance)
            obj2 = obj.copy();
            obj2.path_flre = ['''', obj.path_flre, ''''];
            obj2.path_vac = ['''', obj.path_vac, ''''];
            
            %change arrays to strings for writing
            obj2.gg_fac = sprintf('%.8g ', obj.gg_fac);
            obj2.antenna_fac = sprintf('%.8g ', obj.antenna_fac);
            obj2.iwrite = sprintf('%i ', obj.iwrite);
            obj2.write_formfactors = sprintf('%i ', obj.write_formfactors);
            %transform 0 to .false. and 1 to .true.
            obj2.write_formfactors = strrep(obj2.write_formfactors, '0', '.false.');
            obj2.write_formfactors = strrep(obj2.write_formfactors, '1', '.true.');
            
            %continue with ordinary write of superclass
            write@blueprint(obj2, path_from, path_to);
        end
        
        function obj2 = copy(obj)
            %##############################################################
            %function obj2 = copy(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % creates an exact copy of the class with a new handle
            % (-> new value instance with same properties)
            %##############################################################
            % output:
            %--------------------------------------------------------------
            % obj2  ...  copy of obj
            %##############################################################    
            
            %construct new class
            obj2 = balanceoptions(obj.path_flre, obj.path_vac);
            
            %get metaclass object
        	meta = metaclass(obj);
            %iterate all properties in metaclass
            for k = 1:length(meta.PropertyList)
                %get the current property
                prop = meta.PropertyList(k);
                %skip properties that cant be set
                if strcmp(prop.SetAccess, 'private') || ...
                   strcmp(prop.SetAccess, 'none')
                	continue; 
                end
                %copy property to 2nd object if prop is a property of it
            	if isprop(obj2, prop.Name)
                    obj2.(prop.Name) = obj.(prop.Name);
            	end
            end
        end
    end
end