function c = classprop2cell(obj)
%##############################################################
%function c = classprop2cell(obj)
%##############################################################
% description:
%--------------------------------------------------------------
% converts all class properties into a plain cell array
% *only those which are NOT TRANSIENT
%##############################################################
% input:
%--------------------------------------------------------------
% obj   ... class
%##############################################################
% output:
%--------------------------------------------------------------
% c     ... 1d cell array of class properties
%##############################################################

%author:   Philipp Ulbl
%created:  21.08.2019
%modified: 21.08.2019

     c = {};
     %get metaclass object of class
     prop = metaclass(obj);
     
     %propertylist of metaclass object contains names of all properties
     for k = 1:numel(prop.PropertyList)
         %skip transient properties
         if (prop.PropertyList(k).Transient == true)
             continue; 
         end
         
         %write the property by dynamic indexing .('name')
         c{end+1} = obj.(prop.PropertyList(k).Name);
     end       
end
