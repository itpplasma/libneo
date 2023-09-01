% Make a rotation about z-axis for given groups in coildata.
%
% This function is intended to change the rotation of one or more coil
% groups in a makegrid coil file.
%
% Example usage:
%   coil_data = read_coils('coils.one_period_x16_tf_rmp2');
%   groups = 9:24;
%   coil_data = rotate_coil_group(coil_data, groups, 2*pi/(2*16))
%
% input:
% ------
% data: coil data structure as created by read_coils of matlabVMEC/XGRID.
% groups: array with the indices of the groups to be rotated.
% angle_z: angle in radians about which to rotate.
%
% output:
% -------
% data: modified coil data.
function data = rotate_coil_group(data, groups, angle_z)
  R = [+cos(angle_z), -sin(angle_z), 0;...
       +sin(angle_z), +cos(angle_z), 0;...
                   0,             0, 1];
  for k = 1:size(groups(:),1)
    L = data.vert(5,:) == groups(k);
    data.vert(1:3, L) = R*data.vert(1:3, L);
  end
end
