%##########################################################################
% script_profilerecalc.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script recalcs profiles out of raw profiles, coil files and gfile.
%##########################################################################

%author:   Philipp Ulbl
%created:  12.02.2020

addpath('~/BALANCE/balance/')

shot = 33120;
time = 5500;

profpath = '/temp/ulbl_p/PROFILES/33133_3000/EQ_new/';

neprof = '/temp/ulbl_p/AUG/SHOTS/33133/33133.3000_ne_PED_ulbl_rho_pol.dat';
Teprof = '/temp/ulbl_p/AUG/SHOTS/33133/33133.3000_Te_PED_ulbl_rho_pol.dat';
Tiprof = '/temp/ulbl_p/AUG/SHOTS/33133/33133.3000_Ti_PED_ulbl_rho_pol.dat';
vtprof = '/temp/ulbl_p/AUG/SHOTS/33133/33133.3000_vt_PED_ulbl_rho_pol.dat';

gfile = '/temp/ulbl_p/AUG/SHOTS/33133/g33133.3000_EQH';
pfile = '/itp/MooseFS/kasilov/MESH3D_33133/field.dat';
convexfile = '/proj/plasma/RMP/DATA/ASDEX/convexwall.dat';
fluxdatapath = '/temp/ulbl_p/FLUXDATA/33133/3000/';

r_min = 0.1;
r_max = 70;

mpath = pwd();
cd(profpath);

bal = Balance([profpath, 'bal/'], shot, time, '~/KiLCA_interface/');
bal.FLAG_FORCE_PROFILES = true;
bal.setCoil(pfile);
bal.setEqui(gfile, fluxdatapath);
bal.setProfiles(neprof, Teprof, Tiprof, vtprof);

cd(mpath)