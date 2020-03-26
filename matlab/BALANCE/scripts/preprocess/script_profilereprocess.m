%##########################################################################
% script_profilereprocess.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script recalcs q, Er, Vth and k out of ne, Te, Ti, vt profile
% combinations. e.g.: combinations of Martins and other profiles
%##########################################################################

%author:   Philipp Ulbl
%created:  10.02.2020

addpath('~/BALANCE/balance/')

%profpath = '/temp/ulbl_p/PROFILES/33133_3000/New_ne/';
%profpath = '/temp/ulbl_p/PROFILES/33133_3000/New_Te/';
%profpath = '/temp/ulbl_p/PROFILES/33133_3000/New_ne_Te/';
%profpath = '/temp/ulbl_p/PROFILES/33133_3000/New_vt/';
%profpath = '/temp/ulbl_p/PROFILES/33133_3000/New_Ti/';
%profpath = '/temp/ulbl_p/PROFILES/33133_3000/New_vt_Ti/';
%profpath = '/temp/ulbl_p/PROFILES/33133_3000/New/';
profpath = '/temp/ulbl_p/PROFILES/33133_3000/New/';

gfile = '/proj/plasma/RMP/DATA2017/33133/3.0s/g33133.3000_ed4';
pfile = '/itp/MooseFS/kasilov/MESH3D_33133/field.dat';
convexfile = '/proj/plasma/RMP/DATA/ASDEX/convexwall.dat';
fluxdatapath = '/proj/plasma/RMP/DATA2017/33133/3.0s/PREPROCESS/';

r_min = 0.1;
r_max = 70;

%##########################################################################
% INIT PROFILES
%##########################################################################

mpath = pwd();
cd(profpath);

Te_raw = load('Te.dat');
ne_raw = load('n.dat');
Ti_raw = load('Ti.dat');
vz_raw = load('Vz.dat');

old = load('./../n_old.dat');
r = old(:, 1); %use old quantity

Te = profile_magician('Te', '', '', 1);
Te.r_out = r;
Te.y_out = interp1(Te_raw(:, 1), Te_raw(:, 2), r, 'spline');
Te.write();

ne = profile_magician('n', '', '', 1);
ne.r_out = r;
ne.y_out = interp1(ne_raw(:, 1), ne_raw(:, 2), r, 'spline');
ne.write();

Ti = profile_magician('Ti', '', '', 1);
Ti.r_out = r;
Ti.y_out = interp1(Ti_raw(:, 1), Ti_raw(:, 2), r, 'spline');
Ti.write();

vz = profile_magician('Vz', '', '', 1);
vz.r_out = r;
vz.y_out = interp1(vz_raw(:, 1), vz_raw(:, 2), r, 'spline');
vz.write();

%##########################################################################
% RECALC q, Er, vth, k
%##########################################################################

ProfilesRecalc(profpath, gfile, pfile, convexfile, fluxdatapath, ne, Te, Ti, vz, r_min, r_max);

cd(mpath);