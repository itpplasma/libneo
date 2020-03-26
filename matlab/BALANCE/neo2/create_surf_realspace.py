#!/usr/bin/env python
import sys
import numpy as np
import shutil, errno
import subprocess

surf_dir_name = "s"

s_vec, reff_vec, rbeg_vec, zbeg_vec, T_vec, N_vec, kappa_vec = np.loadtxt('surfaces.dat', unpack=True)

#print "Should I create the surface directories [y,n]?"
#ans = raw_input().lower()

ans = 'y'
if ans == 'y':
    fjobs = open("jobs_list.txt", "wt")
    k = 0
    for s in s_vec:
        dir_name = surf_dir_name + ('%.9e' %(s_vec[k])).replace('e', 'd').replace('d-', 'm')
        shutil.rmtree(dir_name, ignore_errors=True)
        shutil.copytree('TEMPLATE_DIR', dir_name, symlinks=True)
        f = open(dir_name + "/neo2.in", "rt")
        neo2 = f.read()
        f.close()

        neo2 = neo2.replace('<conl_over_mfp>', ('%.9e' %(kappa_vec[k])).replace('e', 'd'))
        neo2 = neo2.replace('<boozer_s>',      ('%.9e' %(s_vec[k])).replace('e', 'd'))
        neo2 = neo2.replace('<rbeg>',         ('%.9e' %(rbeg_vec[k])).replace('e', 'd'))
        neo2 = neo2.replace('<zbeg>',         ('%.9e' %(zbeg_vec[k])).replace('e', 'd'))

        f = open(dir_name + "/neo2.in", "wt")
        f.write(neo2)
        f.close()
        
        print dir_name
        
        fjobs.write(dir_name + "\n")
        
        k = k + 1

    fjobs.close()
