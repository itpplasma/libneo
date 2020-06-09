function bal = singlerun_balance(studyname, shot, time, time_evol, modes, copy, gpecpath)
    %##############################################################
    %function bal = singlerun_balance(studyname, shot, time, time_evol)
    %##############################################################
    % description:
    %--------------------------------------------------------------
    % this function makes a single "standard" run with the balance
    % code. this should be updated as the class evolves. this should
    % NOT be used for any special runs. 
    %##############################################################
    % input:
    %--------------------------------------------------------------
    % studyname ... name of the run
    % shot      ... shot number
    % time      ... shot time
    % time_evol ... flag to set if time evolution is desired
    %##############################################################    
      
    %author:   Philipp Ulbl
    %created:  xx.04.2020

    %default: no time evolution
    if(nargin < 4 || isempty(time_evol))
        time_evol = false;
    end
    
    m = modes(:, 1);
    n = modes(:, 2);

    runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/', num2str(shot), '_', num2str(time),'/'];
    
    %input
    gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
    filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time)];
    cfile  = [filehead,'_coil.dat'];

    dapath = '/temp/ulbl_p/DA/ASTRA/';
    
    neprof = [filehead,'_ne_PED_ULBLP_rho_pol.dat'];
    Teprof = [filehead,'_Te_PED_ULBLP_rho_pol.dat'];
    Tiprof = [filehead,'_Ti_PED_ULBLP_rho_pol.dat'];
    vtprof = [filehead,'_vt_PED_ULBLP_rho_pol.dat'];
    
    %location of field.dat
    pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
    fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

    %##########################################################################
    % RUN BALANCE WITH MATLAB CLASS
    %##########################################################################

    %BALANCE CODE
    bal = Balance(runpath, shot, time, studyname);
    bal.setModes(m, n);
    bal.setCoil(cfile);
    bal.setEqui(gfile, fluxdatapath);
    bal.setTMHDCode('GPEC', gpecpath);
    if(nargin < 6 || isempty(copy))
        bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
    else
        bal.setProfiles(neprof, Teprof, Tiprof, vtprof, copy);
    end
    bal.setKiLCA();
    bal.setDaEstimation(dapath);
    opt = balanceoptions(bal.kil_flre.pathofrun, bal.kil_vacuum.pathofrun);
    opt.stop_time_step=1e-6;
    opt.flag_run_time_evolution = time_evol;
    bal.setOptions(opt);
    bal.write();
    bal.run();

    %##########################################################################
    % EXPORT 2 HDF5 FORMAT
    %##########################################################################

    system(['mkdir -p ~/Balance_Results/', studyname, '/']);
    bal.export2HDF5(['~/Balance_Results/', studyname, '/'], studyname);
end