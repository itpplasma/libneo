function singlerun_balance(studyname, shot, time)

    m = 4:9;
    n = 2 .* ones(size(m));

    runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/', num2str(shot), '_', num2str(time),'/'];
    
    %input
    gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
    filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time)];
    cfile  = [filehead,'_coil.dat'];

    dapath = '/temp/ulbl_p/DA/ASTRA/';
    
    neprof = [filehead,'_ne_PED_ulbl_rho_pol.dat'];
    Teprof = [filehead,'_Te_PED_ulbl_rho_pol.dat'];
    Tiprof = [filehead,'_Ti_PED_ulbl_rho_pol.dat'];
    vtprof = [filehead,'_vt_PED_ulbl_rho_pol.dat'];
    
    %location of field.dat
    pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
    fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

    %##########################################################################
    % RUN BALANCE WITH MATLAB CLASS
    %##########################################################################

    %BALANCE CODE
    bal = Balance(runpath, shot, time, studyname);
    bal.FLAG_FORCE_PROFILES = false;
    bal.setCoil(pfile, cfile);
    bal.setEqui(gfile, fluxdatapath);
    bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
    bal.setDaEstimation(dapath);
    bal.setModes(m, n);
    bal.write();
    bal.run();

    %##########################################################################
    % RESCALE
    %##########################################################################

    bal.rescaleI();
    
    gpecpath = ['/temp/ulbl_p/GPEC/', num2str(shot), '_', num2str(time), '/'];
    gpecfile = [gpecpath, 'gpec_profile_output_n2.nc'];
    
    bal.rescaleD(gpecfile, 'GPEC');

    %##########################################################################
    % EXPORT 2 HDF5 FORMAT
    %##########################################################################

    system(['mkdir -p ~/Balance_Results/', studyname, '/']);
    bal.export2HDF5(['~/Balance_Results/', studyname, '/'], studyname);
end