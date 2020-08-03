# How to run

To set up the Balance code, create a symbolic link to the KiLCA library and compile:

    ln -s /proj/plasma/soft/KiLCA-2.4.2/lib/libKiLCA_Lib_V_2.4.2_MDNO_FPGEN_POLYNOMIAL_Release_64bit.a matlab/BALANCE/template_experimental/lib/libkilca_64bit.a
    cd matlab/BALANCE/template_experimental
    make -f Balance.mk_mpi
    cd -

Next, create a symbolic link to the NEO-2 executable, e.g.,

    ln -s /temp/ulbl_p/NEO-2/NEO-2-PAR/Build-Release/neo_2.x matlab/BALANCE/neo2/TEMPLATE_DIR/neo_2.x

Modify one of the Matlab scripts in `matlab/BALANCE/scripts`, e.g., `script_balance_velscan_dTestudy.m`, so that the variables point to the relevant directories, i.e., data sources, and run in Matlab. If Matlab kicks you out before all runs are completed, modify the corresponding loop to skip the parameters for which you already have results, and start the script again.
