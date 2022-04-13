    iunit=1
    open(iunit,form='unformatted',file='formfactors.dat')
    read (iunit) nmodes_ff,nsqpsi_ff,mpol_min,mpol_max,ntor_min,ntor_max
    close(iunit)
    print *,'nmodes_ff = ',nmodes_ff
    print *,'nsqpsi_ff = ',nsqpsi_ff
    print *,'mpol_min = ',mpol_min
    print *,'mpol_max = ',mpol_max
    print *,'ntor_min = ',ntor_min
    print *,'ntor_max = ',ntor_max
    end
