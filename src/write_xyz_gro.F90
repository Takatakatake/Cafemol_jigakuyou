!write_xyz_gro
!> @brief Outputs coordinate information with GRO style.

subroutine write_xyz_gro(istep, velo_mp)

    use, intrinsic :: iso_fortran_env
    use const_index
    use const_maxsize
    use const_physical
    use var_inp, only: outfile, inperi
    use var_simu, only: tstep
    use var_struct, only: nunit_real, nmp_real, lunit2mp, xyz_mp_rep, &
        ires_mp, cmp2seq, cmp2atom
    use var_replica, only: n_replica_mpi, irep2grep
#ifdef MPI_PAR
    use mpiconst
#endif
    implicit none

    ! --------------------------------------------------------------------
    integer(L_INT), intent(in) :: istep
    real(PREC), intent(in) :: velo_mp(SPACE_DIM, nmp_real, n_replica_mpi)

    ! --------------------------------------------------------------------
    ! local variables
    integer :: ofile
    integer :: imp, iunit
    integer :: irep, grep

    ! --------------------------------------------------------------------
    do irep = 1, n_replica_mpi

        ! Get global replica number
        grep = irep2grep(irep)

        ! Get output file unit number
        ofile = outfile%gro(grep)

        ! Write header line
        write (ofile, "(A,F12.1)") 'System, t=', istep*tstep

        ! Write number of particles
        write (ofile, "(I8)") nmp_real

        ! Write particle positions and velocities
        do iunit = 1, nunit_real
            do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
                write (ofile, "(I5,A5,A5,I5,3(F8.4),3(F8.4))") &
                    ires_mp(imp), &
                    trim(adjustl(cmp2seq(imp))), &
                    trim(adjustl(cmp2atom(imp))), &
                    imp, &
                    xyz_mp_rep(1, imp, irep)/10, &
                    xyz_mp_rep(2, imp, irep)/10, &
                    xyz_mp_rep(3, imp, irep)/10, &
                    velo_mp(1, imp, irep)/10, &
                    velo_mp(2, imp, irep)/10, &
                    velo_mp(3, imp, irep)/10
            end do
        end do

        ! Write box size
        write (ofile, "(3(F15.4))") &
            inperi%psize(1)/10, & ! Box X size in nanometers
            inperi%psize(2)/10, & ! Box Y size in nanometers
            inperi%psize(3)/10    ! Box Z size in nanometers

        ! Write empty line
        write (ofile, "(A)") ""
    enddo

end subroutine write_xyz_gro

