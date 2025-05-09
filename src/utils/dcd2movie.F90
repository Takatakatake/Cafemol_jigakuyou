! dcd2movie
!> @brief Program for coverting dcd to movie format
! ***********************************************************************
program dcd2movie

    implicit none

    ! --------------------------------------------------------------------
    integer, parameter :: MXNUM = 100000
    integer :: ihead = 1
    integer :: n
    integer :: iarg, iargc, infile, ioutfile, ireffile
    integer :: iopen_status, input_status
    integer :: iunit, iunit2, nunit, ntitle, nblock_size
    integer :: imodel, iatom, iat, num, inum, inumber, iatref
    integer :: nset, istrt, nsavc, nstep, nver
    integer :: lunit2mp(MXNUM)
    real(4) :: w1, w2, delta
    real(4) :: x(MXNUM), y(MXNUM), z(MXNUM)
    real(8) :: tempk
    character(80) :: cinfile, coutfile, title, creffile
    character(4) :: hdr
    character(26) :: chain
    character(1)  :: chainid
    character(30) :: atmres(MXNUM)

    ! --------------------------------------------------------------------
    ! open files
    ! --------------------------------------------------------------------
    infile = 11
    ioutfile = 13
    ireffile = 14

    iarg = iargc()
    if (iarg /= 3) then
        write (*, *) 'Usage: % PROGRAM [INPUT_FILE(.pdb)] [INPUT_FILE(.dcd)] [OUTPUT_FILE(.movie)]'
        stop
    end if
    call getarg(1, creffile)
    call getarg(2, cinfile)
    call getarg(3, coutfile)
    write (*, *) creffile, cinfile, coutfile

    open (ireffile, file=creffile, status='OLD', &
          action='READ', iostat=iopen_status)
    if (iopen_status > 0) then
        write (*, *) 'Error: cannot open the .ref file'
        stop
    end if

    open (infile, file=cinfile, status='OLD', &
          action='READ', iostat=iopen_status, &
#ifdef UNFORMATTED
          form='unformatted', access='transparent')
#else
!       form='binary')
    form = 'unformatted', access = 'stream')
#endif
    if (iopen_status > 0) then
        write (*, *) 'Error: cannot open the input file'
        stop
    end if
    open (ioutfile, file=coutfile, status='REPLACE', &
          action='WRITE', iostat=iopen_status)
    if (iopen_status > 0) then
        write (*, *) 'Error: cannot open the output file'
        stop
    end if

    ! --------------------------------------------------------------------
    iatref = 1
    do
        read (ireffile, "(a30)", iostat=input_status) atmres(iatref)
        if (input_status /= 0) then
            if (input_status < 0) then
                close (ireffile)
            else
                write (*, *) 'Error: input error in pdb2crd'
            end if
            exit
        end if
        if (atmres(iatref) (1:4) == 'ATOM') then
            iatref = iatref + 1
            if (iat > MXNUM) then
                write (*, *) 'Error: should be increase MXNUM'
                stop
            end if
        end if
    end do
    iatref = iatref

    ! --------------------------------------------------------------------
    imodel = 1
    inum = 1
    iat = 0

    if (ihead == 1) then
        do
            ! ... block-size
            read (infile) nblock_size
            write (*, *) 'block size = ', nblock_size

            ! ... 'CORD' for coordinate, 'VELD' for velocity
            read (infile) hdr
            write (*, *) 'coodinate (CORD) or velocity (VELD) = ', hdr

            ! ... the number of frames
            nset = 1
            read (infile) nset
            write (*, *) 'the number of the frames =', nset

            ! ... starting step number
            read (infile) istrt

            ! ... step interval
            read (infile) nsavc

            ! ... the number of steps
            read (infile) nstep

            ! ... the number of unit
            read (infile) nunit

            read (infile) num
            read (infile) num
            read (infile) num

            ! ... the number of free atoms, where it is set to be 0.
            read (infile) num

            ! ... time-step
            read (infile) delta

            ! ... unit-cell information
            read (infile) num

            ! ... read int x eight times
            read (infile) num
            read (infile) num
            read (infile) num
            read (infile) num
            read (infile) num
            read (infile) num
            read (infile) num
            read (infile) num

            ! version if CHARMm
            read (infile) nver

            ! block-size
            read (infile) nblock_size

            ! block-size
            read (infile) nblock_size

            ! the line number of title lines
            read (infile) ntitle

            ! title text
            read (infile) title
            read (infile) title

            ! read temperature
            read (infile) title
            read (title, *) tempk

            ! read lunit2mp
            do iunit = 1, nunit
                read (infile) title
                read (title, '(i6)') lunit2mp(iunit)
            end do

            ! block-size
            read (infile) nblock_size

            ! block-size
            read (infile) nblock_size

            ! the number of atoms
            read (infile) iat

            ! block-size
            read (infile) nblock_size

            if (iat > MXNUM) then
                write (*, *) 'Error: should be increase MXNUM'
                stop
            end if
            if (iat /= iatref - 1) then
                write (*, *) 'Error: the number of particle are different'
                stop
            end if

            !  End of the dcd header part
            ihead = 2
            exit
        end do
    end if

    do
        read (infile, iostat=input_status) num
        if (input_status /= 0) then
            if (input_status < 0) then

                close (infile)
                close (ioutfile)
                write (*, *) 'PROGRAM STOP: converted ', trim(cinfile), &
                    ' to ', trim(coutfile)
            else
                write (*, *) 'Error: input error in pdb2crd'
            end if

            stop
        end if

        iat = num/4
        if ((imodel == nset) .or. (imodel <= 2)) then
            write (*, *) 'frame number = ', imodel
        else if (mod(imodel, 10**(int(log10(1.0d0*(imodel))))) == 0) then
            write (*, *) 'frame number = ', imodel
        end if

        read (infile) (x(n), n=1, iat)
        read (infile) num
        read (infile) num
        read (infile) (y(n), n=1, iat)
        read (infile) num
        read (infile) num
        read (infile) (z(n), n=1, iat)
        read (infile) num

        inumber = istrt + (imodel - 1)*nsavc

!!   write(ioutfile,'(a4, i12, f10.3)') '<<<<', inumber, tempk
        write (ioutfile, '(a43, f10.3)') &
            'REMARK    GENERATED BY CAFEMOL at tempk=   ', tempk
        write (ioutfile, '(a16, i12)') 'TITLE     step= ', inumber
!!   write(ioutfile, '(a5, i12)') 'MODEL', inumber
        write (ioutfile, '(a5)') 'MODEL'
        write (ioutfile, '(a4)') '<<<<'

        iunit = 1
        iatom = 0
        w1 = 0.0
        w2 = 0.0
        chain = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        do iatom = 1, iat

            if (iatom == 1 .or. iatom > lunit2mp(iunit)) then
                if (iatom /= 1) then
                    iunit = iunit + 1
                end if
                iunit2 = mod(iunit - 1, 26) + 1
                chainid = chain(iunit2:iunit2)
                if (iatom /= 1) then
                    write (ioutfile, '(a2)') '>>'
                end if
                write (ioutfile, '(a11, i3)') '<< protein_', iunit
            end if

            write (ioutfile, '(a4, (2xi5), a9, a2, i4, (4x3f8.3), 2f6.2)') &
                'ATOM', iatom, atmres(iatom) (12:20), chainid, iunit, &
                x(iatom), y(iatom), z(iatom), w1, w2
        end do
        write (ioutfile, '(a2)') '>>'

        write (ioutfile, '(a4)') '>>>>'
!!   write (ioutfile, '(a3)') 'END'
        write (ioutfile, '(a6)') 'ENDMDL'
        write (ioutfile, *) ''

        imodel = imodel + 1

    end do

end program dcd2movie
