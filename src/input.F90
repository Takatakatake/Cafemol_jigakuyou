#include "config.h"

subroutine input()

    use const_maxsize
    use const_index, only: RUN
    use var_inp, only: infile, i_run_mode, flg_rst
    use var_logger
    use var_file, only: input_file, restart_file
    use mod_error
#ifdef _DEBUG
    use var_inp, only: outfile
#endif
#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    type arguments_t
        character(len=128) :: input_file_name
        character(len=128) :: restart_file_name
    end type arguments_t

    ! Internal variables
    type(error_t) :: error
    type(arguments_t) args

    ! Initialize variables
    flg_rst = .false.

#ifdef LOG_DEBUG
    call debug_logger%log_debug("Subroutine input start")
#endif

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

        call parse_command_line_args(args)

        ! Open input file
        call input_file%open_file(file_name=args%input_file_name, mode="r", &
                                  new_file=.false., file_unit=infile%inp, error=error)
        call check_error(error)

        ! Open restart file
        if (len_trim(args%restart_file_name) /= 0) then
            call restart_file%open_file(file_name=args%restart_file_name, mode="r", &
                                        new_file=.false., binary=.true., file_unit=infile%rst, &
                                        error=error)
            call check_error(error)
            flg_rst = .true.
        end if

#ifdef MPI_PAR
    end if
    call MPI_Bcast(flg_rst, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

    ! reading i_run_mode for replica exchange
    call inp_runmode()

    write (log_msg_buffer, "(A,I0)") "i_run_mode = ", i_run_mode
    call logger%log_info(log_msg_buffer)

    ! reading parameters for replica exchange
    call inp_replica_para()

    ! parallelization settings for MPI
    call inp_split_mpi()

    ! making the replica set tables
    if (i_run_mode == RUN%REPLICA) then
        call inp_replica_tables()
    end if

    ! open data file and movie file
    call inp_filename()

    ! reading parameters of job_control
    call inp_job()

    ! reading unit_and_state
    call inp_unitstate()

    ! -----------------------------------------------------------------------
    ! reading energy function
    call inp_energy_func()

    ! open parameter files
    call inp_datafile()

    ! output replica variable table to .rep file
    if (i_run_mode == RUN%REPLICA) then
        call write_rep_table()
    endif

#ifdef _DEBUG
    write (6, *) 'input : infile%inp = ', infile%inp
    write (6, *) 'input : outfile%data = ', outfile%data
    write (6, *) 'input : outfile%ninfo = ', outfile%ninfo
    do i_debug = 1, MXREPLICA
        write (6, *) 'input : outfile%ts(', i_debug, ') = ', outfile%ts(i_debug)
    enddo
    do i_debug = 1, MXREPLICA
        write (6, *) 'input : outfile%movie(', i_debug, ') = ', outfile%movie(i_debug)
    enddo
    do i_debug = 1, MXREPLICA
        write (6, *) 'input : outfile%crd(', i_debug, ') = ', outfile%crd(i_debug)
    enddo
    do i_debug = 1, MXREPLICA
        write (6, *) 'input : outfile%velo(', i_debug, ') = ', outfile%velo(i_debug)
    enddo
    do i_debug = 1, MXREPLICA
        write (6, *) 'input : outfile%dcd(', i_debug, ') = ', outfile%dcd(i_debug)
    enddo
    do i_debug = 1, MXREPLICA
        write (6, *) 'input : outfile%vdcd(', i_debug, ') = ', outfile%vdcd(i_debug)
    enddo
    do i_debug = 1, MXREPLICA
        write (6, *) 'input : outfile%pdb(', i_debug, ') = ', outfile%pdb(i_debug)
    enddo
    write (6, *) 'input : outfile%rep = ', outfile%rep
    write (*, *) '#### end input'
#endif

#ifdef LOG_DEBUG
    call debug_logger%log_debug("Subroutine input end")
#endif

contains

    subroutine parse_command_line_args(args)

        use, intrinsic :: iso_fortran_env
        implicit none

        type(arguments_t), intent(out) :: args

        integer :: state
        integer :: nargs, iarg
        character(len=128) :: arg

        args%input_file_name = ""
        args%restart_file_name = ""

        state = 0

        ! Get number of arguments
        nargs = command_argument_count()

        ! Print help if no argument
        if (nargs == 0) then
            call print_help()
            stop
        end if

        ! Parse arguments
        do iarg = 1, nargs
            call get_command_argument(iarg, arg)

            if (arg(1:1) == "-") then
                call process_option(arg)
            else
                select case (state)
                case (0)
                    args%input_file_name = arg
                    state = 1
                case (1)
                    args%restart_file_name = arg
                    state = 2
                case default
                    write (error_unit, "(A)") "Error: too many arguments"
                    stop
                end select
            end if
        end do

        ! Check if input file is given
        if (state == 0) then
            call print_help()
            stop
        end if

    end subroutine parse_command_line_args

    subroutine process_option(arg)

        use, intrinsic :: iso_fortran_env
        implicit none

        character(len=128), intent(in) :: arg

        select case (arg)
        case ("--help", "-h")
            call print_help()
            stop
        case ("--version", "-v")
            call print_version()
            stop
        case default
            write (error_unit, "(A)") "Error: invalid option "//arg
            stop
        end select
    end subroutine process_option

    subroutine print_help()

        use, intrinsic :: iso_fortran_env
        implicit none

        write (output_unit, "(A)") "CafeMol"
        write (output_unit, "(A)") ""
        write (output_unit, "(A)") "USAGE:"
        write (output_unit, "(A)") "    cafemol [OPTIONS] FILE.inp [FILE.rst]"
        write (output_unit, "(A)") ""
        write (output_unit, "(A)") "OPTIONS:"
        write (output_unit, "(A)") "    -v, --version     Print version information and exit"
        write (output_unit, "(A)") "    -h, --help        Print help information and exit"
        write (output_unit, "(A)") ""
        write (output_unit, "(A)") "ARGS:"
        write (output_unit, "(A)") "    FILE.inp     An input file"
        write (output_unit, "(A)") "    FILE.rst     A restart file [optional]"

    end subroutine print_help

    subroutine print_version()

        use, intrinsic :: iso_fortran_env
        implicit none

        write (output_unit, "(A)") "cafemol "//PROJECT_VERSION

    end subroutine print_version

end subroutine input
