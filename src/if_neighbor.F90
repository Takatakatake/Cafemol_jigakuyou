! if_neighbor
!> @brief This is the interface for the neighborlist.

module if_neighbor
    interface

        subroutine simu_neighbor_pre(xyz_mp, ineigh_unit)
            use const_maxsize
            implicit none
            real(PREC), intent(in)  :: xyz_mp(:, :)
            integer, intent(out) :: ineigh_unit(MXUNIT, MXUNIT)
        endsubroutine simu_neighbor_pre

    endinterface
endmodule if_neighbor
