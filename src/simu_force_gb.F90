! simu_force_gb
!> @brief Calculate the force of gb solvation

subroutine simu_force_gb(irep, force_mp)

    use const_maxsize
    use const_physical
    use const_index
    use var_neighbor_list, only: ele_gb_list, pnl_gb_list
    use var_setp, only: inele, ingb
    use var_struct, only: xyz_mp_rep, imp2type, &
        nbd, nba, ibd2mp, iba2mp, nmp_all
    use var_replica, only: irep2grep
#ifdef MPI_PAR
    use mpiconst
#endif

    implicit none

    ! ------------------------------------------------------------------------
    integer, intent(in)    :: irep
    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

    ! ------------------------------------------------------------------------
    ! local variables
    integer :: imp, imb, ima, imn, ksta, kend, imp1, imp2
    integer :: grep
    real(PREC) :: dist1, dist2, dist6, distr, dista, alpha12
    real(PREC) :: q12, q2
    real(PREC) :: rin, rew, rev, rcdist, alphac
    real(PREC) :: p0, p1, p2, p3, p4, p5
    real(PREC) :: rimp, rimp1, rimp2, rimp12, rimp122, vimp1, vimp2
    real(PREC) :: xcp, xcc, fcp, fco, con, fcon
    real(PREC) :: alpha(MXMP)
    real(PREC) :: v21(SPACE_DIM), for(SPACE_DIM)
    real(PREC) :: dist2bd(1:nbd)
    real(PREC) :: dist2ba(1:nba)
    real(PREC) :: dist2nl(pnl_gb_list(irep)%num_pairs)
    character(CARRAY_MSG_ERROR) :: error_message

#ifdef MPI_PAR
    if (myrank == 0) then
        error_message = 'Error: MPI is not supported in GB'
        call util_error(ERROR%STOP_ALL, error_message)
    endif
#endif

    grep = irep2grep(irep)
    rew = 1/inele%diele_water
    rcdist = 1/inele%cdist(grep)
    rin = 1/ingb%diele_mol
    rev = ingb%coef
    alphac = 0.01
    p1 = ingb%para(1)
    p2 = ingb%para(2)
    p3 = ingb%para(3)
    p4 = ingb%para(4)
    p5 = ingb%para(5)

    do imp = 1, nmp_all
        rimp = ingb%rvdw(imp)
        if (imp2type(imp) == 1) then
            p0 = ingb%pzpro
        else if (imp2type(imp) == 2 .or. imp2type(imp) == 5 .or. imp2type(imp) == 8) then
            p0 = ingb%pznuc
        end if
        alpha(imp) = p0/rimp + p1/rimp/rimp
    end do

    do imb = 1, nbd
        imp1 = ibd2mp(1, imb)
        imp2 = ibd2mp(2, imb)
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
        dist2bd(imb) = dist2
        rimp1 = ingb%rvdw(imp1)
        rimp2 = ingb%rvdw(imp2)
        vimp1 = rimp1*rimp1*rimp1
        vimp2 = rimp2*rimp2*rimp2
        alpha(imp1) = alpha(imp1) + p2*vimp2/dist2/dist2
        alpha(imp2) = alpha(imp2) + p2*vimp1/dist2/dist2
    end do

    do ima = 1, nba
        imp1 = iba2mp(1, ima)
        imp2 = iba2mp(3, ima)
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
        dist2ba(ima) = dist2
        rimp1 = ingb%rvdw(imp1)
        rimp2 = ingb%rvdw(imp2)
        vimp1 = rimp1*rimp1*rimp1
        vimp2 = rimp2*rimp2*rimp2
        alpha(imp1) = alpha(imp1) + p3*vimp2/dist2/dist2
        alpha(imp2) = alpha(imp2) + p3*vimp1/dist2/dist2
    end do

    ksta = 1
    kend = pnl_gb_list(irep)%num_pairs

    do imn = ksta, kend
        imp1 = pnl_gb_list(irep)%pairs(1, imn)
        imp2 = pnl_gb_list(irep)%pairs(2, imn)
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
        dist2nl(imn) = dist2
        rimp1 = ingb%rvdw(imp1)
        rimp2 = ingb%rvdw(imp2)
        rimp12 = rimp1 + rimp2
        rimp122 = rimp12*rimp12
        vimp1 = rimp1*rimp1*rimp1
        vimp2 = rimp2*rimp2*rimp2
        xcp = p5*dist2/rimp122
        if (xcp > 1) then
            fcp = 1
        else
            fco = 1 - cos(F_PI*xcp)
            fcp = fco*fco/4
        end if
        alpha(imp1) = alpha(imp1) + p4*fcp*vimp2/dist2/dist2
        alpha(imp2) = alpha(imp2) + p4*fcp*vimp1/dist2/dist2
    end do

    do imp = 1, nmp_all
        if (alpha(imp) > alphac) then
            alpha(imp) = 1/alpha(imp)
        else
            alpha(imp) = 0
        end if
    end do

!!$omp do private(imp1,imp2,v21,rimp1,rimp2,vimp1,vimp2,dist6,rimp12,&
!!$omp& rimp122,distr,q12,q2,alpha12,xcc,dista,dist1,con,fcon)
    ! do imb=1, nbd
    !   imp1 = ibd2mp(1,imb)
    !   imp2 = ibd2mp(2,imb)
    !   dist2 = dist2bd(imb)
    !   v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
    !   distr = sqrt(dist2)
    !   q12 = ingb%charge(imp1)*ingb%charge(imp2)
    !   for(1:3) = -rev*rin*q12/distr/dist2*v21(1:3)
    !   force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
    !   force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    !   alpha12 = alpha(imp1)*alpha(imp2)
    !   if(alpha12>0) then
    !     rimp1 = ingb%rvdw(imp1)
    !     rimp2 = ingb%rvdw(imp2)
    !     vimp1 = rimp1 * rimp1 * rimp1
    !     vimp2 = rimp2 * rimp2 * rimp2
    !     dist6 = dist2*dist2*dist2
    !     xcc = exp(-dist2/4/alpha12)
    !     dista = dist2 + alpha12*xcc
    !     dist1 = sqrt(dista)
    !     con = (dist2/4-alpha(imp1)*alpha(imp2))*(alpha(imp1)*vimp2+alpha(imp2)*vimp1)/dist6
    !     fcon = q12*(1-xcc/4-xcc*p2*con)/dist1/dista
    !     for(1:3) = rev*(rin-rew)*fcon*v21(1:3)
    !     force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
    !     force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    !     q2 = ingb%charge(imp1)*ingb%charge(imp1)
    !     for(1:3) = rev*(rin-rew)*q2*vimp2*p4*fcp/dist6*v21(1:3)
    !     force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
    !     force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    !     q2 = ingb%charge(imp2)*ingb%charge(imp2)
    !     for(1:3) = rev*(rin-rew)*q2*vimp1*p4*fcp/dist6*v21(1:3)
    !     force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
    !     force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    !   end if
    ! end do
!!$omp end do
!!$omp do private(imp1,imp2,v21,rimp1,rimp2,vimp1,vimp2,dist6,rimp12,&
!!$omp& rimp122,distr,q12,q2,alpha12,xcc,dista,dist1,con,fcon)
    ! do ima=1, nba
    !   imp1 = iba2mp(1,ima)
    !   imp2 = iba2mp(2,ima)
    !   dist2 = dist2ba(ima)
    !   v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
    !   distr = sqrt(dist2)
    !   q12 = ingb%charge(imp1)*ingb%charge(imp2)
    !   for(1:3) = -rev*rin*q12/distr/dist2*v21(1:3)
    !   force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
    !   force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    !   alpha12 = alpha(imp1)*alpha(imp2)
    !   if(alpha12>0) then
    !     rimp1 = ingb%rvdw(imp1)
    !     rimp2 = ingb%rvdw(imp2)
    !     vimp1 = rimp1 * rimp1 * rimp1
    !     vimp2 = rimp2 * rimp2 * rimp2
    !     dist6 = dist2*dist2*dist2
    !     xcc = exp(-dist2/4/alpha12)
    !     dista = dist2 + alpha12*xcc
    !     dist1 = sqrt(dista)
    !     con = (dist2/4-alpha(imp1)*alpha(imp2))*(alpha(imp1)*vimp2+alpha(imp2)*vimp1)/dist6
    !     fcon = q12*(1-xcc/4-xcc*p3*con)/dist1/dista
    !     for(1:3) = rev*(rin-rew)*fcon*v21(1:3)
    !     force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
    !     force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    !     q2 = ingb%charge(imp1)*ingb%charge(imp1)
    !     for(1:3) = rev*(rin-rew)*q2*vimp2*p4*fcp/dist6*v21(1:3)
    !     force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
    !     force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    !    q2 = ingb%charge(imp2)*ingb%charge(imp2)
    !     for(1:3) = rev*(rin-rew)*q2*vimp1*p4*fcp/dist6*v21(1:3)
    !     force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
    !     force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    !   end if
    ! end do
!!$omp end do
    do imn = ksta, kend
        imp1 = pnl_gb_list(irep)%pairs(1, imn)
        imp2 = pnl_gb_list(irep)%pairs(2, imn)
        alpha12 = alpha(imp1)*alpha(imp2)
        if (alpha12 > 0) then
            dist2 = dist2nl(imn)
            v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
            q12 = ingb%charge(imp1)*ingb%charge(imp2)
            rimp1 = ingb%rvdw(imp1)
            rimp2 = ingb%rvdw(imp2)
            vimp1 = rimp1*rimp1*rimp1
            vimp2 = rimp2*rimp2*rimp2
            dist6 = dist2*dist2*dist2
            rimp12 = rimp1 + rimp2
            rimp122 = rimp12*rimp12
            xcp = p5*dist2/rimp122
            if (xcp > 1) then
                fcp = 1
            else
                fco = 1 - cos(F_PI*xcp)
                fcp = (fco*fco - 2*F_PI*xcp*fco*sin(F_PI*xcp))/4
            end if
            xcc = exp(-dist2/4/alpha12)
            dista = dist2 + alpha12*xcc
            dist1 = sqrt(dista)
            con = (dist2/4 - alpha(imp1)*alpha(imp2))*(alpha(imp1)*vimp2 + alpha(imp2)*vimp1)/dist6
            fcon = q12*(1 - xcc/4 - xcc*p4*fcp*con)/dist1/dista
            for(1:3) = rev*(rin - exp(-dist1*rcdist)*rew)*fcon*v21(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
            fcon = q12*(1 - xcc/4 - xcc*p4*fcp*con)/dista
            for(1:3) = -rev*exp(-dist1*rcdist)*rew*rcdist*fcon*v21(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
            q2 = ingb%charge(imp1)*ingb%charge(imp1)
            con = q2*vimp2*p4*fcp/dist6
            for(1:3) = rev*(rin - exp(-alpha(imp1)*rcdist)*rew)*con*v21(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
            for(1:3) = -rev*exp(-alpha(imp1)*rcdist)*rew*rcdist*alpha(imp1)*con*v21(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
            q2 = ingb%charge(imp2)*ingb%charge(imp2)
            con = q2*vimp1*p4*fcp/dist6
            for(1:3) = rev*(rin - exp(-alpha(imp2)*rcdist)*rew)*con*v21(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
            for(1:3) = -rev*exp(-alpha(imp2)*rcdist)*rew*rcdist*alpha(imp2)*con*v21(1:3)
            force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
            force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
        end if
    end do

    ksta = 1
    kend = ele_gb_list(irep)%num_pairs

    do imn = ksta, kend
        imp1 = ele_gb_list(irep)%pairs(1, imn)
        imp2 = ele_gb_list(irep)%pairs(2, imn)
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
        distr = sqrt(dist2)
        q12 = ingb%charge(imp1)*ingb%charge(imp2)
        for(1:3) = -rev*rin*q12/distr/dist2*v21(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) - for(1:3)
    end do

end subroutine simu_force_gb
