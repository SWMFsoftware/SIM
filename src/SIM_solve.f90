! Copyright (C) 2002 Regents of the University of Michigan, portions used with
! permission. For more information, see http://csem.engin.umich.edu/tools/swmf

module IE_ModSolve

    use IE_ModMain
    use IE_ModMethods, ONLY: calc_ddtheta, calc_ddphi

    use ModUtilities,  ONLY: CON_STOP
    use CON_planet,    ONLY: RadiusPlanet

    implicit none

    ! Coefficients for linear system
    real, allocatable :: A1(:, :), A2(:, :), A3(:, :), A4(:, :)
    real, allocatable :: B1(:, :), B2(:, :), B3(:, :), B4(:, :), B5(:, :)

    ! Sparse representations for ModLinearSolver
    integer :: nPoints, m1, m2
    real, allocatable :: b_I(:), x_I(:), d_I(:), e_I(:), f_I(:)
    real, allocatable :: e1_I(:), f1_I(:), e2_I(:), f2_I(:)

    ! Temporary arrays to store gradients
    real, allocatable :: tmp_0PC(:, :)              ! To speed up
    real, allocatable :: tmp_0HC_eps(:, :)          ! the calculation
    real, allocatable :: tmp_deriv(:, : ,:)

    SAVE

    contains
    ! =========================================================================
    subroutine solve_ie
        ! This is similar to the Ridley solver, except SIMPLE uses a
        ! heptadiagonal linear system.
        use ModLinearSolver, ONLY: gmres, bicgstab, prehepta, Uhepta, Lhepta

        ! ---------------------------------------------------------------------
        real :: Residual
        integer :: i, j, iI, nIteration, iError
        logical :: DoTestMe = .true.

        character(len=*), parameter:: NameSub = 'solve_ie'
        ! ---------------------------------------------------------------------
        call construct_linear_system

        ! Preconditioning
        ! Left side preconditioning: U^{-1}.L^{-1}.A.x = U^{-1}.L^{-1}.rhs
        ! rhs'=U^{-1}.L^{-1}.rhs
        if(UsePreconditioner) then
            call prehepta(nPoints,1,m1,m2,-0.5,d_I,e_I,f_I,e1_I,f1_I,e2_I,f2_I)
            call Lhepta(       nPoints,1,m1,m2,b_I,d_I,e_I,e1_I,e2_I)
            call Uhepta(.true.,nPoints,1,m1,m2,b_I    ,f_I,f1_I,f2_I)
        end if

        Residual = Tolerance
        if(.not.UseInitialGuess) x_I = 0.0

        select case(NameSolver)
        case('gmres')
            nIteration = MaxIteration
            call gmres(matvec_ionosphere, b_I, x_I, UseInitialGuess, &
                       nPoints, MaxIteration, Residual, 'rel', &
                       nIteration, iError, DoTestMe)
        case('bicgstab')
            nIteration = 3 * MaxIteration
            call bicgstab(matvec_ionosphere, b_I, x_I, UseInitialGuess, &
                          nPoints, Residual, 'rel', nIteration, &
                          iError, DoTestMe)
        case default
            call CON_STOP(NameSub//': unknown NameSolver='//NameSolver)
        end select

        if(iError /= 0 .and. iError /=3)then
            write(*,*)'IE_ERROR in iono_solve: solver failed !!!'
            if(iError < 0) &
                call CON_stop('IE_ERROR in iono_solve: residual did not decrease')
        end if


        ! Put solution into 2D array.
        do j=1,nPhi; do i=1,nTheta
            iI = (j - 1) * nTheta + i
            Vars_VII(Potential_,i,j) = x_I(iI)
        end do; end do

        ! Apply pole boundary condition
        Vars_VII(Potential_     ,1,:) = sum(Vars_VII(Potential_,2,:)) / nPhi
        Vars_VII(Potential_,nTheta,:) = sum(Vars_VII(Potential_,nTheta-1,:)) / nPhi

        call clean_linear_system

        write(*,*) 'Solver completed.'

    end subroutine solve_ie


    ! =========================================================================
    subroutine construct_linear_system

        ! This subroutine constructs the coefficients for the linear system.
        ! The Poisson equation can be written in the following form
        !
        ! j_R sin^2(Theta) R1^2 =   A1 * d^2/dTheta^2(Potential)
        !                         + A2 * d^2/dPhi^2(Potential)
        !                         + A3 * d/dTheta(Potential)
        !                         + A4 * d/dPhi(Potential)
        !
        ! After calculating A1, A2, A3 and A4 (which involve calculations of
        ! gradients of certain quantities), we construct the pentadiagonal
        ! matrix D and solve the linear system (Dx = b).
        ! ---------------------------------------------------------------------
        integer :: i, j, iI
        real :: sinsqrtheta
        character(len=*), parameter:: NameSub = 'construct_linear_system'
        ! ---------------------------------------------------------------------
        call init_linear_system

        ! Calculate temporary variables --------
        do j=1,nPhi; do i=1,nTheta

            tmp_0PC(i,j) = &
                Vars_VII(Sigma0_,i,j)*Vars_VII(SigmaP_,i,j)/Vars_VII(C_,i,j)

            tmp_0HC_eps(i,j) = &
                Vars_VII(Sigma0_,i,j)*Vars_VII(SigmaH_,i,j)/Vars_VII(C_,i,j) &
                    *cos(Vars_VII(Epsilon_,i,j))

        end do; end do

        ! Calculate derivatives needed for A3-A4
        call calc_ddtheta(tmp_0PC, tmp_deriv(1,:,:))
        call calc_ddphi(tmp_0HC_eps, tmp_deriv(2,:,:))
        call calc_ddtheta(&
            -tmp_0HC_eps/sin(Vars_VII(Theta_,:,:)), tmp_deriv(3,:,:))
        call calc_ddphi(A2, tmp_deriv(4,:,:))

        ! Calculate A1,A2,A3 and A4
        ! Store the coefficients first in a separate variables
        ! B1,B3,B3,B4,B5 as the sparse variables d,e,f,e1,e2,f1
        ! f2 will be modified during the preconditioning?
        ! B1 - i,j      B2 - i+1,j       B3 - i-1,j
        !               B4 - i,j+1       B5 - i,j-1


        do j=1,nPhi; do i = 1,nTheta

            sinsqrtheta = sin(Vars_VII(Theta_,i,j))**2

            A1(i,j) = tmp_0PC(i,j) * sinsqrtheta

            A2(i,j) = Vars_VII(SigmaP_,i,j) + (Vars_VII(SigmaH_,i,j)**2 * &
                sin(Vars_VII(Epsilon_,i,j))**2 / Vars_VII(C_,i,j))

            A3(i,j) = sinsqrtheta * tmp_deriv(1,i,j) + &
                sin(2*Vars_VII(Theta_,i,j)) / 2. * tmp_0PC(i,j) + &
                sin(Vars_VII(Theta_,i,j)) * tmp_deriv(2,i,j)

            A4(i,j) = sinsqrtheta * tmp_deriv(3,i,j) + tmp_deriv(4,i,j) + &
                cos(Vars_VII(Theta_,i,j)) * (-1.) * tmp_0HC_eps(i,j) * &
                sinsqrtheta


            B1(i,j) = -2.0*(A1(i,j)/dTheta**2 + A2(i,j)/dPhi**2)
            B2(i,j) = A1(i,j)/dTheta**2 + A3(i,j)/(2*dTheta)
            B3(i,j) = A1(i,j)/dTheta**2 - A3(i,j)/(2*dTheta)
            B4(i,j) = A2(i,j)/dPhi**2   + A4(i,j)/(2*dPhi)
            B5(i,j) = A2(i,j)/dPhi**2   - A4(i,j)/(2*dPhi)

        end do; end do

        ! ------------------
        ! Convert to sparse representation
        ! In reality there are five terms corresponding to each point -
        ! (i,j), (i+1,j), (i-1,j), (i,j+1), (i,j-1) which would lead one to
        ! assume that the system would be tridiagonal. But, we also need to
        ! stich the ends at phi=1,nPhi. This leads to a sparse Heptadiagonal
        ! matrix.
        ! d  - diagonal
        ! e  - subdiagonal
        ! f  - supdiagonal
        ! e1 - 1st_subdiag
        ! f1 - 1st_supdiag
        ! e2 - 2nd_subdiag
        ! f2 - 2nd_supdiag

        do j=1,nPhi; do i=1,nTheta
            iI = (j - 1) * nTheta + i
            b_I(iI) = Vars_VII(Jr_,i,j) * sin(Vars_VII(Theta_,i,j))**2 * &
                      (Radius * RadiusPlanet)**2
            x_I(iI) = Vars_VII(Potential_,i,j)
        end do; end do

        if (UsePreconditioner) then

            ! There are five types of points depending on which of d,e,f,e1,f1,
            ! e2,f2 are non-zero.

            i = 1; j = 1; iI = 1
            d_I(iI)  = B1(i,j); f_I(iI)  = B2(i,j);
            f1_I(iI) = B4(i,j); f2_I(iI) = B5(i,j)

            i = nTheta; j = nPhi; iI = nPoints
            d_I(iI)  = B1(i,j); e_I(iI)  = B3(i,j);
            e1_I(iI) = B5(i,j); e2_I(iI) = B4(i,j)

            j = 1
            do i = 2,nTheta
                iI = (j - 1) * nTheta + i
                d_I(iI) = B1(i,j)
                e_I(iI) = B3(i,j)
                f_I(iI) = B2(i,j)
                f1_I(iI) = B4(i,j)
                f2_I(iI) = B5(i,j)
            end do

            j = nPhi
            do i = 1, nTheta - 1
                iI = (j - 1) * nTheta + i
                d_I(iI) = B1(i,j)
                e_I(iI) = B3(i,j)
                f_I(iI) = B2(i,j)
                e1_I(iI) = B5(i,j)
                e2_I(iI) = B4(i,j)
            end do

            do j = 2, nPhi - 1; do i = 1, nTheta
                iI = (j - 1) * nTheta + i
                d_I(iI) = B1(i,j)
                e_I(iI) = B3(i,j)
                f_I(iI) = B2(i,j)
                e1_I(iI) = B5(i,j)
                f1_I(iI) = B4(i,j)
            end do; end do

            m1 = nTheta
            m2 = (nPhi - 1) * nTheta

        end if

    end subroutine construct_linear_system

    ! =========================================================================
    subroutine init_linear_system

        character(len=*), parameter:: NameSub = 'init_linear_system'
        ! ---------------------------------------------------------------------

        nPoints = nTheta * nPhi

        if(.not.allocated(A1)) &
            allocate(A1(1:nTheta, 1:nPhi), A2(1:nTheta, 1:nPhi), &
                     A3(1:nTheta, 1:nPhi), A4(1:nTheta, 1:nPhi), &
                     B1(1:nTheta, 1:nPhi), B2(1:nTheta, 1:nPhi), &
                     B3(1:nTheta, 1:nPhi), B4(1:nTheta, 1:nPhi), &
                     B5(1:nTheta, 1:nPhi))

        if (.not.allocated(tmp_0PC)) then
            allocate(tmp_0PC(1:nTheta, 1:nPhi))
            allocate(tmp_0HC_eps(1:nTheta, 1:nPhi))
            allocate(tmp_deriv(1:4, 1:nTheta, 1:nPhi))
        end if

        A1 = 0.; A2 = 0.; A3 = 0.; A4 = 0.;
        tmp_0PC = 0.; tmp_0HC_eps=0.; tmp_deriv=0.

        if(.not.allocated(b_I)) then
            allocate(b_I(nPoints), x_I(nPoints))
            b_I = 0.0; x_I = 0.0;
        end if

        if(.not.allocated(d_I) .and. UsePreconditioner) then
            allocate(d_I(nPoints), &
                     e_I(nPoints),  f_I(nPoints),   &
                     e1_I(nPoints), f1_I(nPoints), &
                     e2_I(nPoints), f2_I(nPoints))

            d_I = 0.0;
            e_I = 0.0;  f_I = 0.0;
            e1_I = 0.0; f1_I = 0.0;
            e2_I = 0.0; f2_I = 0.0
        end if

    end subroutine init_linear_system

    ! =========================================================================
    subroutine clean_linear_system

        character(len=*), parameter:: NameSub = 'clean_linear_system'
        ! ---------------------------------------------------------------------

        if(allocated(A1)) then
            deallocate(A1, A2, A3, A4)
            deallocate(tmp_0PC)
            deallocate(tmp_0HC_eps)
            deallocate(tmp_deriv)
        end if

        if(allocated(b_I)) deallocate(b_I, x_I)
        if(allocated(d_I)) deallocate(d_I, e_I, f_I, e1_I, f1_I, e2_I, f2_I)


    end subroutine clean_linear_system

    ! =========================================================================

    subroutine matvec_ionosphere(x_I, y_I, n)

        use ModLinearSolver, ONLY: Uhepta, Lhepta

        implicit none
        ! ---------------------------------------------------------------------
        integer, intent(in)  :: n        ! Number of unknowns.
        real,    intent(in)  :: x_I(n)   ! Vector of unknowns.
        real,    intent(out) :: y_I(n)   ! y = A.x
        integer :: i, j, iI

        character(len=*), parameter:: NameSub = 'matvec_ionosphere'

        real :: x_G(nTheta, nPhi)
        ! ---------------------------------------------------------------------

        ! Put 1D vector into 2D solution.
        iI = 0
        do j = 1,nPhi; do i = 1,nTheta
            iI = iI + 1
            x_G(i,j) = x_I(iI)
        end do; end do

        ! Pole boundary condition.
        ! x_G(1,:)        = sum(x_G(2,:)) / nPhi
        ! x_G(nTheta,:)   = sum(x_G(nTheta-1, :)) / nPhi

        do j = 2,nPhi-1; do i = 2,nTheta-1

            iI = (j - 1) * nTheta + i
            y_I(iI) = B1(i,j) * x_G(i,j)     + &
                      B2(i,j) * x_G(i+1,j)   + &
                      B3(i,j) * x_G(i-1,j)   + &
                      B4(i,j) * x_G(i,j+1)   + &
                      B5(i,j) * x_G(i,j-1)

        end do; end do

        ! Stitching the ends ....
        j = 1
        do i = 2,nTheta-1
            iI = (j - 1) * nTheta + i
            y_I(iI) = B1(i,j) * x_G(i,j)     + &
                      B2(i,j) * x_G(i+1,j)   + &
                      B3(i,j) * x_G(i-1,j)   + &
                      B4(i,j) * x_G(i,j+1)   + &
                      B5(i,j) * x_G(i,nPhi)
        end do

        j = nPhi
        do i = 2,nTheta-1
            iI = (j - 1) * nTheta + i
            y_I(iI) = B1(i,j) * x_G(i,j)     + &
                      B2(i,j) * x_G(i+1,j)   + &
                      B3(i,j) * x_G(i-1,j)   + &
                      B4(i,j) * x_G(i,1)     + &
                      B5(i,j) * x_G(i,j-1)
        end do

        ! Top and bottom
        i = 1
        do j = 2,nPhi-1
            iI = (j - 1) * nTheta + i
            y_I(iI) = B1(i,j) * x_G(i,j)     + &
                      B2(i,j) * x_G(i+1,j)   + &
                      B4(i,j) * x_G(i,j+1)   + &
                      B5(i,j) * x_G(i,j-1)
        end do

        i = nTheta
        do j = 2,nPhi-1
            iI = (j - 1) * nTheta + i
            y_I(iI) = B1(i,j) * x_G(i,j)     + &
                      B3(i,j) * x_G(i-1,j)   + &
                      B4(i,j) * x_G(i,j+1)   + &
                      B5(i,j) * x_G(i,j-1)
        end do

        ! Corner cells
        i=1; j=1; iI = (j - 1) * nTheta + i
        y_I(iI) = B1(i,j) * x_G(i,j)     + &
                  B2(i,j) * x_G(i+1,j)   + &
                  B4(i,j) * x_G(i,j+1)   + &
                  B5(i,j) * x_G(i,nPhi)

        i=1; j=nPhi; iI = (j - 1) * nTheta + i
        y_I(iI) = B1(i,j) * x_G(i,j)     + &
                  B2(i,j) * x_G(i+1,j)   + &
                  B4(i,j) * x_G(i,1)     + &
                  B5(i,j) * x_G(i,j-1)

        i=nTheta; j=1; iI = (j - 1) * nTheta + i
        y_I(iI) = B1(i,j) * x_G(i,j)     + &
                  B3(i,j) * x_G(i-1,j)   + &
                  B4(i,j) * x_G(i,j+1)   + &
                  B5(i,j) * x_G(i,nPhi)

        i=nTheta; j=nPhi; iI = (j - 1) * nTheta + i
        y_I(iI) = B1(i,j) * x_G(i,j)     + &
                  B3(i,j) * x_G(i-1,j)   + &
                  B4(i,j) * x_G(i,1)     + &
                  B5(i,j) * x_G(i,j-1)


        if(UsePreconditioner) then
            call Lhepta(       n,1,m1,m2,y_I,d_I,e_I,e1_I,e2_I)
            call Uhepta(.true.,n,1,m1,m2,y_I    ,f_I,f1_I,f2_I)
        end if

    end subroutine matvec_ionosphere
    ! =========================================================================

end module IE_ModSolve
