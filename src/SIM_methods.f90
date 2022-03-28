! Copyright (C) 2002 Regents of the University of Michigan, portions used with
! permission. For more information, see http://csem.engin.umich.edu/tools/swmf

module IE_ModMethods

  use IE_ModMain !, ONLY: nTheta, nPhi, dTheta, dPhi, &
  !                       Time_Simulation, IE_CoordSystem, Radius, &
  !                       Theta_, Phi_, Sigma0_, SigmaP_, Epsilon_, C_

  implicit none

  contains
  ! ===========================================================================
  subroutine calc_ddtheta(var_II, derivative_II)

    real, intent(in)  :: var_II(:, :)
    real, intent(out) :: derivative_II(:, :)

    integer :: i, j
    real :: two_dTheta

    two_dTheta = 2. * dTheta

    ! Middle cells first (2nd order accurate central-difference)
    do j = 1,nPhi; do i = 2,nTheta-1
        derivative_II(i,j) = (var_II(i+1,j) - var_II(i-1,j)) / two_dTheta
    end do; end do

    ! Top and bottom cells later (we use one sided 2nd order finite difference)
    i = 1
    do j = 1,nPhi
        derivative_II(i,j) = &
            (-3*var_II(i,j) + 4*var_II(i+1,j) - var_II(i+2,j)) / two_dTheta
    end do

    i = nTheta
    do j = 1,nPhi
        derivative_II(i,j) = &
            ( 3*var_II(i,j) - 4*var_II(i-1,j) + var_II(i-2,j)) / two_dTheta
    end do

  end subroutine calc_ddtheta

  ! ===========================================================================
  subroutine calc_ddphi(var_II, derivative_II)

      real, intent(in)  :: var_II(:, :)
      real, intent(out) :: derivative_II(:, :)

      integer :: i, j
      real :: two_dPhi

      two_dPhi = 2. * dPhi

      ! Middle cells first
      do j = 2,nPhi-1; do i = 1,nTheta
          derivative_II(i,j) = (var_II(i,j+1) - var_II(i,j-1)) / two_dPhi
      end do; end do

      ! j=1 and j=nPhi are separated by space (stitch the periodic boundaries)
      j = 1
      do i = 1, nTheta
          derivative_II(i,j) = (var_II(i,j+1) - var_II(i,nPhi)) / two_dPhi
      end do

      j = nPhi
      do i = 1, nTheta
          derivative_II(i,j) = (var_II(i,1) - var_II(i,j-1)) / two_dPhi
      end do

  end subroutine calc_ddphi

  ! ===========================================================================
  subroutine calc_epsilon_and_C

      ! Epsilon is defined as the angle between the magnetic field and the
      ! radial vector.
      ! epsilon = arccos(Br/B))
      ! C = Sigma0 * cos^2(epsilon) + SigmaP * sin^2(epsilon)
      ! Time_Simulation is also updated in this step.
      ! -----------------------------------------------------------------------
      use CON_planet_field, ONLY: get_planet_field

      integer :: i, j
      real :: r_D(3), B_D(3)

      r_D = 0.
      B_D = 0.

      do j=1,nPhi; do i=1,nTheta

          r_D(1) = sin(Vars_VII(Theta_, i, j)) * cos(Vars_VII(Phi_, i, j)) * Radius
          r_D(2) = sin(Vars_VII(Theta_, i, j)) * sin(Vars_VII(Phi_, i, j)) * Radius
          r_D(3) = cos(Vars_VII(Theta_, i, j)) * Radius

          ! write(*, *) Time_Simulation, r_D, IE_CoordSystem, B_D
          call get_planet_field(Time_Simulation, r_D, IE_CoordSystem, B_D)

          Vars_VII(Epsilon_, i, j) = &
                acos(sum(B_D * r_D) / &
                sqrt(sum(B_D**2) * sum(r_D**2)))

          Vars_VII(C_, i, j) = &
                Vars_VII(Sigma0_, i, j) * cos(Vars_VII(Epsilon_, i, j))**2 + &
                Vars_VII(SigmaP_, i, j) * sin(Vars_VII(Epsilon_, i, j))**2

      end do; end do

  end subroutine calc_epsilon_and_C

end module IE_ModMethods
