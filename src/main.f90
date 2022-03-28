! Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with
! permission. For more information, see http://csem.engin.umich.edu/tools/swmf

program SIM

  ! This is a standalone program which can be used to test the SIM
  ! ionosphere solver.

  use SIM_main
  use SIM_solve,      ONLY: solve_ie
  use SIM_methods
  use SIM_output,     ONLY: write_ie_output

  use ModPlanetConst, ONLY: init_planet_const
  use CON_planet,     ONLY: set_planet_defaults
  use CON_axes,       ONLY: init_axes

  implicit none

  !----------------------------------------------------------------------------
  call init_planet_const
  call set_planet_defaults
  call init_axes(Time_Simulation)

  call construct_test_input
  call solve_ie
  call write_ie_output

contains
  !============================================================================
  subroutine construct_test_input
    ! Construct an idealized input similar to the paper by
    ! Merkin & Lyon (2010), https://doi.org/10.1029/2010JA015461

    integer :: i, j
    real :: j0, Theta0, dTheta0, PhiNow, ThetaNow

    !--------------------------------------------------------------------------
    call init_variables
    call create_grid

    write(*,*) 'Theta range - [',Vars_VII(Theta_,1,1), &
      ',', Vars_VII(Theta_,nTheta,1), ']'
    write(*,*) 'Phi range   - [',Vars_VII(Phi_,1,1), &
      ',', Vars_VII(Phi_,1,nPhi), ']'

    ! Set constant field-aligned and pedersen conductance
    Vars_VII(SigmaP_, :, :) = 10.
    Vars_VII(Sigma0_, :, :) = 1.0e3
    call calc_epsilon_and_C

    ! Specify input
    j0      = 1e-6
    Theta0  = 22. * cPi / 180.
    dTheta0 = 12. * cPi / 180.

    do j=1,nPhi; do i=1,nTheta
      if (Theta0 <= Vars_VII(Theta_,i,j) .and. &
        Vars_VII(Theta_,i,j) <= Theta0 + dTheta0) then

        Vars_VII(Jr_,i,j) = j0 * sin(Vars_VII(Theta_,i,j)) &
          * sin(Vars_VII(Phi_,i,j))

      else if (cPi - Theta0 >= Vars_VII(Theta_,i,j) .and. &
        Vars_VII(Theta_,i,j) >= cPi - Theta0 - dTheta0) then

        Vars_VII(Jr_,i,j) = j0 * sin(Vars_VII(Theta_,i,j)) &
          * sin(Vars_VII(Phi_,i,j))
      end if

    end do; end do

  end subroutine construct_test_input
  !============================================================================

end program
!==============================================================================
