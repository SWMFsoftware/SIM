! Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with
! permission. For more information, see http://csem.engin.umich.edu/tools/swmf

module SIM_main

  use ModNumConst, ONLY: cPi
  use ModIoUnit,   ONLY: STDOUT_

  implicit none
  save

  ! Declare variables to be used in the model.

  ! Size of spherical grid. (nTheta should be nPhi/2 + 1)
  integer :: nTheta  = 91
  integer :: nPhi    = 180
  real    :: dTheta
  real    :: dPhi

  ! Ionosphere height
  real :: Radius = 1.0        ! Eventually to be read in from ModPlanetConst

  ! Ionosphere variables (indices)
  integer, parameter :: Theta_     = 1
  integer, parameter :: Phi_       = 2
  integer, parameter :: Potential_ = 3
  integer, parameter :: Jr_        = 4
  integer, parameter :: Epsilon_   = 5
  integer, parameter :: C_         = 6
  integer, parameter :: Sigma0_    = 7
  integer, parameter :: SigmaP_    = 8
  integer, parameter :: SigmaH_    = 9

  integer, parameter :: nVars = SigmaH_     ! Change as per last number ^

  real, allocatable :: Vars_VII(:, :, :)    ! Array storing IE variables

  ! Parameters for IE Krylov solver
  integer :: MaxIteration = 10000
  real    :: Tolerance = 1e-8
  logical :: DoDebug = .true.
  logical :: UsePreconditioner = .true.
  logical :: UseInitialGuess = .true.
  character(len=10):: NameSolver = 'bicgstab'

  ! Other parameters
  real :: Time_Simulation = 1.0e9                    ! Fetch from CON_physics
  character(len=10) :: IE_CoordSystem = 'SMG NORM'   ! Default

  ! Output directory
  character(len=100) :: NameOutputDir = 'output/'

  ! MPI stuff
  integer :: iComm=-1, iProc=-1, nProc=-1

  ! Input-Output related variables
  integer :: iUnitOut = STDOUT_
  character (len=6) :: StringPrefix = ''

contains
  !============================================================================
  subroutine init_variables

    !--------------------------------------------------------------------------
    if(.not.allocated(Vars_VII)) then
      allocate(Vars_VII(1:nVars, 1:nTheta, 1:nPhi))
      Vars_VII = 0.0
    end if

  end subroutine init_variables
  !============================================================================

  subroutine clean_variables

    !--------------------------------------------------------------------------
    if(allocated(Vars_VII)) deallocate(Vars_VII)

  end subroutine clean_variables
  !============================================================================

  subroutine create_grid

      ! Creates a uniform grid covering the entire spherical surface.

    integer :: i, j
    real :: ThetaNow, PhiNow

    !--------------------------------------------------------------------------
    dPhi   = 360. / nPhi  * cPi / 180
    dTheta = (179.99 - 0.01) / (nTheta - 1) * cPi / 180

    ! Creating the grid.
    PhiNow = 0.
    do j=1,nPhi
      Vars_VII(Phi_,:,j) = PhiNow
      PhiNow = PhiNow + dPhi
    end do

    ThetaNow = 0.01
    do i=1,nTheta
      Vars_VII(Theta_,i,:) = ThetaNow
      ThetaNow = ThetaNow + dTheta
    end do

  end subroutine create_grid
  !============================================================================

end module SIM_main
!==============================================================================
