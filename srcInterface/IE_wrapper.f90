!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IE_wrapper
  ! Wrapper for the SIM Ionosphere Electrodynamics (IE) model

  use SIM_main
  use SIM_solve, ONLY: solve_ie
  use SIM_unused

  implicit none

  private ! except

  ! CON Wrapper
  public :: IE_set_param
  public :: IE_init_session
  public :: IE_run
  public :: IE_finalize

  ! Coupling with GM
  public :: IE_get_for_gm
  public :: IE_put_from_gm

  ! Not implemented with SIM but needed for compilation
  public :: IE_save_restart
  public :: IE_get_for_pw
  public :: IE_get_for_rb
  public :: IE_get_for_ps
  public :: IE_get_for_im
  public :: IE_put_from_im
  public :: IE_put_from_im_complete
  public :: IE_put_from_UA
  public :: IE_get_info_for_ua
  public :: IE_get_for_ua
  public :: IE_setnLats
  public :: IE_setnMlts

contains
  !============================================================================
  subroutine IE_set_param(CompInfo, TypeAction)

    use ModIoUnit,      ONLY: STDOUT_
    use CON_comp_info,  ONLY: put, get, CompInfoType
    use CON_comp_param, ONLY: IE_
    use CON_coupler,    ONLY: set_grid_descriptor

    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do
    character(len=*), parameter:: NameSub = 'IE_set_param'
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
      call put(CompInfo,&
        Use        =.true., &
        NameVersion='SIMPLE', &
        Version    =0.1)

    case('MPI')
      call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

    case('READ', 'CHECK')
      call read_param

    case('STDOUT')
      iUnitOut     = STDOUT_
      StringPrefix = 'IE:'

    case('FILEOUT')
      call get(CompInfo,iUnitOut=iUnitOut)
      StringPrefix=''

    case('GRID')
      call create_grid

      ! This is necessary for GM-IE coupling.
      call set_grid_descriptor(                        &
        IE_,                                        &! component index
        nDim=2,                                     &! dimensionality
        nRootBlock_D=[1,1],                       &! north+south hemispheres
        nCell_D =[nTheta, nPhi],                  &! size of node based grid
        XyzMin_D=[1.,1.],                         &! minimum indexes
        XyzMax_D=[real(nTheta), real(nPhi)],      &! maximum indexes
        TypeCoord=IE_CoordSystem,                   &! coordsystem
        Coord1_I=Vars_VII(Theta_,:,1),              &! colatitudes
        Coord2_I=Vars_VII(Phi_,1,:),                &! longitudes
        Coord3_I=[Radius],                        &! radial size in meters
        iProc_A = [iProc])                           ! processor assigment

    case default
      call CON_stop(NameSub//': IE_ERROR: SIMPLE: Invalid type/action.')
    end select

    contains
    !==========================================================================
    subroutine read_param

      use ModReadParam, ONLY: i_session_read, i_line_read, n_line_read,&
           read_line, read_command, read_var

      character (len=100) :: NameCommand

      logical :: UseStrict = .true.
      !------------------------------------------------------------------------
      select case(TypeAction)
      case('CHECK')
         write(*,*) NameSub, ': CHECK iSession =', i_session_read()
         RETURN

      case('READ')
         write(*,*) NameSub, ': READ iSession =', i_session_read(), &
              ' iLine=', i_line_read(), ' nLine=', n_line_read()

      end select

      ! Read input data via ModReadParam
      do
        if(.not. read_line()) EXIT
        if(.not.read_command(NameCommand)) CYCLE

        select case(NameCommand)
        case('#KRYLOV')
          call read_var('NameSolver',        NameSolver)
          call read_var('UsePreconditioner', UsePreconditioner)
          call read_var('UseInitialGuess',   UseInitialGuess)
          call read_var('Tolerance',         Tolerance)
          call read_var('MaxIteration',      MaxIteration)

        case('#IECOORDSYSTEM')
          call read_var('IE_CoordSystem', IE_CoordSystem)
          
        case('#DODEBUG')
          call read_var('DoDebug', DoDebug)

        case default
          write(*,'(a,i4,a)') NameSub//' IE_ERROR at line ', &
            i_line_read(), ' invalid command '//trim(NameCommand)
          if(UseStrict) call CON_stop('Correct PARAM.in!')
        end select
      end do

    end subroutine read_param
    !==========================================================================

  end subroutine IE_set_param
  !============================================================================

  subroutine IE_init_session(iSession, TimeSimulation)
    ! Initialize the IE module for session iSession

    use CON_physics,   ONLY: get_time, get_planet, get_axes
    use SIM_methods,   ONLY: calc_epsilon_and_C

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'IE_init_session'
    !--------------------------------------------------------------------------
    call init_variables
    call calc_epsilon_and_C

  end subroutine IE_init_session
  !============================================================================

  subroutine IE_finalize(TimeSimulation)

    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'IE_finalize'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_finalize
  !============================================================================

  subroutine IE_run(TimeSimulation,TimeSimulationLimit)

    real, intent(inout) :: TimeSimulation   ! current time of component

    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter:: NameSub = 'IE_run'
    !--------------------------------------------------------------------------
    call solve_ie

  end subroutine IE_run
  !============================================================================

  subroutine IE_get_for_gm(Buffer_IIV, iSize, jSize, nVar, NameVar_I, &
       tSimulation)

    ! When asked for by CON_couple_gm_ie, pass the variables from IE to GM.

    integer,          intent(in) :: iSize, jSize, nVar
    real,             intent(out):: Buffer_IIV(iSize,jSize,nVar)
    character(len=*), intent(in) :: NameVar_I(nVar)
    real,             intent(in) :: tSimulation

    integer :: iVar
    real :: tSimulationCopy

    character(len=*), parameter:: NameSub = 'IE_get_for_gm_sim'
    !--------------------------------------------------------------------------
    if(iSize /= nTheta .or. jSize /= nPhi) then
       write(*, *) NameSub//' incorrect buffer size=',iSize,jSize,&
         'nTheta=',nTheta,'nPhi=',nPhi
       call CON_stop(NameSub//' SWMF_ERROR')
    end if

    tSimulationCopy = tSimulation  ! This seems necessary because tSimulation
                                   ! has intent in.

    call IE_run(tSimulationCopy, tSimulation)

    do iVar = 1, nVar    ! nVar here is the number of buffer variables
      select case(NameVar_I(iVar))
      case('potential')
        Buffer_IIV(:,:,iVar) = Vars_VII(Potential_,:,:)
      case('sigmahall')
        Buffer_IIV(:,:,iVar) = Vars_VII(SigmaH_,:,:)
      case('sigmapedersen')
        Buffer_IIV(:,:,iVar) = Vars_VII(SigmaP_,:,:)
      case default
        call CON_stop(NameSub//': unknown NameVar='//NameVar_I(iVar))
      end select
    end do

  end subroutine IE_get_for_gm
  !============================================================================

  subroutine IE_put_from_gm(Buffer_IIV, iSize, jSize, nVar)

    ! When asked for by CON_couple_gm_ie, pass the variables from GM to IE.

    integer,          intent(in) :: iSize, jSize, nVar
    real,             intent(in) :: Buffer_IIV(iSize,jSize,nVar)

    character(len=*), parameter:: NameSub = 'IE_put_from_gm_sim'
    !--------------------------------------------------------------------------
    Vars_VII(Jr_, :, :) = Buffer_IIV(:, :, 1)

  end subroutine IE_put_from_gm
  !============================================================================

end module IE_wrapper
!==============================================================================

