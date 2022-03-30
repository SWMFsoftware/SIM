! Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with
! permission. For more information, see http://csem.engin.umich.edu/tools/swmf

module SIM_unused

	use ModUtilities, ONLY: CON_stop

	implicit none

	contains
  !============================================================================
  subroutine IE_save_restart(TimeSimulation)

    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'IE_save_restart'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_save_restart
  !============================================================================
  subroutine IE_get_for_pw(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
       tSimulation)

    integer, intent(in)           :: iSize, jSize, nVar
    real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
    character (len=*),intent(in)  :: NameHem
    character (len=*),intent(in)  :: Name_V(nVar)
    real,             intent(in)  :: tSimulation

    character(len=*), parameter:: NameSub = 'IE_get_for_pw'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_pw
  !============================================================================
  subroutine IE_get_for_rb(Buffer_IIV, iSize, jSize, nVar, Name_V, NameHem,&
       tSimulation)

    integer, intent(in)           :: iSize, jSize, nVar
    real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
    character (len=*),intent(in)  :: NameHem
    character (len=*),intent(in)  :: Name_V(nVar)
    real,             intent(in)  :: tSimulation

    character(len=*), parameter:: NameSub = 'IE_get_for_rb'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_rb
  !============================================================================
  subroutine IE_get_for_ps(Buffer_II, iSize, jSize, tSimulation)

    integer, intent(in) :: iSize, jSize
    real, intent(out)   :: Buffer_II(iSize,jSize)
    real, intent(in)    :: tSimulation

    character(len=*), parameter:: NameSub = 'IE_get_for_ps'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_ps
  !============================================================================
  subroutine IE_get_for_im(nPoint,iPointStart,Index,Weight,Buff_V,nVar)
		use CON_router,   ONLY: IndexPtrType, WeightPtrType
    integer,intent(in)            :: nPoint, iPointStart, nVar
    real,intent(out)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight

    character(len=*), parameter:: NameSub = 'IE_get_for_im'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_im
  !============================================================================
  subroutine IE_put_from_UA(Buffer_IIBV, nMLTs, nLats, nVarIn, NameVarUaIn_V)

    integer,          intent(in) :: nMlts, nLats, nVarIn
    character(len=3), intent(in) :: NameVarUaIn_V(nVarIn)
    real,             intent(in) :: Buffer_IIBV(nMlts, nLats, 2, nVarIn)

    character(len=*), parameter:: NameSub = 'IE_put_from_UA'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_put_from_UA
  !============================================================================
  subroutine IE_get_info_for_ua(nVar, NameVar_V)

    integer,          intent(out)           :: nVar
    character(len=*), intent(out), optional :: NameVar_V(:)

    character(len=*), parameter:: NameSub = 'IE_get_info_for_ua'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_info_for_ua
  !============================================================================

  subroutine IE_get_for_ua(Buffer_IIV,iSize,jSize,nVarIn,NameVar_V, &
       iBlock,tSimulation)

    integer,          intent(in)  :: iSize,jSize, nVarIn, iBlock
    real,             intent(out) :: Buffer_IIV(iSize,jSize,nVarIn)
    character (len=*),intent(in)  :: NameVar_V(nVarIn)
    real,             intent(in)  :: tSimulation

    character(len=*), parameter:: NameSub = 'IE_get_for_ua'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_get_for_ua
  !============================================================================
  subroutine IE_setnMlts(iComponent, nMLTsIn, iError)

    integer, intent(in)  :: iComponent, nMLTsIn
    integer, intent(out) :: iError

    character(len=*), parameter:: NameSub = 'IE_setnMlts'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_setnMlts
  !============================================================================
  subroutine IE_setnLats(iComponent, nLatsIn, iError)

    integer, intent(in)  :: iComponent, nLatsIn
    integer, intent(out) :: iError

    character(len=*), parameter:: NameSub = 'IE_setnLats'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_setnLats
  !============================================================================
  subroutine IE_put_from_im(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)
		use CON_router,   ONLY: IndexPtrType, WeightPtrType
    integer,intent(in)            :: nPoint, iPointStart, nVar
    real, intent(in)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    logical,intent(in)            :: DoAdd

    character(len=*), parameter:: NameSub = 'IE_put_from_im'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': IE_ERROR: empty version cannot be used!')

  end subroutine IE_put_from_im
  !============================================================================
  subroutine IE_put_from_im_complete

    !--------------------------------------------------------------------------
    write(*,*)"What is IE_put_from_im_complete really supposed to do?"

  end subroutine IE_put_from_im_complete
  !============================================================================

end module SIM_unused
!==============================================================================
