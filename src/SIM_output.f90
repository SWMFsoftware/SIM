! Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with
! permission. For more information, see http://csem.engin.umich.edu/tools/swmf

module SIM_output

  use SIM_main

  implicit none

  character(len=100) :: NameOutFile = 'output.dat'

contains
  !============================================================================

  subroutine write_ie_output

    use ModIoUnit, ONLY: UnitTmp_

    integer :: i, j
    character(len=50) :: formatcode
    !--------------------------------------------------------------------------
    open(UnitTmp_, file=trim(NameOutputDir)//trim(NameOutFile), &
         action='write')

    write(UnitTmp_, *) 'i, j, Theta, Phi, Potential, Jr, Epsilon,'// &
                       'C, Sigma0, SigmaP, SigmaH'

    write(UnitTmp_, *) 'nTheta = ', nTheta
    write(UnitTmp_, *) 'nPhi   = ', nPhi

    formatcode = "(I, I, F7.2, F7.2, E10.4, E10.4, E10.4, E10.4, E10.4"// &
                 " ,E10.4, E10.4)"

    do j = 1,nPhi; do i = 1,nTheta
      write(UnitTmp_, *) i, j, Vars_VII(:, i, j)
    end do; end do

    close(UnitTmp_)

  end subroutine write_ie_output
  !============================================================================

end module SIM_output
!==============================================================================
