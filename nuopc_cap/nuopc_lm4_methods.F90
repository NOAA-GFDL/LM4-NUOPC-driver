!! Based of of CMEPS nuopc_shr_methods

module nuopc_lm4_methods

  use ESMF         , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError

  implicit none
  private
  public  :: chkerr

!===============================================================================
contains
!===============================================================================

  logical function chkerr(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

end module nuopc_lm4_methods
