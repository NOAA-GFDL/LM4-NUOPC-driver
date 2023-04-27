module lnd_import_export

  use ESMF                    , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                    , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_END_ABORT, ESMF_Finalize
  use ESMF                    , only : ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF                    , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError, ESMF_FAILURE
  use ESMF                    , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                    , only : operator(/=), operator(==)
  use NUOPC                   , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model             , only : NUOPC_ModelGet
  use lm4_kind_mod            , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cs=>shr_kind_cs
  use nuopc_lm4_methods       , only : chkerr

  implicit none
  private ! except

  public  :: advertise_fields
  public  :: realize_fields
  ! public  :: import_fields
  ! public  :: export_fields

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_getfldptr
  ! private :: fldchk

  type fld_list_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToLnd_num = 0
  integer                :: fldsFrLnd_num = 0
  type (fld_list_type)   :: fldsToLnd(fldsMax)
  type (fld_list_type)   :: fldsFrLnd(fldsMax)

  ! from atm->lnd
  integer                :: ndep_nflds       ! number  of nitrogen deposition fields from atm->lnd/ocn

  ! from lnd->atm
  character(len=cx)      :: carma_fields     ! List of CARMA fields from lnd->atm
  integer                :: drydep_nflds     ! number of dry deposition velocity fields lnd-> atm
  integer                :: megan_nflds      ! number of MEGAN voc fields from lnd-> atm
  integer                :: emis_nflds       ! number of fire emission fields from lnd-> atm

  logical                :: flds_co2a        ! use case
  logical                :: flds_co2b        ! use case
  logical                :: flds_co2c        ! use case
  integer                :: glc_nec          ! number of glc elevation classes
  integer, parameter     :: debug = 0        ! internal debug level


  logical :: send_to_atm = .false.

  character(*),parameter :: F01 = "('(lnd_import_export) ',a,i5,2x,i5,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, flds_scalar_name,  rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: cvalue
    integer                :: n, num
    character(len=*), parameter :: subname='(lnd_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    send_to_atm = .true.
    
    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    ! export to atm
    if (send_to_atm) then
       ! TODO: actually set land frac for this field
       call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_lfrin'     )


       
    end if

    ! Now advertise above export fields
    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------


    ! from atm
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_z')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_tbot')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ta')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_tskn')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_pslv')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_prsl')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_pbot')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_shum')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_qa')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_u')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_v')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ua')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_va')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_exner')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ustar')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swdn')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_lwdn')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swnet')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainc')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainl')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rain')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snow')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowc')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowl')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'vfrac')
    call fldlist_add(fldsToLnd_num, fldsToLnd, 'zorl')
    ! needed?
    !call fldlist_add(fldsToLnd_num, fldsToLnd,'Faxa_garea')
    

    ! Now advertise import fields
    do n = 1,fldsToLnd_num
       call NUOPC_Advertise(importState, standardName=fldsToLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine advertise_fields

  !===============================================================================
  subroutine realize_fields(gcomp, mesh, grid, flds_scalar_name, flds_scalar_num, rc)

    use ESMF, only : ESMF_Mesh, ESMF_Grid
    
    ! input/output variables
    type(ESMF_GridComp) , intent(inout)          :: gcomp
    type(ESMF_Mesh)     , optional , intent(in)  :: mesh
    type(ESMF_Grid)     , optional , intent(in)  :: grid
    
    character(len=*)    , intent(in)             :: flds_scalar_name
    integer             , intent(in)             :: flds_scalar_num
    integer             , intent(out)            :: rc

    ! local variables
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    character(len=*), parameter :: subname='(lnd_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    if (present(mesh)) then
       call fldlist_realize( &
            state=ExportState, &
            fldList=fldsFrLnd, &
            numflds=fldsFrLnd_num, &
            flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, &
            tag=subname//':Land Export',&
            mesh=mesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call fldlist_realize( &
            state=importState, &
            fldList=fldsToLnd, &
            numflds=fldsToLnd_num, &
            flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, &
            tag=subname//':Land Import',&
            mesh=mesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
    else if (present(grid)) then
       call fldlist_realize( &
            state=ExportState, &
            fldList=fldsFrLnd, &
            numflds=fldsFrLnd_num, &
            flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, &
            tag=subname//':Land Export',&
            grid=grid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call fldlist_realize( &
            state=importState, &
            fldList=fldsToLnd, &
            numflds=fldsToLnd_num, &
            flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, &
            tag=subname//':Land Import',&
            grid=grid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
  end subroutine realize_fields


  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(lnd_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, grid, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_INDEX_DELOCAL, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_Grid, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh), optional , intent(in)    :: mesh
    type(ESMF_Grid), optional , intent(in)    :: grid
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(lnd_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             ! Create the field
             if (present(mesh)) then
                if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                   field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                        ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                        ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                        gridToFieldMap=(/2/), rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                else
                   field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
                call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                     ESMF_LOGMSG_INFO)
             else if (present(grid)) then
                ! Note no ungridded bounds. Hope this doesn't cause issues.
                   field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=stdname, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using grid", &
                     ESMF_LOGMSG_INFO)
             else
                call ESMF_LogWrite(subname // 'input must be grid or mesh', ESMF_LOGMSG_INFO)
                rc = ESMF_FAILURE
                return
             end if ! mesh or grid
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(lnd_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================
  subroutine state_getimport_1d(state, fldname, lm4data, rc)

    ! fill in lm4 import data for 1d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in)    :: state
    character(len=*) , intent(in)    :: fldname
    real(r8)         , intent(inout) :: lm4data(:)
    integer          , intent(out)   :: rc

    ! local variables
    real(r8), pointer :: fldPtr1d(:)
    integer           :: g
    character(len=*), parameter :: subname='(lnd_import_export:state_getimport_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call state_getfldptr(State, trim(fldname), fldptr1d=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do g = 1,size(lm4data)
       lm4data(g) = fldptr1d(g)
    end do
    !call check_for_nans(lm4data, trim(fldname), 1)

  end subroutine state_getimport_1d

  !===============================================================================
  ! subroutine state_getimport_2d(state, fldname, lm4data, rc)

  !   ! fill in lm4 import data for 2d field

  !   use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
  !   use ESMF, only : ESMF_Finalize

  !   ! input/output variabes
  !   type(ESMF_State) , intent(in)    :: state
  !   character(len=*) , intent(in)    :: fldname
  !   real(r8)         , intent(inout) :: lm4data(:,:)
  !   integer          , intent(out)   :: rc

  !   ! local variables
  !   real(r8), pointer :: fldPtr2d(:,:)
  !   integer           :: g,n
  !   character(len=CS) :: cnum
  !   character(len=*), parameter :: subname='(lnd_import_export:state_getimport_1d)'
  !   ! ----------------------------------------------

  !   rc = ESMF_SUCCESS

  !   call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
  !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !   do n = 1,size(lm4data, dim=2)
  !      write(cnum,'(i0)') n
  !      do g = 1,size(lm4data,dim=1)
  !         lm4data(g,n) = fldptr2d(n,g)
  !      end do
  !      call check_for_nans(lm4data(:,n), trim(fldname)//trim(cnum), 1)
  !   end do

  ! end subroutine state_getimport_2d

  ! !===============================================================================
  ! subroutine state_setexport_1d(state, fldname, lm4data, minus, rc)

  !   ! fill in lm4 export data for 1d field

  !   use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
  !   use ESMF, only : ESMF_Finalize

  !   ! input/output variabes
  !   type(ESMF_State) , intent(in) :: state
  !   character(len=*) , intent(in) :: fldname
  !   real(r8)         , intent(in) :: lm4data(:)
  !   logical, optional, intent(in) :: minus
  !   integer          , intent(out):: rc

  !   ! local variables
  !   real(r8), pointer :: fldPtr1d(:)
  !   integer           :: g
  !   character(len=*), parameter :: subname='(lnd_export_export:state_setexport_1d)'
  !   ! ----------------------------------------------

  !   call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
  !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !   fldptr1d(:) = 0._r8
  !   if (present(minus)) then
  !      do g = 1,size(lm4data)
  !         fldptr1d(g) = -lm4data(g)
  !      end do
  !   else
  !      do g = 1,size(lm4data)
  !         fldptr1d(g) = lm4data(g)
  !      end do
  !   end if
  !   call check_for_nans(lm4data, trim(fldname), 1)

  ! end subroutine state_setexport_1d

  ! !===============================================================================
  ! subroutine state_setexport_2d(state, fldname, lm4data, minus, rc)

  !   ! fill in lm4 export data for 2d field

  !   use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
  !   use ESMF, only : ESMF_Finalize

  !   ! input/output variabes
  !   type(ESMF_State) , intent(in) :: state
  !   character(len=*) , intent(in) :: fldname
  !   real(r8)         , intent(in) :: lm4data(:,:)
  !   logical, optional, intent(in) :: minus
  !   integer          , intent(out):: rc

  !   ! local variables
  !   real(r8), pointer :: fldPtr2d(:,:)
  !   integer           :: g, n
  !   character(len=CS) :: cnum
  !   character(len=*), parameter :: subname='(lnd_export_export:state_setexport_2d)'
  !   ! ----------------------------------------------

  !   rc = ESMF_SUCCESS

  !   call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
  !   if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !   fldptr2d(:,:) = 0._r8
  !   do n = 1,size(lm4data, dim=2)
  !      write(cnum,'(i0)') n
  !      if (present(minus)) then
  !         do g = 1,size(lm4data, dim=1)
  !            fldptr2d(n,g) = -lm4data(g,n)
  !         end do
  !      else
  !         do g = 1,size(lm4data, dim=1)
  !            fldptr2d(n,g) = lm4data(g,n)
  !         end do
  !      end if
  !      call check_for_nans(lm4data(:,n), trim(fldname)//trim(cnum), 1)
  !   end do

  ! end subroutine state_setexport_2d

  !===============================================================================
  subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fldname
    real(R8), pointer, optional , intent(out)   :: fldptr1d(:)
    real(R8), pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (present(fldptr1d)) then
       if (.not.associated(fldptr1d)) then
          write(*,*) 'fldptr1d not associated'
       endif
       call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (present(fldptr2d)) then
       if (.not.associated(fldptr2d)) then
          write(*,*) 'fldptr2d not associated'
       endif
       call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      call ESMF_LogWrite(trim(subname)//": either fldptr1d or fldptr2d must be an input argument", &
      ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)

      call ESMF_Finalize(endflag=ESMF_END_ABORT)

    end if

  end subroutine state_getfldptr

  !===============================================================================
  ! logical function fldchk(state, fldname)
  !   ! ----------------------------------------------
  !   ! Determine if field with fldname is in the input state
  !   ! ----------------------------------------------

  !   ! input/output variables
  !   type(ESMF_State), intent(in)  :: state
  !   character(len=*), intent(in)  :: fldname

  !   ! local variables
  !   type(ESMF_StateItem_Flag)   :: itemFlag
  !   ! ----------------------------------------------
  !   call ESMF_StateGet(state, trim(fldname), itemFlag)
  !   if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
  !      fldchk = .true.
  !   else
  !      fldchk = .false.
  !   endif
  ! end function fldchk

end module lnd_import_export
