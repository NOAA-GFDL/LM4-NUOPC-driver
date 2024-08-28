!!===============================================================================
!! This module contains the import and export routines for the land model, and
!! related helper routines.
!!===============================================================================
module lm4_import_export

   use ESMF,                     only: ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
   use ESMF,                     only: ESMF_Field, ESMF_FieldGet, ESMF_LOGERR_PASSTHRU
   use ESMF,                     only: ESMF_KIND_R8, ESMF_SUCCESS, ESMF_END_ABORT, ESMF_Finalize
   use ESMF,                     only: ESMF_MAXSTR, ESMF_LOGMSG_INFO
   use ESMF,                     only: ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError, ESMF_FAILURE
   use ESMF,                     only: ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
   use ESMF,                     only: operator(/=), operator(==)
   use ESMF,                     only : ESMF_StateItem_Flag, ESMF_STATEITEM_FIELD
   use NUOPC,                    only: NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
   use NUOPC_Model,              only: NUOPC_ModelGet
   use lm4_kind_mod,             only: r8 => shr_kind_r8, cx=>shr_kind_cx, cs=>shr_kind_cs
   use lm4_type_mod,             only: lm4_type
   use land_data_mod,            only: lnd ! global data
   use land_data_mod,            only: land_data_type, atmos_land_boundary_type
   use nuopc_lm4_methods,        only: chkerr
   use mpp_domains_mod,          only : mpp_pass_sg_to_ug

   implicit none
   private ! except

   public  :: advertise_fields
   public  :: realize_fields
   public  :: import_fields
   public  :: export_fields
   public  :: correct_import_fields

   private :: fldlist_add
   private :: fldlist_realize
   private :: state_getfldptr
   ! private :: fldchk

   type fld_list_type
      character(len=128) :: stdname
      integer :: ungridded_lbound = 0
      integer :: ungridded_ubound = 0
      logical :: connected = .false.
   end type fld_list_type

   integer, parameter     :: fldsMax = 100
   integer                :: fldsToLnd_num = 0
   integer                :: fldsFrLnd_num = 0
   type (fld_list_type)   :: fldsToLnd(fldsMax)
   type (fld_list_type)   :: fldsFrLnd(fldsMax)

   integer                :: ie_debug         ! internal debug level

   logical                :: send_to_atm = .true.

   character(*),parameter :: modName =  "(lm4_import_export)"
   character(*),parameter :: F01 = "('(lm4_import_export) ',a,i5,2x,i5,2x,d21.14)"
   character(*),parameter :: u_FILE_u = &
      __FILE__

!===============================================================================
contains

   !===============================================================================
   subroutine advertise_fields(gcomp, rc)

      ! input/output variables
      type(ESMF_GridComp)            :: gcomp
      integer          , intent(out) :: rc

      ! local variables
      type(ESMF_State)       :: importState
      type(ESMF_State)       :: exportState
      character(ESMF_MAXSTR) :: cvalue
      integer                :: n, num
      character(len=*), parameter :: subname='(lm4_import_export:advertise_fields)'
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return


      !--------------------------------
      ! Advertise export fields
      !--------------------------------

      ! export to atm
      if (send_to_atm) then
         ! TODO: actually set land frac for this field
         call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'Sl_lfrin')

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
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_z')       ! atmosphere export - bottom layer height
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_tbot')    ! atmosphere export - bottom layer temperature
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ta')      ! atmosphere export - bottom layer temperature
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_tskn')    ! atmosphere export - sea surface skin temperature
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_pslv')    ! atmosphere export - instantaneous pressure land and sea surface
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_prsl')    ! atmosphere export - pressure at lowest model layer
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_pbot')    ! atmosphere export - pressure at lowest model layer
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_shum')    ! atmosphere export - bottom layer specific humidity
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_qa')      ! atmosphere export - bottom layer specific humidity
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_u')       ! atmosphere export - bottom layer zonal wind
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_v')       ! atmosphere export - bottom layer meridional wind
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_exner')   ! dimensionless exner function at surface adjacent layer
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Sa_ustar')   ! surface friction velocity
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swdn')  ! atmosphere export -  mean downward SW heat flux
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_lwdn')  ! atmosphere export - mean downward LW heat flux
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swnet') ! mean_net_sw_flx
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainc')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainl')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rain')  ! mean_prec_rate
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snow')  ! mean_fprec_rate
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowc')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowl')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'vfrac')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'zorl')
      ! needed?
      !call fldlist_add(fldsToLnd_num, fldsToLnd,'Faxa_garea')

      ! additional provided by CDEPS
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdf') ! atmosphere export - mean surface downward uv+vis diffuse flux
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndf') ! atmosphere export - mean surface downward nir diffuse flux
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdr') ! atmosphere export - mean surface downward uv+visvdirect flux
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndr') ! atmosphere export - mean surface downward nir direct flux

      ! Needed by CMEPS
      call fldlist_add(fldsFrLnd_num, fldsFrlnd, 'cpl_scalars')


      ! Now advertise import fields
      do n = 1,fldsToLnd_num
         call NUOPC_Advertise(importState, standardName=fldsToLnd(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      enddo

   end subroutine advertise_fields

   !===============================================================================
   subroutine realize_fields(gcomp, lm4_model, mesh, grid, rc)

      use ESMF, only : ESMF_Mesh, ESMF_Grid

      ! input/output variables
      type(ESMF_GridComp) , intent(inout)          :: gcomp
      type(lm4_type)      , intent(inout)          :: lm4_model
      type(ESMF_Mesh)     , optional , intent(in)  :: mesh
      type(ESMF_Grid)     , optional , intent(in)  :: grid
      integer             , intent(out)            :: rc

      ! local variables
      type(ESMF_State)     :: importState
      type(ESMF_State)     :: exportState
      character(len=*), parameter :: subname='(lm4_import_export:realize_fields)'
      real(R8)          :: scalardim(3)
      !---------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return


      if (present(mesh)) then
         call fldlist_realize( &
            state=ExportState, &
            fldList=fldsFrLnd, &
            numflds=fldsFrLnd_num, &
            flds_scalar_name=lm4_model%cpl_scalar%flds_scalar_name, &
            flds_scalar_num=lm4_model%cpl_scalar%flds_scalar_num, &
            tag=subname//':Land Export',&
            mesh=mesh, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         call fldlist_realize( &
            state=importState, &
            fldList=fldsToLnd, &
            numflds=fldsToLnd_num, &
            flds_scalar_name=lm4_model%cpl_scalar%flds_scalar_name, &
            flds_scalar_num=lm4_model%cpl_scalar%flds_scalar_num, &
            tag=subname//':Land Import',&
            mesh=mesh, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

      else if (present(grid)) then
         call fldlist_realize( &
            state=ExportState, &
            fldList=fldsFrLnd, &
            numflds=fldsFrLnd_num, &
            flds_scalar_name=lm4_model%cpl_scalar%flds_scalar_name, &
            flds_scalar_num=lm4_model%cpl_scalar%flds_scalar_num, &
            tag=subname//':Land Export',&
            grid=grid, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         call fldlist_realize( &
            state=importState, &
            fldList=fldsToLnd, &
            numflds=fldsToLnd_num, &
            flds_scalar_name=lm4_model%cpl_scalar%flds_scalar_name, &
            flds_scalar_num=lm4_model%cpl_scalar%flds_scalar_num, &
            tag=subname//':Land Import',&
            grid=grid, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      end if

      ! cpl_scalars for export state
      scalardim = 0.0
      scalardim(1) = real(lm4_model%nml%npx,8)
      scalardim(2) = real(lm4_model%nml%npy,8)
      scalardim(3) = real(lm4_model%nml%ntiles,8)

      if (lm4_model%cpl_scalar%flds_scalar_num > 0) then
         ! Set the scalar data into the exportstate
         call State_SetScalar(scalardim(1), lm4_model%cpl_scalar%flds_scalar_index_nx, exportState, &
            lm4_model%cpl_scalar%flds_scalar_name, lm4_model%cpl_scalar%flds_scalar_num, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         call State_SetScalar(scalardim(2), lm4_model%cpl_scalar%flds_scalar_index_ny, exportState, &
            lm4_model%cpl_scalar%flds_scalar_name, lm4_model%cpl_scalar%flds_scalar_num, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         call State_SetScalar(scalardim(3), lm4_model%cpl_scalar%flds_scalar_index_ntile, exportState, &
            lm4_model%cpl_scalar%flds_scalar_name, lm4_model%cpl_scalar%flds_scalar_num, rc)
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
      character(len=*), parameter :: subname='(lm4_import_export:fldlist_add)'
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
      use ESMF  , only : ESMF_MAXSTR, ESMF_State, ESMF_Mesh, ESMF_Grid, ESMF_StateRemove
      use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
      use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR

      ! input/output variables
      type(ESMF_State)    , intent(inout) :: state
      type(fld_list_type) , intent(inout) :: fldList(:)
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
      character(len=*),parameter  :: subname='(lm4_import_export:fldlist_realize)'
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

            !! Set flag for connected fields
            fldList(n)%connected = .true.

         else
            if (stdname /= trim(flds_scalar_name)) then
               call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
               call ESMF_StateRemove(state, (/stdname/), rc=rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            end if
         end if
      end do

   end subroutine fldlist_realize


   subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(lm4_import_export:SetScalarField)'
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


   subroutine State_SetScalar(scalar_value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)

      use ESMF, only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet

      ! input/output arguments
      real(ESMF_KIND_R8), intent(in)   :: scalar_value
      integer,          intent(in)     :: scalar_id
      type(ESMF_State), intent(inout)  :: State
      character(len=*), intent(in)     :: flds_scalar_name
      integer,          intent(in)     :: flds_scalar_num
      integer,          intent(inout)  :: rc

      ! local variables
      integer           :: mytask
      type(ESMF_Field)  :: lfield
      type(ESMF_VM)     :: vm
      real(ESMF_KIND_R8), pointer :: farrayptr(:,:)

      character(len=*), parameter :: subname = ' (lm4_import_export:state_setscalar) '
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      call ESMF_VMGetCurrent(vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_VMGet(vm, localPet=mytask, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=lfield, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if (mytask == 0) then
         call ESMF_FieldGet(lfield, farrayPtr = farrayptr, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
         if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
            call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO)
            rc = ESMF_FAILURE
            return
         endif
         farrayptr(scalar_id,1) = scalar_value
      endif

   end subroutine State_SetScalar


   ! Import fields that do not need to be altered for the land model
   !!===============================================================================
   subroutine import_fields(gcomp, lm4_model, rc)

      ! input/output variables
      type(ESMF_GridComp),              intent(in)    :: gcomp
      type(lm4_type),                   intent(inout) :: lm4_model
      integer,                          intent(out)   :: rc

      ! local variables
      type(ESMF_State)            :: importState
      character(len=*), parameter :: subname=trim(modName)//':(import_fields)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      ie_debug = lm4_model%nml%lm4_debug

      ! Get import state
      call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      !! Atm input fields
      !  At least for data atm, atm imports are 2d (ie, LM4's Structured Grid)
      !  need to convert to LM4's 1d Unstructured Grid
      ! -----------------------

      ! Get Unstructured Grid data
      call state_getimport_2d(importState, 'Sa_z',       lm4data_1d=lm4_model%atm_forc%z_bot, rc=rc)  ! bottom layer height
      call state_getimport_2d(importState, 'Sa_tbot',    lm4data_1d=lm4_model%atm_forc%t_bot, rc=rc)  ! bottom layer temperature
      !call state_getimport_2d(importState, 'Sa_ta',      lm4data_1d=lm4_model%atm_forc%t_bot, rc=rc)  ! bottom layer temperature (active UFS atm)
      ! call state_getimport_2d(importState, 'Sa_tskn' ...                                            ! surface skin temperature
      call state_getimport_2d(importState, 'Sa_pbot',    lm4data_1d=lm4_model%atm_forc%p_bot, rc=rc)  ! bottom layer pressure
      !call state_getimport_2d(importState, 'Sa_prsl',    lm4data_1d=lm4_model%atm_forc%p_bot, rc=rc)  ! bottom layer pressure (active UFS atm)
      call state_getimport_2d(importState, 'Sa_u',       lm4data_1d=lm4_model%atm_forc%u_bot, rc=rc)  ! bottom layer zonal wind
      call state_getimport_2d(importState, 'Sa_v',       lm4data_1d=lm4_model%atm_forc%v_bot, rc=rc)  ! bottom layer meridional wind
      call state_getimport_2d(importState, 'Sa_shum',    lm4data_1d=lm4_model%atm_forc%q_bot, rc=rc)  ! bottom layer specific humidity
      !call state_getimport_2d(importState, 'Sa_qa',      lm4data_1d=lm4_model%atm_forc%q_bot, rc=rc)  ! bottom layer specific humidity (active UFS atm)
      call state_getimport_2d(importState, 'Sa_pslv',    lm4data_1d=lm4_model%atm_forc%p_surf, rc=rc) ! surface pressure
      call state_getimport_2d(importState, 'Faxa_lwdn',  lm4data_1d=lm4_model%atm_forc%flux_lw, rc=rc)
      call state_getimport_2d(importState, 'Faxa_swvdf', lm4data_1d=lm4_model%atm_forc%flux_sw_down_vis_dif, rc=rc) ! mean surface downward uv+vis diffuse flux
      call state_getimport_2d(importState, 'Faxa_swvdr', lm4data_1d=lm4_model%atm_forc%flux_sw_down_vis_dir, rc=rc) ! mean surface downward uv+vis direct flux
      call state_getimport_2d(importState, 'Faxa_swndf', lm4data_1d=lm4_model%atm_forc%flux_sw_down_nir_dif, rc=rc) ! mean surface downward nir diffuse flux
      call state_getimport_2d(importState, 'Faxa_swndr', lm4data_1d=lm4_model%atm_forc%flux_sw_down_nir_dir, rc=rc) ! mean surface downward nir direct flux

      if (ie_debug > 0) then ! Also want Structured Grid data
         call state_getimport_2d(importState, 'Sa_z',       lm4data_2d=lm4_model%atm_forc2d%z_bot,   rc=rc)
         call state_getimport_2d(importState, 'Sa_tbot',    lm4data_2d=lm4_model%atm_forc2d%t_bot,   rc=rc)
         !call state_getimport_2d(importState, 'Sa_ta',      lm4data_2d=lm4_model%atm_forc2d%t_bot,    rc=rc)
         ! call state_getimport_2d(importState, 'Sa_tskn'
         call state_getimport_2d(importState, 'Sa_pbot',    lm4data_2d=lm4_model%atm_forc2d%p_bot,   rc=rc)
         !call state_getimport_2d(importState, 'Sa_prsl',    lm4data_2d=lm4_model%atm_forc2d%p_bot, rc=rc)
         call state_getimport_2d(importState, 'Sa_u',       lm4data_2d=lm4_model%atm_forc2d%u_bot,   rc=rc)
         call state_getimport_2d(importState, 'Sa_v',       lm4data_2d=lm4_model%atm_forc2d%v_bot,   rc=rc)
         call state_getimport_2d(importState, 'Sa_shum',    lm4data_2d=lm4_model%atm_forc2d%q_bot,   rc=rc)
         !call state_getimport_2d(importState, 'Sa_qa',      lm4data_2d=lm4_model%atm_forc2d%q_bot,   rc=rc)
         call state_getimport_2d(importState, 'Sa_pslv',    lm4data_2d=lm4_model%atm_forc2d%p_surf,  rc=rc)
         call state_getimport_2d(importState, 'Faxa_lwdn',  lm4data_2d=lm4_model%atm_forc2d%flux_lw, rc=rc)
         call state_getimport_2d(importState, 'Faxa_swvdf', lm4data_2d=lm4_model%atm_forc2d%flux_sw_down_vis_dif, rc=rc)
         call state_getimport_2d(importState, 'Faxa_swvdr', lm4data_2d=lm4_model%atm_forc2d%flux_sw_down_vis_dir, rc=rc)
         call state_getimport_2d(importState, 'Faxa_swndf', lm4data_2d=lm4_model%atm_forc2d%flux_sw_down_nir_dif, rc=rc)
         call state_getimport_2d(importState, 'Faxa_swndr', lm4data_2d=lm4_model%atm_forc2d%flux_sw_down_nir_dir, rc=rc)
      end if


      ! call state_getimport_2d(importState, 'Faxa_swdn' , cplr2land%swdn_flux, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! call state_getimport_2d(importState, 'Faxa_swnet', forc%flux_sw, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! -----------------------

      ! call state_getimport_2d(importState, 'Faxa_swndf',
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! call state_getimport_2d(importState, 'Faxa_swndr',
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! -----------------------

      ! ! call state_getimport_2d(importState, 'Sa_exner'  , cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! ! call state_getimport_2d(importState, 'Sa_ustar'  , cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! ! call state_getimport_2d(importState, 'vfrac'     , cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! ! call state_getimport_2d(importState, 'zorl'      , cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

   end subroutine import_fields

   !! Imports and "corrects" the fields from the atmosphere that need to be
   !! altered somehow into what the land model expects
   !!=============================================================================
   subroutine correct_import_fields(gcomp, lm4_model, rc)


      ! input/output variables
      type(ESMF_GridComp), intent(in)    :: gcomp
      type(lm4_type),      intent(inout) :: lm4_model
      integer,             intent(out)   :: rc

      ! local variables
      type(ESMF_State)                  :: importState
      character(len=*),       parameter :: subname=trim(modName)//':(correct_import_fields)'
      real(r8), dimension(:), pointer   :: tmp_ug_data ! temp unstructured grid data
      real(r8), dimension(:,:), pointer   :: tmp_sg_data ! TMP DEBUG

      integer :: i ! loop counter TMP DEBUG
      ! ----------------------------------------------

      rc = ESMF_SUCCESS
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      ie_debug = lm4_model%nml%lm4_debug

      ! Get import state
      call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      allocate(tmp_ug_data(lnd%ls:lnd%le))
      if (ie_debug > 0) then
         allocate(tmp_sg_data(lnd%is:lnd%ie,lnd%js:lnd%je))
      endif


      !! precip
      !! ---------------------------------------------------------------------
      ! rain
      if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_rain')) then
         ! have total liquid precip
         call state_getimport_2d(importState, 'Faxa_rain',  lm4data_1d=lm4_model%atm_forc%lprec, rc=rc)
         if (ie_debug > 0) then
            call state_getimport_2d(importState, 'Faxa_rain', lm4data_2d=lm4_model%atm_forc2d%lprec, rc=rc)
         endif
      elseif ( ( check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_rainc') ) .and. &
         ( check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_rainl') ) ) then
         ! have convective and large-scale total precip
         call state_getimport_2d(importState, 'Faxa_rainc',  lm4data_1d=lm4_model%atm_forc%lprec, rc=rc)
         call state_getimport_2d(importState, 'Faxa_rainl',  lm4data_1d=tmp_ug_data, rc=rc)
         lm4_model%atm_forc%lprec = lm4_model%atm_forc%lprec + tmp_ug_data
         if (ie_debug > 0) then
            call state_getimport_2d(importState, 'Faxa_rainc', lm4data_2d=lm4_model%atm_forc2d%lprec, rc=rc)
            call state_getimport_2d(importState, 'Faxa_rainl', lm4data_2d=tmp_sg_data, rc=rc)
            lm4_model%atm_forc2d%lprec = lm4_model%atm_forc2d%lprec + tmp_sg_data
         endif
      else
         call ESMF_LogWrite(trim(subname)//": Don't have any liquid precip fields", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif

      ! snow
      if (check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_snow')) then
         ! have snow precip
         call state_getimport_2d(importState, 'Faxa_snow',  lm4data_1d=lm4_model%atm_forc%fprec, rc=rc)
         if (ie_debug > 0) then
            call state_getimport_2d(importState, 'Faxa_snow', lm4data_2d=lm4_model%atm_forc2d%fprec, rc=rc)
         endif
      elseif ( ( check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_snowc') ) .and. &
         ( check_for_connected(fldsToLnd, fldsToLnd_num, 'Faxa_snowl') ) ) then
         ! have convective and large-scale snow precip
         call state_getimport_2d(importState, 'Faxa_snowc',  lm4data_1d=lm4_model%atm_forc%fprec, rc=rc)
         call state_getimport_2d(importState, 'Faxa_snowl',  lm4data_1d=tmp_ug_data, rc=rc)
         lm4_model%atm_forc%fprec = lm4_model%atm_forc%fprec + tmp_ug_data
         if (ie_debug > 0) then
            call state_getimport_2d(importState, 'Faxa_snowc', lm4data_2d=lm4_model%atm_forc2d%fprec, rc=rc)
            call state_getimport_2d(importState, 'Faxa_snowl', lm4data_2d=tmp_sg_data, rc=rc)
            lm4_model%atm_forc2d%fprec = lm4_model%atm_forc2d%fprec + tmp_sg_data
         endif
      else
         ! TODO: want to have FMS Coupler's precip scaling functionality here?
         call ESMF_LogWrite(trim(subname)//": Don't have any snow precip fields", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif


      deallocate(tmp_ug_data)
      if (ie_debug > 0) then
         deallocate(tmp_sg_data)
      endif

   end subroutine correct_import_fields

   !===============================================================================
   subroutine export_fields(gcomp, rc)

      ! input/output variables
      type(ESMF_GridComp),              intent(in)    :: gcomp
      !type(lm4_type),                   intent(in)    :: lm4_model
      integer,                          intent(out)   :: rc

      ! local variables
      type(ESMF_State)            :: exportState
      character(len=*), parameter :: subname=trim(modName)//':(export_fields)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      ! Get export state
      call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! export to mediator
      call state_setexport_2d(exportState, 'Sl_lfrin', lnd%sg_landfrac,rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return


   end subroutine export_fields

   !===============================================================================
   subroutine state_getimport_1d(state, fldname, lm4data_1d, rc)

      ! fill in lm4 import data for 1d field

      use ESMF, only : ESMF_END_ABORT, ESMF_LogFoundError
      use ESMF, only : ESMF_Finalize

      ! input/output variabes
      type(ESMF_State) , intent(in)    :: state
      character(len=*) , intent(in)    :: fldname
      real(r8)         , intent(inout) :: lm4data_1d(:)
      integer          , intent(out)   :: rc

      ! local variables
      real(r8), pointer :: fldPtr1d(:)
      type(ESMF_StateItem_Flag)   :: itemType
      character(len=*), parameter :: subname='(lm4_import_export:state_getimport_1d)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      call ESMF_StateGet(state, itemName=trim(fldname), itemType=itemType, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (itemType == ESMF_STATEITEM_FIELD) then

         call state_getfldptr(State, trim(fldname), fldptr1d=fldptr1d, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         lm4data_1d(:) = fldptr1d(:)
         !call check_for_nans(lm4data_1d, trim(fldname), 1)
      else
         call ESMF_LogWrite(subname//' '//trim(fldname)//' is not in the state!', ESMF_LOGMSG_ERROR)
         rc = ESMF_FAILURE
         return
      end if


   end subroutine state_getimport_1d

   !===============================================================================
   subroutine state_getimport_2d(state, fldname, lm4data_1d, lm4data_2d, rc)

      ! fill in lm4 import data for 2d field
      ! In this context, 2d is expected to be structured grid data

      use ESMF, only : ESMF_END_ABORT, ESMF_LogFoundError
      use ESMF, only : ESMF_Finalize

      ! input/output variabes
      type(ESMF_State),  intent(in)    :: state
      character(len=*),  intent(in)    :: fldname
      real(r8), optional, intent(inout) :: lm4data_1d(:)      ! 1d, Unstructured Grid output
      real(r8), optional, intent(inout) :: lm4data_2d(:,:)    ! 2d, Structured Grid output
      ! logical,  optional, intent(in)  :: debug_print
      integer,           intent(out)   :: rc

      ! local variables
      real(r8), pointer :: fldPtr2d(:,:)
      type(ESMF_StateItem_Flag)   :: itemType
      character(len=*), parameter :: subname='(lm4_import_export:state_getimport_2d)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      call ESMF_StateGet(state, itemName=trim(fldname), itemType=itemType, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (itemType == ESMF_STATEITEM_FIELD) then

         call state_getfldptr(State, trim(fldname), fldptr2d=fldptr2d, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         ! pass structured grid data to structured grid
         if (present(lm4data_2d)) then
            lm4data_2d(:,:) = fldptr2d(:,:)
         end if

         ! pass structured grid data to unstructured grid
         if (present(lm4data_1d)) then
            call mpp_pass_sg_to_ug(lnd%ug_domain, fldptr2d, lm4data_1d)
         end if

         ! if (present(debug_print)) then
         !    if (debug_print) then
         !       call ESMF_LogWrite(subname//' '//trim(fldname)//' min/max: '// &
         !                          trim(str(fldptr2d%min))//' '//trim(str(fldptr2d%max)), ESMF_LOGMSG_INFO)

         !    end if
         ! end if


      else
         call ESMF_LogWrite(subname//' '//trim(fldname)//' is not in the state!', ESMF_LOGMSG_INFO)
      end if


   end subroutine state_getimport_2d

   !===============================================================================
   subroutine state_setexport_2d(state, fldname, lm4data, minus, rc)
      ! fill in lm4 export data for 2d field (ie, structured grid)
      use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
      use ESMF, only : ESMF_Finalize

      ! input/output variabes
      type(ESMF_State) , intent(in)  :: state
      character(len=*) , intent(in)  :: fldname
      real(r8)         , intent(in)  :: lm4data(:,:)
      logical, optional, intent(in)  :: minus
      integer          , intent(out) :: rc

      ! local variables
      real(r8), pointer :: fldPtr2d(:,:)
      character(len=*), parameter :: subname='(lnd_export_export:state_setexport_1d)'
      ! ----------------------------------------------


      call state_getfldptr(state, trim(fldname), fldptr2d=fldptr2d, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (present(minus)) then
         fldptr2d(:,:) = -lm4data(:,:)
      else
         fldptr2d(:,:) = lm4data(:,:)
      end if

   end subroutine state_setexport_2d


   !===============================================================================
   subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

      ! ----------------------------------------------
      ! Get pointer to a state field
      ! ----------------------------------------------

      use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
      use ESMF , only : ESMF_StateGet, ESMF_MeshGet
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
      character(len=*), parameter :: subname='(lm4_import_export:state_getfldptr)'
      ! ----------------------------------------------


      rc = ESMF_SUCCESS

      call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (present(fldptr1d)) then
         call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      else if (present(fldptr2d)) then
         call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      else
         call ESMF_LogWrite(trim(subname)//": either fldptr1d or fldptr2d must be an input argument", &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)

      end if

   end subroutine state_getfldptr



   !=============================================================================
   logical function check_for_connected(fldList, numflds, fname)

      ! input/output variables
      type(fld_list_type) , intent(inout) :: fldList(:)
      integer             , intent(in)    :: numflds
      character(len=*)    , intent(in)    :: fname

      ! local variables
      integer :: n
      character(len=*), parameter :: subname='(check_for_connected)'
      ! ----------------------------------------------

      check_for_connected = .false.
      do n = 1, numflds
         if (trim(fname) == trim(fldList(n)%stdname)) then
            check_for_connected = fldList(n)%connected
            exit
         end if
      end do

   end function check_for_connected


end module lm4_import_export
