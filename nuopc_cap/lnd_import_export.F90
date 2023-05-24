module lnd_import_export
!-------------------------------------------------------------------------------
! This module contains the import and export routines for the land model, and
! related helper routines.
!-------------------------------------------------------------------------------

   use ESMF,                     only: ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
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
   !public  :: correct_import_fields

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

   logical :: send_to_atm = .true.

   character(*),parameter :: modName =  "(lnd_import_export)"
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
      ! call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainc')
      ! call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rainl')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_rain')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snow')
      ! call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowc')
      ! call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_snowl')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'vfrac')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'zorl')
      ! needed?
      !call fldlist_add(fldsToLnd_num, fldsToLnd,'Faxa_garea')

      ! additional provided by CDEPS
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdf')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndf')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swvdr')
      call fldlist_add(fldsToLnd_num, fldsToLnd, 'Faxa_swndr')



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

            !! Set flag for connected fields
            !fldList(n)%connected = .true.

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
   subroutine import_fields(gcomp, cplr2land, lm4_model, rc)

      ! input/output variables
      type(ESMF_GridComp),              intent(in)    :: gcomp
      type(atmos_land_boundary_type),   intent(inout) :: cplr2land
      type(lm4_type),                   intent(inout) :: lm4_model
      integer,                          intent(out)   :: rc

      ! local variables
      real(r8), dimension(:,:), pointer  :: datar8

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


      allocate(datar8(lnd%is:lnd%ie,lnd%js:lnd%je))

      if (ie_debug > 0) then ! want both U.G. and S.G. data
         call state_getimport_2d(importState, 'Sa_z',       lm4data_1d=lm4_model%atm_forc%z_bot,   lm4data_2d=lm4_model%atm_forc2d%z_bot,   rc=rc)
         !call state_getimport_2d(importState, 'Sa_ta',      lm4data_1d=lm4_model%atm_forc%t_bot,   lm4data_2d=lm4_model%atm_forc2d%t_bot, rc=rc)
         call state_getimport_2d(importState, 'Sa_tbot',    lm4data_1d=lm4_model%atm_forc%t_bot,   lm4data_2d=lm4_model%atm_forc2d%t_bot,   rc=rc)
         ! call state_getimport_2d(importState, 'Sa_tskn'
         ! call state_getimport_2d(importState, 'Sa_prsl',   lm4data_1d=lm4_model%atm_forc%prsl,    lm4data_2d=lm4_model%atm_forc2d%prsl ,rc=rc)
         call state_getimport_2d(importState, 'Sa_pbot',    lm4data_1d=lm4_model%atm_forc%p_bot,   lm4data_2d=lm4_model%atm_forc2d%p_bot,   rc=rc)
         ! call state_getimport_2d(importState, 'Sa_pslv'
         call state_getimport_2d(importState, 'Sa_u',       lm4data_1d=lm4_model%atm_forc%u_bot,   lm4data_2d=lm4_model%atm_forc2d%u_bot,   rc=rc)
         call state_getimport_2d(importState, 'Sa_v',       lm4data_1d=lm4_model%atm_forc%v_bot,   lm4data_2d=lm4_model%atm_forc2d%v_bot,   rc=rc)
         call state_getimport_2d(importState, 'Faxa_rain',  lm4data_1d=lm4_model%atm_forc%tprec,   lm4data_2d=lm4_model%atm_forc2d%tprec,   rc=rc)
         call state_getimport_2d(importState, 'Faxa_snow',  lm4data_1d=lm4_model%atm_forc%fprec,   lm4data_2d=lm4_model%atm_forc2d%fprec,   rc=rc)
         call state_getimport_2d(importState, 'Faxa_lwdn',  lm4data_1d=lm4_model%atm_forc%flux_lw, lm4data_2d=lm4_model%atm_forc2d%flux_lw, rc=rc)
         call state_getimport_2d(importState, 'Faxa_swvdf', lm4data_1d=lm4_model%atm_forc%flux_sw_down_vis_dif, &
            lm4data_2d=lm4_model%atm_forc2d%flux_sw_down_vis_dif, rc=rc)
         call state_getimport_2d(importState, 'Faxa_swvdr', lm4data_1d=lm4_model%atm_forc%flux_sw_down_vis_dir, &
            lm4data_2d=lm4_model%atm_forc2d%flux_sw_down_vis_dir, rc=rc)
      else ! want only Unstructured Grid data
         call state_getimport_2d(importState, 'Sa_z',       lm4data_1d=lm4_model%atm_forc%z_bot, rc=rc)
         !call state_getimport_2d(importState, 'Sa_ta', lm4data_1d=lm4_model%atm_forc%t_bot, rc=rc)
         call state_getimport_2d(importState, 'Sa_tbot',    lm4data_1d=lm4_model%atm_forc%t_bot, rc=rc)
         ! call state_getimport_2d(importState, 'Sa_tskn'
         ! call state_getimport_2d(importState, 'Sa_prsl', lm4data_1d=lm4_model%atm_forc%prsl, rc=rc)
         call state_getimport_2d(importState, 'Sa_pbot',    lm4data_1d=lm4_model%atm_forc%p_bot, rc=rc)
         ! call state_getimport_2d(importState, 'Sa_pslv'
         call state_getimport_2d(importState, 'Sa_u',       lm4data_1d=lm4_model%atm_forc%u_bot, rc=rc)
         call state_getimport_2d(importState, 'Sa_v',       lm4data_1d=lm4_model%atm_forc%v_bot, rc=rc)
         call state_getimport_2d(importState, 'Faxa_rain',  lm4data_1d=lm4_model%atm_forc%tprec, rc=rc)
         call state_getimport_2d(importState, 'Faxa_snow',  lm4data_1d=lm4_model%atm_forc%fprec, rc=rc)
         call state_getimport_2d(importState, 'Faxa_lwdn',  lm4data_1d=lm4_model%atm_forc%flux_lw, rc=rc)
         call state_getimport_2d(importState, 'Faxa_swvdf', lm4data_1d=lm4_model%atm_forc%flux_sw_down_vis_dif, rc=rc)
         call state_getimport_2d(importState, 'Faxa_swvdr', lm4data_1d=lm4_model%atm_forc%flux_sw_down_vis_dir, rc=rc)
      end if

      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! call state_getimport_2d(importState, 'Sa_shum'   , cplr2land%, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! call state_getimport_2d(importState, 'Sa_qa'     , cplr2land%, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
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

      ! ! call state_getimport_2d(importState, 'Faxa_rainc', cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! ! call state_getimport_2d(importState, 'Faxa_rainl', cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! ! call state_getimport_2d(importState, 'Faxa_snowc', cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! ! call state_getimport_2d(importState, 'Faxa_snowl', cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! ! call state_getimport_2d(importState, 'vfrac'     , cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! ! call state_getimport_2d(importState, 'zorl'      , cplr2land%, rc=rc)
      ! ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

      deallocate(datar8)

      call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

   end subroutine import_fields

   !=============================================================================

   ! subroutine correct_import_fields(gcomp, cplr2land, lm4_model, rc)
   !    !! This routine "corrects" the imported fields from the atmosphere into what the land model expects

   !    ! input/output variables
   !    type(ESMF_GridComp),              intent(in)    :: gcomp
   !    type(atmos_land_boundary_type),   intent(inout) :: cplr2land
   !    type(lm4_type),                   intent(inout) :: lm4_model
   !    integer,                          intent(out)   :: rc

   ! end subroutine correct_import_fields

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

      use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
      use ESMF, only : ESMF_Finalize

      ! input/output variabes
      type(ESMF_State) , intent(in)    :: state
      character(len=*) , intent(in)    :: fldname
      real(r8)         , intent(inout) :: lm4data_1d(:)
      integer          , intent(out)   :: rc

      ! local variables
      real(r8), pointer :: fldPtr1d(:)
      type(ESMF_StateItem_Flag)   :: itemType
      character(len=*), parameter :: subname='(lnd_import_export:state_getimport_1d)'
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
         call ESMF_LogWrite(subname//' '//trim(fldname)//' is not in the state!', ESMF_LOGMSG_INFO)
      end if


   end subroutine state_getimport_1d

   !===============================================================================
   subroutine state_getimport_2d(state, fldname, lm4data_1d, lm4data_2d, rc)

      ! fill in lm4 import data for 2d field
      ! In this context, 2d is expected to be structured grid data

      use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
      use ESMF, only : ESMF_Finalize

      ! input/output variabes
      type(ESMF_State),  intent(in)    :: state
      character(len=*),  intent(in)    :: fldname
      real(r8),optional, intent(inout) :: lm4data_1d(:)      ! 1d, Unstructured Grid output
      real(r8),optional, intent(inout) :: lm4data_2d(:,:)    ! 2d, Structured Grid output
      integer,           intent(out)   :: rc

      ! local variables
      real(r8), pointer :: fldPtr2d(:,:)
      type(ESMF_StateItem_Flag)   :: itemType
      character(len=*), parameter :: subname='(lnd_import_export:state_getimport_2d)'
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
      use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
      use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE
      ! TMP DEBUG
      use ESMF , only : ESMF_ArraySpec
      ! END TMP DEBUG

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

      ! TMP DEBUG VARS
      type(ESMF_ArraySpec) :: arrayspec
      ! END TMP DEBUG VARS

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



end module lnd_import_export
