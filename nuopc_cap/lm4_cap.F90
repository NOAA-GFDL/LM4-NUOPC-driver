! model test after CTSM

module lm4_cap_mod

  !-----------------------------------------------------------------------------
  ! LND Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC                  , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                  , only : NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model            , only : model_routine_SS           => SetServices
  use NUOPC_Model            , only : SetVM
  use NUOPC_Model            , only : model_label_Advance        => label_Advance
  use NUOPC_Model            , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model            , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model            , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model            , only : NUOPC_ModelGet

  use shr_kind_mod           , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use nuopc_shr_methods      , only : chkerr
  use lnd_import_export      , only : advertise_fields, realize_fields
  use fms_mod                , only: fms_init
  
  implicit none
  private ! except

  ! Module public routines
  public  :: SetServices
  public  :: SetVM

  ! Module private routines
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize

  character(len=CL)      :: flds_scalar_name = ''
  integer                :: flds_scalar_num = 0

  character(*),parameter :: modName =  "(lm4_cap_mod)"
  character(len=*) , parameter :: u_FILE_u =  __FILE__

  type(ESMF_GeomType_Flag) :: geomtype
  
!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc


    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)


  end subroutine SetServices


  !===============================================================================
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)

  end subroutine InitializeP0

  !===============================================================================

  !===============================================================================
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)      :: vm
    integer            :: lmpicom
    integer            :: ierr
    integer            :: n
    integer            :: localpet
    character(len=CL)  :: cvalue
    character(len=CL)  :: logmsg
    logical            :: isPresent, isSet
    logical            :: cism_evolve
    integer :: mype, ntasks, mpi_comm_land, mpi_comm_land2
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    
    !-------------------------------------------------------------------------------
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get communicator for fms. If smae proc layout as Atm, will be same communicator
    ! that Atm uses in it's fms_init. But it's ok, this fms_init will return without
    ! doing anything if already called on same proc layout
    call ESMF_VMGetCurrent(vm=VM,rc=RC)
    call ESMF_VMGet(vm=VM, localPet=mype, mpiCommunicator=mpi_comm_land, &
         petCount=ntasks, rc=rc)
    if (mype == 0) write(0,*) 'in lnd comp initadvert, ntasks=',ntasks
    !
    call fms_init(mpi_comm_land)

    !----------------------------------------------------------------------------
    ! advertise fields
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else
    !    call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else
    !    call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call advertise_fields(gcomp, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
  end subroutine InitializeAdvertise


  !===============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    use proc_bounds, only : procbounds, control_init_type
    use fms_io_mod,         only: field_exist, read_data

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    ! cube sphere mosaic                                                                                                                              
    type(ESMF_Decomp_Flag)  :: decompflagPTile(2,6)
    character(256)          :: gridfile
    integer                 :: tl
    integer,dimension(2,6)  :: decomptile                  !define delayout for the 6 cubed-sphere tiles                                              
    type(ESMF_Grid)         :: lndGrid
    character(50)           :: gridchoice
    type (control_init_type)::   ctrl_init

    !! tmp debug
    integer :: mype  
    type(ESMF_VM)      :: vm
    !-------------------------------------------------------------------------------   
    
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !! TMP DEBUG
    call ESMF_VMGetCurrent(vm=VM,rc=RC)
    call ESMF_VMGet(vm=VM, localPet=mype, rc=rc)

    call init_driver(ctrl_init)
    
    geomtype = ESMF_GEOMTYPE_GRID

    gridfile = "grid_spec.nc" ! default
    if (field_exist("INPUT/grid_spec.nc", "atm_mosaic_file")) then
       call read_data("INPUT/grid_spec.nc", "atm_mosaic_file", gridfile)
    endif

    ctrl_init%layout = (/2,4/) !!! TMP DEBUG
    do tl=1,6
       decomptile(1,tl) = ctrl_init%layout(1)
       decomptile(2,tl) = ctrl_init%layout(2)
       decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
    enddo

    if (mype == 0) write(0,*) 'JP DEBUG 5' !!!!
    if (mype == 0) write(0,*) 'JP gridfile', trim(gridfile)
    if (mype == 0) write(0,*) 'JP decomptile', decomptile
    !if (mype == 0) write(0,*) 'JP decompflagPTile', decompflagPTile(1,1)

    lndGrid = ESMF_GridCreateMosaic(filename="INPUT/"//trim(gridfile),     &
         regDecompPTile=decomptile,tileFilePath="INPUT/",                   &
         decompflagPTile=decompflagPTile,                                   &
         staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
         name='lnd_grid', rc=rc)

    if (mype == 0) write(0,*) 'JP DEBUG 6' !!!!
    ! ---------------------
    ! Realize the actively coupled fields
    ! ---------------------
    call realize_fields(gcomp, grid=lndGrid, flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mype == 0) write(0,*) 'JP DEBUG 7' !!!!

  end subroutine InitializeRealize

  !===============================================================================
  subroutine ModelAdvance(gcomp, rc)
    
    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! call lm4_finalize()
    
  end subroutine ModelFinalize


end module
