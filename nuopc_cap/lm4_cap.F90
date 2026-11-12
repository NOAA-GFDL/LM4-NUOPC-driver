module lm4_cap_mod

   !-----------------------------------------------------------------------------
   ! LM4 Component
   !-----------------------------------------------------------------------------

   use ESMF,                 only: ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSetEntryPoint
   use ESMF,                 only: ESMF_State, ESMF_Clock, ESMF_ClockGet, ESMF_VM, ESMF_VMGet, ESMF_VMGetCurrent
   use ESMF,                 only: ESMF_Time, ESMF_TimeGet, ESMF_TimeInterval, ESMF_TimeIntervalGet
   use ESMF,                 only: ESMF_Config, ESMF_ConfigCreate, ESMF_ConfigLoadFile, ESMF_ConfigGetAttribute
   use ESMF,                 only: ESMF_Grid, ESMF_GridCreateMosaic, ESMF_Decomp_Flag, ESMF_DECOMP_SYMMEDGEMAX
   use ESMF,                 only: ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, ESMF_METHOD_INITIALIZE
   use ESMF,                 only: ESMF_LogWrite, ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_LOGERR_PASSTHRU
   use ESMF,                 only: ESMF_SUCCESS, operator(-)
   use NUOPC,                only: NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
   use NUOPC,                only: NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet
   use NUOPC_Model,          only: model_routine_SS           => SetServices
   use NUOPC_Model,          only: SetVM
   use NUOPC_Model,          only: model_label_Advance        => label_Advance
   use NUOPC_Model,          only: model_label_Finalize       => label_Finalize
   use NUOPC_Model,          only: NUOPC_ModelGet

   use lm4_kind_mod,         only: r8 => shr_kind_r8, cl=>shr_kind_cl
   use lm4_type_mod,         only: lm4_type, alloc_atmforc, dealloc_atmforc, &
                                   alloc_atmsfc, dealloc_atmsfc

   use nuopc_lm4_methods,    only: chkerr
   use lm4_import_export,    only: advertise_fields, realize_fields, &
                                   import_fields, correct_import_fields, export_fields
   use fms_mod,              only: fms_init, fms_end, uppercase

   use mpp_mod,              only: mpp_error, FATAL
   use diag_manager_mod,     only: diag_manager_init, diag_manager_end, &
                                      diag_manager_set_time_end

   use lm4_driver,           only: lm4_nml_read, init_driver, end_driver, write_int_restart

   use land_model_mod,       only: land_model_init, land_model_end
   use constants_mod,        only: constants_init
   use lm4_surface_flux_mod, only: lm4_surface_flux_init

   use time_manager_mod,     only: set_calendar_type, set_date, set_time, get_time, &
                                   THIRTY_DAY_MONTHS, JULIAN, GREGORIAN, NOLEAP, NO_CALENDAR
   use sat_vapor_pres_mod,   only: sat_vapor_pres_init

   implicit none
   private ! except

   !---- model defined-types ----

   type(lm4_type) :: lm4_model

   integer :: date_init(6)

   ! Module public routines
   public  :: SetServices
   public  :: SetVM

   ! Module private routines
   private :: InitializeP0
   private :: InitializeAdvertise
   private :: InitializeRealize
   private :: ModelAdvance
   private :: ModelFinalize

   character(*),parameter       :: modName =  "(lm4_cap_mod)"
   character(len=*) , parameter :: u_FILE_u =  __FILE__

   ! for connecting ESMF and FMS timestep objects
   type(ESMF_TimeInterval)  :: model_timestep
   integer                  :: timestep_sec
   
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
      type(ESMF_GridComp)         :: gcomp
      type(ESMF_State)            :: importState, exportState
      type(ESMF_Clock)            :: clock  ! model clock

      integer, intent(out)        :: rc


      ! local variables
      type(ESMF_VM)               :: vm

      type(ESMF_Time)             :: CurrTime, StartTime, StopTime
      type(ESMF_TimeInterval)     :: RunDuration
      type(ESMF_Config)           :: cf
      integer                     :: Run_length
      integer,dimension(6)        :: date, date_end

      character(17)               :: calendar='                 '
      integer                     :: calendar_type = -99
      integer                     :: localpet

      character(len=CL)           :: logmsg


      integer :: mpi_comm_land
      character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '


      !-------------------------------------------------------------------------------



      rc = ESMF_SUCCESS

      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! Get communicator for fms. If same proc layout as Atm, will be same communicator
      ! that Atm uses in it's fms_init. But it's ok, this fms_init will return without
      ! doing anything if already called on same proc layout
      call ESMF_VMGetCurrent(vm=VM,rc=RC)
      call ESMF_VMGet(vm=VM, mpiCommunicator=mpi_comm_land, rc=rc)

      call fms_init(mpi_comm_land)

      call constants_init
      call sat_vapor_pres_init

      call lm4_surface_flux_init

      !------------------------------------------------------------------------
      ! get config variables
      !
      CF = ESMF_ConfigCreate(rc=rc)
      call ESMF_ConfigLoadFile(config=CF ,filename='model_configure' ,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      !
      call ESMF_ConfigGetAttribute(config=CF,value=calendar, &
         label ='calendar:', &
         default='gregorian',rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      ! Read lm4 namelist into lm4_model%nml
      call lm4_nml_read(lm4_model)



      !----------------------------------------------------------------------------
      ! Setting up clock, mirroring behavior in UFS FV3 fcst_grid_comp
      !----------------------------------------------------------------------------

      select case( uppercase(trim(calendar)) )
       case( 'JULIAN' )
         calendar_type = JULIAN
       case( 'GREGORIAN' )
         calendar_type = GREGORIAN
       case( 'NOLEAP' )
         calendar_type = NOLEAP
       case( 'THIRTY_DAY' )
         calendar_type = THIRTY_DAY_MONTHS
       case( 'NO_CALENDAR' )
         calendar_type = NO_CALENDAR
       case default
         call mpp_error ( FATAL, 'calendar must be one of '// &
            'JULIAN|GREGORIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
      end select

      call set_calendar_type (calendar_type)

      call ESMF_ClockGet(clock, CurrTime=CurrTime, StartTime=StartTime, &
         StopTime=StopTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      RunDuration = StopTime - CurrTime

      date_init = 0
      call ESMF_TimeGet (StartTime,                      &
         YY=date_init(1), MM=date_init(2), DD=date_init(3), &
         H=date_init(4),  M =date_init(5), S =date_init(6), RC=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if ( date_init(1) == 0 ) date_init = date
      lm4_model%Time_init  = set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))
      write(logmsg,*) date_init
      call ESMF_LogWrite(trim(subname)//'Land StartTime = '//trim(logmsg), ESMF_LOGMSG_INFO)

      date=0
      call ESMF_TimeGet (CurrTime,                           &
         YY=date(1), MM=date(2), DD=date(3), &
         H=date(4),  M =date(5), S =date(6), RC=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      write(logmsg,*) date
      call ESMF_LogWrite(trim(subname)//'Land CurrTime = '//trim(logmsg), ESMF_LOGMSG_INFO)

      lm4_model%Time_land = set_date (date(1), date(2), date(3),  &
         date(4), date(5), date(6))

      date_end=0
      call ESMF_TimeGet (StopTime,                                       &
         YY=date_end(1), MM=date_end(2), DD=date_end(3), &
         H=date_end(4),  M =date_end(5), S =date_end(6), RC=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if ( date_end(1) == 0 ) date_end = date
      lm4_model%Time_end   = set_date (date_end(1), date_end(2), date_end(3),  &
         date_end(4), date_end(5), date_end(6))
      write(logmsg,*) date_end
      call ESMF_LogWrite(trim(subname)//'Land StopTime = '//trim(logmsg), ESMF_LOGMSG_INFO)

      CALL ESMF_TimeIntervalGet(RunDuration, S=Run_length, RC=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call diag_manager_init (TIME_INIT=date)
      call diag_manager_set_time_end(lm4_model%Time_end)

      ! land time step from the ESMF clock (not dt_atmos, which doesn't apply to a data atmosphere)


  

      call ESMF_ClockGet(clock,  timeStep=model_timestep, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_TimeIntervalGet(model_timestep, s=timestep_sec, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return    
      lm4_model%dt_secs = timestep_sec
      write(logmsg,*) lm4_model%dt_secs
      call ESMF_LogWrite(trim(subname)//'init LM4 fast timestep: '//trim(logmsg), ESMF_LOGMSG_INFO)

      write(logmsg,*) lm4_model%nml%dt_lnd_slow
      call ESMF_LogWrite(trim(subname)//'init LM4 slow timestep: '//trim(logmsg), ESMF_LOGMSG_INFO)

      lm4_model%Time_step_land = set_time(timestep_sec,0)
      lm4_model%Time_step_slow = set_time(lm4_model%nml%dt_lnd_slow,0)
   
      
      !----------------------------------------------------------------------------
      ! Initialize model
      !----------------------------------------------------------------------------

      call land_model_init( lm4_model%From_atm, lm4_model%From_lnd, &
      lm4_model%Time_init, lm4_model%Time_land,                &
      lm4_model%Time_step_land, lm4_model%Time_step_slow     )
      
      ! allocate storage for data to/from nuopc-cap atm
      ! (before init_driver, which initializes atm_forc%gust)
      call alloc_atmforc(lm4_model%atm_forc)
      if (lm4_model%nml%cpl2atm) then  ! if have active 2-way coupling with atm
         call alloc_atmsfc(lm4_model%atm_sfc)
      end if

      call init_driver(lm4_model)

      call ESMF_LogWrite('======== COMPLETED land_model_init ==========', ESMF_LOGMSG_INFO)

      !----------------------------------------------------------------------------
      ! advertise fields
      !----------------------------------------------------------------------------

      call advertise_fields(gcomp, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_LogWrite(subname//' finished', ESMF_LOGMSG_INFO)

   end subroutine InitializeAdvertise


   !===============================================================================
   subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

      use fms2_io_mod,        only: read_data, open_file, close_file, FmsNetcdfFile_t, variable_exists

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

      character(len=CL)           :: cvalue
      character(len=CL)           :: logmsg
      logical                     :: isPresent, isSet

      type(FmsNetcdfFile_t)       :: fileobj

      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      gridfile = "grid_spec.nc" ! default

      if (open_file(fileobj, "INPUT/grid_spec.nc", "read")) then
         if (variable_exists(fileobj, "atm_mosaic_file")) then
            call read_data(fileobj, "atm_mosaic_file", gridfile)
         endif
         call close_file(fileobj)
      endif


      do tl=1,6
         decomptile(1,tl) = lm4_model%nml%layout(1)
         decomptile(2,tl) = lm4_model%nml%layout(2)
         decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
      enddo

      lndGrid = ESMF_GridCreateMosaic(filename="INPUT/"//trim(gridfile),     &
         regDecompPTile=decomptile,tileFilePath="INPUT/",                   &
         decompflagPTile=decompflagPTile,                                   &
         staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
         name='lnd_grid', rc=rc)

      ! set cpl_scalars from config. Default to null values
      lm4_model%cpl_scalar%flds_scalar_name = ''
      lm4_model%cpl_scalar%flds_scalar_num = 0
      lm4_model%cpl_scalar%flds_scalar_index_nx = 0
      lm4_model%cpl_scalar%flds_scalar_index_ny = 0
      lm4_model%cpl_scalar%flds_scalar_index_ntile = 0

      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
         lm4_model%cpl_scalar%flds_scalar_name = trim(cvalue)
         call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//lm4_model%cpl_scalar%flds_scalar_name, ESMF_LOGMSG_INFO)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
         read(cvalue, *) lm4_model%cpl_scalar%flds_scalar_num
         write(logmsg,*) lm4_model%cpl_scalar%flds_scalar_num
         call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
         read(cvalue, *) lm4_model%cpl_scalar%flds_scalar_index_nx
         write(logmsg,*) lm4_model%cpl_scalar%flds_scalar_index_nx
         call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
         read(cvalue, *) lm4_model%cpl_scalar%flds_scalar_index_ny
         write(logmsg,*) lm4_model%cpl_scalar%flds_scalar_index_ny
         call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif
      call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNTile", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet) then
         read(cvalue, *) lm4_model%cpl_scalar%flds_scalar_index_ntile
         write(logmsg,*) lm4_model%cpl_scalar%flds_scalar_index_ntile
         call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ntile = '//trim(logmsg), ESMF_LOGMSG_INFO)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
      endif

      call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (isPresent .and. isSet .and. trim(cvalue) == 'continue') then
         lm4_model%control%restart = .true.
      else
         lm4_model%control%restart = .false.
      endif

      ! ------------------------------------
      ! Realize the actively coupled fields
      ! ------------------------------------
      call realize_fields(gcomp, lm4_model, grid=lndGrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      
      
      ! ---------------------
      ! Create export state
      ! ---------------------
      call export_fields(gcomp, lm4_model, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_LogWrite(subname//' finished', ESMF_LOGMSG_INFO)

   end subroutine InitializeRealize

   !===============================================================================
   subroutine ModelAdvance(gcomp, rc)

      use lm4_driver,           only: sfc_boundary_layer, update_atmos_model_down, &
                                      flux_down_from_atmos, flux_up_to_atmos
      use land_model_mod,       only: update_land_model_fast, update_land_model_slow

      ! Arguments
      type(ESMF_GridComp)  :: gcomp
      type(ESMF_State)     :: importState, exportState
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_Clock)           :: clock
      type(ESMF_Clock)           :: dclock ! driver clock
      type(ESMF_Time)            :: currTime
      
      character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
      integer :: sec
      integer,dimension(6)        :: currdate ! for FMS time

      integer                  :: time_sec   ! current time in seconds
      type(ESMF_TimeInterval)  :: model_time ! current time as ESMF_TimeInterval
      character(len=CL)        :: logmsg

      rc = ESMF_SUCCESS
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)


      !-------------------------------------------------------------------------------
      ! Get import and export states
      !-------------------------------------------------------------------------------
      call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, modelClock=clock, driverClock=dclock, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      !-------------------------------------------------------------------------------
      ! Get import fields
      !-------------------------------------------------------------------------------
      call import_fields(gcomp, lm4_model, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call correct_import_fields(gcomp, lm4_model, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return


      !-------------------------------------------------------------------------------
      ! Run fast LM4 calls
      !-------------------------------------------------------------------------------

      call get_time (lm4_model%Time_step_land, sec)        ! get seconds of timestep
      call sfc_boundary_layer(real(sec), lm4_model)
      call update_atmos_model_down(lm4_model)              ! for gust calculation with data atmosphere
      call flux_down_from_atmos(real(sec), lm4_model)      ! explicit coupling only, implicit_atm not functional
      call update_land_model_fast(lm4_model%From_atm,lm4_model%From_lnd)

      if (lm4_model%nml%cpl2atm) then  ! if have active 2-way coupling with atm
         call flux_up_to_atmos(lm4_model)
      end if



      call ESMF_ClockGet(dclock,  CurrTime=CurrTime, currSimTime=model_time, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_TimeIntervalGet(model_time, s=time_sec, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! sync FMS time with ESMF time
      call ESMF_TimeGet (CurrTime,                           &
         YY=currdate(1), MM=currdate(2), DD=currdate(3), &
         H=currdate(4),  M =currdate(5), S =currdate(6), RC=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      write(logmsg,*) currdate
      call ESMF_LogWrite(trim(subname)//'Land CurrTime = '//trim(logmsg), ESMF_LOGMSG_INFO)

      lm4_model%Time_land = set_date (currdate(1), currdate(2), currdate(3),  &
         currdate(4), currdate(5), currdate(6))      

      
      !-------------------------------------------------------------------------------
      ! Run slow timescale LM4 calls
      !-------------------------------------------------------------------------------         
      ! replicates behavior in FMS coupler
      if ( mod(time_sec+timestep_sec ,lm4_model%nml%dt_lnd_slow) == 0 ) then
         call update_land_model_slow(lm4_model%From_atm,lm4_model%From_lnd)
         call ESMF_LogWrite(trim(subname)//'LM4 update_land_model_slow called', ESMF_LOGMSG_INFO)
         
         call write_int_restart(lm4_model)
      endif

      !-------------------------------------------------------------------------------
      ! Send export fields
      !-------------------------------------------------------------------------------      
      call export_fields(gcomp, lm4_model, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return      

      call ESMF_LogWrite(subname//' finished', ESMF_LOGMSG_INFO)

   end subroutine ModelAdvance

   !===============================================================================

   subroutine ModelFinalize(gcomp, rc)
      ! Arguments
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_Clock)     :: dclock
      character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
      type(ESMF_Time)            :: currTime
      
      character(len=CL)        :: logmsg
      integer,dimension(6)        :: currdate ! for FMS time

      rc = ESMF_SUCCESS
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      call NUOPC_ModelGet(gcomp, driverClock=dclock,rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_ClockGet(dclock,  CurrTime=CurrTime, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! sync FMS time with ESMF time
      call ESMF_TimeGet (CurrTime,                           &
         YY=currdate(1), MM=currdate(2), DD=currdate(3), &
         H=currdate(4),  M =currdate(5), S =currdate(6), RC=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      write(logmsg,*) currdate
      call ESMF_LogWrite(trim(subname)//'Land CurrTime = '//trim(logmsg), ESMF_LOGMSG_INFO)

      lm4_model%Time_land = set_date (currdate(1), currdate(2), currdate(3),  &
         currdate(4), currdate(5), currdate(6))      

      call write_int_restart(lm4_model)
      call end_driver()
      call land_model_end(lm4_model%From_atm, lm4_model%From_lnd)

      call diag_manager_end(lm4_model%Time_end)
      call fms_end    

      ! deallocate storage for data to/from nuopc-cap atm
      call dealloc_atmforc(lm4_model%atm_forc)
      
      if (lm4_model%nml%cpl2atm) then  ! if have active 2-way coupling with atm
         call dealloc_atmsfc(lm4_model%atm_sfc)
      end if      

      call ESMF_LogWrite(subname//' finished', ESMF_LOGMSG_INFO)

   end subroutine ModelFinalize

   !
   !----------------------------------------------------------------------------
end module
