module lm4_cap_mod

   !-----------------------------------------------------------------------------
   ! LM4 Component
   !-----------------------------------------------------------------------------

   use ESMF                  ! TODO: limit to only what is needed
   use NUOPC,                only: NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
   use NUOPC,                only: NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
   use NUOPC_Model,          only: model_routine_SS           => SetServices
   use NUOPC_Model,          only: SetVM
   use NUOPC_Model,          only: model_label_Advance        => label_Advance
   use NUOPC_Model,          only: model_label_DataInitialize => label_DataInitialize
   use NUOPC_Model,          only: model_label_SetRunClock    => label_SetRunClock
   use NUOPC_Model,          only: model_label_Finalize       => label_Finalize
   use NUOPC_Model,          only: NUOPC_ModelGet

   use lm4_kind_mod,         only: r8 => shr_kind_r8, cl=>shr_kind_cl
   use lm4_type_mod,         only: lm4_type, alloc_atmforc, dealloc_atmforc
   use lm4_type_mod,         only: alloc_atmforc2d, dealloc_atmforc2d ! TMP DEBUG

   use nuopc_lm4_methods,    only: chkerr
   use lnd_import_export,    only: advertise_fields, realize_fields, &
                                   import_fields, correct_import_fields, export_fields
   use fms_mod,              only: fms_init, fms_end, uppercase
   use fms_io_mod,           only: fms_io_exit

   use mpp_mod,              only: mpp_error,FATAL, WARNING
   use diag_manager_mod,     only: diag_manager_init, diag_manager_end, &
                                      diag_manager_set_time_end

   use lm4_driver,           only: lm4_nml_read, init_driver, end_driver, debug_diag, write_int_restart

   use land_model_mod,       only: land_model_init, land_model_end
   use land_data_mod,        only: land_data_type, atmos_land_boundary_type, lnd
   use constants_mod,        only: constants_init
   !use monin_obukhov_mod,    only: monin_obukhov_init
   use lm4_surface_flux_mod, only: lm4_surface_flux_init

   use time_manager_mod,     only: time_type, set_calendar_type, set_date,    &
                                   set_time, get_time,                        &
                                   THIRTY_DAY_MONTHS, JULIAN, GREGORIAN,      &
                                   NOLEAP, NO_CALENDAR
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

   character(len=CL)      :: flds_scalar_name = ''
   integer                :: flds_scalar_num = 0

   character(*),parameter       :: modName =  "(lm4_cap_mod)"
   character(len=*) , parameter :: u_FILE_u =  __FILE__

   !type(ESMF_GeomType_Flag) :: geomtype

   ! LM4 debug level
   integer :: debug_cap 

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
      type(ESMF_Clock)            :: dclock ! driver clock

      integer, intent(out)        :: rc


      ! local variables
      type(ESMF_VM)               :: vm
      integer                     :: lmpicom

      type(ESMF_Time)             :: CurrTime, StartTime, StopTime
      type(ESMF_TimeInterval)     :: RunDuration
      type(ESMF_Config)           :: cf
      integer                     :: Run_length
      integer,dimension(6)        :: date, date_end
      !integer                     :: dt_atmos

      character(17)               :: calendar='                 '
      integer                     :: calendar_type = -99
      integer                     :: ierr
      integer                     :: n
      integer                     :: localpet
      character(len=CL)           :: cvalue
      character(len=CL)           :: logmsg
      logical                     :: isPresent, isSet
      logical                     :: cism_evolve
      integer :: mype, ntasks, mpi_comm_land, mpi_comm_land2
      character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '

      ! for connecting ESMF and FMS timestep objects
      type(ESMF_TimeInterval)  :: model_timestep
      integer                  :: timestep_sec

      !-------------------------------------------------------------------------------



      rc = ESMF_SUCCESS

      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! Get communicator for fms. If same proc layout as Atm, will be same communicator
      ! that Atm uses in it's fms_init. But it's ok, this fms_init will return without
      ! doing anything if already called on same proc layout
      call ESMF_VMGetCurrent(vm=VM,rc=RC)
      call ESMF_VMGet(vm=VM, localPet=mype, mpiCommunicator=mpi_comm_land, &
         petCount=ntasks, rc=rc)

      call fms_init(mpi_comm_land)

      call constants_init
      call sat_vapor_pres_init

      call lm4_surface_flux_init
      ! originally in surface_flux_init
      !call monin_obukhov_init

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
      debug_cap = lm4_model%nml%lm4_debug

      ! if lm4_model%nml%lm4_debug is set, and > 0, write out namelist variables read in 
      if (mype == 0 .and. debug_cap > 0) then
         write(*,*) 'lm4_model%nml%lm4_debug: ' ,lm4_model%nml%lm4_debug
         write(*,*) 'lm4_model%nml%grid: '      ,lm4_model%nml%grid
         write(*,*) 'lm4_model%nml%npx: '       ,lm4_model%nml%npx
         write(*,*) 'lm4_model%nml%npy: '       ,lm4_model%nml%npy
         write(*,*) 'lm4_model%nml%layout: '    ,lm4_model%nml%layout
         write(*,*) 'lm4_model%nml%ntiles: '    ,lm4_model%nml%ntiles
         write(*,*) 'lm4_model%nml%blocksize: ' ,lm4_model%nml%blocksize
         write(*,*) 'lm4_model%nml%dt_lnd_slow ',lm4_model%nml%dt_lnd_slow
      endif



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
      if(mype==0) write(*,'(A,6I5)') 'Land StartTime=',date_init

      date=0
      call ESMF_TimeGet (CurrTime,                           &
         YY=date(1), MM=date(2), DD=date(3), &
         H=date(4),  M =date(5), S =date(6), RC=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      if(mype==0) write(*,'(A,6I5)') 'Land CurrTime =',date

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
      if(mype==0) write(*,'(A,6I5)') 'Land StopTime =',date_end

      ! call diag_manager_set_time_end(lm4_model%Time_end)

      CALL ESMF_TimeIntervalGet(RunDuration, S=Run_length, RC=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      call diag_manager_init (TIME_INIT=date)
      call diag_manager_set_time_end(lm4_model%Time_end)

      ! ! This will match the active atmosphere component's time step, but not appropriate for data atmosphere
      ! call ESMF_ConfigGetAttribute(config=CF, value=dt_atmos, label ='dt_atmos:',   rc=rc)
      ! lm4_model%Time_step_land = set_time(dt_atmos,0)
      ! if(mype==0) write(*,*) 'Land dt_atmos =',dt_atmos
      ! ! instead, get from ESMF clock


      ! get model and driver clocks
      ! call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return      
      ! call NUOPC_ModelGet(gcomp, driverClock=dclock, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
  

      call ESMF_ClockGet(clock,  timeStep=model_timestep, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_TimeIntervalGet(model_timestep, s=timestep_sec, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return    
      write(logmsg,*) timestep_sec
      call ESMF_LogWrite(trim(subname)//'init LM4 fast timestep: '//trim(logmsg), ESMF_LOGMSG_INFO)

      write(logmsg,*) lm4_model%nml%dt_lnd_slow
      call ESMF_LogWrite(trim(subname)//'init LM4 slow timestep: '//trim(logmsg), ESMF_LOGMSG_INFO)

      ! ! JP TMP DEBUG
      ! call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return       
      ! call NUOPC_ModelGet(gcomp, driverClock=dclock, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! call ESMF_ClockGet(dclock,  timeStep=model_timestep, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! call ESMF_TimeIntervalGet(model_timestep, s=timestep_sec, rc=rc)
      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
      ! write(logmsg,*) timestep_sec
      ! call ESMF_LogWrite(trim(subname)//'init LM4 driver timestep: '//trim(logmsg), ESMF_LOGMSG_INFO)      
      ! ! END TMP DEBUG

      lm4_model%Time_step_land = set_time(timestep_sec,0)
      lm4_model%Time_step_slow = set_time(lm4_model%nml%dt_lnd_slow,0)
   


      
      !----------------------------------------------------------------------------
      ! Initialize model
      !----------------------------------------------------------------------------

      call land_model_init( lm4_model%From_atm, lm4_model%From_lnd, &
      lm4_model%Time_init, lm4_model%Time_land,                &
      lm4_model%Time_step_land, lm4_model%Time_step_slow     )
      
      call init_driver(lm4_model)

      call ESMF_LogWrite('======== COMPLETED land_model_init ==========', ESMF_LOGMSG_INFO)

      ! allocate storage for the atm forc data
      call alloc_atmforc(lm4_model%atm_forc)
      if (debug_cap > 0) then
         call alloc_atmforc2d(lm4_model%atm_forc2d) ! TMP DEBUG
      endif

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

      call ESMF_LogWrite(subname//' finished', ESMF_LOGMSG_INFO)

   end subroutine InitializeAdvertise


   !===============================================================================
   subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

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


      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      !geomtype = ESMF_GEOMTYPE_GRID

      gridfile = "grid_spec.nc" ! default
      if (field_exist("INPUT/grid_spec.nc", "atm_mosaic_file")) then
         call read_data("INPUT/grid_spec.nc", "atm_mosaic_file", gridfile)
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

      !! JP TMP DEBUG
      if (debug_cap > 0) then
         call wrt_fcst_grid(lndGrid, "diagnostic_lndGrid.nc", rc=rc)
      endif


      ! ------------------------------------
      ! Realize the actively coupled fields
      ! ------------------------------------
      call realize_fields(gcomp, grid=lndGrid, flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! ---------------------
      ! Create export state
      ! ---------------------
      call export_fields(gcomp,rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_LogWrite(subname//' finished', ESMF_LOGMSG_INFO)

   end subroutine InitializeRealize

   !===============================================================================
   subroutine ModelAdvance(gcomp, rc)

      use lm4_driver,           only: sfc_boundary_layer, flux_down_from_atmos
      use land_model_mod,       only: update_land_model_fast, update_land_model_slow
      use ESMF, only: ESMF_ClockPrint, ESMF_AlarmIsRinging ! TMP DEBUG

      ! Arguments
      type(ESMF_GridComp)  :: gcomp
      type(ESMF_State)     :: importState, exportState
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_Clock)           :: clock
      type(ESMF_Clock)           :: dclock ! driver clock
      type(ESMF_Time)            :: prevTime
      type(ESMF_Time)            :: currTime
      type(ESMF_Alarm)           :: slowAlarm
      
      character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
      integer :: sec
      integer,dimension(6)        :: currdate ! for FMS time

      ! JP tmp debug
      integer                  :: time_sec   ! current time in seconds
      type(ESMF_TimeInterval)  :: model_time ! current time as ESMF_TimeInterval
      character(len=CL)        :: logmsg
      ! END TMP DEBUG

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


      ! option to write out diag history of imports
      ! if (debug_cap > 0) then
      !    call debug_diag(lm4_model)
      ! endif

      ! TODO, is sec being used anywhere?
      call get_time (lm4_model%Time_step_land, sec)  ! get seconds of timestep
      call sfc_boundary_layer(real(sec), lm4_model)
      call flux_down_from_atmos(real(sec), lm4_model)      ! JP: needs review of implicit coupling
      call update_land_model_fast(lm4_model%From_atm,lm4_model%From_lnd)


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

      ! JP TMP DEBUG
      write(*,*) 'ModelAdvance: clock info:'
      call ESMF_ClockPrint(dclock,rc=rc)
      write(logmsg,*) time_sec
      call ESMF_LogWrite(trim(subname)//'MA LM4 driver currSimTime: '//trim(logmsg), ESMF_LOGMSG_INFO)
      
      ! quick way to only call on slow timestep
      if ( time_sec /= 0 .and. mod(time_sec,lm4_model%nml%dt_lnd_slow) == 0 ) then
         call update_land_model_slow(lm4_model%From_atm,lm4_model%From_lnd)
         call ESMF_LogWrite('MA LM4 update_land_model_slow called', ESMF_LOGMSG_INFO)
         
         call write_int_restart(lm4_model)
      endif

      call ESMF_LogWrite(subname//' finished', ESMF_LOGMSG_INFO)

   end subroutine ModelAdvance

   !===============================================================================

   subroutine ModelFinalize(gcomp, rc)
      type(ESMF_GridComp)  :: gcomp
      integer, intent(out) :: rc
      character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '

      rc = ESMF_SUCCESS
      call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      call end_driver()
      call land_model_end(lm4_model%From_atm, lm4_model%From_lnd)

      call diag_manager_end(lm4_model%Time_end)
      call fms_io_exit  ! TEST
      call fms_end      ! TEST

      ! deallocate storage for the atm forc data
      call dealloc_atmforc(lm4_model%atm_forc)
      if (debug_cap > 0) then
         call dealloc_atmforc2d(lm4_model%atm_forc2d) ! TMP DEBUG
      endif

      call ESMF_LogWrite(subname//' finished', ESMF_LOGMSG_INFO)

   end subroutine ModelFinalize

   !
   !#######################################################################
   !-- TMP DEBUG write grid to NetCDF file for diagnostics
   !
   subroutine wrt_fcst_grid(grid, fileName, relaxedflag, regridArea, rc)
      type(ESMF_Grid), intent(in)                      :: grid
      character(len=*), intent(in), optional           :: fileName
      logical, intent(in), optional                    :: relaxedflag
      logical, intent(in), optional                    :: regridArea
      integer, intent(out)                             :: rc
      !
      !-----------------------------------------------------------------------
      !***  local variables
      !
      logical                     :: ioCapable
      logical                     :: doItFlag
      character(len=64)           :: lfileName
      character(len=64)           :: gridName
      type(ESMF_Array)            :: array
      type(ESMF_ArrayBundle)      :: arraybundle
      logical                     :: isPresent
      integer                     :: stat
      logical                     :: hasCorners
      logical                     :: lRegridArea
      type(ESMF_Field)            :: areaField
      type(ESMF_FieldStatus_Flag) :: areaFieldStatus

      ioCapable = (ESMF_IO_PIO_PRESENT .and. &
         (ESMF_IO_NETCDF_PRESENT .or. ESMF_IO_PNETCDF_PRESENT))
      doItFlag = .true.
      if (present(relaxedFlag)) then
         doItFlag = .not.relaxedflag .or. (relaxedflag.and.ioCapable)
      endif

      if (doItFlag) then
         ! Process optional arguments
         if (present(fileName)) then
            lfileName = trim(fileName)
         else
            call ESMF_GridGet(grid, name=gridName, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            lfileName = trim(gridName)//".nc"
         endif
         if (present(regridArea)) then
            lRegridArea = regridArea
         else
            lRegridArea = .FALSE.
         endif

         ! Create bundle for storing output
         arraybundle = ESMF_ArrayBundleCreate(rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return

         ! -- Centers --
         call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
            isPresent=isPresent, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (isPresent) then
            call ESMF_GridGetCoord(grid, coordDim=1, &
               staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_ArraySet(array, name="lon_center", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_GridGetCoord(grid, coordDim=2, &
               staggerLoc=ESMF_STAGGERLOC_CENTER, array=array, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_ArraySet(array, name="lat_center", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
         endif

         ! -- Corners --
         call ESMF_GridGetCoord(grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
            isPresent=hasCorners, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (hasCorners) then
            call ESMF_GridGetCoord(grid, coordDim=1, &
               staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
            if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
               call ESMF_ArraySet(array, name="lon_corner", rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return
               call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return
            endif
            call ESMF_GridGetCoord(grid, coordDim=2, &
               staggerLoc=ESMF_STAGGERLOC_CORNER, array=array, rc=rc)
            if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
               call ESMF_ArraySet(array, name="lat_corner", rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return
               call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return
            endif
            if (lRegridArea) then
               areaField = ESMF_FieldCreate(grid=grid, &
                  typekind=ESMF_TYPEKIND_R8, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return
               call ESMF_FieldRegridGetArea(areaField, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return
               call ESMF_FieldGet(areaField, array=array, rc=rc)
               if (.not. ESMF_LogFoundError(rc, line=__LINE__, file=__FILE__)) then
                  call ESMF_ArraySet(array, name="regrid_area", rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                     line=__LINE__, file=__FILE__)) return
                  call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                     line=__LINE__, file=__FILE__)) return
               endif
            endif
         endif
         ! -- Mask --
         call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_MASK, &
            staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (isPresent) then
            call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
               itemflag=ESMF_GRIDITEM_MASK, array=array, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_ArraySet(array, name="mask", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
         endif

         ! -- Area --
         call ESMF_GridGetItem(grid, itemflag=ESMF_GRIDITEM_AREA, &
            staggerLoc=ESMF_STAGGERLOC_CENTER, isPresent=isPresent, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
         if (isPresent) then
            call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
               itemflag=ESMF_GRIDITEM_AREA, array=array, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_ArraySet(array, name="area", rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            call ESMF_ArrayBundleAdd(arraybundle,(/array/), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
         endif

         ! Write array bundle to grid file
         ! note: 6-tile not supported yet
         ! call ESMF_ArrayBundleWrite(arraybundle, fileName=trim(lfileName), rc=rc)
         ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         !      line=__LINE__, file=__FILE__)) return

         ! Clean-up
         if (lRegridArea) then
            call ESMF_FieldGet(areaField, status=areaFieldStatus, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) return
            if (areaFieldStatus.eq.ESMF_FIELDSTATUS_COMPLETE) then
               call ESMF_FieldDestroy(areaField, rc=rc)
               if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, file=__FILE__)) return
            endif
         endif
         call ESMF_ArrayBundleDestroy(arraybundle,rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) return
      endif
   end subroutine wrt_fcst_grid
   !
   !----------------------------------------------------------------------------
end module
