!! Routines to prepare and run the land model
!! ============================================================================
module lm4_driver

   use ESMF,                 only: ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_END_ABORT, ESMF_Finalize, &
                                   ESMF_LOGMSG_ERROR, ESMF_LOGMSG_WARNING
   use mpp_domains_mod,      only: domain2d
   use mpp_mod,              only: input_nml_file

   use lm4_type_mod,         only: lm4_type
   use lm4_kind_mod,         only: r8 => shr_kind_r8, cl=>shr_kind_cl
   use land_data_mod,        only: lnd
   use land_tile_mod,        only : land_tile_map, land_tile_enum_type, land_tile_type, first_elmt, &
                                    loop_over_tiles, max_n_tiles
   use land_tracers_mod,     only: isphum, ntcana
   use lm4_surface_flux_mod, only: lm4_surface_flux_1d, virtual_temp, air_density, alt_gustiness

   use time_manager_mod,      only: time_type, date_to_string, increment_date, decrement_date
   use time_manager_mod,      only: operator(>=), operator(<), operator(==)


   implicit none
   private

   type(domain2D), public :: land_domain

   integer :: date_init(6)
   character(len=32)     :: timestamp


   public :: lm4_nml_read
   public :: init_driver, end_driver
   public :: sfc_boundary_layer, update_atmos_model_down, flux_down_from_atmos, flux_up_to_atmos
   public :: write_int_restart


   ! --- namelist of vars originally from flux exchange nml
   real :: z_ref_heat =  2. !< Reference height (meters) for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q)
   real :: z_ref_mom  = 10. !< Reference height (meters) for momentum diagnostics (u_ref, v_ref, del_m)

   ! name kept from FMS so existing GFDL input.nml files work
   namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom

   ! --- namelist of vars originally from atmos_prescr_nml
   character(len=24) :: gust_to_use = 'computed' ! or 'prescribed'
   real    :: gustiness    = 5.0  ! m/s, wind gustiness if gust_to_use = 'prescribed'
   real    :: gust_min     = 0.0  ! m/s, minimum gustiness when gust_to_use = 'computed'
   namelist /atmos_prescr_nml/ gustiness, gust_to_use, gust_min


   logical :: do_forecast = .false.


   ! variables for between subroutines

   integer, allocatable :: tidx(:)     ! land tile indices for current domain
   integer              :: ntiles_lnd  ! number of land tiles in current domain
   integer              :: ntile_types ! number of land tile types in current domain
   integer              :: k           ! local loop index for land tiles
   integer              :: kk          ! global loop index for land tiles
   integer              :: l           ! land gridcell index
   integer              :: i, j        ! land gridcell indices

   type(land_tile_enum_type) :: ce ! tile list enumerator
   type(land_tile_type), pointer :: tile


   real, allocatable, dimension(:) :: &
      ex_flux_t, ex_flux_lw,      &
      ex_dhdt_surf, ex_dedt_surf, &
      ex_drdt_surf,  ex_dhdt_atm, &
      ex_drag_q,     &   !< q drag.coeff.
      ex_cd_t,       &
      ex_cd_m,       &
      ex_b_star,     &
      ex_u_star,     &
      ex_wind,       &
      ex_z_atm,      &
      ex_p_surf,     &
      ex_t_surf,     &
      ex_t_ca,       &
      ex_f_t_delt_n, &   
      ex_e_t_n        

   logical, allocatable, dimension(:) :: &
      ex_avail,     &   !< true where data on exchange grid are available
      ex_land           !< true if exchange grid cell is over land


   ! these originally had a tracer dimension
   real, allocatable, dimension(:,:) :: &
      ex_tr_atm,      &
      ex_tr_surf,     & !< near-surface tracer fields
      ex_flux_tr,     & !< tracer fluxes
      ex_dfdtr_surf,  & !< d(tracer flux)/d(surf tracer)
      ex_dfdtr_atm,   & !< d(tracer flux)/d(atm tracer)
      ex_e_tr_n,      & !< coefficient in implicit scheme
      ex_f_tr_delt_n    !< coefficient in implicit scheme


   !integer :: n_exch_tr !< number of tracers exchanged between models

   character(len=CL) :: logmsg

contains

   !! Read in lm4 namelist for NUOPC-cap related variables
   !! Also read in other namelists that are not read in from other model components,
   !! but are used in this module
   !! ============================================================================
   subroutine lm4_nml_read(lm4_model)

      use fms_mod,             only: check_nml_error
      use mpp_mod,             only: mpp_pe, mpp_root_pe, stdlog


      type(lm4_type),          intent(inout) :: lm4_model ! land model's variable type

      ! namelist variables for lm4 cap
      ! ------------------------------------------
      integer           :: lm4_debug   = 0          ! debug flag for lm4 (0=off, 1=low, 2=high)
      integer           :: npx         = 0, npy = 0
      integer           :: ntiles      = 0
      integer           :: layout(2)   = (/0,0/)
      character(len=64) :: grid        = 'none'
      integer           :: blocksize   = -1
      integer           :: dt_lnd_slow = 86400  ! time step for slow land processes (s)
      integer, dimension(6) :: restart_interval = (/ 0, 0, 0, 0, 0, 0/) !< The time interval that write out intermediate restart file.
      !! The format is (yr,mo,day,hr,min,sec).  When restart_interval
      !! is all zero, no intermediate restart file will be written out
      logical           :: cpl2atm           = .false.  ! coupling to active atmosphere 
      logical           :: kinematic_flux    = .true.   ! some fluxes to atm are expressed as kinematic
      logical           :: implicit_atm      = .false.  ! implicit coupling to atm (not functional, .true. aborts)

      ! for namelist read
      integer :: io, ierr
      namelist /lm4_nml/ grid, npx, npy, layout, ntiles,      &
         blocksize, lm4_debug, dt_lnd_slow, restart_interval, &
         cpl2atm, kinematic_flux, implicit_atm

      ! read in namelists
      ! ------------------------------------------
      ! lm4_nml
      read(input_nml_file, nml=lm4_nml, iostat=io)
      ierr = check_nml_error(io, 'lm4_nml')
      if (mpp_pe() == mpp_root_pe()) write(stdlog(), nml=lm4_nml)

      if (implicit_atm) then
         call ESMF_LogWrite('lm4_nml_read: implicit_atm = .true. is not functional', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif

      read(input_nml_file, nml=atmos_prescr_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_prescr_nml')

      read(input_nml_file, nml=flux_exchange_nml, iostat=io)
      ierr = check_nml_error(io, 'flux_exchange_nml')


      lm4_model%nml%lm4_debug        = lm4_debug
      lm4_model%nml%grid             = grid
      lm4_model%nml%blocksize        = blocksize
      lm4_model%nml%npx              = npx
      lm4_model%nml%npy              = npy
      lm4_model%nml%layout           = layout
      lm4_model%nml%ntiles           = ntiles
      lm4_model%nml%dt_lnd_slow      = dt_lnd_slow
      lm4_model%nml%restart_interval = restart_interval
      lm4_model%nml%cpl2atm          = cpl2atm
      lm4_model%nml%kinematic_flux   = kinematic_flux
      lm4_model%nml%implicit_atm     = implicit_atm

   end subroutine lm4_nml_read

   !! ============================================================================
   subroutine init_driver(lm4_model)

      use mpp_domains_mod,    only: domain2d, mpp_get_compute_domain
      use land_domain_mod,    only: domain_create
      use land_tile_io_mod,   only : gather_tile_index

      type(lm4_type),          intent(inout) :: lm4_model ! land model's variable type

      ! local
      integer :: isc, iec, jsc, jec

      ! FMS domain creation:
      call domain_create(lm4_model%nml, land_domain)

      call mpp_get_compute_domain(land_domain,isc,iec,jsc,jec)

      ! this gives land tiles same compressed indexing size as in LM4 IO
      call gather_tile_index(land_tile_exists, tidx)
      ntiles_lnd = size(tidx)
      write(logmsg, '(A,I4)') 'init_driver number of land tiles = ', ntiles_lnd
      call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)

      allocate( &
         ex_flux_t(ntiles_lnd), ex_flux_lw(ntiles_lnd),   &
         ex_dhdt_surf(ntiles_lnd), ex_dedt_surf(ntiles_lnd), &
         ex_drdt_surf(ntiles_lnd),  ex_dhdt_atm(ntiles_lnd), &
         ex_drag_q(ntiles_lnd),     &   !< q drag.coeff.
         ex_cd_t(ntiles_lnd),       &
         ex_cd_m(ntiles_lnd),       &
         ex_b_star(ntiles_lnd),     &
         ex_u_star(ntiles_lnd),     &
         ex_wind(ntiles_lnd),       &
         ex_z_atm(ntiles_lnd),      &
         ex_p_surf(ntiles_lnd),     &
         ex_t_surf(ntiles_lnd),     &
         ex_t_ca(ntiles_lnd),       &
         ex_f_t_delt_n(ntiles_lnd), &   
         ex_e_t_n(ntiles_lnd)       &        
         )

      ! these originally had a tracer dimension
      allocate( &
         ex_tr_atm(ntiles_lnd,ntcana),      &
         ex_tr_surf(ntiles_lnd,ntcana),     & !< near-surface tracer fields
         ex_flux_tr(ntiles_lnd,ntcana),     & !< tracer fluxes
         ex_dfdtr_surf(ntiles_lnd,ntcana),  & !< d(tracer flux)/d(surf tracer)
         ex_dfdtr_atm(ntiles_lnd,ntcana),   & !< d(tracer flux)/d(atm tracer)
         ex_e_tr_n(ntiles_lnd,ntcana),      & !< coefficient in implicit scheme
         ex_f_tr_delt_n(ntiles_lnd,ntcana), & !< coefficient in implicit scheme
         ex_avail(ntiles_lnd),              & !< true where data on exchange grid are available
         ex_land(ntiles_lnd)                & !< true if exchange grid cell is over land
         )

      ! initialize reals
      ex_flux_t     = 0.0_r8
      ex_flux_lw    = 0.0_r8
      ex_dhdt_surf  = 0.0_r8
      ex_dedt_surf  = 0.0_r8
      ex_drdt_surf  = 0.0_r8
      ex_dhdt_atm   = 0.0_r8
      ex_drag_q     = 0.0_r8
      ex_cd_t       = 0.0_r8
      ex_cd_m       = 0.0_r8
      ex_b_star     = 0.0_r8
      ex_u_star     = 0.0_r8
      ex_wind       = 0.0_r8
      ex_z_atm      = 0.0_r8
      ex_p_surf     = 0.0_r8
      ex_t_surf     = 0.0_r8
      ex_t_ca       = 0.0_r8
      ex_f_t_delt_n = 0.0_r8
      ex_e_t_n      = 0.0_r8

      ex_tr_atm      = 0.0_r8
      ex_tr_surf     = 0.0_r8
      ex_flux_tr     = 0.0_r8
      ex_dfdtr_surf  = 0.0_r8
      ex_dfdtr_atm   = 0.0_r8
      ex_e_tr_n      = 0.0_r8
      ex_f_tr_delt_n = 0.0_r8

      ! initialize ex_avail and ex_land
      ex_avail    = .TRUE.
      ex_land     = .TRUE.

      ! Set restart time
      if (ALL(lm4_model%nml%restart_interval ==0)) then
         lm4_model%Time_restart = increment_date(lm4_model%Time_end, 0, 0, 10, 0, 0, 0)   ! no intermediate restart
      else

         lm4_model%Time_restart = increment_date(lm4_model%Time_land, lm4_model%nml%restart_interval(1), lm4_model%nml%restart_interval(2), &
            lm4_model%nml%restart_interval(3), lm4_model%nml%restart_interval(4), lm4_model%nml%restart_interval(5), lm4_model%nml%restart_interval(6) )

         ! subtract the fast time step in seconds
         timestamp = date_to_string(lm4_model%Time_restart)
         call ESMF_LogWrite('LM4 init_driver: Time_restart before decrement' //trim(timestamp), ESMF_LOGMSG_INFO)
         lm4_model%Time_restart = decrement_date(lm4_model%Time_restart, 0,0,0,0,0, lm4_model%dt_secs)
         timestamp = date_to_string(lm4_model%Time_restart)
         call ESMF_LogWrite('LM4 init_driver: Time_restart after decrement' //trim(timestamp), ESMF_LOGMSG_INFO)


         if (lm4_model%Time_restart < lm4_model%Time_land) then
            call ESMF_LogWrite('The first intermediate restart time is earlier than the current time', &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
            call ESMF_Finalize(endflag=ESMF_END_ABORT)

         endif

      endif

      ! initialize gust
      if (trim(gust_to_use)=='computed') then
         ! can't use with CDEPS atm, since no ustar/bstar provided in atm forcing or restarts
         call ESMF_LogWrite('Computed gustiness reinitializes gust on initialization with this data atmosphere.', &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=__FILE__)
         call ESMF_LogWrite('Restarting is NOT reproducible with continous run', &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=__FILE__)
      elseif (trim(gust_to_use)=='prescribed') then
         call ESMF_LogWrite('Using prescribed gustiness', ESMF_LOGMSG_INFO)
         write(logmsg, '(A,F6.2)') 'gustiness = ', gustiness
         call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)
         write(logmsg, '(A,F6.2)') 'gust_min = ', gust_min
         call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)
      endif
      ! computed: u_star = b_star = 0 here, so gust starts at gust_min
      call compute_gust(lm4_model)

      ! FV3 GFS physics floors its surface wind at 1 m/s, like alt_gustiness with gust_const = 1
      if (lm4_model%nml%cpl2atm .and. .not. alt_gustiness) then
         call ESMF_LogWrite('cpl2atm without alt_gustiness: land and atm use different gustiness', &
            ESMF_LOGMSG_WARNING, line=__LINE__, file=__FILE__)
      endif

      
      if (lm4_model%nml%lm4_debug > 0) then
         ntile_types = max_n_tiles()
         write(logmsg, '(A,I4)') 'LM4 init_driver max number of land tile types = ', ntile_types
         call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)
         call test_aggregate_land_tiles()
      endif

      
   end subroutine init_driver

   !! ============================================================================
   !! Adapted from GFDL coupler, write intermediate restarts
   !! ============================================================================
   subroutine write_int_restart(lm4_model)
      use land_model_mod,          only: land_model_restart

      type(lm4_type), intent(inout) :: lm4_model ! land model's variable type

      type(time_type)       :: Time_restart_stamp ! datetime stamp for restart file

      !--- write out intermediate restart file when needed.
      if (lm4_model%Time_land >= lm4_model%Time_restart) then
         lm4_model%Time_restart = increment_date(lm4_model%Time_land, lm4_model%nml%restart_interval(1), lm4_model%nml%restart_interval(2), &
            lm4_model%nml%restart_interval(3), lm4_model%nml%restart_interval(4), lm4_model%nml%restart_interval(5), lm4_model%nml%restart_interval(6) )

         ! To match behavior of GFDL coupler, advance the restart file timestamp by the fast time step
         ! (datestamp is modeltime needed for restart)
         Time_restart_stamp = increment_date(lm4_model%Time_land, 0,0,0,0,0, lm4_model%dt_secs)
         timestamp = date_to_string(Time_restart_stamp)
         call ESMF_LogWrite('write_int_restart restart is written for '//trim(timestamp), ESMF_LOGMSG_INFO)
         call land_model_restart(timestamp)
      endif

   end subroutine write_int_restart

   !! ============================================================================
   !! Adapted from GFDL atm_land_ice_flux_exchange,
   !! stripped down to be "land only" on
   !! unstructured grid, returns
   !! explicit fluxes as well as derivatives
   !! used to compute an implicit flux correction.
   !! ============================================================================
   subroutine sfc_boundary_layer( dt,lm4_model )


      use monin_obukhov_mod,  only: mo_profile
      use land_tile_mod,      only: max_n_tiles   ! this should be 1 for now

      real,                  intent(in)     :: dt        ! Time step
      type(lm4_type),        intent(inout)  :: lm4_model ! land model's variable type

      ! local variables --------------------------------------

      real    :: zrefm, zrefh

      real, dimension(ntiles_lnd) :: &
         ex_t_atm,              &
         ex_p_atm,              &
         ex_u_atm, ex_v_atm,    &
         ex_q_atm,              &
         ex_gust,               &
         ex_t_surf4,            &
         ex_u_surf,              &
         ex_v_surf,              &
         ex_rough_mom, ex_rough_heat, ex_rough_moist, &
         ex_rough_scale,        &
         ex_q_star,             &
         ex_cd_q,               &
         ex_ref_u, ex_ref_v, ex_u10, &
         ex_del_m,              &
         ex_del_h,              &
         ex_del_q


      !! -- These were originally allocatable:

      logical, dimension(ntiles_lnd) :: &
         ex_seawater       !< true if exchange grid cell is over seawater

      real, dimension(ntiles_lnd) :: &
         ex_flux_u,    &
         ex_flux_v,    &
         ex_dtaudu_atm,&
         ex_dtaudv_atm

      integer :: tr ! tracer index




      ! ---------------------------
      ex_seawater = .FALSE.

      ! these are 0 for land
      ex_u_surf = 0.0
      ex_v_surf = 0.0


      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0  ! global tile index
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))  
         kk = kk + 1      
         ! prefill q_surf with q_bot. In original code, there were options to send/write
         ! before surface_flux call.
         ex_tr_surf(kk,isphum) = lm4_model%atm_forc%q_bot(l)

         do tr = 1,ntcana
            ex_tr_surf(kk,tr) = lm4_model%From_lnd%tr(l,k,tr)
         enddo
      end do

      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0  ! global tile index
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))  
         kk = kk + 1

         ! fill in tile data from gridcell data
         ex_t_atm(kk)  = lm4_model%atm_forc%t_bot(l)
         ex_q_atm(kk)  = lm4_model%atm_forc%q_bot(l)
         ex_u_atm(kk)  = lm4_model%atm_forc%u_bot(l)
         ex_v_atm(kk)  = lm4_model%atm_forc%v_bot(l)
         ex_p_atm(kk)  = lm4_model%atm_forc%p_bot(l)
         ex_z_atm(kk)  = lm4_model%atm_forc%z_bot(l)
         ex_p_surf(kk) = lm4_model%atm_forc%p_surf(l)
         ex_gust(kk)   = lm4_model%atm_forc%gust(l)

         ! this is mimicking the original code for land roughness vars
         ex_rough_mom(kk)   = lm4_model%From_lnd%rough_mom(l,k)
         ex_rough_moist(kk) = lm4_model%From_lnd%rough_heat(l,k)
         ex_rough_heat(kk)  = lm4_model%From_lnd%rough_heat(l,k)
         ex_rough_scale(kk) = lm4_model%From_lnd%rough_scale(l,k)

         ex_t_surf(kk) = lm4_model%From_lnd%t_surf(l,k)
         ex_t_ca(kk)   = lm4_model%From_lnd%t_ca(l,k)         
      end do

      ! initialize other variables to zero
      ex_flux_t         =  0.0
      ex_flux_tr        =  0.0
      ex_flux_lw        =  0.0
      ex_flux_u         =  0.0
      ex_flux_v         =  0.0
      ex_cd_m           =  0.0
      ex_cd_t           =  0.0
      ex_cd_q           =  0.0
      ex_wind           =  0.0
      ex_u_star         =  0.0
      ex_b_star         =  0.0
      ex_q_star         =  0.0
      ex_dhdt_surf      =  0.0
      ex_dedt_surf      =  0.0
      ex_dfdtr_surf     =  0.0
      ex_drdt_surf      =  0.0
      ex_dhdt_atm       =  0.0
      ex_dfdtr_atm      =  0.0
      ex_dtaudu_atm     =  0.0
      ex_dtaudv_atm     =  0.0


      ! blocking not used for now
      call lm4_surface_flux_1d ( &
      ! inputs
         ex_t_atm, ex_q_atm,  & 
         ex_u_atm, ex_v_atm,  ex_p_atm, &
         ex_z_atm, ex_p_surf, ex_t_surf, &
         ex_t_ca, &
      ! inout
         ex_tr_surf(:,isphum),         & 
      ! more inputs
         ex_u_surf, ex_v_surf,             & ! 0s
         ex_rough_mom, ex_rough_heat, &
         ex_rough_moist, ex_rough_scale,   &
         ex_gust,                                                       & ! gustiness
      ! outputs
         ex_flux_t, ex_flux_tr(:,isphum), ex_flux_lw,   &
         ex_flux_u, ex_flux_v, ex_cd_m,   ex_cd_t, &
         ex_cd_q,   ex_wind,   ex_u_star, ex_b_star, &
         ex_q_star, ex_dhdt_surf, ex_dedt_surf, &
         ex_dfdtr_surf(:,isphum),  ex_drdt_surf,  ex_dhdt_atm, &
         ex_dfdtr_atm(:,isphum),  ex_dtaudu_atm, ex_dtaudv_atm,         &
      ! more inputs
         dt,                                                            & !! timestep doesn't seem to be used
         ex_land, ex_seawater, ex_avail                                 & ! Is land, Is seawater, Is ex. avail
         )



      !! ....
      zrefm = z_ref_mom
      zrefh = z_ref_heat
      !      ---- optimize calculation ----
      call mo_profile ( zrefm, zrefh, ex_z_atm, ex_rough_mom, &
         ex_rough_heat, ex_rough_moist,          &
         ex_u_star, ex_b_star, ex_q_star,        &
         ex_del_m, ex_del_h, ex_del_q, ex_avail  )

      ! Former OMP parallel do
      do k = 1, ntiles_lnd
         ex_u10(k) = 0.
         if(ex_avail(k)) then
            ex_ref_u(k) = ex_u_surf(k) + (ex_u_atm(k)-ex_u_surf(k)) * ex_del_m(k)
            ex_ref_v(k) = ex_v_surf(k) + (ex_v_atm(k)-ex_v_surf(k)) * ex_del_m(k)
            ex_u10(k) = sqrt(ex_ref_u(k)**2 + ex_ref_v(k)**2)
         endif
      enddo

      ! FMS also fills derivatives for tracers other than sphum here; not needed while sphum is the
      ! only exchanged tracer

      if (lm4_model%nml%cpl2atm) call surface_layer_exports(lm4_model, ex_del_h, ex_del_q, ex_rough_mom)
      
      do k = 1, ntiles_lnd
         if(ex_avail(k)) ex_drag_q(k) = ex_wind(k)*ex_cd_q(k)
         ! [6] get mean quantities on atmosphere grid
         ! [6.1] compute t surf for radiation
         ex_t_surf4(k) = ex_t_surf(k) ** 4
      enddo


      if (lm4_model%nml%cpl2atm ) then 
         ! accumulate tile-fraction-weighted albedos for each grid cell l
         ! ignore '_fix' albedos in original code, just send land albedos for export
         ! note SW down from atmosphere is grid cell value, so no need to weight
         ce = first_elmt(land_tile_map, ls=lnd%ls)
         kk = 0  ! global tile index
         lm4_model%atm_sfc%albedo_vis_dir = 0.0
         lm4_model%atm_sfc%albedo_nir_dir = 0.0
         lm4_model%atm_sfc%albedo_vis_dif = 0.0
         lm4_model%atm_sfc%albedo_nir_dif = 0.0

         ! check if associated properly
         if (.not. associated(lm4_model%atm_sfc%albedo_nir_dif)) then
            call ESMF_LogWrite('lm4_model%atm_sfc%albedo_nir_dif not associated in sfc_boundary_layer', ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         end if

         do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))  
            kk = kk + 1
            if (.not. land_tile_exists(tile)) cycle

            lm4_model%atm_sfc%albedo_vis_dir(l) = lm4_model%atm_sfc%albedo_vis_dir(l) + lm4_model%From_lnd%albedo_vis_dir(l,k)*tile%frac
            lm4_model%atm_sfc%albedo_nir_dir(l) = lm4_model%atm_sfc%albedo_nir_dir(l) + lm4_model%From_lnd%albedo_nir_dir(l,k)*tile%frac
            lm4_model%atm_sfc%albedo_vis_dif(l) = lm4_model%atm_sfc%albedo_vis_dif(l) + lm4_model%From_lnd%albedo_vis_dif(l,k)*tile%frac
            lm4_model%atm_sfc%albedo_nir_dif(l) = lm4_model%atm_sfc%albedo_nir_dif(l) + lm4_model%From_lnd%albedo_nir_dif(l,k)*tile%frac
         end do
      end if

   end subroutine sfc_boundary_layer

   !! ============================================================================
   !! Surface-layer quantities exported to the UFS atmosphere, aggregated over land tiles:
   !! 2 m temperature/humidity (as in GFDL flux_exchange t_ref/q_ref diagnostics),
   !! exchange coefficients cmm/chh and zvfun (as in FV3 sfc_diff / NoahMP)
   !! ============================================================================
   subroutine surface_layer_exports(lm4_model, ex_del_h, ex_del_q, ex_rough_mom)

      type(lm4_type),        intent(inout) :: lm4_model ! land model's variable type
      real, dimension(:),    intent(in)    :: ex_del_h  ! mo_profile weight for heat at z_ref_heat
      real, dimension(:),    intent(in)    :: ex_del_q  ! mo_profile weight for moisture at z_ref_heat
      real, dimension(:),    intent(in)    :: ex_rough_mom ! land roughness length for momentum, m

      real, parameter :: z0_low  = 0.1  ! roughness bounds of zvfun, as in FV3 sfc_diff (m)
      real, parameter :: z0_high = 1.0

      real, dimension(ntiles_lnd) :: ex_t_ref, ex_q_ref, ex_cmm, ex_cd_t_wind, ex_zvfun
      real, dimension(size(lm4_model%atm_sfc%chh)) :: rho
      real :: roughness_function

      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))
         kk = kk + 1
         ex_t_ref(kk) = ex_t_ca(kk) + (lm4_model%atm_forc%t_bot(l) - ex_t_ca(kk)) * ex_del_h(kk)
         ex_q_ref(kk) = ex_tr_surf(kk,isphum) + (lm4_model%atm_forc%q_bot(l) - ex_tr_surf(kk,isphum)) * ex_del_q(kk)
         ex_cmm(kk)       = ex_cd_m(kk) * ex_wind(kk)
         ex_cd_t_wind(kk) = ex_cd_t(kk) * ex_wind(kk)
         roughness_function = min(max((ex_rough_mom(kk) - z0_low) / (z0_high - z0_low), 0.0), 1.0)
         ex_zvfun(kk) = sqrt(roughness_function * max(lm4_model%atm_forc%vfrac(l), 0.1))
      end do

      call aggregate_land_tiles(ex_t_ref,     lm4_model%atm_sfc%t_ref)
      call aggregate_land_tiles(ex_q_ref,     lm4_model%atm_sfc%q_ref)
      call aggregate_land_tiles(ex_cmm,       lm4_model%atm_sfc%cmm)
      call aggregate_land_tiles(ex_cd_t_wind, lm4_model%atm_sfc%chh)
      call aggregate_land_tiles(ex_zvfun,     lm4_model%atm_sfc%zvfun)

      ! chh is a mass flux coefficient: multiply by air density at the lowest atmospheric level
      rho = virtual_temp(lm4_model%atm_forc%t_bot, lm4_model%atm_forc%q_bot)
      rho = air_density(lm4_model%atm_forc%p_bot, rho)
      lm4_model%atm_sfc%chh = rho * lm4_model%atm_sfc%chh

   end subroutine surface_layer_exports

   !! ============================================================================
   !! Adapted from GFDL atm_land_ice_flux_exchange,
   !! stripped down to be "land only" on unstructured grid
   !! In original code, caculates radiation, damping, and vertical diffusion of
   !! momentum, tracers, and downward heat/moisture
   !! For Data Atmosphere, this is only used for gustiness
   !! ============================================================================
   subroutine update_atmos_model_down(lm4_model)

      type(lm4_type),        intent(inout)  :: lm4_model ! land model's variable type

      ! Original code here calculated net shortwave fluxes from downward shortwave fluxes and surface albedos

      call compute_gust(lm4_model)

      !! Original code calculated Atmos%Surf_diff%dtmass and reset Sfc%dt_tr = 0.0 here
      !! If needed by LM4 with active atmosphere, this would need to be added back in

   end subroutine update_atmos_model_down

   !! ============================================================================
   !! Put prescribed atmosphere gust calculation in its own subroutine,
   !! so it can be called both from update_atmos_model_down and at initialization.
   !! Note that:
   !! 1) If alt_gustiness= T from surface_flux_nml, absolute wind in suface
   !!    flux will not use this gustiness.
   !! 2) If gust_to_use= 'computed', there will not be restart prodicibility with
   !!    UFS Data Atmosphere, since u_star and b_star are not saved in restarts.
   !! ============================================================================
   subroutine compute_gust(lm4_model)

      type(lm4_type),        intent(inout)  :: lm4_model ! land model's variable type

      if (trim(gust_to_use)=='computed') then
         !compute gust based on u_star and b_star
         where (ex_b_star > 0.)
            lm4_model%atm_forc%gust = (ex_u_star*ex_b_star*1000)**(1./3.)
         elsewhere
            lm4_model%atm_forc%gust = 0.0
         endwhere
         lm4_model%atm_forc%gust = max(lm4_model%atm_forc%gust, gust_min)
      else if (trim(gust_to_use)=='prescribed') then
         lm4_model%atm_forc%gust = gustiness
      else
         call ESMF_LogWrite('update_atmos_down: illegal value of gust_to_use', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

   end subroutine compute_gust



   !! ============================================================================
   !! Adapted from GFDL atm_land_ice_flux_exchange,
   !! stripped down to be "land only" on unstructured grid
   !! Returns fluxes and derivatives corrected for the implicit treatment of atmospheric
   !! diffusive fluxes, as well as the increments in the temperature and specific humidity
   !! of the lowest atmospheric layer due to all explicit processes as well as the diffusive
   !! fluxes through the top of this layer.
   !! ============================================================================
   subroutine flux_down_from_atmos(dt, lm4_model)

      use constants_mod,    only: cp_air

      real,                  intent(in)     :: dt        ! Time step
      type(lm4_type),        intent(inout)  :: lm4_model ! land model's variable type

      ! local variables

      real, dimension(ntiles_lnd) :: &
         ex_flux_sw, ex_flux_lwd, &
         ex_flux_sw_down_vis_dir, ex_flux_sw_down_total_dir,  &
         ex_flux_sw_down_vis_dif, ex_flux_sw_down_total_dif,  &
         ex_dtmass, ex_gamma, &
         ex_delta_t, ex_dflux_t, ex_e_q_n

      real, dimension(ntiles_lnd,ntcana) ::  &
         ex_dflux_tr, & ! tracer flux change, zero for explicit coupling, needed for implicit coupling
         ex_delta_tr ! tracer tendencies, zero for explicit coupling, needed for implicit coupling

      real :: cp_inv

      integer :: tr

      ! NOTE. Not including here (not needed): 
      ! 1. scale_precip_2d functionality
      ! 2. partition_fprec_from_lprec functionality
      ! 3. sw1way_bug, use_AM3_physics, _USE_LEGACY_LAND_, or SCM functionality
      ! 4. OMP parallelization
      ! 5. Stock changes
      ! 6. data overrides

      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))  
         kk = kk + 1

         ex_flux_sw_down_vis_dir(kk)   = lm4_model%atm_forc%flux_sw_down_vis_dir(l)
         ex_flux_sw_down_vis_dif(kk)   = lm4_model%atm_forc%flux_sw_down_vis_dif(l)
         ex_flux_sw_down_total_dir(kk) = lm4_model%atm_forc%flux_sw_down_vis_dir(l) + lm4_model%atm_forc%flux_sw_down_nir_dir(l)
         ex_flux_sw_down_total_dif(kk) = lm4_model%atm_forc%flux_sw_down_vis_dif(l) + lm4_model%atm_forc%flux_sw_down_nir_dif(l)

         ex_flux_lwd(kk) = lm4_model%atm_forc%flux_lw(l)
      end do

      if (lm4_model%nml%implicit_atm) then
         call ESMF_LogWrite('flux_down_from_atmos: implicit atmosphere coupling not functional', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)

         !ex_dtmass = real(dt)*grav/lm4_model%atm_forc%p_surf

      else  ! ie, explicit atm coupling
         ex_dtmass = 0
         cp_inv = 1.0/cp_air
         ex_delta_tr = 0.0
         ex_dflux_tr = 0.0
         ex_delta_t = 0.0  !
         ex_dflux_t = 0.0  !
      end if         


      do k = 1, ntiles_lnd
         !----- compute net longwave flux (down-up) -----
         ! (note: lw up already in ex_flux_lw)
         ! net LW is not passed to LM4, which uses lwdn_flux, but is kept as in FMS
         ex_flux_lw(k) = ex_flux_lwd(k) - ex_flux_lw(k)

         ! temperature
         ex_gamma(k)      =  1./ (1.0 - ex_dtmass(k)*(ex_dflux_t(k) + ex_dhdt_atm(k)*cp_inv))
         ex_e_t_n(k)      =  ex_dtmass(k)*ex_dhdt_surf(k)*cp_inv*ex_gamma(k)
         ex_f_t_delt_n(k) = (ex_delta_t(k) + ex_dtmass(k) * ex_flux_t(k)*cp_inv) * ex_gamma(k)

         ex_flux_t (k)    =  ex_flux_t(k)        + ex_dhdt_atm(k) * ex_f_t_delt_n(k)
         ex_dhdt_surf(k)  =  ex_dhdt_surf(k)     + ex_dhdt_atm(k) * ex_e_t_n(k)

         ! moisture

         ! moisture vs. surface temperature, assuming saturation
         ex_gamma(k)   =  1.0 / (1.0 - ex_dtmass(k)*(ex_dflux_tr(k,isphum) + ex_dfdtr_atm(k,isphum)))
         ex_e_q_n(k)      =  ex_dtmass(k) * ex_dedt_surf(k) * ex_gamma(k)
         ex_dedt_surf(k)  =  ex_dedt_surf(k) + ex_dfdtr_atm(k,isphum) * ex_e_q_n(k)
         ! original code uses n_exch_tr instead of ntcana
         do tr = 1,ntcana
            ex_gamma(k)   =  1.0 / (1.0 - ex_dtmass(k)*(ex_dflux_tr(k,tr) + ex_dfdtr_atm(k,tr)))

            ex_e_tr_n(k,tr)      =  ex_dtmass(k)*ex_dfdtr_surf(k,tr)*ex_gamma(k)
            ex_f_tr_delt_n(k,tr) = (ex_delta_tr(k,tr)+ex_dtmass(k)*ex_flux_tr(k,tr))*ex_gamma(k)

            ex_flux_tr(k,tr)     =  ex_flux_tr(k,tr) + ex_dfdtr_atm(k,tr)*ex_f_tr_delt_n(k,tr)
            ex_dfdtr_surf(k,tr)  =  ex_dfdtr_surf(k,tr) + ex_dfdtr_atm(k,tr)*ex_e_tr_n(k,tr)
         enddo
      enddo

      ! send to land boundary
      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))  
         kk = kk + 1
         lm4_model%From_atm%t_flux(l,k)                  = ex_flux_t(kk)
         lm4_model%From_atm%sw_flux(l,k)                 = ex_flux_sw(kk)
         lm4_model%From_atm%sw_flux_down_vis_dir(l,k)    = ex_flux_sw_down_vis_dir(kk)
         lm4_model%From_atm%sw_flux_down_total_dir(l,k)  = ex_flux_sw_down_total_dir(kk)
         lm4_model%From_atm%sw_flux_down_vis_dif(l,k)    = ex_flux_sw_down_vis_dif(kk)
         lm4_model%From_atm%sw_flux_down_total_dif(l,k)  = ex_flux_sw_down_total_dif(kk)  

         ! From_atm%lw_flux (net LW) is not used by LM4, see lwdn_flux below
         lm4_model%From_atm%dhdt(l,k)                   = ex_dhdt_surf(kk)
         lm4_model%From_atm%drdt(l,k)                   = ex_drdt_surf(kk)

         lm4_model%From_atm%p_surf(l,k) = ex_p_surf(kk)

         ! for a given grid index, distribute precipitation to its tiles
         lm4_model%From_atm%lprec(l,k)  = lm4_model%atm_forc%lprec(l)
         lm4_model%From_atm%fprec(l,k)  = lm4_model%atm_forc%fprec(l)
         ! set Land's precipitation temperature to atmosphere's temperature
         lm4_model%From_atm%tprec(l,k) = lm4_model%atm_forc%t_bot(l)


         ! These originally had data overrides
         lm4_model%From_atm%lwdn_flux(l,k) = ex_flux_lwd(kk)
         lm4_model%From_atm%ustar(l,k)     = ex_u_star(kk)
         lm4_model%From_atm%drag_q(l,k)    = ex_drag_q(kk)

         !! In original code, also had these, but don't appear needed by LM4:
         !! %cd_m, %cd_t, %wind, %z_bot

      end do
      
      





      lm4_model%From_atm%tr_flux = 0.0
      lm4_model%From_atm%dfdtr = 0.0

      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))  
         kk = kk + 1
         do tr = 1,ntcana
            lm4_model%From_atm%tr_flux(l,k,tr) = ex_flux_tr(kk,tr)
            lm4_model%From_atm%dfdtr(l,k,tr)   = ex_dfdtr_surf(kk,tr)
         enddo
      end do

   end subroutine flux_down_from_atmos

   !! ============================================================================
   !! Adapted from GFDL atm_land_ice_flux_exchange.
   !! Corrects the fluxes for consistency with the new surface temperatures in land
   !! Final increments for temperature and specific humidity in the
   !! lowest atmospheric layer are computed to be returned to the atmospheric model
   !! ============================================================================
   subroutine  flux_up_to_atmos( lm4_model )

      use constants_mod, only : cp_air


      type(lm4_type),        intent(inout)  :: lm4_model ! land model's variable type

      real, dimension(ntiles_lnd) :: &
         ex_t_surf_new, &
         ex_dt_t_surf,  &
         ex_delta_t_n,  &         
         ex_t_ca_new,   &    
         ex_dt_t_ca

      real, dimension(size(lm4_model%atm_sfc%shflx)) :: &
         rho               ! air density at the lowest atmospheric level, per grid cell

      real, dimension(ntiles_lnd,ntcana) ::  &
         ex_tr_surf_new,    & ! updated tracer values at the surface
         ex_dt_tr_surf,     & ! tendency of tracers at the surface
         ex_delta_tr_n

      integer :: tr 

      !----- compute surface temperature change ----- 

      ex_t_surf_new = 200.0 
      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))  
         kk = kk + 1
         ex_t_surf_new(kk) = lm4_model%From_lnd%t_surf(l,k)
         ex_t_ca_new(kk)   = lm4_model%From_lnd%t_ca(l,k)
      end do

      do k = 1,ntiles_lnd
         if(ex_avail(k)) then
            ex_dt_t_ca(k)  = ex_t_ca_new(k)   - ex_t_ca(k)   ! changes in near-surface T
            ex_dt_t_surf(k) = ex_t_surf_new(k) - ex_t_surf(k) ! changes in radiative T
         endif
      enddo

      if (do_forecast) then
         do k = 1,ntiles_lnd
               if(ex_avail(k) .and. (.not.ex_land(k))) then
                  ex_dt_t_ca  (k) = 0.
                  ex_dt_t_surf(k) = 0.
               endif
         enddo
      end if

      !-----------------------------------------------------------------------
      !-----  adjust fluxes and atmospheric increments for
      !-----  implicit dependence on surface temperature -----
      ! original FMS code uses n_exch_tr instead of ntcana
      do tr = 1,ntcana
         ! set up updated surface tracer field so that flux to atmos for absent
         ! tracers is zero
         do k = 1, ntiles_lnd
            if(.not.ex_avail(k)) cycle
            if (ex_dfdtr_surf(k,tr)/=0.0) then
               ex_dt_tr_surf(k,tr) = -ex_flux_tr(k,tr)/ex_dfdtr_surf(k,tr)
            else
               ex_dt_tr_surf(k,tr) = 0.0
            endif
               ex_tr_surf_new(k,tr) = ex_tr_surf(k,tr)+ex_dt_tr_surf(k,tr)
            enddo
      enddo
      
      ! get updated canopy-air tracers from land, per tile (same tile order as ex_ arrays)
      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))
         kk = kk + 1
         ex_tr_surf_new(kk,:) = lm4_model%From_lnd%tr(l,k,:)
      end do

      ! update tracer tendencies in the atmosphere
      do tr = 1,ntcana
         do k = 1,ntiles_lnd
            if(ex_avail(k)) then
               ex_dt_tr_surf(k,tr) = ex_tr_surf_new(k,tr) - ex_tr_surf(k,tr)
               ex_delta_tr_n(k,tr) = ex_f_tr_delt_n(k,tr) + ex_dt_tr_surf(k,tr) * ex_e_tr_n(k,tr)
               ex_flux_tr(k,tr)    = ex_flux_tr(k,tr)     + ex_dt_tr_surf(k,tr) * ex_dfdtr_surf(k,tr)
            endif
         enddo
      enddo

      do k = 1,ntiles_lnd
         ex_delta_t_n(k) = 0.0
         if(ex_avail(k)) then
            ex_flux_t(k)    = ex_flux_t(k)  + ex_dt_t_ca(k)   * ex_dhdt_surf(k)
            ex_flux_lw(k)   = ex_flux_lw(k) - ex_dt_t_surf(k) * ex_drdt_surf(k)
            ex_delta_t_n(k) = ex_f_t_delt_n(k)  + ex_dt_t_ca(k)*ex_e_t_n(k)
         endif
      enddo

      !-----------------------------------------------------------------------
      !---- get mean quantites on atmospheric grid ----
      ! UFS FV3 atm takes fluxes and updated surface states, not the FMS dt_t/dt_tr tendencies

      call aggregate_land_tiles(ex_flux_t, lm4_model%atm_sfc%shflx)
      call aggregate_land_tiles(ex_flux_tr(:,isphum), lm4_model%atm_sfc%lhflx)
      ! not originally to flux_up_to_atmos, but needed by ufs atm
      call aggregate_land_tiles(ex_tr_surf_new(:,isphum), lm4_model%atm_sfc%q_surf)

      if (lm4_model%nml%kinematic_flux) then
         ! convert units from W/m2, using air density at the lowest atmospheric level,
         ! which is what the atmosphere uses to convert these fluxes back
         rho = virtual_temp(lm4_model%atm_forc%t_bot, lm4_model%atm_forc%q_bot)
         rho = air_density(lm4_model%atm_forc%p_bot, rho)
         lm4_model%atm_sfc%shflx = lm4_model%atm_sfc%shflx/(rho*cp_air)  ! SH/(rho*c_p)
         ! lhflx holds the water vapor flux (kg/m2/s), not W/m2, so only divide by rho
         lm4_model%atm_sfc%lhflx = lm4_model%atm_sfc%lhflx/rho           ! E/rho, kg/kg m/s
      endif


   end subroutine flux_up_to_atmos

   !! Wrap up
   !! ===========================================================================
   subroutine end_driver()

      deallocate( &
         ex_flux_t, ex_flux_lw, ex_dhdt_surf, ex_dedt_surf, &
         ex_drdt_surf, ex_dhdt_atm, ex_tr_atm, ex_tr_surf, &
         ex_flux_tr, ex_dfdtr_surf, ex_dfdtr_atm, ex_e_tr_n, &
         ex_f_tr_delt_n, &
         ex_drag_q,     &
         ex_cd_t,       &
         ex_cd_m,       &
         ex_b_star,     &
         ex_u_star,     &
         ex_wind,       &
         ex_z_atm,      &
         ex_p_surf,     &
         ex_avail,      & 
         ex_land,       &
         ex_t_surf,     &
         ex_t_ca,       &
         ex_f_t_delt_n, &   
         ex_e_t_n       &             
              )


   end subroutine end_driver

   ! ============================================================================
   ! tile existence detector: returns a logical value indicating wether component
   ! model tile exists or not
   logical function land_tile_exists(tile)
      use land_tile_mod,      only : land_tile_enum_type, land_tile_type
      type(land_tile_type), pointer :: tile
      land_tile_exists = associated(tile)
   end function land_tile_exists

   ! ============================================================================
   !! aggregate all land tiles in a grid cell
   !! weighted average by land tile fraction

   subroutine aggregate_land_tiles(tile_data, grid_data)
      real, dimension(:), intent(in)  :: tile_data ! data on land tiles
      real, dimension(:), intent(out) :: grid_data ! aggregated data on (unstructured) grid cell

      grid_data = 0.0
      
      ! check dimensions
      if (size(tile_data) /= ntiles_lnd) then
         write(logmsg,'(A,2I0)') 'Error in aggregate_land_tiles: size(tile_data) /= ntiles_lnd', size(tile_data), ntiles_lnd
         call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif
      if (size(grid_data) /= size(lnd%ug_landfrac)) then
         write(logmsg,'(A,2I0)') 'Error in aggregate_land_tiles: size(grid_data) /= size(lnd%ug_landfrac)', size(grid_data), size(lnd%ug_landfrac)
         call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
         call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif

      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))       
         kk = kk + 1
         if (.not. land_tile_exists(tile)) then
            write(logmsg,'(A,I4,A,I2)') 'Not land tile for (l,k) = (', l, ',', k, ')'
            call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)
            cycle
         endif
         grid_data(l) = grid_data(l) + tile_data(kk) * tile%frac

      end do   

   end subroutine aggregate_land_tiles


   subroutine test_aggregate_land_tiles()

      use land_tile_mod,      only : nitems

      real, dimension(ntiles_lnd) :: tile_data
      real, dimension(:), allocatable :: grid_data

      real    :: scaled_frac_data     = 0.0 ! tile data scaled by tile fraction
      real    :: sum_scaled_frac_data = 0.0 ! sum of tile data scaled by tile fraction. SHould be equal to grid data

      real    :: tol                  = 1.0e-6
      integer :: num_tiles            = 0  ! number of land tiles in current grid cell


      allocate(grid_data(size(lnd%ug_landfrac)))


      ! create up some tile data
      call random_number(tile_data)

      grid_data = 0.0

      call aggregate_land_tiles(tile_data, grid_data)

      write(logmsg,'(A)') 'Results of test_aggregate_land_tiles:'
      call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)

      ce = first_elmt(land_tile_map, ls=lnd%ls)
      kk = 0
      do while (loop_over_tiles(ce,tile,i=i,j=j,l=l,k=k))       
         kk = kk + 1

         if (k == 1) then !(first tile of new grid cell)

            num_tiles = nitems(land_tile_map(l))

            write(logmsg,'(A)') '----------------------------------'
            call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)
            write(logmsg,'(A,I4,A,I2,A,F8.4)') 'Grid cell ', l, ' (', num_tiles, ' tiles) = ', grid_data(l)
            call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)

            sum_scaled_frac_data = 0.0

         endif
         
         scaled_frac_data = tile_data(kk) * tile%frac
         write(logmsg,'(A,I4,A,I2,A,I4,A,F8.4,A,F8.4,A,F8.4)') '  Tile(l,k,kk=', l, ',', k, ',', kk, ') ', tile_data(kk), ' *', tile%frac, ' =', scaled_frac_data
         call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)

         sum_scaled_frac_data = sum_scaled_frac_data + scaled_frac_data

         if (k == num_tiles) then  ! on last tile of current grid cell
            write(logmsg,'(A,F8.4,A,F8.4)') '  Sum of scaled tile values:', sum_scaled_frac_data
            call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_INFO)
            if (abs(sum_scaled_frac_data - grid_data(l)) > tol) then
               write(logmsg,'(A,I4,A,F8.4,A,F8.4)') '    MISMATCH FOR GRID CELL ', l, ': ', sum_scaled_frac_data, ' vs ', grid_data(l)
               call ESMF_LogWrite(trim(logmsg), ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
            endif
            sum_scaled_frac_data = 0.0
         endif



      end do

      deallocate(grid_data)

   end subroutine test_aggregate_land_tiles  

end module lm4_driver