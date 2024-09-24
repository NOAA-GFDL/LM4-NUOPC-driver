!! Routines to prepare and run the land model
!! ============================================================================
module lm4_driver

   use ESMF,                 only: ESMF_MethodRemove, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, &
                                   ESMF_FAILURE, ESMF_END_ABORT, ESMF_Finalize, ESMF_LOGMSG_ERROR, ESMF_LOGMSG_WARNING
   use mpp_domains_mod,      only: domain2d
   use mpp_mod,              only: mpp_pe, mpp_root_pe

   use lm4_type_mod,         only: lm4_type
   use lm4_kind_mod,         only: r8 => shr_kind_r8, cl=>shr_kind_cl
   use land_data_mod,        only: land_data_type, atmos_land_boundary_type, lnd
   use land_tracers_mod,     only: isphum, ico2, ntcana
   use lm4_surface_flux_mod, only: lm4_surface_flux_1d

   use time_manager_mod,      only: time_type, date_to_string, increment_date, decrement_date
   use time_manager_mod,      only: operator(>=), operator(<), operator(==)
   use mpp_domains_mod,       only: domainUG
   use mpp_domains_mod,       only: mpp_get_UG_compute_domain, mpp_get_UG_domain_ntiles
   use mpp_domains_mod,       only: mpp_get_UG_domain_grid_index
   use diag_manager_mod,      only: diag_axis_init, register_static_field, diag_field_add_attribute, &
                                    register_diag_field
   use diag_axis_mod,         only: diag_axis_add_attribute
   use land_tile_diag_mod,    only: register_tiled_area_fields, register_tiled_diag_field, set_default_diag_filter, &
                                    send_tile_data
   use tile_diag_buff_mod,    only: diag_buff_type, init_diag_buff


   implicit none
   private

   type(domain2D), public :: land_domain

   integer :: date_init(6)
   character(len=32)     :: timestamp


   public :: lm4_nml_read
   public :: init_driver, end_driver
   public :: sfc_boundary_layer, update_atmos_model_down, flux_down_from_atmos
   public :: write_int_restart
   public :: debug_diag


   ! --- namelist of vars originally from flux exchange nml
   real :: z_ref_heat =  2. !< Reference height (meters) for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q)
   ! TODO: rename this nml?
   namelist /flux_exchange_nml/ z_ref_heat

   ! --- namelist of vars originally from atmos_prescr_nml
   character(len=24) :: gust_to_use = 'computed' ! or 'prescribed'
   real    :: gustiness    = 5.0  ! m/s, wind gustiness if gust_to_use = 'prescribed'
   real    :: gust_min     = 0.0  ! m/s, minimum gustiness when gust_to_use = 'computed'
   namelist /atmos_prescr_nml/ gustiness, gust_to_use, gust_min


   logical :: scale_precip_2d = .false.

   ! variables for between subroutines
   real, allocatable, dimension(:) :: &
      ex_flux_t, ex_flux_lw,      &
      ex_dhdt_surf, ex_dedt_surf, &
      ex_drdt_surf,  ex_dhdt_atm, &
      ex_drag_q,    &   !< q drag.coeff.
      ex_cd_t,      &
      ex_cd_m,      &
      ex_b_star,    &
      ex_u_star,    &
      ex_wind,      &
      ex_z_atm   

   ! these originally had a tracer dimension
   real, allocatable, dimension(:,:) :: &
      ex_tr_atm,     &
      ex_tr_surf,    & !< near-surface tracer fields
      ex_flux_tr,    & !< tracer fluxes
      ex_dfdtr_surf, & !< d(tracer flux)/d(surf tracer)
      ex_dfdtr_atm,  & !< d(tracer flux)/d(atm tracer)
      ex_e_tr_n,     & !< coefficient in implicit scheme
      ex_f_tr_delt_n   !< coefficient in implicit scheme


   !integer :: n_exch_tr !< number of tracers exchanged between models
   integer :: ntile = 1 ! true for now, with no subtiling

   ! integers for diag manager fields (TODO: clean up)
   integer :: id_cellarea
   integer :: id_lon, id_lat, id_lonb, id_latb

   ! fields to be written out
   integer :: &
      id_swdn_vf, id_z_bot, id_t_bot, id_p_bot, &
      id_u_bot, id_v_bot, id_q_bot, id_p_surf,  &
      id_lprec, id_fprec,            &
      id_flux_lw, id_flux_sw_dn_vdf, id_flux_sw_dn_vr

   !! for unstructured grid diagnostics
   integer :: id_band !, id_zfull ! IDs of land diagnostic axes
   integer :: id_ug !<Unstructured axis id.
   ! unstructured grid diag field ids
   integer :: &
      id_landfrac,    &
      id_geolon_t, id_geolat_t,    &
      id_frac, id_area, id_ntiles, &
      iug_q_atm, &
      iug_t_atm, &
      iug_u_atm, &
      iug_v_atm, &
      iug_p_atm, &
      iug_z_atm, &
      iug_p_surf, &
      iug_t_surf, &
      iug_t_ca, &
      iug_q_surf, &
      iug_rough_mom  , &
      iug_rough_heat , &
      iug_rough_moist  , &
      iug_rough_scale  , &
      iug_gust

   character(len=CL) :: logmsg

contains

   !! Read in lm4 namelist for NUOPC-cap related variables
   !! Also read in other namelists that are not read in from other model components,
   !! but are used in this module
   !! ============================================================================
   subroutine lm4_nml_read(lm4_model)

      use fms_mod,             only: check_nml_error, close_file, file_exist
#ifdef INTERNAL_FILE_NML
      use mpp_mod,             only: input_nml_file
#else
      use fms_mod,             only: open_namelist_file
#endif

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

      ! TODO: are all these still needed?

      ! for namelist read
      integer :: unit, io, ierr
      namelist /lm4_nml/ grid, npx, npy, layout, ntiles, &
         blocksize, lm4_debug, dt_lnd_slow, restart_interval

      ! read in namelists
      ! ------------------------------------------
      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
         ! lm4_nml
         read(input_nml_file, nml=lm4_nml, iostat=io)
         ierr = check_nml_error(io, 'lm4_nml')

         read(input_nml_file, nml=atmos_prescr_nml, iostat=io)
         ierr = check_nml_error(io, 'atmos_prescr_nml')

         read(input_nml_file, nml=flux_exchange_nml, iostat=io)
         ierr = check_nml_error(io, 'flux_exchange_nml')
#else
         unit = open_namelist_file ( )
         ierr=1
         do while (ierr /= 0)
            read(unit, nml=lm4_nml, iostat=io)
            ierr = check_nml_error(io,'lm4_nml')
         enddo

         ierr=1
         do while (ierr /= 0)
            read(unit, nml=atmos_prescr_nml, iostat=io)
            ierr = check_nml_error(io,'atmos_prescr_nml')
         enddo        
         
         ierr=1
         do while (ierr /= 0)
            read(unit, nml=flux_exchange_nml, iostat=io)
            ierr = check_nml_error(io,'flux_exchange_nml')
         enddo             
         call close_file(unit)
#endif
      endif

      lm4_model%nml%lm4_debug   = lm4_debug
      lm4_model%nml%grid        = grid
      lm4_model%nml%blocksize   = blocksize
      lm4_model%nml%npx         = npx
      lm4_model%nml%npy         = npy
      lm4_model%nml%layout      = layout
      lm4_model%nml%ntiles      = ntiles
      lm4_model%nml%dt_lnd_slow = dt_lnd_slow
      lm4_model%nml%restart_interval = restart_interval

   end subroutine lm4_nml_read

   !! ============================================================================
   subroutine init_driver(lm4_model)
      !! TODO: cleanup unused code

      use mpp_domains_mod,    only: domain2d, mpp_get_compute_domain
      use mpp_mod,            only: mpp_pe, mpp_root_pe
      use land_domain_mod,    only: domain_create
      use block_control_mod,  only: block_control_type, define_blocks_packed
      !use land_restart_mod,   only: sfc_prop_restart_read, sfc_prop_transfer

      type(lm4_type),          intent(inout) :: lm4_model ! land model's variable type

      ! ---------------
      ! local
      !type (block_control_type), target   :: Lnd_block !  Block container
      !integer                             ::   blocksize
      !logical, save                       :: block_message = .true.


      integer :: isc, iec, jsc, jec

      ! FMS domain creation:
      call domain_create(lm4_model%nml, land_domain)

      call mpp_get_compute_domain(land_domain,isc,iec,jsc,jec)

      allocate( &
         ex_flux_t(lnd%ls:lnd%le), ex_flux_lw(lnd%ls:lnd%le),   &
         ex_dhdt_surf(lnd%ls:lnd%le), ex_dedt_surf(lnd%ls:lnd%le), &
         ex_drdt_surf(lnd%ls:lnd%le),  ex_dhdt_atm(lnd%ls:lnd%le), &
         ex_drag_q(lnd%ls:lnd%le),    &   !< q drag.coeff.
         ex_cd_t(lnd%ls:lnd%le),      &
         ex_cd_m(lnd%ls:lnd%le),      &
         ex_b_star(lnd%ls:lnd%le),    &
         ex_u_star(lnd%ls:lnd%le),    &
         ex_wind(lnd%ls:lnd%le),      &
         ex_z_atm(lnd%ls:lnd%le)      &            
      )

      ! these originally had a tracer dimension
      allocate( &
         ex_tr_atm(lnd%ls:lnd%le,ntcana),  &
         ex_tr_surf(lnd%ls:lnd%le,ntcana),    & !< near-surface tracer fields
         ex_flux_tr(lnd%ls:lnd%le,ntcana),    & !< tracer fluxes
         ex_dfdtr_surf(lnd%ls:lnd%le,ntcana), & !< d(tracer flux)/d(surf tracer)
         ex_dfdtr_atm(lnd%ls:lnd%le,ntcana),  & !< d(tracer flux)/d(atm tracer)
         ex_e_tr_n(lnd%ls:lnd%le,ntcana),     & !< coefficient in implicit scheme
         ex_f_tr_delt_n(lnd%ls:lnd%le,ntcana) &  !< coefficient in implicit scheme      

      )

      ! Set restart time  
      if (ALL(lm4_model%nml%restart_interval ==0)) then
         lm4_model%Time_restart = increment_date(lm4_model%Time_end, 0, 0, 10, 0, 0, 0)   ! no intermediate restart
      else
         
         lm4_model%Time_restart = increment_date(lm4_model%Time_init, lm4_model%nml%restart_interval(1), lm4_model%nml%restart_interval(2), &
            lm4_model%nml%restart_interval(3), lm4_model%nml%restart_interval(4), lm4_model%nml%restart_interval(5), lm4_model%nml%restart_interval(6) )

            ! subtract the slow time step in seconds
            timestamp = date_to_string(lm4_model%Time_restart)
            call ESMF_LogWrite('LM4 init_driver: Time_restart before decrement' //trim(timestamp), ESMF_LOGMSG_INFO)
            lm4_model%Time_restart = decrement_date(lm4_model%Time_restart, 0,0,0,0,0, lm4_model%nml%dt_lnd_slow)
            timestamp = date_to_string(lm4_model%Time_restart)
            call ESMF_LogWrite('LM4 init_driver: Time_restart after decrement' //trim(timestamp), ESMF_LOGMSG_INFO)
         

         if (lm4_model%Time_restart < lm4_model%Time_land) then
            call ESMF_LogWrite('The first intermediate restart time is larger than the start time', &
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
         call compute_gust(lm4_model)
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

         ! To match behavior of GFDL coupler, advance the restart file timestamp by the slow time step
         ! (datestamp is modeltime needed for restart)
         Time_restart_stamp = increment_date(lm4_model%Time_land, 0,0,0,0,0, lm4_model%nml%dt_lnd_slow)
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
   ! subroutine sfc_boundary_layer( dt )


      use sat_vapor_pres_mod, only: compute_qs
      use monin_obukhov_mod,  only: mo_profile
      use land_tile_mod,      only: max_n_tiles   ! this should be 1 for now

      real,                  intent(in)     :: dt        ! Time step
      type(lm4_type),        intent(inout)  :: lm4_model ! land model's variable type

      ! local variables --------------------------------------
      !  !! blocking not used for now
      !  integer :: nblocks = 1
      !  integer :: my_nblocks = 1
      !  integer, allocatable :: block_start(:), block_end(:)

      real    :: zrefm, zrefh

      real, dimension(lnd%ls:lnd%le) :: &
         ex_albedo,             &
         ex_albedo_vis_dir,     &
         ex_albedo_nir_dir,     &
         ex_albedo_vis_dif,     &
         ex_albedo_nir_dif,     &
      !ex_land_frac,          &
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
         ex_ref, ex_ref_u, ex_ref_v, ex_u10, &
         ex_ref2,               &
         ex_t_ref,              &
         ex_qs_ref,             &
         ex_qs_ref_cmip,        &
         ex_del_m,              &
         ex_del_h,              &
         ex_del_q,              &
         ex_frac_open_sea


      !! -- These were originally allocatable:
      !!

      logical, dimension(lnd%ls:lnd%le) :: &
         ex_avail,     &   !< true where data on exchange grid are available
         ex_land,      &   !< true if exchange grid cell is over land
         ex_seawater       !< true if exchange grid cell is over seawater

      real, dimension(lnd%ls:lnd%le) :: &
         ex_t_surf   ,  &
         ex_t_ca     ,  &
         ex_t_surf_miz, &
         ex_p_surf   ,  &
         ex_q_surf  ,  &
      !ex_slp      ,  &
         ex_dqsatdt_surf,  &
         ex_f_t_delt_n, &

      ! MOD these were moved from local ! so they can be passed to flux down
         ex_flux_u,    &
         ex_flux_v,    &
         ex_dtaudu_atm,&
         ex_dtaudv_atm,&

      ! values added for LM3

         ex_e_t_n    ,  &
         !ex_e_q_n    ,  &

      !
         ex_albedo_fix,        &
         ex_albedo_vis_dir_fix,&
         ex_albedo_nir_dir_fix,&
         ex_albedo_vis_dif_fix,&
         ex_albedo_nir_dif_fix

      integer :: tr, n, m ! tracer indices
      integer :: i



      ! ---------------------------

      ex_avail    = .TRUE.
      ex_land     = .TRUE.
      ex_seawater = .FALSE.

      ! these are 0 for land
      ex_u_surf = 0.0
      ex_v_surf = 0.0

      ! prefill q_surf with q_bot. In original code, there were options to send/write
      ! before surface_flux call.
      ex_tr_surf(:,isphum) = lm4_model%atm_forc%q_bot
      ex_q_surf            = lm4_model%atm_forc%q_bot

      ! TODO: review
      do tr = 1,ntcana
         ex_tr_surf(:,tr) = lm4_model%From_lnd%tr(:,ntile,tr)
      enddo

      ex_t_atm  = lm4_model%atm_forc%t_bot
      ex_q_atm  = lm4_model%atm_forc%q_bot
      ex_u_atm  = lm4_model%atm_forc%u_bot
      ex_v_atm  = lm4_model%atm_forc%v_bot
      ex_p_atm  = lm4_model%atm_forc%p_bot
      ex_z_atm  = lm4_model%atm_forc%z_bot
      ex_p_surf = lm4_model%atm_forc%p_surf
      ex_gust   = lm4_model%atm_forc%gust


      ! this is mimicking the original code for land roughness vars
      ! TODO: review and optimize
      ex_rough_mom   = lm4_model%From_lnd%rough_mom(:,ntile)
      ex_rough_moist = lm4_model%From_lnd%rough_heat(:,ntile)
      ex_rough_heat  = lm4_model%From_lnd%rough_heat(:,ntile)
      ex_rough_scale = lm4_model%From_lnd%rough_scale(:,ntile)

      ex_t_surf = lm4_model%From_lnd%t_surf(:,ntile)
      ex_t_ca   = lm4_model%From_lnd%t_ca(:,ntile)


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

      ! Note, these are relying on the send_tile_data_0d_array procedure
      ! to avoid bring land_mdel.F90's tile buffer into calls
      ! Note: only use for data with tile dim.?
      ! call send_tile_data(iug_q_atm      , ex_q_atm      )
      ! call send_tile_data(iug_t_atm      , ex_t_atm      )
      ! call send_tile_data(iug_u_atm      , ex_u_atm      )
      ! call send_tile_data(iug_v_atm      , ex_v_atm      )
      ! call send_tile_data(iug_p_atm      , ex_p_atm      )
      ! call send_tile_data(iug_z_atm      , ex_z_atm      )
      ! call send_tile_data(iug_p_surf     , ex_p_surf     )
      ! call send_tile_data(iug_t_surf     , ex_t_surf     )
      ! call send_tile_data(iug_t_ca       , ex_t_ca       )
      ! call send_tile_data(iug_q_surf     , ex_q_surf     )
      ! call send_tile_data(iug_rough_mom  , ex_rough_mom  )
      ! call send_tile_data(iug_rough_heat , ex_rough_heat )
      ! call send_tile_data(iug_rough_moist, ex_rough_moist)
      ! call send_tile_data(iug_rough_scale, ex_rough_scale)
      ! call send_tile_data(iug_gust       , ex_gust       )


         ex_tr_surf(1,isphum) = ex_q_surf(1)  ! TODO: review this connection 


      !1 TODO: make sure output args that are used outside of this routine have the right scope
      call lm4_surface_flux_1d ( &
         ! inputs
         lm4_model%atm_forc%t_bot, lm4_model%atm_forc%q_bot,  & !! TODO: link q_bot var and tracer field
         lm4_model%atm_forc%u_bot, lm4_model%atm_forc%v_bot,  lm4_model%atm_forc%p_bot, &
         lm4_model%atm_forc%z_bot, lm4_model%atm_forc%p_surf, lm4_model%From_lnd%t_surf(:,ntile), &
         lm4_model%From_lnd%t_ca(:,ntile), &
         ! inout
         ex_tr_surf(:,isphum),         & !! TODO review using q_bot as surface Q (this is inout).
         ! more inputs
         ex_u_surf, ex_v_surf,             & ! 0s
         lm4_model%From_lnd%rough_mom(:,ntile), lm4_model%From_lnd%rough_heat(:,ntile), &
         ex_rough_moist, lm4_model%From_lnd%rough_scale(:,ntile),   &
         lm4_model%atm_forc%gust,                                                       & ! gustiness
         ! outputs
         ex_flux_t, ex_flux_tr(:,isphum), ex_flux_lw,   &
         ex_flux_u, ex_flux_v, ex_cd_m,   ex_cd_t, &
         ex_cd_q,   ex_wind,   ex_u_star, ex_b_star, &
         ex_q_star, ex_dhdt_surf, ex_dedt_surf, &
         ex_dfdtr_surf(:,isphum),  ex_drdt_surf,  ex_dhdt_atm, &
         ex_dfdtr_atm(:,isphum),  ex_dtaudu_atm, ex_dtaudv_atm,         &
         dt,                                                            & !! timestep doesn't seem to be used
         ex_land, ex_seawater, ex_avail                                       & ! Is land, Is seawater, Is ex. avail
         )

         ex_q_surf(1) = ex_tr_surf(1,isphum)  ! TODO: review this connection


      !! ....
      zrefm = 10.0
      zrefh = z_ref_heat
      !      ---- optimize calculation ----
      ! write(*,*) 'DEBUG: calling mo_profile'
      call mo_profile ( zrefm, zrefh, lm4_model%atm_forc%z_bot, ex_rough_mom, &
         ex_rough_heat, ex_rough_moist,          &
         ex_u_star, ex_b_star, ex_q_star,        &
         ex_del_m, ex_del_h, ex_del_q, ex_avail  )
      ! write(*,*) 'DEBUG: done calling mo_profile'

      do i = lnd%ls,lnd%le
         ex_u10(i) = 0.
         if(ex_avail(i)) then
            ex_ref_u(i) = ex_u_surf(i) + (lm4_model%atm_forc%u_bot(i)-ex_u_surf(i)) * ex_del_m(i)
            ex_ref_v(i) = ex_v_surf(i) + (lm4_model%atm_forc%v_bot(i)-ex_v_surf(i)) * ex_del_m(i)
            ex_u10(i) = sqrt(ex_ref_u(i)**2 + ex_ref_v(i)**2)
         endif
      enddo

      !! TODO: is this needed?
      ! do n = 1, ex_gas_fields_atm%num_bcs  !{
      !    if (atm%fields%bc(n)%use_10m_wind_speed) then  !{
      !       if (.not. ex_gas_fields_atm%bc(n)%field(ind_u10)%override) then  !{
      !          do i = lnd%ls,lnd%le
      !             ex_gas_fields_atm%bc(n)%field(ind_u10)%values(i) = ex_u10(i)
      !          enddo
      !       endif  !}
      !    endif  !}
      ! enddo  !} n

      do i = lnd%ls,lnd%le
         if(ex_avail(i)) ex_drag_q(i) = ex_wind(i)*ex_cd_q(i)
         ! [6] get mean quantities on atmosphere grid
         ! [6.1] compute t surf for radiation
         ex_t_surf4(i) = ex_t_surf(i) ** 4
      enddo

      ! [6.3] save atmos albedo fix and old albedo (for downward SW flux calculations)
      ! on exchange grid
      do i = lnd%ls,lnd%le
         ex_albedo_fix(i) = 0.
         ex_albedo_vis_dir_fix(i) = 0.
         ex_albedo_nir_dir_fix(i) = 0.
         ex_albedo_vis_dif_fix(i) = 0.
         ex_albedo_nir_dif_fix(i) = 0.
      enddo


      ! TODO: convert these from  xgrid and Land_Ice_Atmos_Boundary to atmos_land_boundary_type?
      ! [6.2] put relevant quantities onto atmospheric boundary
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%t,         'ATM', ex_t_surf4  ,  xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%frac_open_sea,'ATM',ex_frac_open_sea, xmap_sfc)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo,    'ATM', ex_albedo   ,  xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir,    'ATM',   &
      !      ex_albedo_vis_dir   ,  xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir,    'ATM',   &
      !      ex_albedo_nir_dir   ,  xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif,    'ATM',   &
      !      ex_albedo_vis_dif   ,  xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif,    'ATM',   &
      !      ex_albedo_nir_dif   ,  xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%rough_mom, 'ATM', ex_rough_mom,  xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%land_frac, 'ATM', ex_land_frac,  xmap_sfc, complete=.false.)

      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%u_flux,    'ATM', ex_flux_u,     xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%v_flux,    'ATM', ex_flux_v,     xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudu,    'ATM', ex_dtaudu_atm, xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudv,    'ATM', ex_dtaudv_atm, xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%u_star,    'ATM', ex_u_star    , xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%b_star,    'ATM', ex_b_star    , xmap_sfc, complete=.false.)
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%q_star,    'ATM', ex_q_star    , xmap_sfc, complete=.true.)

      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%u_ref,     'ATM', ex_ref_u     , xmap_sfc, complete=.false.) !bqx
      ! call get_from_xgrid (Land_Ice_Atmos_Boundary%v_ref,     'ATM', ex_ref_v     , xmap_sfc, complete=.true.) !bqx



      ! ! TODO: Is this needed? do these have the right scope?
      ! do i = lnd%ls,lnd%le
      !    ex_albedo_fix(i) = (1.0-ex_albedo(i)) / (1.0-ex_albedo_fix(i))
      !    ex_albedo_vis_dir_fix(i) = (1.0-ex_albedo_vis_dir(i)) / (1.0-ex_albedo_vis_dir_fix(i))
      !    ex_albedo_nir_dir_fix(i) = (1.0-ex_albedo_nir_dir(i)) / (1.0-ex_albedo_nir_dir_fix(i))
      !    ex_albedo_vis_dif_fix(i) = (1.0-ex_albedo_vis_dif(i)) / (1.0-ex_albedo_vis_dif_fix(i))
      !    ex_albedo_nir_dif_fix(i) = (1.0-ex_albedo_nir_dif(i)) / (1.0-ex_albedo_nir_dif_fix(i))
      ! enddo


      !=======================================================================
      ! [7] diagnostics section

      ! call mo_profile ( zrefm, zrefh, lm4_model%atm_forc%z_bot,   ex_rough_mom(:), &
      !    ex_rough_heat(:), ex_rough_moist(:),          &
      !    ex_u_star(:), ex_b_star(:), ex_q_star(:),        &
      !    ex_del_m(:), ex_del_h(:), ex_del_q(:), ex_avail(:)  )

      ! !    ------- reference relative humidity -----------
      ! !cjg     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
      ! !cjg          id_rh_ref_cmip > 0 .or. &
      ! !cjg          id_q_ref > 0 .or. id_q_ref_land >0 ) then
      ! do i = lnd%ls,lnd%le
      !    ex_ref(i) = 1.0e-06
      !    if (ex_avail(i)) &
      !       ex_ref(i)   = ex_tr_surf(i,isphum) + (ex_tr_atm(i,isphum)-ex_tr_surf(i,isphum)) * ex_del_q(i)
      ! enddo

      ! do i = lnd%ls,lnd%le
      !    ex_t_ref(i) = 200.
      !    if(ex_avail(i)) &
      !       ex_t_ref(i) = ex_t_ca(i) + (lm4_model%atm_forc%t_bot(i)-ex_t_ca(i)) * ex_del_h(i)
      ! enddo
      ! call compute_qs (ex_t_ref(:), ex_p_surf(:), ex_qs_ref(:), q = ex_ref(:))
      ! call compute_qs (ex_t_ref(:), ex_p_surf(:), ex_qs_ref_cmip(:),  &
      !    q = ex_ref(:), es_over_liq_and_ice = .true.)
      ! do i = lnd%ls,lnd%le
      !    if(ex_avail(i)) then
      !       ! remove cap on relative humidity -- this mod requested by cjg, ljd
      !       !RSH    ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)
      !       ex_ref2(i)   = 100.*ex_ref(i)/ex_qs_ref_cmip(i)
      !       ex_ref(i)    = 100.*ex_ref(i)/ex_qs_ref(i)
      !    endif
      ! enddo

      ! ! lots of send_data stuff originally here, removed
      ! ! TODO: get diag history write back in

    ! JP DEBUG
      ! if (mpp_pe()== mpp_root_pe() ) then
      !    write(*,*) 'calling write_data at end of sfc_boundary_layer-------'
      !    call write_data(lm4_model)
      ! endif

   end subroutine sfc_boundary_layer

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
      ! TODO: move data atmosphere sw flux calculation from imports here?
      
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
   ! subroutine flux_down_from_atmos()

      use constants_mod,    only: cp_air, grav

      real,                  intent(in)     :: dt        ! Time step
      type(lm4_type),        intent(inout)  :: lm4_model ! land model's variable type

      ! local variables

      real, dimension(lnd%ls:lnd%le) :: &
         ex_flux_sw, ex_flux_lwd, &
         ex_flux_sw_dir,  &
         ex_flux_sw_dif,  &
         ex_flux_sw_down_vis_dir, ex_flux_sw_down_total_dir,  &
         ex_flux_sw_down_vis_dif, ex_flux_sw_down_total_dif,  &
         ex_flux_sw_vis, &
         ex_flux_sw_vis_dir, &
         ex_flux_sw_vis_dif, &
         ! ex_tprec, & ! temperature of precipitation, currently equal to atm T
         ex_dtmass, ex_gamma, ex_e_t_n, ex_f_t_delt_n, &
         ex_delta_t, ex_dflux_t, ex_e_q_n

         real, dimension(lnd%ls:lnd%le,ntcana) ::  &
          ex_dflux_tr, & ! tracer flux change. TODO: NEED TO FETCH
          ex_delta_tr ! tracer tendencies. TODO: NEED TO FETCH

      real :: cp_inv

      integer :: l, tr

      ! NOTE. Not including here: (TODO: REVIEW)
      ! 1. scale_precip_2d functionality
      ! 2. partition_fprec_from_lprec functionality
      ! 3. sw1way_bug, use_AM3_physics, _USE_LEGACY_LAND_, or SCM functionality
      ! 4. OMP parallelization
      ! 5. Stock changes
      ! 6. data overrides


      ! ex_flux_sw_dir            =
      ! ex_flux_sw_vis_dir        =
      ! ex_flux_sw_dif            =
      ! ex_flux_sw_vis_dif        =
      ex_flux_sw_down_vis_dir   = lm4_model%atm_forc%flux_sw_down_vis_dir
      ex_flux_sw_down_vis_dif   = lm4_model%atm_forc%flux_sw_down_vis_dif
      ex_flux_sw_down_total_dir = lm4_model%atm_forc%flux_sw_down_vis_dir + lm4_model%atm_forc%flux_sw_down_nir_dir
      ex_flux_sw_down_total_dif = lm4_model%atm_forc%flux_sw_down_vis_dif + lm4_model%atm_forc%flux_sw_down_nir_dif

      ! TODO: ex_flux_lwd
      ex_flux_lwd = lm4_model%atm_forc%flux_lw

      ! this echos standalone LM4.0 driver
      ex_dtmass = real(dt)*grav/lm4_model%atm_forc%p_surf
      cp_inv = 1.0/cp_air

      ! TODO: review behavior of these.
      ! With data atmosphere, no implicit derivatives. Would need to review
      ! for implicit coupling with active atmosphere.
      ex_delta_tr = 0.0
      ex_dflux_tr = 0.0
      ex_delta_t = 0.0  ! 
      ex_dflux_t = 0.0  !

      do l = lnd%ls,lnd%le
         !----- compute net longwave flux (down-up) -----
         ! (note: lw up already in ex_flux_lw)
         ex_flux_lw(l) = ex_flux_lwd(l) - ex_flux_lw(l)  ! TODO: review if needed 

         ! temperature
         ex_gamma(l)      =  1./ (1.0 - ex_dtmass(l)*(ex_dflux_t(l) + ex_dhdt_atm(l)*cp_inv))
         ex_e_t_n(l)      =  ex_dtmass(l)*ex_dhdt_surf(l)*cp_inv*ex_gamma(l)
         ex_f_t_delt_n(l) = (ex_delta_t(l) + ex_dtmass(l) * ex_flux_t(l)*cp_inv) * ex_gamma(l)

         ex_flux_t (l)    =  ex_flux_t(l)        + ex_dhdt_atm(l) * ex_f_t_delt_n(l)
         ex_dhdt_surf(l)  =  ex_dhdt_surf(l)     + ex_dhdt_atm(l) * ex_e_t_n(l)


         ! moisture

         ! moisture vs. surface temperture, assuming saturation
         ex_gamma(l)   =  1.0 / (1.0 - ex_dtmass(l)*(ex_dflux_tr(l,isphum) + ex_dfdtr_atm(l,isphum)))
         ex_e_q_n(l)      =  ex_dtmass(l) * ex_dedt_surf(l) * ex_gamma(l)
         ex_dedt_surf(l)  =  ex_dedt_surf(l) + ex_dfdtr_atm(l,isphum) * ex_e_q_n(l)
         ! original code uses n_exch_tr instead of ntcana
         do tr = 1,ntcana
            ex_gamma(l)   =  1.0 / (1.0 - ex_dtmass(l)*(ex_dflux_tr(l,tr) + ex_dfdtr_atm(l,tr)))

            ex_e_tr_n(l,tr)      =  ex_dtmass(l)*ex_dfdtr_surf(l,tr)*ex_gamma(l)
            ex_f_tr_delt_n(l,tr) = (ex_delta_tr(l,tr)+ex_dtmass(l)*ex_flux_tr(l,tr))*ex_gamma(l)

            ex_flux_tr(l,tr)     =  ex_flux_tr(l,tr) + ex_dfdtr_atm(l,tr)*ex_f_tr_delt_n(l,tr)
            ex_dfdtr_surf(l,tr)  =  ex_dfdtr_surf(l,tr) + ex_dfdtr_atm(l,tr)*ex_e_tr_n(l,tr)
         enddo
      enddo

      ! send to land boundary

      lm4_model%From_atm%t_flux(:,ntile)                  = ex_flux_t
      lm4_model%From_atm%sw_flux(:,ntile)                 = ex_flux_sw
      lm4_model%From_atm%sw_flux_down_vis_dir(:,ntile)    = ex_flux_sw_down_vis_dir
      lm4_model%From_atm%sw_flux_down_total_dir(:,ntile)  = ex_flux_sw_down_total_dir
      lm4_model%From_atm%sw_flux_down_vis_dif(:,ntile)    = ex_flux_sw_down_vis_dif
      lm4_model%From_atm%sw_flux_down_total_dif(:,ntile)  = ex_flux_sw_down_total_dif
      ! TODO: review if is this net LW needed by land?
      ! lm4_model%From_atm%lw_flux(:,ntile)                 = ex_flux_lw
      lm4_model%From_atm%dhdt(:,ntile)                    = ex_dhdt_surf
      lm4_model%From_atm%drdt(:,ntile)                    = ex_drdt_surf
      ! TODO: review if should replace with wider scope versions of ex_p_surf, ex_lprec, ex_fprec
      lm4_model%From_atm%p_surf(:,ntile) = lm4_model%atm_forc%p_surf !ex_p_surf
      lm4_model%From_atm%lprec(:,ntile)  = lm4_model%atm_forc%lprec !ex_lprec
      lm4_model%From_atm%fprec(:,ntile)  = lm4_model%atm_forc%fprec !ex_fprec

      ! set Land's precipitation temperature to atmosphere's temperature
      lm4_model%From_atm%tprec(:,ntile) = lm4_model%atm_forc%t_bot

      ! TODO: review scope of these
      ! These originally had data overrides
      if(associated(lm4_model%From_atm%drag_q)) then
         lm4_model%From_atm%drag_q(:,ntile) = ex_drag_q
      endif
      if(associated(lm4_model%From_atm%lwdn_flux)) then
         lm4_model%From_atm%lwdn_flux(:,ntile) = ex_flux_lwd
      endif
      if(associated(lm4_model%From_atm%cd_m)) then
         lm4_model%From_atm%cd_m(:,ntile) = ex_cd_m
      endif
      if(associated(lm4_model%From_atm%cd_t)) then
         lm4_model%From_atm%cd_t(:,ntile) = ex_cd_t
      endif
      if(associated(lm4_model%From_atm%bstar)) then
         lm4_model%From_atm%bstar(:,ntile) = ex_b_star
      endif
      if(associated(lm4_model%From_atm%ustar)) then
         lm4_model%From_atm%ustar(:,ntile) = ex_u_star
      endif
      if(associated(lm4_model%From_atm%wind)) then
         lm4_model%From_atm%wind(:,ntile) = ex_wind
      endif
      if(associated(lm4_model%From_atm%z_bot)) then
         lm4_model%From_atm%z_bot(:,ntile) = ex_z_atm
      endif


      lm4_model%From_atm%tr_flux = 0.0
      lm4_model%From_atm%dfdtr = 0.0
      do tr = 1,ntcana
            lm4_model%From_atm%tr_flux(:,ntile,tr) = ex_flux_tr(:,tr)
            lm4_model%From_atm%dfdtr(:,ntile,tr)   = ex_dfdtr_surf(:,tr)
      enddo

   end subroutine flux_down_from_atmos

   !! ============================================================================
   !! flux_up_to_atmos will be completed for 2-way coupling to active atmosphere
   !! ============================================================================
   subroutine  flux_up_to_atmos( Land )

      type(land_data_type),  intent(in)    :: Land !< A derived data type to specify land boundary data

      !   where (Land%mask(:,:,1))
      !      t_surf_new = Land%t_surf(:,:,1)
      !      t_ca_new   = Land%t_ca  (:,:,1)
      !   endwhere

      !   !??????? should this be done in land model ??????
      !   call escomp (t_surf_new, q_surf_new)
      !   where (Land%mask(:,:,1))
      !      q_surf_new = Land%tr(:,:,1,1)
      !   elsewhere
      !      !q_surf_new = d622*q_surf_new/(p_surf-d378*q_surf_new)
      !   endwhere

      !   dt_t_ca   = t_ca_new   - t_ca   ! changes in near-surface T
      !   dt_t_surf = t_surf_new - t_surf ! changes in radiative T
      !   dt_q_surf = q_surf_new - q_surf ! changes in near-surface q

      !   ! adjust fluxes and atmospheric increments for
      !   ! implicit dependence on surface temperature

      !   flux_t        = flux_t      + dt_t_ca  *dhdt_surf
      !   flux_lw       = flux_lw     - dt_t_surf*drdt_surf
      !   Boundary%dt_t = f_t_delt_n  + dt_t_ca  *e_t_n

      !   where (Land%mask(:,:,1))
      !      flux_q                     = flux_q      + dt_q_surf*dedq_surf
      !      Boundary%dt_tr(:,:,isphum) = f_q_delt_n  + dt_q_surf*e_q_n
      !   elsewhere
      !      !flux_q                     = flux_q      + dt_t_surf*dedt_surf
      !      !Boundary%dt_tr(:,:,isphum) = f_q_delt_n  + dt_t_surf*e_q_n
      !   endwhere



   end subroutine flux_up_to_atmos

   !! Write out structured grid diagnostic history
   !! ============================================================================
   subroutine debug_diag(lm4_model)
      ! This is a quick and dirty diagnostic routine to write out some fields

      use mpp_domains_mod,  only : mpp_get_ntile_count
      use diag_manager_mod, only : diag_axis_init, register_static_field, &
         register_diag_field, diag_field_add_attribute, send_data
      use time_manager_mod,     only: time_type


      type(lm4_type), intent(inout) :: lm4_model

      ! local variables
      character(len=32) :: mod_name = 'lm4_dbug_diag'  ! diag module name for history
      logical           :: first_call = .true.
      real              :: missval    = -1.0e+20
      logical           :: used
      integer           :: i


      ! only run if first call,
      if (first_call) then
         first_call = .false.

         ! initialize output on structure grid, with cell_area

         if(mpp_get_ntile_count(lnd%sg_domain)==1) then
            ! grid has just one tile, so we assume that the grid is regular lat-lon
            ! define longitude axes and its edges
            id_lonb = diag_axis_init ('lonb', lnd%coord_glonb, 'degrees_E', 'X', 'longitude edges', &
               set_name=mod_name, domain2=lnd%sg_domain )
            id_lon  = diag_axis_init ('lon',  lnd%coord_glon, 'degrees_E', 'X',   'longitude', &
               set_name=mod_name,  edges=id_lonb, domain2=lnd%sg_domain )

            ! define latitude axes and its edges
            id_latb = diag_axis_init ('latb', lnd%coord_glatb, 'degrees_N', 'Y', 'latitude edges',  &
               set_name=mod_name,  domain2=lnd%sg_domain   )
            id_lat = diag_axis_init ('lat',  lnd%coord_glat, 'degrees_N', 'Y', 'latitude', &
               set_name=mod_name, edges=id_latb, domain2=lnd%sg_domain)
         else
            id_lon = diag_axis_init ( 'grid_xt', [(real(i),i=1,size(lnd%coord_glon))], 'degrees_E', 'X', &
               'T-cell longitude', set_name=mod_name,  domain2=lnd%sg_domain)
            id_lat = diag_axis_init ( 'grid_yt', [(real(i),i=1,size(lnd%coord_glat))], 'degrees_N', 'Y', &
               'T-cell latitude', set_name=mod_name,  domain2=lnd%sg_domain)
         endif

         ! register cell area on structured grid
         id_cellarea = register_static_field( mod_name, 'cell_area', (/id_lon, id_lat/), &
            'total area in grid cell', 'm2', missing_value=-1.0 )
         call diag_field_add_attribute(id_cellarea,'cell_methods','area: sum')

         associate ( &
            axes => (/ id_lon, id_lat /), &
            ltime => lm4_model%Time_land  &
            )

            ! register other fields on structured grid
            id_t_bot   = register_diag_field(mod_name, 't_bot', axes, ltime, 'bottom temperature', 'K', missing_value=missval )
            id_p_bot   = register_diag_field(mod_name, 'p_bot', axes, ltime, 'bottom pressure', 'Pa', missing_value=missval )
            id_z_bot   = register_diag_field(mod_name, 'z_bot', axes, ltime, 'bottom depth', 'm', missing_value=missval )
            id_u_bot   = register_diag_field(mod_name, 'u_bot', axes, ltime, 'bottom u velocity', 'm/s', missing_value=missval )
            id_v_bot   = register_diag_field(mod_name, 'v_bot', axes, ltime, 'bottom v velocity', 'm/s', missing_value=missval )
            id_q_bot   = register_diag_field(mod_name, 'q_bot', axes, ltime, 'bottom specific humidity', 'kg/kg', missing_value=missval )
            id_p_surf  = register_diag_field(mod_name, 'p_surf', axes, ltime, 'surface pressure', 'Pa', missing_value=missval )
            id_lprec   = register_diag_field(mod_name, 'lprec', axes, ltime, 'liquid precipitation', 'kg/m2/s', missing_value=missval )
            id_fprec   = register_diag_field(mod_name, 'fprec', axes, ltime, 'frozen precipitation', 'kg/m2/s', missing_value=missval )
            id_flux_lw = register_diag_field(mod_name, 'flux_lw', axes, ltime, 'longwave flux down', 'W/m2', missing_value=missval )
            id_swdn_vf = register_diag_field(mod_name, 'sw_down_vis_dif', axes, ltime,  'shortwave downwelling vis. diffuse radiation', 'W/m2', missing_value=missval )
            id_flux_sw_dn_vdf = register_diag_field(mod_name, 'flux_sw_down_vis_dif', axes, ltime, 'vis. diff. shortwave flux down', 'W/m2', missing_value=missval )
            id_flux_sw_dn_vr  = register_diag_field(mod_name, 'flux_sw_down_vis_dir', axes, ltime, 'vis. dir. shortwave flux down', 'W/m2', missing_value=missval )


         end associate

         ! send out static data
         if (id_cellarea > 0)       used = send_data(id_cellarea,       lnd%sg_cellarea,                           lm4_model%Time_land)

      endif ! first_call

      ! send out data to be written
      if (id_t_bot > 0)          used = send_data(id_t_bot,          lm4_model%atm_forc2d%t_bot,                lm4_model%Time_land)
      if (id_p_bot > 0)          used = send_data(id_p_bot,          lm4_model%atm_forc2d%p_bot,                lm4_model%Time_land)
      if (id_z_bot > 0)          used = send_data(id_z_bot,          lm4_model%atm_forc2d%z_bot,                lm4_model%Time_land)
      if (id_u_bot > 0)          used = send_data(id_u_bot,          lm4_model%atm_forc2d%u_bot,                lm4_model%Time_land)
      if (id_v_bot > 0)          used = send_data(id_v_bot,          lm4_model%atm_forc2d%v_bot,                lm4_model%Time_land)
      if (id_q_bot > 0)          used = send_data(id_q_bot,          lm4_model%atm_forc2d%q_bot,                lm4_model%Time_land)
      if (id_p_surf > 0)         used = send_data(id_p_surf,         lm4_model%atm_forc2d%p_surf,               lm4_model%Time_land)
      if (id_lprec > 0)          used = send_data(id_lprec,          lm4_model%atm_forc2d%lprec,                lm4_model%Time_land)
      if (id_fprec > 0)          used = send_data(id_fprec,          lm4_model%atm_forc2d%fprec,                lm4_model%Time_land)
      if (id_flux_lw > 0)        used = send_data(id_flux_lw,        lm4_model%atm_forc2d%flux_lw,              lm4_model%Time_land)
      if (id_swdn_vf > 0)        used = send_data(id_swdn_vf,        lm4_model%atm_forc2d%flux_sw_down_vis_dif, lm4_model%Time_land)
      if (id_flux_sw_dn_vdf > 0) used = send_data(id_flux_sw_dn_vdf, lm4_model%atm_forc2d%flux_sw_down_vis_dif, lm4_model%Time_land)
      if (id_flux_sw_dn_vr > 0)  used = send_data(id_flux_sw_dn_vr,  lm4_model%atm_forc2d%flux_sw_down_vis_dir, lm4_model%Time_land)




   end subroutine debug_diag





   !! Set up unstructured grid diagnostics
   ! initialize horizontal axes for land grid so that all sub-modules can use them,
   ! instead of creating their own
   !! ============================================================================
   subroutine land_diag_init(clonb, clatb, clon, clat, time, domain, id_band, id_ug)


      !Inputs/outputs
      real,dimension(:),intent(in) :: clonb   !<longitudes of grid cells vertices
      real,dimension(:),intent(in) :: clatb   !<latitudes of grid cells vertices
      real,dimension(:),intent(in) :: clon    !<Longitude of grid cell centers.
      real,dimension(:),intent(in) :: clat    !<Latitude of grid cell centers
      type(time_type),intent(in)   :: time    !<Initial time for diagnostic fields.
      type(domainUG), intent(in)   :: domain  !<
      integer,intent(out)          :: id_band !<"band" axis id.
      integer,intent(out)          :: id_ug   !<Unstructured axis id.

      ! ---- local vars ----------------------------------------------------------
      character(len=32) :: module_name = 'UG_dbug_diag'  ! diag module name for history

      integer :: nlon, nlat       ! sizes of respective axes
      integer             :: axes(1)        ! Array of axes for 1-D unstructured fields.
      integer             :: ug_dim_size    ! Size of the unstructured axis
      integer,allocatable :: ug_dim_data(:) ! Unstructured axis data.
      ! integer             :: id_lon, id_lonb
      ! integer             :: id_lat, id_latb
      integer :: i
      character(32) :: name       ! tracer name

      ! Register the unstructured axis for the unstructured domain.
      call mpp_get_UG_compute_domain(domain, size=ug_dim_size)
      if (.not. allocated(ug_dim_data)) then
         allocate(ug_dim_data(ug_dim_size))
      endif
      call mpp_get_UG_domain_grid_index(domain, ug_dim_data)
      !--- grid_index needs to be starting from 0.
      ug_dim_data = ug_dim_data - 1
      id_ug = diag_axis_init("grid_index",  real(ug_dim_data), "none", "U", long_name="grid indices", &
         set_name=trim(module_name), DomainU=domain, aux="geolon_t geolat_t")
      if (allocated(ug_dim_data)) then
         deallocate(ug_dim_data)
      endif

      ! Register horizontal axes that are required by the post-processing so that the output
      ! files can be "decompressed": converted from unstructured back to lon-lat or cubic sphere.
      ! The "grid_xt" and "grid_yt" axes should run from 1 to the total number of x- and
      ! y-points on cubic sphere face. It is assumed that all faces tiles contain the same
      ! number of x- and y-points.
      nlon = size(clon)
      nlat = size(clat)
      if(mpp_get_UG_domain_ntiles(lnd%ug_domain)==1) then
         ! grid has just one tile, so we assume that the grid is regular lat-lon
         ! define geographic axes and its edges
         id_lonb = diag_axis_init ('lonb', clonb, 'degrees_E', 'X', 'longitude edges', set_name=trim(module_name))
         id_lon  = diag_axis_init ('lon',  clon,  'degrees_E', 'X', 'longitude', set_name=trim(module_name),  edges=id_lonb)
         id_latb = diag_axis_init ('latb', clatb, 'degrees_N', 'Y', 'latitude edges', set_name=trim(module_name))
         id_lat  = diag_axis_init ('lat',  clat,  'degrees_N', 'Y', 'latitude', set_name=trim(module_name), edges=id_latb)
         ! add "compress" attribute to the unstructured grid axis
         call diag_axis_add_attribute(id_ug, "compress", "lat lon")
      else
         id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,nlon)/), 'degrees_E', 'X', &
            'T-cell longitude', set_name=trim(module_name) )
         id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,nlat)/), 'degrees_N', 'Y', &
            'T-cell latitude', set_name=trim(module_name) )
         ! add "compress" attribute to the unstructured grid axis
         call diag_axis_add_attribute(id_ug, "compress", "grid_yt grid_xt")
      endif

      id_band = diag_axis_init ('band',  (/1.0,2.0/), 'unitless', 'Z', 'spectral band', set_name=trim(module_name) )

      ! Set up an array of axes ids, for convenience.
      axes(1) = id_ug

      ! register auxiliary coordinate variables
      id_geolon_t = register_static_field ( module_name, 'geolon_t', axes, &
         'longitude of grid cell centers', 'degrees_E', missing_value = -1.0e+20 )
      id_geolat_t = register_static_field ( module_name, 'geolat_t', axes, &
         'latitude of grid cell centers', 'degrees_N', missing_value = -1.0e+20 )

      ! register static diagnostic fields
      id_landfrac = register_static_field ( module_name, 'land_frac', axes, &
         'fraction of land in grid cell','unitless', missing_value=-1.0, area=id_cellarea)
      call diag_field_add_attribute(id_landfrac,'ocean_fillvalue',0.0)

      ! register areas and fractions for the rest of the diagnostic fields
      call register_tiled_area_fields(module_name, axes, time, id_area, id_frac)

      ! set the default filter (for area and subsampling) for consequent calls to
      ! register_tiled_diag_field
      !call set_default_diag_filter('land')

      ! register regular (dynamic) diagnostic fields

      id_ntiles = register_tiled_diag_field(module_name,'ntiles', axes,  &
         time, 'number of tiles', 'unitless', missing_value=-1.0, op='sum')


      iug_q_atm       = register_tiled_diag_field(module_name, "q_atm"      , axes, time, "q_atm"      , "kg/kg", missing_value=-1.0e+20)
      iug_t_atm       = register_tiled_diag_field(module_name, "t_atm"      , axes, time, "t_atm"      , "K"    , missing_value=-1.0e+20)
      iug_u_atm       = register_tiled_diag_field(module_name, "u_atm"      , axes, time, "u_atm"      , "m/s"  , missing_value=-1.0e+20)
      iug_v_atm       = register_tiled_diag_field(module_name, "v_atm"      , axes, time, "v_atm"      , "m/s"  , missing_value=-1.0e+20)
      iug_p_atm       = register_tiled_diag_field(module_name, "p_atm"      , axes, time, "p_atm"      , "Pa"   , missing_value=-1.0e+20)
      iug_z_atm       = register_tiled_diag_field(module_name, "z_atm"      , axes, time, "z_atm"      , "m"    , missing_value=-1.0e+20)
      iug_p_surf      = register_tiled_diag_field(module_name, "p_surf"     , axes, time, "p_surf"     , "Pa"   , missing_value=-1.0e+20)
      iug_t_surf      = register_tiled_diag_field(module_name, "t_surf"     , axes, time, "t_surf"     , "K"    , missing_value=-1.0e+20)
      iug_t_ca        = register_tiled_diag_field(module_name, "t_ca"       , axes, time, "t_ca"       , "K"    , missing_value=-1.0e+20)
      iug_q_surf      = register_tiled_diag_field(module_name, "q_surf"     , axes, time, "q_surf"     , "kg/kg", missing_value=-1.0e+20)
      iug_rough_mom   = register_tiled_diag_field(module_name, "rough_mom"  , axes, time, "rough_mom"  , "m"    , missing_value=-1.0e+20)
      iug_rough_heat  = register_tiled_diag_field(module_name, "rough_heat" , axes, time, "rough_heat" , "m"    , missing_value=-1.0e+20)
      iug_rough_moist = register_tiled_diag_field(module_name, "rough_moist", axes, time, "rough_moist", "m"    , missing_value=-1.0e+20)
      iug_rough_scale = register_tiled_diag_field(module_name, "rough_scale", axes, time, "rough_scale", "m"    , missing_value=-1.0e+20)
      iug_gust        = register_tiled_diag_field(module_name, "gust"       , axes, time, "gust"       , "m/s"  , missing_value=-1.0e+20)

   end subroutine land_diag_init

   ! !! Write out the land structured grid diagnostics
   ! !! ============================================================================
   ! subroutine sg_send_data()

   !    type(diag_buff_type), intent(inout) :: sg_diag





   ! end subroutine sg_send_data

   !! Wrap up
   !! ===========================================================================
   subroutine end_driver()

      deallocate( &
         ex_flux_t, ex_flux_lw, ex_dhdt_surf, ex_dedt_surf, &
         ex_drdt_surf, ex_dhdt_atm, ex_tr_atm, ex_tr_surf, &
         ex_flux_tr, ex_dfdtr_surf, ex_dfdtr_atm, ex_e_tr_n, &
         ex_f_tr_delt_n, &
         ex_drag_q,    &   
         ex_cd_t,      &
         ex_cd_m,      &
         ex_b_star,    &
         ex_u_star,    &
         ex_wind,      &
         ex_z_atm               )

      
   end subroutine end_driver


end module lm4_driver
