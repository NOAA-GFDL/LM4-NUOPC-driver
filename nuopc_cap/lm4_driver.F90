!! Routines to prepare and run the land model
!! ============================================================================
module lm4_driver

   use mpp_domains_mod,    only: domain2d
   use lm4_type_mod,       only: lm4_type
   use land_data_mod,      only: land_data_type, atmos_land_boundary_type, lnd
   use land_tracers_mod,   only : isphum

   implicit none
   private


   integer                      :: im         ! horiz dimension
   integer                      :: date_init(6)

   type(lm4_type),           public :: lm4_model
   type(domain2D),           public :: land_domain

   public :: lm4_nml_read
   public :: init_driver, end_driver
   public :: sfc_boundary_layer, flux_down_from_atmos
   public :: debug_diag

   ! ---- namelist with default values ------------------------------------------
   logical :: old_dtaudv            = .false. !< The derivative of surface wind stress with respect to the zonal wind and meridional
   !! wind are approximated by the same tendency
   logical :: use_mixing_ratio      = .false. !< An option to provide capability to run the Manabe Climate form of the surface flux
   !! (coded for legacy purposes).
   logical :: do_simple             = .false.
   logical :: no_neg_q              = .false. !< If a_atm_in (specific humidity) is negative (because of numerical truncation),
   !! then override with 0.0
   logical :: alt_gustiness         = .false. !< An alternaive formulation for gustiness calculation.  A minimum bound on the wind
   !! speed used influx calculations, with the bound equal to gust_const

   real    :: gust_const            =  1.0    !< Constant for alternative gustiness calculation
   real    :: gust_min              =  0.0    !< Minimum gustiness used when alt_gustiness is .FALSE.

   ! --- namelist of vars originally from flux exchange nml
   real :: z_ref_heat =  2. !< Reference height (meters) for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q)

   namelist /surface_flux_nml/  no_neg_q,             &
      alt_gustiness,        &
      gust_const,           &
      gust_min,             &
      old_dtaudv,           &
      use_mixing_ratio,     &
      do_simple

   ! TODO: rename this nml?
   namelist /flux_exchange_nml/ z_ref_heat

   logical :: scale_precip_2d = .false.

   ! variables for between subroutines
   real, allocatable, dimension(:)   :: ex_flux_u

   ! integers for diag manager fields (TODO: clean up)
   integer :: id_cellarea
   ! fields to be written out
   integer :: id_swdn_vf, id_z_bot, id_t_bot, id_p_bot, &
              id_u_bot, id_v_bot, id_lprec, id_fprec, id_totprec, &
              id_flux_lw, id_flux_sw_dn_vdf, id_flux_sw_dn_vr
   ! other vars
   integer :: id_lon, id_lat, id_lonb, id_latb   

contains

   !! Read in lm4 namelist
   !! ============================================================================
   subroutine lm4_nml_read(lm4_model)

      use fms_mod,             only: check_nml_error, close_file, file_exist
      use mpp_mod,             only: mpp_pe, mpp_root_pe
#ifdef INTERNAL_FILE_NML
      use mpp_mod,             only: input_nml_file
#else
      use fms_mod,             only: open_namelist_file
#endif

      type(lm4_type),          intent(inout) :: lm4_model ! land model's variable type

      ! namelist variables for lm4
      ! ------------------------------------------
      integer           :: lm4_debug = 0        ! debug flag for lm4 (0=off, 1=low, 2=high)
      integer           :: npx = 0, npy = 0
      integer           :: ntiles = 0
      integer           :: layout(2) = (/0,0/)
      character(len=64) :: grid      = 'none'
      integer           :: blocksize = -1


      ! for namelist read
      integer :: unit, io, ierr
      namelist /lm4_nml/ grid, npx, npy, layout, ntiles, &
         blocksize, lm4_debug

      ! read in namelist

      if ( file_exist('input.nml')) then
#ifdef INTERNAL_FILE_NML
         read(input_nml_file, nml=lm4_nml, iostat=io)
         ierr = check_nml_error(io, 'lm4_nml')
#else
         unit = open_namelist_file ( )
         ierr=1
         do while (ierr /= 0)
            read(unit, nml=lm4_nml, iostat=io)
            ierr = check_nml_error(io,'lm4_nml')
         enddo
         call close_file(unit)
#endif
      endif

      lm4_model%nml%lm4_debug = lm4_debug
      lm4_model%nml%grid      = grid
      lm4_model%nml%blocksize = blocksize
      lm4_model%nml%npx       = npx
      lm4_model%nml%npy       = npy
      lm4_model%nml%layout    = layout
      lm4_model%nml%ntiles    = ntiles

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

      type (block_control_type), target   :: Lnd_block !  Block container
      integer                      ::   blocksize
      logical, save                :: block_message = .true.


      integer :: isc, iec, jsc, jec

      ! FMS domain creation:
      call domain_create(lm4_model%nml, land_domain)

      call mpp_get_compute_domain(land_domain,isc,iec,jsc,jec)
      ! use LM4's data type, using just a part of land_data_init
      !call mpp_get_compute_domain(lnd%sg_domain, lnd%is,lnd%ie,lnd%js,lnd%je)

      ! isc = lnd%is
      ! iec = lnd%ie
      ! jsc = lnd%js
      ! jec = lnd%je

      !im = (iec-isc+1)*(jec-jsc+1)

      ! ! Create blocks, but again, not currently using
      ! call define_blocks_packed('land_model', Lnd_block, isc, iec, jsc, jec, 1, &
      !    lm4_model%nml%blocksize, block_message)


      ! ! Restart read of sfc_data
      ! call sfc_prop_restart_read(lm4_model, land_domain, .false.)
      ! ! Transfer from sfcprop to model data
      ! call sfc_prop_transfer(lm4_model)


      allocate(ex_flux_u(im))

   end subroutine init_driver

   !! ============================================================================
   !! Adapted from GFDL atm_land_ice_flux_exchange,
   !! stripped down to be "land only" on
   !! unstructured grid, returns
   !! explicit fluxes as well as derivatives 
   !! used to compute an implicit flux correction.
   !! ============================================================================
   subroutine sfc_boundary_layer( dt,lm4_model )


      use sat_vapor_pres_mod, only: compute_qs
      use monin_obukhov_mod,  only: mo_profile

      real,                  intent(in)     :: dt        ! Time step
      type(lm4_type),        intent(inout)  :: lm4_model ! land model's variable type
      !type(land_data_type),  intent(inout)  :: Land ! A derived data type to specify land boundary data
   
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
         !ex_t_atm,              &
         !ex_p_atm,              &
         !ex_u_atm, ex_v_atm,    &
         ex_gust,               &
         ex_t_surf4,            &
         ex_u_surf = 0.0,       &
         ex_v_surf = 0.0,        &
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

      ! these originally had a tracer dimension, but we're not using them
      real, dimension(lnd%ls:lnd%le) :: &
         ex_tr_atm,  &
         ex_tr_surf,    & !< near-surface tracer fields
         ex_flux_tr,    & !< tracer fluxes
         ex_dfdtr_surf, & !< d(tracer flux)/d(surf tracer)
         ex_dfdtr_atm,  & !< d(tracer flux)/d(atm tracer)
         ex_e_tr_n,     & !< coefficient in implicit scheme
         ex_f_tr_delt_n   !< coefficient in implicit scheme
    
      !! -- These were originally allocatable:
      !!

      logical, dimension(lnd%ls:lnd%le) :: &
            ex_avail,     &   !< true where data on exchange grid are available
            ex_land           !< true if exchange grid cell is over land

      real, dimension(lnd%ls:lnd%le) :: &
            ex_t_surf   ,  &
            ex_t_surf_miz, &
            ex_p_surf   ,  &
            !ex_slp      ,  &
            ex_t_ca     ,  &
            ex_dhdt_surf,  &
            ex_dedt_surf,  &
            ex_dqsatdt_surf,  &
            ex_drdt_surf,  &
            ex_dhdt_atm ,  &
            ex_flux_t   ,  &
            ex_flux_lw  ,  &
            ex_drag_q   ,  &
            ex_f_t_delt_n, &

         ! MOD these were moved from local ! so they can be passed to flux down
         !ex_flux_u,    &
            ex_flux_v,    &
            ex_dtaudu_atm,&
            ex_dtaudv_atm,&
            ex_seawater,  &

         ! values added for LM3
            ex_cd_t     ,  &
            ex_cd_m     ,  &
            ex_b_star   ,  &
            ex_u_star   ,  &
            ex_wind     ,  &
            !ex_z_atm    ,  &

            ex_e_t_n    ,  &
            ex_e_q_n    ,  &

         !
            ex_albedo_fix,        &
            ex_albedo_vis_dir_fix,&
            ex_albedo_nir_dir_fix,&
            ex_albedo_vis_dif_fix,&
            ex_albedo_nir_dif_fix

    
      integer :: tr, n, m ! tracer indices
      integer :: i
      !  integer :: is,ie,l,j
      !  integer :: isc,iec,jsc,jec


      ! JP end
    
      ! prefill q_surf with q_bot. In original code, there were options to send/write
      ! before surface_flux call.
      ex_tr_surf(lnd%ls:lnd%le) = lm4_model%From_lnd%q_bot      

      ! this is mimicking the original code for land roughness vars
      ! TODO: review and optimize
      ex_rough_moist(lnd%ls:lnd%le) = Land%rough_heat
      ex_rough_scale(lnd%ls:lnd%le) = ex_rough_mom
      ex_rough_scale(lnd%ls:lnd%le) = Land%rough_scale

    call surface_flux_1d ( &
         ! inputs
         lm4_model%atm_forc%t_bot,                                      &
         lm4_model%atm_forc%q_bot,                                      & !! TODO: link q_bot var and tracer field
         lm4_model%atm_forc%u_bot, lm4_model%atm_forc%v_bot,            &
         lm4_model%atm_forc%p_bot, lm4_model%atm_forc%z_bot,            &
         lm4_model%atm_forc%p_surf, lm4_model%From_lnd%t_surf,          & !! surface T and P
         lm4_model%From_lnd%t_a,                                        & !! air temp at the canopy
         ! inout         
         ex_tr_surf(lnd%ls:lnd%le),                                     & !! TODO review surface Q (this is INOUT)
         ! more inputs         
         ex_u_surf(lnd%ls:lnd%le), ex_v_surf(lnd%ls:lnd%le),            & !! TODO review surface wind (= 0)
         Land%rough_mom, Land%rough_heat,                               & !! rough lengths
         ex_rough_moist(lnd%ls:lnd%le), ex_rough_scale(lnd%ls:lnd%le),  & !! moisture roughness length and scale factor
         ex_gust(lnd%ls:lnd%le),                                        & !! TODO gustiness 
         ! outputs
         ex_flux_t(lnd%ls:lnd%le), ex_flux_tr(lnd%ls:lnd%le,isphum),    & !! SH and Evap (Q) fluxes
         ex_flux_lw(lnd%ls:lnd%le), &
         ex_flux_u(lnd%ls:lnd%le), ex_flux_v(lnd%ls:lnd%le), ex_cd_m(lnd%ls:lnd%le),   ex_cd_t(lnd%ls:lnd%le), &
         ex_cd_q(lnd%ls:lnd%le),   ex_wind(lnd%ls:lnd%le),   ex_u_star(lnd%ls:lnd%le), ex_b_star(lnd%ls:lnd%le), &
         ex_q_star(lnd%ls:lnd%le), ex_dhdt_surf(lnd%ls:lnd%le), ex_dedt_surf(lnd%ls:lnd%le), &
         ex_dfdtr_surf(lnd%ls:lnd%le,isphum),  ex_drdt_surf(lnd%ls:lnd%le),  ex_dhdt_atm(lnd%ls:lnd%le), &
         ex_dfdtr_atm(lnd%ls:lnd%le,isphum),  ex_dtaudu_atm(lnd%ls:lnd%le), ex_dtaudv_atm(lnd%ls:lnd%le),       &
         dt,                                                            & !! timestep doesn't seem to be used
         .TRUE., .FALSE., .TRUE.                                        & ! Is land, Is seawater, Is ex. avail
         ) 


   !! ....
   zrefm = 10.0
   zrefh = z_ref_heat


   call mo_profile ( zrefm, zrefh, lm4_model%atm_forc%z_bot, ex_rough_mom(lnd%ls:lnd%le), &
      ex_rough_heat(lnd%ls:lnd%le), ex_rough_moist(lnd%ls:lnd%le),          &
      ex_u_star(lnd%ls:lnd%le), ex_b_star(lnd%ls:lnd%le), ex_q_star(lnd%ls:lnd%le),        &
      ex_del_m(lnd%ls:lnd%le), ex_del_h(lnd%ls:lnd%le), ex_del_q(lnd%ls:lnd%le), ex_avail(lnd%ls:lnd%le)  )

   do i = lnd%ls,lnd%le
      ex_u10(i) = 0.
      if(ex_avail(i)) then
         ex_ref_u(i) = ex_u_surf(i) + (lm4_model%atm_forc%u_bot(i)-ex_u_surf(i)) * ex_del_m(i)
         ex_ref_v(i) = ex_v_surf(i) + (lm4_model%atm_forc%v_bot(i)-ex_v_surf(i)) * ex_del_m(i)
         ex_u10(i) = sqrt(ex_ref_u(i)**2 + ex_ref_v(i)**2)
      endif
   enddo


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


   do i = lnd%ls,lnd%le
      ex_albedo_fix(i) = (1.0-ex_albedo(i)) / (1.0-ex_albedo_fix(i))
      ex_albedo_vis_dir_fix(i) = (1.0-ex_albedo_vis_dir(i)) / (1.0-ex_albedo_vis_dir_fix(i))
      ex_albedo_nir_dir_fix(i) = (1.0-ex_albedo_nir_dir(i)) / (1.0-ex_albedo_nir_dir_fix(i))
      ex_albedo_vis_dif_fix(i) = (1.0-ex_albedo_vis_dif(i)) / (1.0-ex_albedo_vis_dif_fix(i))
      ex_albedo_nir_dif_fix(i) = (1.0-ex_albedo_nir_dif(i)) / (1.0-ex_albedo_nir_dif_fix(i))
   enddo

      !=======================================================================
      ! [7] diagnostics section

   call mo_profile ( zrefm, zrefh, lm4_model%atm_forc%z_bot,   ex_rough_mom(lnd%ls:lnd%le), &
      ex_rough_heat(lnd%ls:lnd%le), ex_rough_moist(lnd%ls:lnd%le),          &
      ex_u_star(lnd%ls:lnd%le), ex_b_star(lnd%ls:lnd%le), ex_q_star(lnd%ls:lnd%le),        &
      ex_del_m(lnd%ls:lnd%le), ex_del_h(lnd%ls:lnd%le), ex_del_q(lnd%ls:lnd%le), ex_avail(lnd%ls:lnd%le)  )

   !    ------- reference relative humidity -----------
   !cjg     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
   !cjg          id_rh_ref_cmip > 0 .or. &
   !cjg          id_q_ref > 0 .or. id_q_ref_land >0 ) then
   do i = lnd%ls,lnd%le
      ex_ref(i) = 1.0e-06
      if (ex_avail(i)) &
         ex_ref(i)   = ex_tr_surf(i,isphum) + (ex_tr_atm(i,isphum)-ex_tr_surf(i,isphum)) * ex_del_q(i)
   enddo

   do i = lnd%ls,lnd%le
      ex_t_ref(i) = 200.
      if(ex_avail(i)) &
         ex_t_ref(i) = ex_t_ca(i) + (lm4_model%atm_forc%t_bot(i)-ex_t_ca(i)) * ex_del_h(i)
   enddo
   call compute_qs (ex_t_ref(lnd%ls:lnd%le), ex_p_surf(lnd%ls:lnd%le), ex_qs_ref(lnd%ls:lnd%le), q = ex_ref(lnd%ls:lnd%le))
   call compute_qs (ex_t_ref(lnd%ls:lnd%le), ex_p_surf(lnd%ls:lnd%le), ex_qs_ref_cmip(lnd%ls:lnd%le),  &
      q = ex_ref(lnd%ls:lnd%le), es_over_liq_and_ice = .true.)
   do i = lnd%ls,lnd%le
      if(ex_avail(i)) then
         ! remove cap on relative humidity -- this mod requested by cjg, ljd
         !RSH    ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)
         ex_ref2(i)   = 100.*ex_ref(i)/ex_qs_ref_cmip(i)
         ex_ref(i)    = 100.*ex_ref(i)/ex_qs_ref(i)
      endif
   enddo

   ! lots of send_data stuff originally here, removed
   ! TODO: get diag history write back in

   end subroutine sfc_boundary_layer

   ! adapted land-only surface_flux_1d from FMS coupler
   ! ============================================================================
   subroutine surface_flux_1d (                                           &
      t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
      p_surf,    t_surf,     t_ca,      q_surf,                         &
      u_surf,    v_surf,                                                &
      rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
      flux_t, flux_q, flux_r, flux_u, flux_v,                           &
      cd_m,      cd_t,       cd_q,                                      &
      w_atm,     u_star,     b_star,     q_star,                        &
      dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
      dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
      dt,        land,      seawater,     avail  )

      use constants_mod, only : rdgas, rvgas, kappa, cp_air, stefan
      ! note if using -DGFS_PHYS or not
      use sat_vapor_pres_mod, only: escomp
      use   monin_obukhov_mod, only: mo_drag


      ! ---- arguments -----------------------------------------------------------
      logical, intent(in), dimension(:) :: land, & !< Indicates where land exists (.TRUE. if exchange cell is on land
         seawater, & !< Indicates where liquid ocean water exists (.TRUE. if exchange cell is on liquid ocean water)
         avail !< .TRUE. where the exchange cell is active
      real, intent(in),  dimension(:) :: t_atm, & !< Air temp lowest atmospheric level.
         q_atm_in, & !< Mixing ratio at lowest atmospheric level (kg/kg).
         u_atm, & !< Zonal wind velocity at lowest atmospheric level.
         v_atm, & !< Meridional wind velocity at lowest atmospheric level.
         p_atm, & !< Pressure lowest atmospheric level.
         z_atm, & !< Height lowest atmospheric level.
         t_ca, & !< Air temp at the canopy
         p_surf, & !< Pressure at the Earth's surface
         t_surf, & !< Temp at the Earth's surface
         u_surf, & !< Zonal wind velocity at the Earth's surface
         v_surf, & !< Meridional wind velocity at the Earth's surface
         rough_mom, & !< Momentum roughness length
         rough_heat, & !< Heat roughness length
         rough_moist, & !< Moisture roughness length
         rough_scale, & !< Scale factor used to topographic roughness calculation
         gust !< Gustiness factor
      real, intent(out), dimension(:) :: flux_t, & !< Sensible heat flux
         flux_q, & !< Evaporative water flux
         flux_r, & !< Radiative energy flux
         flux_u, & !< Zonal momentum flux
         flux_v, & !< Meridional momentum flux
         dhdt_surf, & !< Sensible heat flux temperature sensitivity
         dedt_surf, & !< Moisture flux temperature sensitivity
         dedq_surf, & !< Moisture flux humidity sensitivity
         drdt_surf, & !< Radiative energy flux temperature sensitivity
         dhdt_atm, & !< Derivative of sensible heat flux over temp at the lowest atmos level
         dedq_atm, & !< Derivative of water vapor flux over temp at the lowest atmos level
         dtaudu_atm, & !< Derivative of zonal wind stress with respect to the lowest level zonal wind speed of the atmos
         dtaudv_atm, & !< Derivative of meridional wind stress with respect to the lowest level meridional wind speed of the atmos
         w_atm, & !< Absolute wind at the lowest atmospheric level
         u_star, & !< Turbulent velocity scale
         b_star, & !< Turbulent buoyant scale
         q_star, & !< Turbulent moisture scale
         cd_m, & !< Momentum exchange coefficient
         cd_t, & ! Heat exchange coefficient
         cd_q !< Moisture exchange coefficient
      real, intent(inout), dimension(:) :: q_surf !< Mixing ratio at the Earth's surface (kg/kg)
      real, intent(in) :: dt !< Time step (it is not used presently)

      ! ---- local constants -----------------------------------------------------
      ! temperature increment and its reciprocal value for comp. of derivatives
      real, parameter:: del_temp=0.1, del_temp_inv=1.0/del_temp

      ! ---- local vars ----------------------------------------------------------
      real, dimension(size(t_atm(:))) ::                          &
         thv_atm,  th_atm,   tv_atm,    thv_surf,            &
         e_sat,    e_sat1,   q_sat,     q_sat1,    p_ratio,  &
         t_surf0,  t_surf1,  u_dif,     v_dif,               &
         rho_drag, drag_t,   drag_m,    drag_q,    rho,      &
         q_atm,    q_surf0,  dw_atmdu,  dw_atmdv,  w_gust,   &
         zu,       zt,       zq

      integer :: i, nbad

      real, parameter :: d622 = rdgas/rvgas
      real, parameter :: d378   = 1.-d622
      real, parameter :: d608   = d378/d622
      !real, parameter :: STEFAN  = 5.6734e-8_r8_kind !< Stefan-Boltzmann constant [W/m^2/deg^4]

      ! if (.not. module_is_initialized) &
      !      call mpp_error(FATAL, "surface_flux_1d: surface_flux_init is not called")

      !---- use local value of surf temp ----

      t_surf0 = 200.   !  avoids out-of-bounds in es lookup
      where (avail)
         where (land)
            t_surf0 = t_ca
         elsewhere
            t_surf0 = t_surf
         endwhere
      endwhere

      t_surf1 = t_surf0 + del_temp

      call escomp ( t_surf0, e_sat  )  ! saturation vapor pressure
      call escomp ( t_surf1, e_sat1 )  ! perturbed  vapor pressure

      if(use_mixing_ratio) then
         ! surface mixing ratio at saturation
         q_sat   = d622*e_sat /(p_surf-e_sat )
         q_sat1  = d622*e_sat1/(p_surf-e_sat1)
      elseif(do_simple) then                  !rif:(09/02/09)
         q_sat   = d622*e_sat / p_surf
         q_sat1  = d622*e_sat1/ p_surf
      else
         ! surface specific humidity at saturation
         q_sat   = d622*e_sat /(p_surf-d378*e_sat )
         q_sat1  = d622*e_sat1/(p_surf-d378*e_sat1)
      endif

      ! initilaize surface air humidity according to surface type
      where (land)
         q_surf0 = q_surf ! land calculates it
      elsewhere
         q_surf0 = q_sat  ! everything else assumes saturated sfc humidity
      endwhere

      !! Seawater, not needed for land
      !if (raoult_sat_vap) where (seawater) q_surf0 = 0.98 * q_surf0

      ! check for negative atmospheric humidities
      where(avail) q_atm = q_atm_in
      if(no_neg_q) then
         where(avail .and. q_atm_in < 0.0) q_atm = 0.0
      endif

      ! generate information needed by monin_obukhov
      where (avail)
         p_ratio = (p_surf/p_atm)**kappa

         tv_atm  = t_atm  * (1.0 + d608*q_atm)     ! virtual temperature
         th_atm  = t_atm  * p_ratio                ! potential T, using p_surf as refernce
         thv_atm = tv_atm * p_ratio                ! virt. potential T, using p_surf as reference
         thv_surf= t_surf0 * (1.0 + d608*q_surf0 ) ! surface virtual (potential) T
         !     thv_surf= t_surf0                        ! surface virtual (potential) T -- just for testing tun off the q_surf

         u_dif = u_surf - u_atm                    ! velocity components relative to surface
         v_dif = v_surf - v_atm
      endwhere

      if(alt_gustiness) then
         do i = 1, size(avail)
            if (.not.avail(i)) cycle
            w_atm(i) = max(sqrt(u_dif(i)**2 + v_dif(i)**2), gust_const)
            ! derivatives of surface wind w.r.t. atm. wind components
            if(w_atm(i) > gust_const) then
               dw_atmdu(i) = u_dif(i)/w_atm(i)
               dw_atmdv(i) = v_dif(i)/w_atm(i)
            else
               dw_atmdu(i) = 0.0
               dw_atmdv(i) = 0.0
            endif
         enddo
      else
         if (gust_min > 0.0) then
            where(avail)
               w_gust = max(gust,gust_min) ! minimum gustiness
            end where
         else
            where(avail)
               w_gust = gust
            end where
         endif

         where(avail)
            w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + w_gust*w_gust)
            ! derivatives of surface wind w.r.t. atm. wind components
            dw_atmdu = u_dif/w_atm
            dw_atmdv = v_dif/w_atm
         endwhere
      endif

      !  monin-obukhov similarity theory
      call mo_drag (thv_atm, thv_surf, z_atm,                  &
         rough_mom, rough_heat, rough_moist, w_atm,          &
         cd_m, cd_t, cd_q, u_star, b_star, avail             )

      !! removed ocean flux calculation here

      where (avail)
         ! scale momentum drag coefficient on orographic roughness
         cd_m = cd_m*(log(z_atm/rough_mom+1)/log(z_atm/rough_scale+1))**2
         ! surface layer drag coefficients
         drag_t = cd_t * w_atm
         drag_q = cd_q * w_atm
         drag_m = cd_m * w_atm

         ! density
         rho = p_atm / (rdgas * tv_atm)

         ! sensible heat flux
         rho_drag = cp_air * drag_t * rho
         flux_t = rho_drag * (t_surf0 - th_atm)  ! flux of sensible heat (W/m**2)
         dhdt_surf =  rho_drag                   ! d(sensible heat flux)/d(surface temperature)
         dhdt_atm  = -rho_drag*p_ratio           ! d(sensible heat flux)/d(atmos temperature)

         ! evaporation
         rho_drag  =  drag_q * rho
         flux_q    =  rho_drag * (q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))

         where (land)
            dedq_surf = rho_drag
            dedt_surf = 0
         elsewhere
            dedq_surf = 0
            dedt_surf =  rho_drag * (q_sat1 - q_sat) *del_temp_inv
         endwhere

         dedq_atm  = -rho_drag   ! d(latent heat flux)/d(atmospheric mixing ratio)

         q_star = flux_q / (u_star * rho)             ! moisture scale
         ! ask Chris and Steve K if we still want to keep this for diagnostics
         q_surf = q_atm + flux_q / (rho*cd_q*w_atm)   ! surface specific humidity

         ! upward long wave radiation
         flux_r    =   stefan*t_surf**4               ! (W/m**2)
         drdt_surf = 4*stefan*t_surf**3               ! d(upward longwave)/d(surface temperature)

         ! stresses
         rho_drag   = drag_m * rho
         flux_u     = rho_drag * u_dif   ! zonal      component of stress (Nt/m**2)
         flux_v     = rho_drag * v_dif   ! meridional component of stress

      elsewhere
         ! zero-out un-available data in output only fields
         flux_t     = 0.0
         flux_q     = 0.0
         flux_r     = 0.0
         flux_u     = 0.0
         flux_v     = 0.0
         dhdt_surf  = 0.0
         dedt_surf  = 0.0
         dedq_surf  = 0.0
         drdt_surf  = 0.0
         dhdt_atm   = 0.0
         dedq_atm   = 0.0
         u_star     = 0.0
         b_star     = 0.0
         q_star     = 0.0
         q_surf     = 0.0
         w_atm      = 0.0
      endwhere

      ! calculate d(stress component)/d(atmos wind component)
      dtaudu_atm = 0.0
      dtaudv_atm = 0.0
      if (old_dtaudv) then
         where(avail)
            dtaudv_atm = -rho_drag
            dtaudu_atm = -rho_drag
         endwhere
      else
         where(avail)
            dtaudu_atm = -cd_m*rho*(dw_atmdu*u_dif + w_atm)
            dtaudv_atm = -cd_m*rho*(dw_atmdv*v_dif + w_atm)
         endwhere
      endif

   end subroutine surface_flux_1d

   ! ----------------------------------------

   subroutine flux_down_from_atmos(Land)

      !type(atmos_data_type), intent(inout)  :: Atm  !< A derived data type to specify atmosphere boundary data
      type(land_data_type),  intent(inout)  :: Land !< A derived data type to specify land boundary data

      real, dimension(im) :: ex_flux_sw, ex_flux_lwd, &
         ex_flux_sw_dir,  &
         ex_flux_sw_dif,  &
         ex_flux_sw_down_vis_dir, ex_flux_sw_down_total_dir,  &
         ex_flux_sw_down_vis_dif, ex_flux_sw_down_total_dif,  &
         ex_flux_sw_vis, &
         ex_flux_sw_vis_dir, &
         ex_flux_sw_vis_dif

      integer :: is, ie, l, i,j

      !! JP disabled for now. TODO: Fix this
      ! ! Assume land grid is same as atm grid
      ! if (scale_precip_2d) then
      !    !call mpp_get_compute_domain(Atm%Domain, is_atm, ie_atm, js_atm, je_atm)
      !    !call data_override ('ATM', 'precip_scale2d',    frac_precip,   Time)
      !    do j = lnd%js,lnd%je
      !       do i = lnd%is,lnd%ie
      !          Atm%lprec(i,j) = Atm%lprec(i,j)*frac_precip(i,j)
      !       enddo
      !    enddo
      ! endif


      !! JP disabled for now. TODO: Fix this
      ! ! Assume land grid is same as atm grid
      ! !if (partition_fprec_from_lprec .and. Atm%pe) then
      ! if (partition_fprec_from_lprec) then
      !    !call mpp_get_compute_domain(Atm%Domain, is_atm, ie_atm, js_atm, je_atm)
      !    do j = lnd%js,lnd%je
      !       do i = lnd%is,lnd%ie
      !          if (Atm%t_bot(i,j) < tfreeze) then
      !             Atm%fprec(i,j) = Atm%lprec(i,j)
      !             Atm%lprec(i,j) = 0.0
      !          endif
      !       enddo
      !    enddo
      ! endif

      ! MOD update stresses using atmos delta's but derivatives on exchange grid
      !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_flux_u,ex_delta_u, &
      !$OMP                                  ex_dtaudu_atm,ex_dtaudv_atm,ex_flux_v,ex_delta_v )     &
      !$OMP                          private(is,ie)
      !do l = 1, my_nblocks
      ! is=block_start(l)
      ! ie=block_end(l)
      ! do i = 1, im
      !    ex_flux_u(i) = ex_flux_u(i) + ex_delta_u(i)*ex_dtaudu_atm(i)
      !    ex_flux_v(i) = ex_flux_v(i) + ex_delta_v(i)*ex_dtaudv_atm(i)
      ! enddo
      !enddo
      !! JP: DO I NEED THIS? ^^

      !! JP: IGNORE TEMP FOR NOW
      ! !---- adjust sw flux for albedo variations on exch grid ----
      ! !---- adjust 4 categories (vis/nir dir/dif) separately  ----
      ! ! do l = 1, my_nblocks
      ! !    is=block_start(l)
      ! !    ie=block_end(l)
      !    do i = 1, im
      !       ex_flux_sw_dir(i) = ex_flux_sw_dir(i) - ex_flux_sw_vis_dir(i)     ! temporarily nir/dir
      !       ex_flux_sw_dir(i) = ex_flux_sw_dir(i) * ex_albedo_nir_dir_fix(i)  ! fix nir/dir
      !       ex_flux_sw_vis_dir(i) = ex_flux_sw_vis_dir(i) * ex_albedo_vis_dir_fix(i) ! fix vis/dir
      !       ex_flux_sw_dir(i) = ex_flux_sw_dir(i) + ex_flux_sw_vis_dir(i)     ! back to total dir

      !       ex_flux_sw_dif(i) = ex_flux_sw_dif(i) - ex_flux_sw_vis_dif(i)     ! temporarily nir/dif
      !       ex_flux_sw_dif(i) = ex_flux_sw_dif(i) * ex_albedo_nir_dif_fix(i)  ! fix nir/dif
      !       ex_flux_sw_vis_dif(i) = ex_flux_sw_vis_dif(i) * ex_albedo_vis_dif_fix(i) ! fix vis/dif
      !       ex_flux_sw_dif(i) = ex_flux_sw_dif(i) + ex_flux_sw_vis_dif(i)     ! back to total dif

      !       ex_flux_sw_vis(i) = ex_flux_sw_vis_dir(i) + ex_flux_sw_vis_dif(i) ! legacy, remove later
      !       ex_flux_sw(i)     = ex_flux_sw_dir(i)     + ex_flux_sw_dif(i)     ! legacy, remove later
      !    enddo
      ! ! enddo

      !----- adjust fluxes for implicit dependence on atmosphere ----
      ! is this needed to copy over?

      ! JP TEMP DISABLE FOR NOW
      ! do i = 1, im
      !       !----- compute net longwave flux (down-up) -----
      !       ! (note: lw up already in ex_flux_lw)
      !       ex_flux_lw(i) = ex_flux_lwd(i) - ex_flux_lw(i)
      !       if (ex_avail(i) ) then

      !          ! temperature
      !          ex_gamma(i)      =  1./ (1.0 - ex_dtmass(i)*(ex_dflux_t(i) + ex_dhdt_atm(i)*cp_inv))
      !          ex_e_t_n(i)      =  ex_dtmass(i)*ex_dhdt_surf(i)*cp_inv*ex_gamma(i)
      !          ex_f_t_delt_n(i) = (ex_delta_t(i) + ex_dtmass(i) * ex_flux_t(i)*cp_inv) * ex_gamma(i)

      !          ex_flux_t (i)    =  ex_flux_t(i)        + ex_dhdt_atm(i) * ex_f_t_delt_n(i)
      !          ex_dhdt_surf(i)  =  ex_dhdt_surf(i)     + ex_dhdt_atm(i) * ex_e_t_n(i)
      !          ! moisture
      !          !     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_q + ex_dedq_atm))
      !          ! here it looks like two derivatives with different units are added together,
      !          ! but in fact they are not: ex_dedt_surf and ex_dedq_surf defined in complimentary
      !          ! regions of exchange grid, so that if one of them is not zero the other is, and
      !          ! vice versa.
      !          !     ex_e_q_n      =  ex_dtmass*(ex_dedt_surf+ex_dedq_surf) * ex_gamma
      !          !     ex_f_q_delt_n = (ex_delta_q  + ex_dtmass * ex_flux_q) * ex_gamma
      !          !     ex_flux_q     =  ex_flux_q    + ex_dedq_atm * ex_f_q_delt_n
      !          !     ex_dedt_surf  =  ex_dedt_surf + ex_dedq_atm * ex_e_q_n
      !          !     ex_dedq_surf  =  ex_dedq_surf + ex_dedq_atm * ex_e_q_n
      !          ! moisture vs. surface temperture, assuming saturation
      !          ex_gamma(i)   =  1.0 / (1.0 - ex_dtmass(i)*(ex_dflux_tr(i,isphum) + ex_dfdtr_atm(i,isphum)))
      !          ex_e_q_n(i)      =  ex_dtmass(i) * ex_dedt_surf(i) * ex_gamma(i)
      !          ex_dedt_surf(i)  =  ex_dedt_surf(i) + ex_dfdtr_atm(i,isphum) * ex_e_q_n(i)
      !          do tr = 1,n_exch_tr
      !             ex_gamma(i)   =  1.0 / (1.0 - ex_dtmass(i)*(ex_dflux_tr(i,tr) + ex_dfdtr_atm(i,tr)))

      !             ex_e_tr_n(i,tr)      =  ex_dtmass(i)*ex_dfdtr_surf(i,tr)*ex_gamma(i)
      !             ex_f_tr_delt_n(i,tr) = (ex_delta_tr(i,tr)+ex_dtmass(i)*ex_flux_tr(i,tr))*ex_gamma(i)

      !             ex_flux_tr(i,tr)     =  ex_flux_tr(i,tr) + ex_dfdtr_atm(i,tr)*ex_f_tr_delt_n(i,tr)
      !             ex_dfdtr_surf(i,tr)  =  ex_dfdtr_surf(i,tr) + ex_dfdtr_atm(i,tr)*ex_e_tr_n(i,tr)
      !          enddo
      !       endif
      !    enddo ! i

   end subroutine flux_down_from_atmos

   ! ! ----------------------------------------

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

   ! ----------------------------------------
   ! Write out diagnostic history
   ! ----------------------------------------
   subroutine debug_diag(lm4_model)
      ! This is a quick and dirty diagnostic routine to write out some fields

      use mpp_domains_mod,  only : mpp_get_ntile_count
      use diag_manager_mod, only : diag_axis_init, register_static_field, &
         register_diag_field, diag_field_add_attribute, send_data
      use time_manager_mod,     only: time_type


      type(lm4_type), intent(inout) :: lm4_model
      !integer, intent(out) :: id_cellarea

      integer :: i
      logical :: used
      logical :: first_call = .true.

      ! only run if first call, 
      if (first_call) then
         first_call = .false.

         ! initialize output on structure grid, with cell_area

         if(mpp_get_ntile_count(lnd%sg_domain)==1) then
            ! grid has just one tile, so we assume that the grid is regular lat-lon
            ! define longitude axes and its edges
            id_lonb = diag_axis_init ('lonb', lnd%coord_glonb, 'degrees_E', 'X', 'longitude edges', &
               set_name='lm4_dbug_diag', domain2=lnd%sg_domain )
            id_lon  = diag_axis_init ('lon',  lnd%coord_glon, 'degrees_E', 'X',   'longitude', &
               set_name='lm4_dbug_diag',  edges=id_lonb, domain2=lnd%sg_domain )

            ! define latitude axes and its edges
            id_latb = diag_axis_init ('latb', lnd%coord_glatb, 'degrees_N', 'Y', 'latitude edges',  &
               set_name='lm4_dbug_diag',  domain2=lnd%sg_domain   )
            id_lat = diag_axis_init ('lat',  lnd%coord_glat, 'degrees_N', 'Y', 'latitude', &
               set_name='lm4_dbug_diag', edges=id_latb, domain2=lnd%sg_domain)
         else
            id_lon = diag_axis_init ( 'grid_xt', [(real(i),i=1,size(lnd%coord_glon))], 'degrees_E', 'X', &
               'T-cell longitude', set_name='lm4_dbug_diag',  domain2=lnd%sg_domain)
            id_lat = diag_axis_init ( 'grid_yt', [(real(i),i=1,size(lnd%coord_glat))], 'degrees_N', 'Y', &
               'T-cell latitude', set_name='lm4_dbug_diag',  domain2=lnd%sg_domain)
         endif

         ! register cell area on structured grid
         id_cellarea = register_static_field( 'lm4_dbug_diag', 'cell_area', (/id_lon, id_lat/), &
            'total area in grid cell', 'm2', missing_value=-1.0 )
         call diag_field_add_attribute(id_cellarea,'cell_methods','area: sum')

         ! register other fields on structured grid
         id_swdn_vf = register_diag_field( 'lm4_dbug_diag', 'sw_down_vis_dif', (/id_lon, id_lat/), lm4_model%Time_land,  'shortwave downwelling vis. diffuse radiation', 'W/m2', missing_value=-1.0e+20 )
         id_t_bot = register_diag_field( 'lm4_dbug_diag', 't_bot', (/id_lon, id_lat/), lm4_model%Time_land, 'bottom temperature', 'K', missing_value=-1.0e+20 )
         id_p_bot = register_diag_field( 'lm4_dbug_diag', 'p_bot', (/id_lon, id_lat/), lm4_model%Time_land, 'bottom pressure', 'Pa', missing_value=-1.0e+20 )
         id_z_bot = register_diag_field( 'lm4_dbug_diag', 'z_bot', (/id_lon, id_lat/), lm4_model%Time_land, 'bottom depth', 'm', missing_value=-1.0e+20 )
         id_u_bot = register_diag_field( 'lm4_dbug_diag', 'u_bot', (/id_lon, id_lat/), lm4_model%Time_land, 'bottom u velocity', 'm/s', missing_value=-1.0e+20 )
         id_v_bot = register_diag_field( 'lm4_dbug_diag', 'v_bot', (/id_lon, id_lat/), lm4_model%Time_land, 'bottom v velocity', 'm/s', missing_value=-1.0e+20 )
         id_totprec = register_diag_field( 'lm4_dbug_diag', 'totprec', (/id_lon, id_lat/), lm4_model%Time_land, 'total precipitation', 'kg/m2/s', missing_value=-1.0e+20 )
         id_lprec = register_diag_field( 'lm4_dbug_diag', 'lprec', (/id_lon, id_lat/), lm4_model%Time_land, 'liquid precipitation', 'kg/m2/s', missing_value=-1.0e+20 )
         id_fprec = register_diag_field( 'lm4_dbug_diag', 'fprec', (/id_lon, id_lat/), lm4_model%Time_land, 'frozen precipitation', 'kg/m2/s', missing_value=-1.0e+20 )
         id_flux_lw = register_diag_field( 'lm4_dbug_diag', 'flux_lw', (/id_lon, id_lat/), lm4_model%Time_land, 'longwave flux down', 'W/m2', missing_value=-1.0e+20 )
         id_flux_sw_dn_vdf = register_diag_field( 'lm4_dbug_diag', 'flux_sw_down_vis_dif', (/id_lon, id_lat/), lm4_model%Time_land, 'vis. diff. shortwave flux down', 'W/m2', missing_value=-1.0e+20 )
         id_flux_sw_dn_vr = register_diag_field( 'lm4_dbug_diag', 'flux_sw_down_vis_dir', (/id_lon, id_lat/), lm4_model%Time_land, 'vis. dir. shortwave flux down', 'W/m2', missing_value=-1.0e+20 )
      endif ! first_call

         ! send out data to be written ----------------------------
         if (id_cellarea > 0)       used = send_data(id_cellarea,       lnd%sg_cellarea,                           lm4_model%Time_land)
         if (id_swdn_vf > 0)        used = send_data(id_swdn_vf,        lm4_model%atm_forc2d%flux_sw_down_vis_dif, lm4_model%Time_land)
         if (id_t_bot > 0)          used = send_data(id_t_bot,          lm4_model%atm_forc2d%t_bot,                lm4_model%Time_land)
         if (id_p_bot > 0)          used = send_data(id_p_bot,          lm4_model%atm_forc2d%p_bot,                lm4_model%Time_land)
         if (id_z_bot > 0)          used = send_data(id_z_bot,          lm4_model%atm_forc2d%z_bot,                lm4_model%Time_land)
         if (id_u_bot > 0)          used = send_data(id_u_bot,          lm4_model%atm_forc2d%u_bot,                lm4_model%Time_land)
         if (id_v_bot > 0)          used = send_data(id_v_bot,          lm4_model%atm_forc2d%v_bot,                lm4_model%Time_land)
         if (id_totprec > 0)        used = send_data(id_totprec,        lm4_model%atm_forc2d%totprec,              lm4_model%Time_land)
         if (id_lprec > 0)          used = send_data(id_lprec,          lm4_model%atm_forc2d%lprec,                lm4_model%Time_land)
         if (id_fprec > 0)          used = send_data(id_fprec,          lm4_model%atm_forc2d%fprec,                lm4_model%Time_land)
         if (id_flux_lw > 0)        used = send_data(id_flux_lw,        lm4_model%atm_forc2d%flux_lw,              lm4_model%Time_land)
         if (id_flux_sw_dn_vdf > 0) used = send_data(id_flux_sw_dn_vdf, lm4_model%atm_forc2d%flux_sw_down_vis_dif, lm4_model%Time_land)
         if (id_flux_sw_dn_vr > 0)  used = send_data(id_flux_sw_dn_vr,  lm4_model%atm_forc2d%flux_sw_down_vis_dir, lm4_model%Time_land)




   end subroutine debug_diag


   subroutine end_driver()

      deallocate(ex_flux_u)

   end subroutine end_driver


end module lm4_driver
