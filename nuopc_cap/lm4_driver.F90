module lm4_driver

  !use machine, only: kind_phys

  use proc_bounds,        only: procbounds_type, control_init_type
  use mpp_domains_mod,    only: domain2d
  use lm4_type_mod
  use land_model_mod,          only: land_data_type

  implicit none
  private

  type(lm4_type),           public :: lm4_model
  type(domain2D),           public :: land_domain
  type(control_init_type),  public :: ctrl_init

  public :: init_driver, run_driver

  ! ---- namelist with default values ------------------------------------------
  logical :: use_mixing_ratio      = .false. !< An option to provide capability to run the Manabe Climate form of the surface flux
  !! (coded for legacy purposes).
  logical :: do_simple             = .false.
  logical :: no_neg_q              = .false. !< If a_atm_in (specific humidity) is negative (because of numerical truncation),
  !! then override with 0.0
  logical :: alt_gustiness         = .false. !< An alternaive formulation for gustiness calculation.  A minimum bound on the wind
                                           !! speed used influx calculations, with the bound equal to gust_const

  real    :: gust_const            =  1.0 !< Constant for alternative gustiness calculation
  real    :: gust_min              =  0.0 !< Minimum gustiness used when alt_gustiness is .FALSE.
  
  namelist /surface_flux_nml/  no_neg_q,             &
                               alt_gustiness,        &
                               gust_const,           &
                               gust_min,             &
                               use_mixing_ratio,     &
                               do_simple

  
contains

  !subroutine init_driver(procbounds)
  subroutine init_driver(ctrl_init)

    use mpp_domains_mod,    only: domain2d, mpp_get_compute_domain
    use mpp_mod,            only: mpp_pe, mpp_root_pe
    use land_domain_mod,    only: domain_create
    use block_control_mod,  only: block_control_type, define_blocks_packed
    !use land_restart_mod,   only: sfc_prop_restart_read, sfc_prop_transfer
    type(control_init_type), intent(out)  ::   ctrl_init

    ! ---------------
    ! local

    type (block_control_type), target   :: Lnd_block !  Block container
    integer                      ::   blocksize
    logical, save                :: block_message = .true.


    integer                      ::   im         ! horiz dimension
    integer :: isc, iec, jsc, jec


    call ctrl_init%init()

    if (mpp_pe() == mpp_root_pe()) then
       write(*,*) 'ctrl_init%grid: '     ,ctrl_init%grid
       write(*,*) 'ctrl_init%npx: '      ,ctrl_init%npx
       write(*,*) 'ctrl_init%npy: '      ,ctrl_init%npy
       write(*,*) 'ctrl_init%layout: '   ,ctrl_init%layout
       write(*,*) 'ctrl_init%ntiles: '   ,ctrl_init%ntiles
       write(*,*) 'ctrl_init%blocksize: ',ctrl_init%blocksize
    end if

    ! FMS domain creation:
    call domain_create(ctrl_init, land_domain)

    ! Create blocking a la FV3, but not currently using
    call mpp_get_compute_domain(land_domain,isc,iec,jsc,jec)

    im = (iec-isc+1)*(jec-jsc+1)

    ! Create blocks, but again, not currently using
    call define_blocks_packed('land_model', Lnd_block, isc, iec, jsc, jec, 1, &
         ctrl_init%blocksize, block_message)

    lm4_model%control%isc = isc
    lm4_model%control%iec = iec
    lm4_model%control%jsc = jsc
    lm4_model%control%jec = jec
    lm4_model%static%im   = im
    call lm4_model%Create(im)

    ! ! Restart read of sfc_data
    ! call sfc_prop_restart_read(lm4_model, land_domain, .false.)
    ! ! Transfer from sfcprop to model data
    ! call sfc_prop_transfer(lm4_model)


  end subroutine init_driver

  ! ---------------------------------------
  subroutine run_driver(lm4_model)

    type(lm4_type),        intent(inout) :: lm4_model(:) ! land model's variable type

    ! ! local
    real                   :: dt   ! Timestep
    type  (land_data_type) :: Land ! GFDL model dt
    !real(kind_phys)         :: foodata(noah_model%static%im)
    ! !

    !associate(foodata => noah_model%model%foo_atm2lndfield  &
    !     )
    !!
    
    call sfc_boundary_layer(dt,Land)
    call flux_down_from_atmos(Land )

    ! ! Actually run land model
    ! call update_land_model_fast ( Atmos_land_boundary, Land )

    !end associate
  end subroutine run_driver

  ! ---------------------------------------
  subroutine sfc_boundary_layer( dt,Land )
    real,                  intent(in)     :: dt   !< Time step
    type(land_data_type),  intent(inout)  :: Land !< A derived data type to specify land boundary data

    !! blocking not used for now
    integer :: nblocks = 1
    integer :: my_nblocks = 1

    integer, allocatable :: block_start(:), block_end(:)

    real    :: zrefm, zrefh
    
    real, dimension(im) :: &
         ex_albedo,     &
         ex_albedo_vis_dir,     &
         ex_albedo_nir_dir,     &
         ex_albedo_vis_dif,     &
         ex_albedo_nir_dif,     &
         ex_land_frac,  &
         ex_t_atm,      &
         ex_p_atm,      &
         ex_u_atm, ex_v_atm,    &
         ex_gust,       &
         ex_t_surf4,    &
         ex_u_surf, ex_v_surf,  &
         ex_rough_mom, ex_rough_heat, ex_rough_moist, &
         ex_rough_scale,&
         ex_q_star,     &
         ex_cd_q,       &
         ex_ref, ex_ref_u, ex_ref_v, ex_u10, &
         ex_ref2,       &
         ex_t_ref,      &
         ex_qs_ref,     &
         ex_qs_ref_cmip,     &
         ex_del_m,      &
         ex_del_h,      &
         ex_del_q,      &
         ex_frac_open_sea

    
    real, dimension(im,1) :: &
         ex_tr_atm,  &
         ex_tr_surf,    & !< near-surface tracer fields
         ex_flux_tr,    & !< tracer fluxes
         ex_dfdtr_surf, & !< d(tracer flux)/d(surf tracer)
         ex_dfdtr_atm,  & !< d(tracer flux)/d(atm tracer)
         ex_e_tr_n,     & !< coefficient in implicit scheme
         ex_f_tr_delt_n   !< coefficient in implicit scheme
    
    !! -- These were orignally allocatable:
    !!

    logical, dimension(im) :: &
         ex_avail,     &   !< true where data on exchange grid are available
         ex_land           !< true if exchange grid cell is over land

    real, dimension(im) :: &
         ex_t_surf   ,  &
         ex_t_surf_miz, &
         ex_p_surf   ,  &
         ex_slp      ,  &
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
         ex_flux_u,    &
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
         ex_z_atm    ,  &

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
    integer :: is,ie,l,j
    integer :: isc,iec,jsc,jec

    integer :: isphum = 1       !< index of specific humidity tracer in tracer table 

    ! IS this smart?
    is = 1
    ie = im
    allocate(block_start(nblocks), block_end(nblocks))
    
    block_start = is
    block_end   = im
    
    call surface_flux_1d (&
         ex_t_atm(is:ie), ex_tr_atm(is:ie,isphum),  ex_u_atm(is:ie), ex_v_atm(is:ie),  ex_p_atm(is:ie),  ex_z_atm(is:ie),  &
         ex_p_surf(is:ie),ex_t_surf(is:ie), ex_t_ca(is:ie),  ex_tr_surf(is:ie,isphum),                       &
         ex_u_surf(is:ie), ex_v_surf(is:ie),                                           &
         ex_rough_mom(is:ie), ex_rough_heat(is:ie), ex_rough_moist(is:ie), ex_rough_scale(is:ie),    &
         ex_gust(is:ie),                                                        &
         ex_flux_t(is:ie), ex_flux_tr(is:ie,isphum), ex_flux_lw(is:ie), ex_flux_u(is:ie), ex_flux_v(is:ie),         &
         ex_cd_m(is:ie),   ex_cd_t(is:ie), ex_cd_q(is:ie),                                    &
         ex_wind(is:ie),   ex_u_star(is:ie), ex_b_star(is:ie), ex_q_star(is:ie),                     &
         ex_dhdt_surf(is:ie), ex_dedt_surf(is:ie), ex_dfdtr_surf(is:ie,isphum),  ex_drdt_surf(is:ie),        &
         ex_dhdt_atm(is:ie),  ex_dfdtr_atm(is:ie,isphum),  ex_dtaudu_atm(is:ie), ex_dtaudv_atm(is:ie),       &
         dt,                                                             &
         ex_land(is:ie), ex_seawater(is:ie) .gt. 0.0,  ex_avail(is:ie)            )


    !! ....
    zrefm = 10.0
    zrefh = z_ref_heat

    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       call mo_profile ( zrefm, zrefh, ex_z_atm(is:ie), ex_rough_mom(is:ie), &
            ex_rough_heat(is:ie), ex_rough_moist(is:ie),          &
            ex_u_star(is:ie), ex_b_star(is:ie), ex_q_star(is:ie),        &
            ex_del_m(is:ie), ex_del_h(is:ie), ex_del_q(is:ie), ex_avail(is:ie)  )
       do i = is,ie
          ex_u10(i) = 0.
          if(ex_avail(i)) then
             ex_ref_u(i) = ex_u_surf(i) + (ex_u_atm(i)-ex_u_surf(i)) * ex_del_m(i)
             ex_ref_v(i) = ex_v_surf(i) + (ex_v_atm(i)-ex_v_surf(i)) * ex_del_m(i)
             ex_u10(i) = sqrt(ex_ref_u(i)**2 + ex_ref_v(i)**2)
          endif
       enddo
    enddo ! end of block loop

    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          if(ex_avail(i)) ex_drag_q(i) = ex_wind(i)*ex_cd_q(i)
          ! [6] get mean quantities on atmosphere grid
          ! [6.1] compute t surf for radiation
          ex_t_surf4(i) = ex_t_surf(i) ** 4
       enddo
    enddo

    ! [6.3] save atmos albedo fix and old albedo (for downward SW flux calculations)
    ! on exchange grid
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_albedo_fix(i) = 0.
          ex_albedo_vis_dir_fix(i) = 0.
          ex_albedo_nir_dir_fix(i) = 0.
          ex_albedo_vis_dif_fix(i) = 0.
          ex_albedo_nir_dif_fix(i) = 0.
       enddo
    enddo

    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_albedo_fix(i) = (1.0-ex_albedo(i)) / (1.0-ex_albedo_fix(i))
          ex_albedo_vis_dir_fix(i) = (1.0-ex_albedo_vis_dir(i)) / (1.0-ex_albedo_vis_dir_fix(i))
          ex_albedo_nir_dir_fix(i) = (1.0-ex_albedo_nir_dir(i)) / (1.0-ex_albedo_nir_dir_fix(i))
          ex_albedo_vis_dif_fix(i) = (1.0-ex_albedo_vis_dif(i)) / (1.0-ex_albedo_vis_dif_fix(i))
          ex_albedo_nir_dif_fix(i) = (1.0-ex_albedo_nir_dif(i)) / (1.0-ex_albedo_nir_dif_fix(i))
       enddo
    enddo

    !=======================================================================
    ! [7] diagnostics section

    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       call mo_profile ( zrefm, zrefh, ex_z_atm(is:ie),   ex_rough_mom(is:ie), &
            ex_rough_heat(is:ie), ex_rough_moist(is:ie),          &
            ex_u_star(is:ie), ex_b_star(is:ie), ex_q_star(is:ie),        &
            ex_del_m(is:ie), ex_del_h(is:ie), ex_del_q(is:ie), ex_avail(is:ie)  )

       !    ------- reference relative humidity -----------
       !cjg     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
       !cjg          id_rh_ref_cmip > 0 .or. &
       !cjg          id_q_ref > 0 .or. id_q_ref_land >0 ) then
       do i = is,ie
          ex_ref(i) = 1.0e-06
          if (ex_avail(i)) &
               ex_ref(i)   = ex_tr_surf(i,isphum) + (ex_tr_atm(i,isphum)-ex_tr_surf(i,isphum)) * ex_del_q(i)
       enddo
    enddo

    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is,ie
          ex_t_ref(i) = 200.
          if(ex_avail(i)) &
               ex_t_ref(i) = ex_t_ca(i) + (ex_t_atm(i)-ex_t_ca(i)) * ex_del_h(i)
       enddo
       call compute_qs (ex_t_ref(is:ie), ex_p_surf(is:ie), ex_qs_ref(is:ie), q = ex_ref(is:ie))
       call compute_qs (ex_t_ref(is:ie), ex_p_surf(is:ie), ex_qs_ref_cmip(is:ie),  &
            q = ex_ref(is:ie), es_over_liq_and_ice = .true.)
       do i = is,ie
          if(ex_avail(i)) then
             ! remove cap on relative humidity -- this mod requested by cjg, ljd
             !RSH    ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)
             ex_ref2(i)   = 100.*ex_ref(i)/ex_qs_ref_cmip(i)
             ex_ref(i)    = 100.*ex_ref(i)/ex_qs_ref(i)
          endif
       enddo
    enddo

    ! lots of send_data stuff originally here, removed

  end subroutine sfc_boundary_layer

  ! adapted land-only surface_flux_1d from FMS couple
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


    !use constants_mod, only : rdgas, rvgas
    ! not GFS's RDGAS = 287.05
    real,         public, parameter :: RDGAS  = 287.04_r8_kind           !< Gas constant for dry air [J/kg/deg]
    real,         public, parameter :: RVGAS  = 461.50_r8_kind           !< Gas constant for water vapor [J/kg/deg]
    real,         public, parameter :: KAPPA  = 2.0_r8_kind/7.0_r8_kind  !< RDGAS / CP_AIR [dimensionless]
    real,         public, parameter :: CP_AIR = RDGAS/KAPPA              !< Specific heat capacity of dry air at constant pressure [J/kg/deg]

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

    type(land_data_type),  intent(inout)  :: Land !< A derived data type to specify land boundary data 

    if (scale_precip_2d) then
       call mpp_get_compute_domain(Atm%Domain, is_atm, ie_atm, js_atm, je_atm)
       call data_override ('ATM', 'precip_scale2d',    frac_precip,   Time)
       do j=js_atm,je_atm
          do i=is_atm, ie_atm
             Atm%lprec(i,j) = Atm%lprec(i,j)*frac_precip(i,j)
          enddo
       enddo
    endif

    if (partition_fprec_from_lprec .and. Atm%pe) then
       call mpp_get_compute_domain(Atm%Domain, is_atm, ie_atm, js_atm, je_atm)
       do j=js_atm,je_atm
          do i=is_atm, ie_atm
             if (Atm%t_bot(i,j) < tfreeze) then
                Atm%fprec(i,j) = Atm%lprec(i,j)
                Atm%lprec(i,j) = 0.0
             endif
          enddo
       enddo
    endif

    ! MOD update stresses using atmos delta's but derivatives on exchange grid
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_flux_u,ex_delta_u, &
    !$OMP                                  ex_dtaudu_atm,ex_dtaudv_atm,ex_flux_v,ex_delta_v )     &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_flux_u(i) = ex_flux_u(i) + ex_delta_u(i)*ex_dtaudu_atm(i)
          ex_flux_v(i) = ex_flux_v(i) + ex_delta_v(i)*ex_dtaudv_atm(i)
       enddo
    enddo

    !---- adjust sw flux for albedo variations on exch grid ----
    !---- adjust 4 categories (vis/nir dir/dif) separately  ----
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_flux_sw_dir(i) = ex_flux_sw_dir(i) - ex_flux_sw_vis_dir(i)     ! temporarily nir/dir
          ex_flux_sw_dir(i) = ex_flux_sw_dir(i) * ex_albedo_nir_dir_fix(i)  ! fix nir/dir
          ex_flux_sw_vis_dir(i) = ex_flux_sw_vis_dir(i) * ex_albedo_vis_dir_fix(i) ! fix vis/dir
          ex_flux_sw_dir(i) = ex_flux_sw_dir(i) + ex_flux_sw_vis_dir(i)     ! back to total dir

          ex_flux_sw_dif(i) = ex_flux_sw_dif(i) - ex_flux_sw_vis_dif(i)     ! temporarily nir/dif
          ex_flux_sw_dif(i) = ex_flux_sw_dif(i) * ex_albedo_nir_dif_fix(i)  ! fix nir/dif
          ex_flux_sw_vis_dif(i) = ex_flux_sw_vis_dif(i) * ex_albedo_vis_dif_fix(i) ! fix vis/dif
          ex_flux_sw_dif(i) = ex_flux_sw_dif(i) + ex_flux_sw_vis_dif(i)     ! back to total dif

          ex_flux_sw_vis(i) = ex_flux_sw_vis_dir(i) + ex_flux_sw_vis_dif(i) ! legacy, remove later
          ex_flux_sw(i)     = ex_flux_sw_dir(i)     + ex_flux_sw_dif(i)     ! legacy, remove later
       enddo
    enddo

    !----- adjust fluxes for implicit dependence on atmosphere ----
    ! is this needed to copy over?


    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          !----- compute net longwave flux (down-up) -----
          ! (note: lw up already in ex_flux_lw)
          ex_flux_lw(i) = ex_flux_lwd(i) - ex_flux_lw(i)
          if (ex_avail(i) ) then

             ! temperature
             ex_gamma(i)      =  1./ (1.0 - ex_dtmass(i)*(ex_dflux_t(i) + ex_dhdt_atm(i)*cp_inv))
             ex_e_t_n(i)      =  ex_dtmass(i)*ex_dhdt_surf(i)*cp_inv*ex_gamma(i)
             ex_f_t_delt_n(i) = (ex_delta_t(i) + ex_dtmass(i) * ex_flux_t(i)*cp_inv) * ex_gamma(i)

             ex_flux_t (i)    =  ex_flux_t(i)        + ex_dhdt_atm(i) * ex_f_t_delt_n(i)
             ex_dhdt_surf(i)  =  ex_dhdt_surf(i)     + ex_dhdt_atm(i) * ex_e_t_n(i)
             ! moisture
             !     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_q + ex_dedq_atm))
             ! here it looks like two derivatives with different units are added together,
             ! but in fact they are not: ex_dedt_surf and ex_dedq_surf defined in complimentary
             ! regions of exchange grid, so that if one of them is not zero the other is, and
             ! vice versa.
             !     ex_e_q_n      =  ex_dtmass*(ex_dedt_surf+ex_dedq_surf) * ex_gamma
             !     ex_f_q_delt_n = (ex_delta_q  + ex_dtmass * ex_flux_q) * ex_gamma
             !     ex_flux_q     =  ex_flux_q    + ex_dedq_atm * ex_f_q_delt_n
             !     ex_dedt_surf  =  ex_dedt_surf + ex_dedq_atm * ex_e_q_n
             !     ex_dedq_surf  =  ex_dedq_surf + ex_dedq_atm * ex_e_q_n
             ! moisture vs. surface temperture, assuming saturation
             ex_gamma(i)   =  1.0 / (1.0 - ex_dtmass(i)*(ex_dflux_tr(i,isphum) + ex_dfdtr_atm(i,isphum)))
             ex_e_q_n(i)      =  ex_dtmass(i) * ex_dedt_surf(i) * ex_gamma(i)
             ex_dedt_surf(i)  =  ex_dedt_surf(i) + ex_dfdtr_atm(i,isphum) * ex_e_q_n(i)
             do tr = 1,n_exch_tr
                ex_gamma(i)   =  1.0 / (1.0 - ex_dtmass(i)*(ex_dflux_tr(i,tr) + ex_dfdtr_atm(i,tr)))

                ex_e_tr_n(i,tr)      =  ex_dtmass(i)*ex_dfdtr_surf(i,tr)*ex_gamma(i)
                ex_f_tr_delt_n(i,tr) = (ex_delta_tr(i,tr)+ex_dtmass(i)*ex_flux_tr(i,tr))*ex_gamma(i)

                ex_flux_tr(i,tr)     =  ex_flux_tr(i,tr) + ex_dfdtr_atm(i,tr)*ex_f_tr_delt_n(i,tr)
                ex_dfdtr_surf(i,tr)  =  ex_dfdtr_surf(i,tr) + ex_dfdtr_atm(i,tr)*ex_e_tr_n(i,tr)
             enddo
          endif
       enddo ! i = is, ie
    enddo !  l = 1, my_nblocks


  end subroutine flux_down_from_atmos

  ! ----------------------------------------

  subroutine  flux_up_to_atmos( Land )

    type(land_data_type),  intent(in)    :: Land !< A derived data type to specify land boundary data
    
    where (Land%mask(:,:,1))
       t_surf_new = Land%t_surf(:,:,1)
       t_ca_new   = Land%t_ca  (:,:,1)
    endwhere

    !??????? should this be done in land model ??????
    call escomp (t_surf_new, q_surf_new)
    where (Land%mask(:,:,1))
       q_surf_new = Land%tr(:,:,1,1)
    elsewhere
       !q_surf_new = d622*q_surf_new/(p_surf-d378*q_surf_new)
    endwhere

    dt_t_ca   = t_ca_new   - t_ca   ! changes in near-surface T
    dt_t_surf = t_surf_new - t_surf ! changes in radiative T
    dt_q_surf = q_surf_new - q_surf ! changes in near-surface q

    ! adjust fluxes and atmospheric increments for
    ! implicit dependence on surface temperature

    flux_t        = flux_t      + dt_t_ca  *dhdt_surf
    flux_lw       = flux_lw     - dt_t_surf*drdt_surf
    Boundary%dt_t = f_t_delt_n  + dt_t_ca  *e_t_n

    where (Land%mask(:,:,1))
       flux_q                     = flux_q      + dt_q_surf*dedq_surf
       Boundary%dt_tr(:,:,isphum) = f_q_delt_n  + dt_q_surf*e_q_n
    elsewhere
       !flux_q                     = flux_q      + dt_t_surf*dedt_surf
       !Boundary%dt_tr(:,:,isphum) = f_q_delt_n  + dt_t_surf*e_q_n
    endwhere

  end subroutine flux_up_to_atmos

end module lm4_driver
