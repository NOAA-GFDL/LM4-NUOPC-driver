!! This is a port of some of FMS Coupler's surface_flux_mod.F90 to LM4-UFS
!! ============================================================================
module lm4_surface_flux_mod

   use             fms_mod, only: close_file, mpp_pe, mpp_root_pe, write_version_number
   use             fms_mod, only: file_exist, check_nml_error, open_namelist_file, stdlog
   use   monin_obukhov_mod, only: mo_drag, mo_profile, monin_obukhov_init
   use  sat_vapor_pres_mod, only: escomp, descomp
   use       constants_mod, only: cp_air, hlv, stefan, rdgas, rvgas, grav, vonkarm
   use             mpp_mod, only: input_nml_file, FATAL, mpp_error

   implicit none
   private

! ==== public interface ======================================================
   public  lm4_surface_flux_1d, lm4_surface_flux_init
! ==== end of public interface ===============================================


   logical :: module_is_initialized = .false.

   real, parameter :: d622   = rdgas/rvgas
   real, parameter :: d378   = 1.-d622
   real, parameter :: hlars  = hlv/rvgas
   real, parameter :: gcp    = grav/cp_air
   real, parameter :: kappa  = rdgas/cp_air
   real            :: d608   = d378/d622
   ! d608 set to zero at initialization if the use of
   ! virtual temperatures is turned off in namelist


! ---- namelist with default values ------------------------------------------
   logical :: no_neg_q              = .false. !< If a_atm_in (specific humidity) is negative (because of numerical truncation),
   !! then override with 0.0
   logical :: use_virtual_temp      = .true.  !< If .TRUE., use virtual potential temp to calculate the stability of the surface
   !! layer.  If .FALSE., use potential temp.
   logical :: alt_gustiness         = .false. !< An alternaive formulation for gustiness calculation.  A minimum bound on the wind
   !! speed used influx calculations, with the bound equal to gust_const
   logical :: old_dtaudv            = .false. !< The derivative of surface wind stress with respect to the zonal wind and meridional
   !! wind are approximated by the same tendency
   logical :: use_mixing_ratio      = .false. !< An option to provide capability to run the Manabe Climate form of the surface flux
   !! (coded for legacy purposes).
   real    :: gust_const            =  1.0 !< Constant for alternative gustiness calculation
   real    :: gust_min              =  0.0 !< Minimum gustiness used when alt_gustiness is .FALSE.
   logical :: ncar_ocean_flux       = .false. !< Use NCAR climate model turbulent flux calculation described by Large and Yeager,
   !! NCAR Technical Document, 2004
   logical :: ncar_ocean_flux_orig  = .false. !< Use NCAR climate model turbulent flux calculation described by Large and Yeager,
   !! NCAR Technical Document, 2004, using the original GFDL implementation, which
   !! contains a bug in the specification of the exchange coefficient for the sensible
   !! heat.  This option is available for legacy purposes, and is not recommended for
   !! new experiments.
   logical :: raoult_sat_vap        = .false. !< Reduce saturation vapor pressure to account for seawater
   logical :: do_simple             = .false.


   namelist /surface_flux_nml/ no_neg_q,             &
      use_virtual_temp,     &
      alt_gustiness,        &
      gust_const,           &
      gust_min,             &
      old_dtaudv,           &
      use_mixing_ratio,     &
      ncar_ocean_flux,      &
      ncar_ocean_flux_orig, &
      raoult_sat_vap,       &
      do_simple



contains


! ============================================================================
! Initialization of the surface flux module--reads the nml.
   subroutine lm4_surface_flux_init

      ! ---- local vars ----------------------------------------------------------
      integer :: unit, ierr, io

      ! read namelist
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, surface_flux_nml, iostat=io)
      ierr = check_nml_error(io,'surface_flux_nml')
#else
      if ( file_exist('input.nml')) then
         unit = open_namelist_file ()
         ierr=1;
         do while (ierr /= 0)
            read  (unit, nml=surface_flux_nml, iostat=io)
            ierr = check_nml_error(io,'surface_flux_nml')
         enddo
           call close_file (unit)
      endif
#endif

      unit = stdlog()
      if ( mpp_pe() == mpp_root_pe() )  write (unit, nml=surface_flux_nml)

      if(.not. use_virtual_temp) d608 = 0.0

      call monin_obukhov_init()

      module_is_initialized = .true.

   end subroutine lm4_surface_flux_init



! ============================================================================
   subroutine lm4_surface_flux_1d (                                     &
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
         rho_drag, drag_t,    drag_m,   drag_q,    rho,      &
         q_atm,    q_surf0,  dw_atmdu,  dw_atmdv,  w_gust

      integer :: i, nbad


      if (.not. module_is_initialized) &
         call mpp_error(FATAL, "surface_flux_1d: surface_flux_init is not called")

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

      if (raoult_sat_vap) where (seawater) q_surf0 = 0.98 * q_surf0

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

      !! DEBUG TMP check inputs to mo_drag
      write(*,*) 'thv_atm min/max:' , minval(thv_atm), maxval(thv_atm)
      write(*,*) 'thv_surf min/max:', minval(thv_surf), maxval(thv_surf)
      write(*,*) 'z_atm min/max:'   , minval(z_atm), maxval(z_atm)
      write(*,*) 'rough_mom min/max:', minval(rough_mom), maxval(rough_mom)
      write(*,*) 'rough_heat min/max:', minval(rough_heat), maxval(rough_heat)
      write(*,*) 'rough_moist min/max:', minval(rough_moist), maxval(rough_moist)
      write(*,*) 'w_atm min/max:'   , minval(w_atm), maxval(w_atm)

      if (any(isnan(thv_atm))) then
         write(*,*) 'thv_atm contains NaNs'
      end if
      
      if (any(isnan(thv_surf))) then
         write(*,*) 'thv_surf contains NaNs'
      end if
      
      if (any(isnan(z_atm))) then
         write(*,*) 'z_atm contains NaNs'
      end if
      
      if (any(isnan(rough_mom))) then
         write(*,*) 'rough_mom contains NaNs'
      end if
      
      if (any(isnan(rough_heat))) then
         write(*,*) 'rough_heat contains NaNs'
      end if
      
      if (any(isnan(rough_moist))) then
         write(*,*) 'rough_moist contains NaNs'
      end if
      
      if (any(isnan(w_atm))) then
         write(*,*) 'w_atm contains NaNs'
      end if

      ! Check that thv_atm and thv_surf are positive
      if (any(thv_atm <= 200.0)) then
         write(*,*) 'thv_atm is below 200 K'
      end if
      if (any(thv_atm > 330.0)) then
         write(*,*) 'thv_atm is above 330 K'
      end if
      
      if (any(thv_surf <= 200.0)) then
         write(*,*) 'thv_surf is below 200 K'
      end if
      if (any(thv_surf > 330.0)) then
         write(*,*) 'thv_surf is above 330 K'
      end if
      
      ! Check that z_atm is positive
      if (any(z_atm <= 0.0)) then
         write(*,*) 'z_atm contains non-positive values'
      end if
      
      ! Check that rough_mom, rough_heat, and rough_moist are non-negative
      if (any(rough_mom < 0.0)) then
         write(*,*) 'rough_mom contains negative values'
      end if
      
      if (any(rough_heat < 0.0)) then
         write(*,*) 'rough_heat contains negative values'
      end if
      
      if (any(rough_moist < 0.0)) then
         write(*,*) 'rough_moist contains negative values'
      end if
      
      ! Check that w_atm is within a reasonable range
      if (any(abs(w_atm) > 50.0)) then
         write(*,*) 'w_atm is above 50 m/s'
      end if
      
      !! END DEBUG TMP


      !  monin-obukhov similarity theory
      call mo_drag (thv_atm, thv_surf, z_atm,                &
         rough_mom, rough_heat, rough_moist, w_atm,          &
         cd_m, cd_t, cd_q, u_star, b_star, avail             )

    !   ! override with ocean fluxes from NCAR calculation
    !   if (ncar_ocean_flux .or. ncar_ocean_flux_orig) then
    !      call  ncar_ocean_fluxes (w_atm, th_atm, t_surf0, q_atm, q_surf0, z_atm, &
    !         seawater, cd_m, cd_t, cd_q, u_star, b_star     )
    !   end if

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

   end subroutine lm4_surface_flux_1d


end module lm4_surface_flux_mod
