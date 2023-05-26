module lm4_type_mod

   use land_data_mod,    only: land_data_type, atmos_land_boundary_type
   use time_manager_mod, only: time_type

   public

   type :: lm4_nml_type
      integer           :: lm4_debug    ! debug flag for lm4 (0=off, 1=low, 2=high)
      ! grid, domain, and blocking
      integer           :: npx, npy     
      integer           :: ntiles
      integer           :: layout(2)
      character(len=64) :: grid
      integer           :: blocksize
   end type lm4_nml_type

   type :: lm4_control_type
      logical   :: first_time  ! flag for first time step
      integer   :: mype
      integer   :: nblks, blksz, isc, iec, jsc, jec
   end type lm4_control_type

   ! type for atmospheric forcing data, based off atmos_solo_land's atmos_data_type
   type, public :: atm_forc_type
      real, pointer, dimension(:) ::  &
         t_bot       => NULL(), &   ! temperature at the atm. bottom, degK
         z_bot       => NULL(), &   ! altitude of the atm. bottom above sfc., m
         p_bot       => NULL(), &   ! pressure at the atm. bottom, N/m2
         u_bot       => NULL(), &   ! zonal wind at the atm. bottom, m/s
         v_bot       => NULL(), &   ! merid. wind at the atm. bottom, m/s
         q_bot       => NULL(), &   ! specific humidity at the atm. bottom, kg/kg
         p_surf      => NULL(), &   ! surface pressure, N/m2
         ! LM4 doesn't need slp
         !slp       => NULL(), &   ! sea level pressure, N/m2 
         gust        => NULL(), &   ! gustiness, m/s
         coszen      => NULL(), &   ! cosine of zenith angle
         totprec     => NULL(), &   ! total precipitation,  kg/m2/s
         lprec       => NULL(), &   ! liquid precipitation, kg/m2/s
         fprec       => NULL(), &   ! frozen precipitation, kg/m2/s          
         ! LM4 doesn't need net SW
         !flux_sw   => NULL(), &   ! SW radiation flux (net), W/m2
         flux_lw     => NULL(), &   ! LW radiation flux (down), W/m2
         flux_sw_dir => NULL(), &
         flux_sw_dif => NULL(), &        
         flux_sw_down_vis_dir   => NULL(), &
         flux_sw_down_vis_dif   => NULL(), &
         flux_sw_down_total_dir => NULL(), &
         flux_sw_down_total_dif => NULL(), &
         flux_sw_vis            => NULL(), &
         flux_sw_vis_dir        => NULL(), &
         flux_sw_vis_dif        => NULL()

   end type atm_forc_type

   ! TMP DEBUG
   type, public :: atm_forc2d_type
      real, pointer, dimension(:,:) ::  &
         t_bot       => NULL(), &   ! temperature at the atm. bottom, degK
         z_bot       => NULL(), &   ! altitude of the atm. bottom above sfc., m
         p_bot       => NULL(), &   ! pressure at the atm. bottom, N/m2
         u_bot       => NULL(), &   ! zonal wind at the atm. bottom, m/s
         v_bot       => NULL(), &   ! merid. wind at the atm. bottom, m/s
         q_bot       => NULL(), &   ! specific humidity at the atm. bottom, kg/kg
         p_surf      => NULL(), &   ! surface pressure, N/m2
      ! LM4 doesn't need slp
      !slp       => NULL(), &   ! sea level pressure, N/m2
         totprec     => NULL(), &   ! total precipitation,  kg/m2/s
         lprec       => NULL(), &   ! liquid precipitation, kg/m2/s
         fprec       => NULL(), &   ! frozen precipitation, kg/m2/s          
         gust        => NULL(), &   ! gustiness, m/s
         coszen      => NULL(), &   ! cosine of zenith angle
      ! LM4 doesn't need
      !flux_sw   => NULL(), &   ! SW radiation flux (net), W/m2
         flux_lw   => NULL(), &   ! LW radiation flux (down), W/m2
         flux_sw_dir => NULL(), &
         flux_sw_dif => NULL(), &
         flux_sw_down_vis_dir   => NULL(), &
         flux_sw_down_vis_dif   => NULL(), &
         flux_sw_down_total_dir => NULL(), &
         flux_sw_down_total_dif => NULL(), &
         flux_sw_vis            => NULL(), &
         flux_sw_vis_dir        => NULL(), &
         flux_sw_vis_dif        => NULL()
   end type atm_forc2d_type
   ! END TMP DEBUG




   type :: lm4_type
      type(lm4_nml_type)             :: nml        ! namelist
      type(lm4_control_type)         :: control
      type(atm_forc_type)            :: atm_forc   ! data from atm 
      type(atm_forc2d_type)          :: atm_forc2d ! TMP DEBUG
      ! these are passed to the land model's routines:
      type(land_data_type)           :: From_lnd   ! data from land
      type(atmos_land_boundary_type) :: From_atm   ! data from atm      
      type(time_type)                :: Time_land, Time_init, Time_end,  &
                                        Time_step_land, Time_step_ocean, &
                                        Time_restart, Time_step_restart, &
                                        Time_atstart      
      !contains
      !   procedure, public  :: Create

   end type lm4_type

contains

   subroutine dealloc_atmforc(bnd)
      !! for every variable in atm_forc_type, if associated, deallocate
      
      type(atm_forc_type), intent(inout) :: bnd

      if (associated(bnd%t_bot)) deallocate(bnd%t_bot)
      if (associated(bnd%z_bot)) deallocate(bnd%z_bot)
      if (associated(bnd%p_bot)) deallocate(bnd%p_bot)
      if (associated(bnd%u_bot)) deallocate(bnd%u_bot)
      if (associated(bnd%v_bot)) deallocate(bnd%v_bot)
      if (associated(bnd%q_bot)) deallocate(bnd%q_bot)
      if (associated(bnd%p_surf)) deallocate(bnd%p_surf)
      !if (associated(bnd%slp)) deallocate(bnd%slp)
      if (associated(bnd%gust)) deallocate(bnd%gust)
      if (associated(bnd%coszen)) deallocate(bnd%coszen)
      !if (associated(bnd%flux_sw)) deallocate(bnd%flux_sw)
      if (associated(bnd%flux_lw)) deallocate(bnd%flux_lw)
      if (associated(bnd%flux_sw_dir)) deallocate(bnd%flux_sw_dir)
      if (associated(bnd%flux_sw_dif)) deallocate(bnd%flux_sw_dif)
      if (associated(bnd%flux_sw_down_vis_dir)) deallocate(bnd%flux_sw_down_vis_dir)
      if (associated(bnd%flux_sw_down_vis_dif)) deallocate(bnd%flux_sw_down_vis_dif)
      if (associated(bnd%flux_sw_down_total_dir)) deallocate(bnd%flux_sw_down_total_dir)
      if (associated(bnd%flux_sw_down_total_dif)) deallocate(bnd%flux_sw_down_total_dif)
      if (associated(bnd%flux_sw_vis)) deallocate(bnd%flux_sw_vis)
      if (associated(bnd%flux_sw_vis_dir)) deallocate(bnd%flux_sw_vis_dir)
      if (associated(bnd%flux_sw_vis_dif)) deallocate(bnd%flux_sw_vis_dif)
      if (associated(bnd%totprec)) deallocate(bnd%totprec)
      if (associated(bnd%lprec)) deallocate(bnd%lprec)
      if (associated(bnd%fprec)) deallocate(bnd%fprec)

   end subroutine dealloc_atmforc

   subroutine alloc_atmforc(bnd)
      ! must be called after land_data_init

      use land_data_mod, only : lnd

      type(atm_forc_type), intent(inout) :: bnd

      call dealloc_atmforc(bnd)

      allocate( bnd%t_bot(lnd%ls:lnd%le) )
      allocate( bnd%z_bot(lnd%ls:lnd%le) )
      allocate( bnd%p_bot(lnd%ls:lnd%le) )
      allocate( bnd%u_bot(lnd%ls:lnd%le) )
      allocate( bnd%v_bot(lnd%ls:lnd%le) )
      allocate( bnd%q_bot(lnd%ls:lnd%le) )
      allocate( bnd%p_surf(lnd%ls:lnd%le) )
      !allocate( bnd%slp(lnd%ls:lnd%le) )
      allocate( bnd%gust(lnd%ls:lnd%le) )
      allocate( bnd%coszen(lnd%ls:lnd%le) )
      !allocate( bnd%flux_sw(lnd%ls:lnd%le) )
      allocate( bnd%flux_lw(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_dir(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_dif(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_down_vis_dir(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_down_vis_dif(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_down_total_dir(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_down_total_dif(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_vis(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_vis_dir(lnd%ls:lnd%le) )
      allocate( bnd%flux_sw_vis_dif(lnd%ls:lnd%le) )
      allocate( bnd%totprec(lnd%ls:lnd%le) )
      allocate( bnd%lprec(lnd%ls:lnd%le) )
      allocate( bnd%fprec(lnd%ls:lnd%le) )

   end subroutine alloc_atmforc

   ! TMP DEBUG
   subroutine dealloc_atmforc2d(bnd)
      use land_data_mod, only : lnd

      type(atm_forc2d_type), intent(inout) :: bnd
      ! for every variable in atm_forc_type, if associated, deallocate
      if (associated(bnd%t_bot)) deallocate(bnd%t_bot)
      if (associated(bnd%z_bot)) deallocate(bnd%z_bot)
      if (associated(bnd%p_bot)) deallocate(bnd%p_bot)
      if (associated(bnd%u_bot)) deallocate(bnd%u_bot)
      if (associated(bnd%v_bot)) deallocate(bnd%v_bot)
      if (associated(bnd%q_bot)) deallocate(bnd%q_bot)
      if (associated(bnd%p_surf)) deallocate(bnd%p_surf)
      !if (associated(bnd%slp)) deallocate(bnd%slp)
      if (associated(bnd%gust)) deallocate(bnd%gust)
      if (associated(bnd%coszen)) deallocate(bnd%coszen)
      !if (associated(bnd%flux_sw)) deallocate(bnd%flux_sw)
      if (associated(bnd%flux_lw)) deallocate(bnd%flux_lw)
      if (associated(bnd%flux_sw_dir)) deallocate(bnd%flux_sw_dir)
      if (associated(bnd%flux_sw_dif)) deallocate(bnd%flux_sw_dif)
      if (associated(bnd%flux_sw_down_vis_dir)) deallocate(bnd%flux_sw_down_vis_dir)
      if (associated(bnd%flux_sw_down_vis_dif)) deallocate(bnd%flux_sw_down_vis_dif)
      if (associated(bnd%flux_sw_down_total_dir)) deallocate(bnd%flux_sw_down_total_dir)
      if (associated(bnd%flux_sw_down_total_dif)) deallocate(bnd%flux_sw_down_total_dif)
      if (associated(bnd%flux_sw_vis)) deallocate(bnd%flux_sw_vis)
      if (associated(bnd%flux_sw_vis_dir)) deallocate(bnd%flux_sw_vis_dir)
      if (associated(bnd%flux_sw_vis_dif)) deallocate(bnd%flux_sw_vis_dif)
      if (associated(bnd%totprec)) deallocate(bnd%totprec)
      if (associated(bnd%lprec)) deallocate(bnd%lprec)
      if (associated(bnd%fprec)) deallocate(bnd%fprec)
   end subroutine dealloc_atmforc2d

   subroutine alloc_atmforc2d(bnd)
      ! must be called after land_data_init
      use land_data_mod, only : lnd

      type(atm_forc2d_type), intent(inout) :: bnd

      call dealloc_atmforc2d(bnd)

      allocate( bnd%t_bot(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%z_bot(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%p_bot(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%u_bot(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%v_bot(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%q_bot(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%p_surf(lnd%is:lnd%ie,lnd%js:lnd%je) )
      !allocate( bnd%slp(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%gust(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%coszen(lnd%is:lnd%ie,lnd%js:lnd%je) )
      !allocate( bnd%flux_sw(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_lw(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_dir(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_dif(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_down_vis_dir(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_down_vis_dif(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_down_total_dir(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_down_total_dif(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_vis(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_vis_dir(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%flux_sw_vis_dif(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%totprec(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%lprec(lnd%is:lnd%ie,lnd%js:lnd%je) )
      allocate( bnd%fprec(lnd%is:lnd%ie,lnd%js:lnd%je) )
   end subroutine alloc_atmforc2d
	! TMP DEBUG

end module lm4_type_mod
