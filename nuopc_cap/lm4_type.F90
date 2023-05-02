module lm4_type_mod
	! TODO: clean up
	public

	type :: lm4_control_type
		logical   :: first_time  ! flag for first time step
		integer   :: mype
		integer   :: nblks, blksz, isc, iec, jsc, jec
	end type lm4_control_type

	! type for atmospheric forcing data, based off atmos_solo_land's atmos_data_type
	type, public :: atm_forc_type
		real, pointer, dimension(:,:) ::  &
		t_bot     => NULL(), &   ! temperature at the atm. bottom, degK                                                                                                                                              
		z_bot     => NULL(), &   ! altitude of the atm. bottom above sfc., m                                                                                                                                         
		p_bot     => NULL(), &   ! pressure at the atm. bottom, N/m2                                                                                                                                                 
		u_bot     => NULL(), &   ! zonal wind at the atm. bottom, m/s                                                                                                                                                
		v_bot     => NULL(), &   ! merid. wind at the atm. bottom, m/s                                                                                                                                               
		p_surf    => NULL(), &   ! surface pressure, N/m2                                                                                                                                                            
		slp       => NULL(), &   ! sea level pressure, N/m2                                                                                                                                                          
		gust      => NULL(), &   ! gustiness, m/s                                                                                                                                                                    
		coszen    => NULL(), &   ! cosine of zenith angle                                                                                                                                                            
		flux_sw   => NULL(), &   ! SW radiation flux (net), W/m2                                                                                                                                                     
		flux_lw   => NULL(), &   ! LW radiation flux (down), W/m2                                                                                                                                                    
		flux_sw_dir => NULL(),&
		flux_sw_dif => NULL(),&
		flux_sw_down_vis_dir   => NULL(), &
		flux_sw_down_vis_dif   => NULL(), &
		flux_sw_down_total_dir => NULL(), &
		flux_sw_down_total_dif => NULL(), &
		flux_sw_vis            => NULL(), &
		flux_sw_vis_dir        => NULL(), &
		flux_sw_vis_dif        => NULL(), &
		lprec     => NULL(), &   ! liquid precipitation, kg/m2/s                                                                                                                                                     
		fprec     => NULL()      ! frozen precipitation, kg/m2/s  
	end type atm_forc_type
    


  type :: lm4_type
     type(lm4_control_type) :: control
   !contains

   !   procedure, public  :: Create

  end type lm4_type

contains

	subroutine dealloc_atmforc(bnd)

		type(atm_forc_type), intent(inout) :: bnd

		! for every variable in atm_forc_type, if associated, deallocate
		if (associated(bnd%t_bot)) deallocate(bnd%t_bot)
		if (associated(bnd%z_bot)) deallocate(bnd%z_bot)
		if (associated(bnd%p_bot)) deallocate(bnd%p_bot)
		if (associated(bnd%u_bot)) deallocate(bnd%u_bot)
		if (associated(bnd%v_bot)) deallocate(bnd%v_bot)
		if (associated(bnd%p_surf)) deallocate(bnd%p_surf)
		if (associated(bnd%slp)) deallocate(bnd%slp)
		if (associated(bnd%gust)) deallocate(bnd%gust)
		if (associated(bnd%coszen)) deallocate(bnd%coszen)
		if (associated(bnd%flux_sw)) deallocate(bnd%flux_sw)
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
		allocate( bnd%p_surf(lnd%ls:lnd%le) )
		allocate( bnd%slp(lnd%ls:lnd%le) ) 
		allocate( bnd%gust(lnd%ls:lnd%le) )
		allocate( bnd%coszen(lnd%ls:lnd%le) )
		allocate( bnd%flux_sw(lnd%ls:lnd%le) ) 
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
		allocate( bnd%lprec(lnd%ls:lnd%le) )
		allocate( bnd%fprec(lnd%ls:lnd%le) )

	end subroutine alloc_atmforc

end module lm4_type_mod
