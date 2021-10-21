module lm4_driver

  use machine, only: kind_phys

  use proc_bounds,        only: procbounds_type, control_init_type
  use mpp_domains_mod,    only: domain2d


  implicit none
  private

  type(domain2D),           public :: land_domain
  type(control_init_type),  public :: ctrl_init

  public :: init_driver 

contains

  !subroutine init_driver(procbounds)
  subroutine init_driver(ctrl_init)

    use mpp_domains_mod,    only: domain2d, mpp_get_compute_domain
    use mpp_mod,            only: mpp_pe, mpp_root_pe
    use land_domain_mod,    only: domain_create
    use block_control_mod,  only: block_control_type, define_blocks_packed
    use land_restart_mod,   only: sfc_prop_restart_read, sfc_prop_transfer
    type(control_init_type), intent(out)  ::   ctrl_init
    
    ! ---------------
    ! local

    type (block_control_type), target   :: Lnd_block !  Block container
    integer :: isc, iec, jsc, jec
    
    
    call ctrl_init%init()

    if (mpp_pe() == mpp_root_pe()) then
       write(*,*) 'ctrl_init%grid: '     ,ctrl_init%grid
       write(*,*) 'ctrl_init%npx: '      ,ctrl_init%npx
       write(*,*) 'ctrl_init%npy: '      ,ctrl_init%npy
       write(*,*) 'ctrl_init%layout: '   ,ctrl_init%layout
       write(*,*) 'ctrl_init%ntiles: '   ,ctrl_init%ntiles
       write(*,*) 'ctrl_init%blocksize: ',ctrl_init%blocksize
       write(*,*) 'ctrl_init%ivegsrc: '  ,ctrl_init%ivegsrc
       write(*,*) 'ctrl_init%isot: '     ,ctrl_init%isot
    end if

    ! FMS domain creation:
    call domain_create(ctrl_init, land_domain)

    ! Create blocking a la FV3, but not currently using
    call mpp_get_compute_domain(land_domain,isc,iec,jsc,jec)

    im = (iec-isc+1)*(jec-jsc+1)   
    
    ! Create blocks, but again, not currently using
    call define_blocks_packed('land_model', Lnd_block, isc, iec, jsc, jec, 1, &
         ctrl_init%blocksize, block_message)

    ! lm4_model%control%isc = isc
    ! lm4_model%control%iec = iec
    ! lm4_model%control%jsc = jsc
    ! lm4_model%control%jec = jec
    ! lm4_model%static%im  = im
    ! call lm4_model%Create(im)

    ! ! Restart read of sfc_data
    ! call sfc_prop_restart_read(lm4_model, land_domain, .false.)
    ! ! Transfer from sfcprop to model data
    ! call sfc_prop_transfer(lm4_model) 


  end subroutine init_driver

  
end module lm4_driver
