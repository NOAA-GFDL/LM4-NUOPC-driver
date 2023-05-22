!!! TODO: This should probably be cleaned up and/or moved

module proc_bounds

  ! proc and gridcell bounds for simple  mesh with regular decomp
  ! also, control and init type
  
  implicit none

  public procbounds_type, control_init_type
  
  private

  type procbounds_type
     integer :: gridbeg, gridend ! local gridcell range
     integer :: de               ! local de
     integer :: im               ! # gridcells on de
  end type procbounds_type

  type control_init_type
     ! modelled after FV3 GFS init and control types

     ! namelist variables
     ! ------------------------------------------

     integer           :: lm4_debug    ! debug flag for lm4 (0=off, 1=low, 2=high)
     ! grid, domain, and blocking
     integer           :: npx, npy     
     integer           :: ntiles
     integer           :: layout(2)
     character(len=64) :: grid
     integer           :: blocksize

     ! run
     logical   :: first_time  ! flag for first time step

     ! MPI stuff
     integer :: me                                ! my MPI-rank
     integer :: master                            ! master MPI-rank


  contains
    procedure :: init  => control_initialize
  end type control_init_type
  
  type(procbounds_type), public :: procbounds

contains

  subroutine control_initialize(Model)

    use fms_mod,             only: check_nml_error, close_file, file_exist
    use mpp_mod,             only: mpp_pe, mpp_root_pe
#ifdef INTERNAL_FILE_NML
    use mpp_mod,             only: input_nml_file
#else
    use fms_mod,             only: open_namelist_file
#endif

    
    implicit none

    class(control_init_type)            :: Model
    ! namelist variables
    ! ------------------------------------------
    ! grid, domain, and blocking
    integer           :: lm4_debug = 0
    integer           :: npx = 0, npy = 0, ntiles = 0
    integer           :: layout(2) = (/0,0/)
    character(len=64) :: grid      = 'none'
    integer           :: blocksize = -1


    ! for namelist read
    integer :: unit, io, ierr
    namelist /lm4_nml/ grid, npx, npy, layout, ntiles, &
         blocksize, lm4_debug
    
    ! -------------------------------------------
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

    Model%lm4_debug = lm4_debug
    Model%grid      = grid
    Model%blocksize = blocksize
    Model%npx       = npx
    Model%npy       = npy
    Model%layout    = layout
    Model%ntiles    = ntiles
    !--- MPI parameters
    Model%me        =  mpp_pe()
    Model%master    =  mpp_root_pe()
!     Model%fn_nml           = fn_nml


  end subroutine control_initialize

end module proc_bounds
