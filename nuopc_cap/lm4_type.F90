

module lm4_type_mod

  !use machine, only: kind_phys

  implicit none
  save
  private

  !--- originally from kind_phys
  integer, parameter :: kind_phys = 8
  
  !--- unit number
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6
  
  !--- parameter constants used for default initializations
  real(kind_phys), parameter :: zero      = 0.0_kind_phys
  real(kind_phys), parameter :: clear_val = -9999_kind_phys

  type :: lm4_control_type
     logical   :: first_time  ! flag for first time step
     integer   :: mype
     integer   :: nblks, blksz, isc, iec, jsc, jec
  end type lm4_control_type

  type :: lm4_static_type

     integer            ::   ltile      ! subgrid land tile
     integer            ::   im         ! horiz dimension and num of used pts         1
     integer            ::   km         ! vertical soil layer dimension               1
     real(kind_phys)    ::   grav       ! constant added to call in ccpp
     real(kind_phys)    ::   cp         ! constant added to call in ccpp
     real(kind_phys)    ::   hvap       ! constant added to call in ccpp
     real(kind_phys)    ::   rd         ! constant added to call in ccpp
     real(kind_phys)    ::   eps        ! constant added to call in ccpp
     real(kind_phys)    ::   epsm1      ! constant added to call in ccpp
     real(kind_phys)    ::   rvrdm1     ! constant added to call in ccpp
     real(kind_phys)    ::   delt       ! time interval (second)                      1
     integer            ::   isot       ! sfc soil type data source zobler or statsgo
     integer            ::   ivegsrc    ! sfc veg type data source umd or igbp
     logical            ::   lheatstrg  ! flag for canopy heat storage parameterization  1
     character(len=128) ::   errmsg     ! error messaging added to ccpp
     integer            ::   errflg     ! error messaging added to ccpp
     real(kind_phys)    ::   pertvegf
     logical            ::   thsfc_loc ! this should be changed to match same FV3 var

  end type lm4_static_type
  
  type sfcprop_type
     ! dims (im,ltile)
     real(kind_phys), allocatable  :: landfrac(:,:)
     real(kind_phys), allocatable  :: slmsk(:,:)
     real(kind_phys), allocatable  :: tsfcl(:,:)  
     real(kind_phys), allocatable  :: weasd(:,:)  ! aka sheleg in sfc file
     real(kind_phys), allocatable  :: tg3(:,:)
     real(kind_phys), allocatable  :: zorll(:,:)  ! note, z0rl over land
     real(kind_phys), allocatable  :: alvsf(:,:)
     real(kind_phys), allocatable  :: alvwf(:,:)
     real(kind_phys), allocatable  :: alnsf(:,:)
     real(kind_phys), allocatable  :: alnwf(:,:)
     real(kind_phys), allocatable  :: facsf(:,:)
     real(kind_phys), allocatable  :: facwf(:,:)
     real(kind_phys), allocatable  :: vfrac(:,:)
     real(kind_phys), allocatable  :: canopy(:,:)
     real(kind_phys), allocatable  :: f10m(:,:)
     real(kind_phys), allocatable  :: t2m(:,:)
     real(kind_phys), allocatable  :: q2m(:,:)
     real(kind_phys), allocatable  :: vtype(:,:)
     real(kind_phys), allocatable  :: stype(:,:)
     real(kind_phys), allocatable  :: uustar(:,:)
     real(kind_phys), allocatable  :: ffmm(:,:)
     real(kind_phys), allocatable  :: ffhh(:,:)
     real(kind_phys), allocatable  :: hice(:,:)
     real(kind_phys), allocatable  :: fice(:,:)
     real(kind_phys), allocatable  :: tisfc(:,:)
     real(kind_phys), allocatable  :: tprcp(:,:)
     real(kind_phys), allocatable  :: srflag(:,:)
     real(kind_phys), allocatable  :: snowd(:,:)  ! aka snwdph in sfc file
     real(kind_phys), allocatable  :: shdmin(:,:)
     real(kind_phys), allocatable  :: shdmax(:,:)
     real(kind_phys), allocatable  :: slope(:,:)
     real(kind_phys), allocatable  :: snoalb(:,:)
     real(kind_phys), allocatable  :: sncovr(:,:)

     ! JP TODO: allocate these properly
     real(kind_phys), allocatable  :: stc(:,:)
     real(kind_phys), allocatable  :: smc(:,:)
     real(kind_phys), allocatable  :: slc(:,:)     
  end type sfcprop_type

  type  :: lm4_model_type

     real(kind_phys), allocatable :: foo_atm2lndfield(:,:)
     ! from ufs-land-driver

     real(kind_phys), allocatable :: ps        (:,:) ! surface pressure (pa)                       im
     real(kind_phys), allocatable :: t1        (:,:) ! surface layer mean temperature (k)          im
     real(kind_phys), allocatable :: q1        (:,:) ! surface layer mean specific humidity        im
     integer        , allocatable :: soiltyp   (:,:) ! soil type (integer index)                   im
     integer        , allocatable :: vegtype   (:,:) ! vegetation type (integer index)             im
     real(kind_phys), allocatable :: sigmaf    (:,:) ! areal fractional cover of green vegetation  im
     real(kind_phys), allocatable :: sfcemis   (:,:) ! sfc lw emissivity ( fraction )              im
     real(kind_phys), allocatable :: dlwflx    (:,:) ! total sky sfc downward lw flux ( w/m**2 )   im
     real(kind_phys), allocatable :: dswsfc    (:,:) ! total sky sfc downward sw flux ( w/m**2 )   im
     real(kind_phys), allocatable :: dswsfci   (:,:) ! inst  sky sfc downward sw flux ( w/m**2 )   im
     real(kind_phys), allocatable :: snet      (:,:) ! total sky sfc netsw flx into ground(w/m**2) im
     real(kind_phys), allocatable :: tg3       (:,:) ! deep soil temperature (k)                   im
     real(kind_phys), allocatable :: cm        (:,:) ! surface exchange coeff for momentum (m/s)   im
     real(kind_phys), allocatable :: ch        (:,:) ! surface exchange coeff heat & moisture(m/s) im
     real(kind_phys), allocatable :: prsl1     (:,:) ! sfc layer 1 mean pressure (pa)              im
     real(kind_phys), allocatable :: prslki    (:,:) !                                             im
     real(kind_phys), allocatable :: zf        (:,:) ! height of bottom layer (m)                  im
     logical        , allocatable :: land      (:,:) ! = T if a point with any land                im
     real(kind_phys), allocatable :: wind      (:,:) ! wind speed (m/s)                            im
     integer        , allocatable :: slopetyp  (:,:) ! class of sfc slope (integer index)          im
     real(kind_phys), allocatable :: shdmin    (:,:) ! min fractional coverage of green veg        im
     real(kind_phys), allocatable :: shdmax    (:,:) ! max fractnl cover of green veg (not used)   im
     real(kind_phys), allocatable :: snoalb    (:,:) ! upper bound on max albedo over deep snow    im
     real(kind_phys), allocatable :: sfalb     (:,:) ! mean sfc diffused sw albedo (fractional)    im
     logical        , allocatable :: flag_iter (:,:) !                                             im
     logical        , allocatable :: flag_guess(:,:) !                                             im
     real(kind_phys), allocatable :: bexppert  (:,:)
     real(kind_phys), allocatable :: xlaipert  (:,:)
     real(kind_phys), allocatable :: vegfpert  (:,:)
     real(kind_phys), allocatable :: weasd     (:,:) ! water equivalent accumulated snow depth(mm) im
     real(kind_phys), allocatable :: snwdph    (:,:) ! snow depth (water equiv) over land          im
     real(kind_phys), allocatable :: tskin     (:,:) ! ground surface skin temperature ( k )       im
     real(kind_phys), allocatable :: tprcp     (:,:) ! total precipitation                         im
     real(kind_phys), allocatable :: srflag    (:,:) ! snow/rain flag for precipitation            im
     real(kind_phys), allocatable :: canopy    (:,:) ! canopy moisture content (m)                 im
     real(kind_phys), allocatable :: trans     (:,:) ! total plant transpiration (m/s)             im
     real(kind_phys), allocatable :: tsurf     (:,:) ! surface skin temperature (after iteration)  im
     real(kind_phys), allocatable :: z0rl      (:,:) ! surface roughness                           im
     real(kind_phys), allocatable :: sncovr1   (:,:) ! snow cover over land (fractional)            im
     real(kind_phys), allocatable :: qsurf     (:,:) ! specific humidity at sfc                     im
     real(kind_phys), allocatable :: gflux     (:,:) ! soil heat flux (w/m**2)                      im
     real(kind_phys), allocatable :: drain     (:,:) ! subsurface runoff (mm/s)                     im
     real(kind_phys), allocatable :: evap      (:,:) ! evaperation from latent heat flux            im
     real(kind_phys), allocatable :: hflx      (:,:) ! sensible heat flux                           im
     real(kind_phys), allocatable :: ep        (:,:) ! potential evaporation                        im
     real(kind_phys), allocatable :: runoff    (:,:) ! surface runoff (m/s)                         im
     real(kind_phys), allocatable :: cmm       (:,:) !                                              im
     real(kind_phys), allocatable :: chh       (:,:) !                                              im
     real(kind_phys), allocatable :: evbs      (:,:) ! direct soil evaporation (m/s)                im
     real(kind_phys), allocatable :: evcw      (:,:) ! canopy water evaporation (m/s)               im
     real(kind_phys), allocatable :: sbsno     (:,:) ! sublimation/deposit from snopack (m/s)       im
     real(kind_phys), allocatable :: snowc     (:,:) ! fractional snow cover                        im
     real(kind_phys), allocatable :: stm       (:,:) ! total soil column moisture content (m)       im
     real(kind_phys), allocatable :: snohf     (:,:) ! snow/freezing-rain latent heat flux (w/m**2) im
     real(kind_phys), allocatable :: smcwlt2   (:,:) ! dry soil moisture threshold                  im
     real(kind_phys), allocatable :: smcref2   (:,:) ! soil moisture threshold                      im
     real(kind_phys), allocatable :: wet1      (:,:) ! normalized soil wetness                      im
     real(kind_phys), allocatable :: prslk1    (:,:)

     ! JP TODO: allocate these properly
     real(kind_phys), allocatable :: smc(:,:) ! total soil moisture content (fractional)   im,km
     real(kind_phys), allocatable :: stc(:,:) ! soil temp (k)                              im,km
     real(kind_phys), allocatable :: slc(:,:) ! liquid soil moisture                       im,km

     ! rad
     real(kind_phys), allocatable :: albdvis_lnd (:,:)
     real(kind_phys), allocatable :: albdnir_lnd (:,:)
     real(kind_phys), allocatable :: albivis_lnd (:,:)
     real(kind_phys), allocatable :: albinir_lnd (:,:)
     real(kind_phys), allocatable :: adjvisbmd   (:,:)
     real(kind_phys), allocatable :: adjnirbmd   (:,:)
     real(kind_phys), allocatable :: adjvisdfd   (:,:)
     real(kind_phys), allocatable :: adjnirdfd   (:,:)

     ! from sfc_diff
     real(kind_phys), allocatable :: rb_lnd   (:,:)
     real(kind_phys), allocatable :: fm_lnd   (:,:)
     real(kind_phys), allocatable :: fh_lnd   (:,:)
     real(kind_phys), allocatable :: fm10_lnd (:,:)
     real(kind_phys), allocatable :: fh2_lnd  (:,:)
     real(kind_phys), allocatable :: stress   (:,:)  
     real(kind_phys), allocatable :: ustar    (:,:)
     real(kind_phys), allocatable :: garea    (:,:)                   

  end type lm4_model_type


  type, public :: lm4_type
     type(lm4_static_type)  :: static
     type(lm4_model_type)   :: model
     type(lm4_control_type) :: control
     type(sfcprop_type)      :: sfcprop
   contains

     procedure, public  :: Create

  end type lm4_type

contains

  subroutine Create(lm, im)

    implicit none

    class(lm4_type)     :: lm
    integer, intent(in) :: im

    integer,parameter   :: km = 4 ! tmp for testing. This should come from nml
    integer,parameter   :: ltile = 1
    
    ! --------------------------------------------
    lm%control%first_time = .true.
    lm%control%mype       = clear_val

    ! --------------------------------------------
    !lm%static%im         = clear_val
    lm%static%ltile      = ltile
    lm%static%km         = km 
    lm%static%grav       = clear_val
    lm%static%cp         = clear_val
    lm%static%hvap       = clear_val
    lm%static%rd         = clear_val
    lm%static%eps        = clear_val
    lm%static%epsm1      = clear_val
    lm%static%rvrdm1     = clear_val
    lm%static%delt       = clear_val
    lm%static%isot       = 1  ! TMP for testing. TODO: read in
    lm%static%ivegsrc    = 1  ! TMP for testing. TODO: read in
    lm%static%lheatstrg  = .false.
    lm%static%errmsg     = ""
    lm%static%errflg     = clear_val
    lm%static%pertvegf   = clear_val
    lm%static%thsfc_loc  = .true.
    ! --------------------------------------------
    ! --------------------------------------------    
    allocate(lm%model%foo_atm2lndfield        (im,ltile))
    allocate(lm%model%ps            (im,ltile))
    allocate(lm%model%t1            (im,ltile))
    allocate(lm%model%q1            (im,ltile))
    allocate(lm%model%soiltyp       (im,ltile))
    allocate(lm%model%vegtype       (im,ltile))
    allocate(lm%model%sigmaf        (im,ltile))
    allocate(lm%model%sfcemis       (im,ltile))
    allocate(lm%model%dlwflx        (im,ltile))
    allocate(lm%model%dswsfc        (im,ltile))
    allocate(lm%model%dswsfci       (im,ltile))
    allocate(lm%model%snet          (im,ltile))
    allocate(lm%model%tg3           (im,ltile))
    allocate(lm%model%cm            (im,ltile))
    allocate(lm%model%ch            (im,ltile))
    allocate(lm%model%prsl1         (im,ltile))
    allocate(lm%model%prslki        (im,ltile))
    allocate(lm%model%zf            (im,ltile))
    allocate(lm%model%land          (im,ltile))
    allocate(lm%model%wind          (im,ltile))
    allocate(lm%model%slopetyp      (im,ltile))
    allocate(lm%model%shdmin        (im,ltile))
    allocate(lm%model%shdmax        (im,ltile))
    allocate(lm%model%snoalb        (im,ltile))
    allocate(lm%model%sfalb         (im,ltile))
    allocate(lm%model%flag_iter     (im,ltile))
    allocate(lm%model%flag_guess    (im,ltile))
    allocate(lm%model%bexppert      (im,ltile))
    allocate(lm%model%xlaipert      (im,ltile))
    allocate(lm%model%vegfpert      (im,ltile))
    allocate(lm%model%weasd         (im,ltile))
    allocate(lm%model%snwdph        (im,ltile))
    allocate(lm%model%tskin         (im,ltile))
    allocate(lm%model%tprcp         (im,ltile))
    allocate(lm%model%srflag        (im,ltile))
    allocate(lm%model%canopy        (im,ltile))
    allocate(lm%model%trans         (im,ltile))
    allocate(lm%model%tsurf         (im,ltile))
    allocate(lm%model%z0rl          (im,ltile))
    allocate(lm%model%sncovr1       (im,ltile))
    allocate(lm%model%qsurf         (im,ltile))
    allocate(lm%model%gflux         (im,ltile))
    allocate(lm%model%drain         (im,ltile))
    allocate(lm%model%evap          (im,ltile))
    allocate(lm%model%hflx          (im,ltile))
    allocate(lm%model%ep            (im,ltile))
    allocate(lm%model%runoff        (im,ltile))
    allocate(lm%model%cmm           (im,ltile))
    allocate(lm%model%chh           (im,ltile))
    allocate(lm%model%evbs          (im,ltile))
    allocate(lm%model%evcw          (im,ltile))
    allocate(lm%model%sbsno         (im,ltile))
    allocate(lm%model%snowc         (im,ltile))
    allocate(lm%model%stm           (im,ltile))
    allocate(lm%model%snohf         (im,ltile))
    allocate(lm%model%smcwlt2       (im,ltile))
    allocate(lm%model%smcref2       (im,ltile))
    allocate(lm%model%wet1          (im,ltile))
    allocate(lm%model%albdvis_lnd   (im,ltile))
    allocate(lm%model%albdnir_lnd   (im,ltile))
    allocate(lm%model%albivis_lnd   (im,ltile))
    allocate(lm%model%albinir_lnd   (im,ltile))
    allocate(lm%model%adjvisbmd     (im,ltile))
    allocate(lm%model%adjnirbmd     (im,ltile))
    allocate(lm%model%adjvisdfd     (im,ltile))
    allocate(lm%model%adjnirdfd     (im,ltile))
    allocate(lm%model%prslk1        (im,ltile))
    allocate(lm%model%smc       (im,km))
    allocate(lm%model%stc       (im,km))
    allocate(lm%model%slc       (im,km))
    ! sfc_diff
    allocate(lm%model%rb_lnd        (im,ltile))
    allocate(lm%model%fm_lnd        (im,ltile))
    allocate(lm%model%fh_lnd        (im,ltile))
    allocate(lm%model%fm10_lnd      (im,ltile))
    allocate(lm%model%fh2_lnd       (im,ltile))
    allocate(lm%model%stress        (im,ltile))
    allocate(lm%model%ustar         (im,ltile))
    allocate(lm%model%garea         (im,ltile))
    
    !! Sfcprop -------------------------
    allocate(lm%sfcprop%landfrac (im,ltile))
    allocate(lm%sfcprop%slmsk    (im,ltile))      
    allocate(lm%sfcprop%tsfcl    (im,ltile))
    allocate(lm%sfcprop%weasd    (im,ltile))
    allocate(lm%sfcprop%tg3      (im,ltile))    
    allocate(lm%sfcprop%zorll    (im,ltile))     
    allocate(lm%sfcprop%alvsf    (im,ltile))      
    allocate(lm%sfcprop%alvwf    (im,ltile))      
    allocate(lm%sfcprop%alnsf    (im,ltile))      
    allocate(lm%sfcprop%alnwf    (im,ltile))      
    allocate(lm%sfcprop%facsf    (im,ltile))      
    allocate(lm%sfcprop%facwf    (im,ltile))      
    allocate(lm%sfcprop%vfrac    (im,ltile))      
    allocate(lm%sfcprop%canopy   (im,ltile))       
    allocate(lm%sfcprop%f10m     (im,ltile))     
    allocate(lm%sfcprop%t2m      (im,ltile))    
    allocate(lm%sfcprop%q2m      (im,ltile))    
    allocate(lm%sfcprop%vtype    (im,ltile))      
    allocate(lm%sfcprop%stype    (im,ltile))      
    allocate(lm%sfcprop%uustar   (im,ltile))       
    allocate(lm%sfcprop%ffmm     (im,ltile))     
    allocate(lm%sfcprop%ffhh     (im,ltile))     
    allocate(lm%sfcprop%hice     (im,ltile))     
    allocate(lm%sfcprop%fice     (im,ltile))     
    allocate(lm%sfcprop%tisfc    (im,ltile))      
    allocate(lm%sfcprop%tprcp    (im,ltile))      
    allocate(lm%sfcprop%srflag   (im,ltile))       
    allocate(lm%sfcprop%snowd    (im,ltile))
    allocate(lm%sfcprop%shdmin   (im,ltile))       
    allocate(lm%sfcprop%shdmax   (im,ltile))       
    allocate(lm%sfcprop%slope    (im,ltile))      
    allocate(lm%sfcprop%snoalb   (im,ltile))       
    allocate(lm%sfcprop%sncovr   (im,ltile))          
    
    allocate(lm%sfcprop%smc      (im,km))
    allocate(lm%sfcprop%stc      (im,km))
    allocate(lm%sfcprop%slc      (im,km))

    ! --------------------------------------------------------

    lm%model%foo_atm2lndfield  = clear_val
    lm%model%ps         = clear_val
    lm%model%t1         = clear_val
    lm%model%q1         = clear_val
    lm%model%soiltyp    = clear_val
    lm%model%vegtype    = clear_val
    lm%model%sigmaf     = clear_val
    lm%model%sfcemis    = clear_val
    lm%model%dlwflx     = clear_val
    lm%model%dswsfc     = clear_val
    lm%model%dswsfci    = clear_val
    lm%model%snet       = clear_val
    lm%model%tg3        = clear_val
    lm%model%cm         = clear_val
    lm%model%ch         = clear_val
    lm%model%prsl1      = clear_val
    lm%model%prslki     = clear_val
    lm%model%zf         = clear_val
    lm%model%land       = .false.
    lm%model%wind       = clear_val
    lm%model%slopetyp   = clear_val
    lm%model%shdmin     = clear_val
    lm%model%shdmax     = clear_val
    lm%model%snoalb     = clear_val
    lm%model%sfalb      = clear_val
    lm%model%flag_iter  = .false.
    lm%model%flag_guess = .false.
    lm%model%bexppert   = clear_val
    lm%model%xlaipert   = clear_val
    lm%model%vegfpert   = clear_val
    lm%model%weasd      = clear_val
    lm%model%snwdph     = clear_val
    lm%model%tskin      = clear_val
    lm%model%tprcp      = clear_val
    lm%model%srflag     = clear_val
    lm%model%canopy     = clear_val
    lm%model%trans      = clear_val
    lm%model%tsurf      = clear_val
    lm%model%z0rl       = clear_val
    lm%model%sncovr1    = clear_val
    lm%model%qsurf      = clear_val
    lm%model%gflux      = clear_val
    lm%model%drain      = clear_val
    lm%model%evap       = clear_val
    lm%model%hflx       = clear_val
    lm%model%ep         = clear_val
    lm%model%runoff     = clear_val
    lm%model%cmm        = clear_val
    lm%model%chh        = clear_val
    lm%model%evbs       = clear_val
    lm%model%evcw       = clear_val
    lm%model%sbsno      = clear_val
    lm%model%snowc      = clear_val
    lm%model%stm        = clear_val
    lm%model%snohf      = clear_val
    lm%model%smcwlt2    = clear_val
    lm%model%smcref2    = clear_val
    lm%model%wet1       = clear_val
    lm%model%smc        = clear_val
    lm%model%stc        = clear_val
    lm%model%slc        = clear_val     

    lm%model%rb_lnd     = clear_val
    lm%model%fm_lnd     = clear_val
    lm%model%fh_lnd     = clear_val
    lm%model%fm10_lnd   = clear_val
    lm%model%fh2_lnd    = clear_val
    lm%model%stress     = clear_val
    lm%model%ustar      = clear_val

    !! Surf Prop
    lm%sfcprop%landfrac    = clear_val
    lm%sfcprop%landfrac    = clear_val
    lm%sfcprop%slmsk       = clear_val
    lm%sfcprop%tsfcl       = clear_val
    lm%sfcprop%weasd       = clear_val
    lm%sfcprop%tg3         = clear_val
    lm%sfcprop%zorll       = clear_val
    lm%sfcprop%alvsf       = clear_val
    lm%sfcprop%alvwf       = clear_val
    lm%sfcprop%alnsf       = clear_val
    lm%sfcprop%alnwf       = clear_val
    lm%sfcprop%facsf       = clear_val
    lm%sfcprop%facwf       = clear_val
    lm%sfcprop%vfrac       = clear_val
    lm%sfcprop%canopy      = clear_val
    lm%sfcprop%f10m        = clear_val
    lm%sfcprop%t2m         = clear_val
    lm%sfcprop%q2m         = clear_val
    lm%sfcprop%vtype       = clear_val
    lm%sfcprop%stype       = clear_val
    lm%sfcprop%uustar      = clear_val
    lm%sfcprop%ffmm        = clear_val
    lm%sfcprop%ffhh        = clear_val
    lm%sfcprop%hice        = clear_val
    lm%sfcprop%fice        = clear_val
    lm%sfcprop%tisfc       = clear_val
    lm%sfcprop%tprcp       = clear_val
    lm%sfcprop%srflag      = clear_val
    lm%sfcprop%snowd       = clear_val
    lm%sfcprop%shdmin      = clear_val
    lm%sfcprop%shdmax      = clear_val
    lm%sfcprop%slope       = clear_val
    lm%sfcprop%snoalb      = clear_val
    lm%sfcprop%sncovr      = clear_val

    lm%sfcprop%smc         = clear_val
    lm%sfcprop%stc         = clear_val
    lm%sfcprop%slc         = clear_val     

    lm%model%garea       = zero
    
  end subroutine Create




end module lm4_type_mod
