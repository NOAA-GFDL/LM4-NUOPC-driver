# LM4-NUOPC-driver
This is the repository for the NUOPC driver of the GFDL Land Model version 4 (LM4), for use in the Unified Forecast System


--------------------------------------------------------
LM4 itself is a submodule of this repository.
More information and the LM4 project code is made available through GitHub at https://github.com/NOAA-GFDL/LM4,
but is managed by NOAA-GFDL at https://gitlab.gfdl.noaa.gov.

--------------------------------------------------------


Supported grids and resolutions
--------------------------------------------------------

Currently in development, global C48 and C96 grids are supported

Namelist Options
--------------------------------------------------------

namelist options for tests are set `tests/parm/input_datm_lm4.nml.IN`

It adds the namelist `lm4_nml`, along with several LM4 and a couple FMS coupler / flux exchange namelists

Defaults for the `lm4_nml` if not specified in the run directory `input.nml` are:

      lm4_debug   = 0           ! debug flag for lm4 (0=off, 1=low, 2=high)
      npx         = 0, npy = 0  ! for FMS domain
      ntiles      = 0           ! cubed sphere tiles
      layout(2)   = (/0,0/)     ! layout for decompostition
      grid        = 'none'
      blocksize   = -1          ! for FMS domain blocking
      dt_lnd_slow = 86400       ! time step for slow land processes (s)
      restart_interval = (/ 0, 0, 0, 0, 0, 0/) ! The time interval that write out intermediate restart file.
                                               ! The format is (yr,mo,day,hr,min,sec).  When restart_interval
                                               ! is all zero, no intermediate restart file will be written out

Additional, `flux_exchange_nml` and `atmos_prescr_nml` are existing namelists in `input.nml` at GFDL, and are needed by LM4 in UFS.
They have been brought back in, only for the options used in LM4, and are read in by the land model.

- `flux_exchange_nml` has the option `z_ref_heat`, the reference height for for temperature and relative humidity diagnostics (default 2 m )
- `atmos_prescr_nml` has the options:
  - `gust_to_use` computed or prescribed gust for LM4 boundary layer ( default computed)
  - `gustiness`   wind gustiness (default 5.0 m/s)
  - `gust_min`  minimum gustiness used in computed gust calculation (default 0.0 m/s)

- Finally `surface_flux_nml` is used by LM4 UFS. It's options are:

```fortran
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
```


In `tests/parm/input_datm_lm4.nml.IN`, these options are set to: 

    &lm4_nml
          lm4_debug = 1
          grid = 'CS'
          layout = @[INPES],@[JNPES]
          npx = @[NPX]
          npy = @[NPY]
          ntiles = 6
          blocksize = -1
          dt_lnd_slow = 3600
          restart_interval = 0,0,0,6,0,0
    /
    
    &surface_flux_nml
          alt_gustiness = .FALSE.
          gust_min = 1.0
    /



Restarts
--------------------------------------------------------

LM4 uses FMS to write and read restarts.

- The following restart files are needed generated in the `RESTART` directory, and needed in the `INPUT` directory for a warm start

  - cana.res.tile{1..6}.nc 
  - glac.res.tile{1..6}.nc 
  - lake.res.tile{1..6}.nc 
  - land.res.tile{1..6}.nc 
  - snow.res.tile{1..6}.nc 
  - soil.res.tile{1..6}.nc 
  - vegn1.res.tile{1..6}.nc
  - vegn2.res.tile{1..6}.nc
  - landuse.res	  



Intermediate restarts can be generate during the model run
with the `lm4_nml` name list option:

 `restart_interval = Y,M,D,h,m,s` in `input.nml`

The default is `restart_interval = 0,0,0,6,0,0`, ie, 6 hourly

Note that LM4.0 will not restart accurately unless restarts are at time 00z

Fast and Slow Land Timesteps
--------------------------------------------------------

LM4 has fast and slow processes. 
- The fast timestep is obtained from the NUOPC run sequence.
- The slow timstep is currently set by a namelist option in `lm4_nml`. The default is:
  `dt_lnd_slow = 3600` 



Running with regression tests
--------------------------------------------------------

There are currently 2 Data-atmosphere tests:

- `datm_cdeps_lm4_c48_gswp3`
  - 48 hr test from 2000-01-01 00:00 with a cold start
  - using CDEPS with GSWP3 forcing  

- `datm_cdeps_lm4_c48_gswp3_rst` 
   - a 24 hr restart test from 2000-01-02 00:00 
   - warm start using restarts from `datm_cdeps_lm4_c48_gswp3`

For both, the LM4 restart files listed above are used for comparisons 

These tests use the following unique files:

    export UFS_CONFIGURE="ufs.configure.atm_lm4.IN"
    export FV3_RUN="lm4_run.IN"
    export DIAG_TABLE="diag_table_datm_lm4"
    export FIELD_TABLE_ADDITIONAL=field_table_lm4
    export INPUT_NML="input_datm_lm4.nml.IN"


Surface boundary layer and connection to Land
--------------------------------------------------------

TODO: brief overview of call order:

```fortran
      call sfc_boundary_layer(real(sec), lm4_model)
      call update_atmos_model_down(lm4_model)              ! for gust calculation with data atmosphere
      call flux_down_from_atmos(real(sec), lm4_model)      
```