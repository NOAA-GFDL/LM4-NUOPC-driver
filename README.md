# LM4-NUOPC-driver
This is the repository for the NUOPC driver of the GFDL Land Model version 4 (LM4), for use in the Unified Forecast System


--------------------------------------------------------
LM4 itself is a submodule of this repository.
More information and the LM4 project code is made available through GitHub at https://github.com/NOAA-GFDL/LM4,
but is managed by NOAA-GFDL at https://gitlab.gfdl.noaa.gov.

--------------------------------------------------------


Supported grids and resolutions
--------------------------------------------------------

Global cubed-sphere C96, C48, and C24 grids are supported. LM4 runs either with the CDEPS data atmosphere (`-DAPP=LND-LM4`) or coupled to the active FV3 atmosphere (`-DAPP=ATML-LM4`), through the CMEPS mediator.

Namelist Options
--------------------------------------------------------

Namelist options for tests are set in `tests/parm/input_datm_lm4.nml.IN` (data atmosphere) and `tests/parm/input_atm_lm4.nml.IN` (active atmosphere)

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
      cpl2atm     = .false.     ! coupling to active atmosphere
      kinematic_flux = .true.   ! some fluxes to atm are expressed as kinematic
      implicit_atm   = .false.  ! implicit coupling to atm (not functional, .true. aborts)

Additional, `flux_exchange_nml` and `atmos_prescr_nml` are existing namelists in `input.nml` at GFDL, and are needed by LM4 in UFS.
They have been brought back in, only for the options used in LM4, and are read in by the land model.

- `flux_exchange_nml` has the options `z_ref_heat`, the reference height for temperature and relative humidity diagnostics (default 2 m), and `z_ref_mom`, the reference height for momentum diagnostics (default 10 m)
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
          cpl2atm = .false.
    /
    
    &surface_flux_nml
          alt_gustiness = .TRUE.
          !gust_min      = 0.0    ! m/s, min. gustiness used when alt_gustiness is .FALSE.
    /

    &atmos_prescr_nml
      gust_to_use  = 'prescribed'
          !gustiness    = 5.0     ! m/s, wind gustiness if gust_to_use = 'prescribed'
          !gust_min     = 0.01    ! m/s, min. gustiness when gust_to_use = 'computed'
    /



Restarts
--------------------------------------------------------

LM4 uses FMS to write and read restarts.

- The following restart files are needed generated in the `RESTART` directory, and needed in the `INPUT` directory for a warm start
```
  cana.res.tile{1..6}.nc 
  glac.res.tile{1..6}.nc 
  lake.res.tile{1..6}.nc 
  land.res.tile{1..6}.nc 
  snow.res.tile{1..6}.nc 
  soil.res.tile{1..6}.nc 
  vegn1.res.tile{1..6}.nc
  vegn2.res.tile{1..6}.nc
  river.res.tile{1..6}.nc   (only after a river slow time step, river_nml dt_slow)
  landuse.res	  
```


Intermediate restarts can be generate during the model run
with the `lm4_nml` name list option:

 `restart_interval = Y,M,D,h,m,s` in `input.nml`

In the regression tests this is `LM4_RESTART_INTERVAL="0,0,0,6,0,0"`, ie, 6 hourly. Intermediate restarts are prefixed `YYYYMMDD.HHMMSS.`

Note that LM4.0 will not restart accurately unless restarts are at time 00z

Fast and Slow Land Timesteps
--------------------------------------------------------

LM4 has fast and slow processes. 
- The fast timestep is obtained from the NUOPC run sequence.
- The slow timstep is currently set by a namelist option in `lm4_nml`. The default is:
  `dt_lnd_slow = 3600` 



Running with regression tests
--------------------------------------------------------

Data atmosphere (`LND-LM4`, C96/C48/C24):

- `datm_cdeps_lm4_gswp3{,_c48,_c24}`: 48 hr cold start from 2000-01-01 00:00 using CDEPS with GSWP3 forcing
- `datm_cdeps_lm4_gswp3_rst{,_c48,_c24}`: restart test using restarts from the control

Active atmosphere (`ATML-LM4`, GFS v17 P8 physics, C96/C48):

- `control_p8_atmlm4{,_c48}`: 48 hr cold start from 2021-03-22 00:00
- `control_restart_p8_atmlm4{,_c48}`: restart test using restarts from the control

LM4 restart files (and river restarts, when written) are used for comparisons.

These tests use the following unique files:

    export UFS_CONFIGURE="ufs.configure.atm_lm4.IN"
    export FV3_RUN="lm4_common.IN lm4_datm_cdeps_gswp.IN"   # data atmosphere
    export FV3_RUN="control_run.IN lm4_common.IN"           # active atmosphere
    export DIAG_TABLE="diag_table_datm_lm4"                 # data atmosphere
    export DIAG_TABLE_ADDITIONAL=diag_additional_cpld_lm4.IN # active atmosphere
    export FIELD_TABLE_ADDITIONAL=field_table_lm4
    export INPUT_NML="input_datm_lm4.nml.IN"                # or input_atm_lm4.nml.IN; C24/C48 tests call export_lm4_lowres

LM4 input data is read from `INPUTDATA_LM4` (`${INPUTDATA_ROOT}/LM4_input_data`): C96 grid and hydrography from `c96_grid/` and `c96_LM4/`, and lower resolution exchange grids from `a<LNDRES>o<OCNRES>/`.

Imports and Exports
--------------------------------------------------------

The following are import and export variables. Imports are similar for both atmospheres.

| Import      | Description                                      | Variable Name       |
|-------------|--------------------------------------------------|---------------------|
| Sa_z        | Bottom layer height                              | lm4_model%atm_forc%z_bot |
| Sa_tbot     | Bottom layer temperature                         | lm4_model%atm_forc%t_bot |
| Sa_pbot     | Bottom layer pressure                            | lm4_model%atm_forc%p_bot |
| Sa_u        | Bottom layer zonal wind                          | lm4_model%atm_forc%u_bot |
| Sa_v        | Bottom layer meridional wind                     | lm4_model%atm_forc%v_bot |
| Sa_shum     | Bottom layer specific humidity                   | lm4_model%atm_forc%q_bot |
| Sa_pslv     | Surface pressure                                 | lm4_model%atm_forc%p_surf |
| Faxa_lwdn   | Mean surface downward longwave radiation         | lm4_model%atm_forc%flux_lw |
| Faxa_swvdf  | Mean surface downward UV+visible diffuse flux    | lm4_model%atm_forc%flux_sw_down_vis_dif |
| Faxa_swvdr  | Mean surface downward UV+visible direct flux     | lm4_model%atm_forc%flux_sw_down_vis_dir |
| Faxa_swndf  | Mean surface downward NIR diffuse flux           | lm4_model%atm_forc%flux_sw_down_nir_dif |
| Faxa_swndr  | Mean surface downward NIR direct flux            | lm4_model%atm_forc%flux_sw_down_nir_dir |

| Export      | Description                                      | Variable Name       |
|-------------|--------------------------------------------------|---------------------|      
| Sl_lfrin    | Land fraction for mediator                       | lnd%sg_landfrac     |

With an active atmosphere (`cpl2atm = .true.`), the active-atmosphere imports (`Sa_ta`, `Sa_qa`, `Sa_prsl`, `Sa_exner`, `Sa_ustar`, `Sa_tskn`, `Sa_zorl`, `Sa_vfrac`, `Faxa_swdn`, `Faxa_swnet`, `Faxa_rain`, `Faxa_snow`) are also used, and LM4 computes the land surface fluxes and exports:

| Export      | Description                                      |
|-------------|--------------------------------------------------|
| Sl_sfrac    | Snow fraction                                    |
| Fall_lat    | Latent heat flux                                 |
| Fall_sen    | Sensible heat flux                               |
| Fall_evap   | Evaporation                                      |
| Fall_gflx   | Ground heat flux                                 |
| Fall_roff   | Subsurface runoff                                |
| Fall_soff   | Surface runoff                                   |
| Sl_q        | Canopy-air (surface) specific humidity           |
| Sl_tref, Sl_qref | Reference height temperature and humidity   |
| Sl_cmm, Sl_chh   | Momentum and heat exchange coefficients     |
| Sl_zvfun    | Vegetation function for FV3 surface layer        |
| Sl_avsdr, Sl_anidr, Sl_avsdf, Sl_anidf | Visible/NIR direct/diffuse albedos |


Surface boundary layer and coupling with the atmosphere
--------------------------------------------------------

Each fast (coupling) time step, `ModelAdvance` in `lm4_cap.F90` follows the GFDL FMS coupler sequence (adapted from `atm_land_ice_flux_exchange`, stripped down to land only on the LM4 unstructured grid):

```fortran
      call sfc_boundary_layer(real(sec), lm4_model)        ! surface fluxes and their derivatives
      call update_atmos_model_down(lm4_model)              ! gustiness
      call flux_down_from_atmos(real(sec), lm4_model)      ! fluxes and precip passed down to land
      call update_land_model_fast(lm4_model%From_atm, lm4_model%From_lnd)
      if (lm4_model%nml%cpl2atm) then                      ! active 2-way coupling with atm
         call flux_up_to_atmos(lm4_model)                  ! fluxes and states passed up to atm
      end if
      ...
      if (mod(time_sec+timestep_sec, dt_lnd_slow) == 0) call update_land_model_slow(...)
      call export_fields(gcomp, lm4_model, rc)
```

- `sfc_boundary_layer` computes surface-layer exchange fluxes, and their derivatives on each land tile. With `cpl2atm`, it among other things, aggregates the tile-fraction-weighted visible/NIR, direct/diffuse land albedos
- `update_atmos_model_down` only computes gustiness. For an active atmosphere, `alt_gustiness = .true.` 
- `flux_down_from_atmos` passes radiation, precipitation, and fluxes to land with explicit coupling. `implicit_atm = .true.` is not functional yet.
- `flux_up_to_atmos` (only with `cpl2atm`) corrects the fluxes for the new land surface and canopy-air temperatures and tracers after `update_land_model_fast`, then aggregates over land tiles to the grid cell:
  - sensible heat flux (`Fall_sen`) and water vapor flux (`Fall_lat`)
  - canopy-air specific humidity (`Sl_q`)
  - with `kinematic_flux = .true.` (default), sensible heat is converted to `SH/(rho*cp)` and water vapor flux to `E/rho`, using air density at the lowest atmospheric level, matching what the UFS atmosphere expects
- `export_fields` sends `Sl_lfrin` and, with `cpl2atm`, the fluxes, states, albedos, and surface-layer quantities above to CMEPS, which passes them to FV3. On the atmosphere side, `cpllnd = cpllnd2atm = .true.` in `gfs_physics_nml` replaces the land surface model over land with these imports.

`Sl_sfrac`, `Fall_gflx`, `Fall_roff`, `Fall_soff`, and `Fall_evap` are advertised but currently exported as 0.
