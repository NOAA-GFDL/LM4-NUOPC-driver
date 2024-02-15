# LM4-NUOPC-driver
This is the repository for the NUOPC driver of the GFDL Land Model version 4 (LM4), for use in the Unified Forecast System


--------------------------------------------------------
LM4 itself is a submodule of this repository.
More information and the LM4 project code is made available through GitHub at https://github.com/NOAA-GFDL/LM4,
but is managed by NOAA-GFDL at https://gitlab.gfdl.noaa.gov.

--------------------------------------------------------


Supported grids and resolutions
--------------------------------------------------------

Currently in development, only global C48 and C96 grids supported

Namelist Options
--------------------------------------------------------

namelist options are in `tests/parm/input_datm_lm4.nml.IN`

It adds the namelist `lm4_nml`, along with several LM4 and a couple FMS coupler / flux exchange namelists


[TODO: detail key options]

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
- The slow timstep is currently set by a namelist option, the default is:
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
