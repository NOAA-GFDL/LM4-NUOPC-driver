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

> cana.res.tile1.nc 
> cana.res.tile2.nc 
> cana.res.tile3.nc 
> cana.res.tile4.nc 
> cana.res.tile5.nc 
> cana.res.tile6.nc 
> glac.res.tile1.nc 
> glac.res.tile2.nc 
> glac.res.tile3.nc 
> glac.res.tile4.nc 
> glac.res.tile5.nc 
> glac.res.tile6.nc 
> lake.res.tile1.nc 
> lake.res.tile2.nc 
> lake.res.tile3.nc 
> lake.res.tile4.nc 
> lake.res.tile5.nc 
> lake.res.tile6.nc 
> land.res.tile1.nc 
> land.res.tile2.nc 
> land.res.tile3.nc 
> land.res.tile4.nc 
> land.res.tile5.nc 
> land.res.tile6.nc 
> landuse.res	  
> snow.res.tile1.nc 
> snow.res.tile2.nc 
> snow.res.tile3.nc 
> snow.res.tile4.nc 
> snow.res.tile5.nc 
> snow.res.tile6.nc 
> soil.res.tile1.nc 
> soil.res.tile2.nc 
> soil.res.tile3.nc 
> soil.res.tile4.nc 
> soil.res.tile5.nc 
> soil.res.tile6.nc 
> vegn1.res.tile1.nc
> vegn1.res.tile2.nc
> vegn1.res.tile3.nc
> vegn1.res.tile4.nc
> vegn1.res.tile5.nc
> vegn1.res.tile6.nc
> vegn2.res.tile1.nc
> vegn2.res.tile2.nc
> vegn2.res.tile3.nc
> vegn2.res.tile4.nc
> vegn2.res.tile5.nc
> vegn2.res.tile6.nc


Intermediate restarts can be generate during the model run with the name list option `restart_interval = Y,M,D,h,m,s` in `lm4_nml`

-  default is `restart_interval = 0,0,0,6,0,0`

Note that LM4.0 will not restart accurately unless restarts are at time 00z

Fast and Slow Land Timesteps
--------------------------------------------------------

Running with regression tests
--------------------------------------------------------

There are currently 2 Data-atmopshere tests:

- `datm_cdeps_lm4_c48_gswp3`
  - 48 hr test from 2000-01-01 00:00 with a cold start
  - using CDEPS with GSWP3 forcing  
- `datm_cdeps_lm4_c48_gswp3_rst` 
   - a 24 hr restart test from 2000-01-02 00:00 
   - warm start using restarts from `datm_cdeps_lm4_c48_gswp3`

