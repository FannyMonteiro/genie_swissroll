      character in_sst_filename*200
      parameter(in_sst_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'amipbc_sst_1x1_clim.nc')
      character in_sst_lonname*200
      parameter(in_sst_lonname=
     &  'longitude')
      character in_sst_latname*200
      parameter(in_sst_latname=
     &  'latitude')
      character in_sst_varname*200
      parameter(in_sst_varname=
     &  'SSTamip')

      character in_precip_filename*200
      parameter(in_precip_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/prate.mon.ltm.nc')
      character in_precip_lonname*200
      parameter(in_precip_lonname=
     &  'lon')
      character in_precip_latname*200
      parameter(in_precip_latname=
     &  'lat')
      character in_precip_varname*200
      parameter(in_precip_varname=
     &  'prate')

      character in_netsolar_filename*200
      parameter(in_netsolar_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/nswrs.mon.ltm.nc')
      character in_netsolar_lonname*200
      parameter(in_netsolar_lonname=
     &  'lon')
      character in_netsolar_latname*200
      parameter(in_netsolar_latname=
     &  'lat')
      character in_netsolar_varname*200
      parameter(in_netsolar_varname=
     &  'nswrs')

      character in_netlong_filename*200
      parameter(in_netlong_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/nlwrs.mon.ltm.nc')
      character in_netlong_lonname*200
      parameter(in_netlong_lonname=
     &  'lon')
      character in_netlong_latname*200
      parameter(in_netlong_latname=
     &  'lat')
      character in_netlong_varname*200
      parameter(in_netlong_varname=
     &  'nlwrs')

      character in_sensible_filename*200
      parameter(in_sensible_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/shtfl.mon.ltm.nc')
      character in_sensible_lonname*200
      parameter(in_sensible_lonname=
     &  'lon')
      character in_sensible_latname*200
      parameter(in_sensible_latname=
     &  'lat')
      character in_sensible_varname*200
      parameter(in_sensible_varname=
     &  'shtfl')

      character in_latent_filename*200
      parameter(in_latent_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/lhtfl.mon.ltm.nc')
      character in_latent_lonname*200
      parameter(in_latent_lonname=
     &  'lon')
      character in_latent_latname*200
      parameter(in_latent_latname=
     &  'lat')
      character in_latent_varname*200
      parameter(in_latent_varname=
     &  'lhtfl')

      character in_stressx_filename*200
      parameter(in_stressx_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/uflx.mon.ltm.nc')
      character in_stressx_lonname*200
      parameter(in_stressx_lonname=
     &  'lon')
      character in_stressx_latname*200
      parameter(in_stressx_latname=
     &  'lat')
      character in_stressx_varname*200
      parameter(in_stressx_varname=
     &  'uflx')

      character in_stressy_filename*200
      parameter(in_stressy_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/vflx.mon.ltm.nc')
      character in_stressy_lonname*200
      parameter(in_stressy_lonname=
     &  'lon')
      character in_stressy_latname*200
      parameter(in_stressy_latname=
     &  'lat')
      character in_stressy_varname*200
      parameter(in_stressy_varname=
     &  'vflx')

      character in_runoff_filename*200
      parameter(in_runoff_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/runof.mon.ltm.nc')
      character in_runoff_lonname*200
      parameter(in_runoff_lonname=
     &  'lon')
      character in_runoff_latname*200
      parameter(in_runoff_latname=
     &  'lat')
      character in_runoff_varname*200
      parameter(in_runoff_varname=
     &  'runof')

      character in_data_maskname*200
      parameter(in_data_maskname=
     &  '/home/ggdjl/genie/'//
     &  'genie-igcm3/data/input'//
     &  '/landmask_std_t21.nc')
      character in_mask_varname*200
      parameter(in_mask_varname=
     &  'lsm')

      character in_data_indexname*200
      parameter(in_data_indexname=
     &  '/home/ggdjl/genie/'//
     &  'genie-utilities/data/output'//
     &  '/NOCVS_runoff_new_t21.nc')
      character in_index_varname*200
      parameter(in_index_varname=
     &  'destination_indices')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'NOCVS_limit_fix_ncp.nc')
