      character in_surftemp_filename*200
      parameter(in_surftemp_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_surftemp_lonname*200
      parameter(in_surftemp_lonname=
     &  'longitude')
      character in_surftemp_latname*200
      parameter(in_surftemp_latname=
     &  'latitude')
      character in_surftemp_varname*200
      parameter(in_surftemp_varname=
     &  'temp_1')

      character in_sst_filename*200
      parameter(in_sst_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_sst_lonname*200
      parameter(in_sst_lonname=
     &  'longitude')
      character in_sst_latname*200
      parameter(in_sst_latname=
     &  'latitude')
      character in_sst_varname*200
      parameter(in_sst_varname=
     &  'temp_1')

      character in_precip_filename*200
      parameter(in_precip_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_precip_lonname*200
      parameter(in_precip_lonname=
     &  'longitude')
      character in_precip_latname*200
      parameter(in_precip_latname=
     &  'latitude')
      character in_precip_varname*200
      parameter(in_precip_varname=
     &  'precip')

      character in_netsolar_filename*200
      parameter(in_netsolar_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_netsolar_lonname*200
      parameter(in_netsolar_lonname=
     &  'longitude')
      character in_netsolar_latname*200
      parameter(in_netsolar_latname=
     &  'latitude')
      character in_netsolar_varname*200
      parameter(in_netsolar_varname=
     &  'field186')

      character in_netlong_filename*200
      parameter(in_netlong_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_netlong_lonname*200
      parameter(in_netlong_lonname=
     &  'longitude')
      character in_netlong_latname*200
      parameter(in_netlong_latname=
     &  'latitude')
      character in_netlong_varname*200
      parameter(in_netlong_varname=
     &  'field187')

      character in_sensible_filename*200
      parameter(in_sensible_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_sensible_lonname*200
      parameter(in_sensible_lonname=
     &  'longitude')
      character in_sensible_latname*200
      parameter(in_sensible_latname=
     &  'latitude')
      character in_sensible_varname*200
      parameter(in_sensible_varname=
     &  'sh')

      character in_latent_filename*200
      parameter(in_latent_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_latent_lonname*200
      parameter(in_latent_lonname=
     &  'longitude')
      character in_latent_latname*200
      parameter(in_latent_latname=
     &  'latitude')
      character in_latent_varname*200
      parameter(in_latent_varname=
     &  'lh')

      character in_stressx_filename*200
      parameter(in_stressx_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_stressx_lonname*200
      parameter(in_stressx_lonname=
     &  'longitude_1')
      character in_stressx_latname*200
      parameter(in_stressx_latname=
     &  'latitude_1')
      character in_stressx_varname*200
      parameter(in_stressx_varname=
     &  'field61')

      character in_stressy_filename*200
      parameter(in_stressy_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_stressy_lonname*200
      parameter(in_stressy_lonname=
     &  'longitude_1')
      character in_stressy_latname*200
      parameter(in_stressy_latname=
     &  'latitude_1')
      character in_stressy_varname*200
      parameter(in_stressy_varname=
     &  'field62')

      character in_runofff_filename*200
      parameter(in_runofff_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_runofff_lonname*200
      parameter(in_runofff_lonname=
     &  'longitude')
      character in_runofff_latname*200
      parameter(in_runofff_latname=
     &  'latitude')
      character in_runofff_varname*200
      parameter(in_runofff_varname=
     &  'field111')

      character in_runoffs_filename*200
      parameter(in_runoffs_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'xajpa_climate.nc')
      character in_runoffs_lonname*200
      parameter(in_runoffs_lonname=
     &  'longitude')
      character in_runoffs_latname*200
      parameter(in_runoffs_latname=
     &  'latitude')
      character in_runoffs_varname*200
      parameter(in_runoffs_varname=
     &  'field112')

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
     &                   'NOCVS_climate_xajpa_t21.nc')

      character out_sst_filename*200
      parameter(out_sst_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'NOCVS_sst_xajpa.nc')
