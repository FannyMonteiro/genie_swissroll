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

      character in_seaicefrac_filename*200
      parameter(in_seaicefrac_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &  'AMIP.nc')
      character in_seaicefrac_lonname*200
      parameter(in_seaicefrac_lonname=
     &  'longitude')
      character in_seaicefrac_latname*200
      parameter(in_seaicefrac_latname=
     &  'latitude')
      character in_seaicefrac_varname*200
      parameter(in_seaicefrac_varname=
     &  'SEA_ICE')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'sst_amip.nc')
