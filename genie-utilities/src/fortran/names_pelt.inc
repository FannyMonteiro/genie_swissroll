      character in_relief_filename*200
      parameter(in_relief_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_orog')
      character in_relief_lonname*200
      parameter(in_relief_lonname=
     &  'longitude')
      character in_relief_latname*200
      parameter(in_relief_latname=
     &  'latitude')
      character in_relief_varname*200
      parameter(in_relief_varname=
     &  'orog')

      character in_landice_filename*200
      parameter(in_landice_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_icemask')
      character in_landice_lonname*200
      parameter(in_landice_lonname=
     &  'longitude')
      character in_landice_latname*200
      parameter(in_landice_latname=
     &  'latitude')
      character in_landice_varname*200
      parameter(in_landice_varname=
     &  'icemask')

      character in_vege_filename*200
      parameter(in_vege_filename=
     &  '/home/ggdjl/genie_utilities/datain/cnt/'//
     &                   'vegetation_new_t21.nc')
      character in_vege_lonname*200
      parameter(in_vege_lonname=
     &  'longitude')
      character in_vege_latname*200
      parameter(in_vege_latname=
     &  'latitude')
      character in_vege_varname*200
      parameter(in_vege_varname=
     &  'svege')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'limit_pelt_ice4g_new.nc')
