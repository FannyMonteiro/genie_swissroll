      character in_relief_filename*200
      parameter(in_relief_filename=
     &  '/home/ggdjl/genie_utilities/datain/cnt/'//
     &                   'Relief.nc')
      character in_relief_lonname*200
      parameter(in_relief_lonname=
     &  'longitude')
      character in_relief_latname*200
      parameter(in_relief_latname=
     &  'latitude')
      character in_relief_varname*200
      parameter(in_relief_varname=
     &  'RELIEF')


      character in_sst_filename*200
      parameter(in_sst_filename=
     &  '/home/ggdjl/genie_utilities/datain/cnt/'//
     &                   'amipbc_sst_1x1_clim.nc')
      character in_sst_lonname*200
      parameter(in_sst_lonname=
     &  'longitude')
      character in_sst_latname*200
      parameter(in_sst_latname=
     &  'latitude')
      character in_sst_varname*200
      parameter(in_sst_varname=
     &  'SSTamip')



      character in_seaice_filename*200
      parameter(in_seaice_filename=
     &  '/home/ggdjl/genie_utilities/datain/cnt/'//
     &                   'AMIP.nc')
      character in_seaice_lonname*200
      parameter(in_seaice_lonname=
     &  'longitude')
      character in_seaice_latname*200
      parameter(in_seaice_latname=
     &  'latitude')
      character in_seaice_varname*200
      parameter(in_seaice_varname=
     &  'SEA_ICE')


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


      character in_landice_filename*200
      parameter(in_landice_filename=
     &  '/home/ggdjl/genie_utilities/datain/cnt/'//
     &                   'landiceref.nc')
      character in_landice_lonname*200
      parameter(in_landice_lonname=
     &  'longitude')
      character in_landice_latname*200
      parameter(in_landice_latname=
     &  'latitude')
      character in_landice_varname*200
      parameter(in_landice_varname=
     &  'landice')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie_utilities/dataout/'//
     &                   'limit_cnt_new.nc')
