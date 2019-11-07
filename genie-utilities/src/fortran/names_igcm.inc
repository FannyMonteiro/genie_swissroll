      character in_data_filename*200
      parameter(in_data_filename=
     &  '/home/ggdjl/genie_output/'//
     &  'genie_igcmfixed_t21_l7_makeclim'//
     &  '/igcm/igcm_cl_avg.nc')

      character in_sst_varname*200
      parameter(in_sst_varname=
     &  'surftemp')

      character in_precip1_varname*200
      parameter(in_precip1_varname=
     &  'convrain')

      character in_precip2_varname*200
      parameter(in_precip2_varname=
     &  'lscalerain')

      character in_sensible_varname*200
      parameter(in_sensible_varname=
     &  'sensheat')

      character in_latent_varname*200
      parameter(in_latent_varname=
     &  'latentheat')

      character in_netsolar_varname*200
      parameter(in_netsolar_varname=
     &  'netsolar')

      character in_netlong_varname*200
      parameter(in_netlong_varname=
     &  'netlong')

      character in_stressx_varname*200
      parameter(in_stressx_varname=
     &  'stressx')

      character in_stressy_varname*200
      parameter(in_stressy_varname=
     &  'stressy')

      character in_runoff_varname*200
      parameter(in_runoff_varname=
     &  'runoff')



      character in_data_maskname*200
      parameter(in_data_maskname=
     &  '/home/ggdjl/genie/'//
     &  'genie-igcm3/inputdata'//
     &  '/landmask_std_t21.nc')

      character in_mask_varname*200
      parameter(in_mask_varname=
     &  'lsm')




      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie_utilities/dataout/'//
     &                   'igcm_climate.nc')

