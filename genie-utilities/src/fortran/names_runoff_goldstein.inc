      character in_data_maskname*200
      parameter(in_data_maskname=
     &  '/home/ggdjl/genie/'//
     &  'genie-utilities/data/output'//
     &  '/limit_goldstein.nc')
      character in_mask_varname*200
      parameter(in_mask_varname=
     &  'lsm')

      character in_data_runoffname*200
      parameter(in_data_runoffname=
     &  '/home/ggdjl/genie/'//
     &  'genie-land/data/input'//
     &  '/runoff_mask_igcm64x32.nc')
      character in_runoff_varname*200
      parameter(in_runoff_varname=
     &  'destination_indices')

      character out_data_filename_1*200
      parameter(out_data_filename_1=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'runoff_goldstein_ob_t21.nc')

      character out_data_filename_2*200
      parameter(out_data_filename_2=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'runoff_goldstein_nn_t21.nc')


