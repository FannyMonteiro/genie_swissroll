      character in_data_maskname*200
      parameter(in_data_maskname=
     &  '/home/ggdjl/genie/'//
     &  'genie-igcm3/data/input'//
     &  '/landmask_std_t21.nc')

      character in_mask_varname*200
      parameter(in_mask_varname=
     &  'lsm')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie/genie-utilities/'//
     &                   'data/output/'//
     &                   'NOCVS_igcm_fluxadjust.nc')
