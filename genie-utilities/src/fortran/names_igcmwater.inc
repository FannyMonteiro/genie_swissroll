      character in_data_filename*200
      parameter(in_data_filename=
c     &  '/export/allerod/allerod-02/work/ggdjl/'//
c     &  'genie_output/'//
c     &  'genie_ig_go_sl_control'//
c     &  '/main/genie_2449_12_30_atm.nc')
c     &  '/export/stage3/stage3-01/work/ggdjl/'//
c     &  'genie_output/'//
c     &  'genie_ig_go_sl_gawater'//
c     &  '/main/genie_2449_12_30_atm.nc')
c     &  '/export/lgm/lgm-02/work/ggdjl/'//
c     &  'genie_output/'//
c     &  'genie_ig_fi_fi'//
c     &  '/main/igcm_cl_ann_avg.nc')
c     &  '/export/lgm/lgm-02/work/ggdjl/'//
c     &  'genie_output/'//
c     &  'genie_ig_go_sl_gaalbedo'//
c     &  '/main/genie_2499_12_30_atm.nc')
     &  '/export/eemian/eemian-08/work/ggdjl/'//
     &  'genie_output/'//
     &  'genie_ig_go_sl_gaalbedofluxcorr4'//
     &  '/main/genie_2499_12_30_atm.nc')

      character in_precip_varname*200
      parameter(in_precip_varname=
     &  'precip')

      character in_evap_varname*200
      parameter(in_evap_varname=
     &  'evap')

      character in_runoff_varname*200
      parameter(in_runoff_varname=
     &  'runoff')


      character in_data_maskname*200
      parameter(in_data_maskname=
     &  '/home/ggdjl/genie/'//
     &  'genie-fixedicesheet/data/input'//
     &  '/landmask_std_t21.nc')

      character in_mask_varname*200
      parameter(in_mask_varname=
     &  'lsm')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie/genie-utilities/'//
     &                   'data/output/'//
     &                   'NOCVS_igcm_weights.nc')
