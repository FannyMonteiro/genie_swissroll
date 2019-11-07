


      integer nexp
      parameter(nexp=22)

      CHARACTER (len=2) :: my_name
      CHARACTER(LEN=20), PARAMETER :: FMT1 = "(I2.2)"

      character(len=200),parameter :: 
     &  in_lsm_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_bathymetry.nc'
      character in_lsm_varname*200
      parameter(in_lsm_varname=
     &  'bathymetry')

      character(len=200),parameter :: 
     &  in_orog_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_orog_'
      character in_orog_varname*200
      parameter(in_orog_varname=
     &  'orog')

      character(len=200),parameter :: 
     &  in_icefrac_filename =
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_icemask_'
      character in_icefrac_varname*200
      parameter(in_icefrac_varname=
     &  'icemask')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'icebcond_sealevel_21k.nc')

      character out_ascii_orogname*200
      parameter(out_ascii_orogname=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'orography_sealevel_21k.dat')

      character out_ascii_icefracname*200
      parameter(out_ascii_icefracname=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'icemask_sealevel_21k.dat')

      character out_data_filename_new*200
      parameter(out_data_filename_new=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'icebcond_sealevel_125k.nc')

      character out_ascii_orogname_new*200
      parameter(out_ascii_orogname_new=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'orography_sealevel_125k.dat')

      character out_ascii_icefracname_new*200
      parameter(out_ascii_icefracname_new=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'icemask_sealevel_125k.dat')
