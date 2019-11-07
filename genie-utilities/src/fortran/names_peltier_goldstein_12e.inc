
      integer nexp
      parameter(nexp=1)
      real,parameter,dimension(nexp) ::
     &  expname=(/
     &  12
     &  /)

      character(len=200),parameter ::
     &  in_lsm_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_bathymetry.nc'
      character in_lsm_varname*200
      parameter(in_lsm_varname=
     &  'bathymetry')

      character(len=200), parameter :: in_orog_filename(nexp)=(/
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_orog_12.nc'
     &                              /)
      character in_orog_varname*200
      parameter(in_orog_varname=
     &  'orog')

      character(len=200), parameter :: in_icefrac_filename(nexp)=(/
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_icemask_12.nc'
     &                              /)
      character in_icefrac_varname*200
      parameter(in_icefrac_varname=
     &  'icemask')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'peltier_12e_goldstein.nc')

      character out_ascii_orogname*200
      parameter(out_ascii_orogname=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'orography_12e.dat')

      character out_ascii_icefracname*200
      parameter(out_ascii_icefracname=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'icemask_12e.dat')
