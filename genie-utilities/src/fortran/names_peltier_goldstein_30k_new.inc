


      integer nexp
      parameter(nexp=301)

      real,dimension(nexp) :: expname

      integer, parameter :: ndigit=3
      CHARACTER(LEN=ndigit) :: my_name
      CHARACTER(LEN=20), PARAMETER :: FMT1 = "(I3.3)"

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
     &                   'NOCVS_peltint_ice4g_orog_new_'
      character in_orog_varname*200
      parameter(in_orog_varname=
     &  'orog')

      character(len=200),parameter :: 
     &  in_icefrac_filename =
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_peltint_ice4g_icemask_new_'
      character in_icefrac_varname*200
      parameter(in_icefrac_varname=
     &  'icemask')

      character out_data_filename*200
      parameter(out_data_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'peltier_30k_new_goldstein.nc')

      character out_ascii_orogname*200
      parameter(out_ascii_orogname=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'orography_30k_new.dat')

      character out_ascii_icefracname*200
      parameter(out_ascii_icefracname=
     &  '/home/ggdjl/genie/genie-utilities/data/output/'//
     &                   'icemask_30k_new.dat')
