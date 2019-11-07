      character in_orog_filename*200
      parameter(in_orog_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_xakef.qrparm.orog.nc')
      character in_orog_varname*200
      parameter(in_orog_varname=
     &  'ht')

      character in_icefrac_filename*200
      parameter(in_icefrac_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_xakef.qrparm.veg.nc')
      character in_icefrac_varname*200
      parameter(in_icefrac_varname=
     &  'field322')



      character in_icefraclgm_filename*200
      parameter(in_icefraclgm_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_icemask_21.nc')
      character in_icefraclgm_varname*200
      parameter(in_icefraclgm_varname=
     &  'icemask')

      character in_oroglgm_filename*200
      parameter(in_oroglgm_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_orog_21.nc')
      character in_oroglgm_varname*200
      parameter(in_oroglgm_varname=
     &  'orog')



      character out_orog_filename*200
      parameter(out_orog_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_peltst3_ice4g_orog_30.nc')

      character out_icefrac_filename*200
      parameter(out_icefrac_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_peltst3_ice4g_icemask_30.nc')


      integer nexp
      parameter(nexp=10)
      integer, parameter :: ni=21
      integer ndigit
      parameter(ndigit=2)
      CHARACTER(LEN=20), PARAMETER :: FMT1 = "(I2.2)"

      CHARACTER(LEN=ndigit) :: my_name

      character(len=200),parameter :: 
     &  out_orogtime_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_peltlg3_ice4g_orog_'
      character out_orogtime_varname*200
      parameter(out_orogtime_varname=
     &  'orog')

      character(len=200),parameter :: 
     &  out_icefractime_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_peltlg3_ice4g_icemask_'
      character out_icefractime_varname*200
      parameter(out_icefractime_varname=
     &  'icemask')

      character(len=200),parameter :: 
     &  in_d18o_filename=
     &    '/home/ggdjl/genie/genie-utilities/data/input/'//
     &              'd18o_30k_highres.dat'

      integer nd18o
      parameter(nd18o=31)
