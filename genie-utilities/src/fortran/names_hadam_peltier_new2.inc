
      integer nexpin
      parameter(nexpin=22)

      character in_icefraclgm_filename*200
      parameter(in_icefraclgm_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_icemask_')
      character in_icefraclgm_varname*200
      parameter(in_icefraclgm_varname=
     &  'icemask')

      character in_oroglgm_filename*200
      parameter(in_oroglgm_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'peltier_ice4g_orog_')
      character in_oroglgm_varname*200
      parameter(in_oroglgm_varname=
     &  'orog')

      integer ndigit2
      parameter(ndigit2=2)
      CHARACTER(LEN=20), PARAMETER :: FMT2 = "(I2.2)"
      CHARACTER(LEN=ndigit2) :: my_name_2

      integer nexp
      parameter(nexp=211)
      integer, parameter :: ni=0
      integer ndigit
      parameter(ndigit=3)
      CHARACTER(LEN=20), PARAMETER :: FMT1 = "(I3.3)"
      CHARACTER(LEN=ndigit) :: my_name

      character(len=200),parameter :: 
     &  out_orogtime_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_peltint_ice4g_orog_new_'
      character out_orogtime_varname*200
      parameter(out_orogtime_varname=
     &  'orog')

      character(len=200),parameter :: 
     &  out_icefractime_filename=
     &  '/home/ggdjl/genie/genie-utilities/data/input/'//
     &                   'NOCVS_peltint_ice4g_icemask_new_'
      character out_icefractime_varname*200
      parameter(out_icefractime_varname=
     &  'icemask')

      character(len=200),parameter :: 
     &  in_d18o_filename=
     &    '/home/ggdjl/genie/genie-utilities/data/input/'//
     &              'd18o_30k_highres.dat'

      integer nd18o
      parameter(nd18o=301)
