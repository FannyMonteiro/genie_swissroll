! ******************************************************************************************************************************** !
! gem_cmn.f90
! GEochemistry Model common parameter and variable module
! This Fortran90 module contains parameters and variables common to GENIE biogeochemistry.
! ******************************************************************************************************************************** !


MODULE gem_cmn


  IMPLICIT NONE
  SAVE


  ! ****************************************************************************************************************************** !
  ! MODEL CONFIGURATION CONSTANTS - TRACER ARRAY DIMENSIONS
  ! ****************************************************************************************************************************** !


  ! *** tracer array dimensions ***
  ! NOTE: these definitions must come FIRST, because the namelist arrays are dimensioned by these parameters ...
  INTEGER,PARAMETER::n_atm = 19
  INTEGER,PARAMETER::n_ocn = 77
  INTEGER,PARAMETER::n_sed = 56


  ! ****************************************************************************************************************************** !
  ! *** NAMELIST DEFINITIONS ***************************************************************************************************** !
  ! ****************************************************************************************************************************** !


  ! #### EDIT ADD AND/OR EXTEND NAME-LIST PARAMETER AND CONTROL OPTIONS ########################################################## !
  ! ------------------- TRACER SELECTION ----------------------------------------------------------------------------------------- !
  LOGICAL,DIMENSION(n_atm)::atm_select                                  !
  LOGICAL,DIMENSION(n_ocn)::ocn_select                                  !
  LOGICAL,DIMENSION(n_sed)::sed_select                                  !
  NAMELIST /ini_gem_nml/atm_select,ocn_select,sed_select
  ! ------------------- MISC CONTROLS -------------------------------------------------------------------------------------------- !
  real::par_grid_lon_offset                                             ! assumed lon grid offset (w.r.t. Prime Meridian)
  NAMELIST /ini_gem_nml/par_grid_lon_offset
  CHARACTER(len=63)::par_carbconstset_name                              ! carbonate dissociation constants set
  NAMELIST /ini_gem_nml/par_carbconstset_name
  real::par_carbchem_pH_tolerance                                       ! pH solution tolerance
  integer::par_carbchem_pH_iterationmax                                 ! pH solution maximum number of iterations
  NAMELIST /ini_gem_nml/par_carbchem_pH_tolerance,par_carbchem_pH_iterationmax
  logical::ctrl_carbchem_fail                                           ! Exit upon pH solution failure?
  NAMELIST /ini_gem_nml/ctrl_carbchem_fail
  integer::ctrl_debug_init,ctrl_debug_loop,ctrl_debug_end               ! 
  NAMELIST /ini_gem_nml/ctrl_debug_init,ctrl_debug_loop,ctrl_debug_end
  ! ------------------- I/O: DIRECTORY DEFINITIONS ------------------------------------------------------------------------------- !
  CHARACTER(len=127)::par_gem_indir_name                                ! 
  NAMELIST /ini_gem_nml/par_gem_indir_name
  CHARACTER(len=04)::string_results_ext                                 ! 
  NAMELIST /ini_gem_nml/string_results_ext
  ! ############################################################################################################################## !


  ! ****************************************************************************************************************************** !
  ! *** MODEL CONFIGURATION CONSTANTS ******************************************************************************************** !
  ! ****************************************************************************************************************************** !

  ! *** array dimensions ***
  ! main biogeochem ocean array dimensions 
  INTEGER,PARAMETER::n_carb                               = 10          ! number of ocean box chemistry descriptors
  INTEGER,PARAMETER::n_carbconst                          = 17          ! number of ocean box chemistry constants descriptors
  INTEGER,PARAMETER::n_carbalk                            = 13          ! number of alkalinty chemistry descriptors
  INTEGER,PARAMETER::n_carbisor                           = 08          ! number of carbonate isotopic ratio descriptors
  ! *** array index values ***
  ! ocean tracer IDs
  INTEGER,PARAMETER::io_T                                 = 01    ! 
  INTEGER,PARAMETER::io_S                                 = 02    ! 
  INTEGER,PARAMETER::io_DIC                               = 03    ! 
  INTEGER,PARAMETER::io_DIC_13C                           = 04    ! 
  INTEGER,PARAMETER::io_DIC_14C                           = 05    ! 
  INTEGER,PARAMETER::io_NO3                               = 06    ! 
  INTEGER,PARAMETER::io_NO3_15N                           = 07    ! 
  INTEGER,PARAMETER::io_PO4                               = 08    ! 
  INTEGER,PARAMETER::io_Fe                                = 09    ! 
  INTEGER,PARAMETER::io_O2                                = 10    ! 
  INTEGER,PARAMETER::io_O2_18O                            = 11    ! 
  INTEGER,PARAMETER::io_ALK                               = 12    ! 
  INTEGER,PARAMETER::io_SiO2                              = 13    ! 
  INTEGER,PARAMETER::io_SiO2_30Si                         = 14    ! 
  INTEGER,PARAMETER::io_DOM_C                             = 15    ! 
  INTEGER,PARAMETER::io_DOM_C_13C                         = 16    ! 
  INTEGER,PARAMETER::io_DOM_C_14C                         = 17    ! 
  INTEGER,PARAMETER::io_DOM_N                             = 18    ! 
  INTEGER,PARAMETER::io_DOM_N_15N                         = 19    ! 
  INTEGER,PARAMETER::io_DOM_P                             = 20    ! 
  INTEGER,PARAMETER::io_DOM_Cd                            = 21    ! 
  INTEGER,PARAMETER::io_DOM_Cd_114Cd                      = 52    ! 
  INTEGER,PARAMETER::io_DOM_Fe                            = 22    ! 
  INTEGER,PARAMETER::io_RDOM_C                            = 67    ! 
  INTEGER,PARAMETER::io_RDOM_C_13C                        = 68    ! 
  INTEGER,PARAMETER::io_RDOM_C_14C                        = 69    ! 
  INTEGER,PARAMETER::io_RDOM_N                            = 70    ! 
  INTEGER,PARAMETER::io_RDOM_N_15N                        = 71    ! 
  INTEGER,PARAMETER::io_RDOM_P                            = 72    ! 
  INTEGER,PARAMETER::io_RDOM_Cd                           = 73    ! 
  INTEGER,PARAMETER::io_RDOM_Cd_114Cd                     = 74    ! 
  INTEGER,PARAMETER::io_RDOM_Fe                           = 75    ! 
  INTEGER,PARAMETER::io_FeL                               = 23    ! 
  INTEGER,PARAMETER::io_L                                 = 24    ! 
  INTEGER,PARAMETER::io_CH4                               = 25    ! 
  INTEGER,PARAMETER::io_CH4_13C                           = 26    ! 
  INTEGER,PARAMETER::io_CH4_14C                           = 27    ! 
  INTEGER,PARAMETER::io_NH4                               = 28    ! 
  INTEGER,PARAMETER::io_NH4_15N                           = 29    ! 
  INTEGER,PARAMETER::io_N2                                = 30    ! 
  INTEGER,PARAMETER::io_N2_15N                            = 31    !
  INTEGER,PARAMETER::io_N2O                               = 32    ! 
  INTEGER,PARAMETER::io_N2O_15N                           = 33    !
  INTEGER,PARAMETER::io_Cd                                = 34    !
  INTEGER,PARAMETER::io_Cd_114Cd                          = 51    ! 
  INTEGER,PARAMETER::io_Ca                                = 35    ! 
  INTEGER,PARAMETER::io_Mg                                = 50    ! 
  INTEGER,PARAMETER::io_Li                                = 53    ! 
  INTEGER,PARAMETER::io_Li_7Li                            = 54    ! 
  INTEGER,PARAMETER::io_Nd                                = 55    ! 
  INTEGER,PARAMETER::io_Nd_144Nd                          = 56    ! 
  INTEGER,PARAMETER::io_B                                 = 36    ! 
  INTEGER,PARAMETER::io_F                                 = 37    ! 
  INTEGER,PARAMETER::io_SO4                               = 38    ! 
  INTEGER,PARAMETER::io_SO4_34S                           = 39    ! 
  INTEGER,PARAMETER::io_H2S                               = 40    !  
  INTEGER,PARAMETER::io_H2S_34S                           = 41    ! 
  INTEGER,PARAMETER::io_Ge                                = 42    !  
  INTEGER,PARAMETER::io_231Pa                             = 43    !  
  INTEGER,PARAMETER::io_230Th                             = 44    !  
  INTEGER,PARAMETER::io_CFC11                             = 45    ! 
  INTEGER,PARAMETER::io_CFC12                             = 46    ! 
  INTEGER,PARAMETER::io_SF6                               = 47    !  
  INTEGER,PARAMETER::io_colr                              = 48          ! RED numerical (color) tracer
  INTEGER,PARAMETER::io_colb                              = 49          ! BLUE numerical (color) tracer
  INTEGER,PARAMETER::io_col0                              = 57          ! numerical (color) tracer #0
  INTEGER,PARAMETER::io_col1                              = 58          ! numerical (color) tracer #1
  INTEGER,PARAMETER::io_col2                              = 59          ! numerical (color) tracer #2
  INTEGER,PARAMETER::io_col3                              = 60          ! numerical (color) tracer #3
  INTEGER,PARAMETER::io_col4                              = 61          ! numerical (color) tracer #4
  INTEGER,PARAMETER::io_col5                              = 62          ! numerical (color) tracer #5
  INTEGER,PARAMETER::io_col6                              = 63          ! numerical (color) tracer #6
  INTEGER,PARAMETER::io_col7                              = 64          ! numerical (color) tracer #7
  INTEGER,PARAMETER::io_col8                              = 65          ! numerical (color) tracer #8
  INTEGER,PARAMETER::io_col9                              = 66          ! numerical (color) tracer #9
  ! New NH4 tracer to estimate new production supported by NH4 - Fanny Jun 2015
  INTEGER,PARAMETER::io_NH4_new                           = 76    !
  ! New NO3 tracer to estimate new production supported by NO3 - Fanny Sep 2015
  INTEGER,PARAMETER::io_NO3_new                           = 77    !
  ! atmospheric tracer indices
  INTEGER,PARAMETER::ia_T                                 = 01    ! temperature
  INTEGER,PARAMETER::ia_q                                 = 02    ! specific humidity
  INTEGER,PARAMETER::ia_pCO2                              = 03    ! pCO2
  INTEGER,PARAMETER::ia_pCO2_13C                          = 04    ! 13C (pCO2)
  INTEGER,PARAMETER::ia_pCO2_14C                          = 05    ! 14C (pCO2)
  INTEGER,PARAMETER::ia_pO2                               = 06    ! pO2
  INTEGER,PARAMETER::ia_pO2_18O                           = 07    ! 18O (pO2)
  INTEGER,PARAMETER::ia_pN2                               = 08    ! pN2
  INTEGER,PARAMETER::ia_pN2_15N                           = 09    ! 15N (pN2)
  INTEGER,PARAMETER::ia_pCH4                              = 10    ! pCH4
  INTEGER,PARAMETER::ia_pCH4_13C                          = 11    ! 13C (pCH4)
  INTEGER,PARAMETER::ia_pCH4_14C                          = 12    ! 14C (pCH4)
  INTEGER,PARAMETER::ia_pSF6                              = 13    ! halo-carbon
  INTEGER,PARAMETER::ia_pN2O                              = 14    ! pN2
  INTEGER,PARAMETER::ia_pN2O_15N                          = 15    ! 15N (pN2)
  INTEGER,PARAMETER::ia_pH2S                              = 16    ! pH2S
  INTEGER,PARAMETER::ia_pH2S_34S                          = 17    ! pH2S
  INTEGER,PARAMETER::ia_pCFC11                            = 18    ! halo-carbon
  INTEGER,PARAMETER::ia_pCFC12                            = 19    ! halo-carbon
  ! sediment tracer indices
  INTEGER,PARAMETER::is_NULL1                             = 01    ! 
  INTEGER,PARAMETER::is_NULL2                             = 02    ! 
  INTEGER,PARAMETER::is_POC                               = 03    ! 
  INTEGER,PARAMETER::is_POC_13C                           = 04    ! 
  INTEGER,PARAMETER::is_POC_14C                           = 05    ! 
  INTEGER,PARAMETER::is_PON                               = 06    ! 
  INTEGER,PARAMETER::is_PON_15N                           = 07    ! 
  INTEGER,PARAMETER::is_POP                               = 08    ! 
  INTEGER,PARAMETER::is_POCd                              = 09    ! 
  INTEGER,PARAMETER::is_POCd_114Cd                        = 43    ! 
  INTEGER,PARAMETER::is_POFe                              = 10    ! 
  INTEGER,PARAMETER::is_POM_231Pa                         = 11    ! 
  INTEGER,PARAMETER::is_POM_230Th                         = 12    ! 
  INTEGER,PARAMETER::is_POM_Fe                            = 13    !  
  INTEGER,PARAMETER::is_POM_Nd                            = 47    ! 
  INTEGER,PARAMETER::is_POM_Nd_144Nd                      = 48    ! 
  INTEGER,PARAMETER::is_CaCO3                             = 14    ! 
  INTEGER,PARAMETER::is_CaCO3_13C                         = 15    ! 
  INTEGER,PARAMETER::is_CaCO3_14C                         = 16    ! 
  INTEGER,PARAMETER::is_CaCO3_18O                         = 17    !
  INTEGER,PARAMETER::is_CdCO3                             = 18    ! 
  INTEGER,PARAMETER::is_CdCO3_114Cd                       = 44    !
  INTEGER,PARAMETER::is_LiCO3                             = 45    ! 
  INTEGER,PARAMETER::is_LiCO3_7Li                         = 46    !
  INTEGER,PARAMETER::is_CaCO3_231Pa                       = 19    ! 
  INTEGER,PARAMETER::is_CaCO3_230Th                       = 20    !  
  INTEGER,PARAMETER::is_CaCO3_Fe                          = 21    !  
  INTEGER,PARAMETER::is_CaCO3_Nd                          = 49    ! 
  INTEGER,PARAMETER::is_CaCO3_Nd_144Nd                    = 50    !   
  INTEGER,PARAMETER::is_det                               = 22    ! 
  INTEGER,PARAMETER::is_det_231Pa                         = 23    ! 
  INTEGER,PARAMETER::is_det_230Th                         = 24    ! 
  INTEGER,PARAMETER::is_det_Fe                            = 25    ! 
  INTEGER,PARAMETER::is_det_Nd                            = 51    ! 
  INTEGER,PARAMETER::is_det_Nd_144Nd                      = 52    ! 
  INTEGER,PARAMETER::is_detLi                             = 55    ! 
  INTEGER,PARAMETER::is_detLi_7Li                         = 56    !     
  INTEGER,PARAMETER::is_opal                              = 26    ! 
  INTEGER,PARAMETER::is_opal_30Si                         = 27    ! 
  INTEGER,PARAMETER::is_opalGe                            = 28    ! 
  INTEGER,PARAMETER::is_opal_231Pa                        = 29    ! 
  INTEGER,PARAMETER::is_opal_230Th                        = 30    ! 
  INTEGER,PARAMETER::is_opal_Fe                           = 31    ! 
  INTEGER,PARAMETER::is_opal_Nd                           = 53    ! 
  INTEGER,PARAMETER::is_opal_Nd_144Nd                     = 54    ! 
  INTEGER,PARAMETER::is_ash                               = 32    ! 
  INTEGER,PARAMETER::is_POC_frac2                         = 33    ! 
  INTEGER,PARAMETER::is_CaCO3_frac2                       = 34    ! 
  INTEGER,PARAMETER::is_opal_frac2                        = 35    ! 
  INTEGER,PARAMETER::is_CaCO3_age                         = 36    ! 
  INTEGER,PARAMETER::is_foram_p_13C                       = 37    ! 
  INTEGER,PARAMETER::is_foram_p_14C                       = 38    ! 
  INTEGER,PARAMETER::is_foram_p_18O                       = 39    ! 
  INTEGER,PARAMETER::is_foram_b_13C                       = 40    ! 
  INTEGER,PARAMETER::is_foram_b_14C                       = 41    ! 
  INTEGER,PARAMETER::is_foram_b_18O                       = 42    ! 
  ! (carbonate) chemistry descriptors array indices
  INTEGER,PARAMETER::ic_H                                 = 01    ! H+ concentration
  INTEGER,PARAMETER::ic_fug_CO2                           = 02    ! CO2 fugacity
  INTEGER,PARAMETER::ic_conc_CO2                          = 03    ! CO2(aq) concentration
  INTEGER,PARAMETER::ic_conc_CO3                          = 04    ! CO32- concentration
  INTEGER,PARAMETER::ic_conc_HCO3                         = 05    ! HCO3- concentration
  INTEGER,PARAMETER::ic_ohm_cal                           = 06    ! ohmega(calcite)
  INTEGER,PARAMETER::ic_ohm_arg                           = 07    ! ohmega(aragonite)
  INTEGER,PARAMETER::ic_dCO3_cal                          = 08    ! degree of over-saturation [CO32-] w.r.t. calcite
  INTEGER,PARAMETER::ic_dCO3_arg                          = 09    ! degree of over-saturation [CO32-] w.r.t. aragonite
  INTEGER,PARAMETER::ic_RF0                               = 10    ! Revelle factor
  ! (carbonate) chemistry descriptors array indices
  INTEGER,PARAMETER::icc_k                                = 01    ! 
  INTEGER,PARAMETER::icc_k1                               = 02    ! 
  INTEGER,PARAMETER::icc_k2                               = 03    ! 
  INTEGER,PARAMETER::icc_kB                               = 04    ! 
  INTEGER,PARAMETER::icc_kW                               = 05    ! 
  INTEGER,PARAMETER::icc_kSi                              = 06    ! 
  INTEGER,PARAMETER::icc_kHF                              = 07    ! 
  INTEGER,PARAMETER::icc_kHSO4                            = 08    ! 
  INTEGER,PARAMETER::icc_kP1                              = 09    ! 
  INTEGER,PARAMETER::icc_kP2                              = 10    ! 
  INTEGER,PARAMETER::icc_kP3                              = 11    ! 
  INTEGER,PARAMETER::icc_kcal                             = 12    ! 
  INTEGER,PARAMETER::icc_karg                             = 13    ! 
  INTEGER,PARAMETER::icc_QCO2                             = 14    ! 
  INTEGER,PARAMETER::icc_QO2                              = 15    ! 
  INTEGER,PARAMETER::icc_kH2S                             = 16    ! 
  INTEGER,PARAMETER::icc_kNH4                             = 17    ! 
  ! (carbonate) alkalinty chemistry descriptor indices
  INTEGER,PARAMETER::ica_HCO3                             = 01    !
  INTEGER,PARAMETER::ica_CO3                              = 02    !
  INTEGER,PARAMETER::ica_H4BO4                            = 03    !
  INTEGER,PARAMETER::ica_OH                               = 04    !
  INTEGER,PARAMETER::ica_HPO4                             = 05    !
  INTEGER,PARAMETER::ica_PO4                              = 06    !
  INTEGER,PARAMETER::ica_H3SiO4                           = 07    !
  INTEGER,PARAMETER::ica_NH3                              = 08    !
  INTEGER,PARAMETER::ica_HS                               = 09    !
  INTEGER,PARAMETER::ica_H                                = 10    !
  INTEGER,PARAMETER::ica_HSO4                             = 11    !
  INTEGER,PARAMETER::ica_HF                               = 12    !
  INTEGER,PARAMETER::ica_H3PO4                            = 13    !
  ! (carbonate) isotope descriptor indices
  INTEGER,PARAMETER::ici_DIC_r13C                         = 01    !
  INTEGER,PARAMETER::ici_CO2_r13C                         = 02    !
  INTEGER,PARAMETER::ici_CO3_r13C                         = 03    !
  INTEGER,PARAMETER::ici_HCO3_r13C                        = 04    !
  INTEGER,PARAMETER::ici_DIC_r14C                         = 05    !
  INTEGER,PARAMETER::ici_CO2_r14C                         = 06    !
  INTEGER,PARAMETER::ici_CO3_r14C                         = 07    !
  INTEGER,PARAMETER::ici_HCO3_r14C                        = 08    !
  ! *** tracer descriptors ***
  ! sediment tracer types
  INTEGER,PARAMETER::par_sed_type_bio                     = 01    ! 
  INTEGER,PARAMETER::par_sed_type_abio                    = 02    ! 
  INTEGER,PARAMETER::par_sed_type_POM                     = 03    ! 
  INTEGER,PARAMETER::par_sed_type_CaCO3                   = 04    ! 
  INTEGER,PARAMETER::par_sed_type_opal                    = 05    !  
  INTEGER,PARAMETER::par_sed_type_det                     = 06    ! 
  INTEGER,PARAMETER::par_sed_type_scavenged               = 07    ! 
  INTEGER,PARAMETER::par_sed_type_age                     = 08    ! 
  INTEGER,PARAMETER::par_sed_type_frac                    = 09    ! 
!!!  INTEGER,PARAMETER::par_sed_type_indepsinking            = 21    !

  ! *** tracer arrays ***
  ! tracer description - 'type'
  integer,DIMENSION(n_atm)::atm_type
  integer,DIMENSION(n_ocn)::ocn_type
  integer,DIMENSION(n_sed)::sed_type
  ! tracer description - 'dependency'
  integer,DIMENSION(n_atm)::atm_dep
  integer,DIMENSION(n_ocn)::ocn_dep
  integer,DIMENSION(n_sed)::sed_dep
  ! tracer short names
  CHARACTER(len=16),DIMENSION(n_ocn)::string_ocn
  CHARACTER(len=16),DIMENSION(n_atm)::string_atm
  CHARACTER(len=16),DIMENSION(n_sed)::string_sed
  ! tracer long names (i.e., full description)
  CHARACTER(len=128),DIMENSION(n_ocn)::string_longname_ocn
  CHARACTER(len=128),DIMENSION(n_atm)::string_longname_atm
  CHARACTER(len=128),DIMENSION(n_sed)::string_longname_sed !
  ! tracer descriptions (for netCDF)
  CHARACTER(len=16),DIMENSION(n_atm)::string_atm_tname       ! names of active atm tracers
  CHARACTER(len=128),DIMENSION(n_atm)::string_atm_tlname     ! longnames of active atm tracers
  CHARACTER(len=12),DIMENSION(n_atm)::string_atm_unit        ! main units of active atm tracers
  REAL,DIMENSION(n_atm,2)::atm_mima                          ! atm tracer min and max (for netcdf file)
  CHARACTER(len=16),DIMENSION(n_ocn)::string_ocn_tname       ! names of active ocn tracers
  CHARACTER(len=128),DIMENSION(n_ocn)::string_ocn_tlname     ! longnames of active ocn tracers
  CHARACTER(len=12),DIMENSION(n_ocn)::string_ocn_unit        ! main units of active ocn tracers
  REAL,DIMENSION(n_ocn,2)::ocn_mima                          ! tracer min and max (for netcdf file)
  CHARACTER(len=16),DIMENSION(n_sed)::string_sed_tname       ! names of active sed tracers
  CHARACTER(len=128),DIMENSION(n_sed)::string_sed_tlname     ! longnames of active sed tracers
  CHARACTER(len=12),DIMENSION(n_sed)::string_sed_unit        ! main units of active sed tracers
  REAL,DIMENSION(n_sed,2)::sed_mima                          ! sed tracer min and max (for netcdf file)
  ! number of included (selected) tracers
  integer::n_l_atm
  integer::n_l_ocn
  integer::n_l_sed
  ! conversion of selected tracer index -> absolute index
  INTEGER,ALLOCATABLE,DIMENSION(:)::conv_iselected_ia
  INTEGER,ALLOCATABLE,DIMENSION(:)::conv_iselected_io
  INTEGER,ALLOCATABLE,DIMENSION(:)::conv_iselected_is
  ! conversion of absolute index -> selected tracer index
  INTEGER,DIMENSION(n_atm)::conv_ia_lselected
  INTEGER,DIMENSION(n_ocn)::conv_io_lselected
  INTEGER,DIMENSION(n_sed)::conv_is_lselected
  ! tracer conversion - transformation ratios
  real,DIMENSION(n_sed,n_ocn)::conv_ocn_sed
  real,DIMENSION(n_ocn,n_sed)::conv_sed_ocn
  real,DIMENSION(n_atm,n_ocn)::conv_ocn_atm
  real,DIMENSION(n_ocn,n_atm)::conv_atm_ocn
  real,DIMENSION(n_sed,n_ocn)::conv_DOM_POM
  real,DIMENSION(n_ocn,n_sed)::conv_POM_DOM
  real,DIMENSION(n_sed,n_ocn)::conv_RDOM_POM
  real,DIMENSION(n_ocn,n_sed)::conv_POM_RDOM
  ! tracer conversion - indices for non-zero transformation ratio values
  ! NOTE: the zero index place in the array is used in algorithms identifying null relationships (or something)
  integer,DIMENSION(0:n_sed,0:n_ocn)::conv_ocn_sed_i
  integer,DIMENSION(0:n_ocn,0:n_sed)::conv_sed_ocn_i
  integer,DIMENSION(0:n_atm,0:n_ocn)::conv_ocn_atm_i
  integer,DIMENSION(0:n_ocn,0:n_atm)::conv_atm_ocn_i
  integer,DIMENSION(0:n_sed,0:n_ocn)::conv_DOM_POM_i
  integer,DIMENSION(0:n_ocn,0:n_sed)::conv_POM_DOM_i
  integer,DIMENSION(0:n_sed,0:n_ocn)::conv_RDOM_POM_i
  integer,DIMENSION(0:n_ocn,0:n_sed)::conv_POM_RDOM_i
  ! carbonate chemistry
  CHARACTER(len=16),DIMENSION(n_carb),PARAMETER::string_carb = (/ &
       & 'H               ', &
       & 'fug_CO2         ', &
       & 'conc_CO2        ', &
       & 'conc_CO3        ', &
       & 'conc_HCO3       ', &
       & 'ohm_cal         ', &
       & 'ohm_arg         ', &
       & 'dCO3_cal        ', &
       & 'dCO3_arg        ', &
       & 'RF0             ' /)
  ! carbonate chemistry dissociation constants
  CHARACTER(len=16),DIMENSION(n_carbconst),PARAMETER::string_carbconst = (/ &
       & 'k               ', &
       & 'k1              ', &
       & 'k2              ', &
       & 'kB              ', &
       & 'kW              ', &
       & 'kSi             ', &
       & 'kHF             ', &
       & 'kHSO4           ', &
       & 'kP1             ', &
       & 'kP2             ', &
       & 'kP3             ', &
       & 'kcal            ', &
       & 'karg            ', &
       & 'QCO2            ', &
       & 'QO2             ', &
       & 'kH2S            ', &
       & 'kNH4            ' /)
  ! alkalinty chemistry contributions
  CHARACTER(len=16),DIMENSION(n_carbalk),PARAMETER::string_carbalk = (/ &
       & 'ALK_HCO3        ', &
       & 'ALK_CO3         ', &
       & 'ALK_H4BO3       ', &
       & 'ALK_OH          ', &
       & 'ALK_HPO4        ', &
       & 'ALK_PO4         ', &
       & 'ALK_H3SiO4      ', &
       & 'ALK_NH3         ', &
       & 'ALK_HS          ', &
       & 'ALK_H           ', &
       & 'ALK_HSO4        ', &
       & 'ALK_HF          ', &
       & 'ALK_H3PO4       ' /)
  ! carbon isotopes 
  CHARACTER(len=16),DIMENSION(n_carbisor),PARAMETER::string_carbisor = (/ &
       & 'DIC_r13C        ', &
       & 'CO2_r13C        ', &
       & 'CO3_r13C        ', &
       & 'HCO3_r13C       ', &
       & 'DIC_r14C        ', &
       & 'CO2_r14C        ', &
       & 'CO3_r14C        ', &
       & 'HCO3_r14C       ' /)

  ! *** I/O ***
  ! default I/O parameters
  INTEGER,PARAMETER::in                                   = 12
  INTEGER,PARAMETER::out                                  = 13
  ! array allocation errors
  INTEGER::error,alloc_error,dealloc_error
  ! flag to allow model to end gracefully ... (i.e., with some data saving)
  logical::error_stop     = .FALSE.
  logical::error_carbchem = .FALSE.
  ! file extensions
  ! GHC 10/06/09 string_results_ext now defined in namelist to aid processing with non-matlab programs.
  !CHARACTER(len=04),PARAMETER::string_results_ext         = '.res'
  CHARACTER(len=04),PARAMETER::string_data_ext            = '.dat'

  ! *** conversion factors ***
  ! NOTE: many taken virtually unaltered from 'SUE' [Ridgwell, 2001]
  ! NOTE: conv_atm_mol is calculated from vol / (conv_Pa_atm*const_R_SI*T)
  !       -> value is taken directly from that calculated in atchem.f90 (with T = 273.15K)
  !       (1) assuming an atmospheric thickness of 8000 m => conv_atm_mol = 1.81994E+20
  !           => this gives an equivanence for modern climate of: 0.5059433E+17 mol <-> 0.2780000E-03 atm
  !              or 1 ppm CO2 = 2.1839 PgC
  !       (2) if require the CDIAC recommendation of 1 ppm = 2.13 PgC,
  !           then atmospheric volume must be adjusted by: 2.13/2.1839 so that 1 ppm does nto represent quite so much CO2 mass
  !           => atmospheric thickness = 2.13/2.1839 * 8000 = 7803 m
  !           => conv_atm_mol = 2.13/2.1839 * 1.81994E+20 = 1.7750e+020
  !       (3) if require the OCMIP recommendation of 1 ppm = 2.123 PgC,
  !           then atmospheric volume must be adjusted by: 2.123/2.1839
  !           => atmospheric thickness = 2.123/2.1839 * 8000 = 7777 m
  !           => conv_atm_mol = 2.123/2.1839 * 1.81994E+20 = 1.7692e+020
  ! primary
  REAL,PARAMETER::conv_atm_mol                            = 1.7692e+020
  REAL,PARAMETER::conv_m3_kg                              = 1027.649 ! from Winton and Sarachik [1993] @ 34.7o/oo,0'C
  REAL,PARAMETER::conv_yr_d                               = 365.25 !360.00 !365.0
  REAL,PARAMETER::conv_yr_hr                              = 24.0 * conv_yr_d
  REAL,PARAMETER::conv_yr_s                               = 3600.0 * conv_yr_hr
  REAL,PARAMETER::conv_kyr_yr                             = 1.0E+03
  REAL,PARAMETER::conv_m_cm                               = 1.0E+02
  REAL,PARAMETER::conv_m2_cm2                             = 1.0E+04
  REAL,PARAMETER::conv_m3_cm3                             = 1.0E+06
  REAL,PARAMETER::conv_m3_l                               = 1.0E+03
  REAL,PARAMETER::conv_kg_g                               = 1.0E+03
  REAL,PARAMETER::conv_kg_mg                              = 1.0E+06
  REAL,PARAMETER::conv_mol_mmol                           = 1.0E+03
  REAL,PARAMETER::conv_mol_umol                           = 1.0E+06
  REAL,PARAMETER::conv_mol_nmol                           = 1.0E+09
  REAL,PARAMETER::conv_mol_pmol                           = 1.0E+12
  REAL,PARAMETER::conv_g_mg                               = 1.0E+03
  REAL,PARAMETER::conv_atm_Pa                             = 1.01325E+05
  ! derived
  REAL,PARAMETER::conv_mol_atm                            = 1.0 / conv_atm_mol ! 6.024E-21
  REAL,PARAMETER::conv_kg_m3                              = 1.0 / conv_m3_kg
  REAL,PARAMETER::conv_d_yr                               = 1.0 / conv_yr_d
  REAL,PARAMETER::conv_hr_yr                              = 1.0 / conv_yr_hr
  REAL,PARAMETER::conv_s_yr                               = 1.0 / conv_yr_s
  REAL,PARAMETER::conv_yr_kyr                             = 1.0 / conv_kyr_yr
  REAL,PARAMETER::conv_cm_m                               = 1.0 / conv_m_cm
  REAL,PARAMETER::conv_cm2_m2                             = 1.0 / conv_m2_cm2
  REAL,PARAMETER::conv_cm3_m3                             = 1.0 / conv_m3_cm3
  REAL,PARAMETER::conv_l_m3                               = 1.0 / conv_m3_l
  REAL,PARAMETER::conv_g_kg                               = 1.0 / conv_kg_g
  REAL,PARAMETER::conv_mg_kg                              = 1.0 / conv_kg_mg
  REAL,PARAMETER::conv_mmol_mol                           = 1.0 / conv_mol_mmol
  REAL,PARAMETER::conv_umol_mol                           = 1.0 / conv_mol_umol
  REAL,PARAMETER::conv_nmol_mol                           = 1.0 / conv_mol_nmol
  REAL,PARAMETER::conv_pmol_mol                           = 1.0 / conv_mol_pmol
  REAL,PARAMETER::conv_mg_g                               = 1.0 / conv_g_mg
  REAL,PARAMETER::conv_Pa_atm                             = 1.0 / conv_atm_Pa
  ! other
  REAL,PARAMETER::conv_cm3_kg                             = conv_m3_kg / conv_m3_cm3
  REAL,PARAMETER::conv_kg_cm3                             = 1.0 / conv_cm3_kg   
  REAL,PARAMETER::conv_l_kg                               = conv_m3_kg / conv_m3_l
  REAL,PARAMETER::conv_kg_l                               = 1.0 / conv_l_kg 
  ! density of calcite
  ! NOTE: relative molecular mass of calcite is 100.0, density of pure calcite is approximately 2.7
  REAL,PARAMETER::conv_cal_cm3_g                          = 2.70
  REAL,PARAMETER::conv_cal_g_cm3                          = 1.0 / conv_cal_cm3_g
  REAL,PARAMETER::conv_cal_cm3_mol                        = conv_cal_cm3_g / 100.0
  REAL,PARAMETER::conv_cal_mol_cm3                        = 1.0 / conv_cal_cm3_mol
  REAL,PARAMETER::conv_cal_g_mol                          = conv_cal_g_cm3*conv_cal_cm3_mol
  REAL,PARAMETER::conv_cal_mol_g                          = conv_cal_mol_cm3*conv_cal_cm3_g
  ! density of opal
  ! NOTE: relative molecular mass of SiO2 is 60.0, density of opal is in the range 2.0 - 2.5
  REAL,PARAMETER::conv_opal_cm3_g                         = 2.25
  REAL,PARAMETER::conv_opal_g_cm3                         = 1.0 / conv_opal_cm3_g
  REAL,PARAMETER::conv_opal_cm3_mol                       = conv_opal_cm3_g / 60.0
  REAL,PARAMETER::conv_opal_mol_cm3                       = 1.0 / conv_opal_cm3_mol
  REAL,PARAMETER::conv_opal_g_mol                         = conv_opal_g_cm3*conv_opal_cm3_mol
  REAL,PARAMETER::conv_opal_mol_g                         = conv_opal_mol_cm3*conv_opal_cm3_g
  ! density of refractory material
  ! NOTE: assume average density of refractory material (SiO2) as 3.0 (g cm-3)
  REAL,PARAMETER::conv_det_cm3_g                          = 3.00
  REAL,PARAMETER::conv_det_g_cm3                          = 1.0 / conv_det_cm3_g
  REAL,PARAMETER::conv_det_cm3_mol                        = conv_det_cm3_g / 60.0
  REAL,PARAMETER::conv_det_mol_cm3                        = 1.0 / conv_det_cm3_mol
  REAL,PARAMETER::conv_det_g_mol                          = conv_det_g_cm3*conv_det_cm3_mol
  REAL,PARAMETER::conv_det_mol_g                          = conv_det_mol_cm3*conv_det_cm3_g
  ! density of organic matter
  ! NOTE: assume average density of particulate organic material as 1.0 (g cm-3)
  REAL,PARAMETER::conv_POC_cm3_g                          = 1.00
  REAL,PARAMETER::conv_POC_g_cm3                          = 1.0 / conv_POC_cm3_g    
  REAL,PARAMETER::conv_POC_cm3_mol                        = conv_POC_cm3_g / 12.0
  REAL,PARAMETER::conv_POC_mol_cm3                        = 1.0 / conv_POC_cm3_mol
  REAL,PARAMETER::conv_POC_g_mol                          = conv_POC_g_cm3*conv_POC_cm3_mol
  REAL,PARAMETER::conv_POC_mol_g                          = conv_POC_mol_cm3*conv_POC_cm3_g
  ! moles of carbon per kg
  REAL,PARAMETER::conv_C_kg_mol                           = 83.33
  REAL,PARAMETER::conv_C_mol_kg                           = 1.0 / conv_C_kg_mol
  REAL,PARAMETER::conv_CaCO3_mol_kgC                      = conv_g_kg * 12.0 
  ! moles of Fe per kg
  REAL,PARAMETER::conv_Fe_kg_mol                          = 17.86
  REAL,PARAMETER::conv_Fe_g_mol                           = conv_g_kg * conv_Fe_kg_mol
  REAL,PARAMETER::conv_Fe_mol_kg                          = 1.0 / conv_Fe_kg_mol
  ! moles of SiO2 per kg
  REAL,PARAMETER::conv_SiO2_kg_mol                        = 16.67
  REAL,PARAMETER::conv_SiO2_g_mol                         = conv_g_kg * conv_SiO2_kg_mol
  REAL,PARAMETER::conv_SiO2_mol_kg                        = 1.0 / conv_SiO2_kg_mol

  ! *** isotopes and fractionation ***
  ! 18O:(18O+17O+16O) [estimated from % natural abundance data]
  REAL,PARAMETER::const_stnd_18O_O          = 0.002004
  ! isotopic standard array
  ! NOTE: these ratios must calculated be as a fraction of the abundance of the lighter (stable) member;
  !       13R(C) = 13C / 12C
  !       to be consistent with the form of the equation used to calculate the delta notation values (fun_calc_isotope_delta)
  ! NOTE: these array positions are hard-wired in and must match the tracer config files
  ! NOTE: natural abundances [Zeebe and Wlolf-Gladrow, 2001] are;
  !       12C: 98.8922%
  !       13C:  1.1078%
  !       (14C standard is 0.95*AOx (oxalic acid standard) = 0.95*1.176E-12 = 1.117E-12)
  !       16O: 99.7628%
  !       17O:  0.0372%
  !       18O:  0.20004%
  ! NOTE: 30R(Si) = 30Si/28Si
  ! NOTE: natural abundances and NBS 28 standard (Coplen et al., 2002, Pure Appl. Chem.)
  !       28Si: 92.2223%
  !       29Si:  4.6853%
  !       30Si:  3.0924%
  !       Standards: NBS 28 silica sand defined with delta30Si=0 o/oo
  !                  30R_(NBS 28) = n_30Si/n_28Si = 0.030924/0.922223 = 0.033532
  !       6Li:   7.5%
  !       7Li:  92.5%
  !          => 7R(Li) = 92.5/7.5 = 12.33333
  !       144Nd: 
  !       114Cd: 114Cd/110Cd = 0.2873/0.1249 = 2.3002
  REAL,PARAMETER,DIMENSION(11:19)::const_standards = (/ &
       & 0.011202,  & ! TYPE 11; 13C ! OLD: 0.011237 ! NEW: 0.011202
       & 1.176E-12, & ! TYPE 12; 14C ! OLD: 1.117E-12 ! NEW: 1.176E-12
       & 0.002005,  & ! TYPE 13; 18O
       & 0.003660,  & ! TYPE 14; 15N
       & 0.000000,  & ! TYPE 15; 34S
       & 0.033532,  & ! TYPE 16; 30Si
       & 2.3002,    & ! TYPE 17; 114Cd
       & 12.33333,  & ! TYPE 18; 7Li
       & 0.000000 /)  ! TYPE 19; 144Nd

  ! #### THIS COULD BE RATIONALIZED MUCH SIMPLER ################################################################################# !
  ! *** radioactive decay ***
  ! ln(2) is used for the conversion between half-life and e-folding time of decay
  REAL,PARAMETER::const_ln2                               = 6.9314718055994529e-1 ! log(2.0)
  ! radiocarbon (14C):
  ! NOTE: half-life = 5730.0 years [Orr, 2002] (GOSAC final report)
  REAL,PARAMETER::const_lambda_14C           = 1./8267.0 ! decay constant (years^{-1})
  ! lambda for 14C (years^{-1}) using Libby half-life
  ! NOTE: half-life = 5568.0 years [Stuiver and Polach, 1977]
  REAL,PARAMETER::const_lambda_14C_libby     = 1./8033.0 ! decay constant (years^{-1})
  ! NOTE: the constants const_lamda_14C, const_lamda_14C_libby, and const_fracdecay_14C are retained here for
  !       compatibility reasons
  ! lamda for 14C (yrs)
  ! NOTE: half-life = 5730.0 years [Orr, 2002] (GOSAC final report)
  REAL,PARAMETER::const_lamda_14C            = 8267.0    ! e-folding time of radiocarbon decay (years)
  ! lambda for 14C (yrs) using Libby half-life
  ! NOTE: half-life = 5568.0 years [Stuiver and Polach, 1977]
  REAL,PARAMETER::const_lamda_14C_libby      = 8033.0
  ! yearly fractional reduction factor for 14C ( = EXP[-1.0 / const_lamda_14C] )
  REAL,PARAMETER::const_fracdecay_14C        = 0.9998790
  ! 230Th and 231Pa:
  ! lambda for 230Th (years^{-1}): half-life of 75.2e3 years used by Marchal et al. (2000)
  !                                and Siddall et al. (2005)
  REAL,PARAMETER::const_lambda_230Th         = const_ln2/75200.
  ! lambda for 231Pa (years^{-1}): half-life of 32.5e3 years used by Marchal et al. (2000)
  !                                and Siddall et al. (2005)
  REAL,PARAMETER::const_lambda_231Pa         = const_ln2/32500.
  ! lambda for 234U (years^{-1}): half-life of 245250. (+/- 490) years (Chen et al., 2000)
  ! lambda for 235U (years^{-1}): 9.8485e-10 (Steiger and Jaeger, 1977)
  REAL,PARAMETER::const_lambda_234U          = const_ln2/245250.
  REAL,PARAMETER::const_lambda_235U          = 9.8485e-10
  ! Refernces:
  ! H. Chen, R.L. Edwards, J. Hoff, C.D. Gallup, D.A. Richards, Y. Asmerom (2000), "The half-lives of uranium-234 and thorium-230", Chemical Geology, 169, pp. 17-33
  ! O. Marchal, R. Francois, T.F. Stocker, F. Joos (2000), "Ocean thermohaline circulation and sedimentary 231Pa/230Th ratio", Paleoceanography, 15, pp. 625-641
  ! M. Siddall, G.M. Henderson, N.R. Edwards, M. Frank, S.A. Mueller, T.F. Stocker, F. Joos (2005), "231Pa/230Th fractionation by ocean transport, biogenic particle flux and particle type", Earth Planet. Sci. Lett., 237, pp. 135-155
  ! R.H. Steiger, E. Jaeger (1977), "Subcomission on geochronology: Convention on the use of decay constants in geo- and cosmochronology", Earth Planet. Sci. Lett., 36, pp. 359-362
  ! const_lambda_atm: decay constants for atmospheric tracer (year^{-1})
  REAL,PARAMETER,DIMENSION(n_atm)::const_lambda_atm = (/ &
       & 0.000000, &            ! surface air temperature
       & 0.000000, &            ! specific humidity
       & 0.000000, &            ! CO2
       & 0.000000, &            ! d13C of CO2
       & const_lambda_14C, &    ! d14C of CO2
       & 0.000000, &             ! O2
       & 0.000000, & ! d18O of O2
       & 0.000000, & ! N2
       & 0.000000, & ! d15N of N2
       & 0.000000, & ! CH4
       & 0.000000, & ! d13C of CH4
       & const_lambda_14C, &    ! d14C of CH4
       & 0.000000, & ! SF6
       & 0.000000, & ! N2O
       & 0.000000, & ! d15N of N2O
       & 0.000000, & ! H2S
       & 0.000000, & ! d32S of H2S
       & 0.000000, & ! CFC-11
       & 0.000000 /) ! CFC-12
  ! const_lambda_ocn: decay constants for oceanic tracer (year^{-1})
  REAL,PARAMETER,DIMENSION(n_ocn)::const_lambda_ocn = (/ &
       & 0.000000, & ! temperature
       & 0.000000, & ! salinity
       & 0.000000, & ! DIC
       & 0.000000, & ! d13C of DIC
       & const_lambda_14C, &    ! d14C of DIC
       & 0.000000, & ! NO3
       & 0.000000, & ! d15N of NO3
       & 0.000000, & ! PO4
       & 0.000000, & ! Fe
       & 0.000000, & ! O2
       & 0.000000, & ! d18O of O2
       & 0.000000, & ! ALK
       & 0.000000, & ! H4SiO4
       & 0.000000, & ! d30Si of H4SiO4
       & 0.000000, & ! DOM-C
       & 0.000000, & ! d13C of DOM-C
       & const_lambda_14C, &    ! d14C of DOM-C
       & 0.000000, & ! DOM-N
       & 0.000000, & ! d15N of DOM-N
       & 0.000000, & ! DOM-P
       & 0.000000, & ! DOM-Cd
       & 0.000000, & ! DOM-Fe
       & 0.000000, & ! ligand-bound Fe
       & 0.000000, & ! free ligand
       & 0.000000, & ! CH4
       & 0.000000, & ! d13C of CH4
       & const_lambda_14C, &    ! d14C of CH4
       & 0.000000, & ! NH4
       & 0.000000, & ! d15N of NH3
       & 0.000000, & ! N2
       & 0.000000, & ! d15N of N2
       & 0.000000, & ! N2O
       & 0.000000, & ! d15N of N2O
       & 0.000000, & ! Cd
       & 0.000000, & ! Ca
       & 0.000000, & ! B
       & 0.000000, & ! F
       & 0.000000, & ! SO4
       & 0.000000, & ! 32S of SO4
       & 0.000000, & ! H2S
       & 0.000000, & ! 32S of H2S
       & 0.000000, & ! Ge
       & const_lambda_231Pa, & ! 231Pa
       & const_lambda_230Th, & ! 230Th
       & 0.000000, & ! CFC-11
       & 0.000000, & ! CFC-12
       & 0.000000, & ! SF6
       & 0.000000, & ! RED color tracer
       & 0.000000, & ! BLUE color tracer
       & 0.000000, & ! Mg
       & 0.000000, & ! d114Cd of Cd
       & 0.000000, & ! d114Cd of DOM-Cd
       & 0.000000, & ! Li
       & 0.000000, & ! 7Li of Li
       & 0.000000, & ! Nd
       & 0.000000, & ! d144Nd
       & 0.0, & ! color #0
       & 0.0, & ! color #1
       & 0.0, & ! color #2
       & 0.0, & ! color #3
       & 0.0, & ! color #4
       & 0.0, & ! color #5
       & 0.0, & ! color #6
       & 0.0, & ! color #7
       & 0.0, & ! color #8
       & 0.0, & ! color #9
       & 0.000000, & ! RDOM-C
       & 0.000000, & ! d13C of RDOM-C
       & const_lambda_14C, &    ! d14C of RDOM-C
       & 0.000000, & ! RDOM-N
       & 0.000000, & ! d15N of RDOM-N
       & 0.000000, & ! RDOM-P
       & 0.000000, & ! RDOM-Cd
       & 0.000000, & ! d114Cd of RDOM-Cd
       & 0.000000, & ! RDOM-Fe
       & 0.000000, & ! NH4_new - Fanny (June 2015)
       & 0.000000 /) ! NO3_new - Fanny (Sep 2015)
  ! const_lambda_sed: decay constants for sediment tracer (year^{-1})
  REAL,PARAMETER,DIMENSION(n_sed)::const_lambda_sed = (/ &
       & 0.000000, & ! dummy
       & 0.000000, & ! dummy
       & 0.000000, & ! POC
       & 0.000000, & ! d13C of POM-C
       & const_lambda_14C, &    ! d14C of POM-C
       & 0.000000, & ! POM-N
       & 0.000000, & ! d15N of POM-N
       & 0.000000, & ! POM-P
       & 0.000000, & ! POM-Cd
       & 0.000000, & ! POM-Fe
       & const_lambda_231Pa, & ! POM scavenged 231Pa
       & const_lambda_230Th, & ! POM scavenged 230Th
       & 0.000000, & ! POM scavenged Fe
       & 0.000000, & ! CaCO3
       & 0.000000, & ! d13C of CaCO3
       & const_lambda_14C, &    ! d14C of CaCO3
       & 0.000000, & ! d18O of CaCO3
       & 0.000000, & ! CaCO3 incorporated Cd
       & const_lambda_231Pa, & ! CaCO3 scavenged 231Pa
       & const_lambda_230Th, & ! CaCO3 scavenged 230Th
       & 0.000000, & ! CaCO3 scavenged Fe
       & 0.000000, & ! detrital (refractory) material
       & const_lambda_231Pa, & ! detrital scavenged 231Pa
       & const_lambda_230Th, & ! detrital scavenged 230Th
       & 0.000000, & ! detrital scavenged Fe
       & 0.000000, & ! opal
       & 0.000000, & ! d30Si of opal
       & 0.000000, & ! opal incorporated Ge
       & const_lambda_231Pa, & ! opal scavenged 231Pa
       & const_lambda_230Th, & ! opal scavenged 230Th
       & 0.000000, & ! opal scavenged Fe
       & 0.000000, & ! ash
       & 0.000000, & ! n/a
       & 0.000000, & ! n/a
       & 0.000000, & ! n/a
       & 0.000000, & ! CaCO3 numerical age tracer
       & 0.000000, & ! planktic foraminiferal CaCO3 d13C
       & const_lambda_14C, &    ! planktic foraminiferal CaCO3 d14C
       & 0.000000, & ! planktic foraminiferal CaCO3 d18O
       & 0.000000, & ! benthic foraminiferal CaCO3 d13C
       & const_lambda_14C, &    ! benthic foraminiferal CaCO3 d14C
       & 0.0, & ! benthic foraminiferal CaCO3 d18O
       & 0.0, & ! d114Cd of POCd
       & 0.0, & ! d114Cd of CdCO3
       & 0.0, & ! CaCO3 incorporated Li
       & 0.0, & ! d7Li of incorporated Li
       & 0.0, & ! POM scavenged Nd
       & 0.0, & ! d144Nd of POM scavenged Nd
       & 0.0, & ! CaCO3 scavenged Nd
       & 0.0, & ! d144Nd of CaCO3 scavenged Nd
       & 0.0, & ! det scavenged Nd
       & 0.0, & ! d144Nd of det scavenged Nd
       & 0.0, & ! opal scavenged Nd
       & 0.0, & ! d144Nd of opal scavenged Nd
       & 0.0, & ! det Li
       & 0.0 /) ! d7Li of det Li
  ! ############################################################################################################################## !

  ! *** atmospheric chemistry ***
  ! CH4 oxidation [Osborn and Wigley, 1994]
  ! NOTE: assume reference CH4 conc of 1700 ppb
  ! (not explicitly specified in Wigley [1994], but based on a 1992 calculation relative to 'present-day' CH4 ...)
  REAL,PARAMETER::const_pCH4_oxidation_tau0 = 8.3
  REAL,PARAMETER::const_pCH4_oxidation_C0   = 1700.0E-09
  REAL,PARAMETER::const_pCH4_oxidation_N    = 0.238

  ! *** miscellaneous ***
  ! PI
  REAL,PARAMETER::const_pi                                = 3.141592653589793 ! 4.0*atan(1.0) in GOLDSTEIN
  ! zero degree centigrade in Kelvin
  REAL,PARAMETER::const_zeroC                             = 273.15
  ! gas constant R (bar cm3 mol-1 K-1)
  REAL,PARAMETER::const_R                                 = 83.145
  ! gas constant R in SI units (J mol-1 K-1)
  REAL,PARAMETER::const_R_SI                              = 8.3145
  ! Radius of the Earth (m)
  REAL,PARAMETER::const_rEarth                            = 6.37E+06 ! 6.371E+06
  ! gas molar volume at STP
  REAL,PARAMETER::const_V                                 = 0.022414
  ! H2S oxidation coefficient [Zhang and Millero, 1993] (mM-2 O2 hr-1)
  REAL,PARAMETER::const_oxidation_coeff_H2S               = 1.25
  ! modern mean ocean [Mg2+] (mol kg-1)
  ! NOTE: w.r.t. S = 35 PSU [Zeebe and Wolf-Gladrow, 2001; DOE, 1994]
  REAL,PARAMETER::const_conc_Mg                           = 0.05282    
  ! modern mean ocean [Mg2+]/[Ca2+]
  ! NOTE: value is 5.14 in Tyrrell and Zeebe [2004]
  REAL,PARAMETER::const_conc_MgtoCa                       = 5.155
  ! representative ice density (kg/m**3) [from: initialise_seaice.f]
  REAL,PARAMETER::const_rho_seaice                        = 913.0
  ! Schmidt Number coefficients
  real,dimension(4,n_atm)::par_Sc_coef                                  ! 
  !  Bunsen Solubility Coefficient coefficients
  real,dimension(6,n_atm)::par_bunsen_coef                              ! 

  ! *** miscellaneous - dummy values ***
  REAL,PARAMETER::const_real_null       = -0.999999E+19                 ! 
  REAL,PARAMETER::const_real_nullhigh   = +0.999999E+19                 ! 
  REAL,PARAMETER::const_real_nullsmall  = +0.999999E-19                 ! 
  REAL,PARAMETER::const_real_zero       = +0.000000E+00                 ! 
  REAL,PARAMETER::const_real_one        = +1.000000E+00                 ! 
  integer,PARAMETER::const_integer_zero = 0                             ! 
  integer,PARAMETER::const_integer_one  = 1                             ! 


  ! ****************************************************************************************************************************** !
  ! *** END ********************************************************************************************************************** !
  ! ****************************************************************************************************************************** !

  
END MODULE gem_cmn

