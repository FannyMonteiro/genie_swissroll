! ******************************************************************************************************************************** !
! atchem_lib.f90
! Atmosphere Chemistry
! LIBRARY MODULE
! ******************************************************************************************************************************** !


MODULE atchem_lib


  use genie_control
  use gem_util
  use gem_carbchem
  IMPLICIT NONE
  SAVE


  ! ****************************************************************************************************************************** !
  ! *** NAMELIST DEFINITIONS ***************************************************************************************************** !
  ! ****************************************************************************************************************************** !

  ! ### EDIT ADD AND/OR EXTEND NAME-LIST PARAMETER AND CONTROL OPTIONS ########################################################### !
  ! ------------------- TRACER INITIALIZATION ------------------------------------------------------------------------------------ !
  REAL,DIMENSION(n_atm)::atm_init                                ! atmosphere tracer array initial values
  NAMELIST /ini_atchem_nml/atm_init
  ! ------------------- COSMOGENIC & RADIOGENIC PRODUCTION  ---------------------------------------------------------------------- !
  real::par_atm_F14C                                             ! Global cosmogenic production rate of 14C (mol yr-1)
  NAMELIST /ini_atchem_nml/par_atm_F14C
  ! ------------------- EMISSIONS-TO-ATMOSPHERE ---------------------------------------------------------------------------------- !
  real::par_atm_wetlands_FCH4                                    ! Wetlands CH4 flux (mol yr-1) 
  real::par_atm_wetlands_FCH4_d13C                               ! Wetlands CH4 d13C (o/oo) 
  NAMELIST /ini_atchem_nml/par_atm_wetlands_FCH4,par_atm_wetlands_FCH4_d13C
  ! ------------------- RUN CONTROL ---------------------------------------------------------------------------------------------- !
  logical::ctrl_continuing                                       ! continuing run?
  NAMELIST /ini_atchem_nml/ctrl_continuing
  ! ------------------- I/O DIRECTORY DEFINITIONS -------------------------------------------------------------------------------- !
  CHARACTER(len=255)::par_indir_name                             ! 
  CHARACTER(len=255)::par_outdir_name                            ! 
  CHARACTER(len=255)::par_rstdir_name                            ! 
  NAMELIST /ini_atchem_nml/par_indir_name,par_outdir_name,par_rstdir_name
  CHARACTER(len=127)::par_infile_name,par_outfile_name           ! 
  NAMELIST /ini_atchem_nml/par_infile_name,par_outfile_name
  ! ############################################################################################################################## !


  ! ****************************************************************************************************************************** !
  ! *** MODEL CONFIGURATION CONSTANTS ******************************************************************************************** !
  ! ****************************************************************************************************************************** !

  ! *** array dimensions ***
  ! grid dimensions
  INTEGER,PARAMETER::n_i                                  = ilon1_atm ! 
  INTEGER,PARAMETER::n_j                                  = ilat1_atm ! 
  ! grid properties array dimensions 
  INTEGER,PARAMETER::n_phys_atm                           = 13    ! number of grid properties descriptors

  ! *** array index values ***
  ! atmosperhic 'physics' properties array indices
  INTEGER,PARAMETER::ipa_lat                              = 01    ! latitude (degrees) [mid-point]
  INTEGER,PARAMETER::ipa_lon                              = 02    ! longitude (degrees) [mid-point]
  INTEGER,PARAMETER::ipa_dlat                             = 03    ! latitude (degrees) [width]
  INTEGER,PARAMETER::ipa_dlon                             = 04    ! longitude (degrees) [width]
  INTEGER,PARAMETER::ipa_hmid                             = 05    ! height (m) [mid-point]
  INTEGER,PARAMETER::ipa_dh                               = 06    ! height (m) [thickness]
  INTEGER,PARAMETER::ipa_hbot                             = 07    ! height (m) [bottom]
  INTEGER,PARAMETER::ipa_htop                             = 08    ! height (m) [top]
  INTEGER,PARAMETER::ipa_A                                = 09    ! area (m2)
  INTEGER,PARAMETER::ipa_rA                               = 10    ! reciprocal area (to speed up numerics)
  INTEGER,PARAMETER::ipa_V                                = 11    ! atmospheric box volume (m3)
  INTEGER,PARAMETER::ipa_rV                               = 12    ! reciprocal volume (to speed up numerics)
  INTEGER,PARAMETER::ipa_P                                = 13    ! pressure (atm)

  ! *** array index names ***
  ! atmosphere interface 'physics'
  CHARACTER(len=16),DIMENSION(n_phys_atm),PARAMETER::string_phys_atm = (/ &
       & 'lat             ', &
       & 'lon             ', &
       & 'dlat            ', &
       & 'dlon            ', &
       & 'hmid            ', &
       & 'dh              ', &
       & 'hbot            ', &
       & 'htop            ', &
       & 'A               ', &
       & 'rA              ', &
       & 'V               ', &
       & 'rV              ', &
       & 'P               ' /)

  ! *** miscellaneous ***
  ! effective thickness of atmosphere (m) in the case of a 1-cell thick atmosphere
  ! NOTE: was 8000.0 m in Ridgwell et al. [2007]
  REAL,parameter::par_atm_th = 7777.0

  ! *********************************************************
  ! *** GLOBAL VARIABLE AND RUN-TIME SET PARAMETER ARRAYS ***
  ! *********************************************************

  ! *** PRIMARY ATCHEM ARRAYS ***
  real,dimension(n_atm,n_i,n_j)::atm           ! 
  real,dimension(n_atm,n_i,n_j)::fatm          ! 
  real,dimension(n_phys_atm,n_i,n_j)::phys_atm ! 


END MODULE atchem_lib

