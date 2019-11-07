! ******************************************************************************************************************************** !
! biogem_data.f90
! BioGeM
! DATA LOADING/SAVING ROUTINES
! ******************************************************************************************************************************** !


MODULE biogem_data


  USE biogem_lib
  USE biogem_box
  USE biogem_data_netCDF
  IMPLICIT NONE
  SAVE


CONTAINS


  ! ****************************************************************************************************************************** !
  ! DATA LOADING ROUTINES
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! LOAD BioGeM 'goin' FILE OPTIONS
  SUBROUTINE sub_load_goin_biogem()
    USE genie_util, ONLY: check_unit,check_iostat
    ! local variables
    integer::l,io,ia                                                    ! tracer counter
    integer::ios                                                        !
    ! read data_BIOGEM file
    call check_unit(in,__LINE__,__FILE__)
    open(unit=in,file='data_BIOGEM',status='old',action='read',iostat=ios)
    if (ios /= 0) then
       print*,'ERROR: could not open BIOGEM initialisation namelist file'
       stop
    end if
    ! read in namelist and close data_BIOGEM file
    read(UNIT=in,NML=ini_biogem_nml,IOSTAT=ios)
    if (ios /= 0) then
       print*,'ERROR: could not read BIOGEM namelist'
       stop
    else
       close(unit=in,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
    end if
    ! set and report namelist data
    par_indir_name = trim(par_indir_name)//'/'
    par_outdir_name = trim(par_outdir_name)//'/'
    par_rstdir_name = trim(par_rstdir_name)//'/'
    par_fordir_name = trim(par_fordir_name)//'/'
    if (ctrl_debug_init > 0) then
       ! --- TRACER INITIALIZATION  ---------------------------------------------------------------------------------------------- !
       print*,'--- INITIALIZATION ---------------------------------'
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          print*,'ocn tracer initial value: ',trim(string_ocn(io)),' = ',ocn_init(io)
          print*,'ocn tracer perturbation : ',trim(string_ocn(io)),' = ',ocn_dinit(io)
       end do
       ! --- RUN CONTROL --------------------------------------------------------------------------------------------------------- !
       print*,'--- BIOGEM TIME CONTROL ----------------------------'
       print*,'Continuing run?                                     : ',ctrl_continuing
       print*,'Simulation start year                               : ',par_misc_t_start
       print*,'Simulation run length (yr)                          : ',par_misc_t_runtime
       print*,'Time as Years Before Present?                       : ',ctrl_misc_t_BP
       ! --- MISC CONTROL -------------------------------------------------------------------------------------------------------- !
       print*,'--- MISC CONTROL -----------------------------------'
       print*,'Salanity normalization?                             : ',ctrl_misc_Snorm
       print*,'No salanity normalization?                          : ',ctrl_misc_noSnorm
       print*,'No biological update (and transformations)?         : ',ctrl_misc_nobioupdate
       print*,'Sea-ice brine rejection fraction                    : ',par_misc_brinerejection_frac
       print*,'Max j for sea-ice brine rejection                   : ',par_misc_brinerejection_jmax
       print*,'Include biogeochem in Sea-ice brine rejection?      : ',ctrl_misc_brinerejection_bgc
       ! --- BOUNDARY CONDITIONS ------------------------------------------------------------------------------------------------- !
       print*,'--- BOUNDARY CONDITIONS ----------------------------'
       print*,'Set dissolution flux = rain flux to close system?   : ',ctrl_force_sed_closedsystem
       print*,'Allow temperature / salinity forcing of climate?    : ',ctrl_force_GOLDSTEInTS
       print*,'Replace internal fractional sea-ice cover field?    : ',ctrl_force_seaice
       print*,'Replace internal wind-speed field?                  : ',ctrl_force_windspeed
       print*,'Replace internal CaCO3:POC export rain ratio?       : ',ctrl_force_CaCO3toPOCrainratio
       print*,'Replace internal POCd:POC export rain ratio?        : ',ctrl_force_POCdtoPOCrainratio
       print*,'Replace internal [Cd/P]POM/[Cd/P]SW alpha?          : ',ctrl_force_Cd_alpha
       print*,'Replace internal POC flux for 230Th/231Pa scav.     : ',ctrl_force_scav_fpart_POC
       print*,'Replace internal CaCO3 flux for 230Th/231Pa scav.   : ',ctrl_force_scav_fpart_CaCO3
       print*,'Replace internal opal flux for 230Th/231Pa scav.    : ',ctrl_force_scav_fpart_opal
       print*,'Replace internal det flux for 230Th/231Pa scav.     : ',ctrl_force_scav_fpart_det
       print*,'Value of Wanninkhof [1992] gas transfer coeff (a)   : ',par_gastransfer_a
       print*,'Filename for imposed seaice                         : ',trim(par_seaice_file)
       print*,'Filename for imposed windspeed                      : ',trim(par_windspeed_file)
       print*,'Filename for imposed CaCO3toPOCrainratio_file       : ',trim(par_CaCO3toPOCrainratio_file)
       print*,'Filename for imposed POCdtoPOCrainratio             : ',trim(par_POCdtoPOCrainratio_file)
       print*,'Filename for imposed Cd_alpha                       : ',trim(par_Cd_alpha_file)
       print*,'Filename for imposed scavenging POC flux            : ',trim(par_scav_fpart_POC_file)
       print*,'Filename for imposed scavenging CaCO3 flux          : ',trim(par_scav_fpart_CaCO3_file)
       print*,'Filename for imposed scavenging opal flux           : ',trim(par_scav_fpart_opal_file)
       print*,'Filename for imposed scavenging det flux            : ',trim(par_scav_fpart_det_file)
       print*,'Replace solar constant?                             : ',ctrl_force_solconst
       print*,'Use old tracer forcing file format?                 : ',ctrl_force_oldformat
       print*,'Forcings name                                       : ',trim(par_forcing_name)
       ! --- BIOLOGICAL NEW PRODUCTION ------------------------------------------------------------------------------------------- !
       print*,'--- BIOLOGICAL NEW PRODUCTION ----------------------'
       print*,'Biological scheme ID string                         : ',par_bio_prodopt
       print*,'Base [PO4] uptake rate (mol kg-1 yr-1)              : ',par_bio_k0_PO4
       print*,'Base [NO3] uptake rate (mol kg-1 yr-1)              : ',par_bio_k0_NO3
       print*,'[PO4] M-M half-sat value (mol kg-1)                 : ',par_bio_c0_PO4
       print*,'[PO4] M-M half-sat value [SP]                       : ',par_bio_c0_PO4_sp
       print*,'[PO4] M-M half-sat value [NSP]                      : ',par_bio_c0_PO4_nsp
       print*,'[NO3] M-M half-sat value (mol kg-1)                 : ',par_bio_c0_NO3
       print*,'[NO3]+[NH4] M-M half-sat value (mol kg-1)           : ',par_bio_c0_N
       print*,'[Fe] M-M half-sat value (mol kg-1)                  : ',par_bio_c0_Fe
       print*,'[Fe] M-M half-sat value for diazotrophs (mol kg-1)  : ',par_bio_c0_Fe_Diaz
       print*,'[Fe] M-M half-sat value [SP]                        : ',par_bio_c0_Fe_sp
       print*,'[Fe] M-M half-sat value [NSP]                       : ',par_bio_c0_Fe_nsp
       print*,'[H4SiO4] M-M half-sat value (mol kg-1)              : ',par_bio_c0_SiO2
       print*,'[H4SiO4] M-M half-sat value [SP]                    : ',par_bio_c0_SiO2_sp
       print*,'[H4SiO4] M-M half-sat value [NSP]                   : ',par_bio_c0_SiO2_nsp
       print*,'Biological production zone depth (m) (OCMIP-2)      : ',par_bio_zc
       print*,'Biological production time-scale (days) (OCMIP-2)   : ',par_bio_tau
       print*,'Biological production time-scale -- siliceous plank : ',par_bio_tau_sp
       print*,'Biological production time-scale -- non-siliceous   : ',par_bio_tau_nsp
       print*,'Fract. prod. of si. phytop. in Si/Fe-replete cond.  : ',par_bio_relprod_sp
       print*,'Light e-folding depth (m) (OCMIP-2)                 : ',par_bio_I_eL
       print*,'Coefficient for T-dep. uptake rate modifier         : ',par_bio_kT0
       print*,'e-folding temp. (K) for T-dep. uptake rate modifier : ',par_bio_kT_eT
       ! --- ORGANIC MATTER EXPORT RATIOS ---------------------------------------------------------------------------------------- !
       print*,'--- ORGANIC MATTER EXPORT RATIOS -------------------'
       print*,'N/P organic matter Redfield ratio                   : ',par_bio_red_POP_PON
       print*,'C/P organic matter Redfield ratio                   : ',par_bio_red_POP_POC
       print*,'O2/P organic matter pseudo-Redfield ratio           : ',par_bio_red_POP_PO2
       print*,'ALK/N alkalinty correction factor                   : ',par_bio_red_PON_ALK
       print*,'Production fraction of dissolved organic matter     : ',par_bio_red_DOMfrac
       print*,'Production fraction of R-dissolved organic matter   : ',par_bio_red_RDOMfrac
       print*,'P:C fractionation during POM->DOM production        : ',par_bio_red_rP_POM_DOM
       print*,'N:C fractionation during POM->DOM production        : ',par_bio_red_rN_POM_DOM
       print*,'P:C fractionation during POM->RDOM production       : ',par_bio_red_rP_POM_RDOM
       print*,'N:C fractionation during POM->RDOM production       : ',par_bio_red_rN_POM_RDOM
       ! --- INORGANIC MATTER EXPORT RATIOS -------------------------------------------------------------------------------------- !
       print*,'--- INORGANIC MATTER EXPORT RATIOS -----------------'
       print*,'CaCO3:POC rain ratio option ID string               : ',opt_bio_CaCO3toPOCrainratio
       print*,'Base CaCO3:POC export ratio                         : ',par_bio_red_POC_CaCO3
       print*,'Exponent for modifier of CaCO3:POC export ratio     : ',par_bio_red_POC_CaCO3_pP
       print*,'Ohmega half-sat constant [Gehlen et al., 2007]      : ',par_bio_red_POC_CaCO3_Kmax
       print*,'Heinze [2004] CO2aq reference conc (umol kg-1)      : ',par_bio_red_POC_CaCO3_CO2aqREF
       print*,'Barker et al. [2003] CO3 reference conc (umol kg-1) : ',par_bio_red_POC_CaCO3_CO3REF
       print*,'Base opal:POC export ratio                          : ',par_bio_red_POC_opal
       ! --- REMINERALIZATION ---------------------------------------------------------------------------------------------------- !
       print*,'--- REMINERALIZATION -------------------------------'
       print*,'Fraction of POM remin concverted to RDOM            : ',par_bio_remin_RDOMfrac
       print*,'DOM lifetime (yrs)                                  : ',par_bio_remin_DOMlifetime
       print*,'RDOM lifetime (yrs)                                 : ',par_bio_remin_RDOMlifetime
       print*,'RDOM degradation by (surface) photolysis only?      : ',ctrl_bio_remin_RDOM_photolysis
       print*,'Specific CH4 oxidation rate (d-1)                   : ',par_bio_remin_CH4rate
       print*,'Apply fixed-profile for POM remineralization?       : ',ctrl_bio_remin_POC_fixed
       print*,'Ballasting parameterization?                        : ',ctrl_bio_remin_POC_ballast
       print*,'Initial fractional abundance of POC component #2    : ',par_bio_remin_POC_frac2
       print*,'Remineralization length #1 for POC                  : ',par_bio_remin_POC_eL1
       print*,'Remineralization length #2 for POC                  : ',par_bio_remin_POC_eL2
       print*,'Apply fixed-profile for CaCO3 remineralization?     : ',ctrl_bio_remin_CaCO3_fixed
       print*,'Initial fractional abundance of CaCO3 component #2  : ',par_bio_remin_CaCO3_frac2
       print*,'Remineralization length #1 for CaCO3                : ',par_bio_remin_CaCO3_eL1
       print*,'Remineralization length #2 for CaCO3                : ',par_bio_remin_CaCO3_eL2
       print*,'Apply fixed-profile for opal remineralization?      : ',ctrl_bio_remin_opal_fixed
       print*,'Initial fractional abundance of opal component #2   : ',par_bio_remin_opal_frac2
       print*,'Remineralization length #1 for opal                 : ',par_bio_remin_opal_eL1
       print*,'Remineralization length #2 for opal                 : ',par_bio_remin_opal_eL2
       print*,'Prescribed particle sinking rate (m d-1)            : ',par_bio_remin_sinkingrate
       print*,'Prescribed scavenging sinking rate (m d-1)          : ',par_bio_remin_sinkingrate_scav
       print*,'Organic matter carrying capacity of CaCO3           : ',par_bio_remin_ballast_kc
       print*,'Organic matter carrying capacity of opal            : ',par_bio_remin_ballast_ko
       print*,'Organic matter carrying capacity of lithogenics     : ',par_bio_remin_ballast_kl
       print*,'Aerobic remineralization of OM -> NH4 (not NO3)?    : ',ctrl_bio_remin_ONtoNH4
       print*,'Denitrification [O2] threshold (mol kg-1)           : ',par_bio_remin_denitrO2thresh
       print*,'Catch rapidly-oxidizing species going < 0.0?        : ',ctrl_bio_remin_reminfix
       ! ------------------- ISOTOPIC FRACTIONATION ------------------------------------------------------------------------------ !
       print*,'fractionation for intercellular C fixation          : ',par_d13C_DIC_Corg_ef
       print*,'fract. for intercell. C fixation of si. phytop.     : ',par_d13C_DIC_Corg_ef_sp
       print*,'fract. for intercell. C fixation of non-si. phytop. : ',par_d13C_DIC_Corg_ef_nsp
       print*,'30/28Si fractionation between H4SiO4 and opal       : ',par_d30Si_opal_epsilon
       print*,'*** d114Cd = 1.0006 ***                             : ',par_d114Cd_POCd_epsilon
       print*,'114/???Cd fractionation between Cd and CdCO3        : ',par_d114Cd_CdCO3_epsilon
       print*,'7/6Li fractionation between Li and LiCO3            : ',par_d7Li_LiCO3_epsilon
       print*,'Planktic foram 13C fractionation scheme ID string   : ',opt_bio_foram_p_13C_delta
       ! --- IRON CYCLING -------------------------------------------------------------------------------------------------------- !
       print*,'--- IRON CYCLING -----------------------------------'
       print*,'Aeolian Fe solubility                               : ',par_det_Fe_sol
       print*,'Exponent for aeolian Fe solubility                  : ',par_det_Fe_sol_exp
       print*,'Fixed cellular Fe:C ratio?                          : ',ctrl_bio_red_fixedFetoC
       print*,'C/Fe organic matter ratio                           : ',par_bio_red_POFe_POC
       print*,'Fixed scavening rate (if not: Parekh scheme)?       : ',ctrl_bio_Fe_fixedKscav
       print*,'Fixed Fe scavenging rate (d-1)                      : ',par_scav_Fe_Ks
       print*,'Parekh Fe scavenging rate scale factor: POC         : ',par_scav_Fe_sf_POC
       print*,'Parekh Fe scavenging rate scale factor: CaCO3       : ',par_scav_Fe_sf_CaCO3
       print*,'Parekh Fe scavenging rate scale factor: opal        : ',par_scav_Fe_sf_opal
       print*,'Parekh Fe scavenging rate scale factor: det         : ',par_scav_Fe_sf_det
       print*,'Fraction of scavenged Fe that can be remineralized  : ',par_scav_fremin
       print*,'Prevent return of Fe from the sediments?            : ',ctrl_bio_NO_fsedFe
       print*,'log10 of Fe ligand stability constant (K`(FeL))     : ',par_K_FeL_pP
       print*,'[FeT] dependent Fe:C ratio -- power                 : ',par_bio_FetoC_pP
       print*,'[FeT] dependent Fe:C ratio -- scaling               : ',par_bio_FetoC_K
       print*,'[FeT] dependent Fe:C ratio -- constant              : ',par_bio_FetoC_C
       ! --- SILICA CYCLING ------------------------------------------------------------------------------------------------------ !
       print*,'--- SILICA CYCLING ---------------------------------'
       ! --- NITROGEN CYCLING ---------------------------------------------------------------------------------------------------- !
       print*,'--- NITROGEN CYCLING -------------------------------'
       print*,'mu-1 max rate of export production (yr-1)           : ',par_bio_mu1
       print*,'mu-2 max rate of export from N2-fixation (yr-1)     : ',par_bio_mu2
       print*,'threshold NO3+NH4 for N2 fixation (mol kg-1)        : ',par_bio_N2fixthresh
       print*,'N:P ratio of diazotrophs                            : ',par_bio_NPdiaz
       print*,'constant for dynamical threshold for N2 fixation    : ',par_bio_N2fixdyn
       print*,'N* offset (mol kg-1)                                : ',par_bio_Nstar_offset
       ! --- TRACE METAL CYCLING ------------------------------------------------------------------------------------------------- !
       print*,'--- TRACE METAL CYCLING ----------------------------'
       print*,'Default cellular C:Cd (Cd/C) ratio                  : ',par_bio_red_POC_POCd
       print*,'[Cd/P]POM/[Cd/P]SW partition coefficient (alpha)    : ',par_bio_red_POC_POCd_alpha
       print*,'Fe-limitation dependent Cd:C uptake ratio?          : ',ctrl_bio_red_CdtoC_Felim
       print*,'Minimum (Fe replete) Cd:C uptake ratio              : ',par_bio_red_CdtoC_Felim_min
       print*,'Maximum (Fe limited) Cd:C uptake ratio              : ',par_bio_red_CdtoC_Felim_max
       print*,'Default CaCO3 Ca:Li ratio                           : ',par_bio_red_CaCO3_LiCO3
       print*,'partition coefficient (alpha)                       : ',par_bio_red_CaCO3_LiCO3_alpha
       print*,'Default CaCO3 Ca:Cd ratio                           : ',par_bio_red_CaCO3_CdCO3
       print*,'partition coefficient (alpha)                       : ',par_bio_red_CaCO3_CdCO3_alpha
       ! --- ABIOTIC PRECIPITATION ----------------------------------------------------------------------------------------------- !
       print*,'--- ABIOTIC PRECIPITATION --------------------------'
       print*,'Allow abiotic CaCO3 precipitation?                  : ',ctrl_bio_CaCO3precip
       print*,'Restrict precipitation to surface layer?            : ',ctrl_bio_CaCO3precip_sur
       print*,'Precipitate as calcite (otherwise aragonite)        : ',par_bio_CaCO3precip_calcite
       print*,'Minimum ohmega threshold for precip                 : ',par_bio_CaCO3precip_abioticohm_min
       print*,'Scale factor for CaCO3 precipitation                : ',par_bio_CaCO3precip_sf
       print*,'Rate law power for CaCO3 precipitation              : ',par_bio_CaCO3precip_exp
       ! --- I/O DIRECTORY DEFINITIONS ------------------------------------------------------------------------------------------- !
       print*,'--- I/O DIRECTORY DEFINITIONS ----------------------'
       print*,'Input dir. name                                     : ',trim(par_indir_name)
       print*,'Output dir. name                                    : ',trim(par_outdir_name)
       print*,'Restart (input) dir. name                           : ',trim(par_rstdir_name)
       print*,'Forcings (input) dir. name                          : ',trim(par_fordir_name)
       print*,'Filename for restart input                          : ',trim(par_infile_name)
       print*,'Filename for restart output                         : ',trim(par_outfile_name)
       ! --- DATA SAVING: TIME-SLICES -------------------------------------------------------------------------------------------- !
       print*,'--- BIOGEM DATA SAVING: TIME-SLICES ----------------'
       print*,'Atmospheric (interface) composition (2D)?           : ',ctrl_data_save_slice_ocnatm
       print*,'Ocean composition (3D)?                             : ',ctrl_data_save_slice_ocn
       print*,'Sediment (interface) composition (2D)?              : ',ctrl_data_save_slice_ocnsed
       print*,'Export flux?                                        : ',ctrl_data_save_sig_fexport
       print*,'Air-sea gas exchange flux (2D)?                     : ',ctrl_data_save_slice_fairsea
       print*,'Ocean-sediment flux (2D)?                           : ',ctrl_data_save_slice_focnsed
       print*,'Sediment-ocean flux (2D)?                           : ',ctrl_data_save_slice_fsedocn
       print*,'Biological fluxes (3D)?                             : ',ctrl_data_save_slice_bio
       print*,'Aqueous carbonate system properties (3D)?           : ',ctrl_data_save_slice_carb
       print*,'Aqueous carbonate system constants (3D)?            : ',ctrl_data_save_slice_carbconst
       print*,'Atmospheric physical properties (2D)?               : ',ctrl_data_save_slice_phys_atm
       print*,'Ocean physical properties (3D)?                     : ',ctrl_data_save_slice_phys_ocn
       print*,'Miscellaneous properties (-)?                       : ',ctrl_data_save_slice_misc
       print*,'Biogeochemical diagnostics (3D)?                    : ',ctrl_data_save_slice_diag
       print*,'Integration interval (yr)                           : ',par_data_save_slice_dt
       print*,'Filename for time-slice definition input            : ',trim(par_infile_slice_name)
       print*,'Number of timesteps in sub-inteval saving           : ',par_data_save_slice_n
       ! --- DATA SAVING: TIME-SERIES -------------------------------------------------------------------------------------------- !
       print*,'--- BIOGEM DATA SAVING: TIME-SERIES ----------------'
       print*,'Atmospheric (interface) composition?                : ',ctrl_data_save_sig_ocnatm
       print*,'Oceanic composition?                                : ',ctrl_data_save_sig_ocn
       print*,'Export flux?                                        : ',ctrl_data_save_sig_fexport
       print*,'Air-sea gas exchange flux ?                         : ',ctrl_data_save_sig_fairsea
       print*,'Sediment (interface) composition?                   : ',ctrl_data_save_sig_ocnsed
       print*,'Ocean->atmosphere flux?                             : ',ctrl_data_save_sig_focnatm
       print*,'Ocean->sediment flux?                               : ',ctrl_data_save_sig_focnsed
       print*,'Sediment->ocean flux/                               : ',ctrl_data_save_sig_fsedocn
       print*,'Ocean surface tracers?                              : ',ctrl_data_save_sig_ocn_sur
       print*,'Ocean surface carbonate chemistry?                  : ',ctrl_data_save_sig_carb_sur
       print*,'Miscellaneous properties?                           : ',ctrl_data_save_sig_misc
       print*,'Biogeochemical diagnostics?                         : ',ctrl_data_save_sig_diag
       print*,'Integration interval (yr)                           : ',par_data_save_sig_dt
       print*,'Filename for time-series definition input           : ',trim(par_infile_sig_name)
       ! --- DATA SAVING: MISC --------------------------------------------------------------------------------------------------- !
       print*,'--- BIOGEM DATA SAVING: MISC -----------------------'
       print*,'Save derived data (e.g., S-normalized tracers)?     : ',ctrl_data_save_derived
       print*,'Save global diagnostics (at time-slice intervals)?  : ',ctrl_data_save_GLOBAL
       print*,'Save time-slice data in ASCII format?               : ',ctrl_data_save_slice_ascii
       print*,'Save time-series data in ASCII format?              : ',ctrl_data_save_sig_ascii
       print*,'append data to output files on restart              : ',opt_append_data
       print*,'Minimum depth for benthic average (m)               : ',par_data_save_ben_Dmin
       print*,'Integration count for snap-shot saving              : ',par_t_sig_count
       print*,'Generic N j value (for time-series data saving)     : ',par_sig_j_N
       print*,'Generic S j value (for time-series data saving)     : ',par_sig_j_S
       ! --- TRACER AUDITING AND DEBUGGING OPTIONS ------------------------------------------------------------------------------- !
       print*,'--- TRACER AUDITING AND DEBUGGING OPTIONS ----------'
       print*,'Audit tracer inventory?                             : ',ctrl_audit
       print*,'Halt on audit fail?                                 : ',ctrl_audit_fatal
       print*,'Max allowed relative tracer inventory change        : ',par_misc_audit_relerr
       print*,'Report all run-time warnings?                       : ',ctrl_debug_reportwarnings
       print*,'Report level #0 debug?                              : ',ctrl_debug_lvl0
       print*,'Report level #1 debug?                              : ',ctrl_debug_lvl1
       print*,'Report level #2 debug?                              : ',ctrl_debug_lvl2
       ! --- TRACER FORCING ------------------------------------------------------------------------------------------------------ !
       print*,'--- TRACER FORCING ---------------------------------'
       DO l=1,n_l_atm
          ia = conv_iselected_ia(l)
          print*,'atm tracer forcing time scale factor  : ',trim(string_atm(ia)),' = ',par_atm_force_scale_time(ia)
          print*,'atm tracer forcing value scale factor : ',trim(string_atm(ia)),' = ',par_atm_force_scale_val(ia)
       end do
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          print*,'ocn tracer forcing time scale factor  : ',trim(string_ocn(io)),' = ',par_ocn_force_scale_time(io)
          print*,'ocn tracer forcing value scale factor : ',trim(string_ocn(io)),' = ',par_ocn_force_scale_val(io)
       end do
       print*,'i coordinate of point forcing (0 = DISABLED)        : ',par_force_point_i
       print*,'j coordinate of point forcing (0 = DISABLED)        : ',par_force_point_j
       print*,'k coordinate of point forcing (0 = DISABLED)        : ',par_force_point_k
       ! #### INSERT CODE TO LOAD ADDITIONAL PARAMETERS ########################################################################## !
       !
       ! ######################################################################################################################### !
    end if
    ! filter CaCO3:POC rain ratio options for backwards compatability
    if (ctrl_force_CaCO3toPOCrainratio) opt_bio_CaCO3toPOCrainratio = 'prescribed'
    ! ### TO BE CONVERTED TO NAMELIST ITEMS ###################################################################################### !
    par_misc_t_err = 3600.0*1.0/conv_yr_s ! time-stepping error == 1hr
    opt_data(iopt_data_save_timeslice_fnint) = .FALSE.
    opt_data(iopt_data_save_config) = .FALSE.
    opt_misc(iopt_misc_debugij) = .FALSE.
    par_misc_debug_i = 1
    par_misc_debug_j = 3
    opt_force(iopt_force_freshwater) = .FALSE.
    par_bio_c0_I = 20.0 ! half saturatin value for light (W m-2) [Doney et al., 2006] (30.0 in Parekth et al. [2005])
    par_det_Fe_frac = 0.035 ! mass fraction of Fe in dust
    par_K_FeL = 10**par_K_FeL_pP ! conditional stability constant of ligand-bound Fe [Parekth et al., 2005]
    par_scav_Fe_exp = 0.58 ! (see: Parekth et al. [2005])
    par_scav_Fe_k0  = 0.079 ! (see: Parekth et al. [2005])
    par_scav_Fe_k0 = par_scav_Fe_k0/conv_d_yr ! adjust units of scavening rate constant (d-1 -> yr-1)
    par_part_red_FeTmin = 0.125E-9 ! (see: Ridgwell [2001])
    par_part_red_FetoCmax = 250000.0 !
    ! ############################################################################################################################ !

    ! *** COMPLETE PATHS ***
    par_fordir_name = trim(par_fordir_name)//trim(par_forcing_name)//'/'

    ! *** adjust units ***
    ! adjust units of nutrient update time-scale from days to years
    par_bio_tau     = conv_d_yr*par_bio_tau
    par_bio_tau_sp  = conv_d_yr*par_bio_tau_sp
    par_bio_tau_nsp = conv_d_yr*par_bio_tau_nsp
    ! adjust units of scavening rate constant (d-1 -> yr-1)
    par_scav_Fe_ks = par_scav_Fe_ks/conv_d_yr
    ! adjust units of prescribed particulates sinking rate (m d-1 -> m yr-1)
    par_bio_remin_sinkingrate = par_bio_remin_sinkingrate/conv_d_yr
    par_bio_remin_sinkingrate_scav = par_bio_remin_sinkingrate_scav/conv_d_yr
    ! adjust units of CH4 oxidation (d-1 -> yr-1)
    par_bio_remin_CH4rate = par_bio_remin_CH4rate/conv_d_yr
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! ballast coefficients (g POC m-2 yr-1 (g ballast m-2 yr-1)-1 -> mol POC m-2 yr-1 (mol ballast m-2 yr-1)-1)
    par_bio_remin_ballast_kc = (conv_POC_cm3_mol*conv_POC_g_cm3/(conv_cal_cm3_mol*conv_cal_g_cm3))*par_bio_remin_ballast_kc
    par_bio_remin_ballast_ko = (conv_POC_cm3_mol*conv_POC_g_cm3/(conv_opal_cm3_mol*conv_opal_g_cm3))*par_bio_remin_ballast_ko
    par_bio_remin_ballast_kl = (conv_POC_cm3_mol*conv_POC_g_cm3/(conv_det_cm3_mol*conv_det_g_cm3))*par_bio_remin_ballast_kl
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END SUBROUTINE sub_load_goin_biogem
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! LOAD BioGeM RESTART DATA
  SUBROUTINE sub_load_biogem_restart()
    USE biogem_lib
    USE genie_util, ONLY:check_unit,check_iostat
    ! local variables
    integer::l,io,ios
    CHARACTER(len=255)::loc_filename                       ! 
    integer::loc_n_l_ocn,loc_n_l_sed                       ! number of selected tracers in the re-start file
    integer,DIMENSION(n_ocn)::loc_conv_iselected_io        ! 
    integer,DIMENSION(n_sed)::loc_conv_iselected_is        ! 
    ! initialize local variables
    loc_filename = TRIM(par_rstdir_name)//trim(par_infile_name)
    ! retrieve restart data
    call check_unit(in,__LINE__,__FILE__)
    OPEN(unit=in,status='old',file=loc_filename,form='unformatted',action='read',IOSTAT=ios)
    If (ios /= 0) then
       CALL sub_report_error( &
            & 'biogem_data','sub_load_biogem_restart', &
            & 'You have requested a CONTINUING run, but restart file <'//trim(loc_filename)//'> does not exist', &
            & 'SKIPPING - using default initial values (FILE: gem_config_ocn.par)', &
            & (/const_real_null/),.false. &
            & )
    else
       read(unit=in,iostat=ios)                                          &
            & loc_n_l_ocn,                                               &
            & (loc_conv_iselected_io(l),l=1,loc_n_l_ocn),                &
            & (ocn(loc_conv_iselected_io(l),:,:,:),l=1,loc_n_l_ocn),     &
            & loc_n_l_sed,                                               &
            & (loc_conv_iselected_is(l),l=1,loc_n_l_sed),                &
            & (bio_part(loc_conv_iselected_is(l),:,:,:),l=1,loc_n_l_sed)
       call check_iostat(ios,__LINE__,__FILE__)
    end if
    close(unit=in,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! adjust restart data
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       ocn(io,:,:,:) = ocn(io,:,:,:) + ocn_dinit(io)*phys_ocn(ipo_mask_ocn,:,:,:)
    end do
  end SUBROUTINE sub_load_biogem_restart
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! DATA INITIALIZATION ROUTINES
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE 'BIOLOGICAL' PARAMETERS AND VARIABLES
  SUBROUTINE sub_init_bio()
    ! local variables
    CHARACTER(len=255)::loc_filename

    ! *** initialize global arrays ***
    bio_part(:,:,:,:)     = 0.0
    bio_remin(:,:,:,:)    = 0.0
    bio_settle(:,:,:,:)   = 0.0
    bio_part_red(:,:,:,:) = 0.0

    ! *** set default 'Redfield' ratios ***
    ! trivial self-relationships(!)
    bio_part_red(is_POP,is_POP,:,:)     = 1.0
    bio_part_red(is_POC,is_POC,:,:)     = 1.0
    bio_part_red(is_PON,is_PON,:,:)     = 1.0
    bio_part_red(is_CaCO3,is_CaCO3,:,:) = 1.0
    bio_part_red(is_opal,is_opal,:,:)   = 1.0
    bio_part_red(is_POCd,is_POCd,:,:)   = 1.0
    bio_part_red(is_POFe,is_POFe,:,:)   = 1.0
    ! set values and derived values
    ! NOTE: relate everything to carbon units where it is not already
    IF (abs(par_bio_red_POP_POC) > const_real_nullsmall) then
       bio_part_red(is_POP,is_POC,:,:) = par_bio_red_POP_POC
       bio_part_red(is_POC,is_POP,:,:) = 1.0/bio_part_red(is_POP,is_POC,:,:)
    end if
    IF (abs(par_bio_red_POP_PON) > const_real_nullsmall) then
       bio_part_red(is_POP,is_PON,:,:) = par_bio_red_POP_PON
       bio_part_red(is_POC,is_PON,:,:) = bio_part_red(is_POC,is_POP,:,:)*bio_part_red(is_POP,is_PON,:,:)
       bio_part_red(is_PON,is_POC,:,:) = 1.0/bio_part_red(is_POC,is_PON,:,:)
    end if
    if (abs(par_bio_red_POC_CaCO3) > const_real_nullsmall) then
       bio_part_red(is_POC,is_CaCO3,:,:) = par_bio_red_POC_CaCO3
       bio_part_red(is_CaCO3,is_POC,:,:) = 1.0/bio_part_red(is_POC,is_CaCO3,:,:)
    end if
    if (abs(par_bio_red_POC_opal) > const_real_nullsmall) then
       bio_part_red(is_POC,is_opal,:,:) = par_bio_red_POC_opal
       bio_part_red(is_opal,is_POC,:,:) = 1.0/bio_part_red(is_POC,is_opal,:,:)
    end if
    IF (abs(par_bio_red_POC_POCd) > const_real_nullsmall) then
       bio_part_red(is_POC,is_POCd,:,:) = par_bio_red_POC_POCd
       bio_part_red(is_POCd,is_POC,:,:) = 1.0/bio_part_red(is_POC,is_POCd,:,:)
    end if
    IF (abs(par_bio_red_POFe_POC) > const_real_nullsmall) then
       bio_part_red(is_POFe,is_POC,:,:) = par_bio_red_POFe_POC
       bio_part_red(is_POC,is_POFe,:,:) = 1.0/bio_part_red(is_POFe,is_POC,:,:)
    end if
    ! denifrification and sulphate reduction
    if (par_bio_red_POP_PO2 == -138.0 ) then
       par_bio_red_O2_H2SO4 = 53.0/(-par_bio_red_POP_PO2)
       par_bio_red_O2_NO3 = 84.8/(-par_bio_red_POP_PO2)
    elseif (par_bio_red_POP_PO2 == -150.0 ) then
       par_bio_red_O2_H2SO4 = 59.0/(-par_bio_red_POP_PO2)
       par_bio_red_O2_NO3 = 104.0/(-par_bio_red_POP_PO2)
    else
       par_bio_red_O2_H2SO4 = 0.0
       par_bio_red_O2_NO3 = 0.0
    end if

    ! *** load prescribed CaCO3:POC field (if requested) ***
    if (ctrl_force_CaCO3toPOCrainratio) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_CaCO3toPOCrainratio_file)
       CALL sub_load_data_ij(loc_filename,n_i,n_j,par_bio_CaCO3toPOCrainratio(:,:))
    end if

    ! *** load prescribed POCd:POC field (if requested) ***
    if (ctrl_force_POCdtoPOCrainratio) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_POCdtoPOCrainratio_file)
       CALL sub_load_data_ij(loc_filename,n_i,n_j,par_bio_POCdtoPOCrainratio(:,:))
    end if

    ! *** load prescribed [Cd/P]POM/[Cd/P]SW partition coefficient field (if requested) ***
    if (ctrl_force_Cd_alpha) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_Cd_alpha_file)
       CALL sub_load_data_ij(loc_filename,n_i,n_j,par_bio_Cd_alpha(:,:))
    end if

    ! *** load prescribed POC scavenging coefficient field (if requested) ***
    if (ctrl_force_scav_fpart_POC) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_scav_fpart_POC_file)
       CALL sub_load_data_ijk(loc_filename,n_i,n_j,n_k,par_scav_fpart_POC(:,:,:))
    end if

    ! *** load prescribed CaCO3 scavenging coefficient field (if requested) ***
    if (ctrl_force_scav_fpart_CaCO3) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_scav_fpart_CaCO3_file)
       CALL sub_load_data_ijk(loc_filename,n_i,n_j,n_k,par_scav_fpart_CaCO3(:,:,:))
    end if

    ! *** load prescribed opal scavenging coefficient field (if requested) ***
    if (ctrl_force_scav_fpart_opal) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_scav_fpart_opal_file)
       CALL sub_load_data_ijk(loc_filename,n_i,n_j,n_k,par_scav_fpart_opal(:,:,:))
    end if

    ! *** load prescribed det scavenging coefficient field (if requested) ***
    if (ctrl_force_scav_fpart_det) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_scav_fpart_det_file)
       CALL sub_load_data_ijk(loc_filename,n_i,n_j,n_k,par_scav_fpart_det(:,:,:))
    end if

  END SUBROUTINE sub_init_bio
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! UPDATE RELATIONSHIPS BETWEEN TRACERS
  SUBROUTINE sub_update_tracerrelationships()
    IF (par_bio_prodopt /= 'NONE') then
       ! if NO3 is employed;
       ! calculate alkalnity corrections associated with the formation and destruction of organic matter from NO3
       ! otherwise, convert PO4 units to NO3 via the P:N Redfield ratio and then calculate the ALK correction from NO3
       ! NOTE: ensure that both corrections are mutually exclusive (i.e., make sure that there can be no double ALK correction)
       ! NOTE: catch incidence of par_bio_red_PON_ALK set to 0.0
       if (abs(par_bio_red_PON_ALK) > const_real_nullsmall) then
          if (ocn_select(io_NO3)) then
             conv_sed_ocn(io_ALK,is_PON) = par_bio_red_PON_ALK
             conv_ocn_sed(is_PON,io_ALK) = 1.0/conv_sed_ocn(io_ALK,is_PON)
             conv_sed_ocn(io_ALK,is_POP) = 0.0
             conv_ocn_sed(is_POP,io_ALK) = 0.0
          else
             conv_sed_ocn(io_ALK,is_PON) = 0.0
             conv_ocn_sed(is_PON,io_ALK) = 0.0
             conv_sed_ocn(io_ALK,is_POP) = par_bio_red_PON_ALK*par_bio_red_POP_PON
             conv_ocn_sed(is_POP,io_ALK) = 1.0/conv_sed_ocn(io_ALK,is_POP)
          end if
       else
          conv_sed_ocn(io_ALK,is_PON) = 0.0
          conv_ocn_sed(is_PON,io_ALK) = 0.0
          conv_sed_ocn(io_ALK,is_POP) = 0.0
          conv_ocn_sed(is_POP,io_ALK) = 0.0
       end if
       ! update O2 demand assicated with organic matter (taken as the carbon component)
       if (abs(par_bio_red_POP_POC*par_bio_red_POP_PO2) > const_real_nullsmall) then
          conv_sed_ocn(io_O2,is_POC) = par_bio_red_POP_PO2/par_bio_red_POP_POC
          conv_ocn_sed(is_POC,io_O2) = 1.0/conv_sed_ocn(io_O2,is_POC)
       else
          conv_sed_ocn(io_O2,is_POC) = 0.0
          conv_ocn_sed(is_POC,io_O2) = 0.0
       end if
    end if
  END SUBROUTINE sub_update_tracerrelationships
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE INTEGRATED TIME-SLICE VALUE ARRAYS
  SUBROUTINE sub_init_int_timeslice()
    ! initialize integrated time
    int_t_timeslice = 0.0
    int_t_timeslice_count = 0
    ! initialize time-slice data - ocn
    int_ocn_timeslice(:,:,:,:)        = 0.0
    int_bio_part_timeslice(:,:,:,:)   = 0.0
    int_bio_settle_timeslice(:,:,:,:) = 0.0
    int_bio_remin_timeslice(:,:,:,:)  = 0.0
    int_phys_ocn_timeslice(:,:,:,:)   = 0.0
    int_carb_timeslice(:,:,:,:)       = 0.0
    int_carbconst_timeslice(:,:,:,:)  = 0.0
    int_carbisor_timeslice(:,:,:,:)   = 0.0
    ! initialize time-slice data - ocn-atm
    int_sfcatm1_timeslice(:,:,:)     = 0.0
    int_focnatm_timeslice(:,:,:)     = 0.0
    int_phys_ocnatm_timeslice(:,:,:) = 0.0
    ! initialize time-slice data - ocn-sed
    int_sfcsed1_timeslice(:,:,:) = 0.0
    int_focnsed_timeslice(:,:,:) = 0.0
    int_fsedocn_timeslice(:,:,:) = 0.0
    ! initialize time-slice data - GOLDSTEIn
    int_opsi_timeslice(:,:)  = 0.0
    int_opsia_timeslice(:,:) = 0.0
    int_opsip_timeslice(:,:) = 0.0
    int_zpsi_timeslice(:,:)  = 0.0
    int_u_timeslice(:,:,:,:) = 0.0
    ! integrated time slice storage arrays - diagnostics
    int_diag_bio_timeslice(:,:,:)       = 0.0
    int_diag_geochem_timeslice(:,:,:,:) = 0.0
    int_diag_weather_timeslice(:,:,:)   = 0.0
    int_diag_airsea_timeslice(:,:,:)    = 0.0
    ! ### ADD ADDITIONAL TIME-SLICE ARRAY INITIALIZATIONS HERE ################################################################### !
    !
    ! ############################################################################################################################ !
  END SUBROUTINE sub_init_int_timeslice
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE INTEGRATED TIME-SERIES VALUE ARRAYS
  SUBROUTINE sub_init_int_timeseries()
    ! initialize integrated time
    int_t_sig       = 0.0
    int_t_sig_count = 0
    ! initialize time-series data
    int_misc_gemlite_sig    = 0.0
    int_ocn_tot_M_sig       = 0.0
    int_ocn_tot_M_sur_sig   = 0.0
    int_ocn_tot_V_sig       = 0.0
    int_ocn_sig(:)          = 0.0
    int_fexport_sig(:)      = 0.0
    int_ocnatm_sig(:)       = 0.0
    int_focnatm_sig(:)      = 0.0
    int_focnsed_sig(:)      = 0.0
    int_fsedocn_sig(:)      = 0.0
    int_ocn_sur_sig(:)      = 0.0
    int_ocn_ben_sig(:)      = 0.0
    int_carb_sur_sig(:)     = 0.0
    int_carb_ben_sig(:)     = 0.0
    int_misc_seaice_sig     = 0.0
    int_misc_seaice_sig_th  = 0.0
    int_misc_seaice_sig_vol = 0.0
    int_misc_opsi_min_sig   = 0.0
    int_misc_opsi_max_sig   = 0.0
    int_misc_opsia_min_sig  = 0.0
    int_misc_opsia_max_sig  = 0.0
    int_misc_SLT_sig        = 0.0
    int_misc_det_Fe_tot_sig = 0.0
    int_misc_det_Fe_dis_sig = 0.0
    int_misc_ocn_solfor_sig = 0.0
    int_misc_ocn_fxsw_sig   = 0.0
    int_ocnsed_sig(:)       = 0.0
    int_diag_bio_sig(:)     = 0.0
    int_diag_geochem_sig(:) = 0.0
    int_diag_weather_sig(:) = 0.0
    int_diag_airsea_sig(:)  = 0.0
    int_diag_misc_2D_sig(:) = 0.0
    ! ### ADD ADDITIONAL TIME-SERIES ARRAY INITIALIZATIONS HERE ################################################################## !
    !
    ! ############################################################################################################################ !
  END SUBROUTINE sub_init_int_timeseries
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE FORCING ARRAYS
  SUBROUTINE sub_init_force()
    force_restore_ocn(:,:,:,:)    = 0.0
    force_restore_ocn_I(:,:,:,:)  = 0.0
    force_restore_ocn_II(:,:,:,:) = 0.0
    force_restore_ocn_sig(:,:,:)  = 0.0
    force_restore_ocn_sig_x(:)    = 0.0
    force_restore_ocn_sig_i(:,:)  = 0
    force_restore_ocn_tconst      = 0.0
    force_restore_ocn_select(:)   = .FALSE.
    force_restore_ocn_sur(:)      = .FALSE.
    force_restore_atm(:,:,:)      = 0.0
    force_restore_atm_I(:,:,:)    = 0.0
    force_restore_atm_II(:,:,:)   = 0.0
    force_restore_atm_sig(:,:,:)  = 0.0
    force_restore_atm_sig_x(:)    = 0.0
    force_restore_atm_sig_i(:,:)  = 0
    force_restore_atm_tconst      = 0.0
    force_restore_atm_select(:)   = .FALSE.
    !force_restore_sed(:,:,:)      = 0.0
    !force_restore_sed_I(:,:,:)    = 0.0
    !force_restore_sed_II(:,:,:)   = 0.0
    !force_restore_sed_sig(:,:,:)  = 0.0
    !force_restore_sed_sig_x(:)    = 0.0
    !force_restore_sed_sig_i(:,:)  = 0
    !force_restore_sed_tconst      = 0.0
    !force_restore_sed_select(:)   = .FALSE.
    force_flux_ocn(:,:,:,:)       = 0.0
    force_flux_ocn_I(:,:,:,:)     = 0.0
    force_flux_ocn_II(:,:,:,:)    = 0.0
    force_flux_ocn_sig(:,:,:)     = 0.0
    force_flux_ocn_sig_x(:)       = 0.0
    force_flux_ocn_sig_i(:,:)     = 0
    force_flux_ocn_select(:)      = .FALSE.
    force_flux_ocn_scale(:)       = .FALSE.
    force_flux_atm(:,:,:)         = 0.0
    force_flux_atm_I(:,:,:)       = 0.0
    force_flux_atm_II(:,:,:)      = 0.0
    force_flux_atm_sig(:,:,:)     = 0.0
    force_flux_atm_sig_x(:)       = 0.0
    force_flux_atm_sig_i(:,:)     = 0
    force_flux_atm_select(:)      = .FALSE.
    force_flux_sed_scale(:)       = .FALSE.
    force_flux_sed(:,:,:)         = 0.0
    force_flux_sed_I(:,:,:)       = 0.0
    force_flux_sed_II(:,:,:)      = 0.0
    force_flux_sed_sig(:,:,:)     = 0.0
    force_flux_sed_sig_x(:)       = 0.0
    force_flux_sed_sig_i(:,:)     = 0
    force_flux_sed_select(:)      = .FALSE.
    force_flux_sed_scale(:)       = .FALSE.
    ! misc 
    force_solconst_sig(:,:)       = 0.0
    force_restore_docn_nuts(:)    = 0.0
    force_atm_uniform(:)          = 2
    force_ocn_uniform(:)          = 2
    force_sed_uniform(:)          = 2
    force_atm_point_i(:)          = 01
    force_ocn_point_i(:)          = 01
    force_sed_point_i(:)          = 01
    force_atm_point_j(:)          = 01
    force_ocn_point_j(:)          = 01
    force_sed_point_j(:)          = 01
    force_ocn_point_k(:)          = 01
    ! ### ADD ADDITIONAL FORCINGS ARRAY INITIALIZATIONS HERE ##################################################################### !
    ! 
    ! ############################################################################################################################ !
  END SUBROUTINE sub_init_force
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE AUDIT INVENTORY ARRAYS
  SUBROUTINE sub_init_audit()
    audit_ocn_init(:)       = 0.0
    audit_ocn_old(:)        = 0.0
    audit_ocn_new(:)        = 0.0
    audit_ocn_delta(:)      = 0.0
  END SUBROUTINE sub_init_audit
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE DIAGNOSTICS ARRAYS
  SUBROUTINE sub_init_diag()
    diag_bio(:,:,:)       = 0.0
    diag_geochem(:,:,:,:) = 0.0
 !!!   diag_weather(:,:,:)   = 0.0
    diag_airsea(:,:,:)    = 0.0
  END SUBROUTINE sub_init_diag
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE 'PHYSICS' - OCEAN
  SUBROUTINE sub_init_phys_ocn()
    ! local variables
    INTEGER::i,j,k
    REAL,DIMENSION(0:n_k+1)::loc_grid_dz,loc_grid_dza
    ! initialize local variables
    loc_grid_dz(0:n_k+1)  = 0.0
    loc_grid_dz(1:n_k)    = goldstein_dz(:)
    loc_grid_dza(0:n_k+1) = 0.0
    loc_grid_dza(1:n_k)   = goldstein_dza(:); loc_grid_dza(n_k) = loc_grid_dz(n_k)/2.0
    ! zero array
    phys_ocn(:,:,:,:) = 0.0
    ! initialize array values
    ! NOTE: initialize basic grid structure values for the (i,j,k) grid, not just ocean-only points
    ! NOTE: depth in in unit of m BELOW sealevel (i.e., a +ve scale)
    ! NOTE: set default rho
    DO i=1,n_i
       DO j=1,n_j
          DO k=1,n_k
             phys_ocn(ipo_lat,i,j,k)      = (180.0/const_pi)*ASIN(goldstein_s(j))
             phys_ocn(ipo_lon,i,j,k)      = (360.0/n_i)*(real(i)-0.5) + par_grid_lon_offset
             phys_ocn(ipo_dlat,i,j,k)     = (180.0/const_pi)*(ASIN(goldstein_sv(j)) - ASIN(goldstein_sv(j-1)))
             phys_ocn(ipo_dlon,i,j,k)     = (360.0/n_i)
             phys_ocn(ipo_latn,i,j,k)     = (180.0/const_pi)*ASIN(goldstein_sv(j))
             phys_ocn(ipo_lone,i,j,k)     = (360.0/n_i)*real(i) + par_grid_lon_offset
             phys_ocn(ipo_Dmid,i,j,k)     = SUM(goldstein_dsc*loc_grid_dza(k:n_k))
             phys_ocn(ipo_dD,i,j,k)       = goldstein_dsc*loc_grid_dz(k)
             phys_ocn(ipo_Dbot,i,j,k)     = SUM(goldstein_dsc*loc_grid_dz(k:n_k))
             phys_ocn(ipo_Dtop,i,j,k)     = SUM(goldstein_dsc*loc_grid_dz(k+1:n_k+1))
          end do
          DO k=goldstein_k1(i,j),n_k
             phys_ocn(ipo_A,i,j,k)        = 2.0*const_pi*(const_rEarth**2)*(1.0/n_i)*(goldstein_sv(j) - goldstein_sv(j-1))
             phys_ocn(ipo_rA,i,j,k)       = 1.0 / phys_ocn(ipo_A,i,j,k)
             phys_ocn(ipo_V,i,j,k)        = phys_ocn(ipo_dD,i,j,k)*phys_ocn(ipo_A,i,j,k)
             phys_ocn(ipo_M,i,j,k)        = conv_m3_kg*phys_ocn(ipo_V,i,j,k)
             phys_ocn(ipo_rM,i,j,k)       = 1.0 / phys_ocn(ipo_M,i,j,k)
             phys_ocn(ipo_mask_ocn,i,j,k) = 1.0
             phys_ocn(ipo_rho,i,j,k)      = conv_m3_kg
          END DO
       END DO
    END DO
  END SUBROUTINE sub_init_phys_ocn
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE 'PHYSICS' - OCEAN-ATMOSPHERE INTERFACE
  SUBROUTINE sub_init_phys_ocnatm()
    ! local variables
    INTEGER::i,j
    CHARACTER(len=255)::loc_filename
    ! zero array
    phys_ocnatm(:,:,:) = 0.0
    ! initialize array values
    DO i=1,n_i
       DO j=1,n_j
          phys_ocnatm(ipoa_lat,i,j)  = (180.0/const_pi)*ASIN(goldstein_s(j))
          phys_ocnatm(ipoa_lon,i,j)  = (360.0/n_i)*(real(i)-0.5) + par_grid_lon_offset
          phys_ocnatm(ipoa_dlat,i,j) = (180.0/const_pi)*(ASIN(goldstein_sv(j)) - ASIN(goldstein_sv(j-1)))
          phys_ocnatm(ipoa_dlon,i,j) = (360.0/n_i)
          phys_ocnatm(ipoa_A,i,j)    = 2.0*const_pi*(const_rEarth**2)*(1.0/n_i)*(goldstein_sv(j) - goldstein_sv(j-1))
          phys_ocnatm(ipoa_rA,i,j)   = 1.0/ phys_ocnatm(ipoa_A,i,j)
          IF (n_k >= goldstein_k1(i,j)) THEN
             phys_ocnatm(ipoa_seaice,i,j) = 0.0
             phys_ocnatm(ipoa_u,i,j)      = 0.0
             phys_ocnatm(ipoa_mask_ocn,i,j) = 1.0
          END IF
       END DO
    END DO
    ! load prescribed sea-ice cover (if requested)
    ! NOTE: convert from %cover to fractional cover
    if (ctrl_force_seaice) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_seaice_file)
       CALL sub_load_data_ij(loc_filename,n_i,n_j,par_phys_seaice(:,:))
       par_phys_seaice(:,:) = par_phys_seaice(:,:)/100.0
    end if
    ! load prescribed wind-speed (if requested)
    ! NOTE: (m s-1)
    if (ctrl_force_windspeed) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_windspeed_file)
       CALL sub_load_data_ij(loc_filename,n_i,n_j,par_phys_windspeed(:,:))
    end if
  END SUBROUTINE sub_init_phys_ocnatm
  ! ****************************************************************************************************************************** !


  ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
  ! INITIALIZE 'PHYSICS' - OCEAN
  SUBROUTINE sub_data_init_phys_ocn()
    ! local variables
    INTEGER::i,j,k
    integer::loc_n,loc_k1
    REAL,DIMENSION(0:n_k+1)::loc_grid_dz,loc_grid_dza
    ! initialize local variables
    loc_grid_dz(0:n_k+1)  = 0.0
    loc_grid_dz(1:n_k)    = goldstein_dz(:)
    loc_grid_dza(0:n_k+1) = 0.0
    loc_grid_dza(1:n_k)   = goldstein_dza(:); loc_grid_dza(n_k) = loc_grid_dz(n_k)/2.0
    ! initialize array values
    ! NOTE: depth in in unit of m BELOW sealevel (i.e., a +ve scale)
    ! NOTE: set default rho
    loc_n = 0
    DO i=1,n_i
       DO j=1,n_j
          loc_k1 = goldstein_k1(i,j)
          IF (n_k >= loc_k1) THEN
             loc_n = loc_n + 1
             vphys_ocn(loc_n)%i = i
             vphys_ocn(loc_n)%j = j
             vphys_ocn(loc_n)%k1 = loc_k1
             ! initialize, becasue not all 'k' depths are valid
             vphys_ocn(loc_n)%mk(:,:) = 0.0
             DO k=n_k,loc_k1,-1
                vphys_ocn(loc_n)%mk(ipo_lat,k)      = (180.0/const_pi)*ASIN(goldstein_s(j))
                vphys_ocn(loc_n)%mk(ipo_lon,k)      = (360.0/n_i)*(real(i)-0.5) + par_grid_lon_offset
                vphys_ocn(loc_n)%mk(ipo_dlat,k)     = (180.0/const_pi)*(ASIN(goldstein_sv(j)) - ASIN(goldstein_sv(j-1)))
                vphys_ocn(loc_n)%mk(ipo_dlon,k)     = (360.0/n_i)
                vphys_ocn(loc_n)%mk(ipo_latn,k)     = (180.0/const_pi)*ASIN(goldstein_sv(j))
                vphys_ocn(loc_n)%mk(ipo_lone,k)     = (360.0/n_i)*real(i) + par_grid_lon_offset
                vphys_ocn(loc_n)%mk(ipo_Dmid,k)     = SUM(goldstein_dsc*loc_grid_dza(k:n_k))
                vphys_ocn(loc_n)%mk(ipo_dD,k)       = goldstein_dsc*loc_grid_dz(k)
                vphys_ocn(loc_n)%mk(ipo_Dbot,k)     = SUM(goldstein_dsc*loc_grid_dz(k:n_k))
                vphys_ocn(loc_n)%mk(ipo_Dtop,k)     = SUM(goldstein_dsc*loc_grid_dz(k+1:n_k+1))
                vphys_ocn(loc_n)%mk(ipo_A,k)        = 2.0*const_pi*(const_rEarth**2)*(1.0/real(n_i))* &
                     & (goldstein_sv(j) - goldstein_sv(j-1))
                vphys_ocn(loc_n)%mk(ipo_rA,k)       = 1.0 / vphys_ocn(loc_n)%mk(ipo_A,k)
                vphys_ocn(loc_n)%mk(ipo_V,k)        = vphys_ocn(loc_n)%mk(ipo_dD,k)*vphys_ocn(loc_n)%mk(ipo_A,k)
                vphys_ocn(loc_n)%mk(ipo_M,k)        = conv_m3_kg*vphys_ocn(loc_n)%mk(ipo_V,k)
                vphys_ocn(loc_n)%mk(ipo_rM,k)       = 1.0 / vphys_ocn(loc_n)%mk(ipo_M,k)
                vphys_ocn(loc_n)%mk(ipo_mask_ocn,k) = 1.0
                vphys_ocn(loc_n)%mk(ipo_rho,k)      = conv_m3_kg
             end DO
          end if
       end do
    end do
  END SUBROUTINE sub_data_init_phys_ocn
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !


  ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
  ! INITIALIZE 'PHYSICS' - OCEAN-ATMOSPHERE INTERFACE
  SUBROUTINE sub_data_init_phys_ocnatm()
    ! local variables
    INTEGER::i,j
    integer::loc_n
    CHARACTER(len=255)::loc_filename
    loc_n = 0
    DO i=1,n_i
       DO j=1,n_j
          loc_n = loc_n + 1
          vphys_ocnatm(loc_n)%i = i
          vphys_ocnatm(loc_n)%j = j
          ! initialize, becasue not all 'k' depths are valid
          vphys_ocnatm(loc_n)%m(:) = 0.0
          vphys_ocnatm(loc_n)%m(ipoa_lat)  = (180.0/const_pi)*ASIN(goldstein_s(j))
          vphys_ocnatm(loc_n)%m(ipoa_lon)  = (360.0/n_i)*(real(i)-0.5) + par_grid_lon_offset
          vphys_ocnatm(loc_n)%m(ipoa_dlat) = (180.0/const_pi)*(ASIN(goldstein_sv(j)) - ASIN(goldstein_sv(j-1)))
          vphys_ocnatm(loc_n)%m(ipoa_dlon) = (360.0/n_i)
          vphys_ocnatm(loc_n)%m(ipoa_A)    = 2.0*const_pi*(const_rEarth**2)*(1.0/n_i)*(goldstein_sv(j) - goldstein_sv(j-1))
          vphys_ocnatm(loc_n)%m(ipoa_rA)   = 1.0/ vphys_ocnatm(loc_n)%m(ipoa_A)
          IF (n_k >= goldstein_k1(i,j)) THEN
             vphys_ocnatm(loc_n)%m(ipoa_seaice) = 0.0
             vphys_ocnatm(loc_n)%m(ipoa_u)      = 0.0
             vphys_ocnatm(loc_n)%m(ipoa_mask_ocn) = 1.0
          END IF
       end do
    end do
    ! load prescribed sea-ice cover (if requested)
    ! NOTE: convert from %cover to fractional cover
    if (ctrl_force_seaice) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_seaice_file)
       CALL sub_load_data_ij(loc_filename,n_i,n_j,par_phys_seaice(:,:))
       par_phys_seaice(:,:) = par_phys_seaice(:,:)/100.0
    end if
    ! load prescribed wind-speed (if requested)
    ! NOTE: (m s-1)
    if (ctrl_force_windspeed) then
       loc_filename = TRIM(par_indir_name)//TRIM(par_windspeed_file)
       CALL sub_load_data_ij(loc_filename,n_i,n_j,par_phys_windspeed(:,:))
    end if
  END SUBROUTINE sub_data_init_phys_ocnatm
  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !


  ! ****************************************************************************************************************************** !
  ! CONFIGURE AND INITIALIZE TRACER COMPOSITION - OCEAN
  SUBROUTINE sub_init_tracer_ocn_comp()
    ! local variables
    INTEGER::i,j,k,io
    real::loc_tot,loc_frac,loc_standard
    ! initialize global arrays
    ocn(:,:,:,:) = 0.0
    ! set <ocn> array
    DO i=1,n_i
       DO j=1,n_j
          DO k=goldstein_k1(i,j),n_k
             DO io=1,n_ocn
                IF (ocn_select(io)) THEN
                   SELECT CASE (ocn_type(io))
                   CASE (1)
                      ocn(io,i,j,k) = ocn_init(io)
                   case (11:20)
                      loc_tot  = ocn_init(ocn_dep(io))
                      loc_standard = const_standards(ocn_type(io))
                      loc_frac = fun_calc_isotope_fraction(ocn_init(io),loc_standard)
                      ocn(io,i,j,k) = loc_frac*loc_tot
                   END SELECT
                end IF
             end DO
          END DO
       END DO
    END DO
    ! close file pipe
    CLOSE(unit=in)
  END SUBROUTINE sub_init_tracer_ocn_comp
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! CONFIGURE AND INITIALIZE TRACER FORCING - ATMOSPHERE
  SUBROUTINE sub_init_tracer_forcing_atm()
    USE genie_util, ONLY:check_unit,check_iostat
    ! local variables
    INTEGER::n,ia,ios
    INTEGER::loc_n_elements,loc_n_start
    REAL::loc_force_restore_tconst
    LOGICAL::loc_force_restore_select,loc_force_flux_select,loc_force_flux_scale
    logical::loc_airsea_eqm
    integer::loc_ia
    integer::loc_force_uniform
    integer::loc_force_point_i,loc_force_point_j
    CHARACTER(len=255)::loc_filename
    ! initialize global variables.
    force_restore_atm_select(:)   = .FALSE.
    force_flux_atm_select(:)      = .FALSE.
    force_flux_atm_scale(:)       = .FALSE.
    ocnatm_airsea_eqm(:)          = .FALSE.
    ! check file format
    loc_filename = TRIM(par_fordir_name)//'configure_forcings_atm.dat'
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! open file pipe
    call check_unit(in,__LINE__,__FILE__)
    OPEN(unit=in,file=loc_filename,action='read',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
    END DO
    ! read in default (uniform) atmopshere tracer values
    DO n = 1,loc_n_elements
       if (ctrl_force_oldformat) then
          READ(unit=in,FMT=*,iostat=ios)   &
               & loc_ia,                   & ! COLUMN #00: TRACER NUMBER
               & loc_force_restore_select, & ! COLUMN #01: include restoring forcing of tracer?
               & loc_force_restore_tconst, & ! COLUMN #02: time constant of restoring forcing (years)
               & loc_force_flux_select,    & ! COLUMN #03: include flux forcing of tracer?
               & loc_force_flux_scale,     & ! COLUMN #04: scale flux forcing of tracer?
               & loc_airsea_eqm              ! COLUMN #05: assume ocean in equilibrium with atmosphere?
          call check_iostat(ios,__LINE__,__FILE__)
          loc_force_uniform = -99
          loc_force_point_i = 0
          loc_force_point_j = 0
       else
          READ(unit=in,FMT=*,iostat=ios)   &
               & loc_ia,                   & ! COLUMN #00: TRACER NUMBER
               & loc_force_restore_select, & ! COLUMN #01: include restoring forcing of tracer?
               & loc_force_restore_tconst, & ! COLUMN #02: time constant of restoring forcing (years)
               & loc_force_flux_select,    & ! COLUMN #03: include flux forcing of tracer?
               & loc_force_flux_scale,     & ! COLUMN #04: scale flux forcing of tracer?
               & loc_airsea_eqm,           & ! COLUMN #05: assume ocean in equilibrium with atmosphere?
               & loc_force_uniform,        & ! COLUMN #06: make forcing uniform over this dimension
               & loc_force_point_i,        & ! COLUMN #07: i grid location of point forcing
               & loc_force_point_j           ! COLUMN #08: j grid location of point forcing
          call check_iostat(ios,__LINE__,__FILE__)
       end if
       ia = loc_ia
       force_restore_atm_select(ia) = loc_force_restore_select
       force_restore_atm_tconst(ia) = loc_force_restore_tconst
       force_flux_atm_select(ia)    = loc_force_flux_select
       force_flux_atm_scale(ia)     = loc_force_flux_scale
       ocnatm_airsea_eqm(ia)        = loc_airsea_eqm
       force_atm_uniform(ia)        = loc_force_uniform
       force_atm_point_i(ia)        = loc_force_point_i
       force_atm_point_j(ia)        = loc_force_point_j
       if (force_atm_uniform(ia) > 0) force_flux_atm_scale(ia) = .true.
       if (loc_force_restore_select .AND. (loc_force_restore_tconst < const_real_nullsmall)) then
          CALL sub_report_error( &
               & 'biogem_data','sub_init_atm', &
               & 'Please do not set elected tracer restoring constants to zero (gem_config_atm.par) - '// &
               & 'it can only lead to much unpleasantness later on', &
               & 'STOPPING', &
               & (/const_real_null/),.TRUE. &
               & )
       end if
       IF (loc_force_restore_select .AND. loc_force_flux_select) then
          CALL sub_report_error( &
               & 'biogem_data','init_atm', &
               & 'You are being greedy ... and have both flux AND restoring atmospheric forcing selected'// &
               & '(gem_config_atm.par) - Is this really what you intended?', &
               & 'CONTINUING', &
               & (/const_real_null/),.false. &
               & )
       end if
    END DO
    ! close file pipe
    CLOSE(unit=in,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! blanket namelist over-ride of forcing point source
    IF ((par_force_point_i > 0) .AND. (par_force_point_j > 0)) then
       force_atm_point_i(:) = par_force_point_i
       force_atm_point_j(:) = par_force_point_j
    end IF
  END SUBROUTINE sub_init_tracer_forcing_atm
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! CONFIGURE AND INITIALIZE TRACER FORCING - OCEAN
  SUBROUTINE sub_init_tracer_forcing_ocn()
    USE genie_util, ONLY:check_unit,check_iostat
    ! local variables
    INTEGER::n,io,ios
    INTEGER::loc_n_elements,loc_n_start
    REAL::loc_force_restore_tconst
    LOGICAL::loc_force_restore_select,loc_force_restore_sur
    LOGICAL::loc_force_flux_select,loc_force_flux_scale
    integer::loc_io
    integer::loc_force_uniform
    integer::loc_force_point_i,loc_force_point_j,loc_force_point_k
    CHARACTER(len=255)::loc_filename
    ! initialize global arrays
    force_restore_ocn_select(:) = .FALSE.
    force_flux_ocn_select(:)    = .FALSE.
    force_flux_ocn_scale(:)     = .FALSE.
    ! check file format
    loc_filename = TRIM(par_fordir_name)//'configure_forcings_ocn.dat'
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! open file pipe
    call check_unit(in,__LINE__,__FILE__)
    OPEN(unit=in,file=loc_filename,action='read',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
    END DO
    ! 
    DO n = 1,loc_n_elements
       if (ctrl_force_oldformat) then
          READ(unit=in,FMT=*,iostat=ios)   &
               & loc_io,                   & ! COLUMN #00: TRACER NUMBER
               & loc_force_restore_select, & ! COLUMN #01: include restoring forcing of tracer?
               & loc_force_restore_sur,    & ! COLUMN #02: restrict restoring forcing to surface?
               & loc_force_restore_tconst, & ! COLUMN #03: time constant of restoring forcing (years)
               & loc_force_flux_select,    & ! COLUMN #04: include flux forcing of tracer? 
               & loc_force_flux_scale        ! COLUMN #05: scale flux forcing of tracer?
          call check_iostat(ios,__LINE__,__FILE__)
          loc_force_uniform = -99
          loc_force_point_i = 0
          loc_force_point_j = 0
          loc_force_point_k = 0
       else
          READ(unit=in,FMT=*,iostat=ios)   &
               & loc_io,                   & ! COLUMN #00: TRACER NUMBER
               & loc_force_restore_select, & ! COLUMN #01: include restoring forcing of tracer?
               & loc_force_restore_sur,    & ! COLUMN #02: restrict restoring forcing to surface?
               & loc_force_restore_tconst, & ! COLUMN #03: time constant of restoring forcing (years)
               & loc_force_flux_select,    & ! COLUMN #04: include flux forcing of tracer? 
               & loc_force_flux_scale,     & ! COLUMN #05: scale flux forcing of tracer?
               & loc_force_uniform,        & ! COLUMN #06: make forcing uniform over this dimension
               & loc_force_point_i,        & ! COLUMN #07: i grid location of point forcing
               & loc_force_point_j,        & ! COLUMN #08: j grid location of point forcing
               & loc_force_point_k           ! COLUMN #09: k grid location of point forcing
          call check_iostat(ios,__LINE__,__FILE__)
       endif
       io = loc_io
       force_restore_ocn_select(io) = loc_force_restore_select
       force_restore_ocn_sur(io)    = loc_force_restore_sur
       force_restore_ocn_tconst(io) = loc_force_restore_tconst
       force_flux_ocn_select(io)    = loc_force_flux_select
       force_flux_ocn_scale(io)     = loc_force_flux_scale
       force_ocn_uniform(io)        = loc_force_uniform
       force_ocn_point_i(io)        = loc_force_point_i
       force_ocn_point_j(io)        = loc_force_point_j
       force_ocn_point_k(io)        = loc_force_point_k
       if (force_ocn_uniform(io) > 0) force_flux_ocn_scale(io) = .true.
       ! set local depth limit for ocean boundary conditions (i.e., as as to achieve a surface-only forcing)
       ! NOTE: ensure that land-surface information is preserved (i.e, 'k > n_k')
       ! NOTE: there is currently no restriction of flux forcing to the ocean surface only
       IF (force_restore_ocn_sur(io)) THEN
          force_restore_ocn_k1(io,:,:) = MAX(goldstein_k1(:,:),n_k)
       ELSE
          force_restore_ocn_k1(io,:,:) = goldstein_k1(:,:)
       END IF
       force_flux_ocn_k1(io,:,:) = goldstein_k1(:,:)
       if (loc_force_restore_select .AND. (loc_force_restore_tconst < const_real_nullsmall)) then
          CALL sub_report_error( &
               & 'biogem_data','sub_init_ocn', &
               & 'Please do not set selected tracer restoring constants to zero (gem_config_ocn.par) - '// &
               & 'it can only lead to much unpleasantness later on', &
               & 'STOPPING', &
               & (/const_real_null/),.TRUE. &
               & )
       end if
       IF (loc_force_restore_select .AND. loc_force_flux_select) then
          CALL sub_report_error( &
               & 'biogem_data','init_atm', &
               & 'You are being greedy ... and have both flux AND restoring atmospheric forcing selected'// &
               & '(gem_config_atm.par) - Is this really what you intended?', &
               & 'CONTINUING', &
               & (/const_real_null/),.false. &
               & )
       end if
    end DO
    ! close file pipe
    CLOSE(unit=in,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! blanket namelist over-ride of forcing point source
    IF ((par_force_point_i > 0) .AND. (par_force_point_j > 0) .AND. (par_force_point_k > 0)) then
       force_ocn_point_i(:) = par_force_point_i
       force_ocn_point_j(:) = par_force_point_j
       force_ocn_point_k(:) = par_force_point_k
    end IF

  END SUBROUTINE sub_init_tracer_forcing_ocn
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! CONFIGURE AND INITIALIZE TRACER FORCING - SEDIMENTS
  SUBROUTINE sub_init_tracer_forcing_sed()
    USE genie_util, ONLY:check_unit,check_iostat
    ! local variables
    INTEGER::n,is,ios
    INTEGER::loc_n_elements,loc_n_start
    LOGICAL::loc_force_flux_select,loc_force_flux_scale
    integer::loc_is
    integer::loc_force_uniform
    integer::loc_force_point_i,loc_force_point_j
    CHARACTER(len=255)::loc_filename
    ! initialize global variables
    force_flux_sed_select(:) = .FALSE.
    force_flux_sed_scale(:)  = .FALSE.
    ! check file format
    loc_filename = TRIM(par_fordir_name)//'configure_forcings_sed.dat'
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! open file pipe
    call check_unit(in,__LINE__,__FILE__)
    OPEN(unit=in,file=loc_filename,action='read',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
    END DO
    ! read in default sediment tracer info
    DO n = 1,loc_n_elements
       if (ctrl_force_oldformat) then
          READ(unit=in,FMT=*,iostat=ios) &
               & loc_is,                 & ! COLUMN #00: TRACER NUMBER
               & loc_force_flux_select,  & ! COLUMN #01: include flux forcing of tracer?
               & loc_force_flux_scale      ! COLUMN #02: scale flux forcing of tracer?
          call check_iostat(ios,__LINE__,__FILE__)
          loc_force_uniform = -99
          loc_force_point_i = 0
          loc_force_point_j = 0
       else
          READ(unit=in,FMT=*,iostat=ios) &
               & loc_is,                 & ! COLUMN #00: TRACER NUMBER
               & loc_force_flux_select,  & ! COLUMN #01: include flux forcing of tracer?
               & loc_force_flux_scale,   & ! COLUMN #02: scale flux forcing of tracer?
               & loc_force_uniform,      & ! COLUMN #03: make forcing uniform over this dimension
               & loc_force_point_i,      & ! COLUMN #04: i grid location of point forcing
               & loc_force_point_j         ! COLUMN #05: j grid location of point forcing
          call check_iostat(ios,__LINE__,__FILE__)
       end if
       is = loc_is
       force_flux_sed_select(is) = loc_force_flux_select
       force_flux_sed_scale(is)  = loc_force_flux_scale
       force_sed_uniform(is)     = loc_force_uniform
       force_sed_point_i(is)     = loc_force_point_i
       force_sed_point_j(is)     = loc_force_point_j
       if (force_sed_uniform(is) > 0) force_flux_sed_scale(is) = .true.
    END DO
    ! close file pipe
    CLOSE(unit=in,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! blanket namelist over-ride of forcing point source
    IF ((par_force_point_i > 0) .AND. (par_force_point_j > 0)) then
       force_sed_point_i(:) = par_force_point_i
       force_sed_point_j(:) = par_force_point_j
    end IF
  END SUBROUTINE sub_init_tracer_forcing_sed
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! META-OPTION SETUP AND PARAMETER VALUE CONSISTENCY CHECK
  SUBROUTINE sub_check_par_biogem()
    ! local variables
    LOGICAL::loc_flag
    integer::loc_i,loc_tot_i
    CHARACTER(len=255)::loc_string
    CHARACTER(len=255)::loc_string1,loc_string2
    integer::l,io,ia,is

    ! *** set-up ***
    ! initialize variables
    loc_flag = .FALSE.
    opt_select(:) = .FALSE.
    ! set derived tracer selection options
    opt_select(iopt_select_carbchem)   = ocn_select(io_DIC) .AND. ocn_select(io_ALK)
    opt_select(iopt_select_ocnatm_CO2) = opt_select(iopt_select_carbchem) .AND. atm_select(ia_pCO2)

    ! *** parameter consistency check - biological productivity ***
    ! check first-order consistency between biologial option, and selected dissolved and sedimentary tracers
    ! NOTE: only the existence of inconsistency will be highlighted, not exactly what the problem is ...
    SELECT CASE (par_bio_prodopt)
    CASE (                        &
         & '1N1T_PO4restore',     &
         & '1N1T_PO4restoreLL',   &
         & '1N1T_PO4MM',          &
         & '1N1T_PO4MM_Tdep',     &
         & '2N1T_PO4MM_SiO2',     &
         & '2N2T_PN_Tdep',        &
         & '3N2T_PNFe_Tdep',      &
         & 'bio_PFe',             &
         & 'bio_PFeSi',           &
         & 'bio_PFeSi_Ridgwell02' &
         & )
       IF (.NOT. ocn_select(io_PO4)) loc_flag = .TRUE.
       IF (.NOT. sed_select(is_POP)) loc_flag = .TRUE.
    end select
    SELECT CASE (par_bio_prodopt)
    CASE (                        &
         & '2N1T_PO4MM_SiO2',     &
         & 'bio_PFeSi',           &
         & 'bio_PFeSi_Ridgwell02' &
         & )
       IF (.NOT. ocn_select(io_SiO2)) loc_flag = .TRUE.
       IF (.NOT. sed_select(is_opal)) loc_flag = .TRUE.
    end select
    SELECT CASE (par_bio_prodopt)
    case (                    &
         & '2N2T_PO4MM_NO3',   &
         & '2N2T_PN_Tdep',     &
         & '3N2T_PNFe_Tdep'    &
         & )
       IF (.NOT. ocn_select(io_NO3)) loc_flag = .TRUE.
       IF (.NOT. ocn_select(io_N2)) loc_flag = .TRUE.
       IF (.NOT. ocn_select(io_NH4)) loc_flag = .TRUE.
       IF (.NOT. sed_select(is_PON)) loc_flag = .TRUE.
    end select
    SELECT CASE (par_bio_prodopt)
    case (                    &
         & '3N2T_PNFe_Tdep'   &
         & )
       IF (.NOT. ocn_select(io_Fe)) loc_flag = .TRUE.
       IF (.NOT. ocn_select(io_FeL)) loc_flag = .TRUE.
       IF (.NOT. ocn_select(io_L)) loc_flag = .TRUE.
       IF (.NOT. sed_select(is_POFe)) loc_flag = .TRUE.
       IF (.NOT. sed_select(is_POM_Fe)) loc_flag = .TRUE.
    end select
    if (loc_flag) then
       CALL sub_report_error( &
            & 'biogem_data','sub_check_par', &
            & 'Your chosen biological option '//trim(par_bio_prodopt)// &
            & ' is not consistent with the selected ocean (gem_config_ocn.par) and/or sediment (gem_config_sed) tracers. '// &
            & 'Go double-check your selected options, because frankly, I cant be bothered to do your job for you.', &
            & 'STOPPING', &
            & (/const_real_null/),.true. &
            & )
       loc_flag = .FALSE.
    end IF
    If (par_bio_prodopt == 'NONE') then
       If (ctrl_data_save_sig_diag) then
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'Selected data saving is redundant in the event of no biological scheme being activated.', &
               & '[ctrl_data_save_sig_diag] HAS BEEN DE-SELECTED; CONTINUING', &
               & (/const_real_null/),.false. &
               & )
          ctrl_data_save_sig_diag = .FALSE.
       end IF
       If (ctrl_data_save_slice_diag) then
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'Selected data saving is redundant in the event of no biological scheme being activated', &
               & '[ctrl_data_save_slice_diag] HAS BEEN DE-SELECTED; CONTINUING', &
               & (/const_real_null/),.false. &
               & )
          ctrl_data_save_slice_diag = .FALSE.
       end IF
    end if
    ! #### ADD CHECKS OF ADDITIONAL BIOLOGICAL OPTIONS HERE ###################################################################### !
    !
    ! ############################################################################################################################ !
    ! check nutrient restoring tracer self-consistency
    SELECT CASE (par_bio_prodopt)
    CASE (                     &
         & '1N1T_PO4restore',  &
         & '1N1T_PO4restoreLL' &
         & )
       IF (.NOT. force_restore_ocn_select(io_PO4)) THEN
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'PO4 restoring MUST be enabled (FILE: configure_forcings_ocn.dat) in conjunction with the '//&
               & '<1 x nutrient, 1 x taxa: PO4 restoring biological production> option', &
               & 'ALTERING INTERNAL PARAMETER VALUE; CONTINUING', &
               & (/const_real_null/),.false. &
               & )
          force_restore_ocn_select(io_PO4) = .TRUE.
          force_flux_ocn_select(io_PO4) = .FALSE.
          force_restore_ocn_sur(io_PO4) = .TRUE.
          force_restore_ocn_k1(io_PO4,:,:) = MAX(goldstein_k1(:,:),n_k)
       END IF
    CASE (               &
         & 'bio_POCflux' &
         & )
       IF (.NOT. force_flux_sed_select(is_POC)) THEN
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'POC flux forcing MUST be enabled (FILE: configure_forcings_sed.dat) in conjunction with the '//&
               & 'POC flux based biological production option', &
               & 'STOPPING', &
               & (/const_real_null/),.true. &
               & )
       END IF
    end select
    ! check that the necessary dissolved organic matter tracers have been selected for each particulate (sed) tracer selected and
    ! de-select all DOM tracers (including dependents) if no DOM production is specified
    if (par_bio_red_DOMfrac > const_real_nullsmall) then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_tot_i = conv_POM_DOM_i(0,is)
          do loc_i=1,loc_tot_i
             io = conv_POM_DOM_i(loc_i,is)
             if (.NOT. ocn_select(io)) THEN
                loc_flag = .TRUE.
                loc_string = string_ocn(io)
             end if
          end do
       end do
       if (loc_flag) then
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'you have rather cheekily set a non-zero fraction of dissolved organic matter production '//&
               & 'but have failed to ensure that you have the necessary DOM tracers selected - '//TRIM(loc_string)// &
               & '[HINT: there must be a corresponding dissolved tracer for each particulate tracer selected '// &
               & 'Sadly, it is too late to automatically select '//TRIM(loc_string)// &
               & ' and a bit risky to set DOM to zero for you :(', &
               & 'STOPPING', &
               & (/const_real_null/),.true. &
               & )
          loc_flag = .FALSE.
       end if
    else
       if (loc_flag) then
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'although the production of dissolved organic matter production is set to zero '//&
               & 'you have carelessly left some DOM tracers selected (FILE: gem_config_ocn.par) '// &
               & 'The model will run quicker by de-selecting all currently selected DOM tracers.', &
               & 'CONTINUING', &
               & (/const_real_null/),.FALSE. &
               & )
          loc_flag = .FALSE.
       end if
    end if

    ! *** parameter consistency check - isotopes, forcings ***
    ! OCEAN TRACERS
    do io=1,n_ocn
       IF (ocn_select(io)) THEN
          if (.not. ocn_select(ocn_dep(io))) then
             loc_string = string_ocn(ocn_dep(io))
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'If an isotope tracer is selected, the associated bulk ocean tracer '//TRIM(loc_string)// &
                  & ' must be selected', &
                  & 'STOPPING', &
                  & (/const_real_null/),.true. &
                  & )
          end if
          if (ocn_select(ocn_dep(io)) .AND. (io /= ocn_dep(io))) then
             if ( &
                  & (force_restore_ocn_select(io) .AND. (.NOT. force_restore_ocn_select(ocn_dep(io)))) &
                  & .OR. &
                  & (.NOT. (force_restore_ocn_select(io)) .AND. force_restore_ocn_select(ocn_dep(io))) &
                  & ) then
                loc_string1 = string_ocn(io)
                loc_string2 = string_ocn(ocn_dep(io))
                CALL sub_report_error( &
                     & 'biogem_data','sub_check_par', &
                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// &
                     & ' have been selected, but a restoring forcing for only one of them has been selected.', &
                     & 'CONTINUING', &
                     & (/const_real_null/),.false. &
                     & )
             end if
             if ( &
                  & (force_flux_ocn_select(io) .AND. (.NOT. force_flux_ocn_select(ocn_dep(io)))) &
                  & .OR. &
                  & (.NOT. (force_flux_ocn_select(io)) .AND. force_flux_ocn_select(ocn_dep(io))) &
                  & ) then
                loc_string1 = string_ocn(io)
                loc_string2 = string_ocn(ocn_dep(io))
                CALL sub_report_error( &
                     & 'biogem_data','sub_check_par', &
                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// &
                     & ' have been selected, but a flux forcing for only one of them has been selected.', &
                     & 'CONTINUING', &
                     & (/const_real_null/),.false. &
                     & )
             end if
          end if
       else
          if (force_restore_ocn_select(io)) then
             loc_string = string_ocn(io)
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'Because the ocean tracer '//TRIM(loc_string)//' has not been selected, '// &
                  & 'a restoring forcing of this tracer cannot be performed', &
                  & 'CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
             force_restore_ocn_select(io) = .FALSE.
          end if
          if (force_flux_ocn_select(io)) then
             loc_string = string_ocn(io)
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'Because the ocean tracer '//TRIM(loc_string)//' has not been selected, '// &
                  & 'a flux forcing of this tracer cannot be performed', &
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
             force_flux_ocn_select(io) = .FALSE.
          end if
       end if
    end do
    ! ATMOSPHERE TRACERS
    do ia=1,n_atm
       IF (atm_select(ia)) THEN
          if (.not. atm_select(atm_dep(ia))) then
             loc_string = string_atm(atm_dep(ia))
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'If an isotope tracer is selected, the associated bulk atmosphere tracer '//TRIM(loc_string)// &
                  & ' must be selected', &
                  & 'STOPPING', &
                  & (/const_real_null/),.true. &
                  & )
          end if
          if (atm_select(atm_dep(ia)) .AND. (io /= atm_type(ia))) then
             if ( &
                  & (force_restore_atm_select(ia) .AND. (.NOT. force_restore_atm_select(atm_dep(ia)))) &
                  & .OR. &
                  & (.NOT. (force_restore_atm_select(ia)) .AND. force_restore_atm_select(atm_dep(ia))) &
                  & ) then
                loc_string1 = string_atm(ia)
                loc_string2 = string_atm(atm_dep(ia))
                CALL sub_report_error( &
                     & 'biogem_data','sub_check_par', &
                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// &
                     & ' have been selected, but a restoring forcing for only one of them has been selected.', &
                     & 'CONTINUING', &
                     & (/const_real_null/),.false. &
                     & )
             end if
             if ( &
                  & (force_flux_atm_select(ia) .AND. (.NOT. force_flux_atm_select(atm_dep(ia)))) &
                  & .OR. &
                  & (.NOT. (force_flux_atm_select(ia)) .AND. force_flux_atm_select(atm_dep(ia))) &
                  & ) then
                loc_string1 = string_atm(ia)
                loc_string2 = string_atm(atm_dep(ia))
                If ((ia == ia_pCO2_14C) .AND. (atm_dep(ia) == ia_pCO2)) then
                   CALL sub_report_error( &
                        & 'biogem_data','sub_check_par', &
                        & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// &
                        & ' have been selected, but a flux forcing for only one of them has been selected.', &
                        & 'CONTINUING', &
                        & (/const_real_null/),.false. &
                        & )
                else
                   CALL sub_report_error( &
                        & 'biogem_data','sub_check_par', &
                        & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// &
                        & ' have been selected, but a flux forcing for only one of them has been selected.', &
                        & 'CONTINUING', &
                        & (/const_real_null/),.false. &
                        & )
                end If
             end if
          end if
       else
          if (force_restore_atm_select(ia)) then
             loc_string = string_atm(ia)
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'Because the atmospheric tracer '//TRIM(loc_string)//' has not been selected, '// &
                  & 'a restoring forcing of this tracer cannot be performed', &
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
             force_restore_atm_select(ia) = .FALSE.
          end if
          if (force_flux_atm_select(ia)) then
             loc_string = string_atm(ia)
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'Because the atmospheric tracer '//TRIM(loc_string)//' has not been selected, '// &
                  & 'a flux forcing of this tracer cannot be performed', &
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
             force_flux_atm_select(ia) = .FALSE.
          end if
       end IF
    end do
    ! SEDIMENT TRACERS
    do is=1,n_sed
       IF (sed_select(is)) THEN
          if (.not. sed_select(sed_dep(is))) then
             loc_string1 = string_sed(is)
             loc_string2 = string_sed(sed_dep(is))
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'If an isotope or other dependent tracer  '//TRIM(loc_string1)//' is selected, '// &
                  & 'the associated bulk sediment tracer '//TRIM(loc_string), &
                  & 'STOPPING', &
                  & (/const_real_null/),.true. &
                  & )
          end if
          if (sed_select(sed_dep(is)) .AND. (is /= sed_dep(is)) .AND. (sed_type(is) /= par_sed_type_scavenged)) then
             if ( &
                  & (force_flux_sed_select(is) .AND. (.NOT. force_flux_sed_select(sed_dep(is)))) &
                  & .OR. &
                  & (.NOT. (force_flux_sed_select(is)) .AND. force_flux_sed_select(sed_dep(is))) &
                  & ) then
                loc_string1 = string_sed(is)
                loc_string2 = string_sed(sed_dep(is))
                CALL sub_report_error( &
                     & 'biogem_data','sub_check_par', &
                     & 'An isotope tracer '//TRIM(loc_string1)//' and associated bulk tracer '//TRIM(loc_string2)// &
                     & ' have been selected, but a flux forcing for only one of them has been selected.', &
                     & 'CONTINUING', &
                     & (/const_real_null/),.false. &
                     & )
             end if
          end if
       else
          if (force_flux_sed_select(is)) then
             loc_string = string_sed(is)
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'Because the sediment tracer '//TRIM(loc_string)//' has not been selected, '// &
                  & 'a flux forcing of this tracer cannot be performed', &
                  & 'RESTORING HAS BEEN DE-SELECTED; CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
             force_flux_sed_select(is) = .FALSE.
          end if
       end IF
    end do

    ! *** FIX UP AND MAKE GENERIC ***
    ! verify ocn-atm carbon cycle option selection
    IF (atm_select(ia_pCO2) .NEQV. ocn_select(io_DIC)) THEN
       CALL sub_report_error( &
            & 'biogem_data','sub_check_par', &
            & 'Do you really mean to select CO2 in the atmosphere but not in the ocean (or vice versa)?', &
            & 'CONTINUING', &
            & (/const_real_null/),.false. &
            & )
    ENDIF

    ! *** parameter consistency check - selected tracers and sediment-sediment option combinations ***
    do is=1,n_sed
       if (sed_type(is) == par_sed_type_frac) then
          if (( .NOT. sed_select(is)) .AND. (sed_select(sed_dep(is)))) then
             loc_string2 = string_sed(is)
             CALL sub_report_error( &
                  & 'biogem_data','sub_check_par', &
                  & 'A frac2 tracer must be selected associated with '//TRIM(loc_string2), &
                  & 'STOPPING', &
                  & (/const_real_null/),.true. &
                  & )
          end if
       end if
    end do
    IF (sed_select(is_CaCO3) .AND. (.NOT. sed_select(is_POC))) THEN
       CALL sub_report_error( &
            & 'biogem_data','sub_check_par','The POC tracer must be selected with CaCO3 in biogem_config_sed.par ', &
            & 'STOPPING', &
            & (/const_real_null/),.true. &
            & )
    ENDIF
    If (sed_select(is_CaCO3_age) .AND. (.NOT. sed_select(is_CaCO3))) then
       CALL sub_report_error( &
            & 'biogem_data','sub_check_par', &
            & 'If the sediment CaCO3 age tracer is requested, then the solid CaCO3 tracer must be selected ', &
            & 'STOPPING', &
            & (/const_real_null/),.true. &
            & )
    ENDIF

    ! *** parameter consistency check - selected sediment-ocean tracer option combinations ***
    if (par_bio_prodopt /= 'NONE') then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          select case (sed_type(is))
          case (par_sed_type_bio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_scavenged,11:20)
             loc_tot_i = conv_sed_ocn_i(0,is)
             do loc_i=1,loc_tot_i
                io = conv_sed_ocn_i(loc_i,is)
                if (abs(conv_sed_ocn(io,is)) > const_real_nullsmall) then
                   if (.NOT. ocn_select(io)) then
                      loc_string1 = string_ocn(io)
                      loc_string2 = string_sed(is)
                      CALL sub_report_error( &
                           & 'biogem_data','sub_check_par', &
                           & 'Particulate tracer '//TRIM(loc_string2)// &
                           & ' does does not have the corresponding ocean tracer '//TRIM(loc_string1)//' selected', &
                           & 'CONTINUING', &
                           & (/const_real_null/),.false. &
                           & )
                   end if
                end if
             end do
          end SELECT
       end DO
    end if

    ! *** parameter consistency check - 14C ***
    ! NOTE: just test ocean tracers

    IF (ocn_select(io_DIC_14C) .AND. (.NOT. ocn_select(io_DIC_13C))) THEN
       CALL sub_report_error( &
            & 'biogem_data','sub_check_par', &
            & 'To select 14C isotope tracers, 13C isotope tracers MUST also be selected', &
            & 'STOPPING', &
            & (/const_real_null/),.true. &
            & )
    end if

    ! *** parameter consistency check - data save options ***
    IF (.NOT. opt_select(iopt_select_carbchem)) THEN
       IF (ctrl_data_save_sig_carb_sur) THEN
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'You do not have sufficent ocean tracers selected for a marine carbon cycle', &
               & '[ctrl_data_save_sig_carb_sur] HAS BEEN DE-SELECTED; CONTINUING', &
               & (/const_real_null/),.false. &
               & )
          ctrl_data_save_sig_carb_sur = .FALSE.
       end if
       If (ctrl_data_save_slice_carb) then
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'You do not have sufficent ocean tracers selected for a marine carbon cycle', &
               & '[ctrl_data_save_slice_carb] HAS BEEN DE-SELECTED; CONTINUING', &
               & (/const_real_null/),.false. &
               & )
          ctrl_data_save_slice_carb = .FALSE.
       end if
       If (ctrl_data_save_slice_carbconst) then
          CALL sub_report_error( &
               & 'biogem_data','sub_check_par', &
               & 'You do not have sufficent ocean tracers selected for a marine carbon cycle', &
               & '[ctrl_data_save_slice_carbconst] HAS BEEN DE-SELECTED; CONTINUING', &
               & (/const_real_null/),.false. &
               & )
          ctrl_data_save_slice_carbconst = .FALSE.
       end IF
    end IF

  END SUBROUTINE sub_check_par_biogem
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE CARBONATE SYSTEM
  SUBROUTINE sub_init_carb()
    ! local variables
    INTEGER::i,j,k
    ! zero arrays
    ! NOTE: leave carb_TSn array at its initialized state
    !       so that a full update of carb constants etc is ALWAYS performed upon the first call to tstep_biogem
    carbconst(:,:,:,:) = 0.0
    carb(:,:,:,:)      = 0.0
    carbalk(:,:,:,:)   = 0.0
    carbisor(:,:,:,:)  = 0.0
    carb_TSn(:,:,:,:)  = 0.0
    ! initialize arrays
    DO i=1,n_i
       DO j=1,n_j
          DO k=goldstein_k1(i,j),n_k
             ! calculate carbonate constants
             CALL sub_calc_carbconst(         &
                  & phys_ocn(ipo_Dmid,i,j,k), &
                  & ocn(io_T,i,j,k),          &
                  & ocn(io_S,i,j,k),          &
                  & carbconst(:,i,j,k)        &
                  & )
             ! adjust carbonate constants
             if (ocn_select(io_Ca) .AND. ocn_select(io_Mg)) then
                call sub_adj_carbconst(   &
                     & ocn(io_Ca,i,j,k),  &
                     & ocn(io_Mg,i,j,k),  &
                     & carbconst(:,i,j,k) &
                     & )
             END if
             IF (opt_select(iopt_select_carbchem)) then
                ! estimate Ca and borate concentrations (if not selected and therefore explicitly treated)
                IF (.NOT. ocn_select(io_Ca))  ocn(io_Ca,i,j,k)  = fun_calc_Ca(ocn(io_S,i,j,k))
                IF (.NOT. ocn_select(io_B))   ocn(io_B,i,j,k)   = fun_calc_Btot(ocn(io_S,i,j,k))
                IF (.NOT. ocn_select(io_SO4)) ocn(io_SO4,i,j,k) = fun_calc_SO4tot(ocn(io_S,i,j,k))
                IF (.NOT. ocn_select(io_F))   ocn(io_F,i,j,k)   = fun_calc_Ftot(ocn(io_S,i,j,k))
                ! seed default initial ocean pH
                carb(ic_H,i,j,k) = 10**(-7.8)
                ! calculate carbonate chemistry
                CALL sub_calc_carb(        &
                     & ocn(io_DIC,i,j,k),  &
                     & ocn(io_ALK,i,j,k),  &
                     & ocn(io_Ca,i,j,k),   &
                     & ocn(io_PO4,i,j,k),  &
                     & ocn(io_SiO2,i,j,k), &
                     & ocn(io_B,i,j,k),    &
                     & ocn(io_SO4,i,j,k),  &
                     & ocn(io_F,i,j,k),    &
                     & ocn(io_H2S,i,j,k),  &
                     & ocn(io_NH4,i,j,k),  &
                     & carbconst(:,i,j,k), & 
                     & carb(:,i,j,k),      &
                     & carbalk(:,i,j,k)    &
                     & )
                ! estimate Revelle factor
                CALL sub_calc_carb_RF0(      &
                     & ocn(io_DIC,i,j,k),  &
                     & ocn(io_ALK,i,j,k),  &
                     & ocn(io_PO4,i,j,k),  &
                     & ocn(io_SiO2,i,j,k), &
                     & ocn(io_B,i,j,k),    &
                     & ocn(io_SO4,i,j,k),  &
                     & ocn(io_F,i,j,k),    &
                     & ocn(io_H2S,i,j,k),  &
                     & ocn(io_NH4,i,j,k),  &
                     & carbconst(:,i,j,k), & 
                     & carb(:,i,j,k)    & 
                     & )
                ! calculate carbonate system isotopic properties
                if (ocn_select(io_DIC_13C)) then
                   call sub_calc_carb_r13C(      &
                        & ocn(io_T,i,j,k),       &
                        & ocn(io_DIC,i,j,k),     &
                        & ocn(io_DIC_13C,i,j,k), &
                        & carb(:,i,j,k),         &
                        & carbisor(:,i,j,k)      &
                        & )
                end IF
                if (ocn_select(io_DIC_14C)) then
                   call sub_calc_carb_r14C(      &
                        & ocn(io_T,i,j,k),       &
                        & ocn(io_DIC,i,j,k),     &
                        & ocn(io_DIC_14C,i,j,k), &
                        & carb(:,i,j,k),         &
                        & carbisor(:,i,j,k)      &
                        & )
                end IF
             end if
          END DO
       END DO
    END DO
  END SUBROUTINE sub_init_carb
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE SOLUBILITY CONSTANTS
  SUBROUTINE sub_init_solconst()
    ! local variables
    INTEGER::i,j
    ! zero arrays
    ocnatm_airsea_solconst(:,:,:) = 0.0
    ! initialize array
    DO i=1,n_i
       DO j=1,n_j
          if (n_k >= goldstein_k1(i,j)) then
             call sub_calc_solconst(i,j)
          end if
       END DO
    END DO
  END SUBROUTINE sub_init_solconst
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE DATA SAVING
  SUBROUTINE sub_init_data_save()
    ! local variables
    INTEGER::n
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_n_elements
    INTEGER::loc_sig_i
    INTEGER::loc_timeslice_i
    real::loc_data_scale

    ! *** set time series data save interval details ***
    ! initialize time series indices
    par_data_save_sig_i = n_data_max
    par_data_save_sig(:) = 0.0
    ! load data
    loc_filename = TRIM(par_indir_name)//TRIM(par_infile_sig_name)
    loc_data_scale = 1.0
    CALL sub_load_data_t1(loc_filename,loc_data_scale,par_data_save_sig,loc_n_elements)
    ! if no elements, populate array with default time interval steps
    IF (loc_n_elements == 0) THEN
       ! limit the time-series integration interval
       if (par_data_save_sig_dt > const_real_nullsmall) then
          loc_n_elements = INT(par_misc_t_runtime/par_data_save_sig_dt + const_real_nullsmall)
          do while (loc_n_elements > n_data_max)
             par_data_save_sig_dt = 10.0*par_data_save_sig_dt
             loc_n_elements = INT(par_misc_t_runtime/par_data_save_sig_dt + const_real_nullsmall)
             CALL sub_report_error( &
                  & 'biogem_data','sub_init_data_save','time-series save interval (biogem_config.par) too short - '// &
                  & 'was [lower value] and is now [upper value] (years)', &
                  & 'CONTINUING', &
                  & (/par_data_save_sig_dt,par_data_save_sig_dt/10.0/),.FALSE. &
                  & )
          end do
          DO n=1,loc_n_elements
             par_data_save_sig(n) = &
                  & real(n - 0.5)*par_data_save_sig_dt + (par_misc_t_runtime - real(loc_n_elements)*par_data_save_sig_dt)
          END DO
       else
          CALL sub_report_error( &
               & 'biogem_data','sub_init_data_save','time-series save interval (biogem_config.par) '// &
               & 'must be non-zero and positive', &
               & 'STOPPING', &
               & (/const_real_null/),.TRUE. &
               & )
       endif
    end IF
    ! find first save time lying within total model run-time
    ! NOTE: <loc_sig_i> will be zero if no valid time points have been requested in the time series input file,
    !       and the array has not been populated automatically
    ! NOTE: ensure that the first identified time-series time is at least a full integration interval (required value)
    !       from the start time of the model run
    loc_sig_i = loc_n_elements
    DO while (loc_sig_i > 0)
       IF ( &
            & par_data_save_sig(loc_sig_i) &
            & < &
            & (par_misc_t_runtime - par_data_save_sig_dt/2.0 + par_misc_t_err) &
            & ) THEN
          EXIT
       ELSE
          loc_sig_i = loc_sig_i - 1
       END IF
    END DO
    par_data_save_sig_i = loc_sig_i

    ! *** set time slice data save details ***
    ! NOTE: DO NOT populate the time-slice array automatically if the data file is empty
    ! initialize time slice indices
    par_data_save_timeslice_i = n_data_max
    par_data_save_timeslice(:) = 0.0
    ! load data
    loc_filename = TRIM(par_indir_name)//TRIM(par_infile_slice_name)
    loc_data_scale = 1.0
    CALL sub_load_data_t1(loc_filename,loc_data_scale,par_data_save_timeslice,loc_n_elements)
    ! find first save time lying within total model run-time
    ! NOTE: <par_data_save_timeslice_i> will be zero if no valid time slices have been requested in the time slice input file
    ! NOTE: ensure that the first identified time-slice time is at least a full integration interval (required value)
    !       from the start time of the model run
    loc_timeslice_i = loc_n_elements
    DO while (loc_timeslice_i > 0)
       IF ( &
            & par_data_save_timeslice(loc_timeslice_i) &
            & < &
            & (par_misc_t_runtime - par_data_save_slice_dt/2.0 + par_misc_t_err) &
            & ) THEN
          EXIT
       ELSE
          loc_timeslice_i = loc_timeslice_i - 1
       END IF
    END DO
    if (par_data_save_timeslice(loc_timeslice_i) < (par_data_save_slice_dt/2.0 - par_misc_t_err)) then
       loc_timeslice_i = 0
    end if
    par_data_save_timeslice_i = loc_timeslice_i
    if (par_data_save_timeslice_i == 0) then
       CALL sub_report_error( &
            & 'biogem_data','sub_init_data_save', &
            & 'No time-slice dates listed in file biogem_save_timeslice.dat fall within the model start and end years', &
            & 'CONTINUING', &
            & (/const_real_null/),.false. &
            & )
    end if
  END SUBROUTINE sub_init_data_save
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE RESTORING FORCING - OCEAN
  SUBROUTINE sub_init_force_restore_ocn()
    ! local variables
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_n_elements
    integer::l,i,j,k,io
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk
    real,DIMENSION(2)::loc_data_scale
    ! LOOP
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       IF (force_restore_ocn_select(io)) THEN
          force_restore_ocn_sig_i(io,:) = n_data_max
          force_restore_ocn_sig(io,:,:) = 0.0
          ! load forcing data array #I
          loc_ijk(:,:,:) = const_real_zero
          if (force_ocn_uniform(io) == 3) then
             loc_ijk(:,:,:) = 0.0
          elseif (force_ocn_uniform(io) == 2) then
             loc_ijk(:,:,n_k) = 0.0
          elseif (force_ocn_uniform(io) == 0) then
             loc_ijk(force_ocn_point_i(io),force_ocn_point_j(io),force_ocn_point_k(io)) = 0.0
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_I'//TRIM(string_data_ext)
             CALL sub_load_data_ijk(loc_filename,n_i,n_j,n_k,loc_ijk(:,:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                DO k=force_restore_ocn_k1(io,i,j),n_k
                   force_restore_ocn_I(io,i,j,k) = loc_ijk(i,j,k)
                end do
             end do
          end DO
          ! load forcing data array #II
          loc_ijk(:,:,:) = const_real_zero
          if (force_ocn_uniform(io) == 3) then
             loc_ijk(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
          elseif (force_ocn_uniform(io) == 2) then
             loc_ijk(:,:,n_k) = phys_ocn(ipo_mask_ocn,:,:,n_k)
          elseif (force_ocn_uniform(io) == 0) then
             loc_ijk(force_ocn_point_i(io),force_ocn_point_j(io),force_ocn_point_k(io)) = 1.0
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_II'//TRIM(string_data_ext)
             CALL sub_load_data_ijk(loc_filename,n_i,n_j,n_k,loc_ijk(:,:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                DO k=force_restore_ocn_k1(io,i,j),n_k
                   force_restore_ocn_II(io,i,j,k) = loc_ijk(i,j,k)
                end do
             end do
          end DO
          ! load forcing time series data
          loc_filename =TRIM(par_fordir_name)//'biogem_force_restore_ocn_'//TRIM(string_ocn(io))//'_sig'//TRIM(string_data_ext)
          loc_data_scale(:) = (/par_ocn_force_scale_time(io), par_ocn_force_scale_val(io)/)
          CALL sub_load_data_t2(loc_filename,loc_data_scale(:),force_restore_ocn_sig(io,:,:),loc_n_elements)
          ! set default forcing index values
          ! NOTE: catch missing time series data
          if (loc_n_elements == 0) THEN
             CALL sub_report_error( &
                  & 'biogem_data','init_force_restore_ocn','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), &
                  & 'STOPPING', &
                  & (/const_real_null/),.TRUE. &
                  & )
          else
             force_restore_ocn_sig_i(io,:) = loc_n_elements
          end if
          ! warn if forcing information appears to be 'incomplete'
          ! NOTE: this will catch both possible mismatches of forcing signal and model integration specification
          !       i.e., for _not_ the BP option;
          !             a signal time start year that is after the model start year
          !             or a signal time year that is before the model end year
          !             (or visa versa for a BP time scale)
          if ( &
               & (minval(force_restore_ocn_sig(io,1,1:loc_n_elements)) > 0.0) &
               & .OR. &
               & (maxval(force_restore_ocn_sig(io,1,1:loc_n_elements)) < par_misc_t_runtime) &
               & ) then
             CALL sub_report_error( &
                  & 'biogem_data','sub_init_force_restore_ocn', &
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// &
                  & '(1) alter the model start time year and/or run length (the goin file) or'// &
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', &
                  & 'CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
          end if
       end IF
    end DO
  END SUBROUTINE sub_init_force_restore_ocn
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE RESTORING FORCING - ATMOSPHERE
  SUBROUTINE sub_init_force_restore_atm()
    ! local variables
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_n_elements
    integer::l,i,j,ia
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(2)::loc_data_scale
    ! LOOP
    DO l=3,n_l_atm
       ia = conv_iselected_ia(l)
       IF (force_restore_atm_select(ia)) THEN
          force_restore_atm_sig_i(ia,:) = n_data_max
          force_restore_atm_sig(ia,:,:) = 0.0
          ! load forcing data array #I
          loc_ij(:,:) = const_real_zero
          if (force_atm_uniform(ia) == 2) then
             loc_ij(:,:) = 0.0
          elseif (force_atm_uniform(ia) == 0) then
             loc_ij(force_atm_point_i(ia),force_atm_point_j(ia)) = 0.0
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_I'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   force_restore_atm_I(ia,i,j) = loc_ij(i,j)
                end IF
             end DO
          end DO
          ! load forcing data array #II
          loc_ij(:,:) = const_real_zero
          if (force_atm_uniform(ia) == 2) then
             loc_ij(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_k)
          elseif (force_atm_uniform(ia) == 0) then
             loc_ij(force_atm_point_i(ia),force_atm_point_j(ia)) = 1.0
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_II'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   force_restore_atm_II(ia,i,j) = loc_ij(i,j)
                end IF
             end DO
          end DO
          ! load forcing time series data
          loc_filename = TRIM(par_fordir_name)//'biogem_force_restore_atm_'//TRIM(string_atm(ia))//'_sig'//TRIM(string_data_ext)
          loc_data_scale(:) = (/par_atm_force_scale_time(ia), par_atm_force_scale_val(ia)/)
          CALL sub_load_data_t2(loc_filename,loc_data_scale(:),force_restore_atm_sig(ia,:,:),loc_n_elements)
          ! set default forcing index values
          ! NOTE: catch missing time series data
          if (loc_n_elements == 0) THEN
             CALL sub_report_error( &
                  & 'biogem_data','init_force_restore_atm','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), &
                  & 'STOPPING', &
                  & (/const_real_null/),.TRUE. &
                  & )
          else
             force_restore_atm_sig_i(ia,:) = loc_n_elements
          end if
          ! warn if forcing information appears to be 'incomplete'
          ! NOTE: this will catch both possible mismatches of forcing signal and model integration specification
          !       i.e., for _not_ the BP option;
          !             a signal time start year that is after the model start year
          !             or a signal time year that is before the model end year
          !             (or visa versa for a BP time scale)
          if ( &
               & (minval(force_restore_atm_sig(ia,1,1:loc_n_elements)) > 0.0) &
               & .OR. &
               & (maxval(force_restore_atm_sig(ia,1,1:loc_n_elements)) < par_misc_t_runtime) &
               & ) then
             CALL sub_report_error( &
                  & 'biogem_data','sub_init_force_restore_atm', &
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// &
                  & '(1) alter the model start time year and/or run length (the goin file) or'// &
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', &
                  & 'CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
          end if
       end IF
    end DO
  END SUBROUTINE sub_init_force_restore_atm
  ! ****************************************************************************************************************************** !


!!$  ! INITIALIZE RESTORING FORCING - SEDIMENTS
!!$  SUBROUTINE sub_init_force_restore_sed(dum_is)
!!$    ! dummy arguments
!!$    INTEGER,INTENT(in)::dum_is
!!$    ! local varisbles
!!$    CHARACTER(len=255)::loc_filename
!!$    INTEGER::loc_n_elements
!!$    ! initislize forcing signal indices
!!$    force_restore_sed_sig_i(dum_is,:) = n_data_max
!!$    force_restore_sed_sig(dum_is,:,:) = 0.0
!!$    ! load forcing data array #I
!!$    loc_filename = TRIM(par_indir_name)//'biogem_force_restore_sed_'//TRIM(string_sed(dum_is))//'_I'//TRIM(string_data_ext)
!!$    CALL sub_load_data_ij(loc_filename,n_i,n_j,force_restore_sed_I(dum_is,:,:))
!!$    ! load forcing data array #II
!!$    loc_filename = TRIM(par_indir_name)//'biogem_force_restore_sed_'//TRIM(string_sed(dum_is))//'_II'//TRIM(string_data_ext)
!!$    CALL sub_load_data_ij(loc_filename,n_i,n_j,force_restore_sed_II(dum_is,:,:))
!!$    ! load forcing time series data
!!$    loc_filename = TRIM(par_indir_name)//'biogem_force_restore_sed_'//TRIM(string_sed(dum_is))//'_sig'//TRIM(string_data_ext)
!!$    CALL sub_load_data_t2(loc_filename,force_restore_sed_sig(dum_is,:,:),loc_n_elements)
!!$    ! set default forcing index values
!!$    ! NOTE: catch missing time series data
!!$    if (loc_n_elements == 0) THEN
!!$       CALL sub_report_error( &
!!$            & 'biogem_data','init_force_restore_sed','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), &
!!$            & 'STOPPING', &
!!$            & (/const_real_null/),.TRUE. &
!!$            & )
!!$    else
!!$       force_restore_sed_sig_i(dum_is,:) = loc_n_elements
!!$    end if
!!$    ! warn if forcing information appears to be 'incomplete'
!!$    ! NOTE: this will catch both possible mismatches of forcing signal and model integration specification
!!$    !       i.e., for _not_ the BP option;
!!$    !             a signal time start year that is after the model start year
!!$    !             or a signal time year that is before the model end year
!!$    !             (or visa versa for a BP time scale)
!!$    if (maxval(force_restore_sed_sig(dum_is,1,1:loc_n_elements)) < par_misc_t_runtime) then
!!$       CALL sub_report_error( &
!!$            & 'biogem_data','sub_init_force_restore_sed', &
!!$            & 'The time interval of forcing data does not fully span the requested model integration, either; '// &
!!$            & '(1) alter the model start time year (the goin file) '// &
!!$            & '(2) alter the forcing signal start time year (FILE: '//TRIM(loc_filename)//') '// &
!!$            & '(3) leave everything well alone '// &
!!$            & '(it is legitamite to start a model forcing part way through a run by defining a partial time signal)', &
!!$            & 'CONTINUING', &
!!$            & (/const_real_null/),.false. &
!!$            & )
!!$    end if
!!$    if (minval(force_restore_sed_sig(dum_is,1,1:loc_n_elements)) > 0.0) then
!!$       CALL sub_report_error( &
!!$            & 'biogem_data','sub_init_force_restore_sed', &
!!$            & 'The time interval of forcing data does not fully span the requested model integration, either; '// &
!!$            & '(1) alter the model start time year and/or run length (the goin file) '// &
!!$            & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') '// &
!!$            & '(3) leave everything well alone '// &
!!$            & '(it is legitamite to stop a model forcing part way through a run by defining a partial time signal)', &
!!$            & 'CONTINUING', &
!!$            & (/const_real_null/),.false. &
!!$            & )
!!$    end if
!!$  END SUBROUTINE sub_init_force_restore_sed


  ! ****************************************************************************************************************************** !
  ! INITIALIZE FLUX FORCING - OCEAN
  SUBROUTINE sub_init_force_flux_ocn()
    ! local variables
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_n_elements
    integer::l,i,j,k,io
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk
    real,DIMENSION(2)::loc_data_scale
    ! LOOP
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       IF (force_flux_ocn_select(io)) THEN
          force_flux_ocn_sig_i(io,:) = n_data_max
          force_flux_ocn_sig(io,:,:) = 0.0
          ! load forcing data array #I
          loc_ijk(:,:,:) = const_real_zero
          if (force_ocn_uniform(io) == 3) then
             loc_ijk(:,:,:) = 0.0
          elseif (force_ocn_uniform(io) == 2) then
             loc_ijk(:,:,n_k) = 0.0
          elseif (force_ocn_uniform(io) == 0) then
             loc_ijk(force_ocn_point_i(:),force_ocn_point_j(:),force_ocn_point_k(:)) = 0.0
          elseif (force_ocn_uniform(io) == -1) then
             loc_ijk(:,:,n_k) = 0.0
          elseif (force_ocn_uniform(io) == -2) then
             DO i=1,n_i
                DO j=1,n_j
                   if (goldstein_k1(i,j) <= n_k) then
                      loc_ijk(:,:,goldstein_k1(i,j)) = 0.0
                   end if
                end DO
             end DO
          elseIF ((par_force_point_i > 0) .AND. (par_force_point_j > 0) .AND. (par_force_point_k > 0)) then
             loc_ijk(force_ocn_point_i(:),force_ocn_point_j(:),force_ocn_point_k(:)) = 0.0
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_I'//TRIM(string_data_ext)
             CALL sub_load_data_ijk(loc_filename,n_i,n_j,n_k,loc_ijk(:,:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                DO k=force_flux_ocn_k1(io,i,j),n_k
                   force_flux_ocn_I(io,i,j,k) = loc_ijk(i,j,k)
                end do
             end do
          end DO
          ! load forcing data array #II
          loc_ijk(:,:,:) = const_real_zero
          if (force_ocn_uniform(io) == 3) then
             loc_ijk(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
          elseif (force_ocn_uniform(io) == 2) then
             loc_ijk(:,:,n_k) = phys_ocn(ipo_mask_ocn,:,:,n_k)
          elseif (force_ocn_uniform(io) == 0) then
             loc_ijk(force_ocn_point_i(io),force_ocn_point_j(io),force_ocn_point_k(io)) = 1.0
          elseif (force_ocn_uniform(io) == -1) then
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_SUR'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
             DO i=1,n_i
                DO j=1,n_j
                   if (goldstein_k1(i,j) <= n_k) then
                      loc_ijk(i,j,n_k) = loc_ij(i,j)
                   end if
                end do
             end DO
          elseif (force_ocn_uniform(io) == -2) then
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_BEN'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
             DO i=1,n_i
                DO j=1,n_j
                   if (goldstein_k1(i,j) <= n_k) then
                      loc_ijk(i,j,goldstein_k1(i,j)) = loc_ij(i,j)
                   end if
                end do
             end DO
          elseIF ((par_force_point_i > 0) .AND. (par_force_point_j > 0) .AND. (par_force_point_k > 0)) then
             loc_ijk(force_ocn_point_i(io),force_ocn_point_j(io),force_ocn_point_k(io)) = 1.0
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_II'//TRIM(string_data_ext)
             CALL sub_load_data_ijk(loc_filename,n_i,n_j,n_k,loc_ijk(:,:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                DO k=force_flux_ocn_k1(io,i,j),n_k
                   force_flux_ocn_II(io,i,j,k) = loc_ijk(i,j,k)
                end do
             end do
          end DO
          ! load forcing time series data
          loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_ocn_'//TRIM(string_ocn(io))//'_sig'//TRIM(string_data_ext)
          loc_data_scale(:) = (/par_ocn_force_scale_time(io), par_ocn_force_scale_val(io)/)
          CALL sub_load_data_t2(loc_filename,loc_data_scale(:),force_flux_ocn_sig(io,:,:),loc_n_elements)
          ! set default forcing index values
          ! NOTE: catch missing time series data
          if (loc_n_elements == 0) THEN
             CALL sub_report_error( &
                  & 'biogem_data','init_force_flux_ocn','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), &
                  & 'STOPPING', &
                  & (/const_real_null/),.TRUE. &
                  & )
          else
             force_flux_ocn_sig_i(io,:) = loc_n_elements
          end if
          ! warn if forcing information appears to be 'incomplete'
          ! NOTE: this will catch both possible mismatches of forcing signal and model integration specification
          !       i.e., for _not_ the BP option;a signal time start year that is after the model start year
          !             or a signal time year that is before the model end year (or visa versa for a BP time scale)
          if ( &
               & (minval(force_flux_ocn_sig(io,1,1:loc_n_elements)) > 0.0) &
               & .OR. &
               & (maxval(force_flux_ocn_sig(io,1,1:loc_n_elements)) < par_misc_t_runtime) &
               & ) then
             CALL sub_report_error( &
                  & 'biogem_data','sub_init_force_flux_ocn', &
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// &
                  & '(1) alter the model start time year and/or run length (the goin file) or '// &
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', &
                  & 'CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
          end if
       end IF
    end DO
  END SUBROUTINE sub_init_force_flux_ocn
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE FLUX FORCING - ATMOSPHERE
  SUBROUTINE sub_init_force_flux_atm()
    ! local variables
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_n_elements
    integer::l,i,j,ia
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(2)::loc_data_scale
    ! LOOP
    DO l=3,n_l_atm
       ia = conv_iselected_ia(l)
       IF (force_flux_atm_select(ia)) THEN
          force_flux_atm_sig_i(ia,:) = n_data_max
          force_flux_atm_sig(ia,:,:) = 0.0
          ! load forcing data array #I
          loc_ij(:,:) = const_real_zero
          if (force_atm_uniform(ia) == 2) then
             loc_ij(:,:) = 0.0
          elseif (force_atm_uniform(ia) == 0) then
             loc_ij(force_atm_point_i(ia),force_atm_point_j(ia)) = 0.0
          elseif (force_atm_uniform(ia) == -1) then
             loc_ij(:,:) = 0.0
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_I'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   force_flux_atm_I(ia,i,j) = loc_ij(i,j)
                end IF
             end DO
          end DO
          ! load forcing data array #II
          loc_ij(:,:) = const_real_zero
          if (force_atm_uniform(ia) == 2) then
             loc_ij(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_k)
          elseif (force_atm_uniform(ia) == 0) then
             loc_ij(force_atm_point_i(ia),force_atm_point_j(ia)) = 1.0
          elseif (force_atm_uniform(ia) == -1) then
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_SUR'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_II'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   force_flux_atm_II(ia,i,j) = loc_ij(i,j)
                end IF
             end DO
          end DO
          ! load forcing time series data
          loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_atm_'//TRIM(string_atm(ia))//'_sig'//TRIM(string_data_ext)
          loc_data_scale(:) = (/par_atm_force_scale_time(ia), par_atm_force_scale_val(ia)/)
          CALL sub_load_data_t2(loc_filename,loc_data_scale(:),force_flux_atm_sig(ia,:,:),loc_n_elements)
          ! set default forcing index values
          if (loc_n_elements == 0) THEN
             CALL sub_report_error( &
                  & 'biogem_data','init_force_flux_atm','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), &
                  & 'STOPPING', &
                  & (/const_real_null/),.TRUE. &
                  & )
          else
             force_flux_atm_sig_i(ia,:) = loc_n_elements
          end if
          ! warn if forcing information appears to be 'incomplete'
          if ( &
               & (minval(force_flux_atm_sig(ia,1,1:loc_n_elements)) > 0.0) &
               & .OR. &
               & (maxval(force_flux_atm_sig(ia,1,1:loc_n_elements)) < par_misc_t_runtime) &
               & ) then
             CALL sub_report_error( &
                  & 'biogem_data','sub_init_force_flux_atm', &
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// &
                  & '(1) alter the model start time year and/or run length (the goin file) or '// &
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', &
                  & 'CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
          end if
       end IF
    end DO
  END SUBROUTINE sub_init_force_flux_atm
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE FLUX FORCING - SEDIMENTS
  SUBROUTINE sub_init_force_flux_sed()
    ! local variables
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_n_elements
    integer::l,i,j,is
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(2)::loc_data_scale
    ! LOOP
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       IF (force_flux_sed_select(is)) THEN
          force_flux_sed_sig_i(is,:) = n_data_max
          force_flux_sed_sig(is,:,:) = 0.0
          ! load forcing data array #I
          loc_ij(:,:) = 0.0
          if (force_sed_uniform(is) == 2) then
             loc_ij(:,:) = 0.0
          elseif (force_sed_uniform(is) == 0) then
             loc_ij(force_sed_point_i(is),force_sed_point_j(is)) = 0.0
          elseif (force_sed_uniform(is) == -1) then
             loc_ij(:,:) = 0.0
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_I'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   force_flux_sed_I(is,i,j) = loc_ij(i,j)
                end IF
             end DO
          end DO
          ! load forcing data array #II
          loc_ij(:,:) = 0.0
          if (force_sed_uniform(is) == 2) then
             loc_ij(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_k)
          elseif (force_sed_uniform(is) == 0) then
             loc_ij(force_sed_point_i(is),force_sed_point_j(is)) = 1.0
          elseif (force_sed_uniform(is) == -1) then
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_SUR'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
          else
             loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_II'//TRIM(string_data_ext)
             CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
          end if
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   force_flux_sed_II(is,i,j) = loc_ij(i,j)
                end IF
             end DO
          end DO
          ! load forcing time series data
          loc_filename = TRIM(par_fordir_name)//'biogem_force_flux_sed_'//TRIM(string_sed(is))//'_sig'//TRIM(string_data_ext)
          loc_data_scale(:) = 1.0
          CALL sub_load_data_t2(loc_filename,loc_data_scale(:),force_flux_sed_sig(is,:,:),loc_n_elements)
          ! set default forcing index values
          if (loc_n_elements == 0) THEN
             CALL sub_report_error( &
                  & 'biogem_data','init_force_flux_sed','PLEASE PUT SOME DATA IN TIME SERIES FILE: '//TRIM(loc_filename), &
                  & 'STOPPING', &
                  & (/const_real_null/),.TRUE. &
                  & )
          else
             force_flux_sed_sig_i(is,:) = loc_n_elements
          end if
          ! warn if forcing information appears to be 'incomplete'
          if ( &
               & (minval(force_flux_sed_sig(is,1,1:loc_n_elements)) > 0.0) &
               & .OR. &
               & (maxval(force_flux_sed_sig(is,1,1:loc_n_elements)) < par_misc_t_runtime) &
               & ) then
             CALL sub_report_error( &
                  & 'biogem_data','sub_init_force_flux_sed', &
                  & 'The time interval of forcing data does not fully span the requested model integration, either; '// &
                  & '(1) alter the model start time year and/or run length (the goin file) or'// &
                  & '(2) alter the forcing signal stop time year (FILE: '//TRIM(loc_filename)//') ', &
                  & 'CONTINUING', &
                  & (/const_real_null/),.false. &
                  & )
          end if
       end IF
    end DO
  END SUBROUTINE sub_init_force_flux_sed
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE DATA SAVING
  SUBROUTINE sub_init_data_save_runtime()
    USE genie_util, ONLY:check_unit,check_iostat
    ! local variables
    INTEGER::l,io,ia,is,ic,ios
    integer::ib,id
    CHARACTER(len=255)::loc_filename
    CHARACTER(len=255)::loc_string
    real::loc_t = 0.0
    ! tracer
    IF (ctrl_data_save_sig_ocn) THEN
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','ocn_'//TRIM(string_ocn(io)),string_results_ext &
               & )
          IF (ctrl_data_save_sig_ocn_sur) THEN
             SELECT CASE (ocn_type(io))
             CASE (0)
                If (io == io_T) loc_string = '% time (yr) / temperature (C) / _surT (C) / _benT (C)'
                If (io == io_S) loc_string = '% time (yr) / salinity (o/oo) / _surS (o/oo) / _benS (o/oo)'
             CASE (1)
                loc_string = '% time (yr) / global ' //TRIM(string_ocn(io))//' (mol) / global ' // &
                     & TRIM(string_ocn(io))//' (mol kg-1)'// &
                     & ' / surface ' //TRIM(string_ocn(io))//' (mol kg-1) / benthic ' //TRIM(string_ocn(io))//' (mol kg-1)'
             CASE (11:20)
                loc_string = '% time (yr) / global '//TRIM(string_ocn(io))//' (mol) / global ' // &
                     & TRIM(string_ocn(io))//' (o/oo)'// &
                     & ' / surface ' //TRIM(string_ocn(io))//' (o/oo) / benthic ' //TRIM(string_ocn(io))//' (o/oo)'
             end SELECT
          else
             SELECT CASE (ocn_type(io))
             CASE (0)
                If (io == io_T) loc_string = '% time (yr) / temperature (K)'
                If (io == io_S) loc_string = '% time (yr) / salinity (o/oo)'
             CASE (1)
                loc_string = '% time (yr) / global '//TRIM(string_ocn(io))//' (mol) / global ' // &
                     & TRIM(string_ocn(io))//' (mol kg-1)'
             CASE (11:20)
                loc_string = '% time (yr) / global '//TRIM(string_ocn(io))//' (mol) / global ' // &
                     & TRIM(string_ocn(io))//' (o/oo)'
             end SELECT
          end IF
          SELECT CASE (ocn_type(io))
          CASE (0,1,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF
    ! atmospheric tracer
    IF (ctrl_data_save_sig_ocnatm) THEN
       DO l=1,n_l_atm
          ia = conv_iselected_ia(l)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','atm_'//TRIM(string_atm(ia)),string_results_ext &
               & )
          SELECT CASE (atm_type(ia))
          CASE (0)
             If (ia == ia_T) loc_string = '% time (yr) / surface air temperature (degrees C)'
             If (ia == ia_q) loc_string = '% time (yr) / surface humidity (???)'
          CASE (1)
             loc_string = '% time (yr) / global '//TRIM(string_atm(ia))//' (mol) / global ' //TRIM(string_atm(ia))//' (atm)'
          CASE (11:20)
             loc_string = '% time (yr) / global '//TRIM(string_atm(ia))//' (mol) / global ' //TRIM(string_atm(ia))//' (o/oo)'
          end SELECT
          SELECT CASE (atm_type(ia))
          CASE (0,1,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF
    ! export flux
    IF (ctrl_data_save_sig_fexport) THEN
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','fexport_'//TRIM(string_sed(is)),string_results_ext &
               & )
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged)
             loc_string = '% time (yr) / global '//TRIM(string_sed(is))//' flux (mol yr-1) / global '// &
                  & TRIM(string_sed(is))//' density (mol m-2 yr-1)'
          CASE (par_sed_type_age)
             loc_string = '% time (yr) / CaCO3 age (yr)'
          CASE (11:20)
             loc_string = '% time (yr) / global '//TRIM(string_sed(is))//' flux (mol yr-1) / global '// &
                  & TRIM(string_sed(is))//' delta (o/oo)'
          end SELECT
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged,par_sed_type_age,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF
    ! air-sea gas exchange
    IF (ctrl_data_save_sig_fairsea) THEN
       DO l=3,n_l_atm
          ia = conv_iselected_ia(l)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','fseaair_'//TRIM(string_atm(ia)),string_results_ext &
               & )
          SELECT CASE (atm_type(ia))
          CASE (1)
             loc_string = '% time (yr) / global '//TRIM(string_atm(ia))// ' sea->air transfer flux (mol yr-1) / '//&
                  & 'global '//TRIM(string_atm(ia))// ' density (mol m-2 yr-1)'
          CASE (11:20)
             loc_string = '% time (yr) / global '//TRIM(string_atm(ia))//' sea->air transfer flux (mol yr-1) / ' //&
                  & 'global '//TRIM(string_atm(ia))//' (o/oo)'
          end SELECT
          SELECT CASE (atm_type(ia))
          CASE (1,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    end IF
    ! ocean-atmosphere interface flux
    IF (ctrl_data_save_sig_focnatm) THEN
       DO l=3,n_l_atm
          ia = conv_iselected_ia(l)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','focnatm_'//TRIM(string_atm(ia)),string_results_ext &
               & )
          SELECT CASE (atm_type(ia))
          CASE (1)
             loc_string = '% time (yr) / global '//TRIM(string_atm(ia))//' flux (mol yr-1) / global '// &
                  & TRIM(string_atm(ia))//' density (mol m-2 yr-1) '//&
                  & ' NOTE: is the atmospheric forcing flux *net* of the sea-air gas exchange flux'// &
                  & ' (subtract off the sea->air flux to recover the atm forcing flux).'
          CASE (11:20)
             loc_string = '% time (yr) / global '//TRIM(string_atm(ia))//' flux (mol yr-1) / global '// &
                  & TRIM(string_atm(ia))//' (o/oo)'//&
                  & ' NOTE: is the atmospheric forcing flux *net* of the sea-air gas exchange flux'// &
                  & ' (subtract off the sea->air flux to recover the atm forcing flux).'
          end SELECT
          SELECT CASE (atm_type(ia))
          CASE (1,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF
    ! ocean->sediment flux
    IF (ctrl_data_save_sig_focnsed) THEN
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','focnsed_'//TRIM(string_sed(is)),string_results_ext &
               & )
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged)
             loc_string = '% time (yr) / global '//TRIM(string_sed(is))//' flux (mol yr-1) / global '// &
                  & TRIM(string_sed(is))//' density (mol m-2 yr-1)'
          CASE (par_sed_type_age)
             loc_string = '% time (yr) / CaCO3 age (yr)'
          CASE (11:20)
             loc_string = '% time (yr) / global '//TRIM(string_sed(is))//' flux (mol yr-1) / global '// &
                  & TRIM(string_sed(is))//' delta (o/oo)'
          end SELECT
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged,par_sed_type_age,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF
    ! sediment->ocean flux
    IF (ctrl_data_save_sig_fsedocn) THEN
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','fsedocn_'//TRIM(string_ocn(io)),string_results_ext &
               & )
          SELECT CASE (ocn_type(io))
          CASE (1)
             loc_string = '% time (yr) / global '//TRIM(string_ocn(io))//' flux (mol yr-1) / global '// &
                  & TRIM(string_ocn(io))//' density (mol m-2 yr-1)'
          CASE (11:20)
             loc_string = '% time (yr) / global '//TRIM(string_ocn(io))//' flux (mol yr-1) / global '// &
                  & TRIM(string_ocn(io))//' delta (o/oo)'
          end SELECT
          SELECT CASE (ocn_type(io))
          CASE (1,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF
    ! ocean surface carbonate chemistry
    IF (ctrl_data_save_sig_carb_sur) THEN
       DO ic=1,n_carb
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','carb_sur_'//TRIM(string_carb(ic)),string_results_ext &
               & )
          SELECT CASE (ic)
          CASE (ic_ohm_cal,ic_ohm_arg)
             loc_string = '% time (yr) / mean saturation state'
          CASE (ic_conc_CO2,ic_conc_HCO3,ic_conc_CO3)
             if (ocn_select(io_DIC_14C)) then
                loc_string = '% time (yr) / surface '//TRIM(string_carb(ic))//' (mol kg-1) / surface '//TRIM(string_carb(ic))// &
                     & ' d13C (o/oo) / surface '//TRIM(string_carb(ic))//' d14C (o/oo)'
             elseif (ocn_select(io_DIC_13C)) then
                loc_string = '% time (yr) / surface '//TRIM(string_carb(ic))//' (mol kg-1) / surface '//TRIM(string_carb(ic))// &
                     & ' d13C (o/oo)'
             else
                loc_string = '% time (yr) / surface '//TRIM(string_carb(ic))//' (mol kg-1)'
             end if
          CASE (ic_fug_CO2)
             loc_string = '% time (yr) / surface '//TRIM(string_carb(ic))//' (atm)'
          case default
             loc_string = '% time (yr) / surface '//TRIM(string_carb(ic))//' (mol kg-1)'
          end SELECT
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       END DO
    end if
    ! core-top sediment composition
    IF (ctrl_data_save_sig_ocnsed) THEN
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','sed_'//TRIM(string_sed(is)),string_results_ext &
               & )
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged)
             loc_string = '% time (yr) / mean '//TRIM(string_sed(is))//' composition (wt%)'
          CASE (par_sed_type_age)
             loc_string = '% time (yr) / mean core-top age (yr)'
          CASE (11:20)
             loc_string = '% time (yr) / mean core-top '//TRIM(string_sed(is))//' delta (o/oo)'
          end SELECT
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged,par_sed_type_age,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF
    ! miscellaneous
    IF (ctrl_data_save_sig_misc) THEN
       loc_filename=fun_data_timeseries_filename( &
            & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_gemlite',string_results_ext)
       loc_string = '% time (yr) / GEMlite weighting)'
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt=*,iostat=ios) trim(loc_string)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       loc_filename=fun_data_timeseries_filename( &
            & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_seaice',string_results_ext)
       loc_string = '% time (yr) / global sea-ice area (m2) / mean sea-ice cover (%) / '// &
            & 'global sea-ice volume (m3) / mean sea-ice thickness (m)'
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt=*,iostat=ios) trim(loc_string)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       loc_filename=fun_data_timeseries_filename( &
            & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_opsi',string_results_ext)
       loc_string = '% time (yr) / global min overturning (Sv) / global max overturning (Sv) / '// &
            & 'Atlantic min overturning (Sv) / Atlantic max overturning (Sv)'
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt=*,iostat=ios) trim(loc_string)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       IF (atm_select(ia_pCO2_14C)) THEN
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_atm_D14C',string_results_ext)
          loc_string = '% time (yr) / mean isotopic composition (o/oo)'
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end IF
       IF (ctrl_data_save_sig_carb_sur) THEN
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_surpH',string_results_ext)
          loc_string = '% time (yr) / mean ocean surface pH'
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end IF
       ! surface land (air) temperature SLT
       loc_filename=fun_data_timeseries_filename( &
            & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_SLT',string_results_ext)
       loc_string = '% time (yr) / mean (land) surface air temperature (degrees C)'
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt=*,iostat=ios) trim(loc_string)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       ! aeolian Fe diagnostics
       IF (ocn_select(io_Fe)) THEN
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_det_Fe_tot',string_results_ext)
          loc_string = '% time (yr) / total aeolian Fe input (mol yr-1)'
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_det_Fe_dis',string_results_ext)
          loc_string = '% time (yr) / dissolved aeolian Fe input (mol yr-1)'
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_det_Fe_sol',string_results_ext)
          loc_string = '% time (yr) / mean aeolian Fe solubility (%)'
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end IF
       ! insolation (wet grid only)
       loc_filename=fun_data_timeseries_filename( &
            & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_ocn_insol',string_results_ext)
       loc_string = '% time (yr) / mean (wet grid) insolation (W m-2) / ' // &
            & 'N hemisphere latitude (mean zonal) snap-shot (W m-2) / S hemisphere latitude (mean zonal) snap-shot (W m-2)'
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt=*,iostat=ios) trim(loc_string)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       loc_filename=fun_data_timeseries_filename( &
            & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_ocn_swflux',string_results_ext)
       loc_string = '% time (yr) / mean annual Sw flux at ocean surface (W m-2)'
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt=*,iostat=ios) trim(loc_string)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
    end if
    ! diagnostics
    IF (ctrl_data_save_sig_diag) THEN
       DO ib=1,n_diag_bio
          loc_filename=fun_data_timeseries_filename(loc_t, &
               & par_outdir_name,trim(par_outfile_name)//'_series_diag_bio',trim(string_diag_bio(ib)),string_results_ext)
          loc_string = '% time (yr) / global rate (mol yr-1) / mean rate (mol kg-1 yr-1)'
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end DO
       DO id=1,n_diag_geochem
          loc_filename=fun_data_timeseries_filename(loc_t, &
               & par_outdir_name,trim(par_outfile_name)//'_series_diag_geochem',trim(string_diag_geochem(id)),string_results_ext)
          loc_string = '% time (yr) / global rate (mol yr-1) / mean rate (mol kg-1 yr-1)'
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end DO
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          loc_filename=fun_data_timeseries_filename(loc_t, &
               & par_outdir_name,trim(par_outfile_name)//'_series_diag_weather',TRIM(string_ocn(io)),string_results_ext &
               & )
          SELECT CASE (ocn_type(io))
          CASE (1)
             loc_string = '% time (yr) / global '//TRIM(string_ocn(io))//' (mol yr-1)'
          CASE (11:20)
             loc_string = '% time (yr) / global '//TRIM(string_ocn(io))//' (o/oo)'
          end SELECT
          SELECT CASE (ocn_type(io))
          CASE (1,11:20)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt=*,iostat=ios) trim(loc_string)
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
       ! 2D misc diagnostics
       IF (force_restore_atm_select(ia_pCO2_13C) .AND. force_flux_atm_select(ia_pCO2_13C)) THEN
          loc_filename = fun_data_timeseries_filename(loc_t, &
               & par_outdir_name,trim(par_outfile_name)//'_series_diag_misc_2D','FpCO2',string_results_ext)
          loc_string = '% time (yr) / global rate (mol yr-1)'
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          write(unit=out,fmt=*,iostat=ios) trim(loc_string)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end IF
    END IF
  END SUBROUTINE sub_init_data_save_runtime
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! DATA SAVING ROUTINES - BioGeM
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE RUN-TIME DATA
  SUBROUTINE sub_data_save_runtime(dum_t)
    USE genie_util, ONLY:check_unit,check_iostat
    ! dummy arguments
    REAL,INTENT(in)::dum_t
    ! local variables
    INTEGER::l,io,ia,is,ic,ios
    integer::ib,id
    REAL::loc_t
    real::loc_opsi_scale
    real::loc_ocn_tot_M,loc_ocn_tot_M_sur,loc_ocn_tot_A
    real::loc_sig,loc_sig_sur,loc_sig_ben
    real::loc_tot,loc_tot_sur,loc_tot_ben
    real::loc_frac,loc_frac_sur,loc_frac_ben,loc_standard
    real::loc_d13C,loc_d14C
    CHARACTER(len=255)::loc_filename
    REAL,DIMENSION(n_carbisor)::loc_carbisor

    ! *** set-up local constants ***
    ! calculate local opsi conversion constant
    loc_opsi_scale = goldstein_dsc*goldstein_usc*const_rEarth*1.0E-6
    ! total ocean mass
    loc_ocn_tot_M = int_ocn_tot_M_sig/int_t_sig
    ! ocean surface mass
    loc_ocn_tot_M_sur = int_ocn_tot_M_sur_sig/int_t_sig
    ! ocean surface area
    loc_ocn_tot_A = sum(phys_ocn(ipo_A,:,:,n_k))
    ! local time
    loc_t = dum_t

    ! *** initialize local arrays
    loc_carbisor(:) = 0.0

    ! *** <sig_ocn_*> ***
    ! write ocean tracer data
    ! NOTE: write data both as the total inventory, and as the equivalent mean concentration
    IF (ctrl_data_save_sig_ocn) THEN
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','ocn_'//TRIM(string_ocn(io)),string_results_ext &
               & )
          SELECT CASE (ocn_type(io))
          CASE (0)
             If (io == io_T) then
                loc_sig = int_ocn_sig(io)/int_t_sig - const_zeroC
             else
                loc_sig = int_ocn_sig(io)/int_t_sig
             end If
             IF (ctrl_data_save_sig_ocn_sur) THEN
                If (io == io_T) then
                   loc_sig_sur = int_ocn_sur_sig(io)/int_t_sig - const_zeroC
                   loc_sig_ben = int_ocn_ben_sig(io)/int_t_sig - const_zeroC
                else
                   loc_sig_sur = int_ocn_sur_sig(io)/int_t_sig
                   loc_sig_ben = int_ocn_ben_sig(io)/int_t_sig
                end If
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
                WRITE(unit=out,fmt='(f12.3,3f9.3)',iostat=ios) &
                     & loc_t,                                  &
                     & loc_sig,                                &
                     & loc_sig_sur,                            &
                     & loc_sig_ben
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
             else
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
                WRITE(unit=out,fmt='(f12.3,f9.3)',iostat=ios) &
                     & loc_t,                                 &
                     & loc_sig
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
             end IF
          CASE (1)
             loc_sig = int_ocn_sig(io)/int_t_sig
             IF (ctrl_data_save_sig_ocn_sur) THEN
                loc_sig_sur = int_ocn_sur_sig(io)/int_t_sig
                loc_sig_ben = int_ocn_ben_sig(io)/int_t_sig
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
                WRITE(unit=out,fmt='(f12.3,4e15.7)',iostat=ios) &
                     & loc_t,                                   &
                     & loc_ocn_tot_M*loc_sig,                   &
                     & loc_sig,                                 &
                     & loc_sig_sur,                             &
                     & loc_sig_ben
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)                
             else
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
                WRITE(unit=out,fmt='(f12.3,2e15.7)',iostat=ios) &
                     & loc_t,                                   &
                     & loc_ocn_tot_M*loc_sig,                   &
                     & loc_sig
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
             end IF
          case (11:20)
             loc_tot      = int_ocn_sig(ocn_dep(io))/int_t_sig
             loc_frac     = int_ocn_sig(io)/int_t_sig
             loc_standard = const_standards(ocn_type(io))
             loc_sig      = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
             IF (ctrl_data_save_sig_ocn_sur) THEN
                loc_standard = const_standards(ocn_type(io))
                loc_tot_sur    = int_ocn_sur_sig(ocn_dep(io))/int_t_sig
                loc_frac_sur   = int_ocn_sur_sig(io)/int_t_sig
                loc_sig_sur    = fun_calc_isotope_delta(loc_tot_sur,loc_frac_sur,loc_standard,.FALSE.)
                loc_tot_ben    = int_ocn_ben_sig(ocn_dep(io))/int_t_sig
                loc_frac_ben   = int_ocn_ben_sig(io)/int_t_sig
                loc_sig_ben    = fun_calc_isotope_delta(loc_tot_ben,loc_frac_ben,loc_standard,.FALSE.)
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
                WRITE(unit=out,fmt='(f12.3,e15.7,3f12.3)',iostat=ios) &
                     & loc_t,                                         &
                     & loc_ocn_tot_M*loc_frac,                        &
                     & loc_sig,                                       &
                     & loc_sig_sur,                                   &
                     & loc_sig_ben
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
             else
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
                WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                     & loc_t,                                        &
                     & loc_ocn_tot_M*loc_frac,                       &
                     & loc_sig
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
             end IF
          END SELECT
       END DO
    END IF

    ! *** <sig_carb_sur_*> ***
    IF (ctrl_data_save_sig_carb_sur) THEN
       ! calculate isotopic properties of CO2(aq), HCO3-, and CO32-
       IF (ctrl_data_save_sig_ocn_sur) THEN
          if (ocn_select(io_DIC_13C)) then
             call sub_calc_carb_r13C( &
                  & int_ocn_sur_sig(io_T)/int_t_sig, &
                  & int_ocn_sur_sig(io_DIC)/int_t_sig, &
                  & int_ocn_sur_sig(io_DIC_13C)/int_t_sig, &
                  & int_carb_sur_sig(:)/int_t_sig, &
                  & loc_carbisor(:) &
                  & )
          end IF
          if (ocn_select(io_DIC_14C)) then
             call sub_calc_carb_r14C( &
                  & int_ocn_sur_sig(io_T)/int_t_sig, &
                  & int_ocn_sur_sig(io_DIC)/int_t_sig, &
                  & int_ocn_sur_sig(io_DIC_14C)/int_t_sig, &
                  & int_carb_sur_sig(:)/int_t_sig, &
                  & loc_carbisor(:) &
                  & )
          end IF
       end IF
       ! write ocean surface carbonate chemistry data
       ! NOTE: also write d13C and d14C isotopic properties of the carbonate species (CO2(aq), HCO3-, and CO32-)
       !       depending on whether either or both of these isotopic tracers have been selected
       DO ic=1,n_carb
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','carb_sur_'//TRIM(string_carb(ic)),string_results_ext &
               & )
          SELECT CASE (ic)
          CASE (ic_conc_CO2,ic_conc_HCO3,ic_conc_CO3)
             if (ocn_select(io_DIC_14C)) then
                loc_sig = int_carb_sur_sig(ic)/int_t_sig
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
                WRITE(unit=out,fmt='(f12.3,e15.7,2f12.3)',iostat=ios) &
                     & loc_t, &
                     & loc_sig, &
                     & fun_calc_isotope_delta(loc_sig,loc_carbisor(ic - 1)*loc_sig,const_standards(11),.FALSE.), &
                     & fun_calc_isotope_delta(loc_sig,loc_carbisor(ic + 3)*loc_sig,const_standards(12),.FALSE.)
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
             elseif (ocn_select(io_DIC_13C)) then
                loc_sig = int_carb_sur_sig(ic)/int_t_sig
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                     & loc_t, &
                     & loc_sig, &
                     & fun_calc_isotope_delta(loc_sig,loc_carbisor(ic - 1)*loc_sig,const_standards(11),.FALSE.)
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
             else
                loc_sig = int_carb_sur_sig(ic)/int_t_sig
                call check_unit(out,__LINE__,__FILE__)
                OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
                WRITE(unit=out,fmt='(f12.3,e15.7)',iostat=ios) &
                     & loc_t, &
                     & loc_sig
                call check_iostat(ios,__LINE__,__FILE__)
                CLOSE(unit=out,iostat=ios)
                call check_iostat(ios,__LINE__,__FILE__)
             end if
          case default
             loc_sig = int_carb_sur_sig(ic)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7)',iostat=ios) &
                  & loc_t, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    end if

    ! *** <sig_ocnatm_*> ***
    ! write atmosphere tracer data
    ! NOTE: write data both as the total inventory, and as the equivalent mean partial pressure
    ! NOTE: simple conversion factor from atm to mol is used
    IF (ctrl_data_save_sig_ocnatm) THEN
       DO l=1,n_l_atm
          ia = conv_iselected_ia(l)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','atm_'//TRIM(string_atm(ia)),string_results_ext &
               & )
          SELECT CASE (atm_type(ia))
          CASE (0)
             loc_sig = int_ocnatm_sig(ia)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,f9.3)',iostat=ios) &
                  & loc_t, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          CASE (1)
             loc_sig = int_ocnatm_sig(ia)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,2e15.7)',iostat=ios) &
                  & loc_t, &
                  & conv_atm_mol*loc_sig, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          case (11:20)
             loc_tot  = int_ocnatm_sig(atm_dep(ia))/int_t_sig
             loc_frac = int_ocnatm_sig(ia)/int_t_sig
             loc_standard = const_standards(atm_type(ia))
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                  & loc_t, &
                  & conv_atm_mol*loc_frac, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF

    ! *** <sig_fexport_*> ***
    ! write export flux data
    ! NOTE: write data both as mole and mass flux
    IF (ctrl_data_save_sig_fexport) THEN
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','fexport_'//TRIM(string_sed(is)),string_results_ext &
               & )
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged)
             loc_sig = int_fexport_sig(is)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,2e15.7)',iostat=ios) &
                  & loc_t, &
                  & loc_sig, &
                  & loc_sig/loc_ocn_tot_A
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          CASE (par_sed_type_age)
             if (int_fexport_sig(sed_dep(is)) > const_real_nullsmall) then
                loc_sig = int_fexport_sig(is)/int_t_sig
             else
                loc_sig = 0.0
             end if
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7)',iostat=ios) &
                  & loc_t, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          case (11:20)
             loc_tot  = int_fexport_sig(sed_dep(is))/int_t_sig
             loc_frac = int_fexport_sig(is)/int_t_sig
             loc_standard = const_standards(sed_type(is))
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                  & loc_t, &
                  & loc_frac, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF

    ! *** <int_diag_fseaair_sig_*> ***
    ! write air-sea has exchange flux data
    ! NOTE: write data both as the total flux, and as the equivalent mean flux density
    ! NOTE: a positive value of the array represents net ocean to atmosphere transfer
    IF (ctrl_data_save_sig_fairsea) THEN
       DO l=3,n_l_atm
          ia = conv_iselected_ia(l)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','fseaair_'//TRIM(string_atm(ia)),string_results_ext &
               & )
          SELECT CASE (atm_type(ia))
          CASE (1)
             loc_sig = int_diag_airsea_sig(ia)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                  & loc_t, &
                  & loc_sig, &
                  & loc_sig/loc_ocn_tot_A
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          case (11:20)
             loc_tot  = int_diag_airsea_sig(atm_dep(ia))/int_t_sig
             loc_frac = int_diag_airsea_sig(ia)/int_t_sig
             loc_standard = const_standards(atm_type(ia))
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                  & loc_t, &
                  & loc_frac, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF

    ! *** <sig_focnatm_*> ***
    ! write ocean-atmopshere interface flux data
    ! NOTE: write data both as the total flux, and as the equivalent mean flux density
    IF (ctrl_data_save_sig_focnatm) THEN
       DO l=3,n_l_atm
          ia = conv_iselected_ia(l)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','focnatm_'//TRIM(string_atm(ia)),string_results_ext &
               & )
          SELECT CASE (atm_type(ia))
          CASE (1)
             loc_sig = int_focnatm_sig(ia)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                  & loc_t, &
                  & loc_sig, &
                  & loc_sig/loc_ocn_tot_A
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          case (11:20)
             loc_tot  = int_focnatm_sig(atm_dep(ia))/int_t_sig
             loc_frac = int_focnatm_sig(ia)/int_t_sig
             loc_standard = const_standards(atm_type(ia))
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                  & loc_t, &
                  & loc_frac, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF

    ! *** <sig_focnsed_*> ***
    ! write ocean-sediment flux data
    ! NOTE: write data both as the total flux, and as the equivalent mean flux density
    ! NOTE: the surface ocean area is used as a proxy for the ocean bottom area
    IF (ctrl_data_save_sig_focnsed) THEN
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','focnsed_'//TRIM(string_sed(is)),string_results_ext &
               & )
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged)
             loc_sig = int_focnsed_sig(is)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,2e15.7)',iostat=ios) &
                  & loc_t, &
                  & loc_sig, &
                  & loc_sig/loc_ocn_tot_A
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          CASE (par_sed_type_age)
             if (int_focnsed_sig(sed_dep(is)) > const_real_nullsmall) then
                loc_sig = int_focnsed_sig(is)/int_t_sig
             else
                loc_sig = 0.0
             end if
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7)',iostat=ios) &
                  & loc_t, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          case (11:20)
             loc_tot  = int_focnsed_sig(sed_dep(is))/int_t_sig
             loc_frac = int_focnsed_sig(is)/int_t_sig
             loc_standard = const_standards(sed_type(is))
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                  & loc_t, &
                  & loc_frac, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF

    ! *** <sig_fsedocn_*> ***
    ! write sediment->ocean flux data
    ! NOTE: write data both as the total flux, and as the equivalent mean flux density
    ! NOTE: the surface ocean area is used as a proxy for the ocean bottom area
    IF (ctrl_data_save_sig_fsedocn) THEN
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','fsedocn_'//TRIM(string_ocn(io)),string_results_ext &
               & )
          SELECT CASE (ocn_type(io))
          CASE (1)
             loc_sig = int_fsedocn_sig(io)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,2e15.7)',iostat=ios) &
                  & loc_t, &
                  & loc_sig, &
                  & loc_sig/loc_ocn_tot_A
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          case (11:20)
             loc_tot  = int_fsedocn_sig(ocn_dep(io))/int_t_sig
             loc_frac = int_fsedocn_sig(io)/int_t_sig
             loc_standard = const_standards(ocn_type(io))
             loc_sig = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7,f12.3)',iostat=ios) &
                  & loc_t, &
                  & loc_frac, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF

    ! *** <sig_ocnsed_*> ***
    ! write sediment (core-top) composition data
    ! NOTE: the call to fun_sed_coretop made in populating <loc_sed_coretop> has already made the necessary type conversions
    !       for solid tracers as wt%, isotopes in per mill, and recovery of the carbonate 'age' value
    IF (ctrl_data_save_sig_ocnsed) THEN
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','sed_'//TRIM(string_sed(is)),string_results_ext &
               & )
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged)
             loc_sig = int_ocnsed_sig(is)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,f12.6)',iostat=ios) &
                  & loc_t, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          CASE (par_sed_type_age,11:20)
             loc_sig = int_ocnsed_sig(is)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,f12.3)',iostat=ios) &
                  & loc_t, &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end SELECT
       END DO
    END IF

    ! *** <sig_misc_*> ***
    ! write miscellaneous data (if requested)
    IF (ctrl_data_save_sig_misc) THEN
       loc_filename=fun_data_timeseries_filename( &
            & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_gemlite',string_results_ext)
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       WRITE(unit=out,fmt='(f12.3,e12.4)',iostat=ios) &
            & loc_t, &
            & (int_misc_gemlite_sig/int_t_sig)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       loc_filename=fun_data_timeseries_filename( &
            & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_seaice',string_results_ext)
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       WRITE(unit=out,fmt='(f12.3,e12.4,f9.3,e12.4,f9.3)',iostat=ios) &
            & loc_t, &
            & (int_misc_seaice_sig/int_t_sig), &
            & 100.0*(1.0/SUM(phys_ocn(ipo_A,:,:,n_k)))*int_misc_seaice_sig/int_t_sig, &
            & (int_misc_seaice_sig_vol/int_t_sig), &
            & (int_misc_seaice_sig_th/int_t_sig)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       loc_filename=fun_data_timeseries_filename( &
            & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_opsi',string_results_ext)
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       WRITE(unit=out,fmt='(f12.3,4f9.3)',iostat=ios)          &
            & loc_t,                                           &
            & loc_opsi_scale*int_misc_opsi_min_sig/int_t_sig,  &
            & loc_opsi_scale*int_misc_opsi_max_sig/int_t_sig,  &
            & loc_opsi_scale*int_misc_opsia_min_sig/int_t_sig, &
            & loc_opsi_scale*int_misc_opsia_max_sig/int_t_sig
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       ! atmospheric CO2 D14C
       IF (atm_select(ia_pCO2_14C)) THEN
          loc_tot  = int_ocnatm_sig(atm_dep(ia_pCO2))/int_t_sig
          loc_frac = int_ocnatm_sig(ia_pCO2_13C)/int_t_sig
          loc_standard = const_standards(atm_type(ia_pCO2_13C))
          loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          loc_frac = int_ocnatm_sig(ia_pCO2_14C)/int_t_sig
          loc_standard = const_standards(atm_type(ia_pCO2_14C))
          loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          loc_sig = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_atm_D14C',string_results_ext)
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          WRITE(unit=out,fmt='(f12.3,f12.3)',iostat=ios) loc_t,loc_sig
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end IF
       ! pH
       IF (ctrl_data_save_sig_carb_sur) THEN
          loc_filename=fun_data_timeseries_filename( &
               & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_surpH',string_results_ext)
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          WRITE(unit=out,fmt='(f12.3,f9.6)',iostat=ios) loc_t,-log10(int_carb_sur_sig(ic_H)/int_t_sig)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end IF
       ! SLT
       loc_filename = fun_data_timeseries_filename( &
            & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_SLT',string_results_ext)
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       WRITE(unit=out,fmt='(f12.3,f12.6)',iostat=ios) loc_t,int_misc_SLT_sig/int_t_sig
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       ! aeolian Fe diagnostics
       IF (ocn_select(io_Fe)) THEN
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_det_Fe_tot',string_results_ext)
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          WRITE(unit=out,fmt='(f12.3,e12.4)',iostat=ios) &
               & loc_t, &
               & (int_misc_det_Fe_tot_sig/int_t_sig)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_det_Fe_dis',string_results_ext)
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          WRITE(unit=out,fmt='(f12.3,e12.4)',iostat=ios) &
               & loc_t, &
               & (int_misc_det_Fe_dis_sig/int_t_sig)
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          if (int_misc_det_Fe_tot_sig > const_real_nullsmall) then
             loc_sig = 100.0*int_misc_det_Fe_dis_sig/int_misc_det_Fe_tot_sig
          else
             loc_sig = 0.0
          end if
          loc_filename=fun_data_timeseries_filename( &
               & loc_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_det_Fe_sol',string_results_ext)
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          WRITE(unit=out,fmt='(f12.3,f9.3)',iostat=ios) &
               & loc_t, &
               & loc_sig
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end IF
       ! insolation (wet grid only)
       loc_filename=fun_data_timeseries_filename( &
            & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_ocn_insol',string_results_ext)
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       WRITE(unit=out,fmt='(f12.3,3e12.4)',iostat=ios) &
            & loc_t, &
            & (int_misc_ocn_solfor_sig/int_t_sig), &
            & snap_misc_ocn_solfor_N_sig, &
            & snap_misc_ocn_solfor_S_sig
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       ! SW flux at surface (wet grid only)
       loc_filename=fun_data_timeseries_filename( &
            & dum_t,par_outdir_name,trim(par_outfile_name)//'_series','misc_ocn_swflux',string_results_ext)
       call check_unit(out,__LINE__,__FILE__)
       OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
       WRITE(unit=out,fmt='(f12.3,e12.4)',iostat=ios) &
            & loc_t, &
            & (int_misc_ocn_fxsw_sig/int_t_sig)
       call check_iostat(ios,__LINE__,__FILE__)
       CLOSE(unit=out,iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
    end if
    ! diagnostic diagnostics
    IF (ctrl_data_save_sig_diag) THEN
       DO ib=1,n_diag_bio
          loc_sig = int_diag_bio_sig(ib)/int_t_sig
          loc_filename=fun_data_timeseries_filename(loc_t, &
               & par_outdir_name,trim(par_outfile_name)//'_series_diag_bio',trim(string_diag_bio(ib)),string_results_ext)
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          WRITE(unit=out,fmt='(f16.3,2e15.6)',iostat=ios) &
               & loc_t,                        &
               & loc_ocn_tot_M_sur*loc_sig,    &
               & loc_sig
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end DO
       DO id=1,n_diag_geochem
          loc_sig = int_diag_geochem_sig(id)/int_t_sig
          loc_filename=fun_data_timeseries_filename(loc_t, &
               & par_outdir_name,trim(par_outfile_name)//'_series_diag_geochem',trim(string_diag_geochem(id)),string_results_ext)
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          WRITE(unit=out,fmt='(f16.3,2e15.6)',iostat=ios) &
               & loc_t,                        &
               & loc_ocn_tot_M*loc_sig,        &
               & loc_sig
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end DO
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          loc_filename=fun_data_timeseries_filename(dum_t, &
               & par_outdir_name,trim(par_outfile_name)//'_series_diag_weather',TRIM(string_ocn(io)),string_results_ext)
          SELECT CASE (ocn_type(io))
          CASE (1)
             loc_sig = int_diag_weather_sig(io)/int_t_sig
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,e15.7)',iostat=ios) &
                  & loc_t,                                  &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          case (11:20)
             loc_tot      = int_diag_weather_sig(ocn_dep(io))/int_t_sig
             loc_frac     = int_diag_weather_sig(io)/int_t_sig
             loc_standard = const_standards(ocn_type(io))
             loc_sig      = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             WRITE(unit=out,fmt='(f12.3,f12.3)',iostat=ios) &
                  & loc_t,                                  &
                  & loc_sig
             call check_iostat(ios,__LINE__,__FILE__)
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          END SELECT
       END DO
       ! 2D misc diagnostics
       IF (force_restore_atm_select(ia_pCO2_13C) .AND. force_flux_atm_select(ia_pCO2_13C)) THEN
          loc_filename = fun_data_timeseries_filename(dum_t, &
               & par_outdir_name,trim(par_outfile_name)//'_series_diag_misc_2D','FpCO2',string_results_ext)
          call check_unit(out,__LINE__,__FILE__)
          OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
          WRITE(unit=out,fmt='(f12.3,e15.7)',iostat=ios) loc_t,int_diag_misc_2D_sig(idiag_misc_2D_FpCO2)/int_t_sig
          call check_iostat(ios,__LINE__,__FILE__)
          CLOSE(unit=out,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       end IF
    end IF

  END SUBROUTINE sub_data_save_runtime
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE TIME-SLICE DATA
  SUBROUTINE sub_data_save_timeslice()
  END SUBROUTINE sub_data_save_timeslice
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE TIME-SLICE DATA
  SUBROUTINE sub_data_save_timeslice_sed()
  end SUBROUTINE sub_data_save_timeslice_sed
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE OCEAN-ATMOSPHERE FLUX DATA
  SUBROUTINE sub_data_save_flux_ocnatm()
  END SUBROUTINE sub_data_save_flux_ocnatm
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE DERIVED 'COLOR' TRACER DATA
  SUBROUTINE sub_data_save_ocn_col_extra()
  END SUBROUTINE sub_data_save_ocn_col_extra
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE GLOBAL DATA
  SUBROUTINE sub_data_save_global_snap(dum_t,dum_sfcatm1)
    USE genie_util, ONLY:check_unit,check_iostat
    ! dummy arguments
    REAL,INTENT(IN)::dum_t
    REAL,DIMENSION(n_atm,n_i,n_j),INTENT(in)::dum_sfcatm1      ! atmosphere composition interface array
    ! local variables
    INTEGER::l,ia,io,is,ios
    real::loc_tot,loc_frac,loc_standard
    real::loc_atm_ave,loc_ocn_ave,loc_sed_ave
    real::loc_ocn_tot_M,loc_ocn_tot_A
    CHARACTER(len=255)::loc_filename

    ! *** set local parameters ***
    loc_filename= &
         & fun_data_timesnap_filename( &
         & dum_t,par_outdir_name,trim(par_outfile_name)//'_year','diag_GLOBAL_SNAP',string_results_ext)
    ! total ocean mass
    loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
    ! ocean surface area
    loc_ocn_tot_A = sum(phys_ocn(ipo_A,:,:,n_k))

    ! *** save data - OPEN FILE ***
    call check_unit(out,__LINE__,__FILE__)
    OPEN(unit=out,file=TRIM(loc_filename),action='write',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)

    ! *** save data - ALL ***
    ! write atmospheric data
    DO l=3,n_l_atm
       ia = conv_iselected_ia(l)
       SELECT CASE (atm_type(ia))
       CASE (1)
          loc_atm_ave = &
               & SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia,:,:))/SUM(phys_ocnatm(ipoa_A,:,:))
          write(unit=out,fmt='(A13,A16,A3,f10.3,A15,A5,e13.7,A4)',iostat=ios) &
               & ' Atmospheric ',string_atm(ia),' : ', &
               & conv_mol_umol*loc_atm_ave, &
               & ' uatm          ', &
               & ' <-> ', &
               & conv_atm_mol*loc_atm_ave, &
               & ' mol'
          call check_iostat(ios,__LINE__,__FILE__)
       case (11:20)
          loc_tot = &
               & SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(atm_dep(ia),:,:))/SUM(phys_ocnatm(ipoa_A,:,:))
          loc_frac =  &
               & SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia,:,:))/SUM(phys_ocnatm(ipoa_A,:,:))
          loc_standard = const_standards(atm_type(ia))
          loc_atm_ave = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          write(unit=out,fmt='(A13,A16,A3,f10.3,A5)',iostat=ios) &
               & ' Atmospheric ',string_atm(ia),' : ', &
               & loc_atm_ave, &
               & ' o/oo'
          call check_iostat(ios,__LINE__,__FILE__)
       end SELECT
    END DO
    ! write ocean data
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       SELECT CASE (ocn_type(io))
       CASE (0)
          loc_ocn_ave = &
               & SUM(phys_ocn(ipo_M,:,:,:)*ocn(io,:,:,:))/loc_ocn_tot_M
          if (io == io_T) then
             write(unit=out,fmt='(A7,A16,A9,f10.3,A10)',iostat=ios) &
                  & ' Ocean ',string_ocn(io),'       : ',           &
                  & loc_ocn_ave - const_zeroC,                      &
                  & ' degrees C'
             call check_iostat(ios,__LINE__,__FILE__)
          else
             write(unit=out,fmt='(A7,A16,A9,f10.3,A4)',iostat=ios) &
                  & ' Ocean ',string_ocn(io),'       : ',          &
                  & loc_ocn_ave,                                   &
                  & ' PSU'
          end if
       CASE (1)
          loc_ocn_ave = &
               & SUM(phys_ocn(ipo_M,:,:,:)*ocn(io,:,:,:))/loc_ocn_tot_M
          write(unit=out,fmt='(A7,A16,A9,f10.3,A15,A5,e13.7,A4)',iostat=ios) &
               & ' Ocean ',string_ocn(io),'       : ', &
               & conv_mol_umol*loc_ocn_ave, &
               & ' umol kg-1     ', &
               & ' <-> ', &
               & loc_ocn_tot_M*loc_ocn_ave, &
               & ' mol'
          call check_iostat(ios,__LINE__,__FILE__)
       case (11:20)
          loc_tot = &
               & SUM(phys_ocn(ipo_M,:,:,:)*ocn(ocn_dep(io),:,:,:))/loc_ocn_tot_M
          loc_frac =  &
               & SUM(phys_ocn(ipo_M,:,:,:)*ocn(io,:,:,:))/loc_ocn_tot_M
          loc_standard = const_standards(ocn_type(io))
          loc_ocn_ave = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          write(unit=out,fmt='(A7,A16,A9,f10.3,A10)',iostat=ios) &
               & ' Ocean ',string_ocn(io),'       : ', &
               & loc_ocn_ave, &
               & ' o/oo     '
          call check_iostat(ios,__LINE__,__FILE__)
       end SELECT
    END DO
    ! write export data
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       SELECT CASE (sed_type(is))
       CASE (1,2,4)
          loc_sed_ave = SUM(int_bio_settle_timeslice(is,:,:,n_k))/int_t_timeslice/loc_ocn_tot_A
          write(unit=out,fmt='(A13,A16,A3,f10.3,A15,A5,e13.7,A9)',iostat=ios) &
               & ' Export flux ',string_sed(is),' : ', &
               & conv_mol_umol*loc_sed_ave/conv_m2_cm2, &
               & ' umol cm-2 yr-1', &
               & ' <-> ', &
               & loc_ocn_tot_A*loc_sed_ave, &
               & ' mol yr-1'
          call check_iostat(ios,__LINE__,__FILE__)
       case (11:20)
          loc_tot  = SUM(int_bio_settle_timeslice(sed_dep(is),:,:,n_k))/int_t_timeslice/loc_ocn_tot_A
          loc_frac = SUM(int_bio_settle_timeslice(is,:,:,n_k))/int_t_timeslice/loc_ocn_tot_A
          loc_standard = const_standards(sed_type(is))
          loc_sed_ave = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          write(unit=out,fmt='(A13,A16,A3,f10.3,A10)',iostat=ios) &
               & ' Export flux ',string_sed(is),' : ', &
               & loc_sed_ave, &
               & ' o/oo     '
          call check_iostat(ios,__LINE__,__FILE__)
       end SELECT
    END DO

    ! *** save data - CLOSE FILE ***
    call check_iostat(ios,__LINE__,__FILE__)
    CLOSE(unit=out)

  END SUBROUTINE sub_data_save_global_snap
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE GLOBAL DATA
  SUBROUTINE sub_data_save_global_av()
    USE genie_util, ONLY:check_unit,check_iostat
    ! local variables
    INTEGER::i,j,k,l,ia,io,is,ios,ic
    integer::loc_k1
    real::loc_t,loc_dt,loc_K
    real::loc_tot,loc_frac,loc_standard
    real::loc_atm_ave,loc_ocn_ave,loc_sed_ave
    real::loc_ocn_tot_M,loc_ocn_tot_A
    CHARACTER(len=255)::loc_filename
    REAL,DIMENSION(n_phys_ocn,n_i,n_j,n_k)::loc_phys_ocn       !
    REAL,DIMENSION(n_ocn,n_i,n_j,n_k)::loc_ocn                 !
    REAL,DIMENSION(n_carbconst,n_i,n_j,n_k)::loc_carbconst     !
    REAL,DIMENSION(n_carb,n_i,n_j,n_k)::loc_carb               !
    REAL,DIMENSION(n_carbalk,n_i,n_j,n_k)::loc_carbalk         ! 

    ! *** initialize local variables ***
    loc_phys_ocn(:,:,:,:)  = 0.0
    loc_ocn(:,:,:,:)       = 0.0
    loc_carbconst(:,:,:,:) = 0.0
    loc_carb(:,:,:,:)      = 0.0
    loc_carbalk(:,:,:,:)   = 0.0

    ! *** set local parameters ***
    loc_dt = int_t_timeslice
    loc_filename= &
         & fun_data_timeslice_filename( &
         & par_outdir_name,trim(par_outfile_name)//'_year','diag_GLOBAL_AVERAGE',string_results_ext)
    IF (ctrl_misc_t_BP) THEN
       loc_t = par_data_save_timeslice(par_data_save_timeslice_i) + par_misc_t_end
    ELSE
       loc_t = par_misc_t_end - par_data_save_timeslice(par_data_save_timeslice_i)
    END IF
    ! ocean physics
    If (ctrl_data_save_slice_phys_ocn) then
       loc_phys_ocn(:,:,:,:) = int_phys_ocn_timeslice(:,:,:,:)/int_t_timeslice
    else
       loc_phys_ocn(:,:,:,:) = phys_ocn(:,:,:,:)
    end If
    ! ocean tracers
    If (ctrl_data_save_slice_ocn) then
       loc_ocn(:,:,:,:) = int_ocn_timeslice(:,:,:,:)/int_t_timeslice
    else
       loc_ocn(:,:,:,:) = ocn(:,:,:,:)
    end if
    ! total ocean mass
    loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
    ! ocean surface area
    loc_ocn_tot_A = sum(phys_ocn(ipo_A,:,:,n_k))

    ! *** solve carbonate system ***
    IF (opt_select(iopt_select_carbchem)) THEN
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             DO k=goldstein_k1(i,j),n_k
                ! calculate carbonate dissociation constants
                CALL sub_calc_carbconst(           &
                     & loc_phys_ocn(ipo_Dmid,i,j,k), &
                     & loc_ocn(io_T,i,j,k),          &
                     & loc_ocn(io_S,i,j,k),          &
                     & loc_carbconst(:,i,j,k)        &
                     & )
                ! adjust carbonate constants
                if (ocn_select(io_Ca) .AND. ocn_select(io_Mg)) then
                   call sub_adj_carbconst(   &
                        & loc_ocn(io_Ca,i,j,k),  &
                        & loc_ocn(io_Mg,i,j,k),  &
                        & loc_carbconst(:,i,j,k) &
                        & )
                end if
                ! re-estimate Ca and borate concentrations from salinity (if not selected and therefore explicitly treated)
                IF (.NOT. ocn_select(io_Ca))  loc_ocn(io_Ca,i,j,n_k)  = fun_calc_Ca(loc_ocn(io_S,i,j,n_k))
                IF (.NOT. ocn_select(io_B))   loc_ocn(io_B,i,j,n_k)   = fun_calc_Btot(loc_ocn(io_S,i,j,n_k))
                IF (.NOT. ocn_select(io_SO4)) loc_ocn(io_SO4,i,j,n_k) = fun_calc_SO4tot(loc_ocn(io_S,i,j,n_k))
                IF (.NOT. ocn_select(io_F))   loc_ocn(io_F,i,j,n_k)   = fun_calc_Ftot(loc_ocn(io_S,i,j,n_k))
                ! seed default initial ocean pH
                loc_carb(ic_H,i,j,k) = 10**(-7.8)
                ! calculate carbonate chemistry
                CALL sub_calc_carb(        &
                     & loc_ocn(io_DIC,i,j,k),  &
                     & loc_ocn(io_ALK,i,j,k),  &
                     & loc_ocn(io_Ca,i,j,k),   &
                     & loc_ocn(io_PO4,i,j,k),  &
                     & loc_ocn(io_SiO2,i,j,k), &
                     & loc_ocn(io_B,i,j,k),    &
                     & loc_ocn(io_SO4,i,j,k),  &
                     & loc_ocn(io_F,i,j,k),    &
                     & loc_ocn(io_H2S,i,j,k),  &
                     & loc_ocn(io_NH4,i,j,k),  &
                     & loc_carbconst(:,i,j,k), & 
                     & loc_carb(:,i,j,k),      & 
                     & loc_carbalk(:,i,j,k)    & 
                     & )
             end do
          end DO
       end DO
    end IF

    ! *** save data - OPEN FILE ***
    call check_unit(out,__LINE__,__FILE__)
    OPEN(unit=out,file=TRIM(loc_filename),action='write',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)

    ! *** save data - OPEN FILE HEADER ***
    ! write header
    Write(unit=out,fmt=*) '=========================='
    Write(unit=out,fmt=*) 'GLOBAL DIAGNOSTICS'
    Write(unit=out,fmt=*) '=========================='
    Write(unit=out,fmt=*) ' '
    write(unit=out,fmt='(A23,f12.3)',iostat=ios) &
         & ' Year ............... : ',              &
         & loc_t
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A23,f12.3,A6)',iostat=ios) &
         & ' Integration interval : ',              &
         & int_t_timeslice,                        &
         & ' yr'
    call check_iostat(ios,__LINE__,__FILE__)

    ! *** save data - MISC / GLOBAL PHYSICAL PROPERTIES ***
    ! write misc data
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--------------------------'
    Write(unit=out,fmt=*) 'MISCELLANEOUS PROPERTIES'
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt='(A49,e13.7,A3)',iostat=ios) &
         & ' Global ocean k = n_k (surface) area ............. : ', &
         & SUM(loc_phys_ocn(ipo_A,:,:,n_k)), &
         & ' m2'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt='(A49,e13.7,A3)',iostat=ios) &
         & ' Global ocean k = (n_k - 1) (base of surface layer) area : ', &
         & SUM(loc_phys_ocn(ipo_A,:,:,n_k - 1)), &
         & ' m2'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt='(A49,e13.7,A3)',iostat=ios) &
         & ' Global ocean volume ......................... : ', &
         & SUM(loc_phys_ocn(ipo_V,:,:,:)), &
         & ' m3'
    call check_iostat(ios,__LINE__,__FILE__)
    loc_K = sum(int_phys_ocnatm_timeslice(ipoa_KCO2,:,:)*(1.0 - int_phys_ocnatm_timeslice(ipoa_seaice,:,:)))/ &
         & (sum(int_phys_ocn_timeslice(ipo_mask_ocn,:,:,n_k)*(1.0 - int_phys_ocnatm_timeslice(ipoa_seaice,:,:))))
    Write(unit=out,fmt='(A49,f8.6,A24)',iostat=ios) &
         & ' Global mean air-sea coefficient, K(CO2) ..... : ', &
         & loc_K, &
         & '     mol m-2 yr-1 uatm-1'
    call check_iostat(ios,__LINE__,__FILE__)

    ! *** save data - ATMOSPHERIC TRACER PROPERTIES ***
    ! write atmospheric data
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--------------------------'
    Write(unit=out,fmt=*) 'ATMOSPHERIC PROPERTIES'
    Write(unit=out,fmt=*) ' '
    DO l=3,n_l_atm
       ia = conv_iselected_ia(l)
       SELECT CASE (atm_type(ia))
       CASE (1)
          loc_atm_ave = &
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(ia,:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:))
          write(unit=out,fmt='(A13,A16,A3,f10.3,A5)',iostat=ios) &
               & ' Atmospheric ',string_atm(ia),' : ', &
               & conv_mol_umol*loc_atm_ave, &
               & ' uatm'
          call check_iostat(ios,__LINE__,__FILE__)
       case (11:20)
          loc_tot = &
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(atm_dep(ia),:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:))
          loc_frac =  &
               & SUM(phys_ocnatm(ipoa_A,:,:)*int_sfcatm1_timeslice(ia,:,:)/int_t_timeslice)/SUM(phys_ocnatm(ipoa_A,:,:))
          loc_standard = const_standards(atm_type(ia))
          loc_atm_ave = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          write(unit=out,fmt='(A13,A16,A3,f10.3,A5)',iostat=ios) &
               & ' Atmospheric ',string_atm(ia),' : ', &
               & loc_atm_ave, &
               & ' o/oo'
          call check_iostat(ios,__LINE__,__FILE__)
       end SELECT
    END DO

    ! *** save data - OCEAN TRACER PROPERTIES ***
    ! write ocean data
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--------------------------'
    Write(unit=out,fmt=*) 'BULK OCEAN PROPERTIES'
    Write(unit=out,fmt=*) ' '
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       SELECT CASE (ocn_type(io))
       CASE (0)
          loc_ocn_ave = &
               & SUM(loc_phys_ocn(ipo_M,:,:,:)*loc_ocn(io,:,:,:))/loc_ocn_tot_M
          if (io == io_T) then
             write(unit=out,fmt='(A7,A16,A9,f10.3,A10)',iostat=ios) &
                  & ' Ocean ',string_ocn(io),' ..... : ',           &
                  & loc_ocn_ave - const_zeroC,                      &
                  & ' degrees C'
             call check_iostat(ios,__LINE__,__FILE__)
          else
             write(unit=out,fmt='(A7,A16,A9,f10.3,A4)',iostat=ios) &
                  & ' Ocean ',string_ocn(io),' ..... : ',          &
                  & loc_ocn_ave,                                   &
                  & ' PSU'
          end if
       CASE (1)
          loc_ocn_ave = &
               & SUM(loc_phys_ocn(ipo_M,:,:,:)*loc_ocn(io,:,:,:))/loc_ocn_tot_M
          write(unit=out,fmt='(A7,A16,A9,f10.3,A10,A5,e13.7,A4)',iostat=ios) &
               & ' Ocean ',string_ocn(io),' ..... : ', &
               & conv_mol_umol*loc_ocn_ave, &
               & ' umol kg-1', &
               & ' <-> ', &
               & loc_ocn_tot_M*loc_ocn_ave, &
               & ' mol'
          call check_iostat(ios,__LINE__,__FILE__)
       case (11:20)
          loc_tot = &
               & SUM(loc_phys_ocn(ipo_M,:,:,:)*loc_ocn(ocn_dep(io),:,:,:))/loc_ocn_tot_M
          loc_frac =  &
               & SUM(loc_phys_ocn(ipo_M,:,:,:)*loc_ocn(io,:,:,:))/loc_ocn_tot_M
          loc_standard = const_standards(ocn_type(io))
          loc_ocn_ave = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          write(unit=out,fmt='(A7,A16,A9,f10.3,A10)',iostat=ios) &
               & ' Ocean ',string_ocn(io),' ..... : ', &
               & loc_ocn_ave, &
               & ' o/oo     '
          call check_iostat(ios,__LINE__,__FILE__)
       end SELECT
    END DO

    ! *** save data - CARBONATE CHEMSITRY ***
    ! write carbonate chemsitry data
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--------------------------'
    Write(unit=out,fmt=*) 'BULK OCEAN CARBONATE CHEMSITRY'
    Write(unit=out,fmt=*) ' '
    DO ic=1,n_carb
       SELECT CASE (ic)
       CASE (ic_conc_CO2,ic_conc_HCO3,ic_conc_CO3)
          loc_ocn_ave = &
               & SUM(loc_phys_ocn(ipo_M,:,:,:)*loc_carb(ic,:,:,:))/loc_ocn_tot_M
          write(unit=out,fmt='(A11,A16,A5,f10.3,A10,A5,e13.7,A4)',iostat=ios) &
               & ' Carb chem ',string_carb(ic),' . : ', &
               & conv_mol_umol*loc_ocn_ave, &
               & ' umol kg-1', &
               & ' <-> ', &
               & loc_ocn_tot_M*loc_ocn_ave, &
               & ' mol'
          call check_iostat(ios,__LINE__,__FILE__)
       end SELECT
    END DO

    ! *** save data - BIOLOGICAL EXPLORT ***
    ! write export data
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '------------------------------'
    Write(unit=out,fmt=*) 'SURFACE EXPORT PRODUCTION'
    Write(unit=out,fmt=*) ' '
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       SELECT CASE (sed_type(is))
       CASE (1,2,4)
          loc_sed_ave = SUM(int_bio_settle_timeslice(is,:,:,n_k))/int_t_timeslice/loc_ocn_tot_A
          write(unit=out,fmt='(A13,A16,A3,f10.3,A15,A5,e13.7,A9)',iostat=ios) &
               & ' Export flux ',string_sed(is),' : ', &
               & conv_mol_umol*loc_sed_ave/conv_m2_cm2, &
               & ' umol cm-2 yr-1', &
               & ' <-> ', &
               & loc_ocn_tot_A*loc_sed_ave, &
               & ' mol yr-1'
          call check_iostat(ios,__LINE__,__FILE__)
       case (11:20)
          loc_tot  = SUM(int_bio_settle_timeslice(sed_dep(is),:,:,n_k))/int_t_timeslice/loc_ocn_tot_A
          loc_frac = SUM(int_bio_settle_timeslice(is,:,:,n_k))/int_t_timeslice/loc_ocn_tot_A
          loc_standard = const_standards(sed_type(is))
          loc_sed_ave = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          write(unit=out,fmt='(A13,A16,A3,f10.3,A10)',iostat=ios) &
               & ' Export flux ',string_sed(is),' : ', &
               & loc_sed_ave, &
               & ' o/oo     '
          call check_iostat(ios,__LINE__,__FILE__)
       end SELECT
    END DO

    ! *** save data - SEDIMENTATION FLUX ***
    ! write sedimentation flux
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '------------------------------'
    Write(unit=out,fmt=*) 'SEDIMENTATION'
    Write(unit=out,fmt=*) ' '
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       SELECT CASE (sed_type(is))
       CASE (1,2,4)
          loc_sed_ave = SUM(int_focnsed_timeslice(is,:,:))/int_t_timeslice/loc_ocn_tot_A
          write(unit=out,fmt='(A13,A16,A3,f10.3,A15,A5,e13.7,A9)',iostat=ios) &
               & ' Export flux ',string_sed(is),' : ', &
               & conv_mol_umol*loc_sed_ave/conv_m2_cm2, &
               & ' umol cm-2 yr-1', &
               & ' <-> ', &
               & loc_ocn_tot_A*loc_sed_ave, &
               & ' mol yr-1'
          call check_iostat(ios,__LINE__,__FILE__)
       case (11:20)
          loc_tot  = SUM(int_focnsed_timeslice(sed_dep(is),:,:))/int_t_timeslice/loc_ocn_tot_A
          loc_frac = SUM(int_focnsed_timeslice(is,:,:))/int_t_timeslice/loc_ocn_tot_A
          loc_standard = const_standards(sed_type(is))
          loc_sed_ave = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
          write(unit=out,fmt='(A13,A16,A3,f10.3,A10)',iostat=ios) &
               & ' Export flux ',string_sed(is),' : ', &
               & loc_sed_ave, &
               & ' o/oo     '
          call check_iostat(ios,__LINE__,__FILE__)
       end SELECT
    END DO
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '------------------------------'
    Write(unit=out,fmt=*) 'CARBON FLUX SUMMARY'
    Write(unit=out,fmt=*) ' '
    write(unit=out,fmt='(A22,e15.7,A12,f7.3,A9)',iostat=ios) &
         & ' Total POC export   : ', &
         & SUM(int_bio_settle_timeslice(is_POC,:,:,n_k))/int_t_timeslice, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*SUM(int_bio_settle_timeslice(is_POC,:,:,n_k))/int_t_timeslice, &
         & ' GtC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A22,e15.7,A12,f7.3,A9)',iostat=ios) &
         & ' Total CaCO3 export : ', &
         & SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_k))/int_t_timeslice, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*SUM(int_bio_settle_timeslice(is_CaCO3,:,:,n_k))/int_t_timeslice, &
         & ' GtC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) ' '
    write(unit=out,fmt='(A22,e15.7,A12,f7.3,A9)',iostat=ios) &
         & ' Total POC rain     : ', &
         & SUM(int_focnsed_timeslice(is_POC,:,:))/int_t_timeslice, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*SUM(int_focnsed_timeslice(is_POC,:,:))/int_t_timeslice, &
         & ' GtC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A22,e15.7,A12,f7.3,A9)',iostat=ios) &
         & ' Total CaCO3 rain   : ', &
         & SUM(int_focnsed_timeslice(is_CaCO3,:,:))/int_t_timeslice, &
         & ' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*SUM(int_focnsed_timeslice(is_CaCO3,:,:))/int_t_timeslice, &
         & ' GtC yr-1'

    ! *** save data - CLOSE FILE ***
    call check_iostat(ios,__LINE__,__FILE__)
    CLOSE(unit=out)

  END SUBROUTINE sub_data_save_global_av
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! DATA SAVING ROUTINES - GOLDSTEIn
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE GRID DATA
  SUBROUTINE sub_data_save_topography()
    ! local variables
    CHARACTER(len=255)::loc_filename
    ! (i,j) topography (max height in m)
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_topography'//TRIM(string_data_ext)
    CALL sub_save_data_ij(loc_filename,n_i,n_j,-maxval(phys_ocn(ipo_mask_ocn,:,:,:)*phys_ocn(ipo_Dbot,:,:,:),3))
    ! grid point centre
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lat_mid'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,phys_ocn(ipo_lat,:,:,:))
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lon_mid'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,phys_ocn(ipo_lon,:,:,:))
    ! grid point limits
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lat_n'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,phys_ocn(ipo_latn,:,:,:))
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lat_s'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,phys_ocn(ipo_latn,:,:,:) - phys_ocn(ipo_dlat,:,:,:))
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lon_e'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,phys_ocn(ipo_lone,:,:,:))
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lon_w'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,phys_ocn(ipo_lone,:,:,:) - phys_ocn(ipo_dlon,:,:,:))
    ! layer height (m)
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lay_top'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,-phys_ocn(ipo_Dtop,:,:,:))
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lay_bot'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,-phys_ocn(ipo_Dbot,:,:,:))
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_lay_mid'//TRIM(string_data_ext)
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,-phys_ocn(ipo_Dmid,:,:,:))
  END SUBROUTINE sub_data_save_topography
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE STREAMFUNCTION DATA
  SUBROUTINE sub_data_save_goldstein_opsi()
    USE genie_util, ONLY:check_unit,check_iostat
    ! local variables
    INTEGER::j,k,ios
    REAL::loc_scale
    REAL,DIMENSION(n_k+1)::loc_grid_dz
    CHARACTER(len=255)::loc_filename
    ! initialize local variables
    loc_grid_dz(:) = 0.0
    loc_scale = goldstein_dsc*goldstein_usc*const_rEarth*1.0E-6
    !
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_opsi_lat'//TRIM(string_data_ext)
    call check_unit(out,__LINE__,__FILE__)
    OPEN(out,file=loc_filename,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    DO k=n_k,0,-1
       WRITE(unit=out,fmt='(999e14.6)',iostat=ios) ((180.0/const_pi) * ASIN(goldstein_sv(j)),j=0,n_j)
       call check_iostat(ios,__LINE__,__FILE__)
    ENDDO
    CLOSE(out,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! 
    loc_grid_dz(1:n_k) = goldstein_dz(:)
    loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)//'_grid_opsi_depth'//TRIM(string_data_ext)
    call check_unit(out,__LINE__,__FILE__)
    OPEN(out,file=loc_filename,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    DO k=n_k,0,-1
       WRITE(unit=out,fmt='(999e14.6)',iostat=ios) (SUM(-goldstein_dsc * loc_grid_dz(k+1:n_k+1)),j=0,n_j)
       call check_iostat(ios,__LINE__,__FILE__)
    ENDDO
    CLOSE(out,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! 
    loc_filename= &
         & fun_data_timeslice_filename( &
         & par_outdir_name,trim(par_outfile_name)//'_slice','misc_goldstein_opsi',string_results_ext)
    call check_unit(out,__LINE__,__FILE__)
    OPEN(out,file=loc_filename,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    DO k=n_k,0,-1
       WRITE(unit=out,fmt='(999e14.6)',iostat=ios) (loc_scale*int_opsi_timeslice(j,k)/int_t_timeslice,j=0,n_j)
       call check_iostat(ios,__LINE__,__FILE__)
    ENDDO
    CLOSE(out,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    loc_filename= &
         & fun_data_timeslice_filename( &
         & par_outdir_name,trim(par_outfile_name)//'_slice','misc_goldstein_opsia',string_results_ext)
    call check_unit(out,__LINE__,__FILE__)
    OPEN(out,file=loc_filename,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    DO k=n_k,0,-1
       WRITE(unit=out,fmt='(999e14.6)',iostat=ios) (loc_scale*int_opsia_timeslice(j,k)/int_t_timeslice,j=0,n_j)
       call check_iostat(ios,__LINE__,__FILE__)
    ENDDO
    CLOSE(out,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    loc_filename= &
         & fun_data_timeslice_filename( &
         & par_outdir_name,trim(par_outfile_name)//'_slice','misc_goldstein_opsip',string_results_ext)
    call check_unit(out,__LINE__,__FILE__)
    OPEN(out,file=loc_filename,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    DO k=n_k,0,-1
       WRITE(unit=out,fmt='(999e14.6)',iostat=ios) (loc_scale*int_opsip_timeslice(j,k)/int_t_timeslice,j=0,n_j)
       call check_iostat(ios,__LINE__,__FILE__)
    ENDDO
    CLOSE(out)

  END SUBROUTINE sub_data_save_goldstein_opsi
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE VELOCITY FIELD DATA
  SUBROUTINE sub_data_save_goldstein_u()
    ! local variables
    CHARACTER(len=255)::loc_filename
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk
    ! save data
    ! NOTE: scale to give velocity components in units of (m s-1);
    !       for the horizontal velocity components, the scale factor is usc (= 0.05) [Edwards and Shepherd, 2002]
    !       for the vertical velocity component, the overall scale factor is usc*dsc/rsc 
    !       (= 0.05*4000.0/6.36e6) [Edwards and Shepherd, 2002]
    loc_filename= fun_data_timeslice_filename( &
         & par_outdir_name,trim(par_outfile_name)//'_slice','misc_goldstein_u_1',string_results_ext)
    loc_ijk(:,:,:) = goldstein_usc*int_u_timeslice(1,:,:,:)/int_t_timeslice
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,loc_ijk(:,:,:))
    loc_filename= fun_data_timeslice_filename( &
         & par_outdir_name,trim(par_outfile_name)//'_slice','misc_goldstein_u_2',string_results_ext)
    loc_ijk(:,:,:) = goldstein_usc*int_u_timeslice(2,:,:,:)/int_t_timeslice
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,loc_ijk(:,:,:))
    loc_filename= fun_data_timeslice_filename( &
         & par_outdir_name,trim(par_outfile_name)//'_slice','misc_goldstein_u_3',string_results_ext)
    loc_ijk(:,:,:) = (goldstein_usc*goldstein_dsc/const_rEarth)*int_u_timeslice(3,:,:,:)/int_t_timeslice
    CALL sub_save_data_ijk(loc_filename,n_i,n_j,n_k,loc_ijk(:,:,:))
  END SUBROUTINE sub_data_save_goldstein_u
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! MISCELLANEOUS ROUTINES
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! RUN-TIME REPORTING
  SUBROUTINE sub_echo_runtime(dum_yr,dum_opsi_scale,dum_opsia_minmax,dum_sfcatm1,dum_gemlite)
    ! dummy arguments
    REAL,INTENT(in)::dum_yr
    REAL,INTENT(in)::dum_opsi_scale
    REAL,DIMENSION(2),INTENT(in)::dum_opsia_minmax
    REAL,DIMENSION(n_atm,n_i,n_j),INTENT(in)::dum_sfcatm1
    logical,intent(in)::dum_gemlite                                     ! in GEMlite phase of cycle?
    ! local variables
    real::loc_tot,loc_frac,loc_standard
    real::loc_ocn_tot_M,loc_ocn_tot_A,loc_ocnatm_tot_A
    real::loc_pCO2
    ! calculate local constants
    ! total ocean mass
    loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
    ! total ocean surface area (ice-free)
    loc_ocn_tot_A = sum((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k))
    ! total ocean-atmosphere interface area
    loc_ocnatm_tot_A = sum(phys_ocnatm(ipoa_A,:,:))
    ! calculate local isotopic variables
    loc_tot  = SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2,:,:))/loc_ocnatm_tot_A
    loc_frac = SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2_13C,:,:))/loc_ocnatm_tot_A
    loc_standard = const_standards(atm_type(ia_pCO2_13C))
    if (loc_frac < const_real_nullsmall) then
       loc_frac = fun_calc_isotope_fraction(0.0,loc_standard)*loc_tot
    end if

    ! *** echo run-time ***
    IF (opt_select(iopt_select_carbchem)) THEN
       ! print header (if necessary)
       if (par_misc_t_echo_header) then
          print*,' '
          ! ### MAKE MODIFICATIONS TO SCREEN PRINTING INFORMATION HERE ########################################################### !
          PRINT'(A5,A11,A3,A11,A9,A3,A9,A8,A8,A8,A3,A11,A11)', &
               & '    ',      &
               & ' model year', &
               & '  *',         &
               & ' pCO2(uatm)', &
               & '   d13CO2',   &
               & '  *',         &
               & '  AMO(Sv)',   &
               & '  ice(%)',    &
               & '   <SST>',    &
               & '   <SSS>',    &
               & '  *',         &
               & '  <DIC>(uM)', &
               & '  <ALK>(uM)'
          print*,' '
          par_misc_t_echo_header = .FALSE.
       end if
       ! calculate local variables
       loc_pCO2 = conv_mol_umol*SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2,:,:))/loc_ocnatm_tot_A
       ! print values
       if (dum_gemlite) then
          PRINT'(A5,F11.2,3X,F11.3,F9.3,3X,F9.3,F8.3,F8.3,F8.3,3X,F11.3,F11.3)', &
               & ' #G# ', &
               & dum_yr, &
               & loc_pCO2, &
               & fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.), &
               & dum_opsi_scale*dum_opsia_minmax(2), &
               & 100.0*(1.0/SUM(phys_ocn(ipo_A,:,:,n_k)))*SUM(phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_seaice,:,:)), &
               & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k)*ocn(io_T,:,:,n_k))/loc_ocn_tot_A - const_zeroC, &
               & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k)*ocn(io_S,:,:,n_k))/loc_ocn_tot_A, &
               & conv_mol_umol*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io_DIC,:,:,:))/loc_ocn_tot_M, &
               & conv_mol_umol*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io_ALK,:,:,:))/loc_ocn_tot_M
       else
          PRINT'(A5,F11.2,3X,F11.3,F9.3,3X,F9.3,F8.3,F8.3,F8.3,3X,F11.3,F11.3)', &
               & ' $N$ ', &
               & dum_yr, &
               & loc_pCO2, &
               & fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.), &
               & dum_opsi_scale*dum_opsia_minmax(2), &
               & 100.0*(1.0/SUM(phys_ocn(ipo_A,:,:,n_k)))*SUM(phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_seaice,:,:)), &
               & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k)*ocn(io_T,:,:,n_k))/loc_ocn_tot_A - const_zeroC, &
               & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k)*ocn(io_S,:,:,n_k))/loc_ocn_tot_A, &
               & conv_mol_umol*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io_DIC,:,:,:))/loc_ocn_tot_M, &
               & conv_mol_umol*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io_ALK,:,:,:))/loc_ocn_tot_M
       endif
       ! ######################################################################################################################### !
    else
       if (par_misc_t_echo_header) then
          print*,' '
          ! ### MAKE MODIFICATIONS TO SCREEN PRINTING INFORMATION HERE ########################################################### !
          PRINT'(A5,A11,A3,A9,A8,A8,A8)', &
               & '    ',      &
               & ' model year', &
               & '  *',         &
               & '  AMO(Sv)',   &
               & '  ice(%)',    &
               & '   <SST>',    &
               & '   <SSS>'
          print*,' '
          par_misc_t_echo_header = .FALSE.
       end if
       PRINT'(A5,F11.2,3X,F9.3,F8.3,F8.3,F8.3)', &
            & ' *** ', &
            & dum_yr, &
            & dum_opsi_scale*dum_opsia_minmax(2), &
            & 100.0*(1.0/SUM(phys_ocn(ipo_A,:,:,n_k)))*SUM(phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_seaice,:,:)), &
            & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k)*ocn(io_T,:,:,n_k))/loc_ocn_tot_A - const_zeroC, &
            & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k)*ocn(io_S,:,:,n_k))/loc_ocn_tot_A
       ! ######################################################################################################################### !
    end if

  END SUBROUTINE sub_echo_runtime
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! RUN-TIME max,min VALUE REPORTING
  SUBROUTINE sub_echo_maxmin()
    ! local variables
    integer::i,j,k
    integer::l,io
    integer::loc_i_min,loc_j_min,loc_k_min
    integer::loc_i_max,loc_j_max,loc_k_max
    real::loc_value_min
    real::loc_value_max
    real::loc_value,loc_tot,loc_frac,loc_standard

    ! *** determine max and min ocean tracer values + location ***
    IF (ctrl_audit) THEN
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          loc_value_min = const_real_nullhigh
          loc_value_max = const_real_null
          DO i=1,n_i
             DO j=1,n_j
                do k=goldstein_k1(i,j),n_k
                   SELECT CASE (ocn_type(io))
                   CASE (0,1)
                      loc_value = ocn(io,i,j,k)
                   case (11:20)
                      loc_tot = ocn(ocn_dep(io),i,j,k)
                      loc_frac = ocn(io,i,j,k)
                      loc_standard = const_standards(ocn_type(io))
                      loc_value = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
                   END SELECT
                   if (loc_value < loc_value_min) then
                      loc_value_min = loc_value
                      loc_i_min = i
                      loc_j_min = j
                      loc_k_min = k
                   end if
                   if (loc_value > loc_value_max) then
                      loc_value_max = loc_value
                      loc_i_max = i
                      loc_j_max = j
                      loc_k_max = k
                   end if
                end do
             end do
          end DO
          PRINT'(A5,A16,A3,A6,E10.4,A2,I2,A1,I2,A1,I2,A4,A6,E10.4,A2,I2,A1,I2,A1,I2,A1)', &
               & '     ', &
               & string_ocn(io), &
               & ' / ', &
               & 'min = ', &
               & loc_value_min, &
               & ' (', &
               & loc_i_min, &
               & ',', &
               & loc_j_min, &
               & ',', &
               & loc_k_min, &
               & ') / ', &
               & 'max = ', &
               & loc_value_max, &
               & ' (', &
               & loc_i_max, &
               & ',', &
               & loc_j_max, &
               & ',', &
               & loc_k_max, &
               & ')'
       end do
    end if

  END SUBROUTINE sub_echo_maxmin
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! OUTPUT AUDIT DIAGNOSTICS
  SUBROUTINE sub_data_audit_diagnostics()
    ! local variables
    INTEGER::l,io
    REAL::loc_ocn_rM
    ! calculate local constants
    loc_ocn_rM = 1.0/SUM(phys_ocn(ipo_M,:,:,:))
    DO l=3,n_l_ocn
       io = conv_iselected_io(l)
       SELECT CASE (io)
          ! \/\/\/ MAKE MODIFICATIONS TO REPORT ADDITIONAL TRACER AUDIT HERE
       CASE (io_DIC,io_NO3,io_PO4,io_Fe,io_O2,io_SiO2,io_ALK,io_Ca,io_B,io_SO4,io_F,io_colr,io_colb)
          PRINT*,'INITIAL / FINAL ',string_ocn(io),' inventory:', &
               & audit_ocn_init(io),audit_ocn_init(io)*loc_ocn_rM, &
               & '/', &
               & audit_ocn_new(io),audit_ocn_new(io)*loc_ocn_rM
       end SELECT
    END DO
  END SUBROUTINE sub_data_audit_diagnostics
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! REPLACE PARAMETER VALUES FOR CALIBRATION
  SUBROUTINE sub_init_bio_calibration()
    USE genie_util, ONLY:check_unit,check_iostat
    USE biogem_lib
    IMPLICIT NONE
    ! locals
    integer::ios
    ! *** LOAD PARAMETERS ***
    call check_unit(in,__LINE__,__FILE__)
    open(unit=in,file='goin_BIOGEM',status='old',action='read',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! skip over previously read in parameter section
    READ(unit=in,fmt='(1X)',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    READ(unit=in,fmt='(1X)',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    READ(unit=in,fmt='(1X)',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    READ(unit=in,fmt='(1X)',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    READ(unit=in,fmt='(1X)',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    READ(unit=in,fmt='(1X)',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    READ(unit=in,fmt='(1X)',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    READ(unit=in,fmt='(1X)',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! organic particulate carbon cycling
    read(unit=in,fmt=*,iostat=ios) par_bio_k0_PO4            ! maximum PO4 uptake rate
    call check_iostat(ios,__LINE__,__FILE__)
    read(unit=in,fmt=*,iostat=ios) par_bio_c0_PO4            ! PO4 half-sat constant
    call check_iostat(ios,__LINE__,__FILE__)
    read(unit=in,fmt=*,iostat=ios) par_bio_remin_POC_frac2   ! partitioning of POC into fraction #2
    call check_iostat(ios,__LINE__,__FILE__)
    read(unit=in,fmt=*,iostat=ios) par_bio_remin_POC_eL1     ! e-folding depth of POC fraction #1
    call check_iostat(ios,__LINE__,__FILE__)
!!$    read(unit=in,fmt=*) par_bio_remin_POC_eL2     ! e-folding depth of POC fraction #2
    ! inorganic particulate carbon cycling
    read(unit=in,fmt=*,iostat=ios) par_bio_red_POC_CaCO3     ! CaCO3:POC export 'rain ratio' scalar
    call check_iostat(ios,__LINE__,__FILE__)
    read(unit=in,fmt=*,iostat=ios) par_bio_red_POC_CaCO3_pP  ! calcification rate power
    call check_iostat(ios,__LINE__,__FILE__)
    read(unit=in,fmt=*,iostat=ios) par_bio_remin_CaCO3_frac2 ! partitioning of CaCO3 into fraction #2
    call check_iostat(ios,__LINE__,__FILE__)
    read(unit=in,fmt=*,iostat=ios) par_bio_remin_CaCO3_eL1   ! e-folding depth of CaCO3 fraction #1
    call check_iostat(ios,__LINE__,__FILE__)
!!$    read(unit=in,fmt=*) par_bio_remin_CaCO3_eL2   ! e-folding depth of CaCO3 fraction #2
    ! inorganic carbon cycling
!!$    read(unit=in,fmt=*) par_bio_red_DOMfrac       ! DOM fraction of export production
!!$    read(unit=in,fmt=*) par_bio_remin_DOMlifetime ! lifetime of DOM
    close(unit=in)
    ! DISPLAY
    print*,' '
    print*,'replace parameter values for calibration:'
    print*,'par_bio_k0_PO4            : ',par_bio_k0_PO4
    print*,'par_bio_c0_PO4            : ',par_bio_c0_PO4
    print*,'par_bio_remin_POC_frac2   : ',par_bio_remin_POC_frac2
    print*,'par_bio_remin_POC_eL1     : ',par_bio_remin_POC_eL1
!!$    print*,'par_bio_remin_POC_eL2     : ',par_bio_remin_POC_eL2
    print*,'par_bio_red_POC_CaCO3     : ',par_bio_red_POC_CaCO3
    print*,'par_bio_red_POC_CaCO3_pP  : ',par_bio_red_POC_CaCO3_pP
    print*,'par_bio_remin_CaCO3_frac2 : ',par_bio_remin_CaCO3_frac2
    print*,'par_bio_remin_CaCO3_eL1   : ',par_bio_remin_CaCO3_eL1
!!$    print*,'par_bio_remin_CaCO3_eL2   : ',par_bio_remin_CaCO3_eL2
!!$    print*,'par_bio_red_DOMfrac       : ',par_bio_red_DOMfrac
!!$    print*,'par_bio_remin_DOMlifetime : ',par_bio_remin_DOMlifetime
    print*,' '
  end SUBROUTINE sub_init_bio_calibration
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE SOLAR CONSTANT FORCING
  SUBROUTINE sub_init_force_solconst()
    ! local variables
    CHARACTER(len=255)::loc_filename
    INTEGER::loc_n_elements
    real,DIMENSION(2)::loc_data_scale
    ! load forcing time series data
    loc_filename = TRIM(par_fordir_name)//'biogem_force_solconst_sig'//TRIM(string_data_ext)
    loc_data_scale(:) = 1.0
    CALL sub_load_data_t2(loc_filename,loc_data_scale(:),force_solconst_sig(:,:),loc_n_elements)
    ! check that time elements are identical to those for time-series data saving
    if (size(force_solconst_sig(1,:)) /= size(par_data_save_sig(:))) THEN
       CALL sub_report_error( &
            & 'biogem_data','sub_init_force_solconst','PLEASE ENSURE THAT THE SAME TIME ELEMENTS ARE PRESENT IN: '&
            & //TRIM(loc_filename)//'AS IN THE TIME-SERIES SPECIFICATION FILE', &
            & 'STOPPING', &
            & (/const_real_null/),.TRUE. &
            & )
    end if
  END SUBROUTINE sub_init_force_solconst
  ! ****************************************************************************************************************************** !



END MODULE biogem_data
