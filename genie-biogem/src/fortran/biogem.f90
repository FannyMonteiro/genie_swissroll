
! ******************************************************************************************************************************** !
! BIOGEM LOOP SUBROUTINE
subroutine biogem(       &
     & dum_dts,          &
     & dum_genie_clock,  &
     & dum_sfcatm1,      &
     & dum_sfxatm1,      &
     & dum_sfcocn1,      &
     & dum_sfxocn1,      &
     & dum_sfcsed1,      &
     & dum_sfxsed1,      &
     & dum_sfxsumrok1    &
     & )
  use omp_lib
  use gem_carbchem
  USE biogem_lib
  USE biogem_box
  USE biogem_data
  IMPLICIT NONE
  ! DUMMY ARGUMENTS
  REAL,INTENT(IN)::dum_dts                                       ! biogem time-step length (seconds)
  integer(kind=8),INTENT(IN)::dum_genie_clock                    ! genie clock (ms since start) NOTE: 8-byte integer
  real,intent(in),dimension(n_atm,n_i,n_j)::dum_sfcatm1          ! atmosphere-interface tracer composition; ocn grid
  real,intent(inout),dimension(n_atm,n_i,n_j)::dum_sfxatm1       ! atmosphere -> surface flux; ocn grid
  real,intent(out),dimension(n_ocn,n_i,n_j)::dum_sfcocn1         ! sediment-interface ocean tracer composition; ocn grid
  real,intent(inout),dimension(n_ocn,n_i,n_j)::dum_sfxocn1       ! sediment -> ocean flux; ocn grid
  real,intent(in),dimension(n_sed,n_i,n_j)::dum_sfcsed1          ! sediment-interface sediment composition; ocn grid
  real,intent(inout),dimension(n_sed,n_i,n_j)::dum_sfxsed1       ! ocean -> sediment flux; ocn grid
  real,intent(inout),dimension(n_ocn,n_i,n_j)::dum_sfxsumrok1    ! coastal (weathering) -> ocean flux; ocn grid

  ! LOCAL VARIABLES
  INTEGER::i,j,k,l,io,ia,is,n                                    ! counting indices
  integer::loc_k1                                                ! local topography
  integer::loc_i,loc_j,loc_tot_i                                 ! 
!!$  integer::loc_m,loc_tot_m                                       ! tracer array conversion indices
  real::loc_dts,loc_dtyr,loc_t,loc_yr                            ! local time and time step BLAH actual year
  real::loc_rdts,loc_rdtyr                                       ! time reciprocals
  logical::loc_debug_ij                                          ! 
  logical,DIMENSION(n_ocn)::locio_mask                           ! 
  REAL,DIMENSION(n_ocn,n_i,n_j,n_k)::locijk_ocn                  ! local ocean tracer array
  REAL,DIMENSION(n_ocn,n_i,n_j,n_k)::locijk_focn                 ! local ocean tracer flux array
  REAL,DIMENSION(n_sed,n_i,n_j,n_k)::locijk_fpart                ! local particulate tracer flux array
  REAL,DIMENSION(n_atm,n_i,n_j)::locij_fatm                      ! local atmosphere tracer flux array
  REAL,DIMENSION(n_atm)::loc_datm_restore                        ! 
  REAL,DIMENSION(n_ocn)::loc_force_flux_dust                     ! 
  real::loc_ocn_tot_M,loc_ocn_rtot_M                             ! ocean mass and reciprocal
  real::loc_ocn_tot_V,loc_ocn_rtot_V                             ! ocean volume  and reciprocal
  real::loc_fS,loc_fT                                            !
  REAL,DIMENSION(n_ocn)::loc_force_restore_ocn_tmod              ! relaxation modifier
  REAL,DIMENSION(n_atm)::loc_force_restore_atm_tmod              ! relaxation modifier
  REAL,DIMENSION(n_atm,n_i,n_j)::locij_focnatm                   ! local ocn->atm flux (atm tracer currency), units (mol yr-1)
  REAL,DIMENSION(n_sed,n_i,n_j)::locij_focnsed                   ! local ocn->sed change (sed tracer currency), units (mol)
  REAL,DIMENSION(n_ocn,n_i,n_j)::locij_fsedocn                   ! local sed->ocean change (ocn tracer currency), units (mol)
  REAL,DIMENSION(n_ocn,n_i,n_j)::locij_frokocn                   ! local rok->ocean change (ocn tracer currency), units (mol)
  REAL,DIMENSION(n_ocn)::loc_ocnsed_audit                        ! temporary ocn->sed auditing array, units (mol)
  REAL,DIMENSION(n_ocn)::loc_fracdecay_ocn                       ! local reduction factor for decaying ocanic tracers
  REAL,DIMENSION(n_sed)::loc_fracdecay_sed                       ! local reduction factor for decaying sediment tracers
  real::loc_det_tot,loc_det_sol_tot,loc_det_Fe_sol_sf            !
  REAL,DIMENSION(n_ocn)::loc_fweather_tot,loc_fseddis_tot        ! local total weathering flux, dissolution flux
  real::loc_fsedpres_tot_CaCO3                                   ! local CaCO3 preservation
  REAL,DIMENSION(n_i,n_j)::locij_fsedpres_CaCO3                  ! local CaCO3 preservation
  real::loc_fsedpres_tot_opal                                    ! local opal preservation
  REAL,DIMENSION(n_i,n_j)::locij_fsedpres_opal                   ! local opal preservation
  real::loc_fsedpres_tot_POP                                     ! local Corg (POP) preservation
  REAL,DIMENSION(n_i,n_j)::locij_fsedpres_POP                    ! local Corg (POP) preservation
  real::loc_frac,loc_standard                                    !
  real::loc_delta_actual,loc_delta_target,loc_delta_source       !
  real::loc_sign  

  ! *** RESET GLOBAL ARRAYS ***
  ! reset remin array
  DO l=1,n_l_ocn
     io = conv_iselected_io(l)
     bio_remin(io,:,:,:) = 0.0
  end do
						  
  ! *** INITIALIZE LOCAL ARRAYS ***
  ! NOTE: only bother initializing for selected tracers
  DO l=3,n_l_atm
     ia = conv_iselected_ia(l)
     locij_fatm(:,:,:)    = 0.0
     locij_focnatm(:,:,:) = 0.0
  end do
  DO l=1,n_l_ocn
     io = conv_iselected_io(l)
     locijk_ocn(io,:,:,:)  = 0.0
     locijk_focn(io,:,:,:) = 0.0
     locij_fsedocn(io,:,:) = 0.0
     locij_frokocn(io,:,:) = 0.0
  end do
  DO l=1,n_l_sed
     is = conv_iselected_is(l)
     locij_focnsed(is,:,:)  = 0.0
     locijk_fpart(is,:,:,:) = 0.0
  end do

  ! *** CALCULATE GEM TIME ***
  ! update model time
  ! NOTE: par_misc_t_runtime is counted DOWN in years
  !       => for BIOGEM, the 'end of the world' occurs when time reaches zero
  loc_t = par_misc_t_runtime - real(dum_genie_clock)/(1000.0*conv_yr_s)
  loc_dts   = dum_dts
  loc_rdts  = 1.0/loc_dts
  loc_dtyr  = loc_dts/conv_yr_s
  loc_rdtyr = 1.0/loc_dtyr
  ! calculate actual year (counting years Before Present or otherwise)
  IF (ctrl_misc_t_BP) THEN
     loc_yr = loc_t + par_misc_t_end
  ELSE
     loc_yr = par_misc_t_end - loc_t
  END IF

  ! *** CALCULATE LOCAL CONSTANTS ***
  ! total ocean mass and recpirocal
  loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
  loc_ocn_rtot_M = 1.0/loc_ocn_tot_M
  ! total ocean volume and recpirocal
  loc_ocn_tot_V = sum(phys_ocn(ipo_V,:,:,:))
  loc_ocn_rtot_V = 1.0/loc_ocn_tot_V
  ! fractional reduction factors for decaying isotopes
  loc_fracdecay_ocn(:) = EXP(-loc_dtyr*const_lambda_ocn(:))
  loc_fracdecay_sed(:) = EXP(-loc_dtyr*const_lambda_sed(:))

  ! *** START-UP REPORTING ***
  if (loc_t > (par_misc_t_runtime - const_real_nullsmall)) then
     print*,' '
     print*,'>>> START BioGeM run-time diagnostics >>>'
     par_misc_t_echo_header = .TRUE.
  END IF

  ! ****************************** !
  ! *** UPDATE BIOGEOCHEMISTRY *** !
  ! ****************************** !
  ! main BioGeM loop - all updating of ocean-atmosphere-sediment biogeochemistry occurs within this,
  ! with ocean biogeochemistry is updated on a grid point by grid point ((i,j) loop) basis
  ! NOTE: update conditional on end of biogeochemical simulation not having been reached
  ! NOTE: tracer loops follow the GOLDSTEIn notation, with indices converted to their BioGeM equivalents
  if_go: IF (par_misc_t_go) THEN

     IF (ctrl_debug_lvl1) print*, '******************************'
     IF (ctrl_debug_lvl1) print*, '*** UPDATE BIOGEOCHEMISTRY ***'
     IF (ctrl_debug_lvl1) print*, '******************************'

     ! *** UPDATE RESTORING FORCING TIME-CONSTANTS ***
     IF (ctrl_debug_lvl1) print*, '*** UPDATE RESTORING FORCING TIME-CONSTANTS ***'
     ! recalculate time-varying restoring and flux forcings of the system & fractional relaxation
     ! NOTE: updating of restoring time-series position calculated seperately (sub: biogem_forcing)
     ! ATMOSPHERIC TRACERS (applied at the ocean-atmospere interface)
     DO l=3,n_l_atm
        ia = conv_iselected_ia(l)
        IF (force_restore_atm_select(ia)) THEN
           loc_force_restore_atm_tmod(ia) = 1.0 - EXP(-loc_dtyr/force_restore_atm_tconst(ia))
        END IF
     END DO
     ! OCEAN TRACERS
     DO l=1,n_l_ocn
        io = conv_iselected_io(l)
        IF (force_restore_ocn_select(io)) THEN
           loc_force_restore_ocn_tmod(io) = 1.0 - EXP(-loc_dtyr/force_restore_ocn_tconst(io))
        END IF
     END DO

     ! *** UPDATE DERIVED FORCING DATA ***
     loc_det_tot               = 0.0
     loc_det_sol_tot           = 0.0
     loc_fweather_tot(:)       = 0.0
     loc_fseddis_tot(:)        = 0.0
     loc_fsedpres_tot_CaCO3    = 0.0
     locij_fsedpres_CaCO3(:,:) = 0.0
     loc_fsedpres_tot_opal     = 0.0
     locij_fsedpres_opal(:,:)  = 0.0
     DO i=1,n_i
        DO j=1,n_j
           loc_k1 = goldstein_k1(i,j)
           IF (n_k >= loc_k1) THEN
              ! Aeolian solubilities
              if (force_flux_sed(is_det,i,j) > const_real_nullsmall) then
                 phys_ocnatm(ipoa_solFe,i,j) = force_flux_sed(is_det,i,j)**(par_det_Fe_sol_exp - 1.0)
                 loc_det_tot = loc_det_tot + force_flux_sed(is_det,i,j)
                 loc_det_sol_tot = loc_det_sol_tot + phys_ocnatm(ipoa_solFe,i,j)*force_flux_sed(is_det,i,j)
              else
                 phys_ocnatm(ipoa_solFe,i,j) = 0.0
              end if
              ! total global weathering (mol per time step)
              DO l=3,n_l_ocn
                 io = conv_iselected_io(l)
                 loc_fweather_tot(io) = loc_fweather_tot(io) + dum_sfxsumrok1(io,i,j)
              end DO
              ! sedimentary dissolution (mol per time step)
              DO l=3,n_l_ocn
                 io = conv_iselected_io(l)
                 loc_fseddis_tot(io) = loc_fseddis_tot(io) + loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io_Ca,i,j)
              end DO
              ! CaCO3 preservation (mol per time step)
              locij_fsedpres_CaCO3(i,j) = &
                   & bio_settle(is_CaCO3,i,j,loc_k1) - loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io_Ca,i,j)
              loc_fsedpres_tot_CaCO3 = loc_fsedpres_tot_CaCO3 + locij_fsedpres_CaCO3(i,j)
              ! opal preservation (mol per time step)
              locij_fsedpres_opal(i,j) = &
                   & bio_settle(is_opal,i,j,loc_k1) - loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io_SiO2,i,j)
              loc_fsedpres_tot_opal = loc_fsedpres_tot_opal + locij_fsedpres_opal(i,j)
              ! Corg (POP) preservation (mol per time step)
              locij_fsedpres_POP(i,j) = &
                   & bio_settle(is_POP,i,j,loc_k1) - loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io_PO4,i,j)
              loc_fsedpres_tot_POP = loc_fsedpres_tot_POP + locij_fsedpres_POP(i,j)
           end IF
        end DO
     end DO
     ! re-scale Aeolian solubilities
     ! => calculate the scale factor required to match the requested mean global solubility (<par_det_Fe_sol>)
     if ((loc_det_tot > const_real_nullsmall) .AND. (par_det_Fe_sol > const_real_nullsmall)) then
        loc_det_Fe_sol_sf = par_det_Fe_sol/(loc_det_sol_tot/loc_det_tot)
     else
        loc_det_Fe_sol_sf = 0.0
     end if
     phys_ocnatm(ipoa_solFe,:,:) = loc_det_Fe_sol_sf*phys_ocnatm(ipoa_solFe,:,:)


!!$     DO i=1,n_i
!!$        DO j=1,n_j
!!$           loc_k1 = goldstein_k1(i,j)
!!$           IF (n_k >= loc_k1) THEN
!!$              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
!!$                   & '*** WATER COLUMN REMINERALIZATION - DISSOLVED ORGANIC MATTER ***'
!!$              ! *** WATER COLUMN REMINERALIZATION - DISSOLVED ORGANIC MATTER ***
!!$              call sub_calc_bio_remin_allDOM(i,j,loc_k1,loc_dtyr)
!!$           end if
!!$        end do
!!$     end do

!!$  ! Fork a team of threads giving them their own copies of variables
!!$  ! *NB* No space between '!' and '$omp' below...
!!$  !$omp parallel private(nthreads, tid)
!!$  ! Obtain thread number
!!$  tid = OMP_get_thread_num()
!!$  !
!!$     do n=1,n_vocn
!!$        call sub_box_bio_remin_allDOM(vocn(n),vbio_remin(n),loc_dtyr)
!!$     end do
!!$  ! All threads join master thread and disband
!!$  !$omp end parallel


     ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
     ! *** (i,j) GRID PT LOOP START *** !
     ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

     block_iloop: DO i=1,n_i
        block_jloop: DO j=1,n_j

           ! *** INITIALIZE LOOP VARIABLES ***
           ! set local depth loop limit
           loc_k1 = goldstein_k1(i,j)
           ! set local debug condition
           if (opt_misc(iopt_misc_debugij)) loc_debug_ij = .TRUE.
           if ((i == par_misc_debug_i) .AND. (j == par_misc_debug_j)) then
              loc_debug_ij = .TRUE.
           else
              loc_debug_ij = .FALSE.
           end if

           IF (ctrl_debug_lvl1 .AND. loc_debug_ij) &
                & print*, '*** >>> START (i,j) GRID POINT'

           ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
           ! *** WET GRID PT CONDITIONALITY START *** !
           ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !

           IF (n_k >= loc_k1) THEN
           
              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** BRINE REJECTION ***'
              ! *** BRINE REJECTION ***
              ! calculate sea-ice brine rejection
              if (par_misc_brinerejection_frac > const_real_nullsmall) then
                 call sub_calc_misc_brinerejection(loc_dtyr,i,j,loc_fT,loc_fS)
                 locijk_focn(io_T,i,j,n_k)    = locijk_focn(io_T,i,j,n_k)    - loc_fT
                 locijk_focn(io_T,i,j,loc_k1) = locijk_focn(io_T,i,j,loc_k1) + loc_fT
                 locijk_focn(io_S,i,j,n_k)    = locijk_focn(io_S,i,j,n_k)    - loc_fS
                 locijk_focn(io_S,i,j,loc_k1) = locijk_focn(io_S,i,j,loc_k1) + loc_fS
              end if

              IF (opt_misc(iopt_misc_debugij)) print*,'(i,j); ',i,j

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** UPDATE AIR-SEA INTERFACE AQUEOUS SYSTEM ***'
              ! *** UPDATE AIR-SEA INTERFACE AQUEOUS SYSTEM ***
              ! re-calculate carbonate dissociation constants
              CALL sub_calc_carbconst(           &
                   & phys_ocn(ipo_Dmid,i,j,n_k), &
                   & ocn(io_T,i,j,n_k),          &
                   & ocn(io_S,i,j,n_k),          &
                   & carbconst(:,i,j,n_k)        &
                   & )
              ! adjust carbonate constants
              if (ocn_select(io_Ca) .AND. ocn_select(io_Mg)) then
                 call sub_adj_carbconst(     &
                      & ocn(io_Ca,i,j,n_k),  &
                      & ocn(io_Mg,i,j,n_k),  &
                      & carbconst(:,i,j,n_k) &
                      & )
              end if
              ! re-calculate gas solubility coefficients
              call sub_calc_solconst(i,j)
              ! re-calculate piston velocities
              call sub_calc_pv(i,j)
              ! (re-)solve surface carbonate equilibrium
              ! NOTE: only if DIC and ALK ocean tracers are both selected (hence the 'iopt_select_carbchem' flag is tested)
              IF (opt_select(iopt_select_carbchem)) THEN
                 ! re-estimate Ca and borate concentrations from salinity (if not selected and therefore explicitly treated)
                 IF (.NOT. ocn_select(io_Ca))  ocn(io_Ca,i,j,n_k)  = fun_calc_Ca(ocn(io_S,i,j,n_k))
                 IF (.NOT. ocn_select(io_B))   ocn(io_B,i,j,n_k)   = fun_calc_Btot(ocn(io_S,i,j,n_k))
                 IF (.NOT. ocn_select(io_SO4)) ocn(io_SO4,i,j,n_k) = fun_calc_SO4tot(ocn(io_S,i,j,n_k))
                 IF (.NOT. ocn_select(io_F))   ocn(io_F,i,j,n_k)   = fun_calc_Ftot(ocn(io_S,i,j,n_k))
                 ! re-calculate surface ocean carbonate chemistry
                 CALL sub_calc_carb(          &
                      & ocn(io_DIC,i,j,n_k),  &
                      & ocn(io_ALK,i,j,n_k),  &
                      & ocn(io_Ca,i,j,n_k),   &
                      & ocn(io_PO4,i,j,n_k),  &
                      & ocn(io_SiO2,i,j,n_k), &
                      & ocn(io_B,i,j,n_k),    &
                      & ocn(io_SO4,i,j,n_k),  &
                      & ocn(io_F,i,j,n_k),    &
                      & ocn(io_H2S,i,j,n_k),  &
                      & ocn(io_NH4,i,j,n_k),  &
                      & carbconst(:,i,j,n_k), & 
                      & carb(:,i,j,n_k),      &  
                      & carbalk(:,i,j,n_k)    & 
                      & )
                 ! estimate Revelle factor
                 CALL sub_calc_carb_RF0(      &
                      & ocn(io_DIC,i,j,n_k),  &
                      & ocn(io_ALK,i,j,n_k),  &
                      & ocn(io_PO4,i,j,n_k),  &
                      & ocn(io_SiO2,i,j,n_k), &
                      & ocn(io_B,i,j,n_k),    &
                      & ocn(io_SO4,i,j,n_k),  &
                      & ocn(io_F,i,j,n_k),    &
                      & ocn(io_H2S,i,j,n_k),  &
                      & ocn(io_NH4,i,j,n_k),  &
                      & carbconst(:,i,j,n_k), & 
                      & carb(:,i,j,n_k)    & 
                      & )
                 ! re-calculate carbonate system isotopic properties
                 if (ocn_select(io_DIC_13C)) then
                    call sub_calc_carb_r13C(        &
                         & ocn(io_T,i,j,n_k),       &
                         & ocn(io_DIC,i,j,n_k),     &
                         & ocn(io_DIC_13C,i,j,n_k), &
                         & carb(:,i,j,n_k),         &
                         & carbisor(:,i,j,n_k)      &
                         & )
                 end IF
                 if (ocn_select(io_DIC_14C)) then
                    call sub_calc_carb_r14C(        &
                         & ocn(io_T,i,j,n_k),       &
                         & ocn(io_DIC,i,j,n_k),     &
                         & ocn(io_DIC_14C,i,j,n_k), &
                         & carb(:,i,j,n_k),         &
                         & carbisor(:,i,j,n_k)      &
                         & )
                 end IF
              END IF

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** DECAY RADIOACTIVE TRACERS ***'
              ! *** DECAY RADIOACTIVE TRACERS ***
              ! NOTE: decay of isotopes in atmosphere is implemented in ATCHEM
              ! NOTE: the decay of ocean tracers is implemented as a (-ve.) flux forcing
              !       because of ensuring mass conservation in the salinity-normalization tracer advection scheme later ...
              !       and to do this, units of mol kg-1 (tracer concentration) must be converted to mol yr-1 (flux)
              ! NOTE: particulate tracers can be adjusted directly
              ! OCEAN TRACERS
              DO l=3,n_l_ocn
                 io = conv_iselected_io(l)
                 IF (abs(const_lambda_ocn(io)).gt.const_real_nullsmall) THEN
                    DO k=loc_k1,n_k
                       locijk_focn(io,i,j,k) = locijk_focn(io,i,j,k) - &
                            & phys_ocn(ipo_M,i,j,k)*(1.0 - loc_fracdecay_ocn(io))*ocn(io,i,j,k)/loc_dtyr
                    END DO
                 end IF
              end do
              ! SEDIMENT TRACERS
              DO l=1,n_l_sed
                 is = conv_iselected_is(l)
                 IF (abs(const_lambda_sed(is)).gt.const_real_nullsmall) THEN
                    DO k=loc_k1,n_k
                       bio_part(is,i,j,k) = loc_fracdecay_sed(is)*bio_part(is,i,j,k)
                    END DO
                 end if
              END DO

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) &
                   & print*,'*** PRODUCTION OF DECAY CHAIN ISOTOPES ***'
              ! PRODUCTION OF DECAY CHAIN ISOTOPES
              ! NOTE: Temporary version with implicit 230Th and 231Pa production;
              !       to be replaced by generic code that calculates production
              ! from the decay of the corresponding parent tracer
              call sub_calc_geochem_decay_chain(locijk_focn(:,i,j,:),i,j,loc_k1)

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** CALCULATE RESTORING BOUNDARY CONDITIONS ***'
              ! *** CALCULATE RESTORING BOUNDARY CONDITIONS ***
              ! calculate tracer changes necessary to meet any imposed atmospheric restoring boundary conditions
              ! ATMOSPHERIC TRACERS
              ! NOTE: divide restoring atmospheric flux by <ocn_dt> to produce a value in units of mol yr-1,
              !       (since in updating the atmosphere subsequently, the flux is multiplied by the timestep 
              !       to get the required increment)
              ! NOTE: use a crude conversion factor for partial pressure (atm) -> mol, BUT
              !       take into account that only 1/(imax x jmax) of the entire atmosphere is being forced under each grid point
              ! NOTE: do not force restore if time is outside the specified restoring interval
              !       (indicated by the time indices being identical)
              DO l=3,n_l_atm
                 ia = conv_iselected_ia(l)
                 IF (force_restore_atm_select(ia)) THEN
                    IF (force_restore_atm_sig_i(ia,1) /= force_restore_atm_sig_i(ia,2)) THEN
                       ! catch missing restoring data (non-isotope tracer values < 0.0) => force restoring flux to zero
                       SELECT CASE (atm_type(ia))
                       CASE (1)
                          if (force_restore_atm(ia,i,j) < 0.0) force_restore_atm(ia,i,j) = dum_sfcatm1(ia,i,j)
                       end select
                       ! set restoring flux
                       loc_datm_restore(ia) = (force_restore_atm(ia,i,j) - dum_sfcatm1(ia,i,j))*loc_force_restore_atm_tmod(ia)
                       locij_fatm(ia,i,j) = (1.0/real(n_i*n_j))*conv_atm_mol*loc_datm_restore(ia)*loc_rdtyr
                    END IF
                 end IF
              END DO
              ! OCEAN TRACERS #1
              ! hack to overwrite surface ocean restoring field and equilibrate with atmosphere if requested
              ! => this by-passes the need for explicit air-sea gas exchange
              ! NOTE: restoring time-constant (set in the tracer config file) will still apply
              ! NOTE: atmospheric inventory is not adjusted
              ! NOTE: assume one-to-one mapping between atm and ocn tracers
              !       => array need not be looped through to fine io corresponding to any chosen ia
              ! NOTE: if no correspondence, the NULL index of the ocean tracer array provides the dustbin
              DO l=3,n_l_atm
                 ia = conv_iselected_ia(l)
                 IF (atm_type(ia) == 1) THEN
                    if (ocnatm_airsea_eqm(ia)) then
                       loc_tot_i = conv_atm_ocn_i(0,ia)
                       do loc_i=1,loc_tot_i
                          io = conv_atm_ocn_i(loc_i,ia)
                          force_restore_ocn(io,i,j,n_k) = ocnatm_airsea_solconst(ia,i,j)*dum_sfcatm1(ia,i,j)
                          ! ### INSERT CODE TO DEAL WITH RELATED ISOTOPE COMPOSITION BOUNDARY CONDITIONS ######################### !
                          !
                          ! ###################################################################################################### !
                       end do
                    end if
                 end IF
              end DO
              ! OCEAN TRACERS #2
              ! calculate tracer changes necessary to meet any imposed oceanic restoring boundary conditions
              ! NOTE: fractional sea ice-covered area is taken into account and used to modify restoring flux at the surface
              ! NOTE: do not force restore if time is outside the specified restoring interval
              !       (indicated by the time indices being identical)
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 IF (force_restore_ocn_select(io)) THEN
                    IF (force_restore_ocn_sig_i(io,1) /= force_restore_ocn_sig_i(io,2)) THEN
                       DO k=force_restore_ocn_k1(io,i,j),n_k
                          ! catch missing restoring data (non-isotope tracer values < 0.0) => force restoring flux to zero
                          SELECT CASE (ocn_type(io))
                          CASE (0,1)
                             if (force_restore_ocn(io,i,j,k) < -const_real_nullsmall) force_restore_ocn(io,i,j,k) = ocn(io,i,j,k)
                          end select
                          ! set restoring flux
                          force_restore_docn_nuts(io) = (force_restore_ocn(io,i,j,k) - ocn(io,i,j,k))*loc_force_restore_ocn_tmod(io)
                          locijk_focn(io,i,j,k) = force_restore_docn_nuts(io)*phys_ocn(ipo_M,i,j,k)*loc_rdtyr
                       END DO
                       ! modify surface layer restoring forcing according to fractional sea-ice cover
                       locijk_focn(io,i,j,n_k) = (1.0 - phys_ocnatm(ipoa_seaice,i,j))*locijk_focn(io,i,j,n_k)
                    end IF
                 END IF
              END DO
              ! OCEAN TRACERS #3
              ! if a nutrient-restoring biological productivity 'biological' option has been set,
              ! => do not directly change nutrient concentrations via <locijk_focn>
              !    (this will instead be done in sub_calc_bio_uptake)
              SELECT CASE (par_bio_prodopt)
              CASE ('1N1T_PO4restore','1N1T_PO4restoreLL')
                 locijk_focn(io_PO4,i,j,n_k) = 0.0
              CASE DEFAULT
                 ! ### INSERT CODE TO INCLUDE ADDITIONAL SPECIAL CASES OF OCEAN TRACER RESTORING FORCING ######################### !
                 !
                 ! ############################################################################################################### !
              end SELECT

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** EXTERNAL FLUX FORCING ***'
              ! *** EXTERNAL FLUX FORCING ***
              ! calculate value of applied flux forcings
              ! NOTE: <force_flux*> in units of (mol yr-1)
              ! ATMOSPHERIC TRACERS (applied at the ocean-atmospere interface)
              DO l=3,n_l_atm
                 ia = conv_iselected_ia(l)
                 IF (force_flux_atm_select(ia)) THEN
                    locij_fatm(ia,i,j) = locij_fatm(ia,i,j) + force_flux_atm(ia,i,j)
                 END IF
              END DO
              ! OCEAN TRACERS
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 IF (force_flux_ocn_select(io)) THEN
                    DO k=loc_k1,n_k
                       locijk_focn(io,i,j,k) = locijk_focn(io,i,j,k) + force_flux_ocn(io,i,j,k)
                    END DO
                 END IF
              END DO
              ! SEDIMENT TRACERS #1
              ! NOTE: currently, fluxes are valid at the ocean surface only
              ! NOTE: addition is made directly to particulate sedimentary tracer array (scaled by time-step and cell mass)
              DO l=1,n_l_sed
                 is = conv_iselected_is(l)
                 IF (force_flux_sed_select(is)) THEN
                    locijk_fpart(is,i,j,n_k) = locijk_fpart(is,i,j,n_k) + force_flux_sed(is,i,j)
                 END IF
              END DO
              ! SEDIMENT TRACERS #2
              ! if a prescribed export flux based biological productivity 'biological' option has been set,
              ! => do not directly change nutrient concentrations via <locijk_focn>
              !    (this will instead be done in sub_calc_bio_uptake)
              ! => hack the value of the nutrient restoring array
              ! NOTE: force_restore_docn_nuts(io_PO4) in units of mol kg-1 (converted from the flux forcing units of mol kg-1)
              ! NOTE: use BIOGEM array conversion of POC -> POP and POP -> PO4 for completeness
              ! NOTE: negative sign of force_restore_docn_nuts for particulate creation ...
              SELECT CASE (par_bio_prodopt)
              CASE ('bio_POCflux')
                 force_restore_docn_nuts(io_PO4) = -phys_ocn(ipo_rM,i,j,n_k)*loc_dtyr* &
                      & conv_sed_ocn(io_PO4,is_POP)*bio_part_red(is_POC,is_POP,i,j)*locijk_fpart(is_POC,i,j,n_k)
                 locijk_fpart(is_POC,i,j,n_k) = 0.0
              CASE DEFAULT
                 ! ### INSERT CODE TO INCLUDE ADDITIONAL SPECIAL CASES OF SEDIMENT FLUX FORCING ################################## !
                 !
                 ! ############################################################################################################### !
              end SELECT

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** FLUX FORCING BELOW A SET THRESHOLD ***'
              ! NOTE: this code has a similar effect of restoring forcing, except that it allows 
              !       the isotopic properties of the flux to be prescribed rather than just the final isotopic state
              ! NOTE: the threshold is set via a prescribed restoring signal
              ! NOTE: this is all a bit of a hack ... :o)
              IF (force_restore_atm_select(ia_pCO2) .AND. force_flux_atm_select(ia_pCO2)) THEN
                 If (dum_sfcatm1(ia_pCO2,i,j) > force_restore_atm(ia_pCO2,i,j)) then
                    locij_fatm(ia_pCO2,i,j)     = 0.0
                    locij_fatm(ia_pCO2_13C,i,j) = 0.0
                 else
                    locij_fatm(ia_pCO2,i,j)     = force_flux_atm(ia_pCO2,i,j)
                    locij_fatm(ia_pCO2_13C,i,j) = force_flux_atm(ia_pCO2_13C,i,j)
                 end If
              end If

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** FLUX FORCING FOLLOWING A PRESCRIBED d13C HISTORY ***'
              ! NOTE: this code has a similar effect of restoring forcing, except that it allows
              !       the isotopic properties of a target signal to be followed via input/loss of carbon with a prescribed d13C
              ! NOTE: this is even more of a hack than before ;)
              IF (force_restore_atm_select(ia_pCO2_13C) .AND. force_flux_atm_select(ia_pCO2_13C)) THEN
                 IF (force_restore_atm_select(ia_pCO2) .AND. force_flux_atm_select(ia_pCO2)) THEN
                    ! calculate local variables
                    loc_frac = force_flux_atm(ia_pCO2_13C,i,j)/force_flux_atm(ia_pCO2,i,j)
                    loc_standard = const_standards(atm_type(ia_pCO2_13C))
                    loc_delta_actual = fun_calc_isotope_delta( &
                         & dum_sfcatm1(ia_pCO2,i,j),dum_sfcatm1(ia_pCO2_13C,i,j),loc_standard,.FALSE. &
                         & )
                    loc_delta_target = fun_calc_isotope_delta( &
                         & force_restore_atm(ia_pCO2,i,j),force_restore_atm(ia_pCO2_13C,i,j),loc_standard,.FALSE. &
                         & )
                    loc_delta_source = fun_calc_isotope_delta( &
                         & force_flux_atm(ia_pCO2,i,j),force_flux_atm(ia_pCO2_13C,i,j),loc_standard,.FALSE. &
                         & )
                    ! calculate the sign of the CO2 input
                    If (loc_delta_target > loc_delta_actual) then
                       if (loc_delta_source > loc_delta_actual) then
                          loc_sign = 1.0
                       else
                          loc_sign = -1.0
                       end If
                    else
                       if (loc_delta_source > loc_delta_actual) then
                          loc_sign = -1.0
                       else
                          loc_sign = 1.0
                       end If
                    end If
                    ! calculate flux of CO2 to atmosphere with specified d13C to approach atmospheric d13C target
                    ! NOTE: units of (mol yr-1)
                    locij_fatm(ia_pCO2,i,j) = loc_sign*force_flux_atm(ia_pCO2,i,j)
                    locij_fatm(ia_pCO2_13C,i,j) = loc_frac*locij_fatm(ia_pCO2,i,j)
                 else
                    locij_fatm(ia_pCO2,i,j)     = 0.0
                    locij_fatm(ia_pCO2_13C,i,j) = 0.0
                 end IF
              else
                 loc_delta_source = 0.0
              end If
              diag_misc_2D(idiag_misc_2D_FpCO2,i,j)     = locij_fatm(ia_pCO2,i,j)
              diag_misc_2D(idiag_misc_2D_FpCO2_13C,i,j) = loc_delta_source

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** OCEAN-ATMOPSHERE EXCHANGE FLUXES ***'
              ! *** OCEAN-ATMOPSHERE EXCHANGE FLUXES ***
              ! calculate ocean-atmosphere exchange
              ! NOTE: a positive value represents net ocean to atmosphere transfer
              ! NOTE: units of (mol yr-1)
              ! NOTE: <locij_fatm> is an array used to update the atmospheric tracer inventory
              ! NOTE: <locij_focnatm> is an array used to store the ocean -> atm gas flux for results reporting;
              !                       => it is not used in the updating of mass balance anywhere
              locij_focnatm(:,i,j) = fun_calc_ocnatm_flux(i,j,dum_sfcatm1(:,i,j),loc_dtyr)
              ! set local flux arrays for the updating of ocean and atmosphere reservoirs
              DO l=3,n_l_atm
                 ia = conv_iselected_ia(l)
                 locij_fatm(ia,i,j) = locij_fatm(ia,i,j) + locij_focnatm(ia,i,j)
                 loc_tot_i = conv_atm_ocn_i(0,ia)
                 do loc_i=1,loc_tot_i
                    io = conv_atm_ocn_i(loc_i,ia)
                    locijk_focn(io,i,j,n_k) = locijk_focn(io,i,j,n_k) - conv_atm_ocn(io,ia)*locij_focnatm(ia,i,j)
                 end do
              end DO
              diag_airsea(:,i,j) = locij_focnatm(:,i,j)
              ! KMM: Add sulfate and negative alkalinity to surface ocean to compensate for H2S loss to atmosphere
              ! NOTE: H2S can escape to the atmosphere, and this code implicitly accounts for the return sulphuric acid rain
              ! NOTE: also account for O2 consumed in the atmosphere in the oxidation of H2S
              ! NOTE: only do anythign if there is O2 in the atmosphere!!!
              if (ocn_select(io_H2S)) then
                 IF (                                                       &
                      & (locij_focnatm(ia_pH2S,i,j) > const_real_nullsmall) &
                      &  .AND.                                              &
                      & (dum_sfcatm1(ia_pO2,i,j) > const_real_nullsmall)    &
                      & ) THEN
                    ! mass balance adjustments
                    locijk_focn(io_ALK,i,j,n_k) = locijk_focn(io_ALK,i,j,n_k) - 2.0*locij_focnatm(ia_pH2S,i,j)
                    locijk_focn(io_SO4,i,j,n_k) = locijk_focn(io_SO4,i,j,n_k) + locij_focnatm(ia_pH2S,i,j)
                    locij_fatm(ia_pO2,i,j)  = locij_fatm(ia_pO2,i,j) - 2.0*locij_focnatm(ia_pH2S,i,j)
                    locij_fatm(ia_pH2S,i,j) = 0.0
                    ! update flux reporting
                    locij_focnatm(ia_pO2,i,j)  = locij_focnatm(ia_pO2,i,j) - 2.0*locij_focnatm(ia_pH2S,i,j)
                    locij_focnatm(ia_pH2S,i,j) = 0.0
                    ! ### INSERT CODE TO UPDATE <locij_fatmocn> ARRAY (NOT YET EVEN CREATED ...) ################################# !
                    !
                    ! ############################################################################################################ !
                 end IF
              END IF
              ! record air-sea gas exchange coefficient for posterity
              if ((1.0 - phys_ocnatm(ipoa_seaice,i,j)) > const_real_nullsmall) then
                 phys_ocnatm(ipoa_KCO2,i,j) = conv_umol_mol*phys_ocn(ipo_rA,i,j,n_k)* &
                      & (locij_focnatm(ia_pCO2,i,j)/(1.0 - phys_ocnatm(ipoa_seaice,i,j)))/ &
                      & ((carb(ic_conc_CO2,i,j,n_k)/ocnatm_airsea_solconst(ia_pCO2,i,j)) - dum_sfcatm1(ia_pCO2,i,j))
              end if

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** AEROSOL DISSOLUTION INPUT ***'
              ! *** AEROSOL DISSOLUTION INPUT ***
              ! dissolved from Fe from dust
              ! NOTE: par_det_frac_Fe is the mass fraction of Fe in dust
              ! NOTE: first convert detrital concentration (mol kg-1) back to mass,
              !       and then from mass of dust to mass of iron in the dust ...
              ! total aeolian Fe input (mol yr-1)
              phys_ocnatm(ipoa_totFe,i,j) = &
                   & conv_Fe_g_mol*par_det_Fe_frac*conv_det_mol_g*phys_ocn(ipo_M,i,j,n_k)*bio_part(is_det,i,j,n_k)/loc_dtyr
              ! dissolved Fe input (mol yr-1)
              loc_force_flux_dust(io_Fe) = phys_ocnatm(ipoa_solFe,i,j)*phys_ocnatm(ipoa_totFe,i,j)
              ! update ocean flux array
              locijk_focn(io_Fe,i,j,n_k) = locijk_focn(io_Fe,i,j,n_k) + loc_force_flux_dust(io_Fe)
              ! ### ADD CODE FOR ADDITIONAL AEOLIAN DISSOLUTION INPUTS ########################################################### !
              ! 
              ! ################################################################################################################## !

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** SEDIMENT DISSOLUTION INPUT ***'
              ! *** SEDIMENT DISSOLUTION INPUT ***
              ! modify remineralization array according to dissolution input from sediments
              ! NOTE: <dum_sfxocn1> in units of (mol m-2 s-1) - needs to be converted to (mol per time step)
              ! NOTE: if the model is configured as a 'closed' system (biogem_config.par),
              !       set dissolution flux equal to rain flux and add to overlying ocean
              ! NOTE: scale returned scavenged Fe according to value of par_scav_fremin
              ! NOTE: allow remineralization of Fe incorporated into organic matter (is_POFe)
              ! NOTE: in the event of a 'closed system', subtract any prescribed weathering from the sed->ocn return flux 
              !       + subtract solute in proportion to the estimated CaCO3 preservation flux
              !       (if not weathering flux is specified and ROKGEM not included, no modification is made)
              locij_fsedocn(:,i,j) = 0.0
              if (ctrl_force_sed_closedsystem) then
                 ! *** carbonate and Ca budget ***
                 If ((loc_fseddis_tot(io_Ca) > const_real_nullsmall) .AND. (loc_fweather_tot(io_Ca) > const_real_nullsmall)) then
                    ! adjust weathering input
                    ! NOTE: loop through sed tracers, and for CaCO3:
                    !       find CaCO3 related sed tracer; find corresponding ocn tracer
                    !       + re-scale ROKGEM weathering array appropriately
                    ! NOTE: specify dissolution flux according to sediment model
                    DO l=1,n_l_sed
                       is = conv_iselected_is(l)
                       if (sed_dep(is) == is_CaCO3) then
                          loc_tot_i = conv_sed_ocn_i(0,is)
                          do loc_i=1,loc_tot_i
                             io = conv_sed_ocn_i(loc_i,is)
                             dum_sfxsumrok1(io,i,j) = (loc_fsedpres_tot_CaCO3/loc_fweather_tot(io_Ca))*dum_sfxsumrok1(io,i,j)
                             locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                  & loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io,i,j)
                          end do
                       end if
                    end do
                 else
                    ! set dissolution flux equal to rain flux
                    DO l=1,n_l_sed
                       is = conv_iselected_is(l)
                       if (sed_dep(is) == is_CaCO3) then
                          loc_tot_i = conv_sed_ocn_i(0,is)
                          do loc_i=1,loc_tot_i
                             io = conv_sed_ocn_i(loc_i,is)
                             if (sed_type(is) == par_sed_type_scavenged) then
                                locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                     & par_scav_fremin*conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)
                             else
                                locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                     & conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)
                             end if
                          end do
                       end if
                    end DO
                 end If
                 ! *** opal and Si budget ***
                 If ((loc_fseddis_tot(io_SiO2) > const_real_nullsmall) .AND. (loc_fweather_tot(io_SiO2) > const_real_nullsmall)) then
                    ! adjust weathering input
                    ! NOTE: loop through sed tracers, and for opal:
                    !       find opal related sed tracer; find corresponding ocn tracer
                    !       + re-scale ROKGEM weathering array appropriately
                    ! NOTE: specify dissolution flux according to sediment model
                    DO l=1,n_l_sed
                       is = conv_iselected_is(l)
                       if (sed_dep(is) == is_opal) then
                          loc_tot_i = conv_sed_ocn_i(0,is)
                          do loc_i=1,loc_tot_i
                             io = conv_sed_ocn_i(loc_i,is)
                             dum_sfxsumrok1(io,i,j) = (loc_fsedpres_tot_opal/loc_fweather_tot(io_SiO2))*dum_sfxsumrok1(io,i,j)
                             locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                  & loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io,i,j)
                          end do
                       end if
                    end do
                 else
                    ! set dissolution flux equal to rain flux
                    DO l=1,n_l_sed
                       is = conv_iselected_is(l)
                       if (sed_dep(is) == is_opal) then
                          loc_tot_i = conv_sed_ocn_i(0,is)
                          do loc_i=1,loc_tot_i
                             io = conv_sed_ocn_i(loc_i,is)
                             if (sed_type(is) == par_sed_type_scavenged) then
                                locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                     & par_scav_fremin*conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)
                             else
                                locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                     & conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)
                             end if
                          end do
                       end if
                    end DO
                 end If
                 ! *** Corg and PO4 budget ***
                 ! NOTE: too complicated, even for just P, as there are e.g. Corg-associated ALK issues etc.
                 !       => omitt attempting to balance organic matter and/or its constituents
                 ! *** remaining solids ***
                 ! set dissolution flux equal to rain flux for remaining solids
                    DO l=1,n_l_sed
                       is = conv_iselected_is(l)
                       if ((sed_dep(is) /= is_CaCO3) .AND. (sed_dep(is) /= is_opal)) then
                          loc_tot_i = conv_sed_ocn_i(0,is)
                          do loc_i=1,loc_tot_i
                             io = conv_sed_ocn_i(loc_i,is)
                             if (sed_type(is) == par_sed_type_scavenged) then
                                locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                     & par_scav_fremin*conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)
                             else
                                locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                                     & conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)
                             end if
                          end do
                       end if
                    end DO
                 ! prevent return of dissolved Fe?
                 if (ctrl_bio_NO_fsedFe) locij_fsedocn(io_Fe,i,j) = 0.0
                 ! prevent return of dissolved 230Th, and 231Pa
                 locij_fsedocn(io_230Th,i,j) = 0.0
                 locij_fsedocn(io_231Pa,i,j) = 0.0
              else
                 ! set dissolution fluxes
                 DO l=3,n_l_ocn
                    io = conv_iselected_io(l)
                    locij_fsedocn(io,i,j) = locij_fsedocn(io,i,j) + &
                         & loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io,i,j)
                 end do
              end if
              DO l=3,n_l_ocn
                 io = conv_iselected_io(l)
                 bio_remin(io,i,j,loc_k1) = bio_remin(io,i,j,loc_k1) + &
                      & phys_ocn(ipo_rM,i,j,loc_k1)*locij_fsedocn(io,i,j)
              end do

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** TERRESTRIAL WEATHERING INPUT ***'
              ! *** TERRESTRIAL WEATHERING INPUT ***
              ! modify remineralization array according to terrestrial weathering input
              ! NOTE: <dum_sfxsumrok1> in units of (mol) (per time-step)
              ! NOTE: clear integrated weathering flux array when values are used
              ! NOTE: no screening for a 'closed system' is made as substractions are made appropriatly from
              !       'SEDIMENT DISSOLUTION INPUT' (below) if a weathering flux exists
              DO l=3,n_l_ocn
                 io = conv_iselected_io(l)
                 locij_frokocn(io,i,j) = locij_frokocn(io,i,j) + dum_sfxsumrok1(io,i,j)
!!!                 dum_sfxsumrok1(io,i,j) = 0.0
              end do
              DO l=3,n_l_ocn
                 io = conv_iselected_io(l)
                 bio_remin(io,i,j,n_k) = bio_remin(io,i,j,n_k) + phys_ocn(ipo_rM,i,j,n_k)*locij_frokocn(io,i,j)
              end do
!!!              diag_weather(:,i,j) = locij_frokocn(:,i,j)

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN GEOCHEMISTRY - H2S OXIDATION ***'
              ! *** WATER COLUMN GEOCHEMISTRY - H2S OXIDATION ***
              if (ocn_select(io_O2) .AND. ocn_select(io_SO4) .AND. ocn_select(io_H2S)) then
                 call sub_calc_geochem_oxidize_H2S(i,j,loc_k1,loc_dtyr)
              end If

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - CH4 OXIDATION ***'
              ! *** WATER COLUMN REMINERALIZATION - CH4 OXIDATION ***
              if (ocn_select(io_O2) .AND. ocn_select(io_CH4)) then
                 call sub_calc_bio_remin_oxidize_CH4(i,j,loc_k1,loc_dtyr)
              end If

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - DISSOLVED ORGANIC MATTER ***'
              ! *** WATER COLUMN REMINERALIZATION - DISSOLVED ORGANIC MATTER ***
              call sub_calc_bio_remin_allDOM(i,j,loc_k1,loc_dtyr)

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - PARTICULATE MATTER ***'
              ! *** WATER COLUMN REMINERALIZATION - PARTICULATE MATTER ***
              ! NOTE: remineralize particulate material left over from PREVIOUS time step,
              !       i.e., do not immediately remineralize material newly created in the current time step
              !       (which is why the subroutine <calc_bio_remin> is called before <calc_bio_uptake>)
              ! NOTE: units of (mol kg-1)
              call sub_calc_bio_remin(i,j,loc_k1,loc_dtyr)

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - CONVERT ON -> NH4 ***'
              ! *** WATER COLUMN REMINERALIZATION - CONVERT ON -> NH4 ***
              if (ocn_select(io_O2) .AND. ocn_select(io_NO3) .AND. ocn_select(io_NH4)) then
                 if(ctrl_bio_remin_ONtoNH4) call sub_calc_bio_remin_NO3toNH4(i,j,loc_k1)
              end If

	      ! Moved after conversion of ON into NH4 to fix bug - Fanny (April 2015)
              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - NH4 OXIDATION ***'
              ! *** WATER COLUMN REMINERALIZATION - NH4 OXIDATION ***
              if (ocn_select(io_O2) .AND. ocn_select(io_NO3) .AND. ocn_select(io_NH4)) then
                 call sub_calc_bio_remin_oxidize_NH4(i,j,loc_k1,loc_dtyr)
              end If

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - NITRATE REDUCTION ***'
              ! *** WATER COLUMN REMINERALIZATION - NITRATE REDUCTION ***
              if (ocn_select(io_O2) .AND. ocn_select(io_NO3) .AND. ocn_select(io_N2) .AND. ocn_select(io_NH4)) then
                 call sub_calc_bio_remin_reduce_NO3(i,j,loc_k1)
              end If

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** WATER COLUMN REMINERALIZATION - SULPHATE REDUCTION ***'
              ! *** WATER COLUMN REMINERALIZATION - SULPHATE REDUCTION ***
              if (ocn_select(io_O2) .AND. ocn_select(io_SO4) .AND. ocn_select(io_H2S)) then
                 if (ocn_select(io_NO3) .AND. ocn_select(io_NH4)) then
                    call sub_calc_bio_remin_reduce_SO4(i,j,loc_k1)
                 else
                    call sub_calc_bio_remin_reduce_SO4_SIMPLE(i,j,loc_k1)                    
                 end if
              end If

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** SURFACE OCEAN BIOLOGICAL PRODUCTIVITY ***'
              ! *** SURFACE OCEAN BIOLOGICAL PRODUCTIVITY ***
              ! NOTE: particulate production and update in units of (mol kg-1)
              call sub_calc_bio(i,j,loc_k1,loc_dtyr)

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) print*, &
                   & '*** OCEAN ABIOTIC PRECIPITATION ***'
              ! *** OCEAN ABIOTIC PRECIPITATION ***
              ! NOTE: units of (mol kg-1)
              if (ctrl_bio_CaCO3precip .AND. sed_select(is_CaCO3)) then
                 call sub_calc_bio_uptake_abio(i,j,loc_k1,loc_dtyr)
              end if

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) &
                   & print*,'*** WATER COLUMN GEOCHEMISTRY - Fe SPECIATION ***'
              ! *** WATER COLUMN GEOCHEMISTRY - Fe SPECIATION ***
              ! NOTE: although <locijk_focn> Fe is added to the remin array within the sub_calc_geochem_Fe subroutine,
              !       the same flux is subtracted again after equilibrium has been calculated,
              !       hence <locijk_focn> Fe later can be added 'as normal' in updating the <ocn> array
              if (sed_select(is_det) .AND. ocn_select(io_Fe)) then
                 call sub_calc_geochem_Fe(i,j,loc_k1,loc_dtyr*phys_ocn(ipo_rM,i,j,:)*locijk_focn(io_Fe,i,j,:))
              end If

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) &
                   & print*, '*** INTERFACE ARRAY UPDATE ***'
              ! *** INTERFACE ARRAY UPDATE ***
              ! (1) set ocn->atm flux
              ! NOTE: convert units from (mol yr-1) to (mol m-2 s-1)
              ! NOTE: update external interface array with TOTAL flux to atmosphere
              !       (i.e., due to both air-sea exchange and any forcing of the atmosphere)
              DO l=3,n_l_atm
                 ia = conv_iselected_ia(l)
                 dum_sfxatm1(ia,i,j) = phys_ocnatm(ipoa_rA,i,j)*conv_s_yr*locij_fatm(ia,i,j)
              end do
              ! (2) set bottom-water tracers
              ! NOTE: estimate Ca and borate concentrations from salinity (if not selected and therefore not explicitly treated)
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 dum_sfcocn1(io,i,j) = ocn(io,i,j,loc_k1) + &
                      & bio_remin(io,i,j,loc_k1) + loc_dtyr*phys_ocn(ipo_rM,i,j,loc_k1)*locijk_focn(io,i,j,loc_k1)
              end do
              IF (.NOT. ocn_select(io_Ca))    dum_sfcocn1(io_Ca,i,j) = fun_calc_Ca(dum_sfcocn1(io_S,i,j))
              IF (.NOT. ocn_select(io_B))     dum_sfcocn1(io_B,i,j)  = fun_calc_Btot(dum_sfcocn1(io_S,i,j))
              IF (.NOT. ocn_select(io_SO4))   dum_sfcocn1(io_SO4,i,j)= fun_calc_SO4tot(dum_sfcocn1(io_S,i,j))
              IF (.NOT. ocn_select(io_F))     dum_sfcocn1(io_F,i,j)  = fun_calc_Ftot(dum_sfcocn1(io_S,i,j))
              ! (3) set ocn->sed flux
              ! NOTE: convert units from (mol per timestep) to (mol m-2 s-1)
              ! NOTE: if 'allow particulate flux to sediments' option in biogem_config is not selected,
              !       the value of the particulate flux to sediments local array <locij_ocnsed> is zero
              DO l=1,n_l_sed
                 is = conv_iselected_is(l)
                 locij_focnsed(is,i,j) = bio_settle(is,i,j,loc_k1)
                 dum_sfxsed1(is,i,j) = phys_ocn(ipo_rA,i,j,loc_k1)*locij_focnsed(is,i,j)*loc_rdts
              end do
              ! (4) set CaCO3 age tracer
              if (sed_select(is_CaCO3_age)) then
                 dum_sfxsed1(is_CaCO3_age,i,j) = loc_t*dum_sfxsed1(is_CaCO3,i,j)
              end if

              IF (ctrl_debug_lvl1 .AND. loc_debug_ij) &
                   & print*, '*** AUDIT FLUX UPDATE ***'
              ! *** AUDIT FLUX UPDATE ***
              IF (ctrl_audit) THEN
                 ! update net integrated fluxes into the ocean for tracer inventory auditing (if selected)
                 ! NOTE: take into account any sediment tracer flux forcing
                 loc_ocnsed_audit(:) = -locij_fsedocn(:,i,j) - locij_frokocn(:,i,j)
                 DO l=1,n_l_sed
                    is = conv_iselected_is(l)
                    loc_tot_i = conv_sed_ocn_i(0,is)
                    do loc_i=1,loc_tot_i
                       io = conv_sed_ocn_i(loc_i,is)
                       loc_ocnsed_audit(io) = loc_ocnsed_audit(io) + &
                            & conv_sed_ocn(io,is)*(locij_focnsed(is,i,j) - loc_dtyr*force_flux_sed(is,i,j))
                       locio_mask(io) = .FALSE.
                    end do
                 end DO
                 DO l=3,n_l_ocn
                    io = conv_iselected_io(l)
                    audit_ocn_delta(io) = audit_ocn_delta(io) + &
                         & loc_dtyr*SUM(locijk_focn(io,i,j,loc_k1:n_k)) - loc_ocnsed_audit(io)
                 END DO
              end IF

              ! Setting NH4_new to NH4 and NO3_new to NO3 below a depth threshold - Fanny (June 2015)
              !DO k=loc_k1,n_k
              !   if (ocn_select(io_NO3).and.ocn_select(io_NH4)) then ! .and. (k<n_k)) then
              !      ocn(io_NH4_new,i,j,k) = ocn(io_NH4,i,j,k)
              !      ocn(io_NO3_new,i,j,k) = ocn(io_NO3,i,j,k)
              !   end if 

                 !Fanny
              !   print*, &
              !   & 'loc_NH4_new =', ocn(io_NH4_new,i,j,k), ocn(io_NH4,i,j,k), i, j,k
              !end DO

           end if

           ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
           ! *** WET GRID PT CONDITIONALITY END *** !
           ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

           IF (ctrl_debug_lvl1 .AND. loc_debug_ij) &
                & print*, '*** <<< END (i,j) GRID POINT'

        END DO block_jloop
     END DO block_iloop

     ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
     ! *** (i,j) GRID PT LOOP END *** !
     ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

     ! *** CALCULATE TRACER ANOMOLY ***
     ! NOTE: also, for now, vectorize <ocn>
     IF (ctrl_debug_lvl1) print*, '*** CALCULATE TRACER ANOMOLY ***'
     do n=1,n_vocn
        loc_i = vocn(n)%i
        loc_j = vocn(n)%j
        loc_k1 = vocn(n)%k1
        vdocn(n)%i = loc_i
        vdocn(n)%j = loc_j
        vdocn(n)%k1 = loc_k1
        DO k=n_k,loc_k1,-1
           DO l=1,n_l_ocn
              io = conv_iselected_io(l)
              vdocn(n)%mk(l,k) = bio_remin(io,loc_i,loc_j,k) + loc_dtyr*vphys_ocn(n)%mk(ipo_rM,k)*locijk_focn(io,i,j,k)
              vocn(n)%mk(l,k) = ocn(io,loc_i,loc_j,k)
           end do
        end DO
     end do
     DO i=1,n_i
        DO j=1,n_j
           DO k=goldstein_k1(i,j),n_k
              !
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 docn(io,i,j,k) = loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
              end DO
              ! adjust particulate tracer concentrations
              DO l=1,n_l_sed
                 is = conv_iselected_is(l)
                 dbio_part(is,i,j,k) = loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_fpart(is,i,j,k)
              end do
           end do
        end do
     end DO

  END IF if_go

  IF (ctrl_debug_lvl1) print*, '*** TEST FOR END-OF-RUN ***'
  ! *** TEST FOR END-OF-RUN ***
  ! if local (BioGeM) model time has surpassed the set model end time, then set the flag updating system biogeochemistry to false
  IF (loc_t < const_real_nullsmall) par_misc_t_go = .FALSE.
  If (error_stop) then
     PRINT*,'*** END BioGeM - ERROR ***'
     PRINT*,' '
     ! reset array values
     CALL sub_init_int_timeslice()
     call diag_biogem_timeslice( &
          dum_dts,               &
          dum_genie_clock,       &
          & dum_sfcatm1,         &
          & dum_sfxatm1,         &
          & dum_sfxocn1,         &
          & dum_sfcsed1,         &
          & dum_sfxsed1          &
          )
     CALL end_biogem()
     STOP
  end If

  return

end subroutine biogem
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! BIOGEM LOOP SUBROUTINE
subroutine biogem_tracercoupling_tmp( &
     & dum_ts,                    &
     & dum_ts1                    &
     & )
  use omp_lib
  USE biogem_lib
  USE biogem_box
  USE genie_util, ONLY: check_iostat
  IMPLICIT NONE
  ! DUMMY ARGUMENTS
  real,intent(inout),dimension(intrac_ocn,n_i,n_j,n_k)::dum_ts   ! NOTE: number of tracers in GOLDSTEIN used in dimension #1
  real,intent(inout),dimension(intrac_ocn,n_i,n_j,n_k)::dum_ts1  ! NOTE: number of tracers in GOLDSTEIN used in dimension #1
  ! LOCAL VARIABLES
  INTEGER::i,j,k,l,io,is                                         ! counting indices
  integer::loc_k1                                                ! local topography
  real::loc_ocn_tot_M,loc_ocn_rtot_M                             ! ocean mass and reciprocal
  real::loc_ocn_tot_V,loc_ocn_rtot_V                             ! ocean volume  and reciprocal
  real::loc_ocn_mean_S_OLD,loc_ocn_rmean_S_OLD                   ! old mean ocean salinity and reciprocal
  real::loc_ocn_mean_S_NEW                                       ! new mean ocean salinity
  real::loc_Sratio,loc_rSratio                                   !
  REAL,DIMENSION(n_ocn)::loc_ocn_tot_OLD,loc_ocn_tot_NEW         !
  REAL,DIMENSION(n_ocn)::loc_ocn_rtot_NEW                        !
  REAL,DIMENSION(n_ocn,n_i,n_j,n_k)::locijk_ocn                  ! local ocean tracer array

  ! *** CALCULATE LOCAL CONSTANTS ***
  ! total ocean mass and recpirocal
  loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
  loc_ocn_rtot_M = 1.0/loc_ocn_tot_M
  ! total ocean volume and recpirocal
  loc_ocn_tot_V = sum(phys_ocn(ipo_V,:,:,:))
  loc_ocn_rtot_V = 1.0/loc_ocn_tot_V

        IF (ctrl_debug_lvl1) print*, '*** OCEAN TRACER UPDATE ***'
        if (ctrl_misc_Snorm) then
           ! [SALINITY NORMALIZED SCHEME]
           ! NOTE: units of aqueous concentration (mol kg-1)
           ! NOTE: loop though only tracers actually selected (i.e., 3:lmax)
           ! calculate mean ocean salinity
           loc_ocn_mean_S_OLD = SUM(ocn(io_S,:,:,:)*phys_ocn(ipo_V,:,:,:))*loc_ocn_rtot_V
           loc_ocn_rmean_S_OLD = 1.0/loc_ocn_mean_S_OLD
           ! calculate unadjusted biogeochem tracer inventory
           DO l=1,n_l_ocn
              io = conv_iselected_io(l)
              if (ocn_type(io) > 0) loc_ocn_tot_OLD(io) = SUM(ocn(io,:,:,:)*phys_ocn(ipo_M,:,:,:))
           end DO
           ! (1) salinity-adjust circulation-updated GOLDSTEIn <ts> tracer array
           ! NOTE: no offset (array: <tstoocn_offset()>) required for biogeochem-only tracers
           DO i=1,n_i
              DO j=1,n_j
                 DO k=goldstein_k1(i,j),n_k
                    DO l=1,n_l_ocn
                       io = conv_iselected_io(l)
                       if (ocn_type(io) > 0) locijk_ocn(io,i,j,k) = dum_ts(l,i,j,k)*(ocn(io_S,i,j,k)*loc_ocn_rmean_S_OLD)
                    end DO
                 end DO
              end do
           end DO
           ! calculate new salinity-adjusted biogeochem tracer inventory
           ! NOTE: be sure to avoid potential divide-by-zero problems
           DO l=1,n_l_ocn
              io = conv_iselected_io(l)
              if (ocn_type(io) > 0) then
                 loc_ocn_tot_NEW(io) = SUM(locijk_ocn(io,:,:,:)*phys_ocn(ipo_M,:,:,:))
                 if (loc_ocn_tot_NEW(io) < const_real_nullsmall) then
                    loc_ocn_rtot_NEW(io) = const_real_zero
                 else
                    loc_ocn_rtot_NEW(io) = 1.0/loc_ocn_tot_NEW(io)
                 end if
              end if
           end do
        elseif (ctrl_misc_noSnorm) then
           ! [NON-SALINITY NORMALIZED SCHEME]
           ! (1) just straight-copy the circulation-updated GOLDSTEIn <ts> tracer array
           ! NOTE: no offset (array: <tstoocn_offset()>) required for biogeochem-only tracers
           DO i=1,n_i
              DO j=1,n_j
                 DO k=goldstein_k1(i,j),n_k
                    DO l=1,n_l_ocn
                       io = conv_iselected_io(l)
                       if (ocn_type(io) > 0) locijk_ocn(io,i,j,k) = dum_ts(l,i,j,k)
                    end DO
                 end DO
              end do
           end DO
        else
           ! [NO OCEAN CIRCULATION]
        endif
        ! (2) update <ocn> T,S
        ! NOTE: includes hack for T,S boundary condition forcing in biogeochemical model only (i.e., not see by GOLDSTEIN)
        if (ctrl_misc_Snorm .OR. ctrl_misc_noSnorm .OR. ctrl_misc_nobioupdate) then
           DO i=1,n_i
              DO j=1,n_j
                 ! set local depth loop limit
                 loc_k1 = goldstein_k1(i,j)
                 DO k=loc_k1,n_k
                    DO l=1,n_l_ocn
                       io = conv_iselected_io(l)
                       if (ocn_type(io) == 0) then
                          If (ctrl_force_GOLDSTEInTS) then
                             dum_ts(l,i,j,k)  = dum_ts(l,i,j,k) + docn(io,i,j,k)
                             dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)
                             ocn(io,i,j,k) = dum_ts(l,i,j,k) + tstoocn_offset(l)
                          else
                             if (force_restore_ocn_select(io)) then
                                ocn(io,i,j,k) = ocn(io,i,j,k) + docn(io,i,j,k)
                             else
                                ocn(io,i,j,k) = dum_ts(l,i,j,k) + tstoocn_offset(l)
                             end if
                          end if
                       end if
                    end DO
                 end do
              end DO
           end DO
        else
           DO i=1,n_i
              DO j=1,n_j
                 DO k=goldstein_k1(i,j),n_k
                    DO l=1,2
                       io = conv_iselected_io(l)
                       ocn(io,i,j,k) = ocn(io,i,j,k) + docn(io,i,j,k)
                    end DO
                 end DO
              end DO
           end DO
        end if
        if (ctrl_misc_Snorm) then
           ! [SALINITY NORMALIZED SCHEME]
           ! re-calculate mean ocean salinity as well as the relative change in global salinity compared to the previous time-step
           ! NOTE: calculate loc_ocn_mean_S_NEW with ocean volume, because it is assumed invarient
           loc_ocn_mean_S_NEW = SUM(ocn(io_S,:,:,:)*phys_ocn(ipo_V,:,:,:))*loc_ocn_rtot_V
           loc_Sratio = loc_ocn_mean_S_NEW/loc_ocn_mean_S_OLD
           loc_rSratio = 1.0/loc_Sratio
           ! (3a) adjust biogeochemal cell masses according to changes in global salinity inventory
           ! (3b)  correct the <ocn> biogeochem tracer inventories and then update with biogeochem fluxes
           !       salinity normalize <ts> biogeochem tracers and copy to <ts1>
           ! NOTE: zero the global remin array after having added its contents to the ocean tracer array
           DO i=1,n_i
              DO j=1,n_j
                 DO k=goldstein_k1(i,j),n_k
                   DO l=1,n_l_ocn
                       io = conv_iselected_io(l)
                       if (ocn_type(io) > 0) then
                          ocn(io,i,j,k) = locijk_ocn(io,i,j,k)*(loc_ocn_tot_OLD(io)*loc_ocn_rtot_NEW(io)) + &
                               & bio_remin(io,i,j,k) + docn(io,i,j,k)
                          ocn(io,i,j,k) = loc_Sratio*ocn(io,i,j,k)
                       end if
                    end do
                    ! adjust particulate tracer concentrations
                    DO l=1,n_l_sed
                       is = conv_iselected_is(l)
                       bio_part(is,i,j,k) = bio_part(is,i,j,k) + dbio_part(is,i,j,k)
                       bio_part(is,i,j,k) = loc_Sratio*bio_part(is,i,j,k)
                    end do

                    ! update GOLDSTEIn tracer arrays and salinity-normalize
                    DO l=1,n_l_ocn
                       io = conv_iselected_io(l)
                       if (ocn_type(io) > 0) then
                          dum_ts(l,i,j,k)  = (loc_ocn_mean_S_NEW/ocn(io_S,i,j,k))*ocn(io,i,j,k)
                          dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)
                       end if
                    end do
                 end DO
              end do
           end DO
           phys_ocn(ipo_M,:,:,:) = loc_rSratio*phys_ocn(ipo_M,:,:,:)
           phys_ocn(ipo_rM,:,:,:) = loc_Sratio*phys_ocn(ipo_rM,:,:,:)
        elseif (ctrl_misc_noSnorm) then
           ! [NON-SALINITY NORMALIZED SCHEME]
!           ! (3) update <ocn> tracers with biogeochem fluxes
!           DO i=1,n_i
!              DO j=1,n_j
!                 DO k=goldstein_k1(i,j),n_k
!                    DO l=1,n_l_ocn
!                       io = conv_iselected_io(l)
!                       if (ocn_type(io) > 0) then
!                          ocn(io,i,j,k) = locijk_ocn(io,i,j,k) + docn(io,i,j,k)
!!!                          ocn(io,i,j,k) = locijk_ocn(io,i,j,k) + &
!!!                               & bio_remin(io,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
!                          !!!! reset remin array
!                          !!!bio_remin(io,i,j,k) = 0.0
!                       end if
!                    end do
!                    ! catch mis-behaving rapidly-oxidizing species going < 0.0
!                    if (ctrl_bio_remin_reminfix) then
!                       if (ocn_select(io_NH4)) call sub_calc_bio_remin_fix_NH4(ocn(:,i,j,k))
!                       if (ocn_select(io_H2S)) call sub_calc_bio_remin_fix_H2S(ocn(:,i,j,k))
!                    end if
!!!                    ! adjust particulate tracer concentrations
!!!                    DO l=1,n_l_sed
!!!                       is = conv_iselected_is(l)
!!!                       bio_part(is,i,j,k) = bio_part(is,i,j,k) + &
!!!                            & loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_fpart(is,i,j,k)
!!!                    end do
!                    ! update GOLDSTEIn tracer arrays
!                    DO l=1,n_l_ocn
!                       io = conv_iselected_io(l)
!                       if (ocn_type(io) > 0) then
!                          dum_ts(l,i,j,k)  = ocn(io,i,j,k)
!                          dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)
!                       end if
!                    end DO
!                 end DO
!              end do
!           end DO
        elseif (ctrl_misc_nobioupdate) then
           ! [NO BIOGEM MODIFICATION OF TRACER FIELD]
!           ! (3) copy  <ocn> tracers 
!          DO i=1,n_i
!              DO j=1,n_j
!                 DO k=goldstein_k1(i,j),n_k
!                    DO l=1,n_l_ocn
!                       io = conv_iselected_io(l)
!                       if (ocn_type(io) > 0) then
!                          ocn(io,i,j,k) = locijk_ocn(io,i,j,k)
!                          ! reset arrays
!                          !!!bio_remin(io,i,j,k)   = 0.0
!!!                          locijk_focn(io,i,j,k) = 0.0
!                       end if
!                    end do
!                    DO l=1,n_l_sed
!                       is = conv_iselected_is(l)
!                       bio_part(is,i,j,k) = 0.0
!                    end DO
!                 end DO
!                 DO l=3,n_l_atm
!                    ia = conv_iselected_ia(l)
!!!                    locij_fatm(ia,i,j) = 0.0
!                 end DO
!              end DO
!           end DO
        else
           ! [NO UPDATE OF OCEAN ARRAY]
!           ! (3) update internal tracer arrays only with biogeochem fluxes
!           DO i=1,n_i
!              DO j=1,n_j
!                 DO k=goldstein_k1(i,j),n_k
!                    DO l=1,n_l_ocn
!                       if (ocn_type(io) > 0) then
!                          io = conv_iselected_io(l)
!                          ocn(io,i,j,k) = ocn(io,i,j,k) + docn(io,i,j,k)
!!!                          ocn(io,i,j,k) = ocn(io,i,j,k) + bio_remin(io,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
!                          !!!! reset remin array
!                          !!!bio_remin(io,i,j,k) = 0.0
!                       end if
!                    end do
!                    ! catch mis-behaving rapidly-oxidizing species going < 0.0
!                    if (ctrl_bio_remin_reminfix) then
!                       if (ocn_select(io_NH4)) call sub_calc_bio_remin_fix_NH4(ocn(:,i,j,k))
!                       if (ocn_select(io_H2S)) call sub_calc_bio_remin_fix_H2S(ocn(:,i,j,k))
!                    end if
!!!                    ! adjust particulate tracer concentrations
!!!                    DO l=1,n_l_sed
!!!                       is = conv_iselected_is(l)
!!!                       bio_part(is,i,j,k) = bio_part(is,i,j,k) + &
!!!                            & loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_fpart(is,i,j,k)
!!!                    end do
!                 end DO
!              end do
!           end DO
        endif

end subroutine biogem_tracercoupling_tmp
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! BIOGEM LOOP SUBROUTINE
! NOTE: <loc_vocn> is the array into which the GOLDSTEIn tracer field <ts> is going to be de-salinity normalized and placed
!       plus a copy of <ts> T,S for completeness
subroutine biogem_tracercoupling( &
     & dum_ts,                    &
     & dum_ts1                    &
     & )
  USE biogem_lib
  USE biogem_box
  USE genie_util, ONLY: check_iostat
  IMPLICIT NONE
  ! DUMMY ARGUMENTS
  real,intent(inout),dimension(intrac_ocn,n_i,n_j,n_k)::dum_ts   ! NOTE: number of tracers in GOLDSTEIN used in dimension #1
  real,intent(inout),dimension(intrac_ocn,n_i,n_j,n_k)::dum_ts1  ! NOTE: number of tracers in GOLDSTEIN used in dimension #1
  ! LOCAL VARIABLES
  INTEGER::k,l,n                                                 ! counting indices
  integer::loc_k1                                                ! local topography
  real::loc_ocn_tot_M,loc_ocn_rtot_M                             ! ocean mass and reciprocal
  real::loc_ocn_tot_V,loc_ocn_rtot_V                             ! ocean volume  and reciprocal
  real::loc_ocn_mean_S_OLD,loc_ocn_rmean_S_OLD                   ! old mean ocean salinity and reciprocal
  real::loc_ocn_mean_S_NEW                                       ! new mean ocean salinity
  real::loc_Sratio,loc_rSratio                                   !
  REAL,DIMENSION(n_l_ocn)::loc_ocn_tot_OLD,loc_ocn_tot_NEW       !
  REAL,DIMENSION(n_l_ocn)::loc_ocn_rtot_NEW                      !
  type(fieldocn),DIMENSION(:),ALLOCATABLE::loc_vocn              ! 
  type(fieldocn),DIMENSION(:),ALLOCATABLE::loc_vts               ! 
		
  ! *** INITIALIZE LOCAL ARRAYS ***
  ! dimension arrays	
  ALLOCATE(loc_vocn(1:n_vocn),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(loc_vts(1:n_vocn),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  do n=1,n_vocn
     allocate(loc_vocn(n)%mk(1:n_l_ocn,1:n_k),STAT=alloc_error)
     call check_iostat(alloc_error,__LINE__,__FILE__)
     allocate(loc_vts(n)%mk(1:n_l_ocn,1:n_k),STAT=alloc_error)
     call check_iostat(alloc_error,__LINE__,__FILE__)
  end do
  ! set grid information and zero tracer values
  loc_vocn(:) = fun_lib_init_vocn()
  loc_vts(:)  = fun_lib_conv_tsTOvocn(dum_ts(:,:,:,:))

  ! *** CALCULATE LOCAL CONSTANTS ***
  ! total ocean mass and recpirocal
  loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
  loc_ocn_rtot_M = 1.0/loc_ocn_tot_M
  ! total ocean volume and recpirocal
  loc_ocn_tot_V = sum(phys_ocn(ipo_V,:,:,:))
  loc_ocn_rtot_V = 1.0/loc_ocn_tot_V

  ! *** OCEAN TRACER UPDATE ***
  If (.NOT. error_stop) then
     IF (ctrl_debug_lvl1) print*, '*** OCEAN TRACER UPDATE ***'
     ! NOTE: units of aqueous concentration (mol kg-1)
     ! (0) ORIGINAL TRACER INVENTORIES 
     !     calculate original (OLD) mean ocean salinity (fom BIOGEM tracer array)
     loc_ocn_mean_S_OLD = 0.0
     do n=1,n_vocn
        loc_k1 = vocn(n)%k1
           loc_ocn_mean_S_OLD = loc_ocn_mean_S_OLD + &
                & sum(vocn(n)%mk(io_S,loc_k1:n_k)*vphys_ocn(n)%mk(ipo_V,loc_k1:n_k))*loc_ocn_rtot_V
     end do
!!$     loc_ocn_mean_S_OLD = SUM(ocn(io_S,:,:,:)*phys_ocn(ipo_V,:,:,:))*loc_ocn_rtot_V
     loc_ocn_rmean_S_OLD = 1.0/loc_ocn_mean_S_OLD
     ! calculate original (OLD) BIOGEM tracer inventories
     do n=1,n_vocn
        loc_k1 = vocn(n)%k1
        DO l=3,n_l_ocn
           loc_ocn_tot_OLD(l) = sum(vocn(n)%mk(l,loc_k1:n_k)*vphys_ocn(n)%mk(ipo_M,loc_k1:n_k))
        end DO
     end do
!!$     DO l=1,n_l_ocn
!!$        io = conv_iselected_io(l)
!!$        if (ocn_type(io) > 0) loc_ocn_tot_OLD(io) = SUM(ocn(io,:,:,:)*phys_ocn(ipo_M,:,:,:))
!!$     end DO
     ! (1) salinity-adjust (circulation-updated) GOLDSTEIn <loc_vts> tracer array
     !     NOTE: copy T,S for completeness
     !     NOTE: no offset (array: <tstoocn_offset()>) required for biogeochem-only tracers
     do n=1,n_vocn
        loc_k1 = vocn(n)%k1
        DO k=n_k,loc_k1,-1
           loc_vocn(n)%mk(1:2,k)       = loc_vts(n)%mk(1:2,k) + tstoocn_offset(1:2)
           loc_vocn(n)%mk(3:n_l_ocn,k) = loc_vts(n)%mk(3:n_l_ocn,k)*vocn(n)%mk(2,k)*loc_ocn_rmean_S_OLD
        end DO
     end do
!!$     DO i=1,n_i
!!$        DO j=1,n_j
!!$           DO k=goldstein_k1(i,j),n_k
!!$              DO l=1,n_l_ocn
!!$                 io = conv_iselected_io(l)
!!$                 if (ocn_type(io) > 0) locijk_ocn(io,i,j,k) = dum_ts(l,i,j,k)*(ocn(io_S,i,j,k)*loc_ocn_rmean_S_OLD)
!!$              end DO
!!$           end DO
!!$        end do
!!$     end DO
     ! calculate new salinity-adjusted biogeochem tracer inventory
     ! NOTE: be sure to avoid potential divide-by-zero problems
     loc_ocn_tot_NEW(:) = 0.0
     do n=1,n_vocn
        loc_k1 = vocn(n)%k1
        DO l=3,n_l_ocn
           loc_ocn_tot_NEW(l) = sum(loc_vocn(n)%mk(l,loc_k1:n_k)*vphys_ocn(n)%mk(ipo_M,loc_k1:n_k))
        end DO
     end do
     DO l=3,n_l_ocn
        if (loc_ocn_tot_NEW(l) < const_real_nullsmall) then
           loc_ocn_rtot_NEW(l) = const_real_zero
        else
           loc_ocn_rtot_NEW(l) = 1.0/loc_ocn_tot_NEW(l)
        end if
     end do
!!$     DO l=1,n_l_ocn
!!$        io = conv_iselected_io(l)
!!$        if (ocn_type(io) > 0) then
!!$           loc_ocn_tot_NEW(io) = SUM(locijk_ocn(io,:,:,:)*phys_ocn(ipo_M,:,:,:))
!!$           if (loc_ocn_tot_NEW(io) < const_real_nullsmall) then
!!$              loc_ocn_rtot_NEW(io) = const_real_zero
!!$           else
!!$              loc_ocn_rtot_NEW(io) = 1.0/loc_ocn_tot_NEW(io)
!!$           end if
!!$        end if
!!$     end do
     ! (2) update <ocn> T,S
     do n=1,n_vocn
        loc_k1 = vocn(n)%k1
        DO k=n_k,loc_k1,-1
           vocn(n)%mk(1:2,k) = loc_vocn(n)%mk(1:2,k) + vdocn(n)%mk(1:2,k)
        end DO
     end do
!!$     DO i=1,n_i
!!$        DO j=1,n_j
!!$           DO k=goldstein_k1(i,j),n_k
!!$              DO l=1,2
!!$                 io = conv_iselected_io(l)
!!$                 ocn(io,i,j,k) = ocn(io,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
!!$              end DO
!!$           end DO
!!$        end DO
!!$     end DO
     ! re-calculate mean ocean salinity as well as the relative change in global salinity compared to the previous time-step
     ! NOTE: calculate loc_ocn_mean_S_NEW with ocean volume, because it is assumed invarient
     loc_ocn_mean_S_NEW = 0.0
     do n=1,n_vocn
        loc_k1 = vocn(n)%k1
           loc_ocn_mean_S_NEW = loc_ocn_mean_S_NEW + &
                & sum(vocn(n)%mk(2,loc_k1:n_k)*vphys_ocn(n)%mk(ipo_V,loc_k1:n_k))*loc_ocn_rtot_V
     end do
!!$     loc_ocn_mean_S_NEW = SUM(ocn(io_S,:,:,:)*phys_ocn(ipo_V,:,:,:))*loc_ocn_rtot_V
     loc_Sratio = loc_ocn_mean_S_NEW/loc_ocn_mean_S_OLD
     loc_rSratio = 1.0/loc_Sratio
     ! (3a) adjust biogeochemal cell masses according to changes in global salinity inventory
     ! (3b)  correct the <ocn> biogeochem tracer inventories and then update with biogeochem fluxes
     !       salinity normalize <ts> biogeochem tracers and copy to <ts1>
     ! NOTE: zero the global remin array after having added its contents to the ocean tracer array
     do n=1,n_vocn
        loc_k1 = vocn(n)%k1
        DO k=n_k,loc_k1,-1
           !
           vocn(n)%mk(3:n_l_ocn,k) = &
                & (loc_ocn_tot_OLD(3:n_l_ocn)*loc_ocn_rtot_NEW(3:n_l_ocn))*vocn(n)%mk(3:n_l_ocn,k) + &
                & vdocn(n)%mk(3:n_l_ocn,k)
           vocn(n)%mk(3:n_l_ocn,k) = loc_Sratio*vocn(n)%mk(3:n_l_ocn,k)

!print*,n,k,vdocn(n)%mk(3:n_l_ocn,k)
!           vocn(n)%mk(3:n_l_ocn,k) = vocn(n)%mk(3:n_l_ocn,k) + vdocn(n)%mk(3:n_l_ocn,k)

!!$!              ! adjust particulate tracer concentrations

           ! update GOLDSTEIn tracer arrays and salinity-normalize
           loc_vts(n)%mk(1:2,k)       = loc_vocn(n)%mk(1:2,k) - tstoocn_offset(1:2)
           loc_vts(n)%mk(3:n_l_ocn,k) = (loc_ocn_mean_S_NEW/vocn(n)%mk(2,k))*vocn(n)%mk(3:n_l_ocn,k)
        end DO
     end do
!!$     DO i=1,n_i
!!$        DO j=1,n_j
!!$           DO k=goldstein_k1(i,j),n_k 
!!$              DO l=1,n_l_ocn
!!$                 io = conv_iselected_io(l)
!!$                 if (ocn_type(io) > 0) then
!!$                    ocn(io,i,j,k) = locijk_ocn(io,i,j,k)*(loc_ocn_tot_OLD(io)*loc_ocn_rtot_NEW(io)) + &
!!$                         & bio_remin(io,i,j,k) + loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_focn(io,i,j,k)
!!$                    ocn(io,i,j,k) = loc_Sratio*ocn(io,i,j,k)
!!$                 end if
!!$              end do
!!$!              ! adjust particulate tracer concentrations
!!$!              DO l=1,n_l_sed
!!$!                 is = conv_iselected_is(l)
!!$!                 bio_part(is,i,j,k) = bio_part(is,i,j,k) + &
!!$!                      & loc_dtyr*phys_ocn(ipo_rM,i,j,k)*locijk_fpart(is,i,j,k)
!!$!                 bio_part(is,i,j,k) = loc_Sratio*bio_part(is,i,j,k)
!!$!              end do
!!$              ! update GOLDSTEIn tracer arrays and salinity-normalize
!!$              DO l=1,n_l_ocn
!!$                 io = conv_iselected_io(l)
!!$                 if (ocn_type(io) > 0) then
!!$                    dum_ts(l,i,j,k)  = (loc_ocn_mean_S_NEW/ocn(io_S,i,j,k))*ocn(io,i,j,k)
!!$                    dum_ts1(l,i,j,k) = dum_ts(l,i,j,k)
!!$                 end if
!!$              end do
!!$           end DO
!!$        end do
!!$     end DO
     do n=1,n_vocn
        loc_k1 = vocn(n)%k1
        vphys_ocn(n)%mk(ipo_M,loc_k1:n_k) = loc_rSratio*vphys_ocn(n)%mk(ipo_M,loc_k1:n_k)
        vphys_ocn(n)%mk(ipo_rM,loc_k1:n_k) = loc_rSratio*vphys_ocn(n)%mk(ipo_rM,loc_k1:n_k)
     end do
!!$     phys_ocn(ipo_M,:,:,:) = loc_rSratio*phys_ocn(ipo_M,:,:,:)
!!$     phys_ocn(ipo_rM,:,:,:) = loc_Sratio*phys_ocn(ipo_rM,:,:,:)

     ! set dummy variable values
  dum_ts(:,:,:,:) = fun_lib_conv_vocnTOts(loc_vts(:))
  dum_ts1(:,:,:,:) = dum_ts(:,:,:,:)
!
!!!  ocn(:,:,:,:) = fun_lib_conv_vocnTOocn(vocn(:))

!!$  DO i=1,n_i
!!$     DO j=1,n_j
!!$        DO k=goldstein_k1(i,j),n_k 
!!$           ! catch mis-behaving rapidly-oxidizing species going < 0.0
!!$           if (ctrl_bio_remin_reminfix) then
!!$              if (ocn_select(io_NH4)) call sub_calc_bio_remin_fix_NH4(ocn(:,i,j,k))
!!$              if (ocn_select(io_H2S)) call sub_calc_bio_remin_fix_H2S(ocn(:,i,j,k))
!!$           end if
!!$        end DO
!!$     end DO
!!$  end DO

  end If

end subroutine biogem_tracercoupling
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! BIOGEM FORCINGS
subroutine biogem_forcing( &
     & dum_genie_clock     &
     & )
  USE biogem_lib
  USE biogem_box
  USE biogem_data
  IMPLICIT NONE
  ! DUMMY ARGUMENTS
  integer(kind=8),INTENT(IN)::dum_genie_clock                    ! genie clock (ms since start) NOTE: 8-byte integer
  ! LOCAL VARIABLES
  INTEGER::l,io,ia,is                                            ! counting indices
  real::loc_t                                                    ! local time

  ! *** CALCULATE GEM TIME ***
  ! update model time
  ! NOTE: par_misc_t_runtime is counted DOWN in years
  !       => for BIOGEM, the 'end of the world' occurs when time reaches zero
  loc_t = par_misc_t_runtime - real(dum_genie_clock)/(1000.0*conv_yr_s)
  ! *** UPDATE FORCING TIME SERIES DATA ***
  IF (ctrl_debug_lvl1) print*, '*** UPDATE FORCING TIME SERIES DATA ***'
  ! recalculate time-varying restoring and flux forcings of the system & fractional relaxation
  ! ATMOSPHERIC TRACERS (applied at the ocean-atmospere interface)
  DO l=3,n_l_atm
     ia = conv_iselected_ia(l)
     IF (force_restore_atm_select(ia)) THEN
        CALL sub_update_force_restore_atm(loc_t,ia)
     END IF
     IF (force_flux_atm_select(ia)) THEN
        CALL sub_update_force_flux_atm(loc_t,ia)
     END IF
  END DO
  ! OCEAN TRACERS
  DO l=1,n_l_ocn
     io = conv_iselected_io(l)
     IF (force_restore_ocn_select(io)) THEN
        CALL sub_update_force_restore_ocn(loc_t,io)
     END IF
     IF (force_flux_ocn_select(io)) THEN
        CALL sub_update_force_flux_ocn(loc_t,io)
     END IF
  END DO
  ! SEDIMENT TRACERS (applied at the ocean surface)
  DO l=1,n_l_sed
     is = conv_iselected_is(l)
     IF (force_flux_sed_select(is)) THEN
        CALL sub_update_force_flux_sed(loc_t,is)
     END IF
  END DO

end subroutine biogem_forcing
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! BIOGEM LOOP SUBROUTINE - CLIMATE STATE UPDATE
     subroutine biogem_climate( &
     & dum_hght_sic,            &
     & dum_frac_sic,            &
     & dum_cost,                &
     & dum_solfor,              &
     & dum_fxsw,                &
     & dum_u,                   &
     & dum_tau,                 &
     & dum_mld,                 &
     & dum_solconst             &
     & )
  USE biogem_lib
  USE biogem_box
  IMPLICIT NONE
  ! DUMMY ARGUMENTS
  real,dimension(n_i,n_j),INTENT(in)::dum_hght_sic               ! sea-ice height (2-D)
  real,dimension(n_i,n_j),INTENT(in)::dum_frac_sic               ! sea-ice fractional cover (2-D)
  real,dimension(n_i,n_j),INTENT(in)::dum_cost                   ! convective frequency something-or-other whatsit thingy (2-D)
  REAL,DIMENSION(n_j),INTENT(in)::dum_solfor                     ! potential solar radiation incident at Earth's surface
  REAL,DIMENSION(n_i,n_j),INTENT(in)::dum_fxsw                   ! actual solar radiation incident at Earth's surface
  REAL,DIMENSION(3,n_i,n_j,n_k),INTENT(in)::dum_u                ! GOLDSTEIN velocity field (3-D)
  REAL,DIMENSION(2,n_i,n_j),INTENT(in)::dum_tau                  ! GOLDSTEIN wind stress field (2-D)
  real,dimension(n_i,n_j),INTENT(in)::dum_mld                    ! mixed layer depth
  REAL,INTENT(inout)::dum_solconst                               ! solar constant

  ! LOCAL VARIABLES
  INTEGER::i,j,k
  integer::loc_k1                                                ! local topography
  real::loc_tau_scale                                            !
  real,dimension(n_i,n_j)::locij_seaice_V_old                    ! sea-ice volume

  ! *** INITIALIZE ***
  ! calculate local variables
  loc_tau_scale = goldstein_rh0sc*goldstein_dsc*goldstein_usc*goldstein_fsc/goldstein_scf
  ! calculate sea-ice volume at previous (BIOGEM) time-step
  DO i=1,n_i
     DO j=1,n_j
        loc_k1 = goldstein_k1(i,j)
        IF (n_k >= loc_k1) THEN
           locij_seaice_V_old(i,j) = phys_ocnatm(ipoa_seaice,i,j)*phys_ocnatm(ipoa_A,i,j)*phys_ocnatm(ipoa_seaice_th,i,j)
        end IF
     end DO
  end DO

  ! *** UPDATE CLIMATE PROPERTIES ***
  ! (i.e., just the physical properties of the ocean and ocean-atmosphere interface that BIOGEM needs to know about)
  DO i=1,n_i
     DO j=1,n_j
        loc_k1 = goldstein_k1(i,j)
        IF (n_k >= loc_k1) THEN
           ! solar insolation
           phys_ocnatm(ipoa_solfor,i,j) = dum_solfor(j)
           phys_ocnatm(ipoa_fxsw,i,j)   = dum_fxsw(i,j)
           ! wind stress
           phys_ocnatm(ipoa_tau_u,i,j) = loc_tau_scale*dum_tau(1,i,j)
           phys_ocnatm(ipoa_tau_v,i,j) = loc_tau_scale*dum_tau(2,i,j)
           ! convective 'cost'
           phys_ocnatm(ipoa_cost,i,j) = dum_cost(i,j)
           ! MLD
           phys_ocnatm(ipoa_mld,i,j) = dum_mld(i,j)
        end IF
        do k = loc_k1,n_k
           ! density
           phys_ocn(ipo_rho,i,j,k) = fun_calc_rho(ocn(io_T,i,j,k),ocn(io_S,i,j,k))
           ! velocity
           phys_ocn(ipo_gu,i,j,k) = dum_u(1,i,j,k)
           phys_ocn(ipo_gv,i,j,k) = dum_u(2,i,j,k) 
           phys_ocn(ipo_gw,i,j,k) = dum_u(3,i,j,k) 
        end do
     end DO
  end DO
  ! force to prescribed fractional sea-ice cover if requested
  if (ctrl_force_seaice) then
     phys_ocnatm(ipoa_seaice,:,:) = par_phys_seaice(:,:)
  else
     phys_ocnatm(ipoa_seaice,:,:) = dum_frac_sic
  end if
  phys_ocnatm(ipoa_seaice_th,:,:) = dum_hght_sic
  ! calculate change in sea-ice volume
  DO i=1,n_i
     DO j=1,n_j
        loc_k1 = goldstein_k1(i,j)
        IF (n_k >= loc_k1) THEN
           phys_ocnatm(ipoa_seaice_dV,i,j) = &
                & phys_ocnatm(ipoa_seaice,i,j)*phys_ocnatm(ipoa_A,i,j)*phys_ocnatm(ipoa_seaice_th,i,j) - &
                & locij_seaice_V_old(i,j)
        end IF
     end DO
  end DO
  ! force to prescribed wind-speed if requested
  ! otherwise, calculate wind speed from wind stress
  if (ctrl_force_windspeed) then
     phys_ocnatm(ipoa_u,:,:) = par_phys_windspeed(:,:)
  else
     phys_ocnatm(ipoa_u,:,:) = fun_calc_u()
  end if
  ! replace solar constant(!)
  if (ctrl_force_solconst) dum_solconst = force_solconst_sig(2,par_data_save_sig_i)
  ! make internal copy of solar constant
  phys_solar_constant = dum_solconst

end subroutine biogem_climate
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! BIOGEM pCO2 diagnostics (for GEMlite)
subroutine diag_biogem_pCO2( &
     & dum_sfcatm1,          &
     & dum_pCO2              &
     & )
  USE biogem_lib
  IMPLICIT NONE                     
  ! DUMMY ARGUMENTS
  REAL,DIMENSION(n_atm,n_i,n_j),INTENT(in)::dum_sfcatm1          ! 
  REAL,INTENT(inout)::dum_pCO2                                   ! 
  ! calculate pCO2
  dum_pCO2 = conv_mol_umol*SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia_pCO2,:,:))/sum(phys_ocnatm(ipoa_A,:,:))
end subroutine diag_biogem_pCO2
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! RESTART BioGeM (save data)
SUBROUTINE rest_biogem()
  USE biogem_lib
  USE genie_util, ONLY:check_unit,check_iostat
  IMPLICIT NONE
  ! local variables
  integer::l                                                     ! tracer counter
  CHARACTER(len=255)::loc_filename                               ! guess ...
  integer::ios                                                   ! file checks
  INTEGER::i,j,io,is                                             ! counting indices
  integer::loc_k1                                                ! local topography
  integer::loc_i,loc_tot_i                                       ! tracer array conversion indices
  REAL,DIMENSION(n_ocn,n_i,n_j,n_k)::loc_docn                    ! local ocean tracer (change) array

  ! *** ENSURE TRACER CONSERVATION ***
  ! 'flush' settling solid sediment interface array
  loc_docn(:,:,:,:) = 0.0
  if (ctrl_force_sed_closedsystem) then
     DO i=1,n_i
        DO j=1,n_j
           loc_k1 = goldstein_k1(i,j)
           IF (n_k >= loc_k1) THEN
              DO l=1,n_l_sed
                 is = conv_iselected_is(l)
                 loc_tot_i = conv_sed_ocn_i(0,is)
                 do loc_i=1,loc_tot_i
                    ! convert settling solids at sediment interface to dissolved constituents
                    io = conv_sed_ocn_i(loc_i,is)
                    loc_docn(io,i,j,loc_k1) = loc_docn(io,i,j,loc_k1) + &
                         & phys_ocn(ipo_rM,i,j,loc_k1)*conv_sed_ocn(io,is)*bio_settle(is,i,j,loc_k1)
                    ! prevent return of dissolved Fe
                    loc_docn(io_Fe,i,j,loc_k1) = 0.0
                 end do
              end DO
           end if
        end DO
     end DO
  end if
  ! correct tracer inventory for restart save
  ocn(:,:,:,:) = ocn(:,:,:,:) + loc_docn(:,:,:,:)
  bio_settle(:,:,:,:) = 0.0
  ! *** DUMP RESTART DATA ***
  ! NOTE: data is saved unformatted for minimal file size
  !       also means that arrays can be written directly to file without needing to loop thought data
  loc_filename = TRIM(par_outdir_name)//trim(par_outfile_name)
  call check_unit(out,__LINE__,__FILE__)
  OPEN(unit=out,status='replace',file=loc_filename,form='unformatted',action='write',iostat=ios)
  call check_iostat(ios,__LINE__,__FILE__)
  WRITE(unit=out,iostat=ios) &
       & n_l_ocn,                                           &
       & (conv_iselected_io(l),l=1,n_l_ocn),                &
       & (ocn(conv_iselected_io(l),:,:,:),l=1,n_l_ocn),     &
       & n_l_sed,                                           &
       & (conv_iselected_is(l),l=1,n_l_sed),                &
       & (bio_part(conv_iselected_is(l),:,:,:),l=1,n_l_sed)
  call check_iostat(ios,__LINE__,__FILE__)
  close(unit=out,iostat=ios)
  call check_iostat(ios,__LINE__,__FILE__)

end SUBROUTINE rest_biogem
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! BIOGEM DIAGNOSTICS
subroutine diag_biogem( &
     & dum_genie_clock, &
     & dum_sfcatm1,     &
     & dum_gemlite      &
     & )
  USE biogem_lib
  USE biogem_box
  USE biogem_data
  IMPLICIT NONE
  ! DUMMY ARGUMENTS
  integer(kind=8),INTENT(IN)::dum_genie_clock                           ! genie clock (milliseconds since start) * 8-byte integer *
  real,intent(in),dimension(n_atm,n_i,n_j)::dum_sfcatm1                 ! atmosphere-surface tracer composition; ocn grid
  logical,intent(in)::dum_gemlite                                       ! in GEMlite phase of cycle?
  ! LOCAL VARIABLES
  real::loc_t,loc_yr                                                    ! local time and time step BLAH actual year
  REAL,DIMENSION(0:n_j,0:n_k)::loc_opsi                                 !
  REAL,DIMENSION(0:n_j,0:n_k)::loc_opsia,loc_opsip,loc_zpsi             !
  real::loc_opsi_scale                                                  !
  REAL,DIMENSION(2)::loc_opsia_minmax,loc_opsip_minmax                  !  

  ! calculate local variables
  loc_t = par_misc_t_runtime - real(dum_genie_clock)/(1000.0*conv_yr_s)
  ! calculate actual year (counting years Before Present or otherwise)
  IF (ctrl_misc_t_BP) THEN
     loc_yr = loc_t + par_misc_t_end
  ELSE
     loc_yr = par_misc_t_end - loc_t
  END IF
  ! calculate local opsi conversion constant
  loc_opsi_scale = goldstein_dsc*goldstein_usc*const_rEarth*1.0E-6

  ! *** RUN-TIME REPORTING ***
  ! run-time reporting
  ! NOTE: carry out an updated tracer inventory audit at this time (if selected)
  ! NOTE: adjusted reporting frequency to every year
  if_report: if ( &
       & ((abs(loc_yr - real(INT(loc_yr))) < conv_s_yr) .AND. ctrl_debug_lvl0) &
       & .OR. &
       & (error_stop) &
       & ) then
     IF (ctrl_debug_lvl1) print*, '*** RUN-TIME REPORTING ***'
     ! run-time data echo-ing
     ! ### UN-COMMENT TO PERIODICALLY RE-PRINT HEADER INFORMATION ############################################################### !
     ! if (mod(par_data_save_sig_i,50) == 0) par_misc_t_echo_header = .TRUE.
     ! ########################################################################################################################## !
     if (error_stop) par_misc_t_echo_header = .TRUE.
     CALL sub_calc_psi(phys_ocn(ipo_gu:ipo_gw,:,:,:),loc_opsi,loc_opsia,loc_opsip,loc_zpsi,loc_opsia_minmax,loc_opsip_minmax)
     call sub_echo_runtime(loc_yr,loc_opsi_scale,loc_opsia_minmax,dum_sfcatm1(:,:,:),dum_gemlite)
     ! carry out tracer audit + echo max,min ocean tracer values (and location)
     ! NOTE: use audit switch to turn on/off
     IF (ctrl_audit) THEN
        CALL sub_echo_maxmin()
        CALL sub_audit_update()
     END IF
  end IF if_report

end subroutine diag_biogem
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! biogem DIAGNOSTICS - TIME-SLICE
SUBROUTINE diag_biogem_timeslice( &
     & dum_dts,                   &
     & dum_genie_clock,           &
     & dum_sfcatm1,               &
     & dum_sfxatm1,               &
     & dum_sfxocn1,               &
     & dum_sfcsed1,               &
     & dum_sfxsed1,               &
     & dum_sfxsumrok1,            &
     & dum_save,                  &
     & dum_gemlite                &
     & )
  USE biogem_lib
  USE biogem_box
  USE biogem_data
  ! dummy arguments
  REAL,INTENT(IN)::dum_dts                                       ! biogem time-step length (seconds)
  integer(kind=8),INTENT(IN)::dum_genie_clock                    ! genie clock (milliseconds since start) NOTE: 8-byte integer
  REAL,DIMENSION(n_atm,n_i,n_j),INTENT(in)::dum_sfcatm1          ! atmosphere composition interface array
  REAL,DIMENSION(n_atm,n_i,n_j),INTENT(in)::dum_sfxatm1          ! atmospheric flux interface array
  REAL,DIMENSION(n_ocn,n_i,n_j),INTENT(in)::dum_sfxocn1          ! sediment dissolution flux interface array
  REAL,DIMENSION(n_sed,n_i,n_j),INTENT(in)::dum_sfcsed1          ! sediment composition interface array
  REAL,DIMENSION(n_sed,n_i,n_j),INTENT(in)::dum_sfxsed1          ! sediment rain flux interface array
  real,dimension(n_ocn,n_i,n_j),intent(in)::dum_sfxsumrok1       ! coastal (weathering) -> ocean flux; ocn grid
  logical,INTENT(IN)::dum_save                                   ! average and save data?
  logical,intent(in)::dum_gemlite                                ! in GEMlite phase of cycle?
  ! local variables
  INTEGER::i,j,k,l,io,ia,is  
  integer::loc_k1                                                ! 
  real::loc_t,loc_dts,loc_dtyr                                   ! 
  real::loc_yr_save                                              ! 
  REAL,DIMENSION(0:n_j,0:n_k)::loc_opsi,loc_zpsi,loc_opsia,loc_opsip !
  REAL,DIMENSION(n_atm,n_i,n_j)::locij_focnatm                   ! local ocn->atm flux (atm tracer currency) (mol yr-1)
  REAL,DIMENSION(n_sed,n_i,n_j)::locij_focnsed                   ! local ocn->sed change (sed tracer currency) (mol)
  REAL,DIMENSION(n_ocn,n_i,n_j)::locij_fsedocn                   ! local sed->ocean change (ocn tracer currency) (mol)
  REAL,DIMENSION(2)::loc_opsia_minmax,loc_opsip_minmax           !

  ! *** TIME-SLICE DATA UPDATE ***
  IF (ctrl_debug_lvl1) print*, '*** TIME-SLICE DATA UPDATE ***'
  ! update time slice data
  ! NOTE: carried out only when the local (BioGeM) time falls between a selected time slice time plus integration time,
  !       and the time slice time itself
  ! NOTE: do not construct and save a time slice if <par_data_save_timeslice_i> == 0,
  !       i.e., no valid (in-range) time slices have been requested in the time slice input file 'biogem_save_timeslice.dat',
  !       or the final requested time slice has been made
  ! calculate local model time and time step length

  loc_t = par_misc_t_runtime - real(dum_genie_clock)/(1000.0*conv_yr_s)
  loc_dts  = dum_dts
  loc_dtyr = loc_dts/conv_yr_s
  if_save1: if (par_data_save_timeslice_i > 0) then
     if_save2: IF ( &
          & ((loc_t - (par_data_save_timeslice(par_data_save_timeslice_i) + par_data_save_slice_dt/2.0)) < -conv_s_yr) &
          & .OR. &
          & (error_stop) &
          & ) THEN

        int_t_timeslice       = int_t_timeslice       + loc_dtyr
        int_t_timeslice_TOT   = int_t_timeslice_TOT   + loc_dtyr
        int_t_timeslice_count = int_t_timeslice_count + 1

        if_save3: if (dum_save) then

           ! reconstruct local interface fluxes and update whole-ocean carbonate equilibrium
           IF (opt_select(iopt_select_carbchem)) THEN
              DO i=1,n_i
                 DO j=1,n_j
                    loc_k1 = goldstein_k1(i,j)
                    IF (n_k >= loc_k1) THEN
                       ! ocn->atm
                       ! NOTE: convert units from (mol m-2 s-1) to (mol yr-1)
                       DO l=3,n_l_atm
                          ia = conv_iselected_ia(l)
                          locij_focnatm(ia,i,j) = conv_yr_s*phys_ocnatm(ipoa_A,i,j)*dum_sfxatm1(ia,i,j)
                       end do
                       ! ocn->sed
                       ! NOTE: convert units from (mol m-2 s-1) to (mol per timestep)
                       DO l=1,n_l_sed
                          is = conv_iselected_is(l)
                          locij_focnsed(is,i,j) = loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxsed1(is,i,j)
                       end DO
                       ! sed->ocn
                       ! NOTE: convert units from (mol m-2 s-1) to (mol per timestep)
                       DO l=3,n_l_ocn
                          io = conv_iselected_io(l)
                          locij_fsedocn(io,i,j) = loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io,i,j)
                       end do
                    end IF
                    DO k=goldstein_k1(i,j),n_k
                       ! calculate carbonate dissociation constants
                       CALL sub_calc_carbconst(           &
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
                       end if
                       ! re-estimate Ca and borate concentrations from salinity (if not selected and therefore explicitly treated)
                       IF (.NOT. ocn_select(io_Ca))  ocn(io_Ca,i,j,k)  = fun_calc_Ca(ocn(io_S,i,j,k))
                       IF (.NOT. ocn_select(io_B))   ocn(io_B,i,j,k)   = fun_calc_Btot(ocn(io_S,i,j,k))
                       IF (.NOT. ocn_select(io_SO4)) ocn(io_SO4,i,j,k) = fun_calc_SO4tot(ocn(io_S,i,j,k))
                       IF (.NOT. ocn_select(io_F))   ocn(io_F,i,j,k)   = fun_calc_Ftot(ocn(io_S,i,j,k))
                       ! re-calculate surface ocean carbonate chemistry
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
                            & ocn(io_DIC,i,j,n_k),  &
                            & ocn(io_ALK,i,j,n_k),  &
                            & ocn(io_PO4,i,j,n_k),  &
                            & ocn(io_SiO2,i,j,n_k), &
                            & ocn(io_B,i,j,n_k),    &
                            & ocn(io_SO4,i,j,n_k),  &
                            & ocn(io_F,i,j,n_k),    &
                            & ocn(io_H2S,i,j,n_k),  &
                            & ocn(io_NH4,i,j,n_k),  &
                            & carbconst(:,i,j,n_k), & 
                            & carb(:,i,j,n_k)       & 
                            & )
                       ! re-calculate carbonate system isotopic properties
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
                    end do
                 end DO
              end DO
           end IF
           ! update time slice data - ocean
           ! NOTE: do not time increment weight quantities such as <int_bio_remin_timeslice> or <int_bio_settle_timeslice>,
           !       because they represent discrete increments rather than a flux or weightable concentration value
           int_ocn_timeslice(:,:,:,:)        = int_ocn_timeslice(:,:,:,:)        + loc_dtyr*ocn(:,:,:,:)
           int_bio_part_timeslice(:,:,:,:)   = int_bio_part_timeslice(:,:,:,:)   + loc_dtyr*bio_part(:,:,:,:)
           int_bio_settle_timeslice(:,:,:,:) = int_bio_settle_timeslice(:,:,:,:) + bio_settle(:,:,:,:)
           int_bio_remin_timeslice(:,:,:,:)  = int_bio_remin_timeslice(:,:,:,:)  + bio_remin(:,:,:,:)
           int_phys_ocn_timeslice(:,:,:,:)   = int_phys_ocn_timeslice(:,:,:,:)   + loc_dtyr*phys_ocn(:,:,:,:)
           int_carb_timeslice(:,:,:,:)       = int_carb_timeslice(:,:,:,:)       + loc_dtyr*carb(:,:,:,:)
           int_carbconst_timeslice(:,:,:,:)  = int_carbconst_timeslice(:,:,:,:)  + loc_dtyr*carbconst(:,:,:,:)
           int_carbisor_timeslice(:,:,:,:)   = int_carbisor_timeslice(:,:,:,:)   + loc_dtyr*carbisor(:,:,:,:)
           ! update time slice data - ocean-atmosphere interface
           int_sfcatm1_timeslice(:,:,:)      = int_sfcatm1_timeslice(:,:,:)     + loc_dtyr*dum_sfcatm1(:,:,:)
           int_focnatm_timeslice(:,:,:)      = int_focnatm_timeslice(:,:,:)     + loc_dtyr*locij_focnatm(:,:,:)
           int_phys_ocnatm_timeslice(:,:,:)  = int_phys_ocnatm_timeslice(:,:,:) + loc_dtyr*phys_ocnatm(:,:,:)
           ! update time slice data - ocean-sediment interface
           int_sfcsed1_timeslice(:,:,:)  = int_sfcsed1_timeslice(:,:,:) + loc_dtyr*dum_sfcsed1(:,:,:)
           int_focnsed_timeslice(:,:,:)  = int_focnsed_timeslice(:,:,:) + locij_focnsed(:,:,:)
           int_fsedocn_timeslice(:,:,:)  = int_fsedocn_timeslice(:,:,:) + locij_fsedocn(:,:,:)
           ! update time-slice data - GOLDSTEIn
           CALL sub_calc_psi(phys_ocn(ipo_gu:ipo_gw,:,:,:),loc_opsi,loc_opsia,loc_opsip,loc_zpsi,loc_opsia_minmax,loc_opsip_minmax)
           int_opsi_timeslice(:,:)  = int_opsi_timeslice(:,:)  + loc_dtyr*loc_opsi(:,:)
           int_opsia_timeslice(:,:) = int_opsia_timeslice(:,:) + loc_dtyr*loc_opsia(:,:)
           int_opsip_timeslice(:,:) = int_opsip_timeslice(:,:) + loc_dtyr*loc_opsip(:,:)
           int_zpsi_timeslice(:,:)  = int_zpsi_timeslice(:,:)  + loc_dtyr*loc_zpsi(:,:)
           int_u_timeslice(:,:,:,:) = int_u_timeslice(:,:,:,:) + loc_dtyr*phys_ocn(ipo_gu:ipo_gw,:,:,:)
           ! update time-slice data - diagnostics
           int_diag_bio_timeslice(:,:,:)       = int_diag_bio_timeslice(:,:,:)       + diag_bio(:,:,:)
           int_diag_geochem_timeslice(:,:,:,:) = int_diag_geochem_timeslice(:,:,:,:) + diag_geochem(:,:,:,:)
                 if (dum_gemlite) then
           int_diag_weather_timeslice(:,:,:)   = int_diag_weather_timeslice(:,:,:)   + loc_dtyr*dum_sfxsumrok1(:,:,:)
                 else
           int_diag_weather_timeslice(:,:,:)   = int_diag_weather_timeslice(:,:,:)   + dum_sfxsumrok1(:,:,:)
                 end if
           int_diag_airsea_timeslice(:,:,:)    = int_diag_airsea_timeslice(:,:,:)    + diag_airsea(:,:,:)

        end if if_save3

        ! write time-slice data and re-set integration
        if ( ((par_data_save_slice_dt - int_t_timeslice) < conv_s_yr) &
             & .OR.                                                                &
             & (error_stop)                                                        &
             & .OR.                                                                &
             & (int_t_timeslice_count == par_data_save_slice_n)                    &
             & ) then

           if_save3b: if (dum_save) then

           ! set save time
           if (int_t_timeslice_count == par_data_save_slice_n) then
              if (ctrl_misc_t_BP) then
                 loc_yr_save = loc_t + int_t_timeslice/2.0 - par_misc_t_end
              else
                 loc_yr_save = par_misc_t_end - loc_t - int_t_timeslice/2.0
              end if
           elseif ((par_data_save_slice_dt - int_t_timeslice) < conv_s_yr) then
              if (ctrl_misc_t_BP) then
                 loc_yr_save = loc_t + int_t_timeslice/2.0 - par_misc_t_end
              else
                 loc_yr_save = par_misc_t_end - loc_t - int_t_timeslice/2.0
              end if
           elseif (error_stop) then
              if (ctrl_misc_t_BP) then
                 loc_yr_save = loc_t - par_misc_t_end
              else
                 loc_yr_save = par_misc_t_end - loc_t
              end if
           else
              ! NOTHING
           end if
           ! clean up year
           ! NOTE: original 'rounding' equation resulted in integer time > 2 Myr ended up over the precision limit
           !             (-2147483.648 ended up being reported)
           loc_yr_save = &
                & real(int(loc_yr_save)) + real(int(1000.0*(loc_yr_save - real(int(loc_yr_save)) + 0.0005)))/1000.0
			   
              ! reporting
              if (int_t_timeslice_count == par_data_save_slice_n) then
           WRITE(unit=6,fmt='(A57,f12.3)') &
               & ' >>> SAVING BIOGEM MONTHLY DATA & year :                 ', &
               & loc_yr_save
              elseif ((par_data_save_slice_dt - int_t_timeslice) < conv_s_yr) then
           WRITE(unit=6,fmt='(A57,f12.3)') &
               & ' >>> SAVING BIOGEM TIME-SLICE AVERAGE CENTERED @ year  : ', &
               & loc_yr_save
              elseif (error_stop) then
           WRITE(unit=6,fmt='(A57,f12.3)') &
               & ' >>> SAVING FATAL ERROR DATA DUMP @ year :               ', &
               & loc_yr_save
              else
                 ! NOTHING
              end if
              ! save time-slice data
              IF (ctrl_data_save_slice_ascii) THEN
                 CALL sub_data_save_timeslice()
                 CALL sub_data_save_timeslice_sed()
              END IF
              ! re-open netcdf file
              call sub_save_netcdf(loc_yr_save,2)
              call sub_save_netcdf(loc_yr_save,3)
              CALL sub_save_netcdf_2d()
              CALL sub_save_netcdf_3d()
              CALL sub_save_netcdf_2d_sed()
              ! close netcdf file and update record number
              call sub_closefile(ncout2d_iou)
              call sub_closefile(ncout3d_iou)
              ncout2d_ntrec = ncout2d_ntrec + 1
              ncout3d_ntrec = ncout3d_ntrec + 1
              ! save global diagnostics
              If (ctrl_data_save_GLOBAL) call sub_data_save_global_av()
              If (ctrl_data_save_GLOBAL) call sub_data_save_global_snap(loc_t,dum_sfcatm1(:,:,:))
              ! reset array values
              CALL sub_init_int_timeslice()

           end if if_save3b

              ! update time slice index at end of primary integration interval (not the end of monthly or seasonal averages)
              if ((par_data_save_slice_dt - int_t_timeslice_TOT) < conv_s_yr) then
                 int_t_timeslice_TOT = 0.0
                 par_data_save_timeslice_i = par_data_save_timeslice_i - 1
              end if

        END IF

     end if if_save2
  end if if_save1

end SUBROUTINE diag_biogem_timeslice
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! BioGeM DIAGNOSTICS - TIME-SERIES
SUBROUTINE diag_biogem_timeseries( &
     & dum_dts,                    &
     & dum_genie_clock,            &
     & dum_sfcatm1,                &
     & dum_sfxatm1,                &
     & dum_sfxocn1,                &
     & dum_sfcsed1,                &
     & dum_sfxsed1,                &
     & dum_sfxsumrok1,             &
     & dum_save,                   &
     & dum_gemlite                 &
     & )
  USE biogem_lib
  USE biogem_box
  USE biogem_data
  ! dummy arguments
  REAL,INTENT(IN)::dum_dts                                       ! biogem time-step length (seconds)
  integer(kind=8),INTENT(IN)::dum_genie_clock                    ! genie clock (milliseconds since start) NOTE: 8-byte integer
  REAL,DIMENSION(n_atm,n_i,n_j),INTENT(in)::dum_sfcatm1          ! atmosphere composition interface array
  REAL,DIMENSION(n_atm,n_i,n_j),INTENT(in)::dum_sfxatm1          ! atmospheric flux interface array
  REAL,DIMENSION(n_ocn,n_i,n_j),INTENT(in)::dum_sfxocn1          ! sediment dissolution flux interface array
  REAL,DIMENSION(n_sed,n_i,n_j),INTENT(in)::dum_sfcsed1          ! sediment composition interface array
  REAL,DIMENSION(n_sed,n_i,n_j),INTENT(in)::dum_sfxsed1          ! sediment rain flux interface array
  real,dimension(n_ocn,n_i,n_j),intent(in)::dum_sfxsumrok1       ! coastal (weathering) -> ocean flux; ocn grid
  logical,INTENT(IN)::dum_save                                   ! average and save data?
  logical,intent(in)::dum_gemlite                                ! in GEMlite phase of cycle?
  ! local variables
  INTEGER::i,j,l,io,ia,is,ic 
  integer::ib,id                                                 ! counting variables
  integer::loc_k1                                                ! 
  real::loc_t,loc_dts,loc_dtyr                                   !
  real::loc_yr,loc_yr_save                                       !
  REAL,DIMENSION(0:n_j,0:n_k)::loc_opsi,loc_zpsi                 !
  REAL,DIMENSION(0:n_j,0:n_k)::loc_opsia,loc_opsip               !
  REAL,DIMENSION(n_atm,n_i,n_j)::locij_focnatm                   ! local ocn->atm flux (atm tracer currency) (mol yr-1)
  REAL,DIMENSION(n_sed,n_i,n_j)::locij_focnsed                   ! local ocn->sed change (sed tracer currency) (mol)
  REAL,DIMENSION(n_ocn,n_i,n_j)::locij_fsedocn                   ! local sed->ocean change (ocn tracer currency) (mol)
  REAL,DIMENSION(n_ocn,n_i,n_j)::locij_ocn_ben                   ! local benthic ocean composition
  REAL,DIMENSION(n_i,n_j)::locij_mask_ben                        ! benthic save mask
  real::loc_ocn_tot_M,loc_ocn_tot_A,loc_ocnatm_tot_A             ! 
  real::loc_ocn_rtot_M,loc_ocn_rtot_A,loc_ocnatm_rtot_A          ! 
  real::loc_ocnsed_tot_A,loc_ocnsed_tot_A_ben                    ! 
  real::loc_ocnsed_rtot_A,loc_ocnsed_rtot_A_ben                  ! 
  real::loc_tot_A                                                ! 
  real::loc_sig                                                  ! 
  REAL,DIMENSION(2)::loc_opsia_minmax,loc_opsip_minmax                  !
  real::loc_opsi_scale                                                  !

  ! *** TIME-SERIES DATA UPDATE ***
  IF (ctrl_debug_lvl1) print*, '*** RUN-TIME DATA UPDATE ***'
  ! update time slice data
  ! NOTE: carried out only when the local (BioGeM) time falls between a selected time slice time plus integration time,
  !       and the time slice time itself
  ! NOTE: do not construct and save a time slice if <par_data_save_timeslice_i> == 0,
  !       i.e., no valid (in-range) time slices have been requested in the time slice input file 'biogem_save_timeslice.dat',
  !       or the final requested time slice has been made
  ! calculate local model time and time step length
  loc_t = par_misc_t_runtime - real(dum_genie_clock)/(1000.0*conv_yr_s)
  loc_dts  = dum_dts
  loc_dtyr = loc_dts/conv_yr_s
  if_save1: if( par_data_save_sig_i > 0 ) then
     if_save2: IF ((loc_t - (par_data_save_sig(par_data_save_sig_i) + par_data_save_sig_dt/2.0)) < -conv_s_yr) THEN

        if_save3: if (dum_save) then

           ! calculate local constants
           ! total ocean mass and recpirocal
           loc_ocn_tot_M = sum(phys_ocn(ipo_M,:,:,:))
           loc_ocn_rtot_M = 1.0/loc_ocn_tot_M
           ! total ocean-atmosphere interface area
           loc_ocnatm_tot_A = sum(phys_ocnatm(ipoa_A,:,:))
           loc_ocnatm_rtot_A = 1.0/loc_ocnatm_tot_A
           ! total ocean surface area (ice-free)
           loc_ocn_tot_A = sum((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k))
           loc_ocn_rtot_A = 1.0/loc_ocn_tot_A
           ! total ocean-sediment interface area
           loc_ocnsed_tot_A = sum(phys_ocn(ipo_A,:,:,n_k))
           loc_ocnsed_rtot_A = 1.0/loc_ocnsed_tot_A
           ! reconstruct local interface fluxes
           ! also: benthic mask
           DO i=1,n_i
              DO j=1,n_j
                 loc_k1 = goldstein_k1(i,j)
                 IF (n_k >= loc_k1) THEN
                    ! ocn->atm
                    ! NOTE: convert units from (mol m-2 s-1) to (mol yr-1)
                    DO l=3,n_l_atm
                       ia = conv_iselected_ia(l)
                       locij_focnatm(ia,i,j) = conv_yr_s*phys_ocnatm(ipoa_A,i,j)*dum_sfxatm1(ia,i,j)
                    end do
                    ! ocn->sed
                    ! NOTE: convert units from (mol m-2 s-1) to (mol per timestep)
                    DO l=1,n_l_sed
                       is = conv_iselected_is(l)
                       locij_focnsed(is,i,j) = loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxsed1(is,i,j)
                    end DO
                    ! sed->ocn
                    ! NOTE: convert units from (mol m-2 s-1) to (mol per timestep)
                    DO l=3,n_l_ocn
                       io = conv_iselected_io(l)
                       locij_fsedocn(io,i,j) = loc_dts*phys_ocn(ipo_A,i,j,loc_k1)*dum_sfxocn1(io,i,j)
                    end do
                    ! interface
                    if (phys_ocn(ipo_Dbot,i,j,loc_k1) > par_data_save_ben_Dmin) then
                       DO l=1,n_l_ocn
                          io = conv_iselected_io(l)
                          locij_ocn_ben(io,i,j) = ocn(io,i,j,loc_k1)
                       end do
                       locij_mask_ben(i,j) = 1.0
                    else
                       locij_mask_ben(i,j) = 0.0
                       locij_ocn_ben(:,i,j) = 0.0
                    end if
                 end IF
              end DO
           end DO
           ! benthic mask area
           loc_ocnsed_tot_A_ben = sum(locij_mask_ben(:,:)*phys_ocn(ipo_A,:,:,n_k))
           loc_ocnsed_rtot_A_ben = 1.0/loc_ocnsed_tot_A_ben
           ! update signal-averaging ingetrated tracers
           ! NOTE: mass-weight ocean tracers, and surface area-weight atmospheric tracers
           ! NOTE: for sea surface (_sur) tracer values weight by fractional ice-free cover
           ! NOTE: for core-top sediment properties filter out grid points with no solid material from isotopic value averaging
           ! NOTE: for calculating mean of sedimentary core-top isotopic composition,
           !       additioanlly weight the integrated area used to normalize the sum by (wt%) CaCO3 abundance
           int_ocn_tot_M_sig     = int_ocn_tot_M_sig     + loc_dtyr*SUM(phys_ocn(ipo_M,:,:,:))
           int_ocn_tot_M_sur_sig = int_ocn_tot_M_sur_sig + loc_dtyr*SUM(phys_ocn(ipo_M,:,:,n_k))
           int_ocn_tot_V_sig     = int_ocn_tot_M_sig     + loc_dtyr*SUM(phys_ocn(ipo_V,:,:,:))
           ! main time-series
           IF (ctrl_data_save_sig_ocn) THEN
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 int_ocn_sig(io) = int_ocn_sig(io) + &
                      & loc_dtyr*SUM(phys_ocn(ipo_M,:,:,:)*ocn(io,:,:,:))*loc_ocn_rtot_M
              END DO
           end if
           IF (ctrl_data_save_sig_ocnatm) THEN
              DO l=1,n_l_atm
                 ia = conv_iselected_ia(l)
                 int_ocnatm_sig(ia) = int_ocnatm_sig(ia) + &
                      & loc_dtyr*SUM(phys_ocnatm(ipoa_A,:,:)*dum_sfcatm1(ia,:,:))*loc_ocnatm_rtot_A
              END DO
           end if
           IF (ctrl_data_save_sig_fexport) THEN
              DO l=1,n_l_sed
                 is = conv_iselected_is(l)
                 int_fexport_sig(is) = int_fexport_sig(is) + &
                      & SUM(bio_settle(is,:,:,n_k))
              END DO
           end if
           IF (ctrl_data_save_sig_focnatm) THEN
              DO l=3,n_l_atm
                 ia = conv_iselected_ia(l)
                 int_focnatm_sig(ia) = int_focnatm_sig(ia) + &
                      & loc_dtyr*SUM(locij_focnatm(ia,:,:))
              END DO
           end if
           IF (ctrl_data_save_sig_focnsed) THEN
              DO l=1,n_l_sed
                 is = conv_iselected_is(l)
                 int_focnsed_sig(is) = int_focnsed_sig(is) + &
                      & SUM(locij_focnsed(is,:,:))
              END DO
           end if
           IF (ctrl_data_save_sig_fsedocn) THEN
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 int_fsedocn_sig(io) = int_fsedocn_sig(io) + loc_dts*&
                      & SUM(phys_ocn(ipo_A,:,:,n_k)*dum_sfxocn1(io,:,:))
              END DO
           end if
           IF (ctrl_data_save_sig_ocn_sur) THEN
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 int_ocn_sur_sig(io) = int_ocn_sur_sig(io) + loc_dtyr*&
                      & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k)*ocn(io,:,:,n_k))*loc_ocn_rtot_A
              END DO
           end if
           IF (ctrl_data_save_sig_carb_sur) THEN
              DO ic=1,n_carb
                 int_carb_sur_sig(ic) = int_carb_sur_sig(ic) + loc_dtyr*&
                      & SUM((1.0 - phys_ocnatm(ipoa_seaice,:,:))*phys_ocn(ipo_A,:,:,n_k)*carb(ic,:,:,n_k))*loc_ocn_rtot_A
              END DO
           end if
           IF (ctrl_data_save_sig_ocn_sur) THEN
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 int_ocn_ben_sig(io) = int_ocn_ben_sig(io) + loc_dtyr*&
                      & SUM(locij_mask_ben(:,:)*phys_ocn(ipo_A,:,:,n_k)*locij_ocn_ben(io,:,:))*loc_ocnsed_rtot_A_ben
              END DO
           end if
           IF (ctrl_data_save_sig_misc) THEN
              ! record GEMlite phase
              if (dum_gemlite) int_misc_gemlite_sig = int_misc_gemlite_sig + loc_dtyr
              ! sea-ice
              int_misc_seaice_sig = int_misc_seaice_sig + &
                   & loc_dtyr*SUM(phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_seaice,:,:))
              If (sum(phys_ocnatm(ipoa_seaice,:,:)) > const_real_nullsmall) then
                 int_misc_seaice_sig_th = int_misc_seaice_sig_th + &
                      & loc_dtyr* &
                      & SUM(phys_ocnatm(ipoa_seaice_th,:,:)*phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_seaice,:,:))/ &
                      & SUM(phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_seaice,:,:))
              end If
              int_misc_seaice_sig_vol = int_misc_seaice_sig_vol + &
                   & loc_dtyr*SUM(phys_ocnatm(ipoa_seaice_th,:,:)*phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_seaice,:,:))
              ! overturning streamfunction
              CALL sub_calc_psi( &
                   & phys_ocn(ipo_gu:ipo_gw,:,:,:),loc_opsi,loc_opsia,loc_opsip,loc_zpsi,loc_opsia_minmax,loc_opsip_minmax &
                   & )
              int_misc_opsi_min_sig = int_misc_opsi_min_sig + loc_dtyr*minval(loc_opsi(:,:))
              int_misc_opsi_max_sig = int_misc_opsi_max_sig + loc_dtyr*maxval(loc_opsi(:,:))
              int_misc_opsia_min_sig = int_misc_opsia_min_sig + loc_dtyr*loc_opsia_minmax(1)
              int_misc_opsia_max_sig = int_misc_opsia_max_sig + loc_dtyr*loc_opsia_minmax(2)
              ! calculate current mean surface land (air) temperature SLT (degrees C)
              loc_sig = 0.0
              loc_tot_A = 0.0
              DO i=1,n_i
                 DO j=1,n_j
                    loc_k1 = goldstein_k1(i,j)
                    IF (n_k < loc_k1) THEN
                       loc_sig = loc_sig + phys_ocnatm(ipoa_A,i,j)*dum_sfcatm1(ia_T,i,j)
                       loc_tot_A = loc_tot_A + phys_ocnatm(ipoa_A,i,j)
                    end IF
                 end DO
              end DO
              if (loc_tot_A > const_real_nullsmall) then
                 int_misc_SLT_sig = int_misc_SLT_sig + loc_dtyr*loc_sig/loc_tot_A
              else
                 int_misc_SLT_sig = 0.0
              end if
              ! solar insolation (and orbitally-related information)
              ! NOTE: apply ocean mask (@ surface)
              ! (1) mean global properties
              int_misc_ocn_solfor_sig = int_misc_ocn_solfor_sig + &
                   & loc_dtyr*loc_ocn_rtot_A*sum(phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_solfor,:,:))
              int_misc_ocn_fxsw_sig = int_misc_ocn_fxsw_sig + &
                   & loc_dtyr*loc_ocn_rtot_A*sum(phys_ocn(ipo_A,:,:,n_k)*phys_ocnatm(ipoa_fxsw,:,:))
              ! (2) latitudinal/seasonal properties
              !     NOTE: done very crudely and taking values from a single specified time-step of the averaging only
              if (int_t_sig_count == par_t_sig_count) then
                 snap_misc_ocn_solfor_N_sig = &
                      & sum(phys_ocnatm(ipoa_A,:,par_sig_j_N)*phys_ocnatm(ipoa_solfor,:,par_sig_j_N)) / &
                      & sum(phys_ocnatm(ipoa_A,:,par_sig_j_N))
                 snap_misc_ocn_solfor_S_sig = &
                      & sum(phys_ocnatm(ipoa_A,:,par_sig_j_S)*phys_ocnatm(ipoa_solfor,:,par_sig_j_S)) / &
                      & sum(phys_ocnatm(ipoa_A,:,par_sig_j_S))
              end if
           end if
           ! calcualte mean global sediement properties
           ! NOTE: for isotopic properties, set thresholds for including sediments in the global mean
           !       (to avoid possible NaNs or non-sensicle values associated with vanishigly small bulk abundances):
           !        0.1% for wt% POM
           !       10.0% for wt% CaCO3
           !        1.0% for wt% opal
           !        0.1% for wt% det
           IF (ctrl_data_save_sig_ocnsed) THEN
              DO l=1,n_l_sed
                 is = conv_iselected_is(l)
                 SELECT CASE (sed_type(is))
                 CASE (par_sed_type_bio,par_sed_type_abio, &
                      & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
                      & par_sed_type_scavenged)
                    loc_sig = SUM(phys_ocn(ipo_A,:,:,n_k)*dum_sfcsed1(is,:,:))
                    loc_tot_A = loc_ocnsed_tot_A
                 case (par_sed_type_age,11:20)
                    loc_sig = 0.0
                    loc_tot_A = 0.0
                    DO i=1,n_i
                       DO j=1,n_j
                          if ((sed_type(sed_dep(is)) == par_sed_type_POM) .OR. (sed_dep(is) == is_POC)) then
                             If (dum_sfcsed1(is_POC,i,j) > 0.1) then
                                loc_sig = loc_sig + phys_ocn(ipo_A,i,j,n_k)*dum_sfcsed1(sed_dep(is),i,j)*dum_sfcsed1(is,i,j)
                                loc_tot_A = loc_tot_A + phys_ocn(ipo_A,i,j,n_k)*dum_sfcsed1(sed_dep(is),i,j)
                             end if
                          end if
                          if ((sed_type(sed_dep(is)) == par_sed_type_CaCO3) .OR. (sed_dep(is) == is_CaCO3)) then
                             If (dum_sfcsed1(is_CaCO3,i,j) > 10.0) then
                                loc_sig = loc_sig + phys_ocn(ipo_A,i,j,n_k)*dum_sfcsed1(sed_dep(is),i,j)*dum_sfcsed1(is,i,j)
                                loc_tot_A = loc_tot_A + phys_ocn(ipo_A,i,j,n_k)*dum_sfcsed1(sed_dep(is),i,j)
                             end if
                          end if
                          if ((sed_type(sed_dep(is)) == par_sed_type_opal) .OR. (sed_dep(is) == is_opal)) then
                             If (dum_sfcsed1(is_opal,i,j) > 1.0) then
                                loc_sig = loc_sig + phys_ocn(ipo_A,i,j,n_k)*dum_sfcsed1(sed_dep(is),i,j)*dum_sfcsed1(is,i,j)
                                loc_tot_A = loc_tot_A + phys_ocn(ipo_A,i,j,n_k)*dum_sfcsed1(sed_dep(is),i,j)
                             end if
                          end if
                          if ((sed_type(sed_dep(is)) == par_sed_type_det) .OR. (sed_dep(is) == is_det)) then
                             If (dum_sfcsed1(is_det,i,j) > 0.1) then
                                loc_sig = loc_sig + phys_ocn(ipo_A,i,j,n_k)*dum_sfcsed1(sed_dep(is),i,j)*dum_sfcsed1(is,i,j)
                                loc_tot_A = loc_tot_A + phys_ocn(ipo_A,i,j,n_k)*dum_sfcsed1(sed_dep(is),i,j)
                             end if
                          end if
                       end DO
                    end DO
                 case default
                    loc_sig = 0.0
                    loc_tot_A = 0.0
                 end SELECT
                 if (loc_tot_A > const_real_nullsmall) then
                    int_ocnsed_sig(is) = int_ocnsed_sig(is) + loc_dtyr*loc_sig/loc_tot_A
                 else
                    int_ocnsed_sig(is) = 0.0
                 end if
              END DO
           end if
           ! aeolian Fe input and solubility etc
           int_misc_det_Fe_tot_sig = int_misc_det_Fe_tot_sig + SUM(phys_ocnatm(ipoa_totFe,:,:))
           int_misc_det_Fe_dis_sig = int_misc_det_Fe_dis_sig + SUM(phys_ocnatm(ipoa_solFe,:,:)*phys_ocnatm(ipoa_totFe,:,:))
           ! diagnostic diagnostics
           ! NOTE: the area used in calculating a mean global biological uptake rates must be total ocean surface area
           !       (i.e., including sea-ice covered area)
           loc_tot_A = sum(phys_ocn(ipo_A,:,:,n_k))
           IF (ctrl_data_save_sig_diag) THEN
              DO ib=1,n_diag_bio
                 int_diag_bio_sig(ib) = int_diag_bio_sig(ib) + &
                      & SUM(phys_ocn(ipo_A,:,:,n_k)*diag_bio(ib,:,:))/loc_tot_A
              END DO
              DO id=1,n_diag_geochem
                 int_diag_geochem_sig(id) = int_diag_geochem_sig(id) + &
                      & SUM(phys_ocn(ipo_M,:,:,:)*diag_geochem(id,:,:,:))*loc_ocn_rtot_M
              END DO
              DO l=1,n_l_ocn
                 io = conv_iselected_io(l)
                 if (dum_gemlite) then
                    int_diag_weather_sig(io) = int_diag_weather_sig(io) + loc_dtyr*SUM(dum_sfxsumrok1(io,:,:))
                 else
                    int_diag_weather_sig(io) = int_diag_weather_sig(io) + SUM(dum_sfxsumrok1(io,:,:))
                 end if
              END DO
              DO l=1,n_l_atm
                 ia = conv_iselected_ia(l)
                 int_diag_airsea_sig(ia) = int_diag_airsea_sig(ia) + loc_dtyr*SUM(diag_airsea(ia,:,:))
              END DO
           end if
           ! misc - 2D diagnostics
           IF (force_restore_atm_select(ia_pCO2_13C) .AND. force_flux_atm_select(ia_pCO2_13C)) THEN
              int_diag_misc_2D_sig(idiag_misc_2D_FpCO2) = &
                   & int_diag_misc_2D_sig(idiag_misc_2D_FpCO2) + loc_dtyr*SUM(diag_misc_2D(idiag_misc_2D_FpCO2,:,:))
              int_diag_misc_2D_sig(idiag_misc_2D_FpCO2_13C) = &
                   & int_diag_misc_2D_sig(idiag_misc_2D_FpCO2_13C) + loc_dtyr*SUM(diag_misc_2D(idiag_misc_2D_FpCO2_13C,:,:))
           end IF
           ! ### ADD ADDITIONAL TIME-SERIES UPDATES HERE ######################################################################### !
           !
           ! ##################################################################################################################### !

        end if if_save3

        ! update signal-averaging ingetrated time
        int_t_sig = int_t_sig + loc_dtyr
        int_t_sig_count = int_t_sig_count + 1

     end IF if_save2

     ! set save time
	 ! 28/12/2011: changed 'rounding' equation because integer time > 2 Myr ended up over the precision limit
	 !             (-2147483.648 ended up being reported)
	 !             => instead simply 'bumped' the year value to avoid the occurrence of e.g. x.999999 reported values
     if ((par_data_save_sig_dt - int_t_sig) < conv_s_yr) then
        if (ctrl_misc_t_BP) then
           loc_yr_save = loc_t + int_t_sig/2.0 - par_misc_t_end
        else
           loc_yr_save = par_misc_t_end - loc_t - int_t_sig/2.0
        end if
           ! clean up year
           ! NOTE: original 'rounding' equation resulted in integer time > 2 Myr ended up over the precision limit
           !             (-2147483.648 ended up being reported)
           loc_yr_save = &
                & real(int(loc_yr_save)) + real(int(1000.0*(loc_yr_save - real(int(loc_yr_save)) + 0.0005)))/1000.00

        if (dum_save) then

           WRITE(unit=6,fmt='(A57,f12.3)') &
               & '>>> SAVING BIOGEM TIME-SERIES AVERAGE CENTERED @ year : ', &
               & loc_yr_save
           ! check that there is no chance of dividing-by-zero ...
           IF (int_t_sig > const_real_nullsmall) then
              IF (ctrl_data_save_sig_ascii) then
                 CALL sub_data_save_runtime(loc_yr_save)
              else
                 CALL sub_save_netcdf_runtime(loc_yr_save)
              end IF
              ! ### OPTIONAL CODE ################################################################################################ !
              ! NOTE: netCDF time-series saving is disabled by default, partly because updating has lagged behind the ASCII version
              !       (and partly because personally, I never used the data files ...)
              !CALL sub_save_netcdf_runtime(loc_yr_save)
              ! ################################################################################################################## !
              if (.NOT. ctrl_debug_lvl0) then
                 IF (ctrl_misc_t_BP) THEN
                    loc_yr = loc_t + par_misc_t_end
                 ELSE
                    loc_yr = par_misc_t_end - loc_t
                 END IF
                 loc_opsi_scale = goldstein_dsc*goldstein_usc*const_rEarth*1.0E-6
                 CALL sub_calc_psi( &
                      & phys_ocn(ipo_gu:ipo_gw,:,:,:),loc_opsi,loc_opsia,loc_opsip,loc_zpsi,loc_opsia_minmax,loc_opsip_minmax &
                      & )
                 call sub_echo_runtime(loc_yr ,loc_opsi_scale,loc_opsia_minmax,dum_sfcatm1(:,:,:),dum_gemlite)
              endif
           end if

        else
           PRINT*,'>>> *SKIPPING* SAVING BIOGEM TIME-SERIES AVERAGE CENTERED @ year : ',loc_yr_save
        end if

        par_data_save_sig_i = par_data_save_sig_i - 1
        ! reset array values
        CALL sub_init_int_timeseries()
     END IF
  end IF if_save1

end SUBROUTINE diag_biogem_timeseries
! ******************************************************************************************************************************** !
