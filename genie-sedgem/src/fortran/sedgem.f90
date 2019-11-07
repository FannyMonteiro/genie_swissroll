
! ******************************************************************************************************************************** !
! TIME-STEP SEDGEM
SUBROUTINE sedgem(          &
     & dum_dts,             &
     & dum_sfxsumsed,       &
     & dum_sfcsumocn,       &
     & dum_sfcsed,          &
     & dum_sfxocn,          &
     & dum_reinit_sfxsumsed &
     & )
  USE sedgem_lib
  USE sedgem_box
  USE sedgem_data
  IMPLICIT NONE
  ! dummy arguments
  REAL,INTENT(in)::dum_dts                                              ! time-step
  real,DIMENSION(n_sed,n_i,n_j),intent(inout)::dum_sfxsumsed            ! sediment rain flux interface array
  real,DIMENSION(n_ocn,n_i,n_j),intent(in)::dum_sfcsumocn               ! ocean composition interface array
  real,DIMENSION(n_sed,n_i,n_j),intent(out)::dum_sfcsed                 ! sediment composition interface array
  real,DIMENSION(n_ocn,n_i,n_j),intent(out)::dum_sfxocn                 ! sediment dissolution flux interface array
  logical::dum_reinit_sfxsumsed                                         ! reinitialize sedimentation array?
  ! local variables
  integer::i,j,l,io,is                                         ! COUNTING AND ARRAY INDEX VARIABLES
  integer::loc_i,loc_tot_i                                     ! array index conversion variables
  real::loc_dtyr                                               ! local time step (in years)
  real::loc_dts                                                ! local time step (in seconds)
  real::loc_tot_A                                              ! local total area
  logical::loc_flag_save                                       ! local flag
  REAL,DIMENSION(n_sed)::loc_fracdecay_sed                     ! local reduction factor for decaying sediment tracers
  real,DIMENSION(n_ocn)::loc_fhydrothermal                     ! local dissolved tracer array for hydrothermal input
  real,DIMENSION(n_ocn)::loc_flowTalteration                   ! local dissolved tracer array for low T alteration sink
  real,DIMENSION(n_i,n_j)::loc_phys_sed_mask_deepsea           ! 
  real::loc_tot,loc_standard                                   ! 
  real::loc_r7Li                                               ! 
  real::loc_alpha,loc_R,loc_delta                              ! 
  real::loc_fsed                                               ! 
  real,DIMENSION(n_sed,n_i,n_j)::loc_sfxsumsed_OLD                      ! sediment rain flux interface array (COPY)

  ! *** STORE PREVIOUS ITERATION DATA ***
  sed_fsed_OLD(:,:,:) = sed_fsed(:,:,:) 
  sed_fdis_OLD(:,:,:) = sed_fdis(:,:,:)
  ! copy current (passed) sediemnt flux
  loc_sfxsumsed_OLD(:,:,:) = dum_sfxsumsed(:,:,:)

  ! *** INITIALIZE RESULTS ARRAYS ***
  dum_sfxocn(:,:,:)  = 0.0     ! 
  sed_fdis(:,:,:)    = 0.0     ! 
  sedocn_fnet(:,:,:) = 0.0     !

  ! *** INITIALIZE LOCAL ARRAYS ***
  loc_fhydrothermal(:)   = 0.0
  loc_flowTalteration(:) = 0.0
  loc_phys_sed_mask_deepsea(:,:) = 0.0

  ! *** CALCULATE SEDGEM TIME ***
  IF (ctrl_misc_debug4) print*,'*** CALCULATE SEDGEM TIME ***'
  ! calculate sediment model time step length
  ! NOTE: convert between time units of BioGeM (years) and GOLDSTEIn (use <tsc> scale factor to convert to seconds)
  loc_dts  = dum_dts
  loc_dtyr = loc_dts/conv_yr_s

  ! *** DECAY RADIOACTIVE TRACERS ***
  ! calculate fractional reduction factors for decaying isotopes
  loc_fracdecay_sed(:) = EXP(-loc_dtyr*const_lambda_sed(:))
  ! decay radioactive tracers
  DO l=1,n_l_sed
     is = conv_iselected_is(l)
     IF (abs(const_lambda_sed(is)).gt.const_real_nullsmall) THEN
        sed_top(is,:,:) = loc_fracdecay_sed(is)*sed_top(is,:,:)
        sed(is,:,:,:)   = loc_fracdecay_sed(is)*sed(is,:,:,:)
     END if
  end do

  ! *** UPDATE MASKS ***
  IF (ctrl_misc_debug4) print*,'*** UPDATE MASKS ***'
  DO i=1,n_i
     DO j=1,n_j
        ! catch a zero salinity as an indication that there is no valid corresponding ocean grid point
        ! => permanently amend sediment mask
        IF (sed_mask(i,j) .AND. (dum_sfcsumocn(io_S,i,j) < const_real_nullsmall)) then
           sed_mask(i,j)      = .FALSE.
           sed_mask_reef(i,j) = .FALSE.
           sed_mask_muds(i,j) = .FALSE.
           phys_sed(ips_mask_sed,i,j)      = 0.0
           phys_sed(ips_mask_sed_reef,i,j) = 0.0
           phys_sed(ips_mask_sed_muds,i,j) = 0.0
           sed_save_mask(i,j) = .FALSE.
        end IF
     end DO
  end DO
  ! set local deep-sea mask & calculate total area
  loc_phys_sed_mask_deepsea(:,:) = phys_sed(ips_mask_sed,:,:) - phys_sed(ips_mask_sed_reef,:,:) - phys_sed(ips_mask_sed_muds,:,:)
  loc_tot_A = sum(loc_phys_sed_mask_deepsea(:,:)*phys_sed(ips_A,:,:))

  ! *** UPDATE CARBONATE CHEMSITRY ***
  DO i=1,n_i
     DO j=1,n_j
        IF (sed_mask(i,j)) then
           call sub_calc_carbconst(        &
                & phys_sed(ips_D,i,j),     &
                & dum_sfcsumocn(io_T,i,j), &
                & dum_sfcsumocn(io_S,i,j), &
                & sed_carbconst(:,i,j)     &
                & )
           if (ocn_select(io_Ca) .AND. ocn_select(io_Mg)) then
              call sub_adj_carbconst(          &
                   & dum_sfcsumocn(io_Ca,i,j), &
                   & dum_sfcsumocn(io_Mg,i,j), &
                   & sed_carbconst(:,i,j)      &
                   & )
           end if
           call sub_calc_carb(                &
                & dum_sfcsumocn(io_DIC,i,j),  &
                & dum_sfcsumocn(io_ALK,i,j),  &
                & dum_sfcsumocn(io_Ca,i,j),   &
                & dum_sfcsumocn(io_PO4,i,j),  &
                & dum_sfcsumocn(io_SiO2,i,j), &
                & dum_sfcsumocn(io_B,i,j),    &
                & dum_sfcsumocn(io_SO4,i,j),  &
                & dum_sfcsumocn(io_F,i,j),    &
                & dum_sfcsumocn(io_H2S,i,j),  &
                & dum_sfcsumocn(io_NH4,i,j),  &
                & sed_carbconst(:,i,j),       &
                & sed_carb(:,i,j),            &
                & sed_carbalk(:,i,j)          &
                & )
           ! re-calculate carbonate system isotopic properties
           if (ocn_select(io_DIC_13C)) then
              call sub_calc_carb_r13C(              &
                   & dum_sfcsumocn(io_T,i,j),       &
                   & dum_sfcsumocn(io_DIC,i,j),     &
                   & dum_sfcsumocn(io_DIC_13C,i,j), &
                   & sed_carb(:,i,j),               &
                   & sed_carbisor(:,i,j)            &
                   & )
           end IF
           if (ocn_select(io_DIC_14C)) then
              call sub_calc_carb_r14C(              &
                   & dum_sfcsumocn(io_T,i,j),       &
                   & dum_sfcsumocn(io_DIC,i,j),     &
                   & dum_sfcsumocn(io_DIC_14C,i,j), &
                   & sed_carb(:,i,j),               &
                   & sed_carbisor(:,i,j)            &
                   & )
           end IF
        end IF
     end DO
  end DO

  ! *** EARLY DIAGENESIS PROCESSES ***
  ! NOTE: <dum_sfxsumsed> in units of (mol m-2 per time-step)
  ! NOTE: dum_sfxocn(io,:,:) in units of (mol m-2 s-1)
  DO i=1,n_i
     DO j=1,n_j
        IF (sed_mask(i,j)) THEN
           ! ammend sediment rain flux according to prescribed detrital input
           ! NOTE: convert units from (g cm-2 kyr-1) to (mol m-2 (per time-step))
           if (sed_select(is_det)) then
              dum_sfxsumsed(is_det,i,j) = dum_sfxsumsed(is_det,i,j) + &
                   & conv_m2_cm2*conv_det_g_mol*(conv_yr_kyr*loc_dtyr)*par_sed_fdet
           endif
           ! add ash layer (if selected)
           if (sed_select(is_ash)) then
              if (par_sed_ashevent) then
                 dum_sfxsumsed(is_ash,i,j) = dum_sfxsumsed(is_ash,i,j) + &
                      & conv_m2_cm2*conv_det_g_mol*(conv_yr_kyr*loc_dtyr)*par_sed_ashevent_fash
              end if
           endif
           ! account for clay formation
           If (ocn_select(io_Li) .AND. ocn_select(io_Ca)) then
              loc_fsed = par_sed_clay_fLi_alpha*dum_sfxsumsed(is_det,i,j)*(dum_sfcsumocn(io_Li,i,j)/dum_sfcsumocn(io_Ca,i,j))
              dum_sfxsumsed(is_detLi,i,j) = dum_sfxsumsed(is_detLi,i,j) + loc_fsed
              dum_sfxocn(io_Li,i,j) = dum_sfxocn(io_Li,i,j) - loc_fsed/dum_dts
              if (ocn_select(io_Li_7Li)) then
                 loc_standard = const_standards(ocn_type(io_Li_7Li))
                 if (dum_sfcsumocn(io_Li,i,j) > const_real_nullsmall) then
                    loc_r7Li = dum_sfcsumocn(io_Li_7Li,i,j)/dum_sfcsumocn(io_Li,i,j)
                 else
                    loc_r7Li = 0.0
                 end if
                 loc_alpha = 1.0 + par_sed_clay_7Li_epsilon/1000.0
                 loc_R = loc_r7Li/(1.0 - loc_r7Li)
                 loc_fsed = (loc_alpha*loc_R/(1.0 + loc_alpha*loc_R))*loc_fsed
                 dum_sfxsumsed(is_detLi_7Li,i,j) = dum_sfxsumsed(is_detLi_7Li,i,j) + loc_fsed
                 dum_sfxocn(io_Li_7Li,i,j) = dum_sfxocn(io_Li_7Li,i,j) - loc_fsed/dum_dts
              end if
           end if
        end IF
     end DO
  end DO
  ! deselect ash fall
  if (sed_select(is_ash)) then
     if (par_sed_ashevent) par_sed_ashevent = .false.
  endif

  ! *** FORAM TRACERS ***
  ! 
  DO i=1,n_i
     DO j=1,n_j
        IF (sed_mask(i,j)) THEN
           ! add foram tracers
           if (sed_select(is_CaCO3) .AND. sed_select(is_foram_b_13C)) then
              ! calculate 13C/12C fractionation between DIC and CaCO3
              SELECT CASE (opt_sed_foram_b_13C_delta)
              CASE ('NONE')
                 loc_delta = 0.0
              case ('SPERO')
                 ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                 loc_delta = 0.0
                 ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              end SELECT
              loc_alpha = 1.0 + loc_delta/1000.0
              loc_R = sed_carbisor(ici_HCO3_r13C,i,j)/(1.0 - sed_carbisor(ici_HCO3_r13C,i,j))
              loc_fsed = (loc_alpha*loc_R/(1.0 + loc_alpha*loc_R))*dum_sfxsumsed(is_CaCO3,i,j)
              dum_sfxsumsed(is_foram_b_13C,i,j) = loc_fsed
           end if
        end IF
     end DO
  end DO

  ! *** UPDATE SEDIMENTS ***
  IF (ctrl_misc_debug4) print*,'*** UPDATE SEDIMENTS ***'
  DO i=1,n_i
     DO j=1,n_j
        ! calculate sediment rain flux from sediment->ocean flux (convert units)
        ! NOTE: if the sediment grid point lies outside of the sediment mask, then
        !       dissolve all sediment tracers and set the ocean tracer dissolution flux equal to this
        ! convert units of sedimentation flux
        ! NOTE: <dum_sfxsumsed> in units of (mol m-2)
        ! NOTE: <sed_fsed> in units of (mol cm-2)
        sed_fsed(:,i,j) = conv_cm2_m2*dum_sfxsumsed(:,i,j)
        IF (sed_mask(i,j)) THEN
           ! call sediment composition update
           ! NOTE: the values in both <sed_fsed> and <ocnsed_fnet> are updated by this routine
           if (sed_mask_reef(i,j)) then
              IF (ctrl_misc_debug3) print*,'> UPDATE SED: reef'
              call sub_update_sed_reef(    &
                   & loc_dtyr,             &
                   & i,j,                  &
                   & phys_sed(ips_D,i,j),  &
                   & dum_sfcsumocn(:,i,j), &
                   & loc_flag_save         &
                   & )
           elseif (sed_mask_muds(i,j)) then
              IF (ctrl_misc_debug3) print*,'> UPDATE SED: mud'
              call sub_update_sed_mud(     &
                   & i,j,                  &
                   & loc_flag_save         &
                   & )
           else
              IF (ctrl_misc_debug3) print*,'> UPDATE SED: (normal)'
              call sub_update_sed(         &
                   & loc_dtyr,             &
                   & i,j,                  &
                   & phys_sed(ips_D,i,j),  &
                   & dum_sfcsumocn(:,i,j), &
                   & loc_flag_save         &
                   & )
           end if
           ! save sedcore environmental conditions
           if (sed_save_mask(i,j) .AND. loc_flag_save) then
              call sub_sedgem_save_sedcoreenv( &
                   & loc_dtyr,                 &
                   & i,j,                      &
                   & sed_top(:,i,j),           &
                   & sed_fsed(:,i,j),          &
                   & sed_fdis(:,i,j),          &
                   & dum_sfcsumocn(:,i,j),     &
                   & sed_carb(:,i,j)           &
                   & )
           end if
        else
           ! NO SEDIMENTS HERE
           ! set dissolution flux (as sediment solids)
           sed_fdis(:,i,j) = sed_fsed(:,i,j)
           ! calculate equivalent ocean tracer flux
           DO l=1,n_l_sed
              is = conv_iselected_is(l)
              loc_tot_i = conv_sed_ocn_i(0,is)
              do loc_i=1,loc_tot_i
                 io = conv_sed_ocn_i(loc_i,is)
                 sedocn_fnet(io,i,j) = sedocn_fnet(io,i,j) + conv_sed_ocn(io,is)*sed_fsed(is,i,j)
              end do
           end DO
        end if
     end do
  end do
  
  ! *** HYDROTHERMAL / SEAFLOOR ALTERATION ***
  IF (ctrl_misc_debug4) print*,'*** HYDROTHERMAL / SEAFLOOR ALTERATION ***'
  ! calculate hydrothermal source fluxes (inputs)
  ! NOTE: loc_fhydrothermal_input(io) in units of (mol yr-1)
  ! NOTE: dum_sfxocn(io,:,:) in units of (mol m-2 s-1)
  If (ocn_select(io_Li)) then
     loc_fhydrothermal(io_Li) = par_sed_hydroip_fLi
     if (ocn_select(io_Li_7Li)) then
        loc_tot = loc_fhydrothermal(io_Li)
        loc_standard = const_standards(ocn_type(io_Li_7Li))
        loc_fhydrothermal(io_Li_7Li) = fun_calc_isotope_fraction(par_sed_hydroip_fLi_d7Li,loc_standard)*loc_tot
     end if
  end if
  ! re-scale global total and apread over *valid* grid points
  DO i=1,n_i
     DO j=1,n_j
        if (loc_phys_sed_mask_deepsea(i,j) > const_real_nullsmall) then
           DO l=1,n_l_ocn
              io = conv_iselected_io(l)
              dum_sfxocn(io,i,j) = dum_sfxocn(io,i,j) + loc_fhydrothermal(io)/loc_tot_A/conv_yr_s
           end DO
        end if
     end DO
  end DO
  ! calculate weathering/alteration fluxes (net sink)
  ! NOTE: dum_sfxocn(io,:,:) in units of (mol m-2 s-1)
  ! NOTE: the value of par_sed_lowTalt_fLi_alpha is calculated on the basis of a sink of x mol yr-1 (e.f. 1.0E10 mol yr-1)
  !       => [Ca/Li] * x / Atot / (365.25*24*3600)
  !          where Atot is the total sediment area and [Ca] and [Li] are assumed to be 0.01025 mol Ca, 26 umol Li
  If (ocn_select(io_Li) .AND. ocn_select(io_Ca)) then
     DO i=1,n_i
        DO j=1,n_j
           if (loc_phys_sed_mask_deepsea(i,j) > const_real_nullsmall) then
              loc_flowTalteration(io_Li) = par_sed_lowTalt_fLi_alpha*(dum_sfcsumocn(io_Li,i,j)/dum_sfcsumocn(io_Ca,i,j))
              dum_sfxocn(io_Li,i,j) = dum_sfxocn(io_Li,i,j) - loc_flowTalteration(io_Li)
              if (ocn_select(io_Li_7Li)) then
                 loc_standard = const_standards(ocn_type(io_Li_7Li))
                 if (dum_sfcsumocn(io_Li,i,j) > const_real_nullsmall) then
                    loc_r7Li = dum_sfcsumocn(io_Li_7Li,i,j)/dum_sfcsumocn(io_Li,i,j)
                 else
                    loc_r7Li = 0.0
                 end if
                 loc_alpha = 1.0 + par_sed_lowTalt_7Li_epsilon/1000.0
                 loc_R = loc_r7Li/(1.0 - loc_r7Li)
                 loc_flowTalteration(io_Li_7Li) = (loc_alpha*loc_R/(1.0 + loc_alpha*loc_R))*loc_flowTalteration(io_Li)
                 dum_sfxocn(io_Li_7Li,i,j) = dum_sfxocn(io_Li_7Li,i,j) - loc_flowTalteration(io_Li_7Li)
              end if
           end if
        end DO
     end DO
  end if

  ! *** UPDATE INTERFACE ***
  IF (ctrl_misc_debug4) print*,'*** UPDATE INTERFACE ***'
  ! update update sed->ocn interface
  ! NOTE: <sedocn_fnet> in units of (mol cm-2)
  ! NOTE: <dum_sfxocn> in units of (mol m-2 s-1)
  DO l=1,n_l_ocn
     io = conv_iselected_io(l)
     dum_sfxocn(io,:,:) = dum_sfxocn(io,:,:) + conv_m2_cm2*sedocn_fnet(io,:,:)/loc_dts
  end DO
  ! re-initialize the interfacing integrated sediment rain flux array
  ! NOTE: in normal operation, dum_reinit_sfxsumsed is .true., hence once sediment rain has been added to the sediments
  !       the flux array is zero-ed
  !       (during GEMlite phases, there is no BIOGEM updating of the ocean rain flux, hence need to preserve the flux value)
  if (dum_reinit_sfxsumsed) then
     dum_sfxsumsed(:,:,:) = 0.0
  else
     dum_sfxsumsed(:,:,:) = loc_sfxsumsed_OLD(:,:,:)
  end if
  ! update sediment interface composition data
  dum_sfcsed(:,:,:) = fun_sed_coretop()

  ! *** DEBUG ***
  ! print some debugging info if 'ctrl_misc_debug1' option is selected
  IF (ctrl_misc_debug1) THEN
     i = par_misc_debug_i
     j = par_misc_debug_j
     print*,''
     print*,'--- DEBUG ---'
     print*,'> SEDGEM LOCATION (i,j): ',i,j
     print*,'> TIME, TIME-STEP: ',loc_dts,loc_dtyr
     print*, &
          & phys_sed(ips_D,i,j),               &
          & dum_sfcsumocn(io_T,i,j),           &
          & dum_sfcsumocn(io_S,i,j)
     print*, &
          & dum_sfcsumocn(io_DIC,i,j),         &
          & dum_sfcsumocn(io_ALK,i,j),         &
          & 1.0E+06*sed_carb(ic_dCO3_cal,i,j)
     print*, &
          & 100.0*sed_top(is_CaCO3,i,j),       &
          & 100.0*sed_top(is_opal,i,j)
     print*, &
          & 1.0E+06*sed_fsed(is_CaCO3,i,j)/loc_dtyr, &
          & 1.0E+06*sed_fsed(is_POC,i,j)/loc_dtyr,   &
          & 1.0E+06*sed_fdis(is_CaCO3,i,j)/loc_dtyr, &
          & 1.0E+06*sed_fdis(is_POC,i,j)/loc_dtyr
     print*, &
          & sum(sed_fsed(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/loc_dtyr, &
          & sum(sed_fsed(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/loc_dtyr,   &
          & sum(sed_fdis(is_CaCO3,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/loc_dtyr, &
          & sum(sed_fdis(is_POC,:,:)*conv_m2_cm2*phys_sed(ips_A,:,:))/loc_dtyr
     print*,'---'
     print*,''
  end if
  ! catch any carbonate chemsitry errors arising in sub_update_sed
  if (error_stop) then
     call end_sedgem(     &
          & dum_dts,      &
          & dum_sfcsumocn &
          & )
  end if

  ! *** RUN-TIME OUTPUT ***
  ! GHC 20/05/09 - Save time-series output
  IF (ctrl_timeseries_output) THEN
     ! increment timestep counter  
     tstep_count = tstep_count + 1  
     ! if output due then change year  
     CALL sub_output_year()  
     IF (tstep_count.eq.output_tsteps_0d(output_counter_0d)) THEN
        call sub_data_save_seddiag_GLOBAL(loc_dtyr,dum_sfcsumocn)  
     ENDIF
     IF (tstep_count.eq.output_tsteps_2d(output_counter_2d)) THEN
        ! save requested sediment cores as ASCII     
        ! call sub_sedgem_save_sedcore()
        ! save oecan-sediment interface properties
        !if (ctrl_data_save_ascii) call sub_data_save_seddiag_2D(loc_dtyr,dum_sfcsumocn)
        call sub_save_netcdf(year)
        call sub_save_netcdf_sed2d(loc_dtyr,dum_sfcsumocn)
        call sub_closefile(ntrec_siou)
        ntrec_sout = ntrec_sout + 1  
     ENDIF
     ! if output then increment output counter  
     CALL sub_output_counters()
  ENDIF

end SUBROUTINE sedgem
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! AGE SEDIMENTS
SUBROUTINE sedgem_dsedage(  &
     & dum_dts,             &
     & dum_sfxsumsed        &
     & )
  USE sedgem_lib
  IMPLICIT NONE
  ! dummy arguments
  REAL,INTENT(in)::dum_dts                                              ! time-step
  real,DIMENSION(n_sed,n_i,n_j),intent(inout)::dum_sfxsumsed            ! sediment rain flux interface array
  ! local variables
  integer::i,j                                                 ! COUNTING AND ARRAY INDEX VARIABLES
  real::loc_dtyr                                               ! 
  ! set age decrement (years)
  loc_dtyr = dum_dts/conv_yr_s
  ! decrement (CaCO3) sediment age tracer
  DO i=1,n_i
     DO j=1,n_j
        IF (sed_mask(i,j)) THEN
           dum_sfxsumsed(is_CaCO3_age,i,j) = dum_sfxsumsed(is_CaCO3_age,i,j) - loc_dtyr*dum_sfxsumsed(is_CaCO3,i,j)
        end IF
     end DO
  end DO
end SUBROUTINE sedgem_dsedage
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! RESTART SEDGEM (RESTART DATA DUMP)
SUBROUTINE rest_sedgem()
  USE sedgem_lib
  USE genie_util, ONLY: check_unit, check_iostat
  IMPLICIT NONE
  ! local variables
  integer::l
  integer::ios  ! for file open, close, read & write checks
  CHARACTER(len=255)::loc_filename
  ! initialize local variables
  loc_filename = TRIM(par_outdir_name)//TRIM(par_outfile_name)
  ! dump data
  call check_unit(out,__LINE__,__FILE__)
  OPEN(unit=out,status="replace",file=loc_filename,form="unformatted",action="write",iostat=ios)
  call check_iostat(ios,__LINE__,__FILE__)
  WRITE(unit=out,iostat=ios)                              &
       & n_l_sed,                                         &
       & (conv_iselected_is(l),l=1,n_l_sed),              &
       & (sed(conv_iselected_is(l),:,:,:),l=1,n_l_sed),   &
       & (sed_top(conv_iselected_is(l),:,:),l=1,n_l_sed), &
       & sed_top_h(:,:)
  call check_iostat(ios,__LINE__,__FILE__)
  close(unit=out,iostat=ios)
  call check_iostat(ios,__LINE__,__FILE__)
end SUBROUTINE rest_sedgem
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! ENKF DATA DUMP
SUBROUTINE enkf_sedgem_datadump()
  USE sedgem_lib
  USE sedgem_box
  USE genie_util, ONLY: check_unit, check_iostat
  IMPLICIT NONE
  ! local variables
  CHARACTER(len=255)::loc_filename
  REAL,DIMENSION(n_sed,n_i,n_j)::loc_sed_coretop
  integer::ios  ! for file open, close, read & write checks
  ! calculate core-top sediment composition data
  loc_sed_coretop(:,:,:) = fun_sed_coretop()
  ! dump data
  ! NOTE: data is saved unformatted for minimal file size
  !       also means that arrays can be written directly to file without needing to loop thought data
  loc_filename = TRIM(par_outdir_name)//trim(par_outfile_name)//'.EnKF'
  call check_unit(out,__LINE__,__FILE__)
  OPEN(unit=out,status='replace',file=loc_filename,form='unformatted',action='write',iostat=ios)
  call check_iostat(ios,__LINE__,__FILE__)
  WRITE(unit=out,iostat=ios) loc_sed_coretop(is_CaCO3,:,:)
  call check_iostat(ios,__LINE__,__FILE__)
  close(unit=out,iostat=ios)
  call check_iostat(ios,__LINE__,__FILE__)
end SUBROUTINE enkf_sedgem_datadump
! ******************************************************************************************************************************** !

