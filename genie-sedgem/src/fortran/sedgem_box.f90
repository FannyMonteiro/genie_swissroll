! ******************************************************************************************************************************** !
! sedgem_box.f90
! Integral sediment system routines
! ******************************************************************************************************************************** !


MODULE sedgem_box


  use genie_control
  USE gem_carbchem
  USE sedgem_lib
  USE sedgem_nnutils
  use sedgem_box_archer1991_sedflx
  use sedgem_box_ridgwelletal2003_sedflx
  IMPLICIT NONE
  SAVE


CONTAINS


  ! ****************************************************************************************************************************** !
  ! UPDATE SEDIMENT COMPOSITION
  ! this scheme is called form within a loop of all sediment grid point, with the following sequence of events;
  ! (a) initialize variables and calculate local constants
  ! (b) calculate new sedimenting material to be added to the sediment top ('well mixed') layer, and
  !     calculate the thickness of material that this represents
  ! (c) estimate dissolution from sediment top layer
  ! (d) update sediment top and sediment stack:
  !     - add new sedimenting material to sediment top layer
  !     - remove remineralized material
  !     - calculate potential thickness of sediment top layer including new sedimenting material
  !     - remove material to the sediment stack below if the thickness is >= par_sed_top_th cm, or
  !       add material from the sediment stack below if the thickness is < par_sed_top_th cm
  !     - update sediment stack height
  ! (e) mix the sediment stack if this is required
  ! (f) check the thickness of sediment stack:
  !     if it is within 1.0 cm of the maximum thickness, then remove the bottom n_sedcor_tot cm 
  !     entirely, and move the remaining sedimentary layers down to start at the bottom of the stack
  ! (g) calculate sedimment fluxes to ocean
  ! (h) sort out consumption of O2/NO3/SO4 resulting from organic matter remineralization
  ! NOTE: the surface ('top') sediment layer has a thickness of par_sed_top_th (cm) and
  !       has porosity of par_sed_poros_top (cm3 cm-3) (i.e., the solids actually occupy 1.0 - par_sed_poros_top cm3 cm-3)
  SUBROUTINE sub_update_sed( &
       & dum_dtyr,           &
       & dum_i,dum_j,        &
       & dum_D,              &
       & dum_sfcsumocn,      &
       & dum_flag_save       &
       & )
    IMPLICIT NONE
    ! dummy arguments
    REAL,INTENT(in)::dum_dtyr                                  ! time-step
    integer,INTENT(in)::dum_i,dum_j                            ! grid point (i,j)
    REAL,INTENT(in)::dum_D                                     ! depth
    real,DIMENSION(n_ocn),intent(in)::dum_sfcsumocn            ! ocean composition interface array
    logical,intent(inout)::dum_flag_save                       ! save sedcore environmental conditions?
    ! local variables
    INTEGER::l,is,io                                           ! tracer index counters
    integer::loc_i,loc_tot_i                                   ! array index conversion variables
    INTEGER::loc_n_sed_stack_top                               ! sediment stack top (incomplete) layer number
    REAL::loc_new_sed_vol                                      ! new sediment volume (as SOLID material)
    REAL::loc_dis_sed_vol                                      ! dis sediment volume (as SOLID material)
    REAL::loc_sed_top_dth                                      ! potential change in sediment surface ('top' layer) thickness
    REAL::loc_exe_sed_th                                       ! exchanged sediment thickness (w.r.t. surface sediment porosity)
    REAL::loc_sed_stack_top_th                                 ! sediment stack top thickness (i.e., of the incomplete sub-layer)
    real::loc_sed_dis_frac                                     ! (organic matter) fraction which is remineralized (<-> dissolution)
    real::loc_sed_dis_frac_P                                   ! (organic matter P) fraction which is remineralized (<-> dissolution)
    real::loc_sed_diagen_fracC                                 ! fraction of organic matter available for (CaCO3) diagenesis
    real::loc_sed_poros                                        ! 
    real::loc_sed_poros_top                                    ! 
    real::loc_r_sed_por                                        ! thickness ratio due to porosity differences (sediment stack / surface layer)
    real::loc_frac_CaCO3                                       ! 
    real::loc_frac_CaCO3_top                                   ! 
    REAL,DIMENSION(n_sed)::loc_new_sed                         ! new (sedimenting) top layer material
    REAL,DIMENSION(n_sed)::loc_dis_sed                         ! remineralized top layer material
    REAL,DIMENSION(n_sed)::loc_exe_sed                         ! top layer material to be exchanged with stack
    logical::loc_flag_stackgrow,loc_flag_stackshrink           ! growing or shrinking of stack occurs (by a 1 cm layer)

    ! *** (a) initialize variables
    IF (ctrl_misc_debug3) print*,'(a) initialize variables'
    ! zero local sediment arrays
    loc_new_sed(:) = 0.0
    loc_dis_sed(:) = 0.0
    loc_exe_sed(:) = 0.0
    ! initialize relevant location in global sediment dissolution results array
    sed_fdis(:,dum_i,dum_j) = 0.0
    ! initialize flags recording growing or shrinking of sediment stack (i.e., new layers being added or removed, respectively)
    loc_flag_stackgrow   = .FALSE.
    loc_flag_stackshrink = .FALSE.

    IF (ctrl_misc_debug3) print*,'(b) calculate Archer et al. [2002] bioturbation intensity profile (if selected)'
    ! *** (b) calculate bioturbation intensity
    IF (ctrl_sed_bioturb) then
       if (ctrl_sed_bioturb_Archer) then
          call sub_calc_sed_mix(               &
               & dum_dtyr,                     &
               & sed_fsed(is_POC,dum_i,dum_j), &
               & dum_sfcsumocn(io_O2),         &
               & par_sed_mix_k(:)              &
               & )
       end if
       phys_sed(ips_mix_k0,dum_i,dum_j) = par_sed_mix_k(par_n_sed_mix)
    else
       phys_sed(ips_mix_k0,dum_i,dum_j) = par_sed_mix_k_sur_min
    end if

    IF (ctrl_misc_debug3) print*,'(c) calculate new sedimenting material to be added to the sediment top layer'
    ! *** (c) calculate new sedimenting material to be added to the sediment top layer
    !         NOTE: sedimentary material is represented as SOILDS (i.e., zero porosity)
    !         NOTE: convert unts if new particulate matter added to sediments; from (mol cm-2) to (cm3 cm-2)
    ! calculate rain input of new solid material
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_new_sed(is) = conv_sed_mol_cm3(is)*sed_fsed(is,dum_i,dum_j)
    end do
    ! calculate volume of added material (as SOILD matter. i.e., zero porosity), in units of cm3 (cm-2)
    ! NOTE: nutrients associated with organic carbon (POP, PON, POFe) have only a 'virtual volume' and so are not counted
    ! NOTE: ditto for isotope tracers
    ! NOTE: check that rain thickness of sedimentating material is >= 0.0 cm
    loc_new_sed_vol = fun_calc_sed_vol(loc_new_sed(:))
    IF (loc_new_sed_vol < -const_real_nullsmall) THEN
       CALL sub_report_error( &
            & 'sedgem_box','sub_update_sed','sediment input < 0.0 cm3', &
            & 'STOPPING', &
            & (/real(dum_i), real(dum_j), loc_new_sed_vol, &
            & loc_new_sed(is_POC),   &
            & loc_new_sed(is_CaCO3), &
            & loc_new_sed(is_opal),  &
            & loc_new_sed(is_det)    &
            & /),.TRUE. &
            & )
    END IF

    IF (ctrl_misc_debug3) print*,'(d) estimate dissolution/remineralization from sediment top (mixed) sedimentary layer'
    ! *** (d) estimate dissolution/remineralization from sediment top ('well mixed') sedimentary layer
    !         NOTE: material is represented and stored as cm3 of SOILDs (i.e., zero porosity) in the sediment layers

    IF (ctrl_misc_debug4) print*,'*** diagenesis - organic matter remineralization ***'
    ! *** diagenesis - organic matter remineralization ***
    !     NOTE: particulate fluxes are in units of (mol cm-2 yr-1)
    select case (par_sed_diagen_Corgopt)
    case ('archer2002muds')
       ! MUDS
       ! ### <INSERT CODE> ####################################################################################################### !
       ! 
       ! ######################################################################################################################### !
    case default
       ! the only alternative to a wonderful diagenetic scheme with all possible redox whistles and bells that I can cook up,
       ! is a VERY crude way of estimating preservation of organic matter depending in the occurrence of anoxic/euxinic conditions
       If (ocn_select(io_O2) .AND. ocn_select(io_H2S)) then
          If (dum_sfcsumocn(io_H2S) > const_real_nullsmall) then
             loc_sed_dis_frac     = 1.0 - par_sed_diagen_fracCpres_eux
             loc_sed_dis_frac_P   = 1.0 - par_sed_diagen_fracC2Ppres_eux*par_sed_diagen_fracCpres_eux
             loc_sed_diagen_fracC = 1.0 - par_sed_diagen_fracCpres_eux
          elseif (dum_sfcsumocn(io_O2) < const_real_nullsmall) then
             loc_sed_dis_frac     = 1.0 - par_sed_diagen_fracCpres_anox
             loc_sed_dis_frac_P   = 1.0 - par_sed_diagen_fracC2Ppres_anox*par_sed_diagen_fracCpres_anox
             loc_sed_diagen_fracC = 1.0 - par_sed_diagen_fracCpres_anox
          else
             loc_sed_dis_frac     = 1.0 - par_sed_diagen_fracCpres_ox
             loc_sed_dis_frac_P   = 1.0 - par_sed_diagen_fracC2Ppres_ox*par_sed_diagen_fracCpres_ox
             loc_sed_diagen_fracC = 1.0 - par_sed_diagen_fracCpres_ox
          end If
       else
          loc_sed_dis_frac     = 1.0
          loc_sed_dis_frac_P   = 1.0
          loc_sed_diagen_fracC = 1.0
       end if
       ! calculate the return rain flux back to ocean
       ! NOTE: apply estimated fractional preservation
       ! NOTE: particle-reactive elements (e.g., 231Pa) remain in the sediments
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ((sed_dep(is) == is_POC) .OR. (sed_type(is) == par_sed_type_POM) .OR. (sed_dep(is) == is_PON)) then
             if (sed_type(is) == par_sed_type_scavenged) then
                loc_dis_sed(is) = 0.0
                ! deal with how particle-reactive elements are left in the sediments (i.e., what do they stick on?) ...
                ! ### <INSERT CODE> ############################################################################################## !
                ! 
                ! ################################################################################################################ !
             else
                select case (is)
                case (is_POP)
                   loc_dis_sed(is) = loc_sed_dis_frac_P*loc_new_sed(is)
                case default
                   loc_dis_sed(is) = loc_sed_dis_frac*loc_new_sed(is)
                end select
             end if
          end if
       end DO
    end select
    ! error-catching of negative dissoluiton: return rain flux back to ocean
    If (loc_dis_sed(is_POC) < -const_real_nullsmall) then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ((sed_dep(is) == is_POC) .OR. (sed_type(is) == par_sed_type_POM) .OR. (sed_dep(is) == is_PON)) then
             loc_dis_sed(is) = loc_new_sed(is)
          end if
       end DO
       CALL sub_report_error( &
            & 'sedgem_box','sub_update_sed','loc_dis_sed(is_POC) < 0.0 cm3', &
            & 'CONTINUING', &
            & (/real(dum_i), real(dum_j), loc_dis_sed(is_POC) &
            & /),.FALSE. &
            & )
    end if

    IF (ctrl_misc_debug4) print*,'*** diagenesis - CaCO3 dissolution ***'
    ! *** diagenesis - CaCO3 dissolution ***
    select case (par_sed_diagen_CaCO3opt)
    case (                       &
         & 'archer1991explicit', &
         & 'ridgwell2001lookup', &
         & 'ridgwell2001nn',     &
         & 'archer2002nn'        &
         & )
       CALL sub_calc_sed_dis_CaCO3(                                                  &
            & dum_dtyr,                                                              &
            & dum_D,sed_carb(ic_dCO3_cal,dum_i,dum_j),loc_sed_diagen_fracC,          &
            & dum_sfcsumocn(:),sed_carbconst(:,dum_i,dum_j),sed_carb(:,dum_i,dum_j), &
            & loc_dis_sed(:),loc_new_sed(:),sed_top(:,dum_i,dum_j),                  &
            & phys_sed(ips_mix_k0,dum_i,dum_j)                                       &
            & )
       if (error_Archer .AND. ctrl_misc_report_err) then
          CALL sub_report_error(                                                                                           &
               & 'sedgem_box','sub_update_sed','Failure of Archer [1991] sediment scheme calculation (singular matrix)',   &
               & 'CONTINUING',                                                                                             &
               & (/real(dum_i),real(dum_j),dum_D,                                                                          &
               & dum_sfcsumocn(io_T),dum_sfcsumocn(io_S),dum_sfcsumocn(io_DIC),dum_sfcsumocn(io_ALK),dum_sfcsumocn(io_O2), &
               & sed_carb(ic_dCO3_cal,dum_i,dum_j),                                                                        &
               & sed_top(is_CaCO3,dum_i,dum_j),sed_top(is_opal,dum_i,dum_j),sed_top(is_det,dum_i,dum_j),                   &
               & loc_new_sed(is_CaCO3),loc_new_sed(is_POC),loc_new_sed(is_opal),loc_new_sed(is_det),                       &
               & loc_dis_sed(is_CaCO3)                                                                                     &
               & /),.false.                                                                                                &
               & )
       end if
       error_Archer = .FALSE.
    case default
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ((sed_dep(is) == is_CaCO3) .OR. (sed_type(is) == par_sed_type_CaCO3)) then
             if (sed_type(is) == par_sed_type_scavenged) then
                loc_dis_sed(is) = 0.0
                ! deal with how particle-reactive elements are left in the sediments (i.e., what do they stick on?) ...
                ! ### <INSERT CODE> ############################################################################################## !
                ! 
                ! ################################################################################################################ !
             else
                loc_dis_sed(is) = loc_new_sed(is)
             end if
          end if
       end DO
    end select
    ! error-catching of negative dissoluiton: return rain flux back to ocean
    If (loc_dis_sed(is_CaCO3) < -const_real_nullsmall) then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ((sed_dep(is) == is_CaCO3) .OR. (sed_type(is) == par_sed_type_CaCO3)) then
             loc_dis_sed(is) = loc_new_sed(is)
          end if
       end DO
       CALL sub_report_error( &
            & 'sedgem_box','sub_update_sed','loc_dis_sed(is_CaCO3) < 0.0 cm3', &
            & 'CONTINUING', &
            & (/real(dum_i), real(dum_j), loc_dis_sed(is_CaCO3) &
            & /),.FALSE. &
            & )
    end if
    ! *** FORCING OF 'RAPID' (HA HA) OCEAN EQUAILIBRIUM WITH CaCO3
    ! NOTE: a flux to the ocean of 0.01 mol (CaCO3 or Ca2+) cm-2 per time step is approximately equal to:
    !       1 g = 0.37 cm3 solid CaCO3 or ~ 1.5 cm equivalent of erosion (of porous 100 wt% CaCO3 sedimentary material)
    if (par_sed_ohmegamin > const_real_nullsmall) then
       if (sed_carb(ic_ohm_cal,dum_i,dum_j) < par_sed_ohmegamin) then
          sedocn_fnet(io_Ca,dum_i,dum_j)  = par_sed_ohmegamin_flux
          sedocn_fnet(io_ALK,dum_i,dum_j) = 2.0*sedocn_fnet(io_Ca,dum_i,dum_j)
          IF (.NOT. ctrl_sed_forcedohmega_ca) then
             ! adjust chemistry by carbonate-only rather than silicate-only weathering
             sedocn_fnet(io_DIC,dum_i,dum_j) = sedocn_fnet(io_Ca,dum_i,dum_j)
             sedocn_fnet(io_DIC_13C,dum_i,dum_j) = &
                  & fun_calc_isotope_fraction(0.0,const_standards(sed_type(is_CaCO3_13C)))*sedocn_fnet(io_DIC,dum_i,dum_j)
          end IF
       end if
    end if

    IF (ctrl_misc_debug4) print*,'*** diagenesis - opal dissolution ***'
    ! *** diagenesis - opal dissolution ***
    ! select opal daigenesis scheme
    select case (par_sed_diagen_opalopt)
    case (                      &
         & 'ridgwelletal2003lookup' &
         & )
       CALL calc_sed_dis_opal(        &
            & dum_dtyr,               &           
            & dum_sfcsumocn(io_T),    &
            & dum_sfcsumocn(io_SiO2), &
            & loc_dis_sed(:),         &
            & loc_new_sed(:),         &
            & sed_top(:,dum_i,dum_j)  &
            & )
    case default
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ((sed_dep(is) == is_opal) .OR. (sed_type(is) == par_sed_type_opal)) then
             if (sed_type(is) == par_sed_type_scavenged) then
                loc_dis_sed(is) = 0.0
                ! deal with how particle-reactive elements are left in the sediments (i.e., what do they stick on?) ...
                ! ### <INSERT CODE> ############################################################################################## !
                ! 
                ! ################################################################################################################ !
             else
                loc_dis_sed(is) = loc_new_sed(is)
             end if
          end if
       end DO
    end select
    ! default and error-catching of negative dissoluiton: return rain flux back to ocean
    If (loc_dis_sed(is_opal) < -const_real_nullsmall) then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ((sed_dep(is) == is_opal) .OR. (sed_type(is) == par_sed_type_opal)) then
             loc_dis_sed(is) = loc_new_sed(is)
          end if
       end DO
       CALL sub_report_error( &
            & 'sedgem_box','sub_update_sed','loc_dis_sed(is_opal) < 0.0 cm3', &
            & 'CONTINUING', &
            & (/real(dum_i), real(dum_j), loc_dis_sed(is_opal) &
            & /),.FALSE. &
            & )
    end if

    ! calculate volume of removed material (as SOILD matter. i.e., zero porosity), in units of cm3 (cm-2)
    ! NOTE: nutrients associated with organic carbon (POP, PON, POFe) have only a 'virtual volume' and so are not counted
    ! NOTE: ditto for isotope tracers
    loc_dis_sed_vol = fun_calc_sed_vol(loc_dis_sed(:))
    ! catch negative sediment dissolution
    ! NOTE: checks introduced individually for POC, CaCO3, opal should make this impossible. But you never know!!!
    IF (loc_dis_sed_vol < -const_real_nullsmall) THEN
       CALL sub_report_error( &
            & 'sedgem_box','sub_update_sed','sediment dissolution < 0.0 cm3', &
            & 'STOPPING', &
            & (/real(dum_i), real(dum_j), loc_dis_sed_vol, &
            & loc_dis_sed(is_POC),   &
            & loc_dis_sed(is_CaCO3), &
            & loc_dis_sed(is_opal),  &
            & loc_dis_sed(is_det)    &
            & /),.TRUE. &
            & )
    END IF

    IF (ctrl_misc_debug3) print*,'(d) update sediment stack'
    ! *** (d) update sediment stack
    !         add the new sediment to the top sediment layer, and deduct the calculated dissolved material
    !         NOTE: all sediment volume (cm3) is of SOLID material (i.e., as if porosity was zero)
    !         NOTE: all sediment thickness (cm) is actual thickness (taking into account the porosity of the sediments)
    ! calculate maximum potential sediment thickness change (taking into account porosity),
    ! and test if the net (rain - dis) thickness of sedimentating material is > 1.0 cm yr-1, or
    ! net (dis - rain) > 1.0 cm yr-1 (i.e., not about to try and remove too much)
    ! if so - take the simplest response - reject all sediment input and set dissolution = rain
    ! NOTE: the 1 cm limit arises because of the way in which excess sedimentary material 
    !        is removed from the top layer and added to the sediment stack, which has layers of thickness 1.0 cm
    loc_sed_top_dth = (loc_new_sed_vol - loc_dis_sed_vol)/(1.0 - max(par_sed_poros_det,par_sed_poros_CaCO3))
    IF ((loc_sed_top_dth > 1.0) .OR. (loc_sed_top_dth < -1.0)) THEN
       loc_dis_sed(:) = loc_new_sed(:)
       loc_dis_sed_vol = loc_new_sed_vol
       loc_sed_top_dth = 0.0
    end IF
    ! update surface ('top') mixed layer sediment composition and calculate temporary surface layer volume
    ! (i.e., after net rain input minus dissolution, before any exchange with underlying sediments)
    sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) + loc_new_sed(:) - loc_dis_sed(:)
    ! calculate surface layer porosity
    ! NOTE: calculate porosity as a function of the VOLUME (not weight) fraction of CaCO3
    loc_frac_CaCO3_top = sed_top(is_CaCO3,dum_i,dum_j)/fun_calc_sed_vol(sed_top(:,dum_i,dum_j))
    loc_sed_poros_top = fun_calc_sed_poros_nsur(loc_frac_CaCO3_top,par_sed_top_th)
    ! calculate thickness of material to be exchanged
    loc_sed_top_dth = fun_calc_sed_vol(sed_top(:,dum_i,dum_j))/(1.0 - loc_sed_poros_top) - par_sed_top_th
    loc_exe_sed_th = ABS(loc_sed_top_dth)
    ! calculate sub-layer number and thickness of top (incomplete) sub-layer of sediment stack
    loc_n_sed_stack_top  = INT(sed_top_h(dum_i,dum_j)) + 1
    loc_sed_stack_top_th = sed_top_h(dum_i,dum_j) - REAL(loc_n_sed_stack_top - 1)
    ! set porosity ratio (to convert the thickness of stack material into an equivalent thickness of surface material)
    loc_r_sed_por = fun_calc_r_sed_por(loc_frac_CaCO3_top,par_sed_top_th)
    ! keep thickness of top layer = par_sed_top_th by transfer to/from sediment stack
    ! (by calculating what sedimentary material needs to be exchanged exchanged - exe_sed(:))
    !   => remove material to the sediment stack below if loc_sed_top_dth > 0.0 cm, or
    !   => add material from the sediment stack below if loc_sed_top_dth < 0.0 cm
    ! (and do nothing if there is no net change in surface sediment layer thickness)
    ! NOTE: take no action if loc_sed_top_dth is zero (to avoid potential 'divide-by-zero' problems)
    IF (loc_sed_top_dth > const_real_nullsmall) THEN
       ! ADD material to the sediment stack
       ! test thickness of top (incomplete) sub-layer compared with thickness of material to be added;
       !   => if exchange th < remaining (unfilled) thickness of top sub-layer of sediment stack, then
       !      add required material to top sub-layer only 
       !   => if exchange th >= remaining (unfilled) thickness of top sub-layer of sediment stack, then
       !      add sufficient material to top sub-layer to fill it plus additional material to next layer up
       ! update sediment surface layer
       loc_frac_CaCO3_top = sed_top(is_CaCO3,dum_i,dum_j)/fun_calc_sed_vol(sed_top(:,dum_i,dum_j))
       loc_r_sed_por = fun_calc_r_sed_por(loc_frac_CaCO3_top,par_sed_top_th)
       IF (loc_exe_sed_th < (loc_r_sed_por*(1.0 - loc_sed_stack_top_th))) THEN
          ! calculate material to be exchanged
          loc_exe_sed(:) = (loc_exe_sed_th/(par_sed_top_th + loc_exe_sed_th))*sed_top(:,dum_i,dum_j)
          ! update sediment stack
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) + loc_exe_sed(:)
       ELSE
          ! calculate material to be exchanged and update sediment stack
          loc_exe_sed(:) = (loc_exe_sed_th/(par_sed_top_th + loc_exe_sed_th))*sed_top(:,dum_i,dum_j)
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) + &
               & (loc_r_sed_por*(1.0 - loc_sed_stack_top_th)/loc_exe_sed_th)*loc_exe_sed(:)
          sed(:,dum_i,dum_j,(loc_n_sed_stack_top + 1)) = sed(:,dum_i,dum_j,(loc_n_sed_stack_top + 1)) + &
               & (1.0 - (loc_r_sed_por*(1.0 - loc_sed_stack_top_th)/loc_exe_sed_th))*loc_exe_sed(:)
          loc_flag_stackgrow = .TRUE.
       ENDIF
       ! deduct sediment material from the sediment surface ('top') layer
       sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) - loc_exe_sed(:)
       ! update sediment stack height
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) + loc_exe_sed_th/loc_r_sed_por
    elseif (loc_sed_top_dth < -const_real_nullsmall) then
       ! REMOVE material from the sediment stack
       ! test thickness of top (incomplete) sub-layer compared with thickness of material to be removed;
       ! - if exchange vol <= thickness of top sub-layer of sediment stack, then
       !   remove required material from the top stack sub-layer only 
       ! - if exchange vol > thickness of top sub-layer of sediment stack, then
       !   remove all material in the top stack sub-layer, plus additional material from next layer down
       ! NOTE: the porosity of the sediment stack must be taken into account
       ! NOTE: a small error in the transfer of material from the sediment stack to the surface layer during times of erosion
       !       will occur due to the use of a single porosity ratio factor,
       !       based on the CaCO3 fraction in the uppermost (incomplete) stack layer rather than the first full layer below it
       !       (but this simplifies everything considerably ...)
       !       (any small error in surface layer thickness will implicitly be 'corrected' at the next time-step)
       ! check whether there is any material in the uppermost (incomplete) stack layer BEFORE trying to calculate porosity ...
       if (loc_sed_stack_top_th < const_real_nullsmall) then
          loc_frac_CaCO3 = &
               & sed(is_CaCO3,dum_i,dum_j,loc_n_sed_stack_top - 1)/fun_calc_sed_vol(sed(:,dum_i,dum_j,loc_n_sed_stack_top - 1))
       else
          loc_frac_CaCO3 = &
               & sed(is_CaCO3,dum_i,dum_j,loc_n_sed_stack_top)/fun_calc_sed_vol(sed(:,dum_i,dum_j,loc_n_sed_stack_top))
       end if
       loc_r_sed_por = fun_calc_r_sed_por(loc_frac_CaCO3,par_sed_top_th)
       IF (loc_exe_sed_th <= (loc_r_sed_por*loc_sed_stack_top_th)) THEN
          ! calculate material to be exchanged
          loc_exe_sed(:) = (loc_exe_sed_th/(loc_r_sed_por*loc_sed_stack_top_th))*sed(:,dum_i,dum_j,loc_n_sed_stack_top)
          ! update sediment stack
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) - loc_exe_sed(:)
       ELSE
          ! calculate material to be exchanged and update sediment stack
          ! material to be exchanged will be equal to ALL the material in the top (incomplete) sediment stack sub-layer,
          ! plus a proportion of the material in the sub-layer immediately below
          loc_exe_sed(:) = ((loc_exe_sed_th - (loc_r_sed_por*loc_sed_stack_top_th))/loc_r_sed_por)* &
               & sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1))
          sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) = sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) - loc_exe_sed(:)
          loc_exe_sed(:) = loc_exe_sed(:) + sed(:,dum_i,dum_j,loc_n_sed_stack_top)
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = 0.0
          loc_flag_stackshrink = .TRUE.
       ENDIF
       ! add eroded sediment material to the sediment surface ('top') layer
       sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) + loc_exe_sed(:)
       ! update sediment stack height
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) - loc_exe_sed_th/loc_r_sed_por
    else
       ! *** DO NOTHING ***
    ENDIF
    ! store surface porosity
    phys_sed(ips_poros,dum_i,dum_j) = loc_sed_poros_top
    ! update accumulated sediment change
    sed_top_INTdth(dum_i,dum_j) = sed_top_INTdth(dum_i,dum_j) + abs(loc_sed_top_dth)
    ! update local variables of sub-layer number and thickness of top (incomplete) sub-layer of sediment stack
    loc_n_sed_stack_top = INT(sed_top_h(dum_i,dum_j)) + 1
    loc_sed_stack_top_th = sed_top_h(dum_i,dum_j) - REAL(loc_n_sed_stack_top - 1)
    if (loc_flag_stackgrow .OR. loc_flag_stackshrink) then
       dum_flag_save = .TRUE.
    elseif (sed_top_INTdth(dum_i,dum_j) > par_sed_save_top_INTdth) then
       sed_top_INTdth(dum_i,dum_j) = 0.0
       dum_flag_save = .TRUE.
    else
       dum_flag_save = .FALSE.
    endif

    IF (ctrl_misc_debug3) print*,'(e) mix the sediment stack'
    ! *** (e) vertically mix the sediment stack (if bioturbation is selected as an option)
    !        NOTE: although the entire mixable portion of the sediment stack array is passed,
    !              only the activited tracers are mixed in sub_sed_mix
    IF (ctrl_sed_bioturb) THEN
       CALL sub_sed_mix(                                                                    &
            & sed(:,dum_i,dum_j,(loc_n_sed_stack_top - par_n_sed_mix):loc_n_sed_stack_top), &
            & sed_top(:,dum_i,dum_j),                                                       &
            & par_sed_mix_k(:),                                                             &
            & loc_sed_stack_top_th                                                          &
            & )
    ENDIF

    IF (ctrl_misc_debug3) print*,'(f) check the thickness of sediment stack'
    ! *** (f) check the thickness of sediment stack
    !         NOTE: if the sediment stack has reached the last available sub-layer, 
    !               then remove a number of sub-layers equal to 'n_sed_tot_drop' from the bottom, 
    !               and re-index the remaining sublayers starting from the bottom
    IF (loc_n_sed_stack_top == n_sed_tot) THEN
       ! shift sediment down the stack
       sed(:,dum_i,dum_j,1:(n_sed_tot - n_sed_tot_drop)) = sed(:,dum_i,dum_j,(n_sed_tot_drop + 1):n_sed_tot)
       sed(:,dum_i,dum_j,(n_sed_tot - n_sed_tot_drop + 1):n_sed_tot) = 0.0
       ! update sediment height and top sub-layer number
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) - REAL(n_sed_tot_drop)
       loc_n_sed_stack_top = INT(sed_top_h(dum_i,dum_j)) + 1
    ENDIF

    IF (ctrl_misc_debug3) print*,'(g) calculate sediment dissolution flux to ocean'
    ! *** (g) calculate sediment dissolution flux to ocean
    !         NOTE: first, convert flux units from cm3 cm-2 to mol cm-2
    sed_fdis(:,dum_i,dum_j) = conv_sed_cm3_mol(:)*loc_dis_sed(:)
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_tot_i = conv_sed_ocn_i(0,is)
       do loc_i=1,loc_tot_i
          io = conv_sed_ocn_i(loc_i,is)
          sedocn_fnet(io,dum_i,dum_j) = sedocn_fnet(io,dum_i,dum_j) + conv_sed_ocn(io,is)*sed_fdis(is,dum_i,dum_j)
       end do
    end DO

!!$    ! ############################################################################################################################ !
!!$    ! ### FIX THIS UP!!! ######################################################################################################### !
!!$    ! ############################################################################################################################ !
!!$    IF (ctrl_misc_debug3) print*,'(g) deal with low bottom-water [O2]'
!!$    ! *** (g) deal with low bottom-water [O2]
!!$    !         NOTE: because SEDGEM knows nothing about the geometry of the overlying ocean,
!!$    !               the only criterion that can be applied in deciding whether there is sufficient O2 for oxidizing organic matter
!!$    !               is to test whether [O2] is above zero or not.
!!$    !               => if [O2] is zero (or rather, very close to it) then do NO3 then SO4 reduction
!!$    !               => replacement of O2 (oxidation) deficit must be done in its entirety, by either NO3 OR SO4,
!!$    !                  because there is no way of knowing how much NO3 is available in the ocean befoer SO4 has to be used ...
!!$    !         NOTE: because consumption of NO3 of SO4 diffusing into the sediments must be complete,
!!$    !               no fractionation w.r.t. 15N or 34S can occur
!!$    !         NOTE: all fluxes in units of mol cm-2
!!$    If (ocn_select(io_O2) .AND. (dum_sfcsumocn(io_O2) < const_real_nullsmall)) then
!!$       loc_potO2def = -sedocn_fnet(io_O2,dum_i,dum_j)
!!$       if (ocn_select(io_NO3) .AND. (dum_sfcsumocn(io_NO3) > const_real_nullsmall)) then
!!$          loc_r15N = dum_sfcsumocn(io_NO3_15N)/dum_sfcsumocn(io_NO3)
!!$          sedocn_fnet(io_NO3,dum_i,dum_j) = -(2.0/3.0)*loc_potO2def
!!$          sedocn_fnet(io_N2,dum_i,dum_j)  = 0.5*(2.0/3.0)*loc_potO2def
!!$          sedocn_fnet(io_O2,dum_i,dum_j)  = 0.0
!!$          sedocn_fnet(io_NO3_15N,dum_i,dum_j) = -loc_r15N*(2.0/3.0)*loc_potO2def
!!$          sedocn_fnet(io_N2_15N,dum_i,dum_j)  = loc_r15N*0.5*(2.0/3.0)*loc_potO2def
!!$       else if (ocn_select(io_SO4) .AND. (dum_sfcsumocn(io_SO4) > const_real_nullsmall)) then
!!$          loc_r34S = dum_sfcsumocn(io_SO4_34S)/dum_sfcsumocn(io_SO4)
!!$          sedocn_fnet(io_SO4,dum_i,dum_j) = -0.5*loc_potO2def
!!$          sedocn_fnet(io_H2S,dum_i,dum_j) = 0.5*loc_potO2def
!!$          sedocn_fnet(io_O2,dum_i,dum_j)  = 0.0
!!$          sedocn_fnet(io_SO4_34S,dum_i,dum_j) = -loc_r34S*0.5*loc_potO2def
!!$          sedocn_fnet(io_H2S_34S,dum_i,dum_j) = loc_r34S*0.5*loc_potO2def
!!$       end if
!!$    end If
!!$    ! ############################################################################################################################ !
!!$    ! ############################################################################################################################ !
!!$    ! ############################################################################################################################ !

    ! *** DEBUG ***
    ! finally ... print some dull debugging info if 'iopt_sed_debug2' option is selected
    IF (ctrl_misc_debug2) THEN
       if ((dum_i == par_misc_debug_i) .AND. (dum_j == par_misc_debug_j)) then
          PRINT*,'---'
          PRINT*,'dum_i,dum_j'
          PRINT*,dum_i,dum_j
          print*,'D,T,S,DIC,PO4,SiO2,ALK,B,Ca,SO4,F'
          print*,dum_D,dum_sfcsumocn(io_T),dum_sfcsumocn(io_S), &
               & dum_sfcsumocn(io_DIC),dum_sfcsumocn(io_PO4),dum_sfcsumocn(io_SiO2),dum_sfcsumocn(io_ALK), &
               & dum_sfcsumocn(io_B),dum_sfcsumocn(io_Ca),dum_sfcsumocn(io_SO4),dum_sfcsumocn(io_F)
          print*,'sed_carb(:)'
          print*,sed_carb(:,dum_i,dum_j)
          PRINT*,'new_sed(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),loc_new_sed(is)
          end do
          PRINT*,'dis_sed(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),loc_dis_sed(is)
          end do
          PRINT*,'sed_top(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),sed_top(is,par_misc_debug_i,par_misc_debug_j)
          end do
          PRINT*,'new sed age = ',loc_new_sed(is_CaCO3_age) / (loc_new_sed(is_CaCO3) + 1.0E-14)
          PRINT*,'new sed thickness = ',loc_new_sed_vol / (1.0 - loc_sed_poros_top)
          PRINT*,'actual sedtop thickness = ',                           &
               & (                                                       &
               &   sed_top(is_CaCO3,par_misc_debug_i,par_misc_debug_j) + &
               &   sed_top(is_opal,par_misc_debug_i,par_misc_debug_j)  + &
               &   sed_top(is_det,par_misc_debug_i,par_misc_debug_j)   + &
               &   sed_top(is_POC,par_misc_debug_i,par_misc_debug_j)     &
               & ) / (1.0 - loc_sed_poros_top)
          PRINT*,'n_sed_stack_top  = ',INT(sed_top_h(par_misc_debug_i,par_misc_debug_j)) + 1
          PRINT*,'sed_stack_top_th = ',sed_top_h(par_misc_debug_i,par_misc_debug_j) - &
               & REAL(((INT(sed_top_h(par_misc_debug_i,par_misc_debug_j)) + 1) - 1))
          PRINT*,'---'
       end if
       IF ((MINVAL(loc_new_sed(:)) < 0.0) .OR. (MINVAL(loc_dis_sed(:)) < 0.0)) then
          PRINT*,dum_i,dum_j
          print*,'MINVAL(loc_new_sed(:)) < 0.0 or MINVAL(loc_dis_sed(:)) < 0.0'
          print*,'loc_new_sed(:); ',loc_new_sed(:)
          print*,'loc_dis_sed(:); ',loc_dis_sed(:)
          PRINT*,'======='
          STOP
       ENDIF
    endif
    ! Ooo - alright then, more debug it is ...
    IF (ctrl_misc_debug1) then
       IF (loc_sed_top_dth > 1.0) THEN
          loc_frac_CaCO3 = (loc_new_sed(is_CaCO3) - loc_dis_sed(is_CaCO3))/fun_calc_sed_vol(loc_new_sed(:) - loc_dis_sed(:))
          loc_sed_poros = fun_calc_sed_poros_nsur(loc_frac_CaCO3,par_sed_top_th)
          CALL sub_report_error(                                                                 &
               & 'sedgem_box','sub_update_sed','magnitude of net sediment vol change > 1.0 cm3', &
               & 'CONTINUING',                                                                   &
               & (/real(dum_i),real(dum_j),loc_sed_top_dth,                                      &
               & (loc_new_sed(is_POC)-loc_dis_sed(is_POC))     / (1.0 - loc_sed_poros),          &
               & (loc_new_sed(is_CaCO3)-loc_dis_sed(is_CaCO3)) / (1.0 - loc_sed_poros),          &
               & (loc_new_sed(is_opal)-loc_dis_sed(is_opal))   / (1.0 - loc_sed_poros),          &
               & (loc_new_sed(is_det)-loc_dis_sed(is_det))     / (1.0 - loc_sed_poros)           &
               & /),.false.                                                                      &
               & )
       END IF
       IF (loc_sed_top_dth < -1.0) THEN
          CALL sub_report_error(                                                                  &
               & 'sedgem_box','sub_update_sed','magnitude of net sediment vol change < -1.0 cm3', &
               & 'CONTINUING',                                                                    &
               & (/real(dum_i),real(dum_j),loc_sed_top_dth,                                       &
               & (loc_new_sed(is_POC)-loc_dis_sed(is_POC))     / (1.0 - loc_sed_poros),           &
               & (loc_new_sed(is_CaCO3)-loc_dis_sed(is_CaCO3)) / (1.0 - loc_sed_poros),           &
               & (loc_new_sed(is_opal)-loc_dis_sed(is_opal))   / (1.0 - loc_sed_poros),           &
               & (loc_new_sed(is_det)-loc_dis_sed(is_det))     / (1.0 - loc_sed_poros)            &
               & /),.false.                                                                       &
               & )
       end IF
    end IF

  END SUBROUTINE sub_update_sed
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! UPDATE SEDIMENT COMPOSITION - SHALLOW WATER CARBONATES (REEFS)
  ! [scheme modelled vaguely on sub_update_sed (above)]
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_update_sed_reef( &
       & dum_dtyr,                &
       & dum_i,dum_j,             &
       & dum_D,                   &
       & dum_sfcsumocn,           &
       & dum_flag_save            &
       & )
    IMPLICIT NONE
    ! dummy arguments
    REAL,INTENT(in)::dum_dtyr                                    ! time-step
    integer,INTENT(in)::dum_i,dum_j                              ! grid point (i,j)
    REAL,INTENT(in)::dum_D                                       ! depth
    real,DIMENSION(n_ocn),intent(in)::dum_sfcsumocn              ! ocean composition interface array
    logical,intent(inout)::dum_flag_save                         ! save sedcore environmental conditions?
    ! local variables
    INTEGER::l,is,io                                         ! tracer index counters
    integer::loc_i,loc_tot_i                                 ! array index conversion variables
    INTEGER::loc_n_sed_stack_top                             ! sediment stack top (incomplete) layer number
    REAL::loc_new_sed_vol                                    ! new sediment volume (as SOLID material)
    REAL::loc_dis_sed_vol                                    ! dis sediment volume (as SOLID material)
    REAL::loc_sed_top_dth                                    ! potential change in sediment surface ('top' layer) thickness
    REAL::loc_exe_sed_th                                     ! exchanged sediment thickness (w.r.t. surface sediment porosity)
    REAL::loc_sed_stack_top_th                               ! sediment stack top thickness (i.e., of the incomplete sub-layer) 
    real::loc_sed_poros_top                                  ! 
    real::loc_r_sed_por                                      ! thickness ratio due to porosity differences (sediment stack / surface layer)
    real::loc_frac_CaCO3                                     ! 
    real::loc_frac_CaCO3_top                                 ! 
    REAL,DIMENSION(n_sed)::loc_new_sed                       ! new (sedimenting) top layer material
    REAL,DIMENSION(n_sed)::loc_dis_sed                       ! remineralized top layer material
    REAL,DIMENSION(n_sed)::loc_exe_sed                       ! top layer material to be exchanged with stack
    real::loc_delta_CaCO3
    real::loc_alpha                                          ! 
    real::loc_R,loc_r7Li                                     ! local isotope R, local 7Li r
    real::loc_ohm                                            ! 

    ! *** INITIALIZE VARIABLES ****************************************************************************************************
    ! zero local sediment arrays
    loc_new_sed(:) = 0.0
    loc_dis_sed(:) = 0.0
    loc_exe_sed(:) = 0.0
    ! initialize relevant location in global sediment dissolution results array
    sed_fdis(:,dum_i,dum_j) = 0.0

    ! *** CALCULATE SEDIMENT RAIN FLUX ********************************************************************************************
    !     dissolve all particulate rain
    !     NOTE: units of mol cm-2
    !     NOTE: do not set the <loc_dis_sed> array for returning all particulate rain -- set <sedocn_fnet> directly
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_tot_i = conv_sed_ocn_i(0,is)
       do loc_i=1,loc_tot_i
          io = conv_sed_ocn_i(loc_i,is)
          sedocn_fnet(io,dum_i,dum_j) = sedocn_fnet(io,dum_i,dum_j) + conv_sed_ocn(io,is)*sed_fsed(is,dum_i,dum_j)
       end do
       sed_fsed(is,dum_i,dum_j) = 0.0
    end DO

    ! *** CALCULATE DISSOLUTION FLUX **********************************************************************************************
    !     NOTE: assume zero for now ...
    ! set zero volume of dissolving material
    loc_dis_sed(:)  = 0.0
    loc_dis_sed_vol = 0.0

    ! *** CALCULATE CaCO3 PRECIPITATION *******************************************************************************************
    ! 
    ! select for calcite vs. aragonite precipitation
    if (par_sed_reef_calcite) then
       loc_ohm = sed_carb(ic_ohm_cal,dum_i,dum_j)
    else
       loc_ohm = sed_carb(ic_ohm_arg,dum_i,dum_j)
    end if
    ! calculate CaCO3 precipitation
    ! NOTE: CaCO3 precipitation calculated as mol cm-2 per time-step
    !       (precipitation rate scaling constants are mol cm-2 yr-1)
    if (loc_ohm > 1.0) then
       sed_fsed(is_CaCO3,dum_i,dum_j) = &
	   & dum_dtyr*par_sed_reef_CaCO3precip_sf*(loc_ohm - 1.0)**par_sed_reef_CaCO3precip_exp
    elseif (loc_ohm > par_sed_CaCO3_abioticohm_min) then
       sed_fsed(is_CaCO3,dum_i,dum_j) = &
	   & dum_dtyr*par_sed_reef_CaCO3precip_sf*(loc_ohm - 1.0)**par_sed_reef_CaCO3precip_exp + &
	   & dum_dtyr*par_sed_CaCO3precip_sf*(loc_ohm - par_sed_CaCO3_abioticohm_min)**par_sed_CaCO3precip_exp
    else
       sed_fsed(is_CaCO3,dum_i,dum_j) = 0.0
    end if
	
    if (sed_select(is_CaCO3_13C)) then
       ! re-calculate carbonate system isotopic properties
       if (ocn_select(io_DIC_13C)) then
          call sub_calc_carb_r13C(           &
               & dum_sfcsumocn(io_T),        &
               & dum_sfcsumocn(io_DIC),      &
               & dum_sfcsumocn(io_DIC_13C),  &
               & sed_carb(:,dum_i,dum_j),    &
               & sed_carbisor(:,dum_i,dum_j) &
               & )
       end IF
       ! calculate 13C/12C fractionation between DIC and CaCO3
       ! NOTE: T-dependent fractionation for calcite following Mook [1986]
       ! NOTE: CaCO3 fractionation w.r.t. HCO3-
       loc_delta_CaCO3 = 15.10 - 4232.0/dum_sfcsumocn(io_T)
       loc_alpha = 1.0 + loc_delta_CaCO3/1000.0
       loc_R = sed_carbisor(ici_HCO3_r13C,dum_i,dum_j)/(1.0 - sed_carbisor(ici_HCO3_r13C,dum_i,dum_j))
       sed_fsed(is_CaCO3_13C,dum_i,dum_j) = (loc_alpha*loc_R/(1.0 + loc_alpha*loc_R))*sed_fsed(is_CaCO3,dum_i,dum_j)
    end if
    ! calculate Li incorporation and isotopes (if selected)	
    if (ocn_select(io_Li) .AND. ocn_select(io_Ca) .AND. sed_select(is_LiCO3)) then
       sed_fsed(is_LiCO3,dum_i,dum_j) = sed_fsed(is_CaCO3_13C,dum_i,dum_j)* &
            & (par_bio_red_CaCO3_LiCO3 + par_bio_red_CaCO3_LiCO3_alpha*dum_sfcsumocn(io_Li)/dum_sfcsumocn(io_Ca))
       if (sed_select(is_LiCO3_7Li)) then
          ! calculate 7/6Li fractionation between Li and LiCO3
          if (dum_sfcsumocn(io_Li)> const_real_nullsmall) then
             loc_r7Li = dum_sfcsumocn(io_Li_7Li)/dum_sfcsumocn(io_Li)
          else
             loc_r7Li = 0.0
          end if
          loc_alpha = 1.0 + par_d7Li_LiCO3_epsilon/1000.0
          loc_R = loc_r7Li/(1.0 - loc_r7Li)
          sed_fsed(is_LiCO3_7Li,dum_i,dum_j) = sed_fsed(is_LiCO3,dum_i,dum_j)* &
               & (loc_alpha*loc_R/(1.0 + loc_alpha*loc_R))
       end if
    end if
    ! calculate volume of precipitated carbonate
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_new_sed(is) = conv_sed_mol_cm3(is)*sed_fsed(is,dum_i,dum_j)
    end do
    ! calculate volume of added material (as SOILD matter. i.e., assuming zero porosity), in units of cm3 (cm-2)
    loc_new_sed_vol = fun_calc_sed_vol(loc_new_sed(:))
    ! set exchange with ocean to account for removed of solutes via precipitation
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_tot_i = conv_sed_ocn_i(0,is)
       do loc_i=1,loc_tot_i
          io = conv_sed_ocn_i(loc_i,is)
          sedocn_fnet(io,dum_i,dum_j) = sedocn_fnet(io,dum_i,dum_j) - conv_sed_ocn(io,is)*sed_fsed(is,dum_i,dum_j)
       end do
    end DO

    ! *** UPDATE REEF SEDIMENT STACK **********************************************************************************************
    !         add the new sediment to the top sediment layer, and deduct the calculated dissolved material
    !         NOTE: all sediment volume (cm3) is of SOLID material (i.e., as if porosity was zero)
    !         NOTE: all sediment thickness (cm) is actual thickness (taking into account the porosity of the sediments)
    !         BUT ... *** assume a fixed porosity for now (ZERO) ***
    ! calculate maximum potential sediment thickness change (taking into account porosity),
    ! and test if the net (rain - dis) thickness of sedimentating material is > 1.0 cm per time-step, or
    ! net (dis - rain) > 1.0 cm per time-step (i.e., not about to try and remove too much)
    ! if so - take the simplest response - reject all sediment input and set dissolution = rain
    ! NOTE: the 1 cm limit arises because of the way in which excess sedimentary material 
    !        is removed from the top layer and added to the sediment stack, which has layers of thickness 1.0 cm
    loc_sed_top_dth = (loc_new_sed_vol - loc_dis_sed_vol)/(1.0 - par_sed_poros_CaCO3_reef)
    IF (loc_sed_top_dth > 1.0) THEN
       loc_dis_sed(:) = loc_dis_sed(:) + ((loc_sed_top_dth - 1.0)/(loc_new_sed_vol/(1.0 - par_sed_poros_CaCO3_reef)))*loc_new_sed(:)
       loc_dis_sed_vol = loc_dis_sed_vol + (1.0 - par_sed_poros_CaCO3_reef)*(loc_sed_top_dth - 1.0)
       loc_sed_top_dth = 1.0
    elseif (loc_sed_top_dth < -1.0) then
       loc_dis_sed(:) = loc_dis_sed(:) + ((-1.0 - loc_sed_top_dth)/(loc_dis_sed_vol/(1.0 - par_sed_poros_CaCO3_reef)))* &
            & loc_dis_sed(:)
       loc_dis_sed_vol = loc_dis_sed_vol + (1.0 - par_sed_poros_CaCO3_reef)*(-1.0 - loc_sed_top_dth)
       loc_sed_top_dth = -1.0
    end IF
    ! update surface ('top') mixed layer sediment composition and calculate temporary surface layer volume
    ! (i.e., after net rain input minus dissolution, before any exchange with underlying sediments)
    sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) + loc_new_sed(:) - loc_dis_sed(:)
    ! calculate surface layer porosity
    ! NOTE: calculate porosity as a function of the VOLUME (not weight) fraction of CaCO3
    loc_frac_CaCO3_top = sed_top(is_CaCO3,dum_i,dum_j)/fun_calc_sed_vol(sed_top(:,dum_i,dum_j))
    loc_sed_poros_top = par_sed_poros_CaCO3_reef
    ! calculate thickness of material to be exchanged
    loc_sed_top_dth = fun_calc_sed_vol(sed_top(:,dum_i,dum_j))/(1.0 - loc_sed_poros_top) - par_sed_top_th
    loc_exe_sed_th = ABS(loc_sed_top_dth)
    ! calculate sub-layer number and thickness of top (incomplete) sub-layer of sediment stack
    loc_n_sed_stack_top  = INT(sed_top_h(dum_i,dum_j)) + 1
    loc_sed_stack_top_th = sed_top_h(dum_i,dum_j) - REAL(loc_n_sed_stack_top - 1)
    ! set porosity ratio (to convert the thickness of stack material into an equivalent thickness of surface material)
    loc_r_sed_por = 1.0
    ! keep thickness of top layer = par_sed_top_th by transfer to/from sediment stack
    ! (by calculating what sedimentary material needs to be exchanged exchanged - exe_sed(:))
    !   => remove material to the sediment stack below if loc_sed_top_dth > 0.0 cm, or
    !   => add material from the sediment stack below if loc_sed_top_dth < 0.0 cm
    ! (and do nothing if there is no net change in surface sediment layer thickness)
    ! NOTE: take no action if loc_sed_top_dth is zero (to avoid potential 'divide-by-zero' problems)
    IF (loc_sed_top_dth > const_real_nullsmall) THEN
       ! ADD material to the sediment stack
       ! test thickness of top (incomplete) sub-layer compared with thickness of material to be added;
       !   => if exchange th < remaining (unfilled) thickness of top sub-layer of sediment stack, then
       !      add required material to top sub-layer only 
       !   => if exchange th >= remaining (unfilled) thickness of top sub-layer of sediment stack, then
       !      add sufficient material to top sub-layer to fill it plus additional material to next layer up
       ! update sediment surface layer
       loc_frac_CaCO3_top = sed_top(is_CaCO3,dum_i,dum_j)/fun_calc_sed_vol(sed_top(:,dum_i,dum_j))
       loc_r_sed_por = 1.0
       IF (loc_exe_sed_th < (loc_r_sed_por*(1.0 - loc_sed_stack_top_th))) THEN
          ! calculate material to be exchanged
          loc_exe_sed(:) = (loc_exe_sed_th/(par_sed_top_th + loc_exe_sed_th))*sed_top(:,dum_i,dum_j)
          ! update sediment stack
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) + loc_exe_sed(:)
       ELSE
          ! calculate material to be exchanged and update sediment stack
          loc_exe_sed(:) = (loc_exe_sed_th/(par_sed_top_th + loc_exe_sed_th))*sed_top(:,dum_i,dum_j)
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) + &
               & (loc_r_sed_por*(1.0 - loc_sed_stack_top_th)/loc_exe_sed_th)*loc_exe_sed(:)
          sed(:,dum_i,dum_j,(loc_n_sed_stack_top + 1)) = sed(:,dum_i,dum_j,(loc_n_sed_stack_top + 1)) + &
               & (1.0 - (loc_r_sed_por*(1.0 - loc_sed_stack_top_th)/loc_exe_sed_th))*loc_exe_sed(:)
       ENDIF
       ! deduct sediment material from the sediment surface ('top') layer
       sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) - loc_exe_sed(:)
       ! update sediment stack height
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) + loc_exe_sed_th/loc_r_sed_por
    elseif (loc_sed_top_dth < -const_real_nullsmall) then
       ! REMOVE material from the sediment stack
       ! test thickness of top (incomplete) sub-layer compared with thickness of material to be removed;
       ! - if exchange vol <= thickness of top sub-layer of sediment stack, then
       !   remove required material from the top stack sub-layer only 
       ! - if exchange vol > thickness of top sub-layer of sediment stack, then
       !   remove all material in the top stack sub-layer, plus additional material from next layer down
       ! NOTE: the porosity of the sediment stack must be taken into account
       ! NOTE: a small error in the transfer of material from the sediment stack to the surface layer during times of erosion
       !       will occur due to the use of a single porosity ratio factor,
       !       based on the CaCO3 fraction in the uppermost (incomplete) stack layer rather than the first full layer below it
       !       (but this simplifies everything considerably ...)
       !       (any small error in surface layer thickness will implicitly be 'corrected' at the next time-step)
       ! check whether there is any material in the uppermost (incomplete) stack layer BEFORE trying to calculate porosity ...
       if (loc_sed_stack_top_th < const_real_nullsmall) then
          loc_frac_CaCO3 = &
               & sed(is_CaCO3,dum_i,dum_j,loc_n_sed_stack_top - 1)/fun_calc_sed_vol(sed(:,dum_i,dum_j,loc_n_sed_stack_top - 1))
       else
          loc_frac_CaCO3 = &
               & sed(is_CaCO3,dum_i,dum_j,loc_n_sed_stack_top)/fun_calc_sed_vol(sed(:,dum_i,dum_j,loc_n_sed_stack_top))
       end if
       loc_r_sed_por = 1.0
       IF (loc_exe_sed_th <= (loc_r_sed_por*loc_sed_stack_top_th)) THEN
          ! calculate material to be exchanged
          loc_exe_sed(:) = (loc_exe_sed_th/(loc_r_sed_por*loc_sed_stack_top_th))*sed(:,dum_i,dum_j,loc_n_sed_stack_top)
          ! update sediment stack
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) - loc_exe_sed(:)
       ELSE
          ! calculate material to be exchanged and update sediment stack
          ! material to be exchanged will be equal to ALL the material in the top (incomplete) sediment stack sub-layer,
          ! plus a proportion of the material in the sub-layer immediately below
          loc_exe_sed(:) = ((loc_exe_sed_th - (loc_r_sed_por*loc_sed_stack_top_th))/loc_r_sed_por)* &
               & sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1))
          sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) = sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) - loc_exe_sed(:)
          loc_exe_sed(:) = loc_exe_sed(:) + sed(:,dum_i,dum_j,loc_n_sed_stack_top)
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = 0.0
       ENDIF
       ! add eroded sediment material to the sediment surface ('top') layer
       sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) + loc_exe_sed(:)
       ! update sediment stack height
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) - loc_exe_sed_th/loc_r_sed_por
    ENDIF
    ! store surface porosity
    phys_sed(ips_poros,dum_i,dum_j) = loc_sed_poros_top
    ! update accumulated sediment change
    sed_top_INTdth(dum_i,dum_j) = sed_top_INTdth(dum_i,dum_j) + abs(loc_sed_top_dth)
    ! update local variables of sub-layer number and thickness of top (incomplete) sub-layer of sediment stack
    loc_n_sed_stack_top = INT(sed_top_h(dum_i,dum_j)) + 1
    loc_sed_stack_top_th = sed_top_h(dum_i,dum_j) - REAL(loc_n_sed_stack_top - 1)
    if (sed_top_INTdth(dum_i,dum_j) > par_sed_save_top_INTdth) then
       sed_top_INTdth(dum_i,dum_j) = 0.0
       dum_flag_save = .TRUE.
    else
       dum_flag_save = .FALSE.
    endif

    ! *** CHECK HEIGHT OF SEDIMENT STACK ******************************************************************************************
    !     adjust the sediment stack if full
    IF (loc_n_sed_stack_top == n_sed_tot) THEN
       ! shift sediment down the stack
       sed(:,dum_i,dum_j,1:(n_sed_tot - n_sed_tot_drop)) = sed(:,dum_i,dum_j,(n_sed_tot_drop + 1):n_sed_tot)
       sed(:,dum_i,dum_j,(n_sed_tot - n_sed_tot_drop + 1):n_sed_tot) = 0.0
       ! update sediment height and top sub-layer number
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) - REAL(n_sed_tot_drop)
       loc_n_sed_stack_top = INT(sed_top_h(dum_i,dum_j)) + 1
    ENDIF

    ! *** DEBUG *******************************************************************************************************************
    !     finally ... print some dull debugging info if 'iopt_sed_debug2' option is selected
    IF (ctrl_misc_debug2) THEN
       if ((dum_i == par_misc_debug_i) .AND. (dum_j == par_misc_debug_j)) then
          PRINT*,'---'
          PRINT*,'dum_i,dum_j'
          PRINT*,dum_i,dum_j
          print*,'D,T,S,DIC,PO4,SiO2,ALK,B,Ca,SO4,F'
          print*,dum_D,dum_sfcsumocn(io_T),dum_sfcsumocn(io_S), &
               & dum_sfcsumocn(io_DIC),dum_sfcsumocn(io_PO4),dum_sfcsumocn(io_SiO2),dum_sfcsumocn(io_ALK), &
               & dum_sfcsumocn(io_B),dum_sfcsumocn(io_Ca),dum_sfcsumocn(io_SO4),dum_sfcsumocn(io_F)
          print*,'sed_carb(:)'
          print*,sed_carb(:,dum_i,dum_j)
          PRINT*,'new_sed(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),loc_new_sed(is)
          end do
          PRINT*,'dis_sed(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),loc_dis_sed(is)
          end do
          PRINT*,'sed_top(:)'
          do is=1,n_sed
             if (sed_select(is)) print*,is,string_sed(is),sed_top(is,par_misc_debug_i,par_misc_debug_j)
          end do
          PRINT*,'new sed age = ',loc_new_sed(is_CaCO3_age) / (loc_new_sed(is_CaCO3) + 1.0E-14)
          PRINT*,'new sed thickness = ',loc_new_sed_vol / (1.0 - loc_sed_poros_top)
          PRINT*,'actual sedtop thickness = ',                           &
               & (                                                       &
               &   sed_top(is_CaCO3,par_misc_debug_i,par_misc_debug_j) + &
               &   sed_top(is_opal,par_misc_debug_i,par_misc_debug_j)  + &
               &   sed_top(is_det,par_misc_debug_i,par_misc_debug_j)   + &
               &   sed_top(is_POC,par_misc_debug_i,par_misc_debug_j)     &
               & ) / (1.0 - loc_sed_poros_top)
          PRINT*,'n_sed_stack_top  = ',INT(sed_top_h(par_misc_debug_i,par_misc_debug_j)) + 1
          PRINT*,'sed_stack_top_th = ',sed_top_h(par_misc_debug_i,par_misc_debug_j) - &
               & REAL(((INT(sed_top_h(par_misc_debug_i,par_misc_debug_j)) + 1) - 1))
          PRINT*,'---'
       end if
    endif

  END SUBROUTINE sub_update_sed_reef
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! UPDATE SEDIMENT COMPOSITION - MUD!!!
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_update_sed_mud( &
       & dum_i,dum_j,            &
       & dum_flag_save           &
       & )
    IMPLICIT NONE
    ! dummy arguments
    integer,INTENT(in)::dum_i,dum_j                              ! grid point (i,j)
    logical,intent(inout)::dum_flag_save                         ! save sedcore environmental conditions?
    ! local variables
    INTEGER::l,is,io                                             ! tracer index counters
    integer::loc_i,loc_tot_i                                     ! array index conversion variables
    INTEGER::loc_n_sed_stack_top                                 ! sediment stack top (incomplete) layer number
    REAL::loc_new_sed_vol                                        ! new sediment volume (as SOLID material)
    REAL::loc_dis_sed_vol                                        ! dis sediment volume (as SOLID material)
    REAL::loc_sed_top_dth                                        ! potential change in sediment surface ('top' layer) thickness
    REAL::loc_exe_sed_th                                         ! exchanged sediment thickness (w.r.t. surface sediment porosity)
    REAL::loc_sed_stack_top_th                                   ! sediment stack top thickness (i.e., of the incomplete sub-layer) 
    real::loc_sed_poros_top                                      ! 
    real::loc_r_sed_por                                          ! thickness ratio due to porosity differences
    ! (sediment stack / surface layer)
    real::loc_frac_CaCO3                                         ! 
    real::loc_frac_CaCO3_top                                     ! 
    REAL,DIMENSION(n_sed)::loc_new_sed                           ! new (sedimenting) top layer material
    REAL,DIMENSION(n_sed)::loc_dis_sed                           ! remineralized top layer material
    REAL,DIMENSION(n_sed)::loc_exe_sed                           ! top layer material to be exchanged with stack

    ! *** INITIALIZE VARIABLES ****************************************************************************************************
    ! zero local sediment arrays
    loc_new_sed(:) = 0.0
    loc_dis_sed(:) = 0.0
    loc_exe_sed(:) = 0.0
    ! initialize relevant location in global sediment dissolution results array
    sed_fdis(:,dum_i,dum_j) = 0.0

    ! *** CALCULATE SEDIMENT RAIN FLUX ********************************************************************************************
    !     calculate new sedimenting material to be added to the sediment top layer
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_new_sed(is) = conv_sed_mol_cm3(is)*sed_fsed(is,dum_i,dum_j)
    end do
    ! calculate volume of added material (as SOILD matter. i.e., zero porosity), in units of cm3 (cm-2)
    ! NOTE: nutrients associated with organic carbon (POP, PON, POFe) have only a 'virtual volume' and so are not counted
    ! NOTE: ditto for isotope tracers
    ! NOTE: check that rain thickness of sedimentating material is >= 0.0 cm
    loc_new_sed_vol = fun_calc_sed_vol(loc_new_sed(:))
    IF (loc_new_sed_vol < -const_real_nullsmall) THEN
       CALL sub_report_error( &
            & 'sedgem_box','sub_update_sed','sediment input < 0.0 cm3', &
            & 'STOPPING', &
            & (/real(dum_i), real(dum_j), loc_new_sed_vol, &
            & loc_new_sed(is_POC),   &
            & loc_new_sed(is_CaCO3), &
            & loc_new_sed(is_opal),  &
            & loc_new_sed(is_det)    &
            & /),.TRUE. &
            & )
    END IF

    ! *** CALCULATE DISSOLUTION FLUX **********************************************************************************************
    !     NOTE: dissolve all particulate rain
    !     NOTE: units of mol cm-2
    !     NOTE: do not set the <loc_dis_sed> array for returning all particulate rain -- set <sedocn_fnet> directly
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       if ( &
            & ((sed_dep(is) == is_CaCO3) .OR. (sed_type(is) == par_sed_type_CaCO3)) &
            & .OR. &
            & ((sed_dep(is) == is_opal) .OR. (sed_type(is) == par_sed_type_opal)) &
            & .OR. &
            & ((sed_dep(is) == is_POC) .OR. (sed_type(is) == par_sed_type_POM)) &
            & ) then
          loc_dis_sed(is) = loc_new_sed(is)
       end if
    end DO
    ! calculate volume of removed material (as SOILD matter. i.e., zero porosity), in units of cm3 (cm-2)
    ! NOTE: nutrients associated with organic carbon (POP, PON, POFe) have only a 'virtual volume' and so are not counted
    ! NOTE: ditto for isotope tracers
    loc_dis_sed_vol = fun_calc_sed_vol(loc_dis_sed(:))
    IF (loc_dis_sed_vol < -const_real_nullsmall) THEN
       CALL sub_report_error( &
            & 'sedgem_box','sub_update_sed','sediment dissolution < 0.0 cm3', &
            & 'STOPPING', &
            & (/real(dum_i), real(dum_j), loc_dis_sed_vol, &
            & loc_dis_sed(is_POC),   &
            & loc_dis_sed(is_CaCO3), &
            & loc_dis_sed(is_opal),  &
            & loc_dis_sed(is_det)    &
            & /),.TRUE. &
            & )
    END IF

    ! *** UPDATE SEDIMENT STACK ***************************************************************************************************
    !         add the new sediment to the top sediment layer, and deduct the calculated dissolved material
    !         NOTE: all sediment volume (cm3) is of SOLID material (i.e., as if porosity was zero)
    !         NOTE: all sediment thickness (cm) is actual thickness (taking into account the porosity of the sediments)
    !         BUT ... *** assume a fixed porosity for now (ZERO) ***
    ! calculate maximum potential sediment thickness change (taking into account porosity),
    ! and test if the net (rain - dis) thickness of sedimentating material is > 1.0 cm per time-step, or
    ! net (dis - rain) > 1.0 cm per time-step (i.e., not about to try and remove too much)
    ! if so - take the simplest response - reject all sediment input and set dissolution = rain
    ! NOTE: the 1 cm limit arises because of the way in which excess sedimentary material 
    !        is removed from the top layer and added to the sediment stack, which has layers of thickness 1.0 cm
    loc_sed_top_dth = (loc_new_sed_vol - loc_dis_sed_vol)/(1.0 - par_sed_poros_det)
    IF (loc_sed_top_dth > 1.0) THEN
       loc_dis_sed(:) = loc_dis_sed(:) + ((loc_sed_top_dth - 1.0)/(loc_new_sed_vol/(1.0 - par_sed_poros_CaCO3_reef)))*loc_new_sed(:)
       loc_dis_sed_vol = loc_dis_sed_vol + (1.0 - par_sed_poros_CaCO3_reef)*(loc_sed_top_dth - 1.0)
       loc_sed_top_dth = 1.0
    elseif (loc_sed_top_dth < -1.0) then
       loc_dis_sed(:) = loc_dis_sed(:) + ((-1.0 - loc_sed_top_dth)/(loc_dis_sed_vol/(1.0 - par_sed_poros_CaCO3_reef)))* &
            & loc_dis_sed(:)
       loc_dis_sed_vol = loc_dis_sed_vol + (1.0 - par_sed_poros_CaCO3_reef)*(-1.0 - loc_sed_top_dth)
       loc_sed_top_dth = -1.0
    end IF
    ! update surface ('top') mixed layer sediment composition and calculate temporary surface layer volume
    ! (i.e., after net rain input minus dissolution, before any exchange with underlying sediments)
    sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) + loc_new_sed(:) - loc_dis_sed(:)
    ! calculate surface layer porosity
    ! NOTE: calculate porosity as a function of the VOLUME (not weight) fraction of CaCO3
    loc_frac_CaCO3_top = sed_top(is_CaCO3,dum_i,dum_j)/fun_calc_sed_vol(sed_top(:,dum_i,dum_j))
    loc_sed_poros_top = par_sed_poros_CaCO3_reef
    ! calculate thickness of material to be exchanged
    loc_sed_top_dth = fun_calc_sed_vol(sed_top(:,dum_i,dum_j))/(1.0 - loc_sed_poros_top) - par_sed_top_th
    loc_exe_sed_th = ABS(loc_sed_top_dth)
    ! calculate sub-layer number and thickness of top (incomplete) sub-layer of sediment stack
    loc_n_sed_stack_top  = INT(sed_top_h(dum_i,dum_j)) + 1
    loc_sed_stack_top_th = sed_top_h(dum_i,dum_j) - REAL(loc_n_sed_stack_top - 1)
    ! set porosity ratio (to convert the thickness of stack material into an equivalent thickness of surface material)
    loc_r_sed_por = 1.0
    ! keep thickness of top layer = par_sed_top_th by transfer to/from sediment stack
    ! (by calculating what sedimentary material needs to be exchanged exchanged - exe_sed(:))
    !   => remove material to the sediment stack below if loc_sed_top_dth > 0.0 cm, or
    !   => add material from the sediment stack below if loc_sed_top_dth < 0.0 cm
    ! (and do nothing if there is no net change in surface sediment layer thickness)
    ! NOTE: take no action if loc_sed_top_dth is zero (to avoid potential 'divide-by-zero' problems)
    IF (loc_sed_top_dth > const_real_nullsmall) THEN
       ! ADD material to the sediment stack
       ! test thickness of top (incomplete) sub-layer compared with thickness of material to be added;
       !   => if exchange th < remaining (unfilled) thickness of top sub-layer of sediment stack, then
       !      add required material to top sub-layer only 
       !   => if exchange th >= remaining (unfilled) thickness of top sub-layer of sediment stack, then
       !      add sufficient material to top sub-layer to fill it plus additional material to next layer up
       ! update sediment surface layer
       loc_frac_CaCO3_top = sed_top(is_CaCO3,dum_i,dum_j)/fun_calc_sed_vol(sed_top(:,dum_i,dum_j))
       loc_r_sed_por = 1.0
       IF (loc_exe_sed_th < (loc_r_sed_por*(1.0 - loc_sed_stack_top_th))) THEN
          ! calculate material to be exchanged
          loc_exe_sed(:) = (loc_exe_sed_th/(par_sed_top_th + loc_exe_sed_th))*sed_top(:,dum_i,dum_j)
          ! update sediment stack
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) + loc_exe_sed(:)
       ELSE
          ! calculate material to be exchanged and update sediment stack
          loc_exe_sed(:) = (loc_exe_sed_th/(par_sed_top_th + loc_exe_sed_th))*sed_top(:,dum_i,dum_j)
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) + &
               & (loc_r_sed_por*(1.0 - loc_sed_stack_top_th)/loc_exe_sed_th)*loc_exe_sed(:)
          sed(:,dum_i,dum_j,(loc_n_sed_stack_top + 1)) = sed(:,dum_i,dum_j,(loc_n_sed_stack_top + 1)) + &
               & (1.0 - (loc_r_sed_por*(1.0 - loc_sed_stack_top_th)/loc_exe_sed_th))*loc_exe_sed(:)
       ENDIF
       ! deduct sediment material from the sediment surface ('top') layer
       sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) - loc_exe_sed(:)
       ! update sediment stack height
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) + loc_exe_sed_th/loc_r_sed_por
    elseif (loc_sed_top_dth < -const_real_nullsmall) then
       ! REMOVE material from the sediment stack
       ! test thickness of top (incomplete) sub-layer compared with thickness of material to be removed;
       ! - if exchange vol <= thickness of top sub-layer of sediment stack, then
       !   remove required material from the top stack sub-layer only 
       ! - if exchange vol > thickness of top sub-layer of sediment stack, then
       !   remove all material in the top stack sub-layer, plus additional material from next layer down
       ! NOTE: the porosity of the sediment stack must be taken into account
       ! NOTE: a small error in the transfer of material from the sediment stack to the surface layer during times of erosion
       !       will occur due to the use of a single porosity ratio factor,
       !       based on the CaCO3 fraction in the uppermost (incomplete) stack layer rather than the first full layer below it
       !       (but this simplifies everything considerably ...)
       !       (any small error in surface layer thickness will implicitly be 'corrected' at the next time-step)
       ! check whether there is any material in the uppermost (incomplete) stack layer BEFORE trying to calculate porosity ...
       if (loc_sed_stack_top_th < const_real_nullsmall) then
          loc_frac_CaCO3 = &
               & sed(is_CaCO3,dum_i,dum_j,loc_n_sed_stack_top - 1)/fun_calc_sed_vol(sed(:,dum_i,dum_j,loc_n_sed_stack_top - 1))
       else
          loc_frac_CaCO3 = &
               & sed(is_CaCO3,dum_i,dum_j,loc_n_sed_stack_top)/fun_calc_sed_vol(sed(:,dum_i,dum_j,loc_n_sed_stack_top))
       end if
       loc_r_sed_por = 1.0
       IF (loc_exe_sed_th <= (loc_r_sed_por*loc_sed_stack_top_th)) THEN
          ! calculate material to be exchanged
          loc_exe_sed(:) = (loc_exe_sed_th/(loc_r_sed_por*loc_sed_stack_top_th))*sed(:,dum_i,dum_j,loc_n_sed_stack_top)
          ! update sediment stack
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = sed(:,dum_i,dum_j,loc_n_sed_stack_top) - loc_exe_sed(:)
       ELSE
          ! calculate material to be exchanged and update sediment stack
          ! material to be exchanged will be equal to ALL the material in the top (incomplete) sediment stack sub-layer,
          ! plus a proportion of the material in the sub-layer immediately below
          loc_exe_sed(:) = ((loc_exe_sed_th - (loc_r_sed_por*loc_sed_stack_top_th))/loc_r_sed_por)* &
               & sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1))
          sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) = sed(:,dum_i,dum_j,(loc_n_sed_stack_top - 1)) - loc_exe_sed(:)
          loc_exe_sed(:) = loc_exe_sed(:) + sed(:,dum_i,dum_j,loc_n_sed_stack_top)
          sed(:,dum_i,dum_j,loc_n_sed_stack_top) = 0.0
       ENDIF
       ! add eroded sediment material to the sediment surface ('top') layer
       sed_top(:,dum_i,dum_j) = sed_top(:,dum_i,dum_j) + loc_exe_sed(:)
       ! update sediment stack height
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) - loc_exe_sed_th/loc_r_sed_por
    ENDIF
    ! update local variables of sub-layer number and thickness of top (incomplete) sub-layer of sediment stack
    loc_n_sed_stack_top = INT(sed_top_h(dum_i,dum_j)) + 1
    loc_sed_stack_top_th = sed_top_h(dum_i,dum_j) - REAL(loc_n_sed_stack_top - 1)
    ! store surface porosity
    phys_sed(ips_poros,dum_i,dum_j) = loc_sed_poros_top
    ! set no sedenv save
    dum_flag_save = .FALSE.

    ! *** CHECK HEIGHT OF SEDIMENT STACK ******************************************************************************************
    !     adjust the sediment stack if full
    IF (loc_n_sed_stack_top == n_sed_tot) THEN
       ! shift sediment down the stack
       sed(:,dum_i,dum_j,1:(n_sed_tot - n_sed_tot_drop)) = sed(:,dum_i,dum_j,(n_sed_tot_drop + 1):n_sed_tot)
       sed(:,dum_i,dum_j,(n_sed_tot - n_sed_tot_drop + 1):n_sed_tot) = 0.0
       ! update sediment height and top sub-layer number
       sed_top_h(dum_i,dum_j) = sed_top_h(dum_i,dum_j) - REAL(n_sed_tot_drop)
       loc_n_sed_stack_top = INT(sed_top_h(dum_i,dum_j)) + 1
    ENDIF

    ! *** CALCULATE SEDIMENT DISSOLUTION FLUX TO THE OCEAN ************************************************************************
    !     NOTE: will be zero currently, but include for completeness
    !     NOTE: first, convert flux units from cm3 cm-2 to mol cm-2
    sed_fdis(:,dum_i,dum_j) = conv_sed_cm3_mol(:)*loc_dis_sed(:)
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_tot_i = conv_sed_ocn_i(0,is)
       do loc_i=1,loc_tot_i
          io = conv_sed_ocn_i(loc_i,is)
          sedocn_fnet(io,dum_i,dum_j) = sedocn_fnet(io,dum_i,dum_j) + conv_sed_ocn(io,is)*sed_fdis(is,dum_i,dum_j)
       end do
    end DO

  END SUBROUTINE sub_update_sed_mud
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! CALCULATE SEDIMENT CACO3 DIAGENESIS
  subroutine sub_calc_sed_dis_CaCO3( &
       & dum_dtyr,                   &
       & dum_D,                      &
       & dum_dCO3_cal,               &
       & dum_sed_diagen_fPOCfrac,    &
       & dum_ocn,                    &
       & dum_carbconst,              &
       & dum_carb,                   &
       & dum_sed_dis,                &
       & dum_sed_new,                &
       & dum_sed_top,                &
       & dum_sed_mix_k               &
       & )
    IMPLICIT NONE
    ! dummy arguments
    real,intent(in)::dum_dtyr
    real,intent(in)::dum_D,dum_dCO3_cal,dum_sed_diagen_fPOCfrac
    real,intent(in),DIMENSION(n_ocn)::dum_ocn
    REAL,intent(in),DIMENSION(n_carbconst)::dum_carbconst
    REAL,INTENT(in),DIMENSION(n_carb)::dum_carb
    REAL,INTENT(inout),DIMENSION(n_sed)::dum_sed_dis
    REAL,INTENT(in),DIMENSION(n_sed)::dum_sed_new
    REAL,INTENT(in),DIMENSION(n_sed)::dum_sed_top
    real,intent(in)::dum_sed_mix_k
    ! local variables
    INTEGER::is,l              ! COUNTERS
    REAL::loc_frac_CaCO3_top   ! dry mass fraction CaCO3 in top sediment layer
    REAL::loc_fPOC             ! POC rain rate driving diagenesis (and CaCO3 dissolution) (mol cm-2 a-1)
    REAL::loc_sed_wt_top       ! sediment mass (in top sediment)
    REAL::loc_sed_dis_new      ! fraction dissolved from new sediment
    REAL::loc_sed_dis_top      ! fraction dissolved from top sediment
    REAL::loc_dis              ! raw dissolution flux from lookup table
    real::loc_O2               ! local bottom-water oxygen concentration
    !REAL,DIMENSION(par_nn_input)::loc_nn_inp,loc_nn_inpnorm     ! neural network variables

    ! *** USER-DEFINABLE OPTIONS ***
    ! NOTE: settings not included in the run-time configuration files for clarity
    ! ******************************
    par_sed_diagen_fPOCmax = 100.0e-6 ! max Corg flux allowed to drive carbonate dissolution (mol cm-2 yr-1) [50.0]
    ! ******************************

    ! *** calculate sediment diagenesis defining variables ***
    ! NOTE: the units of the Corg flux must be changed from (cm3 cm-2) to (mol cm-2 yr-1),
    !       in order to calculate CaCO3 diagenesis
    ! NOTE: the lookup table is only generated over a range of POC flux up to 50 (umol cm-2 yr-1)
    !       (set by the parameter 'lookup_fCorg_max')
    !       so care must be taken when exceeding this range (hence the purpose of 'par_sed_diagen_fPOCmax')
    !       and enough [O2] may not in actual fact be available to oxidize such a flux in the sediments fully
    ! mass fraction (dry) of calcite in sediments
    loc_sed_wt_top = fun_calc_sed_mass(dum_sed_top(:))
    loc_frac_CaCO3_top = conv_cal_cm3_g*dum_sed_top(is_CaCO3)/loc_sed_wt_top
    ! Corg rain rate
    loc_fPOC = dum_sed_diagen_fPOCfrac*conv_POC_cm3_mol*dum_sed_new(is_POC)/dum_dtyr
    IF (loc_fPOC > par_sed_diagen_fPOCmax) loc_fPOC = par_sed_diagen_fPOCmax
    ! prevent negative oxygen concentrations being passed
    loc_O2 = dum_ocn(io_O2)
    if (loc_O2 < const_real_nullsmall) loc_O2 = 0.0

    ! *** calculate potential calcite dissolution flux ***
    ! NOTE: the dissolution flux is calculated in units of (mol cm-2 yr-1);
    !       this must be converted to cm3 cm-2 per time step
    ! NOTE: catch (sigular matrix) stability failures of the Archer scheme
    !       => assume no CaCO3 preservation (100% dissoution of rain flux)
    select case (par_sed_diagen_CaCO3opt)
    case ('archer1991explicit')
       ! CALCULATE SEDIMENT CACO3 DIAGENESIS EXPLICITLY
       ! NOTE: model of Archer [1991]
       loc_dis = fun_archer1991_sedflx(                                                           &
            & loc_O2,loc_frac_CaCO3_top,loc_fPOC,                                                 &
            & dum_carb(ic_conc_CO2),dum_carb(ic_conc_HCO3),dum_carb(ic_conc_CO3),                 &
            & dum_carbconst(icc_k1),dum_carbconst(icc_k2),dum_carbconst(icc_kcal)/dum_ocn(io_Ca), &
            & dum_sed_mix_k                                                                       &
            & )
       loc_dis = conv_cal_mol_cm3*loc_dis*dum_dtyr
       if (error_Archer) loc_dis = dum_sed_new(is_CaCO3)
    case ('ridgwell2001lookup')
       ! CALCULATE SEDIMENT CACO3 DIAGENESIS VIA A LOOK-UP TABLE [Ridgwell, 2001]
       loc_dis = fun_interp_4D(                                                  &
            & lookup_sed_dis_cal,dum_D,dum_dCO3_cal,loc_frac_CaCO3_top,loc_fPOC, &
            & lookup_D_max,lookup_dCO3_max,lookup_frac_max,lookup_fCorg_max,     &
            & lookup_i_D_min,lookup_i_D_max,                                     &
            & lookup_i_dCO3_min,lookup_i_dCO3_max,                               &
            & lookup_i_frac_min,lookup_i_frac_max,                               &
            & lookup_i_fCorg_min,lookup_i_fCorg_max                              &
            & )
       loc_dis = conv_cal_mol_cm3*loc_dis*dum_dtyr
    case ('ridgwell2001nn')
!!$       ! CALCULATE SEDIMENT CACO3 DIAGENESIS VIA A NEURAL NETWORK (#1)
!!$       ! NOTE: underlying model is Archer [1991]
!!$       loc_nn_inp(1) = dum_D
!!$       loc_nn_inp(2) = dum_dCO3_cal
!!$       loc_nn_inp(3) = 200.0e-6
!!$       loc_nn_inp(4) = loc_frac_CaCO3_top
!!$       loc_nn_inp(5) = loc_fPOC
!!$       ! NOTE: the nn_minp, nn_maxp scaling limits appear to have screwy values in the netCDF file
!!$       !       => use har-coded limits (nn_min, nn_max) set in sedgem_nnutils for now
!!$       loc_nn_inpnorm(:) = fun_nn_parnorm(loc_nn_inp(:),nn_min(:),nn_max(:),par_nn_input)
!!$       call sub_nn_forward_model(loc_nn_inpnorm)
!!$       loc_dis = fun_nn_posmn(nn_a(1),nn_mint,nn_maxt)
!!$       loc_dis = conv_cal_mol_cm3*conv_umol_mol*loc_dis*dum_dtyr
    case ('archer1991nn')
!!$       ! CALCULATE SEDIMENT CACO3 DIAGENESIS VIA A NEURAL NETWORK (#1)
!!$       ! NOTE: underlying model is Archer [1991]
!!$       loc_nn_inp(1) = dum_D
!!$       loc_nn_inp(2) = dum_dCO3_cal
!!$       loc_nn_inp(3) = 200.0e-6
!!$       loc_nn_inp(4) = loc_frac_CaCO3_top
!!$       loc_nn_inp(5) = loc_fPOC
!!$       ! NOTE: the nn_minp, nn_maxp scaling limits appear to have screwy values in the netCDF file
!!$       !       => use har-coded limits (nn_min, nn_max) set in sedgem_nnutils for now
!!$       loc_nn_inpnorm(:) = fun_nn_parnorm(loc_nn_inp(:),nn_min(:),nn_max(:),par_nn_input)
!!$       call sub_nn_forward_model(loc_nn_inpnorm)
!!$       loc_dis = fun_nn_posmn(nn_a(1),nn_mint,nn_maxt)
!!$       loc_dis = conv_cal_mol_cm3*conv_umol_mol*loc_dis*dum_dtyr
    case ('archer2002mudsnn')
       ! CALCULATE SEDIMENT CACO3 DIAGENESIS VIA A NEURAL NETWORK (#2)
       ! NOTE: underlying model is Archer et al. [2002]
       ! *********************
       ! *** <INSERT CODE> ***
       ! *********************
    case default
       ! 
    end select

    ! *** calculate actual dissolution flux ***
    IF (loc_dis < const_real_nullsmall) THEN
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ( &
               & (sed_dep(is) == is_CaCO3) .OR. &
               & (sed_type(is) == par_sed_type_CaCO3) .OR. &
               & (sed_type(sed_dep(is)) == par_sed_type_CaCO3) &
               & ) then
             dum_sed_dis(is) = 0.0
          end if
       end do
    else
       ! cap maximum dissolution
       IF (loc_dis >= (dum_sed_new(is_CaCO3) + dum_sed_top(is_CaCO3))) THEN
          dum_sed_dis(is_CaCO3) = dum_sed_new(is_CaCO3) + dum_sed_top(is_CaCO3)
       ELSE
          dum_sed_dis(is_CaCO3) = loc_dis
       ENDIF
       ! calculate dissolution components from new (rain) and old (core-top) sediments
       ! NOTE: there are TWO possible end-member models for where the carbonate dissolution takes place
       !       (see Ridgwell [2001]; www.seao2.org/publications/ridgwell_thesis.pdf )
       !       the code for both is provided - simply comment out the one that is not wanted (and re-compile) to change over
       ! loc_sed_dis_new       - is the amount of newly arrived carbonate that is dissolved 
       ! loc_sed_dis_top       - is the amount of pre-existing carbonate dissolved from the surface ('top') layer
       ! dum_sed_dis(is_CaCO3) - is the total amount of carbonate that needs to be dissolved
       ! option #1
       ! \/\/\/ CODE FOR INTERFACE DISSOLUTION \/\/\/
       IF (dum_sed_dis(is_CaCO3) >= dum_sed_new(is_CaCO3)) THEN
          loc_sed_dis_new = dum_sed_new(is_CaCO3)
          loc_sed_dis_top = dum_sed_dis(is_CaCO3) - dum_sed_new(is_CaCO3)
       else
          loc_sed_dis_new = dum_sed_dis(is_CaCO3)
          loc_sed_dis_top = 0.0
       end IF
       ! /\/\/\ ****************************** /\/\/\
       !       ! option #2
       !       ! \/\/\/ ALTERNATIVE CODE FOR HOMOGENEOUS DISSOLUTION \/\/\/
       !       ! What it is saying is determine proportion of carbonate that has newly arrived at the sediments (dum_sed_new(is_CaCO3))
       !       ! compared to the total amount of carbonate available for dissolution (dum_sed_new(is_CaCO3) + dum_sed_top(is_CaCO3)).
       !       ! Dissolve carbonate in this proportion - this should have the same effect as if the 'new' carbonate had been 
       !       ! mixed into the surface layer and then dissolution calculated.
       !       ! NOTE: this subroutine is only called if there is some CaCO3 (so dividing-by-zero cannot occur)
       !       loc_sed_new_frac = dum_sed_new(is_CaCO3)/(dum_sed_new(is_CaCO3) + dum_sed_top(is_CaCO3))
       !       loc_sed_dis_new = loc_sed_new_frac*dum_sed_dis(is_CaCO3)
       !       loc_sed_dis_top = (1.0 - loc_sed_new_frac)*dum_sed_dis(is_CaCO3)
       !       ! /\/\/\ ******************************************** /\/\/\
       ! calculate isotope and 'age' dissolution fluxes
       ! NOTE: assume no fractionation associated with dissolution
       ! NOTE: generic routine also includes age and foram tracers
       ! NOTE: ensure that bulk CaCO3 is not re-processed (it has itself as its dependency)
       ! NOTE: assume that particle-reactive elements remain in sediments
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ( &
               & (sed_dep(is) == is_CaCO3) .OR. &
               & (sed_type(is) == par_sed_type_CaCO3) .OR. &
               & (sed_type(sed_dep(is)) == par_sed_type_CaCO3) &
               & ) then
             if ((is /= is_CaCO3) .AND. (sed_type(is) /= par_sed_type_scavenged)) then
                dum_sed_dis(is) = 0.0
                if (dum_sed_new(is_CaCO3) > const_real_nullsmall) then
                   dum_sed_dis(is) = dum_sed_dis(is) + (dum_sed_new(is)/dum_sed_new(is_CaCO3))*loc_sed_dis_new
                end if
                if (dum_sed_top(is_CaCO3) > const_real_nullsmall) then
                   dum_sed_dis(is) = dum_sed_dis(is) + (dum_sed_top(is)/dum_sed_top(is_CaCO3))*loc_sed_dis_top
                end if
             end if
          end if
       end do
    end IF

  END subroutine sub_calc_sed_dis_CaCO3
  ! ****************************************************************************************************************************** !


  ! ********************************************************************************************************************************
  ! CALCULATE SEDIMENT OPAL DIAGENESIS
  SUBROUTINE calc_sed_dis_opal( &
       & dum_dtyr,              &
       & dum_T,                 &
       & dum_SiO2,              &
       & dum_sed_dis,           &
       & dum_sed_new,           &
       & dum_sed_top            &
       & )
    IMPLICIT NONE
    ! dummy arguments
    real,intent(in)::dum_dtyr
    real::dum_T,dum_SiO2
    REAL,INTENT(inout),DIMENSION(n_sed)::dum_sed_dis
    REAL,INTENT(in),DIMENSION(n_sed)::dum_sed_new
    REAL,INTENT(in),DIMENSION(n_sed)::dum_sed_top
    ! local variables
    INTEGER::is,l              ! COUNTERS
    REAL::loc_frac_opal        ! volume fraction of opal in the sediment (cm3 cm-3 solids)
    REAL::loc_KSi              ! estimated opal dissolution constant
    REAL::loc_opaltorefrac     ! %refrac/%opal
    REAL::loc_sed_dis_new      ! fraction dissolved from new sediment
    REAL::loc_sed_dis_top      ! fraction dissolved from top sediment
    REAL::loc_dis              ! raw (opal) dissolution flux from lookup table
    REAL::loc_sed_vol_top      ! sediment solids volume (in top sediment)

    ! *** USER-DEFINABLE OPTIONS ***
    ! NOTE: settings not included in the run-time configuration files for clarity
    ! ******************************
    opt_sed(iopt_sed_diagen_AltoasymSi) = .TRUE. ! asymptotic [Si] dependence on %refrac/%opal?
    opt_sed(iopt_sed_diagen_AltoKSi)    = .TRUE. ! KSi dependence on %refrac/%opal?
    par_sed_opal_Sitoopalmax            = 15.0   ! asymptotic [Si] %refrac/%opal max limit
    par_sed_opal_KSi0                   = 0.1    ! fixed opal KSi value (yr-1)
    ! ******************************

    ! *** calculate sediment diagenesis defining variables ***
    ! volume fraction of opal in sediments
    loc_sed_vol_top = fun_calc_sed_vol(dum_sed_top(:))
    loc_frac_opal = dum_sed_top(is_opal)/loc_sed_vol_top
    ! make asymptotic [Si] dependent on the ratio of sediment (top) refrac to opal (wt%/wt%) (if selected)
    ! NOTE: cap ratio at 'par_sed_opal_Sitoopalmax'
    IF (opt_sed(iopt_sed_diagen_AltoasymSi)) THEN
       if (dum_sed_top(is_opal) < const_real_nullsmall) then
          loc_opaltorefrac = 0.0
       elseIF ((conv_det_cm3_g*dum_sed_top(is_det)) > (par_sed_opal_Sitoopalmax*conv_opal_cm3_g*dum_sed_top(is_opal))) THEN
          loc_opaltorefrac = par_sed_opal_Sitoopalmax
       ELSE
          loc_opaltorefrac = (conv_det_cm3_g*dum_sed_top(is_det))/(conv_opal_cm3_g*dum_sed_top(is_opal))
       END IF
    ELSE
       loc_opaltorefrac = 0.0
    ENDIF
    ! base opal dissolution rate constant
    IF (opt_sed(iopt_sed_diagen_AltoKSi)) THEN
       loc_KSi = (0.0500 + 0.0550/((0.0164 + loc_opaltorefrac)**0.75))/conv_yr_s
    ELSE
       loc_KSi = par_sed_opal_KSi0
    ENDIF

    ! *** calculate potential opal dissolution flux ***
    ! NOTE: the dissolution flux is calculated in units of (mol cm-2 yr-1);
    !       this must be converted to cm3 cm-2
    select case (par_sed_diagen_opalopt)
    case ('ridgwelletal2003explicit')
       ! CALCULATE SEDIMENT OPAL DIAGENESIS EXPLICITLY
       ! NOTE: model of Ridgwell [2001], Ridgwell et al. [2002]
       loc_dis = fun_ridgwelletal2003_sedflx(loc_frac_opal,dum_SiO2,dum_T,loc_KSi,loc_opaltorefrac)
       loc_dis = conv_cal_mol_cm3*loc_dis*dum_dtyr
    case ('ridgwelletal2003lookup')
       ! CALCULATE SEDIMENT OPAL DIAGENESIS VIA A LOOK-UP TABLE [Ridgwell, 2001]
       loc_dis = fun_interp_5D(lookup_sed_dis_opal,loc_frac_opal,dum_SiO2,dum_T,loc_KSi,loc_opaltorefrac, &
            & lookup_opalpc_max,lookup_concSi_max,lookup_T_max,lookup_KSi0_max,lookup_opaltorefrac_max, &
            & lookup_i_opalpc_min,lookup_i_opalpc_max, &
            & lookup_i_concSi_min,lookup_i_concSi_max, &
            & lookup_i_T_min,lookup_i_T_max, &
            & lookup_i_KSi0_min,lookup_i_KSi0_max, &
            & lookup_i_opaltorefrac_min,lookup_i_opaltorefrac_max)
       loc_dis = conv_opal_mol_cm3*loc_dis*dum_dtyr
    case ('archer2002muds')
       ! CALCULATE SEDIMENT OPAL DIAGENESIS VIA A NEURAL NETWORK (#1)
       ! NOTE: underlying model is Ridgwell et al. [2002]
       ! *********************
       ! *** <INSERT CODE> ***
       ! *********************
    case default
       ! CALCULATE SEDIMENT OPAL DIAGENESIS VIA A NEURAL NETWORK (#2)
       ! NOTE: underlying model is Archer et al. [2002]
       ! *********************
       ! *** <INSERT CODE> ***
       ! *********************
    end select

    ! *** calculate actual dissolution flux ***
    ! NOTE: the procedure for calculating 'interface' dissolution of opal
    !       from new sedimenting and top sediment material is:
    !       IF the estimated dissolution flux < 0.0, reset to 0.0
    !           (this is possible in the case of the look-up table,
    !           becuase of the use of extrapolation for values falling outside of the look-up table bounds)
    !       ELSEIF the estimated dissolution flux is greater than or equal to the combined inventories of 
    !           new material and top material,
    !           => cap the dissolution flux at this total
    !       ELSE, dissolution flux is unchanged
    IF (loc_dis < const_real_nullsmall) THEN
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ( &
               & (sed_dep(is) == is_opal) .OR. &
               & (sed_type(is) == par_sed_type_opal) .OR. &
               & (sed_type(sed_dep(is)) == par_sed_type_opal) &
               & ) then
             dum_sed_dis(is) = 0.0
          end if
       end do
    else
       ! cap maximum dissolution
       IF (loc_dis >= (dum_sed_new(is_opal) + dum_sed_top(is_opal))) THEN
          dum_sed_dis(is_opal) = dum_sed_new(is_opal) + dum_sed_top(is_opal)
       ELSE
          dum_sed_dis(is_opal) = loc_dis
       ENDIF
       ! calculate dissolution components from new (rain) and old (core-top) sediments
       IF (dum_sed_dis(is_opal) >= dum_sed_new(is_opal)) THEN
          loc_sed_dis_new = dum_sed_new(is_opal)
          loc_sed_dis_top = dum_sed_dis(is_opal) - dum_sed_new(is_opal)
       else
          loc_sed_dis_new = dum_sed_dis(is_opal)
          loc_sed_dis_top = 0.0
       end IF
       ! calculate isotope and trace metal dissolution fluxes
       ! NOTE: ensure that bulk opal is not re-processed (it has itself as its dependency)
       ! NOTE: assume that particle-reactive elements remain in sediments
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ( &
               & (sed_dep(is) == is_opal) .OR. &
               & (sed_type(is) == par_sed_type_opal) .OR. &
               & (sed_type(sed_dep(is)) == par_sed_type_opal) &
               & ) then
             if ((is /= is_opal) .AND. (sed_type(is) /= par_sed_type_scavenged)) then
                dum_sed_dis(is) = 0.0
                if (dum_sed_new(is_opal) > const_real_nullsmall) then
                   dum_sed_dis(is) = dum_sed_dis(is) + (dum_sed_new(is)/dum_sed_new(is_opal))*loc_sed_dis_new
                end if
                if (dum_sed_top(is_opal) > const_real_nullsmall) then
                   dum_sed_dis(is) = dum_sed_dis(is) + (dum_sed_top(is)/dum_sed_top(is_opal))*loc_sed_dis_top
                end if
             end if
          end if
       end do
    ENDIF

  END SUBROUTINE calc_sed_dis_opal
  ! ********************************************************************************************************************************


  ! ****************************************************************************************************************************** !
  ! CALCULATE SEDIMENT BIOTURBATION MIXING COEFFICIENT PROFILE
  ! NOTE: see: Archer et al. [2002]
  SUBROUTINE sub_calc_sed_mix( &
       & dum_dtyr,             &
       & dum_fPOC,             &
       & dum_O2,               &
       & dum_sed_mix_k         &
       & )
    IMPLICIT NONE
    ! dummy arguments
    real,intent(in)::dum_dtyr
    real,intent(in)::dum_O2
    REAL,INTENT(in)::dum_fPOC
    REAL,INTENT(inout),DIMENSION(0:par_n_sed_mix)::dum_sed_mix_k
    ! local variables
    integer::d                 ! 
    REAL::loc_fPOC             ! POC rain rate driving diagenesis (and CaCO3 dissolution) (umol cm-2 yr-1)
    real::loc_kO2              ! 
    real::loc_mix_k_sur        ! 

    ! #### USER-DEFINABLE OPTIONS ################################################################################################ !
    ! NOTE: settings not included in the run-time configuration files for clarity
    par_sed_mix_zmix  = 5.0        ! depth scale for bioturbation (cm) [Martin and Sayles, 1996]
    !                                NOTE: Archer et al. [2002] gives zmix = 8 cm
    par_sed_mix_c0_O2 = 20.0E-6    ! half saturation constant of O2 for bioturbation (mol kg-1) [Archer et al., 2002]
    ! ############################################################################################################################ !

    ! *** calculate bioturbational mixing profie ***
    ! calculate local variables
    ! NOTE: convert Corg rain rate from (mol cm-2) to (umol cm-2 yr-1)
    loc_fPOC = conv_mol_umol*dum_fPOC/dum_dtyr
    loc_kO2 = dum_O2/(dum_O2 + par_sed_mix_c0_O2)
    ! calculate surface bioturbational mixing rate coefficient
    ! NOTE: restict value to ensure stability
    loc_mix_k_sur = loc_kO2*(0.0232*(loc_fPOC)**0.85)
    if (loc_mix_k_sur < par_sed_mix_k_sur_min) loc_mix_k_sur = par_sed_mix_k_sur_min
    if (loc_mix_k_sur > par_sed_mix_k_sur_max) loc_mix_k_sur = par_sed_mix_k_sur_max
    ! calculate mixing coefficient down the sediment profile
    ! NOTE: d = 0 is the BOTTOM of the mixing profile (not the top!!!)
    DO d = 0,par_n_sed_mix,1
       dum_sed_mix_k(d) = loc_mix_k_sur*exp(-(real(par_n_sed_mix - d)/par_sed_mix_zmix)**2)
    END DO

  END SUBROUTINE sub_calc_sed_mix
  ! ****************************************************************************************************************************** !


  ! ********************************************************************************************************************************
  ! MIX SEDIMENT STACK SUB_LAYERS
  SUBROUTINE sub_sed_mix( &
       & dum_sed, &
       & dum_sed_top, &
       & dum_sed_mix_k, &
       & dum_sed_stack_top_th &
       & )
    IMPLICIT NONE
    ! dummy arguments
    REAL,INTENT(inout),DIMENSION(n_sed,0:par_n_sed_mix)::dum_sed
    REAL,INTENT(inout),DIMENSION(n_sed)::dum_sed_top
    REAL,INTENT(in),DIMENSION(0:par_n_sed_mix)::dum_sed_mix_k
    REAL,INTENT(in)::dum_sed_stack_top_th
    ! local variables
    INTEGER::l,is,d                                              !
    REAL::loc_r_sed_por                                          ! porosity ratio
    REAL::loc_frac_CaCO3                                         ! 
    REAL,DIMENSION(n_l_sed,0:par_n_sed_mix)::loc_sed             ! 
    REAL,DIMENSION(n_l_sed)::loc_sed_top                         ! 
    REAL,DIMENSION(0:par_n_sed_mix)::loc_mix                     ! adapted sediment mixing rate profile
    REAL,DIMENSION(n_l_sed)::loc_exe                             ! exchange sedimentary material
    REAL,DIMENSION(n_l_sed,0:par_n_sed_mix)::loc_dsed            ! sediment composition change
    REAL,DIMENSION(n_l_sed)::loc_dsed_top                        ! sediment top composition change

    ! *** mix the sediment stack ***
    ! NOTE: restrictions on the dimensioning of arrays in f90 requires that 
    !       the partial sediment stack passed into this procedure has 
    !       an index of 'n_sed_mix' for the top (incomplete) sub-layer, and '0' for the bottom
    ! NOTE: the sediment mixing rate profile has a corresponing ordering, with the maximum mixing rate 
    !       (measured at the uppermost surface of the incomplete stack layer) having the highest ('n_sed_mix') index
    ! NOTE: sedimentary layers are sequentially mixed upwardswards in pairs
    ! NOTE: the sediment stack must be mixed in with the top layer, in addition to the sediment stack 
    !       being mixed internally
    ! NOTE: mix only selected tracers
    ! NOTE: mixing exchange (and the sediment exe(:) array) is measured UPWARDS
    ! CONFIGURATION;
    !       (A) sediment mixing profile dum_sed_mix_k(:) has its n_sed_mix index value aligned with the top of the sediment stack
    !           array elements correspond to a seperation of 1 cm (the stack layer resolution) between mixing rate values
    !           (the depths at which the dum_sed_mix_k(:) values are alighed are indicated by '===')
    !       (B) the sediment mixing profile loc_mix(:) is linearly interpolated to the depths seperating each stack layer ('---')
    !           NOTE: strictly, the bottom-most (d = 0) layer should only be partially exchanged with the overlying layer,
    !                 otherwise the total bioturbated depth will vary with the thickness of the uppermost partially-filled layer
    !                 However, because the mixing rate is typically close to zero at the base of the bioturbated zone,
    !                 the error will be negiable
    !       (C) loc_mix(d) therefore defines the mixing between stack layers (d + 1) and d
    !           NOTE: a special case has to be made for d = (n_sed_mix - 1) because the biodiffusion distance is < 1 cm
    !           NOTE: a special case also has to be made for mixing between the d = n_sed_mix stack layer and surface layer
    !         
    !       ------------------------------------------                                                                             !
    !       ******************************************                                                                             !
    !       *** SURFACE SEDIMENT ('TOP') LAYER *******                                                                             !
    !       ******************************************                                                                             !
    !       ------------------------------------------                                                                             !
    !                                                                                                                              !
    !       ------------------------------------------ === sed_mix_k(n_sed_mix)     --- loc_mix(n_sed_mix)                         !
    !       ### PARTIAL STACK LAYER; d = n_sed_mix ### /|\                                                                         !
    !       ------------------------------------------  |                           --- loc_mix(n_sed_mix - 1)                     !
    !       ##########################################  |  1 cm                     /|\                                            !
    !       ### FULL LAYER; d = (n_sed_mix - 1) ###### \|/                           |  1 cm                                       !
    !       ########################################## === sed_mix_k(n_sed_mix - 1) \|/                                            !
    !       ------------------------------------------ /|\                          --- loc_mix(n_sed_mix - 2)                     !
    !       ##########################################  |  1 cm                     /|\                                            !
    !       ### FULL LAYER; d = (n_sed_mix - 2) ###### \|/                           |  1 cm                                       !
    !       ########################################## === sed_mix_k(n_sed_mix - 2) \|/                                            !
    !       ------------------------------------------ /|\                          --- loc_mix(n_sed_mix - 3)                     !
    !       #                                                                       /|\                                            !
    !       #                                          \|/                                                                         !
    !       #                                          ===                          \|/                                            !
    !       ------------------------------------------ /|\                          --- loc_mix(1)                                 !
    !       ##########################################  |  1 cm                     /|\                                            !
    !       ### FULL LAYER; d = 1 #################### \|/                           |  1 cm                                       !
    !       ########################################## === sed_mix_k(1)             \|/                                            !
    !       ------------------------------------------ /|\                          --- loc_mix(0)                                 !
    !       ##########################################  |  1 cm                                                                    !
    !       ### FULL LAYER; d = 0 #################### \|/                                                                         !
    !       ########################################## === sed_mix_k(0)                                                            !
    !       ------------------------------------------                                                                             !
    !       #                                                                                                                      !

    ! (0) copy sediment composition to local variables
    !     -> transform total tracer array to enabled tracer array indices
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_sed(l,:)   = dum_sed(is,:)
       loc_sed_top(l) = dum_sed_top(is)
    END DO
    ! (1) initialize (zero) sediment compositional change variable
    loc_dsed(:,:)   = 0.0
    loc_dsed_top(:) = 0.0
    ! (2) calculate sediment porosities
    loc_frac_CaCO3 = dum_sed_top(is_CaCO3)/fun_calc_sed_vol(dum_sed_top(:))
    loc_r_sed_por = fun_calc_r_sed_por(loc_frac_CaCO3,par_sed_top_th)
    ! (3) re-grid biodiffusion coefficient to boudaries between sediment layers
    !     (d = n_sed_mix being the index corresponding to the boundary between stack and surface ('top') sediment layer)
    d = par_n_sed_mix
    loc_mix(d) = dum_sed_mix_k(par_n_sed_mix)
    DO d = (par_n_sed_mix - 1),0,-1
       loc_mix(d) = (1.0 - dum_sed_stack_top_th)*dum_sed_mix_k(d + 1) + dum_sed_stack_top_th*dum_sed_mix_k(d)
    END DO
    ! (4) calculate mixing exchange between complete stack layers
    !     NOTE: implicit in this calculation is that the distance over which biodiffusion operates is 1 cm
    DO d = (par_n_sed_mix - 2),0,-1
       loc_exe(:) = loc_mix(d)* &
            & (loc_sed(:,d) - loc_sed(:,d + 1))
       loc_dsed(:,d + 1) = loc_dsed(:,d + 1) + loc_exe(:)
       loc_dsed(:,d)     = loc_dsed(:,d)     - loc_exe(:)
    END DO
    ! (5) calculate mixing exchange between uppermost complete stack layer and overlying (incomplete) stack layer
    d = par_n_sed_mix - 1
    loc_exe(:) = (loc_mix(d)/(0.5*(1.0 + dum_sed_stack_top_th)))* &
         & ((dum_sed_stack_top_th*loc_sed(:,d)) - loc_sed(:,d + 1))
    loc_dsed(:,d + 1) = loc_dsed(:,d + 1) + loc_exe(:)
    loc_dsed(:,d)     = loc_dsed(:,d)     - loc_exe(:)
    ! (6) calculate mixing exchange between incomplete stack layer and sediment surface ('top') layer
    !     NOTE: because it is solids that are mixed, rather than the whole sediment (including porosity),
    !           there will be no change in uppermost stack layer and surface sediment layer thicknesses
    !           (because of differing porosities)
    d = par_n_sed_mix
    loc_exe(:) = (loc_mix(d)/(0.5*(par_sed_top_th + dum_sed_stack_top_th)))* &
         & (loc_sed(:,d) - (loc_r_sed_por*dum_sed_stack_top_th/par_sed_top_th)*loc_sed_top(:))
    loc_dsed_top(:) = loc_dsed_top(:) + loc_exe(:)
    loc_dsed(:,d)   = loc_dsed(:,d)   - loc_exe(:)
    ! (7) update mixed sub-layer values
    loc_sed(:,:)   = loc_sed(:,:)   + loc_dsed(:,:)
    loc_sed_top(:) = loc_sed_top(:) + loc_dsed_top(:)
    ! (8) convert tracer array back
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       dum_sed(is,:)   = loc_sed(l,:)
       dum_sed_top(is) = loc_sed_top(l)
    END DO

  END SUBROUTINE sub_sed_mix
  ! ********************************************************************************************************************************


END MODULE sedgem_box
