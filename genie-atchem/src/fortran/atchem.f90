
! ******************************************************************************************************************************** !
! TIMESTEP AtChem
SUBROUTINE atchem(    &
     & dum_dts,       &
     & dum_sfxsumatm, &
     & dum_sfcatm     &
     )
  USE atchem_lib
  USE atchem_box
  IMPLICIT NONE
  ! dummy arguments
  real,intent(in)::dum_dts  
  real,dimension(n_atm,n_i,n_j),intent(inout)::dum_sfxsumatm
  real,dimension(n_atm,n_i,n_j),intent(out)::dum_sfcatm
  ! local variables
  integer::ia,l,i,j                                                     ! 
  real::loc_dtyr                                                        ! 
  REAL::loc_atm_tot_V                                                   ! 
  REAL,DIMENSION(n_atm)::loc_atm_tot                                    !  
  REAL,DIMENSION(n_i,n_j)::loc_conv_atm_mol,loc_conv_mol_atm            !
  REAL,DIMENSION(n_atm,n_i,n_j)::locij_fatm                             ! local flux to atmosphere (mol)
  REAL,DIMENSION(n_atm)::loc_fracdecay_atm                              ! local reduction factor for decaying atmospheric tracers

  ! *** INITIALIZE LOCAL VARIABLES ***
  locij_fatm(:,:,:) = 0.0

  ! *** CALCULATE LOCAL CONSTANTS ***
  ! local constants for converting between partial pressure and molar quantity
  loc_conv_atm_mol(:,:) = phys_atm(ipa_V,:,:)/(conv_Pa_atm*const_R_SI*atm(ia_T,:,:))
  loc_conv_mol_atm(:,:) = 1.0/loc_conv_atm_mol(:,:)
  ! time
  loc_dtyr = dum_dts/conv_yr_s
  ! fractional reduction factors for decaying isotopes
  loc_fracdecay_atm(:) = EXP(-loc_dtyr*const_lambda_atm(:))

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
  ! *** (i,j) GRID PT LOOP START *** !
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
  
  block_iloop: DO i=1,n_i
     block_jloop: DO j=1,n_j 
        
        ! *** DECAY RADIOACTIVE TRACERS ***
        DO l=3,n_l_atm
           ia = conv_iselected_ia(l)
           ! radioactive decay of isotopes
           IF (abs(const_lambda_atm(ia)).gt.const_real_nullsmall) THEN
              atm(ia,i,j) = loc_fracdecay_atm(ia)*atm(ia,i,j)
           END if
        end do
        
        ! *** OXIDIZE CH4 ***
        IF (atm_select(ia_pCH4) .AND. atm_select(ia_pCO2) .AND. atm_select(ia_pO2)) THEN
           CALL sub_calc_oxidize_CH4(i,j,loc_dtyr)
        END IF 

        ! *** ADD CH4 ***
        IF (atm_select(ia_pCH4) .AND. atm_select(ia_pCO2) .AND. atm_select(ia_pO2)) THEN
           CALL sub_calc_wetlands_CH4(loc_dtyr,locij_fatm(:,i,j))
        END IF
        
!!$        ! *** OXIDIZE H2S ***
!!$        IF (atm_select(ia_pH2S)) THEN
!!$           CALL sub_calc_oxidize_H2S(i,j,loc_dtyr)
!!$        END IF
        
        ! *** PRODUCE RADIOACTIVE TRACERS ***
        IF (atm_select(ia_pCO2_14C)) THEN
           CALL sub_calc_generate_14C(loc_dtyr,locij_fatm(:,i,j))
        END IF
        
     END DO block_jloop
  END DO block_iloop
  
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
  ! *** (i,j) GRID PT LOOP END *** !
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

  ! *** UPDATE ATMOSPHERIC COMPOSITION ***
  ! set internal atmospheric flux
  fatm(:,:,:) = dum_sfxsumatm(:,:,:)
  ! NOTE: flux <fatm> in (mol m-2 per timestep)
  ! update atmospheric composition
  ! NOTE: units of partial pressure (atm)
  ! NOTE: carry out at every (i.e, wet + dry) grid point
  DO l=3,n_l_atm
     ia = conv_iselected_ia(l)
     ! update atmospheric tracers
     atm(ia,:,:) = atm(ia,:,:) + loc_conv_mol_atm(:,:)*phys_atm(ipa_A,:,:)*fatm(ia,:,:) + loc_conv_mol_atm(:,:)*locij_fatm(ia,:,:)
     ! <HACK TO HOMOGENIZE ATMOSPHERIC COMPOSITION>
     ! homogenize the partial pressure of tracers in the atmopshere across (all grid points)
     loc_atm_tot(ia) = SUM(loc_conv_atm_mol(:,:)*atm(ia,:,:))
     loc_atm_tot_V = SUM(phys_atm(ipa_V,:,:))
     atm(ia,:,:) = (loc_atm_tot(ia)/loc_atm_tot_V)*conv_Pa_atm*const_R_SI*atm(ia_T,:,:)
  end do

  ! *** UPDATE INTERFACE ARRAYS ***
  ! return new <atm>
  dum_sfcatm(:,:,:) = atm(:,:,:)
  ! reset integrated flux array
  dum_sfxsumatm(:,:,:) = 0.0

END SUBROUTINE atchem
! ******************************************************************************************************************************** !


! ******************************************************************************************************************************** !
! RESTART AtChem (save data)
SUBROUTINE rest_atchem()
  USE atchem_lib
  IMPLICIT NONE
  ! local variables
  integer::l
  CHARACTER(len=255)::loc_filename
  ! dump restart data
  loc_filename = TRIM(par_outdir_name)//trim(par_outfile_name)
  OPEN(unit=out,status='replace',file=loc_filename,form='unformatted',action='write')
  WRITE(unit=out)                                    &
       & n_l_atm,                                    &
       & (conv_iselected_ia(l),l=1,n_l_atm),         &
       & (atm(conv_iselected_ia(l),:,:),l=1,n_l_atm)
  close(unit=out)
END SUBROUTINE rest_atchem
! ******************************************************************************************************************************** !

