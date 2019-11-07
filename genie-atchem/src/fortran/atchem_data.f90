! ******************************************************************************************************************************** !
! atchem_data.f90
! Atmospheric Chemistry
! DATA LOADING/SAVING/INITIALIZATION ROUTINES
! ******************************************************************************************************************************** !


MODULE atchem_data

  
  USE atchem_lib
  IMPLICIT NONE
  SAVE
  
  
CONTAINS
  
  
  ! ****************************************************************************************************************************** !
  ! LOAD AtCheM 'goin' FILE OPTIONS
  SUBROUTINE sub_load_goin_atchem()
    ! local variables
    integer::l,ia                                                ! tracer counter
    integer::ios                                                 !
    ! read data_ATCHEM file
    open(unit=in,file='data_ATCHEM',status='old',action='read',iostat=ios)
    if (ios /= 0) then
       print*,'ERROR: could not open ATCHEM initialisation namelist file'
       stop
    end if
    ! read in namelist and close data_ATCHEM file
    read(UNIT=in,NML=ini_atchem_nml,IOSTAT=ios)
    if (ios /= 0) then
       print*,'ERROR: could not read ATCHEM namelist'
       stop
    else
       close(unit=in)
    end if
    ! set and report namelist data
    par_indir_name = trim(par_indir_name)//'/'
    par_outdir_name = trim(par_outdir_name)//'/'
    par_rstdir_name = trim(par_rstdir_name)//'/'
    if (ctrl_debug_init > 0) then
       ! --- TRACER INITIALIZATION ----------------------------------------------------------------------------------------------- !
       print*,'--- TRACER INITIALIZATION --------------------------'
       DO l=1,n_l_atm
          ia = conv_iselected_ia(l)
          print*,'atm tracer initial value: ',trim(string_atm(ia)),' = ',atm_init(ia)
       end do
       ! --- COSMOGENIC & RADIOGENIC PRODUCTION ---------------------------------------------------------------------------------- !
       print*,'--- COSMOGENIC & RADIOGENIC PRODUCTION -------------'
       print*,'Global cosmogenic production rate of 14C (mol yr-1) : ',par_atm_F14C
       ! --- EMISSIONS-TO-ATMOSPHERE --------------------------------------------------------------------------------------------- !
       print*,'--- EMISSIONS-TO-ATMOSPHERE ------------------------'
       print*,'Wetlands CH4 flux (mol yr-1)                        : ',par_atm_wetlands_FCH4
       print*,'Wetlands CH4 d13C (o/oo)                            : ',par_atm_wetlands_FCH4_d13C
       ! --- RUN CONTROL --------------------------------------------------------------------------------------------------------- !
       print*,'--- RUN CONTROL ------------------------------------'
       print*,'Continuing run?                                     : ',ctrl_continuing
       ! --- I/O DIRECTORY DEFINITIONS ------------------------------------------------------------------------------------------- !
       print*,'--- I/O DIRECTORY DEFINITIONS ----------------------'
       print*,'Input dir. name                                     : ',trim(par_indir_name)
       print*,'Output dir. name                                    : ',trim(par_outdir_name)
       print*,'Restart (input) dir. name                           : ',trim(par_rstdir_name)
       print*,'Filename for restart input                          : ',trim(par_infile_name)
       print*,'Filename for restart output                         : ',trim(par_outfile_name)
       ! #### INSERT CODE TO LOAD ADDITIONAL PARAMETERS ########################################################################## !
       !
       ! ######################################################################################################################### !
    end if
  END SUBROUTINE sub_load_goin_atchem
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! load AtChem restart data
  SUBROUTINE sub_load_atchem_restart()
    IMPLICIT NONE
    ! local variables
    integer::l,ios                                  ! local counting variables
    CHARACTER(len=255)::loc_filename                ! filename string
    integer::loc_n_l_atm                            ! number of selected tracers in the re-start file
    integer,DIMENSION(n_atm)::loc_conv_iselected_ia ! number of selected atmospheric tracers in restart
    ! retrieve restart data
    loc_filename = TRIM(par_rstdir_name)//trim(par_infile_name)
    OPEN(unit=in,status='old',file=loc_filename,form='unformatted',action='read',IOSTAT=ios)
    If (ios /= 0) then
       CALL sub_report_error( &
            & 'atchem_data','sub_load_atchem_restart', &
            & 'You have requested a CONTINUING run, but restart file <'//trim(loc_filename)//'> does not exist', &
            & 'SKIPPING - using default initial values (FILE: gem_config_atm.par)', &
            & (/const_real_null/),.false. &
            & )
    else
       read(unit=in) &
            & loc_n_l_atm,                                          &
            & (loc_conv_iselected_ia(l),l=1,loc_n_l_atm),           &
            & (atm(loc_conv_iselected_ia(l),:,:),l=1,loc_n_l_atm)
    end if
    close(unit=in)
  end SUBROUTINE sub_load_atchem_restart
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! initialize atmosphere grid
  SUBROUTINE sub_init_phys_atm()
    ! local variables
    INTEGER::i,j
    real::loc_th0,loc_th1,loc_s0,loc_s1,loc_ds
    real,dimension(0:n_j)::loc_s,loc_sv
    ! zero array
    phys_atm(:,:,:) = 0.0
    ! calculate local constants
    loc_th0 = -const_pi/2 
    loc_th1 = const_pi/2 
    loc_s0 = sin(loc_th0)    
    loc_s1 = sin(loc_th1)  
    loc_ds = (loc_s1-loc_s0)/real(n_j)
    DO j=0,n_j
       loc_sv(j) = loc_s0 + real(j)*loc_ds
       loc_s(j) = loc_sv(j) - 0.5*loc_ds
    end do
    ! initialize array values
    DO i=1,n_i
       DO j=1,n_j
          phys_atm(ipa_lat,i,j)  = (180.0/const_pi)*ASIN(loc_s(j))
          phys_atm(ipa_lon,i,j)  = (360.0/real(n_i))*(real(i)-0.5) + par_grid_lon_offset
          phys_atm(ipa_dlat,i,j) = (180.0/const_pi)*(ASIN(loc_sv(j)) - ASIN(loc_sv(j-1)))
          phys_atm(ipa_dlon,i,j) = (360.0/real(n_i))
          phys_atm(ipa_dh,i,j)   = par_atm_th
          phys_atm(ipa_A,i,j)    = 2.0*const_pi*(const_rEarth**2)*(1.0/real(n_i))*(loc_sv(j) - loc_sv(j-1))
          phys_atm(ipa_rA,i,j)   = 1.0/phys_atm(ipa_A,i,j)
          phys_atm(ipa_V,i,j)    = phys_atm(ipa_dh,i,j)*phys_atm(ipa_A,i,j)
          phys_atm(ipa_rV,i,j)   = 1.0/phys_atm(ipa_V,i,j)
       END DO
    END DO
  END SUBROUTINE sub_init_phys_atm
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! CONFIGURE AND INITIALIZE TRACER COMPOSITION - ATMOSPHERE
  SUBROUTINE sub_init_tracer_atm_comp()
    ! local variables
    INTEGER::i,j,ia
    real::loc_tot,loc_frac,loc_standard
    ! initialize global arrays
    atm(:,:,:)  = 0.0
    ! set <atm> array
    ! NOTE: need to seed ia_T as temperature is required in order to convert between mole (total) and partial pressure
    DO i=1,n_i
       DO j=1,n_j
          DO ia=1,n_atm
             IF (atm_select(ia)) THEN
                SELECT CASE (atm_type(ia))
                CASE (0)
                   if (ia == ia_T) atm(ia,i,j) = const_zeroC                
                CASE (1)
                   atm(ia,i,j) = atm_init(ia)
                CASE (11,12,13,14)
                   loc_tot  = atm_init(atm_dep(ia))
                   loc_standard = const_standards(atm_type(ia))
                   loc_frac = fun_calc_isotope_fraction(atm_init(ia),loc_standard)
                   atm(ia,i,j) = loc_frac*loc_tot
                END SELECT
             end if
          END DO
       END DO
    END DO
  END SUBROUTINE sub_init_tracer_atm_comp
  ! ****************************************************************************************************************************** !


END MODULE atchem_data
