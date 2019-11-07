
! ******************************************************************************************************************************** !
SUBROUTINE initialise_sedgem( &
     & dum_genie_timestep,    &
     & dum_sfxsumsed,         &
     & dum_sfcsumocn,         &
     & dum_sfcsed,            &
     & dum_sfxocn)
  USE sedgem_lib
  USE sedgem_data
  USE genie_util, ONLY: check_iostat
  IMPLICIT NONE
  ! dummy arguments
  REAL,INTENT(in)::dum_genie_timestep                                   ! genie time-step (in seconds)
  real,dimension(n_sed,n_i,n_j),intent(inout)::dum_sfxsumsed
  real,dimension(n_ocn,n_i,n_j),intent(inout)::dum_sfcsumocn
  real,dimension(n_sed,n_i,n_j),intent(inout)::dum_sfcsed
  real,dimension(n_ocn,n_i,n_j),intent(inout)::dum_sfxocn
  ! local variables
  integer::i,j

  print*,'======================================================='
  print*,' >>> Initialising SEDGEM sediment geochem. module ...'

  ! *** load goin information ***
  call sub_load_goin_sedgem()

  ! *** initialize interface arrays ***
  ! (i.e., set all elements to zero)
  dum_sfxsumsed(:,:,:) = 0.0   ! 
  dum_sfcsumocn(:,:,:) = 0.0   ! 
  dum_sfcsed(:,:,:)    = 0.0   ! 
  dum_sfxocn(:,:,:)    = 0.0   ! 

  ! *** dimension the size of the 2-D sediment grid arrays ***
  ! NOTE: check for problems allocating array space
  ALLOCATE(phys_sed(n_phys_sed,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_mask(n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_mask_reef(n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_mask_muds(n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed(n_sed,n_i,n_j,n_sed_tot),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_top(n_sed,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_top_h(n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_top_INTdth(n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_fsed(n_sed,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_fdis(n_sed,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sedocn_fnet(n_ocn,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_carb(n_carb,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_carbconst(n_carbconst,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_carbalk(n_carbalk,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_carbisor(n_carbisor,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_save_mask(n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_fsed_OLD(n_sed,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)
  ALLOCATE(sed_fdis_OLD(n_sed,n_i,n_j),STAT=alloc_error)
  call check_iostat(alloc_error,__LINE__,__FILE__)

  ! *** INITIALIZE SedGeM ***
  ! initialize dynamically-allocated arrays (those not done elsewhere)
  sed_fsed(:,:,:)      = 0.0   !
  sed_fdis(:,:,:)      = 0.0   !
  sedocn_fnet(:,:,:)   = 0.0   !
  sed_carb(:,:,:)      = 0.0   !
  sed_carbconst(:,:,:) = 0.0   !
  ! setup SedGeM grid
  CALL sub_init_phys_sed()
  ! set meta-options and verify self-consistency of chosen parameters
  call sub_check_par_sedgem()
  ! initialize sediment sub-system
  call sub_init_sed()
  call sub_init_sed_layers_default()
  ! initialize sediment core location data save mask
  call sub_init_sedgem_save_sed_data()
  ! initialize sediment core environmental conditions saving
  call sub_sedgem_init_sedcoresenv()
  ! seed the aqueous carbonate system with an initial value of [H+]
  sed_carb(ic_H,:,:) = 1.0E-8
  ! initialize bioturbation profile
  IF (ctrl_sed_bioturb) THEN
     if (ctrl_sed_bioturb_Archer) then
        ALLOCATE(par_sed_mix_k(0:par_n_sed_mix),STAT=error)
     else
        ! load bioturbation profile data
        call sub_load_sed_mix_k()
     end if
  end IF
  ! seed the integrated ocean->sediment flux interface array with an ash pulse (1 mol m-2 should do)
  ! NOTE: <dum_sfxsumsed> in units of (mol m-2)
  DO i=1,n_i
     DO j=1,n_j
        IF (sed_mask(i,j)) dum_sfxsumsed(is_ash,i,j) = 1.0
     end DO
  end DO

  ! *** INITIALIZE netCDF OUTPUT ***
  string_nctsglob  = TRIM(par_outdir_name)//'ts_sedgem_glob.nc'
  string_nctop     = TRIM(par_outdir_name)//'toplayer_sedgem.nc'
  string_ncout2d   = TRIM(par_outdir_name)//'fields_sedgem_2d.nc'
  string_nccore    = TRIM(par_outdir_name)//'fields_sedgem_3d.nc'
  IF (ctrl_timeseries_output) THEN
     IF (ctrl_continuing.AND.ctrl_append_data) THEN
        OPEN(unit=in,status='old',file=TRIM(par_rstdir_name)//'netcdf_record_numbers',form='formatted',action='read')
        READ(unit=in,fmt='(i6)') ntrec_sout
        close(unit=in)
     ELSE
        ntrec_sout = 0
     ENDIF
     print*, 'netcdf record number: ',ntrec_sout
     print*,'par_outdir_name = par_rstdir_name:',par_outdir_name.eq.par_rstdir_name
  ELSE
     ntrec_sout = 0
  ENDIF

  ! *** INITIALIZE DATA SAVING ***
  IF (ctrl_timeseries_output) THEN
     ! initialize timestep counter
     tstep_count = 0
     tsteps_per_year = conv_yr_s/(dum_genie_timestep*kocn_loop*conv_kocn_ksedgem)
     PRINT*,'timesteps per year                                  :',tsteps_per_year
     ! load in years for output generation
     CALL sub_data_output_years()
     year = min(output_years_0d(output_counter_0d),output_years_2d(output_counter_2d))
  ENDIF

  ! LOAD SEDIMENT RE-START ***
  ! NOTE: modify sediment ages if a continuing run
  IF (ctrl_continuing) then
     call sub_load_sedgem_restart()
     sed(is_CaCO3_age,:,:,:)    = sed(is_CaCO3_age,:,:,:)    + par_misc_t_runtime*sed(is_CaCO3,:,:,:)
     sed_top(is_CaCO3_age,:,:)  = sed_top(is_CaCO3_age,:,:)  + par_misc_t_runtime*sed_top(is_CaCO3,:,:)
     ! update sediment interface composition data arrays with restart 
     dum_sfcsed(:,:,:) = fun_sed_coretop()
     ! set ash event (to mark run start)
     par_sed_ashevent = .true.
  else
     par_sed_ashevent = .false.
  end if

  print*,' <<< Initialisation complete'
  print*,'======================================================='

end SUBROUTINE initialise_sedgem
! ******************************************************************************************************************************** !
