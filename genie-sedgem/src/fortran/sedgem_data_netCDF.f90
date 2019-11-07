! ******************************************************************************************************************************** !
! sedgem_data_netCDF.f90
! DATA LOADING/SAVING ROUTINES
! ******************************************************************************************************************************** !


MODULE sedgem_data_netCDF

  
  use genie_control
  USE gem_cmn
  USE gem_util
  USE gem_netcdf
  USE sedgem_lib
  USE sedgem_box
  IMPLICIT NONE
  SAVE


CONTAINS


  SUBROUTINE sub_save_netcdf (dum_yr)

    real,         intent(in) :: dum_yr
    character(120) :: loc_title, loc_timunit
    real           :: loc_c0, loc_c1, loc_c10, loc_c100, loc_c500, loc_c1e3, loc_c1e6
    real           :: loc_c6e3, loc_dat
    integer        :: loc_it(6), loc_i, loc_id_time, loc_id_lonm, loc_id_misc
    integer        :: loc_id_latm, loc_id_lon_e
    integer        :: loc_id_lat_e
    real,dimension(0:n_i) :: loc_lon_e
    real,dimension(0:n_j) :: loc_lat_e
    real,dimension(n_i,n_j) :: loc_mask
    logical :: loc_defined

    loc_c0 = 0.
    loc_c1 = 1.
    loc_c10 = 10.
    loc_c100 = 100.
    loc_c500 = 500.
    loc_c1e3 = 1.e3
    loc_c6e3 = 6.e3
    loc_c1e6 = 1.e6

    !-----------------------------------------------------------------------
    !     open file and get latest record number
    !-----------------------------------------------------------------------

    loc_defined = .true.
    loc_i = 0
    if (ntrec_sout .eq. 0) then 
       loc_defined = .false.
       loc_i = 1
    end if
    !print*,'filename = ',string_ncout2d
    call sub_opennext (string_ncout2d, dum_yr, loc_i, ntrec_sout, ntrec_siou)

    if (.not. loc_defined ) then
      
     IF (ctrl_timeseries_output) THEN
       !-----------------------------------------------------------------------       
       !       open file        
       !-----------------------------------------------------------------------       
       call sub_opennew (string_ncout2d, ntrec_siou)
     ENDIF
     
       !-----------------------------------------------------------------------
       !       start definitions
       !-----------------------------------------------------------------------
       call sub_redef (ntrec_siou)

       !-----------------------------------------------------------------------
       !       set global attributes
       !-----------------------------------------------------------------------
       loc_title = '2-D surface sediment and bottom-water properties'
       write (loc_timunit,'(a,F12.2)') 'equal_month_year since 0000-01-01 00:00:00'
       call sub_putglobal (ntrec_siou, string_ncout2d, loc_title, string_ncrunid, loc_timunit)

       !-----------------------------------------------------------------------
       !       define dimensions
       !-----------------------------------------------------------------------
       call sub_defdim('time',ntrec_siou,const_integer_zero,loc_id_time)
       call sub_defdim('para',ntrec_siou,const_integer_one,loc_id_misc)
       call sub_defdim ('lon', ntrec_siou, n_i, loc_id_lonm)
       call sub_defdim ('lat', ntrec_siou, n_j, loc_id_latm)
       call sub_defdim ('lon_edges', ntrec_siou, n_i+1, loc_id_lon_e)
       call sub_defdim ('lat_edges', ntrec_siou, n_j+1, loc_id_lat_e)

       loc_it(1) = loc_id_time
       call sub_defvar ('time', ntrec_siou, 1, loc_it, loc_c0, loc_c0, 'T', 'D' &
            &, 'Year', 'time', trim(loc_timunit))
       call sub_defvar ('year', ntrec_siou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','year', ' ',' ')
       !-----------------------------------------------------------------------
       !       define 1d data (x, y or z)
       !-----------------------------------------------------------------------
       loc_it(1) = loc_id_lonm
       call sub_defvar ('lon', ntrec_siou, 1, loc_it, loc_c0, loc_c0, 'X', 'D' , &
            &'longitude of the t grid', 'longitude', 'degrees_east')
       loc_it(1) = loc_id_latm
       call sub_defvar ('lat', ntrec_siou, 1, loc_it, loc_c0, loc_c0, 'Y', 'D' , &
            &'latitude of the t grid', 'latitude', 'degrees_north')
       loc_it(1) = loc_id_lon_e
       call sub_defvar ('lon_edges', ntrec_siou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
            &'longitude of t grid edges', ' ', 'degrees')
       loc_it(1) = loc_id_lat_e
       call sub_defvar ('lat_edges', ntrec_siou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
            &'latitude of t grid edges', ' ', 'degrees')

       loc_it(1) = loc_id_time
       call sub_defvar ('month', ntrec_siou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','month', ' ',' ')
       call sub_defvar ('day', ntrec_siou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','day', ' ',' ')
       call sub_defvar ('hour', ntrec_siou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','hour', ' ',' ')
       call sub_defvar ('minute', ntrec_siou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','minute', ' ',' ')
       call sub_defvar ('second', ntrec_siou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','second', ' ',' ')

!!$       DO ips=1,6
!!$          if(mod(ips,2).eq.0) then
!!$             loc_it(1) = loc_id_lonm
!!$          else
!!$             loc_it(1) = loc_id_latm
!!$          end if
!!$          call sub_defvar (trim(string_phys_sed(ips)), ntrec_siou, 1, loc_it, loc_c0, loc_c0, &
!!$               & ' ', 'D', trim(string_phys_sed(ips)), ' ','1')
!!$       END DO
!!$       loc_it(1) = loc_id_misc
!!$       !     only write A and rA
!!$       DO ips=8,9
!!$          call sub_defvar (trim(string_phys_sed(ips)), ntrec_siou, 1, loc_it, loc_c0, loc_c0, &
!!$               & ' ', 'D', trim(string_phys_sed(ips)), ' ','1')
!!$       END DO

       !-----------------------------------------------------------------------
       !       define 2d data (x,y)
       !-----------------------------------------------------------------------
       loc_it(1) = loc_id_lonm
       loc_it(2) = loc_id_latm
       call sub_defvar ('grid_mask', ntrec_siou, 2, loc_it, loc_c0, loc_c100, ' ', 'D', &
            &'sediment mask', ' ' ,'n/a')
       call sub_defvar ('grid_topo', ntrec_siou, 2, loc_it, -loc_c6e3, loc_c0, ' ', 'D', &
            &'sediment topography', ' ' ,'m')
       !-----------------------------------------------------------------------
       !       end definitions
       !-----------------------------------------------------------------------
       call sub_enddef (ntrec_siou)

       !-----------------------------------------------------------------------
       !       write 1d data (x, y or z)
       !-----------------------------------------------------------------------

       call sub_putvars  ('time', ntrec_siou, ntrec_sout, dum_yr, loc_c1, loc_c0)
       call sub_putvarIs ('year', ntrec_siou, ntrec_sout, floor(dum_yr), loc_c1, loc_c0)
       loc_dat = 1.
       call sub_putvars  ('month', ntrec_siou, ntrec_sout, loc_dat, loc_c1, loc_c0)
       call sub_putvars  ('day', ntrec_siou, ntrec_sout, loc_dat, loc_c1, loc_c0)
       loc_dat = 0.
       call sub_putvars  ('hour', ntrec_siou, ntrec_sout, loc_dat, loc_c1, loc_c0)
       call sub_putvars  ('minute', ntrec_siou, ntrec_sout, loc_dat, loc_c1, loc_c0)
       call sub_putvars  ('second', ntrec_siou, ntrec_sout, loc_dat, loc_c1, loc_c0)
       call sub_putvar1d ('lon', ntrec_siou, n_i, ntrec_sout, n_i, phys_sed(ips_lon,:,1), loc_c1, loc_c0)
       call edge_maker (1, loc_lon_e, phys_sed(ips_lon,:,1), &
            & phys_sed(ips_lone,:,1), phys_sed(ips_dlon,:,1), n_i)
       call sub_putvar1d ('lon_edges', ntrec_siou, n_i+1, ntrec_sout, n_i+1, loc_lon_e, loc_c1, loc_c0)
       call sub_putvar1d ('lat', ntrec_siou, n_j, ntrec_sout, n_j, phys_sed(ips_lat,1,:), loc_c1, loc_c0)
       call edge_maker (1, loc_lat_e, phys_sed(ips_lat,1,:), &
            & phys_sed(ips_latn,1,:), phys_sed(ips_dlat,1,:), n_j)
       call sub_putvar1d ('lat_edges', ntrec_siou, n_j+1, ntrec_sout, n_j+1, loc_lat_e, loc_c1, loc_c0)
!!$       DO ips=1,6
!!$          if(mod(ips,2).eq.0) then
!!$             call sub_putvar1d (trim(string_phys_sed(ips)), ntrec_siou, n_i, ntrec_sout, &
!!$                  & n_i, phys_sed(ips,:,1), loc_c1, loc_c0)
!!$          else
!!$             call sub_putvar1d (trim(string_phys_sed(ips)), ntrec_siou, n_j, ntrec_sout, &
!!$                  & n_j, phys_sed(ips,1,:), loc_c1, loc_c0)
!!$          end if
!!$       END DO
!!$       DO ips=8,9
!!$          call sub_putvar1d (trim(string_phys_sed(ips)), ntrec_siou, 1, ntrec_sout, &
!!$               & 1, phys_sed(ips,1,1), loc_c1, loc_c0)
!!$       END DO
       !-----------------------------------------------------------------------
       !       write 2d data (x,y)
       !-----------------------------------------------------------------------
       loc_mask = phys_sed(ips_mask_sed,:,:)
       call sub_putvar2d ('grid_mask', ntrec_siou, n_i, n_j, ntrec_sout, &
            & phys_sed(ips_mask_sed,:,:), loc_mask)
       call sub_putvar2d ('grid_topo', ntrec_siou, n_i, n_j, ntrec_sout, &
            & -phys_sed(ips_mask_sed,:,:)*phys_sed(ips_D,:,:), loc_mask)
    else
      call sub_putvars  ('time', ntrec_siou, ntrec_sout, dum_yr, loc_c1, loc_c0)
      call sub_putvarIs ('year', ntrec_siou, ntrec_sout, floor(dum_yr), loc_c1, loc_c0)
    end if

    call sub_sync(ntrec_siou)

  END SUBROUTINE sub_save_netcdf


  SUBROUTINE sub_save_netcdf_sed2d(dum_dtyr,dum_sfcsumocn)
    ! dummy valiables
    real,INTENT(in)::dum_dtyr
    real,DIMENSION(n_ocn,n_i,n_j),intent(in)::dum_sfcsumocn
    ! local variables
    INTEGER::i,j,l,io,is,ic
    CHARACTER(len=255)::loc_unitsname
    REAL,DIMENSION(n_sed,n_i,n_j)::loc_sed_coretop
    REAL,DIMENSION(n_sed,n_i,n_j)::loc_sed_preservation
    REAL,DIMENSION(n_i,n_j)::loc_ij,loc_mask
    REAL,DIMENSION(n_i,n_j)::loc_mask_reef,loc_mask_muds       ! 
    REAL,DIMENSION(n_i,n_j)::loc_mask_dsea                     ! 
    real::loc_tot,loc_frac,loc_standard                        ! 
    real::loc_c0,loc_c1,loc_c100,loc_cr1e15,loc_cr1e18         ! 

    ! *** INITIALIZE ***
    ! initialize local variables
    loc_c0 = 0.
    loc_c1 = 1.
    loc_c100 = 100.
    loc_cr1e15 = 1.e-15
    loc_cr1e18 = 1.e-18
    ! initialize local masks
    loc_mask(:,:)      = phys_sed(ips_mask_sed,:,:)
    loc_mask_reef(:,:) = phys_sed(ips_mask_sed_reef,:,:)
    loc_mask_muds(:,:) = phys_sed(ips_mask_sed_muds,:,:)
    loc_mask_dsea(:,:) = phys_sed(ips_mask_sed,:,:)*(1.0 - loc_mask_reef(:,:))*(1.0 - loc_mask_muds(:,:))
    ! calculate core-top sediment composition data
    loc_sed_coretop(:,:,:) = fun_sed_coretop()
    ! calculate local sediment preservation
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       DO i=1,n_i
          DO j=1,n_j
             IF (sed_fsed(is,i,j) > 0.0) THEN
                loc_sed_preservation(is,i,j) = 100.0*(sed_fsed(is,i,j) - sed_fdis(is,i,j))/sed_fsed(is,i,j)
             else
                loc_sed_preservation(is,i,j) = 0.0
             end if
          end do
       end do
    end do

    ! SAVE *ALL* BENTHIC DATA
    ! ocean interface tracer data field
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             SELECT CASE (ocn_type(io))
             CASE (0)
                loc_ij(i,j) = dum_sfcsumocn(io,i,j)
                If (io == io_T) loc_unitsname = 'K'
                If (io == io_S) loc_unitsname = 'o/oo'
             CASE (1)
                loc_ij(i,j) = dum_sfcsumocn(io,i,j)
                loc_unitsname = 'mol kg-1'
             case (11:20)
                loc_tot  = dum_sfcsumocn(ocn_dep(io),i,j)
                loc_frac = dum_sfcsumocn(io,i,j)
                loc_standard = const_standards(ocn_type(io))
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
                loc_unitsname = 'o/oo'
             end SELECT
          end do
       end do
       call sub_adddef_netcdf(ntrec_siou,3,'ocn_'//trim(string_ocn(io)), &
            & 'overlying ocean tracer properties - '//trim(string_ocn(io)),trim(loc_unitsname),ocn_mima(io,1),ocn_mima(io,2))
       call sub_putvar2d('ocn_'//trim(string_ocn(io)),ntrec_siou,n_i, n_j,ntrec_sout,loc_ij(:,:),loc_mask)
    END DO
    ! carbonate chemistry data field
    DO ic=1,n_carb
       SELECT CASE (ocn_type(io))
       CASE (ic_H,ic_ohm_cal,ic_ohm_arg)
          loc_unitsname = 'n/a'
       case default
          loc_unitsname = 'mol kg-1'
       end SELECT
       call sub_adddef_netcdf(ntrec_siou,3,'carb_'//trim(string_carb(ic)), &
            & 'overlying ocean carbonate chemistry - '//trim(string_carb(ic)),trim(loc_unitsname),loc_c0,loc_c0)
       call sub_putvar2d('carb_'//trim(string_carb(ic)),ntrec_siou,n_i,n_j,ntrec_sout,sed_carb(ic,:,:),loc_mask)
    END DO

    ! SAVE *ALL* SEDIMENT FLUXES
    ! interface flux data
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             SELECT CASE (sed_type(is))
             CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det)
                loc_ij(i,j) = sed_fsed(is,i,j)/dum_dtyr
                loc_unitsname = 'mol cm-2 yr-1'
             case (11:20)
                loc_tot  = sed_fsed(sed_dep(is),i,j)
                loc_frac = sed_fsed(is,i,j)
                loc_standard = const_standards(sed_type(is))
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
                loc_unitsname = 'o/oo'
             END SELECT
          end do
       end do
       SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & 11:20)
          call sub_adddef_netcdf(ntrec_siou,3,'fsed_'//trim(string_sed(is)), &
               & 'sediment rain flux - '//trim(string_sed(is)), trim(loc_unitsname),loc_c0,loc_c0)
          call sub_putvar2d('fsed_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
       END SELECT
    END DO
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             SELECT CASE (sed_type(is))
             CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det)
                loc_ij(i,j) = sed_fdis(is,i,j)/dum_dtyr
                loc_unitsname = 'mol cm-2 yr-1'
             case (11:20)
                loc_tot  = sed_fsed(sed_dep(is),i,j)
                loc_frac = sed_fsed(is,i,j)
                loc_standard = const_standards(sed_type(is))
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.)
                loc_unitsname = 'o/oo'
             END SELECT
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & 11:20)
          call sub_adddef_netcdf(ntrec_siou,3,'fdis_'//trim(string_sed(is)), &
               & 'sediment dissolution flux - '//trim(string_sed(is)),trim(loc_unitsname),loc_c0,loc_c0)
          call sub_putvar2d('fdis_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
       END SELECT
    END DO
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             SELECT CASE (sed_type(is))
             CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det)
                loc_ij(i,j) = loc_sed_preservation(is,i,j)
                loc_unitsname = '%'
             END SELECT
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det)
          call sub_adddef_netcdf(ntrec_siou,3,'fpres_'//trim(string_sed(is)), &
               & 'sediment burial (preservation) flux - '//trim(string_sed(is)),trim(loc_unitsname),loc_c0,loc_c0)
          call sub_putvar2d('fpres_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
       end SELECT
    END DO
    ! save interface flux data - dust (log10)
    IF (sed_select(is_det)) THEN
       loc_ij(:,:) = const_real_zero
       loc_unitsname = 'mol cm-2 yr-1'
       ! log10 data
       DO i=1,n_i
          DO j=1,n_j
             IF (sed_fsed(is_det,i,j) > 0.0) THEN
                loc_ij(i,j) = log10(sed_fsed(is_det,i,j)/dum_dtyr)
             else
                loc_ij(i,j) = const_real_null
             end if
          end do
       end do
       call sub_adddef_netcdf(ntrec_siou,3,'fsed_'//trim(string_sed(is_det))//'_log10', &
            & 'sediment rain flux - detrital material (log10)',trim(loc_unitsname),loc_c0,loc_c0)
       call sub_putvar2d('fsed_'//trim(string_sed(is_det))//'_log10',ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
    END IF
    ! CaCO3:POC 'rain ratio'
    IF (sed_select(is_CaCO3) .AND. sed_select(is_POC)) THEN
       loc_unitsname = 'n/a'
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             if (sed_fsed(is_POC,i,j) > const_real_nullsmall) then
                loc_ij(i,j) = sed_fsed(is_CaCO3,i,j)/sed_fsed(is_POC,i,j)
             end if
          END DO
       END DO
       call sub_adddef_netcdf(ntrec_siou,3,'fsed_CaCO3toPOC', &
            & 'sediment rain flux - CaCO3 to POC rain ratio',trim(loc_unitsname),loc_c0,loc_c0)
       call sub_putvar2d('fsed_CaCO3toPOC',ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
    end if

    ! SAVE DEEP-SEA SEDIMENT DATA
    ! core-top data
    DO is=1,n_sed
       IF (sed_select(is)) THEN
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det)
             loc_unitsname = 'wt%'
          CASE (11:20)
             loc_unitsname = 'o/oo'
          CASE (par_sed_type_age)
             loc_unitsname = 'years'
          END SELECT
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio,par_sed_type_age, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & 11:20)
             call sub_adddef_netcdf(ntrec_siou,3,'sed_'//trim(string_sed(is)), &
                  & 'surface sediment composition - '//trim(string_sed(is)),trim(loc_unitsname),sed_mima(is,1),sed_mima(is,2))
             call sub_putvar2d('sed_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_sed_coretop(is,:,:),loc_mask_dsea)
          end SELECT
       END IF
    END DO

    ! SAVE REEF SEDIMENTS
    ! bulk composition
    is = is_CaCO3
    loc_unitsname = 'wt%'
    loc_ij(:,:) = loc_sed_coretop(is,:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'reef_sed_'//trim(string_sed(is)), &
         & 'reef composition - '//trim(string_sed(is)),trim(loc_unitsname),sed_mima(is,1),sed_mima(is,2))
    call sub_putvar2d('reef_sed_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask_reef)
    is = is_CaCO3_13C
    loc_unitsname = 'o/oo'
    loc_ij(:,:) = loc_sed_coretop(is,:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'reef_sed_'//trim(string_sed(is)), &
         & 'reef composition - '//trim(string_sed(is)),trim(loc_unitsname),sed_mima(is,1),sed_mima(is,2))
    call sub_putvar2d('reef_sed_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask_reef)
    ! precipitation
    is = is_CaCO3
    loc_unitsname = 'mol cm-2 yr-1'
    loc_ij(:,:) = sed_fsed(is,:,:)/dum_dtyr
    call sub_adddef_netcdf(ntrec_siou,3,'reef_fsed_'//trim(string_sed(is)), &
         & 'reef precipitation rate - '//trim(string_sed(is)),trim(loc_unitsname),loc_c0,loc_c0)
    call sub_putvar2d('reef_fsed_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask_reef)
    is = is_CaCO3_13C
    loc_unitsname = 'o/oo'
    DO i=1,n_i
       DO j=1,n_j
          loc_tot  = sed_fsed(sed_dep(is),i,j)/dum_dtyr
          loc_frac = sed_fsed(is,i,j)/dum_dtyr
          loc_standard = const_standards(sed_type(is))
          loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
       end do
    end do
    call sub_adddef_netcdf(ntrec_siou,3,'reef_fsed_'//trim(string_sed(is)), &
         & 'reef precipitation rate - '//trim(string_sed(is)),trim(loc_unitsname),sed_mima(is,1),sed_mima(is,2))
    call sub_putvar2d('reef_fsed_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask_reef)

    ! SAVE SHALLOW WATER SEDIMENTS
    ! bulk composition
    is = is_POC
    loc_unitsname = 'wt%'
    loc_ij(:,:) = loc_sed_coretop(is,:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'muds_sed_'//trim(string_sed(is)), &
         & 'shallow water sediment composition - '//trim(string_sed(is)),trim(loc_unitsname),sed_mima(is,1),sed_mima(is,2))
    call sub_putvar2d('muds_sed_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_sed_coretop(is,:,:),loc_mask_muds)
    is = is_POC_13C
    loc_unitsname = 'o/oo'
    loc_ij(:,:) = loc_sed_coretop(is,:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'muds_sed_'//trim(string_sed(is)), &
         & 'shallow water sediment composition - '//trim(string_sed(is)),trim(loc_unitsname),sed_mima(is,1),sed_mima(is,2))
    call sub_putvar2d('muds_sed_'//trim(string_sed(is)),ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(i,j),loc_mask_muds)

    ! MISCELLANEOUS FIELDS
    ! 
    loc_unitsname = 'cm3 cm-3'
    loc_ij(:,:) = phys_sed(ips_poros,:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'misc_sed_porosity', &
         & 'sediment surface porosity',trim(loc_unitsname),loc_c0,loc_c0)
    call sub_putvar2d('misc_sed_porosity',ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
    ! 
    loc_unitsname = 'cm2 yr-1'
    loc_ij(:,:) = phys_sed(ips_mix_k0,:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'misc_sed_max_k', &
         & 'maximum (surface) sediment bioturbation mixing rate',trim(loc_unitsname),loc_c0,loc_c0)
    call sub_putvar2d('misc_sed_max_k',ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
    ! 
    loc_unitsname = 'n/a'
    loc_ij(:,:) = loc_mask_dsea(:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'grid_mask_dsea','deep sea sediments mask',trim(loc_unitsname),loc_c0,loc_c0)
    call sub_putvar2d('grid_mask_dsea',ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
    ! 
    loc_unitsname = 'n/a'
    loc_ij(:,:) = loc_mask_reef(:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'grid_mask_reef','reef area mask',trim(loc_unitsname),loc_c0,loc_c0)
    call sub_putvar2d('grid_mask_reef',ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)
    ! 
    loc_unitsname = 'n/a'
    loc_ij(:,:) = loc_mask_muds(:,:)
    call sub_adddef_netcdf(ntrec_siou,3,'grid_mask_muds','shallow sediments mask',trim(loc_unitsname),loc_c0,loc_c0)
    call sub_putvar2d('grid_mask_muds',ntrec_siou,n_i,n_j,ntrec_sout,loc_ij(:,:),loc_mask)

  end SUBROUTINE sub_save_netcdf_sed2d


SUBROUTINE sub_save_netcdf_sed3d(dum_yr)
  ! dummy valiables
  real, intent(in)  :: dum_yr
  ! local variables
  integer::i,j,o,is, loc_extra, l
  real::loc_tot,loc_frac,loc_standard
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:)::loc_sed_save               ! hold reordered data saving array
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_age_cal         ! sediment age data saving array
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_age_ash         ! sediment age data saving array
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_ash_norm        ! sediment age data saving array
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_age_14C         ! sediment age data saving array
  REAL,ALLOCATABLE,DIMENSION(:,:,:)::loc_sed_save_CaCO3_D14C      ! sediment CaCO3 D14C data saving array
  REAL::loc_sed_tot_wt                                            ! total mass of solid coponents
  REAL::loc_sed_tot_vol                                           ! total volume of solid coponents
  INTEGER::loc_ash_max_o                                          ! running ash volume maximum sub-layer number
  REAL::loc_ash_max                                               ! running ash volume maximum
  REAL::loc_ash_max_depth                                         ! running ash volume maximum down-core depth
  REAL::loc_ash_conv_dbs_age                                      ! convert depth to age using ash stratigraphy
  real::loc_ash_tot
  INTEGER,DIMENSION(n_i,n_j)::loc_n_sed_stack_top ! sediment stack top layer number
  REAL,DIMENSION(n_i,n_j)::loc_sed_stack_top_th   ! sediment stack top layer thickness

  ! *** initialize variables ***
  loc_extra = 6
  ! allocate array for holding sediment data reordered for writing to file
  ALLOCATE(loc_sed_save(n_sed+loc_extra,n_i,n_j,0:n_sed_tot),STAT=alloc_error) 
  ALLOCATE(loc_sed_save_age_cal(n_i,n_j,0:n_sed_tot),STAT=alloc_error)      
  ALLOCATE(loc_sed_save_age_ash(n_i,n_j,0:n_sed_tot),STAT=alloc_error)       
  ALLOCATE(loc_sed_save_ash_norm(n_i,n_j,0:n_sed_tot),STAT=alloc_error)     
  ALLOCATE(loc_sed_save_age_14C(n_i,n_j,0:n_sed_tot),STAT=alloc_error)     
  ALLOCATE(loc_sed_save_CaCO3_D14C(n_i,n_j,0:n_sed_tot),STAT=alloc_error)  
  ! check for problems allocating array space
  IF (error /= 0) THEN
     CALL sub_report_error( &
          & 'sedgem_data','sub_save_netcdf_sed3d', &
          & 'Array space could not be allocated', &
          & 'STOPPING', &
          & (/const_real_null/),.true. &
          & )
  ENDIF
  ! zero local variables
  loc_sed_save(:,:,:,:)          = const_real_zero
  loc_sed_save_age_cal(:,:,:)    = const_real_zero
  loc_sed_save_age_ash(:,:,:)    = const_real_zero
  loc_sed_save_ash_norm(:,:,:)   = const_real_zero
  loc_sed_save_age_14C(:,:,:)    = const_real_zero
  loc_sed_save_CaCO3_D14C(:,:,:) = const_real_zero

  ! *** transform sediment array for saving to file ***

  DO i = 1,n_i
     DO j = 1,n_j

        ! *** (a) calculate local constants
        loc_n_sed_stack_top(i,j)  = INT(sed_top_h(i,j)) + 1
        loc_sed_stack_top_th(i,j) = sed_top_h(i,j) - REAL((loc_n_sed_stack_top(i,j) - 1))

        ! *** (b) copy core top layer data
        loc_sed_save(:,i,j,0) = sed_top(:,i,j)

        ! *** (c) copy core stack sub-layer data
        DO o = loc_n_sed_stack_top(i,j),1,-1
           loc_sed_save(:,i,j,(loc_n_sed_stack_top(i,j) - o + 1)) = sed(:,i,j,o)
        END DO

        ! *** (d) calculate carbonate internal age
        DO o = 0,n_sed_tot
           IF (loc_sed_save(is_CaCO3,i,j,o) > const_real_nullsmall) THEN
              loc_sed_save_age_cal(i,j,o) = loc_sed_save(is_CaCO3_age,i,j,o)/loc_sed_save(is_CaCO3,i,j,o)
           ELSE
              loc_sed_save_age_cal(i,j,o) = 0.0
           ENDIF
        ENDDO

        ! *** (e) normailze solid sediment components to a mass fraction basis (if required)
        DO o = 0,n_sed_tot
           IF (ctrl_data_save_wtfrac) THEN
              loc_sed_tot_wt = fun_calc_sed_mass(loc_sed_save(:,i,j,o))
              IF (loc_sed_tot_wt > const_real_nullsmall) THEN
                 loc_sed_save(:,i,j,o) = conv_sed_cm3_g(:)*loc_sed_save(:,i,j,o)/loc_sed_tot_wt
              end IF
           else
              loc_sed_tot_vol = fun_calc_sed_vol(loc_sed_save(:,i,j,o))
              IF (loc_sed_tot_vol > const_real_nullsmall) THEN
                 loc_sed_save(:,i,j,o) = loc_sed_save(:,i,j,o)/loc_sed_tot_vol
              end IF
           end if
        END DO

        ! *** (f) normalize ash content wt% (or vol%)
        loc_ash_tot = &
             & par_sed_top_th*loc_sed_save(is_ash,i,j,0) + &
             & loc_sed_stack_top_th(i,j)*loc_sed_save(is_ash,i,j,1)+ &
             & sum(loc_sed_save(is_ash,i,j,2:n_sed_tot))
        if (loc_ash_tot > const_real_nullsmall) then
           loc_sed_save_ash_norm(i,j,:) = loc_sed_save(is_ash,i,j,:)/loc_ash_tot
        else
           loc_sed_save_ash_norm(i,j,:) = const_real_zero
        end if

        ! *** (g) produce stratigraphic marker age scale
        loc_ash_max = 0.0
        DO o = 0,n_sed_tot
           IF (loc_sed_save(is_ash,i,j,o) > loc_ash_max + const_real_nullsmall) THEN
              loc_ash_max   = loc_sed_save(is_ash,i,j,o)
              loc_ash_max_o = o
           ENDIF
        END DO
        ! calculate ash maximum depth
        SELECT CASE (loc_ash_max_o)
        CASE (0)
           loc_ash_max_depth = par_sed_top_th/2.0
        CASE (1)
           loc_ash_max_depth = par_sed_top_th + loc_sed_stack_top_th(i,j)/2.0
        CASE default
           loc_ash_max_depth = par_sed_top_th + loc_sed_stack_top_th(i,j) + REAL((loc_ash_max_o - 2)) + 0.5
        END SELECT
        ! calculate linear age-depth relation
        loc_ash_conv_dbs_age = 1.0/loc_ash_max_depth
        ! generate age scale
        o = 0
        loc_sed_save_age_ash(i,j,o) = loc_ash_conv_dbs_age*(par_sed_top_th/2.0)
        o = 1
        loc_sed_save_age_ash(i,j,o) = loc_ash_conv_dbs_age*(par_sed_top_th + loc_sed_stack_top_th(i,j)/2.0)
        DO o = 2,n_sed_tot
           loc_sed_save_age_ash(i,j,o) = loc_ash_conv_dbs_age*(par_sed_top_th + loc_sed_stack_top_th(i,j) + REAL((o - 2)) + 0.5)
        END DO

        ! *** (h) convert mass or volume fraction to wt% or vol%, respectively
        DO l=1,n_l_sed
           is = conv_iselected_is(l)
           SELECT CASE (sed_type(is))
           case (par_sed_type_bio,par_sed_type_abio)
              DO o = 0,n_sed_tot
                 loc_sed_save(is,i,j,o) = 100.0*loc_sed_save(is,i,j,o)
              end DO
           end SELECT
        end do
        
        ! *** (i) calculate isotope per mils
        DO l=1,n_l_sed
           is = conv_iselected_is(l)
           SELECT CASE (sed_type(is))
           case (11:20)
              DO o = 0,n_sed_tot
                 loc_tot  = loc_sed_save(sed_dep(is),i,j,o)
                 loc_frac = loc_sed_save(is,i,j,o)
                 loc_standard = const_standards(sed_type(is))
                 loc_sed_save(is,i,j,o) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.)
              end DO
           end SELECT
        end do

        ! *** (j) calculate D14C and radiocarbon age
        loc_sed_save_CaCO3_D14C(i,j,:) = const_real_zero
        loc_sed_save_age_14C(i,j,:) = const_real_zero
        if (sed_select(is_CaCO3_14C)) then
           DO o = 0,n_sed_tot
              IF (loc_sed_save(is_CaCO3,i,j,o) > const_real_nullsmall) THEN
                 loc_sed_save_CaCO3_D14C(i,j,o) = &
                      & fun_convert_delta14CtoD14C(loc_sed_save(is_CaCO3_13C,i,j,o),loc_sed_save(is_CaCO3_14C,i,j,o))
                 loc_sed_save_age_14C(i,j,o) = &
                      & fun_convert_D14Ctoage(loc_sed_save_CaCO3_D14C(i,j,o))
              end if
           end DO
        end if

     END DO
  END DO

  ! *** save prescribed sediment (code) location data ***
  DO i = 1,n_i
     DO j = 1,n_j
        o = 0
        loc_sed_save(n_sed+1,i,j,o) = par_sed_top_th/2.0
        loc_sed_save(n_sed+2,i,j,o) = loc_sed_save_age_cal(i,j,o)
        loc_sed_save(n_sed+3,i,j,o) = loc_sed_save_age_ash(i,j,o)
        loc_sed_save(n_sed+4,i,j,o) = loc_sed_save_age_14C(i,j,o)
        loc_sed_save(n_sed+5,i,j,o) = loc_sed_save_CaCO3_D14C(i,j,o)
        loc_sed_save(n_sed+6,i,j,o) = loc_sed_save_ash_norm(i,j,o)

        o = 1
        loc_sed_save(n_sed+1,i,j,o) = par_sed_top_th + loc_sed_stack_top_th(i,j)/2.0
        loc_sed_save(n_sed+2,i,j,o) = loc_sed_save_age_cal(i,j,o)
        loc_sed_save(n_sed+3,i,j,o) = loc_sed_save_age_ash(i,j,o)
        loc_sed_save(n_sed+4,i,j,o) = loc_sed_save_age_14C(i,j,o)
        loc_sed_save(n_sed+5,i,j,o) = loc_sed_save_CaCO3_D14C(i,j,o)
        loc_sed_save(n_sed+6,i,j,o) = loc_sed_save_ash_norm(i,j,o)

        do o=2,n_sed_tot
           loc_sed_save(n_sed+1,i,j,o) = par_sed_top_th + loc_sed_stack_top_th(i,j) + REAL((o - 2)) + 0.5
           loc_sed_save(n_sed+2,i,j,o) = loc_sed_save_age_cal(i,j,o)
           loc_sed_save(n_sed+3,i,j,o) = loc_sed_save_age_ash(i,j,o)
           loc_sed_save(n_sed+4,i,j,o) = loc_sed_save_age_14C(i,j,o)
           loc_sed_save(n_sed+5,i,j,o) = loc_sed_save_CaCO3_D14C(i,j,o)
           loc_sed_save(n_sed+6,i,j,o) = loc_sed_save_ash_norm(i,j,o)
        end do
     end DO
  end DO

  call sub_save_netcdf_core(dum_yr,loc_sed_save, loc_extra)

  ! deallocate local arrays
  DEALLOCATE(loc_sed_save,STAT=dealloc_error)
  DEALLOCATE(loc_sed_save_age_cal,STAT=dealloc_error)
  DEALLOCATE(loc_sed_save_age_ash,STAT=dealloc_error)
  DEALLOCATE(loc_sed_save_ash_norm,STAT=dealloc_error)
  DEALLOCATE(loc_sed_save_age_14C,STAT=dealloc_error)
  DEALLOCATE(loc_sed_save_CaCO3_D14C,STAT=dealloc_error)
  ! check for problems de-allocating array space
  IF (dealloc_error /= 0) THEN
     CALL sub_report_error( &
          & 'sedgem_data','sub_save_netcdf_sed3d', &
          & 'Array space could not be deallocated', &
          & 'STOPPING', &
          & (/const_real_null/),.true. &
          & )
  ENDIF

end SUBROUTINE sub_save_netcdf_sed3d


SUBROUTINE sub_save_netcdf_full(dum_dtyr, dum_yr, dum_coretop, dum_sfcsumocn)
  ! dummy arguments
  REAL,INTENT(in)   :: dum_dtyr, dum_yr
  real,DIMENSION(n_sed,n_i,n_j),intent(in)::dum_coretop
  real,DIMENSION(n_ocn,n_i,n_j),intent(in)::dum_sfcsumocn
  ! local variables
  integer        :: loc_ntrec, loc_iou, loc_it(3), loc_id_time
  integer        :: loc_id_lonm, loc_id_latm, loc_id_lon_e, loc_id_lat_e
  character(120) :: loc_title, loc_timunit
  real :: loc_c0, loc_c1, loc_c100, loc_c1e6, loc_dat, loc_cr1e15, loc_cr1e18
  real, DIMENSION(n_i,n_j):: loc_mask, loc_tmp
  real,dimension(0:n_i) :: loc_lon_e
  real,dimension(0:n_j) :: loc_lat_e

  loc_c0 = 0.
  loc_c1 = 1.
  loc_c100 = 100.
  loc_c1e6 = 1.e+6
  loc_cr1e15 = 1.e-15
  loc_cr1e18 = 1.e-18

  !-----------------------------------------------------------------------
  !     open file and get latest record number
  !-----------------------------------------------------------------------

  call sub_opennext (string_nctop, dum_yr, 0, loc_ntrec, loc_iou)

  if ( loc_ntrec .le. 1 ) then

     !-----------------------------------------------------------------------
     !       start definitions
     !-----------------------------------------------------------------------
     call sub_redef (loc_iou)

     !-----------------------------------------------------------------------
     !       set global attributes
     !-----------------------------------------------------------------------
     loc_title = 'Sediment diagnostics'
     write (loc_timunit,'(a,F12.2)') 'equal_month_year since 0000-01-01 00:00:00'
     call sub_putglobal (loc_iou, string_nctop, loc_title, string_ncrunid, loc_timunit)

     !-----------------------------------------------------------------------
     !       define dimensions
     !-----------------------------------------------------------------------
     call sub_defdim ('time', loc_iou, 0, loc_id_time)
     call sub_defdim ('lon', loc_iou, n_i, loc_id_lonm)
     call sub_defdim ('lat', loc_iou, n_j, loc_id_latm)
     call sub_defdim ('lon_edges', loc_iou, n_i+1, loc_id_lon_e)
     call sub_defdim ('lat_edges', loc_iou, n_j+1, loc_id_lat_e)

     !-----------------------------------------------------------------------
     !       define 1d data (t)
     !-----------------------------------------------------------------------
     loc_it(1) = loc_id_time
     call sub_defvar ('time', loc_iou, 1, loc_it, loc_c0, loc_c0, 'T', 'D' &
          &, 'Year', 'time', trim(loc_timunit))
     loc_it(1) = loc_id_lonm
     call sub_defvar ('lon', loc_iou, 1, loc_it, loc_c0, loc_c0, 'X', 'D' , &
          &'longitude of the t grid', 'longitude', 'degrees_east')
     loc_it(1) = loc_id_latm
     call sub_defvar ('lat', loc_iou, 1, loc_it, loc_c0, loc_c0, 'Y', 'D' , &
          &'latitude of the t grid', 'latitude', 'degrees_north')
     loc_it(1) = loc_id_lon_e
     call sub_defvar ('lon_edges', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
          &'longitude of t grid edges', ' ', 'degrees')
     loc_it(1) = loc_id_lat_e
     call sub_defvar ('lat_edges', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
          &'latitude of t grid edges', ' ', 'degrees')

     loc_it(1) = loc_id_time
     call sub_defvar ('year', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','year', ' ',' ')
     call sub_defvar ('month', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','month', ' ',' ')
     call sub_defvar ('day', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','day', ' ',' ')
     call sub_defvar ('hour', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','hour', ' ',' ')
     call sub_defvar ('minute', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','minute', ' ',' ')
     call sub_defvar ('second', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','second', ' ',' ')

     loc_it(1) = loc_id_lonm
     loc_it(2) = loc_id_latm
     loc_it(3) = loc_id_time
     call sub_defvar ('D_m', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'ocean depth', ' ' ,'m')
     call sub_defvar ('T_K', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'temperature', ' ' ,'K')
     call sub_defvar ('S_mil', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'salinity', ' ' ,'psu')
     call sub_defvar ('CO2_uM', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'total dissolved carbon (DIC)', ' ' ,'mu mol')
     call sub_defvar ('ALK_uM', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'alkalinity (ALK)', ' ' ,'mu mol')
     call sub_defvar ('O2_uM', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'dissolved oxygen concentration', ' ' ,'mu mol')
     call sub_defvar ('Ca_uM', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'dissolved calcium concentration', ' ' ,'mu mol')
     call sub_defvar ('CO3_uM',loc_iou,3,loc_it,loc_c0,loc_c0,' ','F','carbonate ion (CO32-) concentration',' ','mu mol')
     call sub_defvar ('ohm', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'saturation state (calcite)', ' ' ,'mu mol')
     call sub_defvar ('dCO3_uM', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'deviation of CO3 from saturation', ' ' ,'mu mol')
     call sub_defvar ('POC', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment POC content', ' ' ,'%')
     call sub_defvar ('cal', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment CaCO3 content', ' ' ,'%')
     call sub_defvar ('opal', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment opal content', ' ' ,'%')
     call sub_defvar ('fsedPOC', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment POC rain flux', ' ' ,'1')
     call sub_defvar ('fsedcal', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment CaCO3 rain flux', ' ' ,'1')
     call sub_defvar ('fsedopal', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment opal rain flux', ' ' ,'1')
     call sub_defvar ('fdisPOC', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment POC dissolution flux', ' ' ,'1')
     call sub_defvar ('fdisopal', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment opal dissolution flux', ' ' ,'1')
     call sub_defvar ('fdiscal', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', 'sediment CaCO3 dissolution flux', ' ' ,'1')

     !-----------------------------------------------------------------------
     !       end definitions
     !-----------------------------------------------------------------------
     call sub_enddef (loc_iou)

     if (loc_ntrec .eq. 0) loc_ntrec = 1

     call sub_putvar1d ('lon', loc_iou, n_i, loc_ntrec, n_i, phys_sed(ips_lon,:,1), loc_c0, loc_c0)
     call edge_maker (1, loc_lon_e, phys_sed(ips_lon,:,1), &
          & phys_sed(ips_lone,:,1), phys_sed(ips_dlon,:,1), n_i)
     call sub_putvar1d ('lon_edges', loc_iou, n_i+1, loc_ntrec, n_i+1, loc_lon_e, loc_c0, loc_c0)
     call sub_putvar1d ('lat', loc_iou, n_j, loc_ntrec, n_j, phys_sed(ips_lat,1,:), loc_c0, loc_c0)
     call edge_maker (1, loc_lat_e, phys_sed(ips_lat,1,:), &
          & phys_sed(ips_latn,1,:), phys_sed(ips_dlat,1,:), n_j)
     call sub_putvar1d ('lat_edges', loc_iou, n_j+1, loc_ntrec, n_j+1, loc_lat_e, loc_c0, loc_c0)       

  endif

  call sub_putvars ('time', loc_iou, loc_ntrec, dum_yr, loc_c0, loc_c0)
  call sub_putvarIs ('year', loc_iou, loc_ntrec, floor(dum_yr), loc_c1, loc_c0)
  loc_dat = 1.
  call sub_putvars  ('month', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
  call sub_putvars  ('day', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
  loc_dat = 0.
  call sub_putvars  ('hour', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
  call sub_putvars  ('minute', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
  call sub_putvars  ('second', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)

  loc_mask= 0.
  where ( sed_mask )
     loc_mask = 1.
  endwhere

  loc_tmp = phys_sed(ips_D,:,:)
  call sub_putvar2d ('D_m', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = dum_sfcsumocn(io_T,:,:)
  call sub_putvar2d ('T_K', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = dum_sfcsumocn(io_S,:,:)
  call sub_putvar2d ('S_mil', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*dum_sfcsumocn(io_DIC,:,:)
  call sub_putvar2d ('DIC_uM', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*dum_sfcsumocn(io_ALK,:,:)
  call sub_putvar2d ('ALK_uM', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*dum_sfcsumocn(io_O2,:,:)
  call sub_putvar2d ('O2_uM', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*dum_sfcsumocn(io_Ca,:,:)
  call sub_putvar2d ('Ca_uM', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*sed_carb(ic_conc_CO3,:,:)
  call sub_putvar2d ('CO3_uM', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = sed_carb(ic_ohm_cal,:,:)
  call sub_putvar2d ('ohm', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*sed_carb(ic_dCO3_cal,:,:)
  call sub_putvar2d ('dCO3_uM', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c100*dum_coretop(is_POC,:,:)
  call sub_putvar2d ('POC', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c100*dum_coretop(is_CaCO3,:,:)
  call sub_putvar2d ('cal', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c100*dum_coretop(is_opal,:,:)
  call sub_putvar2d ('opal', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*sed_fsed(is_POC,:,:)/dum_dtyr
  call sub_putvar2d ('fsedPOC', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*sed_fsed(is_CaCO3,:,:)/dum_dtyr
  call sub_putvar2d ('fsedcal', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*sed_fsed(is_opal,:,:)/dum_dtyr 
  call sub_putvar2d ('fsedopal', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*sed_fdis(is_POC,:,:)/dum_dtyr
  call sub_putvar2d ('fdisPOC', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*sed_fdis(is_CaCO3,:,:)/dum_dtyr
  call sub_putvar2d ('fdiscal', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_c1e6*sed_fdis(is_opal,:,:)/dum_dtyr 
  call sub_putvar2d ('fdisopal', loc_iou, n_i, n_j, loc_ntrec, loc_tmp, loc_mask)

  call sub_closefile(loc_iou)

END SUBROUTINE sub_save_netcdf_full



SUBROUTINE sub_save_netcdf_core(dum_yr, dum_sed, dum_extra)
  ! dummy arguments
  real,   intent(in) :: dum_yr
  integer,intent(in) :: dum_extra
  real, dimension(n_sed+dum_extra,n_i,n_j,0:n_sed_tot),intent(in):: dum_sed
  ! local variables
  integer        :: is, n, l
  integer        :: loc_ntrec, loc_iou, loc_it(4), loc_id_time, loc_id_zt
  integer        :: loc_id_lonm, loc_id_latm, loc_id_lon_e, loc_id_lat_e
  character(120) :: loc_title, loc_timunit
  real :: loc_c0, loc_c1, loc_c100, loc_c1e6, loc_dat, loc_cr1e15, loc_cr1e18
  real, dimension(0:n_i) :: loc_lon_e
  real, dimension(0:n_j) :: loc_lat_e
  real, dimension(n_i,n_j,0:n_sed_tot):: loc_mask, loc_tmp
  real, dimension(n_sed+dum_extra,n_i,n_j,0:n_sed_tot):: loc_sed
  !
  loc_c0 = 0.
  loc_c1 = 1.
  loc_c100 = 100.
  loc_c1e6 = 1.e+6
  loc_cr1e15 = 1.e-15
  loc_cr1e18 = 1.e-18

  !-----------------------------------------------------------------------
  !     open file and get latest record number
  !-----------------------------------------------------------------------
  call sub_opennext (string_nccore, dum_yr, 0, loc_ntrec, loc_iou)

  if ( loc_ntrec .le. 1 ) then

     !-----------------------------------------------------------------------
     !       start definitions
     !-----------------------------------------------------------------------
     call sub_redef (loc_iou)

     !-----------------------------------------------------------------------
     !       set global attributes
     !-----------------------------------------------------------------------
     loc_title = '3-D global sediment stack volume'
     write (loc_timunit,'(a,F12.2)') 'equal_month_year since 0000-01-01 00:00:00'
     call sub_putglobal (loc_iou, string_nccore, loc_title, string_ncrunid, loc_timunit)

     !-----------------------------------------------------------------------
     !       define dimensions
     !-----------------------------------------------------------------------
     call sub_defdim ('time', loc_iou, 0, loc_id_time)
     call sub_defdim ('lon', loc_iou, n_i, loc_id_lonm)
     call sub_defdim ('lat', loc_iou, n_j, loc_id_latm)
     call sub_defdim ('lon_edges', loc_iou, n_i+1, loc_id_lon_e)
     call sub_defdim ('lat_edges', loc_iou, n_j+1, loc_id_lat_e)
     call sub_defdim ('level', loc_iou, n_sed_tot+1, loc_id_zt)

     !-----------------------------------------------------------------------
     !       define 1d data (t)
     !-----------------------------------------------------------------------
     loc_it(1) = loc_id_time
     call sub_defvar ('time', loc_iou, 1, loc_it, loc_c0, loc_c0, 'T', 'D' &
          &, 'time (zero by default)', 'time', trim(loc_timunit))
     loc_it(1) = loc_id_lonm
     call sub_defvar ('lon', loc_iou, 1, loc_it, loc_c0, loc_c0, 'X', 'D' , &
          &'longitude of the t grid', 'longitude', 'degrees_east')
     loc_it(1) = loc_id_latm
     call sub_defvar ('lat', loc_iou, 1, loc_it, loc_c0, loc_c0, 'Y', 'D' , &
          &'latitude of the t grid', 'latitude', 'degrees_north')
     loc_it(1) = loc_id_lon_e
     call sub_defvar ('lon_edges', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
          &'longitude of t grid edges', ' ', 'degrees')
     loc_it(1) = loc_id_lat_e
     call sub_defvar ('lat_edges', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
          &'latitude of t grid edges', ' ', 'degrees')
     loc_it(1) = loc_id_zt
     call sub_defvar ('level', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'D' , &
          &'depth levels', ' ', '1')
     loc_it(1) = loc_id_time
     call sub_defvar ('year', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','year', ' ',' ')
     call sub_defvar ('month', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','month', ' ',' ')
     call sub_defvar ('day', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','day', ' ',' ')
     call sub_defvar ('hour', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','hour', ' ',' ')
     call sub_defvar ('minute', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','minute', ' ',' ')
     call sub_defvar ('second', loc_iou, 1, loc_it, loc_c0, loc_c0, ' ', 'F','second', ' ',' ')

     loc_it(1) = loc_id_lonm
     loc_it(2) = loc_id_latm
     loc_it(3) = loc_id_zt

     call sub_defvar ('depth', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', &
          & 'sediment depth below surface', ' ', 'cm')
     call sub_defvar ('CaCO3_age', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', &
          & 'age of bulk carbonate', ' ', 'yr')
     call sub_defvar ('ash_age', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', &
          & 'age scale from ash peak', ' ', 'yr')
     call sub_defvar ('14C_age', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', &
          & 'radiocarbon age of bulk carbonate', ' ', 'yr')
     call sub_defvar ('D14C', loc_iou, 3, loc_it, loc_c0, loc_c0, ' ', 'F', &
          & 'D14C of bulk carbonate', ' ', 'per mil')
     DO l=1,n_l_sed
        is = conv_iselected_is(l)
        call sub_defvar (trim(string_sed(is)), loc_iou, 3, &
             & loc_it, loc_c0, loc_c0, ' ', 'F', &
             & string_longname_sed(is), ' sediment tracer - '//trim(string_sed(is)), ' ')
     end do

     !-----------------------------------------------------------------------
     !       end definitions
     !-----------------------------------------------------------------------
     call sub_enddef (loc_iou)
     call sub_sync(loc_iou)

     if (loc_ntrec .eq. 0) loc_ntrec = 1

     call sub_putvar1d ('lon', loc_iou, n_i, loc_ntrec, n_i, &
          & phys_sed(ips_lon,:,1), loc_c1, loc_c0)
     call edge_maker (1, loc_lon_e, phys_sed(ips_lon,:,1), &
          & phys_sed(ips_lone,:,1), phys_sed(ips_dlon,:,1), n_i)
     call sub_putvar1d ('lon_edges', loc_iou, n_i+1, loc_ntrec, &
          & n_i+1, loc_lon_e, loc_c1, loc_c0)
     call sub_putvar1d ('lat', loc_iou, n_j, loc_ntrec, n_j, &
          & phys_sed(ips_lat,1,:), loc_c1, loc_c0)
     call edge_maker (1, loc_lat_e, phys_sed(ips_lat,1,:), &
          & phys_sed(ips_latn,1,:), phys_sed(ips_dlat,1,:), n_j)
     call sub_putvar1d ('lat_edges', loc_iou, n_j+1, loc_ntrec, &
          & n_j+1, loc_lat_e, loc_c1, loc_c0)       

  endif

  !-----------------------------------------------------------------------
  !      write everything there is to know about the slime
  !-----------------------------------------------------------------------

    call sub_putvars ('time', loc_iou, loc_ntrec, dum_yr, loc_c1, loc_c0)
    call sub_putvarIs ('year', loc_iou, loc_ntrec, floor(dum_yr), loc_c1, loc_c0)
    loc_dat=1.
    call sub_putvars  ('month', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
    call sub_putvars  ('day', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
    loc_dat = 0.
    call sub_putvars  ('hour', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
    call sub_putvars  ('minute', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)
    call sub_putvars  ('second', loc_iou, loc_ntrec, loc_dat, loc_c1, loc_c0)

  loc_sed = dum_sed
  where ( dum_sed .eq. -1000 )
     loc_sed = 0.0
  endwhere

  loc_mask = 0.
  do n = 0,n_sed_tot
     where (sed_mask)
        loc_mask(:,:,n) = 1.
     endwhere
  enddo
  loc_tmp = loc_sed(n_sed+1,:,:,:)
  call sub_putvar3d ('depth', loc_iou, n_i, n_j, n_sed_tot+1,  &
       & loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_sed(n_sed+2,:,:,:)
  call sub_putvar3d ('CaCO3_age', loc_iou, n_i, n_j, n_sed_tot+1,  &
       & loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_sed(n_sed+3,:,:,:)
  call sub_putvar3d ('ash_age', loc_iou, n_i, n_j, n_sed_tot+1,  &
       & loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_sed(n_sed+4,:,:,:)
  call sub_putvar3d ('14C_age', loc_iou, n_i, n_j, n_sed_tot+1,  &
       & loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_sed(n_sed+5,:,:,:)
  call sub_putvar3d ('D14C', loc_iou, n_i, n_j, n_sed_tot+1,  &
       & loc_ntrec, loc_tmp, loc_mask)
  loc_tmp = loc_sed(n_sed+6,:,:,:)
  DO l=1,n_l_sed
     is = conv_iselected_is(l)
     loc_tmp(:,:,:)= loc_sed(is,:,:,:)
     call sub_putvar3d (trim(string_sed(is)), loc_iou, n_i,  &
          & n_j, n_sed_tot+1, loc_ntrec, loc_tmp(:,:,:), loc_mask)
  end do

  call sub_closefile(loc_iou)

END SUBROUTINE sub_save_netcdf_core


END MODULE sedgem_data_netCDF
