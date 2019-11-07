! ******************************************************************************************************************************** !
! END SEDGEM
SUBROUTINE end_sedgem( &
     & dum_dts,        &
     & dum_sfcsumocn   &
     & )
  use genie_control
  USE gem_cmn
  USE sedgem_lib
  USE sedgem_data
  USE sedgem_data_netCDF
  USE genie_util, ONLY: check_iostat
  IMPLICIT NONE
  ! dummy arguments
  REAL,INTENT(in)::dum_dts
  real,DIMENSION(n_ocn,n_i,n_j),intent(in)::dum_sfcsumocn
  ! local variables
  real::loc_dtyr ! local time step in years
  real::loc_dts  ! local time step in seconds

  print*,'======================================================='
  print*,' >>> Initialising SEDGEM module shutdown ...'

  ! *** SAVE DATA ***
  IF (ctrl_timeseries_output) THEN
     ! data saved during runtime for years specified in output file
     ! (DEFAULT IS genie-sedgem/data/input/sedgem_output_years.dat) 
     !IF (ctrl_debug_lvl2) PRINT*,'saving netcdf record numbers',ntrec_sout
     PRINT*,'saving netcdf record number',ntrec_sout
     OPEN(unit=out,status='replace',file=TRIM(par_outdir_name)//'netcdf_record_numbers',form='formatted',action='write')
     WRITE(unit=out,fmt='(i6)') ntrec_sout                             
     close(unit=out)
  ELSE
     ! calculate sediment model time step length
     loc_dts  = dum_dts
     loc_dtyr = loc_dts/conv_yr_s
     !save diagnostics
     call sub_data_save_seddiag_GLOBAL(loc_dtyr,dum_sfcsumocn)
     ! save requested sediment cores as ASCII
     call sub_sedgem_save_sedcore()
     !save final oecan-sediment interface properties
     ! NOTE: netCDF was set up to be able to save multiple time-slices, but is only being used here to save a single final slice
     if (ctrl_data_save_ascii) call sub_data_save_seddiag_2D(loc_dtyr,dum_sfcsumocn)
     call sub_save_netcdf(const_real_zero)
     call sub_save_netcdf_sed2d(loc_dtyr,dum_sfcsumocn)
     call sub_closefile(ntrec_siou)
     ntrec_sout = ntrec_sout + 1
  ENDIF
  
  ! *** clean up ***
  ! deallocate local arrays
  DEALLOCATE(phys_sed,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_mask,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_mask_reef,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_mask_muds,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_top,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_top_h,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_fsed,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_fdis,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sedocn_fnet,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_carb,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_carbconst,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_carbalk,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_carbisor,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_save_mask,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_fsed_OLD,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)
  DEALLOCATE(sed_fdis_OLD,STAT=dealloc_error)
  call check_iostat(dealloc_error,__LINE__,__FILE__)

  print*,' <<< Shutdown complete'
  print*,'======================================================='

END SUBROUTINE end_sedgem
! ******************************************************************************************************************************** !
