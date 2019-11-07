! 
! File:          ocean_goldstein_Impl.F90
! Symbol:        ocean.goldstein-v0.1.1
! Symbol Type:   class
! Babel Version: 0.11.0
! sidl Created:  20061123 17:47:08 GMT
! Generated:     20061123 17:47:20 GMT
! Description:   Server-side implementation for ocean.goldstein
! 
! WARNING: Automatically generated; only changes within splicers preserved
! 
! babel-version = 0.11.0
! source-line   = 106
! source-url    = file:/homes/sp1003/genie-babel/genie-goldstein/babel/goldstein.sidl
! 


! 
! Symbol "ocean.goldstein" (version 0.1.1)
! 


#include "ocean_goldstein_interface_fAbbrev.h"
#include "sidl_ClassInfo_fAbbrev.h"
#include "ocean_goldstein_fAbbrev.h"
#include "sidl_BaseInterface_fAbbrev.h"
#include "sidl_RuntimeException_fAbbrev.h"
#include "sidl_BaseException_fAbbrev.h"
#include "sidl_BaseClass_fAbbrev.h"
#include "sidl_double_fAbbrev.h"
#include "sidl_int_fAbbrev.h"
! DO-NOT-DELETE splicer.begin(_miscellaneous_code_start)
! Insert-Code-Here {_miscellaneous_code_start} (extra code)
! DO-NOT-DELETE splicer.end(_miscellaneous_code_start)




! 
! Method:  _ctor[]
! Class constructor called when the class is created.
! 

recursive subroutine ocean_goldstein__ctor_mi(self, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use ocean_goldstein
  use ocean_goldstein_impl
  ! DO-NOT-DELETE splicer.begin(ocean.goldstein._ctor.use)
  ! Insert-Code-Here {ocean.goldstein._ctor.use} (use statements)
  ! DO-NOT-DELETE splicer.end(ocean.goldstein._ctor.use)
  implicit none
  type(ocean_goldstein_t) :: self ! in
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(ocean.goldstein._ctor)
! Insert-Code-Here {ocean.goldstein._ctor} (_ctor method)
! DO-NOT-DELETE splicer.end(ocean.goldstein._ctor)
end subroutine ocean_goldstein__ctor_mi


! 
! Method:  _dtor[]
! Class destructor called when the class is deleted.
! 

recursive subroutine ocean_goldstein__dtor_mi(self, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use ocean_goldstein
  use ocean_goldstein_impl
  ! DO-NOT-DELETE splicer.begin(ocean.goldstein._dtor.use)
  ! Insert-Code-Here {ocean.goldstein._dtor.use} (use statements)
  ! DO-NOT-DELETE splicer.end(ocean.goldstein._dtor.use)
  implicit none
  type(ocean_goldstein_t) :: self ! in
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(ocean.goldstein._dtor)
! Insert-Code-Here {ocean.goldstein._dtor} (_dtor method)
! DO-NOT-DELETE splicer.end(ocean.goldstein._dtor)
end subroutine ocean_goldstein__dtor_mi


! 
! Method:  _load[]
! Static class initializer called exactly once before any user-defined method is dispatched
! 

recursive subroutine ocean_goldstein__load_mi(exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use ocean_goldstein
  use ocean_goldstein_impl
  ! DO-NOT-DELETE splicer.begin(ocean.goldstein._load.use)
  ! Insert-Code-Here {ocean.goldstein._load.use} (use statements)
  ! DO-NOT-DELETE splicer.end(ocean.goldstein._load.use)
  implicit none
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(ocean.goldstein._load)
! Insert-Code-Here {ocean.goldstein._load} (_load method)
! DO-NOT-DELETE splicer.end(ocean.goldstein._load)
end subroutine ocean_goldstein__load_mi


! 
! Method:  initialise[]
! 

recursive subroutine ocean_goldstein_initialise_mi(self, alon1_ocn_bbl,        &
  alat1_ocn_bbl, alon2_ocn_bbl, alat2_ocn_bbl, alon3_ocn_bbl, alat3_ocn_bbl,   &
  aboxedge1_lon_ocn_bbl, aboxedge1_lat_ocn_bbl, aboxedge2_lon_ocn_bbl,         &
  aboxedge2_lat_ocn_bbl, aboxedge3_lon_ocn_bbl, aboxedge3_lat_ocn_bbl,         &
  depth1_ocn_bbl, depth2_ocn_bbl, ilandmask1_ocn_bbl, ilandmask2_ocn_bbl,      &
  ilandmask3_ocn_bbl, koverall_total_bbl, tstar_ocn_bbl, sstar_ocn_bbl,        &
  ustar_ocn_bbl, vstar_ocn_bbl, albedo_ocn_bbl, ias_out_bbl, iaf_out_bbl,      &
  ips_out_bbl, ipf_out_bbl, jsf_out_bbl, lrestart_genie_bbl, go_npstp_bbl,     &
  go_iwstp_bbl, go_itstp_bbl, go_ianav_bbl, go_saln0_bbl, go_rhoair_bbl,       &
  go_cd_bbl, go_ds_bbl, go_dphi_bbl, go_ips_bbl, go_ipf_bbl, bg_usc_bbl,       &
  go_rsc_bbl, go_tsc_bbl, go_dsc_bbl, go_fsc_bbl, go_gsc_bbl, go_rh0sc_bbl,    &
  go_rhosc_bbl, go_cpsc_bbl, go_k1_bbl, go_lmax_bbl, go_dz_bbl, go_dza_bbl,    &
  go_ias_bbl, go_iaf_bbl, go_jsf_bbl, go_c_bbl, go_cv_bbl, go_s_bbl,           &
  go_sv_bbl, go_ts_bbl, go_ts1_bbl, go_lin_bbl, go_lout_bbl, retval,           &
  exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use ocean_goldstein
  use sidl_double_array
  use sidl_int_array
  use ocean_goldstein_impl
  ! DO-NOT-DELETE splicer.begin(ocean.goldstein.initialise.use)
  ! Insert-Code-Here {ocean.goldstein.initialise.use} (use statements)
  use genie_control
  use castings
 
  ! DO-NOT-DELETE splicer.end(ocean.goldstein.initialise.use)
  implicit none
  type(ocean_goldstein_t) :: self ! in
  type(sidl_double_1d) :: alon1_ocn_bbl ! out
  type(sidl_double_1d) :: alat1_ocn_bbl ! out
  type(sidl_double_1d) :: alon2_ocn_bbl ! out
  type(sidl_double_1d) :: alat2_ocn_bbl ! out
  type(sidl_double_1d) :: alon3_ocn_bbl ! out
  type(sidl_double_1d) :: alat3_ocn_bbl ! out
  type(sidl_double_1d) :: aboxedge1_lon_ocn_bbl ! out
  type(sidl_double_1d) :: aboxedge1_lat_ocn_bbl ! out
  type(sidl_double_1d) :: aboxedge2_lon_ocn_bbl ! out
  type(sidl_double_1d) :: aboxedge2_lat_ocn_bbl ! out
  type(sidl_double_1d) :: aboxedge3_lon_ocn_bbl ! out
  type(sidl_double_1d) :: aboxedge3_lat_ocn_bbl ! out
  type(sidl_double_1d) :: depth1_ocn_bbl ! out
  type(sidl_double_1d) :: depth2_ocn_bbl ! out
  type(sidl_int_2d) :: ilandmask1_ocn_bbl ! out
  type(sidl_int_2d) :: ilandmask2_ocn_bbl ! out
  type(sidl_int_2d) :: ilandmask3_ocn_bbl ! out
  integer (kind=sidl_int) :: koverall_total_bbl ! in
  type(sidl_double_2d) :: tstar_ocn_bbl ! out
  type(sidl_double_2d) :: sstar_ocn_bbl ! out
  type(sidl_double_2d) :: ustar_ocn_bbl ! out
  type(sidl_double_2d) :: vstar_ocn_bbl ! out
  type(sidl_double_2d) :: albedo_ocn_bbl ! out
  type(sidl_int_1d) :: ias_out_bbl ! out
  type(sidl_int_1d) :: iaf_out_bbl ! out
  type(sidl_int_1d) :: ips_out_bbl ! out
  type(sidl_int_1d) :: ipf_out_bbl ! out
  integer (kind=sidl_int) :: jsf_out_bbl ! out
  logical :: lrestart_genie_bbl ! in
  integer (kind=sidl_int) :: go_npstp_bbl ! out
  integer (kind=sidl_int) :: go_iwstp_bbl ! out
  integer (kind=sidl_int) :: go_itstp_bbl ! out
  integer (kind=sidl_int) :: go_ianav_bbl ! out
  real (kind=sidl_double) :: go_saln0_bbl ! out
  real (kind=sidl_double) :: go_rhoair_bbl ! out
  real (kind=sidl_double) :: go_cd_bbl ! out
  type(sidl_double_1d) :: go_ds_bbl ! out
  real (kind=sidl_double) :: go_dphi_bbl ! out
  type(sidl_int_1d) :: go_ips_bbl ! out
  type(sidl_int_1d) :: go_ipf_bbl ! out
  real (kind=sidl_double) :: bg_usc_bbl ! out
  real (kind=sidl_double) :: go_rsc_bbl ! out
  real (kind=sidl_double) :: go_tsc_bbl ! out
  real (kind=sidl_double) :: go_dsc_bbl ! out
  real (kind=sidl_double) :: go_fsc_bbl ! out
  real (kind=sidl_double) :: go_gsc_bbl ! out
  real (kind=sidl_double) :: go_rh0sc_bbl ! out
  real (kind=sidl_double) :: go_rhosc_bbl ! out
  real (kind=sidl_double) :: go_cpsc_bbl ! out
  type(sidl_int_2d) :: go_k1_bbl ! out
  integer (kind=sidl_int) :: go_lmax_bbl ! out
  type(sidl_double_1d) :: go_dz_bbl ! out
  type(sidl_double_1d) :: go_dza_bbl ! out
  type(sidl_int_1d) :: go_ias_bbl ! out
  type(sidl_int_1d) :: go_iaf_bbl ! out
  integer (kind=sidl_int) :: go_jsf_bbl ! out
  type(sidl_double_1d) :: go_c_bbl ! out
  type(sidl_double_1d) :: go_cv_bbl ! out
  type(sidl_double_1d) :: go_s_bbl ! out
  type(sidl_double_1d) :: go_sv_bbl ! out
  type(sidl_double_4d) :: go_ts_bbl ! out
  type(sidl_double_4d) :: go_ts1_bbl ! out
  character (len=*) :: go_lin_bbl ! out
  character (len=*) :: go_lout_bbl ! out
  integer (kind=sidl_int) :: retval ! out
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(ocean.goldstein.initialise)
! Insert-Code-Here {ocean.goldstein.initialise} (initialise method)

  real, dimension(:), ALLOCATABLE :: alon1_ocn
  real, dimension(:), ALLOCATABLE :: alat1_ocn
  real, dimension(:), ALLOCATABLE :: alon2_ocn
  real, dimension(:), ALLOCATABLE :: alat2_ocn
  real, dimension(:), ALLOCATABLE :: alon3_ocn
  real, dimension(:), ALLOCATABLE :: alat3_ocn
  real, dimension(:), ALLOCATABLE :: aboxedge1_lon_ocn
  real, dimension(:), ALLOCATABLE :: aboxedge1_lat_ocn
  real, dimension(:), ALLOCATABLE :: aboxedge2_lon_ocn
  real, dimension(:), ALLOCATABLE :: aboxedge2_lat_ocn
  real, dimension(:), ALLOCATABLE :: aboxedge3_lon_ocn
  real, dimension(:), ALLOCATABLE :: aboxedge3_lat_ocn
  real, dimension(:), ALLOCATABLE :: depth1_ocn
  real, dimension(:), ALLOCATABLE :: depth2_ocn
  integer, dimension(:,:), ALLOCATABLE :: ilandmask1_ocn
  integer, dimension(:,:), ALLOCATABLE :: ilandmask2_ocn
  integer, dimension(:,:), ALLOCATABLE :: ilandmask3_ocn
  integer koverall_total_int
  real, dimension(:,:), ALLOCATABLE :: tstar_ocn
  real, dimension(:,:), ALLOCATABLE :: sstar_ocn
  real, dimension(:,:), ALLOCATABLE :: ustar_ocn
  real, dimension(:,:), ALLOCATABLE :: vstar_ocn
  real, dimension(:,:), ALLOCATABLE :: albedo_ocn
  
  integer, dimension(:), ALLOCATABLE :: ias_out
  integer, dimension(:), ALLOCATABLE :: iaf_out
  integer, dimension(:), ALLOCATABLE :: ips_out
  integer, dimension(:), ALLOCATABLE :: ipf_out
  integer jsf_out
  logical lrestart_genie
    
  integer go_npstp
  integer go_iwstp
  integer go_itstp
  integer go_ianav
  real go_saln0
  real go_rhoair
  real go_cd
!  real, dimension(:), ALLOCATABLE :: go_ds
  real, dimension(ilat1_ocn) :: go_ds
  real go_dphi
!  integer, dimension(:), ALLOCATABLE :: go_ips
  integer,dimension(ilat1_ocn) :: go_ips
!  integer, dimension(:), ALLOCATABLE :: go_ipf
  integer,dimension(ilat1_ocn) :: go_ipf
  real bg_usc
  real go_rsc
  real go_tsc
  real go_dsc
  real go_fsc
  real go_gsc
  real go_rh0sc
  real go_rhosc
  real go_cpsc
!  integer, dimension(:,:), ALLOCATABLE :: go_k1
  integer,dimension(ilon1_ocn,ilat1_ocn) :: go_k1
  integer go_lmax
!  go_dz
  real,dimension(1:inl1_ocn) :: go_dz
!  go_dza
  real,dimension(1:inl1_ocn) :: go_dza
!  integer, dimension(:), ALLOCATABLE :: go_ias
  integer,dimension(ilat1_ocn) :: go_ias
!  integer, dimension(:), ALLOCATABLE :: go_iaf
  integer,dimension(ilat1_ocn) :: go_iaf
  integer go_jsf
!  go_c
  real,dimension(0:ilat1_ocn) :: go_c
!  go_cv
  real,dimension(0:ilat1_ocn) :: go_cv
!  go_s
  real,dimension(0:ilat1_ocn) :: go_s
!  go_sv
  real,dimension(0:ilat1_ocn) :: go_sv
!  go_ts
  real,dimension(1:intrac_ocn,1:ilon1_ocn,1:ilat1_ocn,1:inl1_ocn)::go_ts
!  go_ts1
  real,dimension(1:intrac_ocn,1:ilon1_ocn,1:ilat1_ocn,1:inl1_ocn)::go_ts1
!  go_lin
  character::go_lin*13
!  go_lout
  character::go_lout*3
  
    allocate(alon1_ocn(ilon1_ocn))
    alon1_ocn(:) = 0.0
    allocate(alat1_ocn(ilat1_ocn))
    alat1_ocn(:) = 0.0
    allocate(alon2_ocn(ilon2_ocn))
    alon2_ocn(:) = 0.0
    allocate(alat2_ocn(ilat2_ocn))
    alat2_ocn(:) = 0.0
    allocate(alon3_ocn(ilon3_ocn))
    alon3_ocn(:) = 0.0
    allocate(alat3_ocn(ilat3_ocn))
    alat3_ocn(:) = 0.0
    allocate(aboxedge1_lon_ocn(ilon1_ocn+1))
    aboxedge1_lon_ocn(:) = 0.0
    allocate(aboxedge1_lat_ocn(ilat1_ocn+1))
    aboxedge1_lat_ocn(:) = 0.0
    allocate(aboxedge2_lon_ocn(ilon2_ocn+1))
    aboxedge2_lon_ocn(:) = 0.0
    allocate(aboxedge2_lat_ocn(ilat2_ocn+1))
    aboxedge2_lat_ocn(:) = 0.0
    allocate(aboxedge3_lon_ocn(ilon3_ocn+1))
    aboxedge3_lon_ocn(:) = 0.0
    allocate(aboxedge3_lat_ocn(ilat3_ocn+1))
    aboxedge3_lat_ocn(:) = 0.0
    
    allocate(depth1_ocn(inl1_ocn))
    allocate(depth2_ocn(inl2_ocn))
    
    allocate(ilandmask1_ocn(ilon1_ocn,ilat1_ocn))
    ilandmask1_ocn(:, :) = 0.0
    allocate(ilandmask2_ocn(ilon2_ocn,ilat2_ocn))
    ilandmask2_ocn(:, :) = 0.0
    allocate(ilandmask3_ocn(ilon3_ocn,ilat3_ocn))
    ilandmask3_ocn(:, :) = 0.0

    koverall_total_int = koverall_total_bbl
    
    allocate(tstar_ocn(ilon1_ocn,ilat1_ocn))
    tstar_ocn(:, :) = 0.0
    allocate(sstar_ocn(ilon1_ocn,ilat1_ocn))
    sstar_ocn(:, :) = 0.0
    allocate(ustar_ocn(ilon1_ocn,ilat1_ocn))
    ustar_ocn(:, :) = 0.0
    allocate(vstar_ocn(ilon1_ocn,ilat1_ocn))
    vstar_ocn(:, :) = 0.0
    allocate(albedo_ocn(ilon1_ocn,ilat1_ocn))
    albedo_ocn(:, :) = 0.0


    allocate(ias_out(ilat1_ocn))
    ias_out(:) = 0
    allocate(iaf_out(ilat1_ocn))
    iaf_out(:) = 0
    allocate(ips_out(ilat1_ocn))
    ips_out(:) = 0
    allocate(ipf_out(ilat1_ocn))
    ipf_out(:) = 0


    lrestart_genie = lrestart_genie_bbl

    
    call initialise_goldstein(alon1_ocn, alat1_ocn, alon2_ocn, alat2_ocn, &
           alon3_ocn, alat3_ocn, aboxedge1_lon_ocn, aboxedge1_lat_ocn, &
           aboxedge2_lon_ocn, aboxedge2_lat_ocn, &
           aboxedge3_lon_ocn, aboxedge3_lat_ocn, &
           depth1_ocn, depth2_ocn, &
           ilandmask1_ocn, ilandmask2_ocn, &
           ilandmask3_ocn, koverall_total_int, &
           tstar_ocn, sstar_ocn, ustar_ocn, vstar_ocn, albedo_ocn, &
           ias_out, iaf_out, ips_out, ipf_out, jsf_out, &
           lrestart_genie, go_npstp, &
           go_iwstp, go_itstp, go_ianav, &
           go_saln0, go_rhoair, &
           go_cd, go_ds, go_dphi, &
           go_ips, go_ipf, bg_usc, &
           go_rsc, go_tsc, go_dsc, &
           go_fsc, go_gsc, go_rh0sc, &
           go_rhosc, go_cpsc, go_k1, &
           go_lmax, go_dz, go_dza, &
           go_ias, go_iaf, go_jsf, &
           go_c, go_cv, go_s, &
           go_sv, go_ts, go_ts1,&
           go_lin, go_lout)
      
    call create1d(ilon1_ocn, alon1_ocn_bbl)
    call create1d(ilat1_ocn, alat1_ocn_bbl)
    call create1d(ilon2_ocn, alon2_ocn_bbl)
    call create1d(ilat2_ocn, alat2_ocn_bbl)
    call create1d(ilon3_ocn, alon3_ocn_bbl)
    call create1d(ilat3_ocn, alat3_ocn_bbl)
    
    call create1d(ilon1_ocn+1, aboxedge1_lon_ocn_bbl)
    call create1d(ilat1_ocn+1, aboxedge1_lat_ocn_bbl)
    call create1d(ilon2_ocn+1, aboxedge2_lon_ocn_bbl)
    call create1d(ilat2_ocn+1, aboxedge2_lat_ocn_bbl)
    call create1d(ilon3_ocn+1, aboxedge3_lon_ocn_bbl)
    call create1d(ilat3_ocn+1, aboxedge3_lat_ocn_bbl)

    call create1d(inl1_ocn, depth1_ocn_bbl)
    call create1d(inl2_ocn, depth2_ocn_bbl)
    
    call create2dcol(ilon1_ocn, ilat1_ocn, ilandmask1_ocn_bbl)
    call create2dcol(ilon2_ocn, ilon2_ocn, ilandmask2_ocn_bbl)
    call create2dcol(ilon3_ocn, ilon3_ocn, ilandmask3_ocn_bbl)

    call create2dcol(ilon1_ocn, ilat1_ocn, tstar_ocn_bbl)
    call create2dcol(ilon1_ocn, ilat1_ocn, sstar_ocn_bbl)
    call create2dcol(ilon1_ocn, ilat1_ocn, ustar_ocn_bbl)
    call create2dcol(ilon1_ocn, ilat1_ocn, vstar_ocn_bbl)
    call create2dcol(ilon1_ocn, ilat1_ocn, albedo_ocn_bbl)
    
   
    call create1d(ilat1_ocn, ias_out_bbl)
    call create1d(ilat1_ocn, iaf_out_bbl)
    call create1d(ilat1_ocn, ips_out_bbl)
    call create1d(ilat1_ocn, ipf_out_bbl)

    call create1d(ilat1_ocn, go_ds_bbl)
    call create1d(ilat1_ocn, go_ips_bbl)
    call create1d(ilat1_ocn, go_ipf_bbl)
    
    call create2dcol(ilon1_ocn,ilat1_ocn, go_k1_bbl)
    
    call create1d(inl1_ocn, go_dz_bbl)
    call create1d(inl1_ocn, go_dza_bbl)

    call create1d(ilat1_ocn, go_ias_bbl)
    call create1d(ilat1_ocn, go_iaf_bbl)

    call createCol([0], [ilat1_ocn], go_c_bbl)
    call createCol([0], [ilat1_ocn], go_cv_bbl)
    call createCol([0], [ilat1_ocn], go_s_bbl)
    call createCol([0], [ilat1_ocn], go_sv_bbl)
        
    call createCol([0, 0, 0, 0], [intrac_ocn-1, ilon1_ocn-1, ilat1_ocn-1, inl1_ocn-1], go_ts_bbl)
    call createCol([0, 0, 0, 0], [intrac_ocn-1, ilon1_ocn-1, ilat1_ocn-1, inl1_ocn-1], go_ts1_bbl)
    
    call castSimple1DtoSidlDouble1D(alon1_ocn, alon1_ocn_bbl)
    call castSimple1DtoSidlDouble1D(alat1_ocn, alat1_ocn_bbl)
    call castSimple1DtoSidlDouble1D(alon2_ocn, alon2_ocn_bbl)
    call castSimple1DtoSidlDouble1D(alat2_ocn, alat2_ocn_bbl)
    call castSimple1DtoSidlDouble1D(alon3_ocn, alon3_ocn_bbl)
    call castSimple1DtoSidlDouble1D(alat3_ocn, alat3_ocn_bbl)
    
    call castSimple1DtoSidlDouble1D(aboxedge1_lon_ocn, aboxedge1_lon_ocn_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge1_lat_ocn, aboxedge1_lat_ocn_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge2_lon_ocn, aboxedge2_lon_ocn_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge2_lat_ocn, aboxedge2_lat_ocn_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge3_lon_ocn, aboxedge3_lon_ocn_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge3_lat_ocn, aboxedge3_lat_ocn_bbl)

    call castSimple1DtoSidlDouble1D(depth1_ocn, depth1_ocn_bbl)
    call castSimple1DtoSidlDouble1D(depth2_ocn, depth2_ocn_bbl)

    call castSimple2DtoSidlInt2D(ilandmask1_ocn, ilandmask1_ocn_bbl)
    call castSimple2DtoSidlInt2D(ilandmask2_ocn, ilandmask2_ocn_bbl)
    call castSimple2DtoSidlInt2D(ilandmask3_ocn, ilandmask3_ocn_bbl)

    call castSimple2DtoSidlDouble2D(tstar_ocn, tstar_ocn_bbl)
    call castSimple2DtoSidlDouble2D(sstar_ocn, sstar_ocn_bbl)
    call castSimple2DtoSidlDouble2D(ustar_ocn, ustar_ocn_bbl)
    call castSimple2DtoSidlDouble2D(vstar_ocn, vstar_ocn_bbl)
    call castSimple2DtoSidlDouble2D(albedo_ocn, albedo_ocn_bbl)
    
    
    call castSimple1DtoSidlInt1D(ias_out, ias_out_bbl)
    call castSimple1DtoSidlInt1D(iaf_out, iaf_out_bbl)
    call castSimple1DtoSidlInt1D(ips_out, ips_out_bbl)
    call castSimple1DtoSidlInt1D(ipf_out, ipf_out_bbl)

    jsf_out_bbl = jsf_out

    go_npstp_bbl = go_npstp
    go_iwstp_bbl = go_iwstp
    go_itstp_bbl = go_itstp
    go_ianav_bbl = go_ianav

    go_saln0_bbl = go_saln0
    go_rhoair_bbl = go_rhoair
    go_cd_bbl = go_cd
 
    call castSimple1DtoSidlDouble1D(go_ds, go_ds_bbl)
    go_dphi_bbl = go_dphi
    call castSimple1DtoSidlInt1D(go_ips, go_ips_bbl)
    call castSimple1DtoSidlInt1D(go_ipf, go_ipf_bbl)
    
    bg_usc_bbl = bg_usc
    go_rsc_bbl = go_rsc
    go_tsc_bbl = go_tsc
    go_dsc_bbl = go_dsc
    go_fsc_bbl = go_fsc
    go_gsc_bbl = go_gsc
    go_rh0sc_bbl = go_rh0sc
    go_rhosc_bbl = go_rhosc
    go_cpsc_bbl = go_cpsc

    call castSimple2DtoSidlInt2D(go_k1, go_k1_bbl)
    go_lmax_bbl = go_lmax

    call castSimple1DtoSidlDouble1D(go_dz, go_dz_bbl)
    call castSimple1DtoSidlDouble1D(go_dza, go_dza_bbl)
    
    call castSimple1DtoSidlInt1D(go_ias, go_ias_bbl)
    call castSimple1DtoSidlInt1D(go_iaf, go_iaf_bbl)
    
    go_jsf_bbl = go_jsf

    call cast0Simple1DtoSidlDouble1D(go_c, go_c_bbl)
    call cast0Simple1DtoSidlDouble1D(go_cv, go_cv_bbl)
    call cast0Simple1DtoSidlDouble1D(go_s, go_s_bbl)
    call cast0Simple1DtoSidlDouble1D(go_sv, go_sv_bbl)
    
    call castSimple4DtoSidlDouble4D(go_ts, go_ts_bbl)
    call castSimple4DtoSidlDouble4D(go_ts1, go_ts1_bbl)
   
    go_lin_bbl = go_lin
    go_lout_bbl = go_lout

! DO-NOT-DELETE splicer.end(ocean.goldstein.initialise)
end subroutine ocean_goldstein_initialise_mi


! 
! Method:  run[]
! 

recursive subroutine ocean_goldstein_run_mi(self, istep_ocn_bbl,               &
  latent_ocn_bbl, sensible_ocn_bbl, netsolar_ocn_bbl, netlong_ocn_bbl,         &
  conductflux_ocn_bbl, evap_ocn_bbl, precip_ocn_bbl, runoff_ocn_bbl,           &
  waterflux_ocn_bbl, ocean_stressx2_ocn_bbl, ocean_stressy2_ocn_bbl,           &
  ocean_stressx3_ocn_bbl, ocean_stressy3_ocn_bbl, tstar_ocn_bbl,               &
  sstar_ocn_bbl, ustar_ocn_bbl, vstar_ocn_bbl, albedo_ocn_bbl,                 &
  test_energy_ocean_bbl, test_water_ocean_bbl, koverall_bbl, go_ts_bbl,        &
  go_ts1_bbl, go_cost_bbl, go_u_bbl, go_ext_bbl, retval, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use ocean_goldstein
  use sidl_double_array
  use ocean_goldstein_impl
  ! DO-NOT-DELETE splicer.begin(ocean.goldstein.run.use)
  ! Insert-Code-Here {ocean.goldstein.run.use} (use statements)
  
  use genie_control
  use castings
  
  ! DO-NOT-DELETE splicer.end(ocean.goldstein.run.use)
  implicit none
  type(ocean_goldstein_t) :: self ! in
  integer (kind=sidl_int) :: istep_ocn_bbl ! in
  type(sidl_double_2d) :: latent_ocn_bbl ! in
  type(sidl_double_2d) :: sensible_ocn_bbl ! in
  type(sidl_double_2d) :: netsolar_ocn_bbl ! in
  type(sidl_double_2d) :: netlong_ocn_bbl ! in
  type(sidl_double_2d) :: conductflux_ocn_bbl ! in
  type(sidl_double_2d) :: evap_ocn_bbl ! in
  type(sidl_double_2d) :: precip_ocn_bbl ! in
  type(sidl_double_2d) :: runoff_ocn_bbl ! in
  type(sidl_double_2d) :: waterflux_ocn_bbl ! in
  type(sidl_double_2d) :: ocean_stressx2_ocn_bbl ! in
  type(sidl_double_2d) :: ocean_stressy2_ocn_bbl ! in
  type(sidl_double_2d) :: ocean_stressx3_ocn_bbl ! in
  type(sidl_double_2d) :: ocean_stressy3_ocn_bbl ! in
  type(sidl_double_2d) :: tstar_ocn_bbl ! out
  type(sidl_double_2d) :: sstar_ocn_bbl ! out
  type(sidl_double_2d) :: ustar_ocn_bbl ! out
  type(sidl_double_2d) :: vstar_ocn_bbl ! out
  type(sidl_double_2d) :: albedo_ocn_bbl ! out
  real (kind=sidl_double) :: test_energy_ocean_bbl ! out
  real (kind=sidl_double) :: test_water_ocean_bbl ! out
  integer (kind=sidl_int) :: koverall_bbl ! in
  type(sidl_double_4d) :: go_ts_bbl ! inout
  type(sidl_double_4d) :: go_ts1_bbl ! out
  type(sidl_double_2d) :: go_cost_bbl ! out
  type(sidl_double_4d) :: go_u_bbl ! out
  character (len=*) :: go_ext_bbl ! out
  integer (kind=sidl_int) :: retval ! out
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(ocean.goldstein.run)
! Insert-Code-Here {ocean.goldstein.run} (run method)
    
    integer istep_ocn

    real latent_ocn(ilon1_ocn, ilat1_ocn)
    real sensible_ocn(ilon1_ocn, ilat1_ocn)
    real netsolar_ocn(ilon1_ocn, ilat1_ocn)
    real netlong_ocn(ilon1_ocn, ilat1_ocn)
    real conductflux_ocn(ilon1_ocn, ilat1_ocn)
    real evap_ocn(ilon1_ocn, ilat1_ocn)
    real precip_ocn(ilon1_ocn, ilat1_ocn)
    real runoff_ocn(ilon1_ocn, ilat1_ocn)
    real waterflux_ocn(ilon1_ocn, ilat1_ocn)
    real ocean_stressx2_ocn(ilon1_ocn, ilat1_ocn)
    real ocean_stressy2_ocn(ilon1_ocn, ilat1_ocn)
    real ocean_stressx3_ocn(ilon1_ocn, ilat1_ocn)
    real ocean_stressy3_ocn(ilon1_ocn, ilat1_ocn)
    real tstar_ocn(ilon1_ocn, ilat1_ocn)
    real sstar_ocn(ilon1_ocn, ilat1_ocn)
    real ustar_ocn(ilon1_ocn, ilat1_ocn)
    real vstar_ocn(ilon1_ocn, ilat1_ocn)
    real albedo_ocn(ilon1_ocn, ilat1_ocn)
        
    real test_energy_ocean
    real test_water_ocean
    integer koverall
    real go_ts(1:intrac_ocn,1:ilon1_ocn,1:ilat1_ocn,1:inl1_ocn)
    real go_ts1(1:intrac_ocn,1:ilon1_ocn,1:ilat1_ocn,1:inl1_ocn)
    real go_cost(1:ilon1_ocn,1:ilat1_ocn)
    real go_u(1:3,1:ilon1_ocn,1:ilat1_ocn,1:inl1_ocn)
    character go_ext*3
       
    istep_ocn = istep_ocn_bbl

    call castSidlDouble2DToSimple2D(latent_ocn_bbl, latent_ocn)
    call castSidlDouble2DtoSimple2D(sensible_ocn_bbl, sensible_ocn)
    call castSidlDouble2DtoSimple2D(netsolar_ocn_bbl, netsolar_ocn)
    call castSidlDouble2DtoSimple2D(netlong_ocn_bbl, netlong_ocn)
    call castSidlDouble2DtoSimple2D(conductflux_ocn_bbl, conductflux_ocn)
    call castSidlDouble2DtoSimple2D(evap_ocn_bbl, evap_ocn)
    call castSidlDouble2DtoSimple2D(precip_ocn_bbl, precip_ocn)
    call castSidlDouble2DtoSimple2D(runoff_ocn_bbl, runoff_ocn)
    call castSidlDouble2DtoSimple2D(waterflux_ocn_bbl, waterflux_ocn)
    call castSidlDouble2DtoSimple2D(ocean_stressx2_ocn_bbl, ocean_stressx2_ocn)
    call castSidlDouble2DtoSimple2D(ocean_stressy2_ocn_bbl, ocean_stressy2_ocn)
    call castSidlDouble2DtoSimple2D(ocean_stressx3_ocn_bbl, ocean_stressx3_ocn)
    call castSidlDouble2DtoSimple2D(ocean_stressy3_ocn_bbl, ocean_stressy3_ocn)
        
    koverall = koverall_bbl

    call goldstein(istep_ocn, latent_ocn, sensible_ocn, &
                   netsolar_ocn, netlong_ocn, conductflux_ocn, evap_ocn,  &
                   precip_ocn, runoff_ocn, waterflux_ocn, ocean_stressx2_ocn, &
                   ocean_stressy2_ocn, ocean_stressx3_ocn, &
                   ocean_stressy3_ocn, tstar_ocn, sstar_ocn, ustar_ocn, vstar_ocn, &
                   albedo_ocn, test_energy_ocean, test_water_ocean, koverall, &
                   go_ts, go_ts1, go_cost, go_u, go_ext)
    
    call create2dcol(ilon1_ocn, ilat1_ocn, tstar_ocn_bbl)
    call create2dcol(ilon1_ocn, ilat1_ocn, sstar_ocn_bbl)
    call create2dcol(ilon1_ocn, ilat1_ocn, ustar_ocn_bbl)
    call create2dcol(ilon1_ocn, ilat1_ocn, vstar_ocn_bbl)
    call create2dcol(ilon1_ocn, ilat1_ocn, albedo_ocn_bbl) 
    
    call castSimple2DtoSidlDouble2D(tstar_ocn, tstar_ocn_bbl)
    call castSimple2DtoSidlDouble2D(sstar_ocn, sstar_ocn_bbl)
    call castSimple2DtoSidlDouble2D(ustar_ocn, ustar_ocn_bbl)
    call castSimple2DtoSidlDouble2D(vstar_ocn, vstar_ocn_bbl)
    call castSimple2DtoSidlDouble2D(albedo_ocn, albedo_ocn_bbl)

    go_ext_bbl = go_ext
  
! DO-NOT-DELETE splicer.end(ocean.goldstein.run)
end subroutine ocean_goldstein_run_mi


! DO-NOT-DELETE splicer.begin(_miscellaneous_code_end)
! Insert-Code-Here {_miscellaneous_code_end} (extra code)
! DO-NOT-DELETE splicer.end(_miscellaneous_code_end)
