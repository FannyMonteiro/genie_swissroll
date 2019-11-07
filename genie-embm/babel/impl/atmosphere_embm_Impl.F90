! 
! File:          atmosphere_embm_Impl.F90
! Symbol:        atmosphere.embm-v0.1.1
! Symbol Type:   class
! Babel Version: 0.11.0
! sidl Created:  20061123 18:01:34 GMT
! Generated:     20061123 18:01:50 GMT
! Description:   Server-side implementation for atmosphere.embm
! 
! WARNING: Automatically generated; only changes within splicers preserved
! 
! babel-version = 0.11.0
! source-line   = 94
! source-url    = file:/homes/sp1003/genie-babel/genie-embm/babel/embm.sidl
! 


! 
! Symbol "atmosphere.embm" (version 0.1.1)
! 


#include "atmosphere_embm_interface_fAbbrev.h"
#include "atmosphere_embm_fAbbrev.h"
#include "sidl_ClassInfo_fAbbrev.h"
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

recursive subroutine atmosphere_embm__ctor_mi(self, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use atmosphere_embm
  use atmosphere_embm_impl
  ! DO-NOT-DELETE splicer.begin(atmosphere.embm._ctor.use)
  ! Insert-Code-Here {atmosphere.embm._ctor.use} (use statements)
  ! DO-NOT-DELETE splicer.end(atmosphere.embm._ctor.use)
  implicit none
  type(atmosphere_embm_t) :: self ! in
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(atmosphere.embm._ctor)
! Insert-Code-Here {atmosphere.embm._ctor} (_ctor method)
! DO-NOT-DELETE splicer.end(atmosphere.embm._ctor)
end subroutine atmosphere_embm__ctor_mi


! 
! Method:  _dtor[]
! Class destructor called when the class is deleted.
! 

recursive subroutine atmosphere_embm__dtor_mi(self, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use atmosphere_embm
  use atmosphere_embm_impl
  ! DO-NOT-DELETE splicer.begin(atmosphere.embm._dtor.use)
  ! Insert-Code-Here {atmosphere.embm._dtor.use} (use statements)
  ! DO-NOT-DELETE splicer.end(atmosphere.embm._dtor.use)
  implicit none
  type(atmosphere_embm_t) :: self ! in
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(atmosphere.embm._dtor)
! Insert-Code-Here {atmosphere.embm._dtor} (_dtor method)
! DO-NOT-DELETE splicer.end(atmosphere.embm._dtor)
end subroutine atmosphere_embm__dtor_mi


! 
! Method:  _load[]
! Static class initializer called exactly once before any user-defined method is dispatched
! 

recursive subroutine atmosphere_embm__load_mi(exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use atmosphere_embm
  use atmosphere_embm_impl
  ! DO-NOT-DELETE splicer.begin(atmosphere.embm._load.use)
  ! Insert-Code-Here {atmosphere.embm._load.use} (use statements)
  ! DO-NOT-DELETE splicer.end(atmosphere.embm._load.use)
  implicit none
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(atmosphere.embm._load)
! Insert-Code-Here {atmosphere.embm._load} (_load method)
! DO-NOT-DELETE splicer.end(atmosphere.embm._load)
end subroutine atmosphere_embm__load_mi


! 
! Method:  initialise[]
! 

recursive subroutine atmosphere_embm_initialise_mi(self, alon1_atm_bbl,        &
  alat1_atm_bbl, alon2_atm_bbl, alat2_atm_bbl, alon3_atm_bbl, alat3_atm_bbl,   &
  aboxedge1_lon_atm_bbl, aboxedge1_lat_atm_bbl, aboxedge2_lon_atm_bbl,         &
  aboxedge2_lat_atm_bbl, aboxedge3_lon_atm_bbl, aboxedge3_lat_atm_bbl,         &
  ilandmask1_atm_bbl, ilandmask2_atm_bbl, ilandmask3_atm_bbl, ias_out_bbl,     &
  iaf_out_bbl, ips_out_bbl, ipf_out_bbl, tstar_ocn_bbl, koverall_total_bbl,    &
  co2_atm_bbl, ocean_stressx2_ocn_bbl, ocean_stressy2_ocn_bbl,                 &
  ocean_stressx3_ocn_bbl, ocean_stressy3_ocn_bbl, tstar_atm_bbl,               &
  surf_qstar_atm_bbl, atmos_dt_tim_bbl, go_solconst_bbl, retval, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use atmosphere_embm
  use sidl_double_array
  use sidl_int_array
  use atmosphere_embm_impl
  ! DO-NOT-DELETE splicer.begin(atmosphere.embm.initialise.use)
  ! Insert-Code-Here {atmosphere.embm.initialise.use} (use statements)
  
  use genie_control
  use castings
  
  ! DO-NOT-DELETE splicer.end(atmosphere.embm.initialise.use)
  implicit none
  type(atmosphere_embm_t) :: self ! in
  type(sidl_double_1d) :: alon1_atm_bbl ! out
  type(sidl_double_1d) :: alat1_atm_bbl ! out
  type(sidl_double_1d) :: alon2_atm_bbl ! out
  type(sidl_double_1d) :: alat2_atm_bbl ! out
  type(sidl_double_1d) :: alon3_atm_bbl ! out
  type(sidl_double_1d) :: alat3_atm_bbl ! out
  type(sidl_double_1d) :: aboxedge1_lon_atm_bbl ! out
  type(sidl_double_1d) :: aboxedge1_lat_atm_bbl ! out
  type(sidl_double_1d) :: aboxedge2_lon_atm_bbl ! out
  type(sidl_double_1d) :: aboxedge2_lat_atm_bbl ! out
  type(sidl_double_1d) :: aboxedge3_lon_atm_bbl ! out
  type(sidl_double_1d) :: aboxedge3_lat_atm_bbl ! out
  type(sidl_int_2d) :: ilandmask1_atm_bbl ! out
  type(sidl_int_2d) :: ilandmask2_atm_bbl ! out
  type(sidl_int_2d) :: ilandmask3_atm_bbl ! out
  type(sidl_int_1d) :: ias_out_bbl ! in
  type(sidl_int_1d) :: iaf_out_bbl ! in
  type(sidl_int_1d) :: ips_out_bbl ! in
  type(sidl_int_1d) :: ipf_out_bbl ! in
  type(sidl_double_2d) :: tstar_ocn_bbl ! in
  integer (kind=sidl_int) :: koverall_total_bbl ! in
  type(sidl_double_2d) :: co2_atm_bbl ! out
  type(sidl_double_2d) :: ocean_stressx2_ocn_bbl ! out
  type(sidl_double_2d) :: ocean_stressy2_ocn_bbl ! out
  type(sidl_double_2d) :: ocean_stressx3_ocn_bbl ! out
  type(sidl_double_2d) :: ocean_stressy3_ocn_bbl ! out
  type(sidl_double_2d) :: tstar_atm_bbl ! out
  type(sidl_double_2d) :: surf_qstar_atm_bbl ! out
  real (kind=sidl_double) :: atmos_dt_tim_bbl ! out
  real (kind=sidl_double) :: go_solconst_bbl ! out
  integer (kind=sidl_int) :: retval ! out
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(atmosphere.embm.initialise)
! Insert-Code-Here {atmosphere.embm.initialise} (initialise method)

    real,dimension(:), ALLOCATABLE :: alon1_atm
    real,dimension(:), ALLOCATABLE :: alat1_atm
    real,dimension(:), ALLOCATABLE :: alon2_atm
    real,dimension(:), ALLOCATABLE :: alat2_atm
    real,dimension(:), ALLOCATABLE :: alon3_atm
    real,dimension(:), ALLOCATABLE :: alat3_atm
    real,dimension(:), ALLOCATABLE :: aboxedge1_lon_atm
    real,dimension(:), ALLOCATABLE :: aboxedge1_lat_atm
    real,dimension(:), ALLOCATABLE :: aboxedge2_lon_atm
    real,dimension(:), ALLOCATABLE :: aboxedge2_lat_atm
    real,dimension(:), ALLOCATABLE :: aboxedge3_lon_atm
    real,dimension(:), ALLOCATABLE :: aboxedge3_lat_atm
    integer, dimension(:,:), ALLOCATABLE :: ilandmask1_atm
    integer, dimension(:,:), ALLOCATABLE :: ilandmask2_atm
    integer, dimension(:,:), ALLOCATABLE :: ilandmask3_atm
    integer, dimension(:), ALLOCATABLE :: ias_out
    integer, dimension(:), ALLOCATABLE :: iaf_out
    integer, dimension(:), ALLOCATABLE :: ips_out
    integer, dimension(:), ALLOCATABLE :: ipf_out
    real,dimension(:,:), ALLOCATABLE :: tstar_ocn
    integer koverall_total_int
    real, dimension(:,:), ALLOCATABLE :: co2_atm

    real,dimension(:,:), ALLOCATABLE :: ocean_stressx2_ocn
    real,dimension(:,:), ALLOCATABLE :: ocean_stressy2_ocn
    real,dimension(:,:), ALLOCATABLE :: ocean_stressx3_ocn
    real,dimension(:,:), ALLOCATABLE :: ocean_stressy3_ocn
    real,dimension(:,:), ALLOCATABLE :: tstar_atm
    real,dimension(:,:), ALLOCATABLE :: surf_qstar_atm
    real atmos_dt_tim
    real go_solconst

    allocate(alon1_atm(ilon1_atm))
    allocate(alat1_atm(ilat1_atm))
    allocate(alon2_atm(ilon2_atm))
    allocate(alat2_atm(ilat2_atm))
    allocate(alon3_atm(ilon3_atm))
    allocate(alat3_atm(ilat3_atm))
    allocate(aboxedge1_lon_atm(ilon1_atm+1))
    allocate(aboxedge1_lat_atm(ilat1_atm+1))
    allocate(aboxedge2_lon_atm(ilon2_atm+1))
    allocate(aboxedge2_lat_atm(ilat2_atm+1))
    allocate(aboxedge3_lon_atm(ilon3_atm+1))
    allocate(aboxedge3_lat_atm(ilat3_atm+1))
    allocate(ilandmask1_atm(ilon1_atm,ilat1_atm))
    allocate(ilandmask2_atm(ilon2_atm,ilat2_atm))
    allocate(ilandmask3_atm(ilon3_atm,ilat3_atm))
    allocate(ias_out(ilat1_ocn))
    allocate(iaf_out(ilat1_ocn))
    allocate(ips_out(ilat1_ocn))
    allocate(ipf_out(ilat1_ocn))
    allocate(tstar_ocn(ilon1_ocn,ilat1_ocn))
    
    koverall_total_int = koverall_total_bbl
    
    allocate(co2_atm(ilon1_atm,ilat1_atm))
    
    allocate(ocean_stressx2_ocn(ilon1_ocn,ilat1_ocn))
    allocate(ocean_stressy2_ocn(ilon1_ocn,ilat1_ocn))
    allocate(ocean_stressx3_ocn(ilon1_ocn,ilat1_ocn))
    allocate(ocean_stressy3_ocn(ilon1_ocn,ilat1_ocn))
    allocate(tstar_atm(ilon1_atm,ilat1_atm))
    allocate(surf_qstar_atm(ilon1_atm,ilat1_atm))

    call castSidlInt1DToSimple1D(ias_out_bbl, ias_out)
    call castSidlInt1DToSimple1D(iaf_out_bbl, iaf_out)
    call castSidlInt1DToSimple1D(ips_out_bbl, ips_out)
    call castSidlInt1DToSimple1D(ipf_out_bbl, ipf_out)
    
    call castSidlDouble2DToSimple2D(tstar_ocn_bbl, tstar_ocn)

    call initialise_embm(alon1_atm, alat1_atm, alon2_atm, alat2_atm, alon3_atm, &
             alat3_atm, aboxedge1_lon_atm, aboxedge1_lat_atm, aboxedge2_lon_atm, &
             aboxedge2_lat_atm, aboxedge3_lon_atm, aboxedge3_lat_atm, ilandmask1_atm, &
             ilandmask2_atm, ilandmask3_atm, ias_out, iaf_out, ips_out, ipf_out, &
             tstar_ocn, koverall_total_int, co2_atm, ocean_stressx2_ocn, &
             ocean_stressy2_ocn, ocean_stressx3_ocn, ocean_stressy3_ocn, &
             tstar_atm, surf_qstar_atm, atmos_dt_tim, go_solconst)
    
    
    call sidl_double__array_create1d_m(ilon1_atm, alon1_atm_bbl)
    call sidl_double__array_create1d_m(ilat1_atm, alat1_atm_bbl)
    call sidl_double__array_create1d_m(ilon2_atm, alon2_atm_bbl)
    call sidl_double__array_create1d_m(ilat2_atm, alat2_atm_bbl)
    call sidl_double__array_create1d_m(ilon3_atm, alon3_atm_bbl)
    call sidl_double__array_create1d_m(ilat3_atm, alat3_atm_bbl)
    call sidl_double__array_create1d_m(ilon1_atm+1, aboxedge1_lon_atm_bbl)
    call sidl_double__array_create1d_m(ilat1_atm+1, aboxedge1_lat_atm_bbl)
    call sidl_double__array_create1d_m(ilon2_atm+1, aboxedge2_lon_atm_bbl)
    call sidl_double__array_create1d_m(ilat2_atm+1, aboxedge2_lat_atm_bbl)
    call sidl_double__array_create1d_m(ilon3_atm+1, aboxedge3_lon_atm_bbl)
    call sidl_double__array_create1d_m(ilat3_atm+1, aboxedge3_lat_atm_bbl)

    call castSimple1DtoSidlDouble1D(alon1_atm, alon1_atm_bbl)
    call castSimple1DtoSidlDouble1D(alat1_atm, alat1_atm_bbl)
    call castSimple1DtoSidlDouble1D(alon2_atm, alon2_atm_bbl)
    call castSimple1DtoSidlDouble1D(alat2_atm, alat2_atm_bbl)
    call castSimple1DtoSidlDouble1D(alon3_atm, alon3_atm_bbl)
    call castSimple1DtoSidlDouble1D(alat3_atm, alat3_atm_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge1_lon_atm, aboxedge1_lon_atm_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge1_lat_atm, aboxedge1_lat_atm_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge2_lon_atm, aboxedge2_lon_atm_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge2_lat_atm, aboxedge2_lat_atm_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge3_lon_atm, aboxedge3_lon_atm_bbl)
    call castSimple1DtoSidlDouble1D(aboxedge3_lat_atm, aboxedge3_lat_atm_bbl)
    
    call create2dCol(ilon1_atm, ilat1_atm, ilandmask1_atm_bbl)
    call create2dCol(ilon2_atm, ilat2_atm, ilandmask2_atm_bbl)
    call create2dCol(ilon3_atm, ilat3_atm, ilandmask3_atm_bbl)

    
    call castSimple2DtoSidlInt2D(ilandmask1_atm, ilandmask1_atm_bbl)
    call castSimple2DtoSidlInt2D(ilandmask2_atm, ilandmask2_atm_bbl)
    call castSimple2DtoSidlInt2D(ilandmask3_atm, ilandmask3_atm_bbl)
    
    call create2dCol(ilon1_atm, ilat1_atm, co2_atm_bbl)
    call create2dCol(ilon1_ocn, ilat1_ocn, ocean_stressx2_ocn_bbl)
    call create2dCol(ilon1_ocn, ilat1_ocn, ocean_stressy2_ocn_bbl)
    call create2dCol(ilon1_ocn, ilat1_ocn, ocean_stressx3_ocn_bbl)
    call create2dCol(ilon1_ocn, ilat1_ocn, ocean_stressy3_ocn_bbl)
    call create2dCol(ilon1_atm, ilat1_atm, tstar_atm_bbl)
    call create2dCol(ilon1_atm, ilat1_atm, surf_qstar_atm_bbl)
    
    call castSimple2DtoSidlDouble2D(co2_atm, co2_atm_bbl)
    call castSimple2DtoSidlDouble2D(ocean_stressx2_ocn, ocean_stressx2_ocn_bbl)
    call castSimple2DtoSidlDouble2D(ocean_stressy2_ocn, ocean_stressy2_ocn_bbl)
    call castSimple2DtoSidlDouble2D(ocean_stressx3_ocn, ocean_stressx3_ocn_bbl)
    call castSimple2DtoSidlDouble2D(ocean_stressy3_ocn, ocean_stressy3_ocn_bbl)
    call castSimple2DtoSidlDouble2D(tstar_atm, tstar_atm_bbl)
    call castSimple2DtoSidlDouble2D(surf_qstar_atm, surf_qstar_atm_bbl)
    
    go_solconst_bbl = go_solconst

! DO-NOT-DELETE splicer.end(atmosphere.embm.initialise)
end subroutine atmosphere_embm_initialise_mi


! 
! Method:  run[]
! 

recursive subroutine atmosphere_embm_run_mi(self, istep_atm_bbl,               &
  surf_latent_atm_bbl, surf_sensible_atm_bbl, netsolar_atm_bbl,                &
  netlong_atm_bbl, evap_atm_bbl, precip_atm_bbl, ocean_stressx2_ocn_bbl,       &
  ocean_stressy2_ocn_bbl, ocean_stressx3_ocn_bbl, ocean_stressy3_ocn_bbl,      &
  tstar_atm_bbl, surf_qstar_atm_bbl, koverall_bbl, retval, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use atmosphere_embm
  use sidl_double_array
  use atmosphere_embm_impl
  ! DO-NOT-DELETE splicer.begin(atmosphere.embm.run.use)
  
  use genie_control
  use castings
  
  ! Insert-Code-Here {atmosphere.embm.run.use} (use statements)
  ! DO-NOT-DELETE splicer.end(atmosphere.embm.run.use)
  implicit none
  type(atmosphere_embm_t) :: self ! in
  integer (kind=sidl_int) :: istep_atm_bbl ! in
  type(sidl_double_2d) :: surf_latent_atm_bbl ! in
  type(sidl_double_2d) :: surf_sensible_atm_bbl ! in
  type(sidl_double_2d) :: netsolar_atm_bbl ! in
  type(sidl_double_2d) :: netlong_atm_bbl ! in
  type(sidl_double_2d) :: evap_atm_bbl ! in
  type(sidl_double_2d) :: precip_atm_bbl ! in
  type(sidl_double_2d) :: ocean_stressx2_ocn_bbl ! out
  type(sidl_double_2d) :: ocean_stressy2_ocn_bbl ! out
  type(sidl_double_2d) :: ocean_stressx3_ocn_bbl ! out
  type(sidl_double_2d) :: ocean_stressy3_ocn_bbl ! out
  type(sidl_double_2d) :: tstar_atm_bbl ! out
  type(sidl_double_2d) :: surf_qstar_atm_bbl ! out
  integer (kind=sidl_int) :: koverall_bbl ! in
  integer (kind=sidl_int) :: retval ! out
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(atmosphere.embm.run)
! Insert-Code-Here {atmosphere.embm.run} (run method)

  integer istep_atm
  real surf_latent_atm(ilon1_atm,ilat1_atm)
  real surf_sensible_atm(ilon1_atm,ilat1_atm)
  real netsolar_atm(ilon1_atm,ilat1_atm)
  real netlong_atm(ilon1_atm,ilat1_atm)
  real evap_atm(ilon1_atm,ilat1_atm)
  real precip_atm(ilon1_atm,ilat1_atm)
  real ocean_stressx2_ocn(ilon1_ocn, ilat1_ocn)
  real ocean_stressy2_ocn(ilon1_ocn, ilat1_ocn)
  real ocean_stressx3_ocn(ilon1_ocn, ilat1_ocn)
  real ocean_stressy3_ocn(ilon1_ocn, ilat1_ocn)
  real tstar_atm(ilon1_atm,ilat1_atm)
  real surf_qstar_atm(ilon1_atm,ilat1_atm)
  integer koverall
  
  istep_atm = istep_atm_bbl

  call castSidlDouble2DToSimple2D(surf_latent_atm_bbl, surf_latent_atm)
  call castSidlDouble2DToSimple2D(surf_sensible_atm_bbl, surf_sensible_atm)
  call castSidlDouble2DToSimple2D(netsolar_atm_bbl, netsolar_atm)
  call castSidlDouble2DToSimple2D(netlong_atm_bbl, netlong_atm)
  call castSidlDouble2DToSimple2D(evap_atm_bbl, evap_atm)
  call castSidlDouble2DToSimple2D(precip_atm_bbl, precip_atm)
    
  koverall = koverall_bbl
  
  call embm(istep_atm, surf_latent_atm, surf_sensible_atm, &
           netsolar_atm, netlong_atm, evap_atm, precip_atm, &
           ocean_stressx2_ocn, ocean_stressy2_ocn, ocean_stressx3_ocn, &
           ocean_stressy3_ocn, tstar_atm, surf_qstar_atm, koverall)

  call create2dcol(ilon1_ocn, ilat1_ocn, ocean_stressx2_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, ocean_stressy2_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, ocean_stressx3_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, ocean_stressy3_ocn_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, tstar_atm_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, surf_qstar_atm_bbl)

  call castSimple2DtoSidlDouble2D(ocean_stressx2_ocn, ocean_stressx2_ocn_bbl)
  call castSimple2DtoSidlDouble2D(ocean_stressy2_ocn, ocean_stressy2_ocn_bbl)
  call castSimple2DtoSidlDouble2D(ocean_stressx3_ocn, ocean_stressx3_ocn_bbl)
  call castSimple2DtoSidlDouble2D(ocean_stressy3_ocn, ocean_stressy3_ocn_bbl)
  call castSimple2DtoSidlDouble2D(tstar_atm, tstar_atm_bbl)
  call castSimple2DtoSidlDouble2D(surf_qstar_atm, surf_qstar_atm_bbl)

! DO-NOT-DELETE splicer.end(atmosphere.embm.run)
end subroutine atmosphere_embm_run_mi


! 
! Method:  run_surflux[]
! 

recursive subroutine atmosphere_embm_run_surflux_mi(self, istep_ocn_bbl,       &
  tstar_ocn_bbl, sstar_ocn_bbl, tstar_atm_bbl, surf_qstar_atm_bbl,             &
  hght_sic_bbl, frac_sic_bbl, temp_sic_bbl, albd_sic_bbl,                      &
  ocean_stressx2_ocn_bbl, ocean_stressy2_ocn_bbl, ocean_stressx3_ocn_bbl,      &
  ocean_stressy3_ocn_bbl, co2_atm_bbl, albedo_ocn_bbl, latent_ocn_bbl,         &
  sensible_ocn_bbl, netsolar_ocn_bbl, netlong_ocn_bbl, evap_ocn_bbl,           &
  precip_ocn_bbl, runoff_ocn_bbl, surf_latent_atm_bbl, surf_sensible_atm_bbl,  &
  netsolar_atm_bbl, netlong_atm_bbl, evap_atm_bbl, precip_atm_bbl,             &
  dhght_sic_bbl, dfrac_sic_bbl, atmos_lowestlh_atm_bbl, go_solfor_bbl,         &
  go_tau_bbl, retval, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use atmosphere_embm
  use sidl_double_array
  use atmosphere_embm_impl
  ! DO-NOT-DELETE splicer.begin(atmosphere.embm.run_surflux.use)
  ! Insert-Code-Here {atmosphere.embm.run_surflux.use} (use statements)
  
  use genie_control
  use castings
  
  ! DO-NOT-DELETE splicer.end(atmosphere.embm.run_surflux.use)
  implicit none
  type(atmosphere_embm_t) :: self ! in
  integer (kind=sidl_int) :: istep_ocn_bbl ! in
  type(sidl_double_2d) :: tstar_ocn_bbl ! in
  type(sidl_double_2d) :: sstar_ocn_bbl ! in
  type(sidl_double_2d) :: tstar_atm_bbl ! in
  type(sidl_double_2d) :: surf_qstar_atm_bbl ! in
  type(sidl_double_2d) :: hght_sic_bbl ! in
  type(sidl_double_2d) :: frac_sic_bbl ! in
  type(sidl_double_2d) :: temp_sic_bbl ! in
  type(sidl_double_2d) :: albd_sic_bbl ! in
  type(sidl_double_2d) :: ocean_stressx2_ocn_bbl ! in
  type(sidl_double_2d) :: ocean_stressy2_ocn_bbl ! in
  type(sidl_double_2d) :: ocean_stressx3_ocn_bbl ! in
  type(sidl_double_2d) :: ocean_stressy3_ocn_bbl ! in
  type(sidl_double_2d) :: co2_atm_bbl ! in
  type(sidl_double_2d) :: albedo_ocn_bbl ! out
  type(sidl_double_2d) :: latent_ocn_bbl ! out
  type(sidl_double_2d) :: sensible_ocn_bbl ! out
  type(sidl_double_2d) :: netsolar_ocn_bbl ! out
  type(sidl_double_2d) :: netlong_ocn_bbl ! out
  type(sidl_double_2d) :: evap_ocn_bbl ! out
  type(sidl_double_2d) :: precip_ocn_bbl ! out
  type(sidl_double_2d) :: runoff_ocn_bbl ! out
  type(sidl_double_2d) :: surf_latent_atm_bbl ! out
  type(sidl_double_2d) :: surf_sensible_atm_bbl ! out
  type(sidl_double_2d) :: netsolar_atm_bbl ! out
  type(sidl_double_2d) :: netlong_atm_bbl ! out
  type(sidl_double_2d) :: evap_atm_bbl ! out
  type(sidl_double_2d) :: precip_atm_bbl ! out
  type(sidl_double_2d) :: dhght_sic_bbl ! out
  type(sidl_double_2d) :: dfrac_sic_bbl ! out
  type(sidl_double_2d) :: atmos_lowestlh_atm_bbl ! out
  type(sidl_double_2d) :: go_solfor_bbl ! out
  type(sidl_double_3d) :: go_tau_bbl ! out
  integer (kind=sidl_int) :: retval ! out
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(atmosphere.embm.run_surflux)
! Insert-Code-Here {atmosphere.embm.run_surflux} (run_surflux method)

  integer istep_ocn
  real tstar_ocn(ilon1_ocn,ilat1_ocn)
  real sstar_ocn(ilon1_ocn,ilat1_ocn)
  real tstar_atm(ilon1_atm,ilat1_atm)
  real surf_qstar_atm(ilon1_atm,ilat1_atm)
  real hght_sic(ilon1_sic,ilat1_sic)
  real frac_sic(ilon1_sic,ilat1_sic)
  real temp_sic(ilon1_sic,ilat1_sic)
  real albd_sic(ilon1_sic,ilat1_sic)
  real ocean_stressx2_ocn(ilon1_ocn, ilat1_ocn)
  real ocean_stressy2_ocn(ilon1_ocn, ilat1_ocn)
  real ocean_stressx3_ocn(ilon1_ocn, ilat1_ocn)
  real ocean_stressy3_ocn(ilon1_ocn, ilat1_ocn)
  real co2_atm(ilon1_atm,ilat1_atm)
  real albedo_ocn(ilon1_ocn, ilat1_ocn)
  real latent_ocn(ilon1_ocn,ilat1_ocn)
  real sensible_ocn(ilon1_ocn, ilat1_ocn)
  real netsolar_ocn(ilon1_ocn, ilat1_ocn)
  real netlong_ocn(ilon1_ocn, ilat1_ocn)
  real evap_ocn(ilon1_ocn, ilat1_ocn)
  real precip_ocn(ilon1_ocn, ilat1_ocn)
  real runoff_ocn(ilon1_ocn, ilat1_ocn)
  real surf_latent_atm(ilon1_atm,ilat1_atm)
  real surf_sensible_atm(ilon1_atm,ilat1_atm)
  real netsolar_atm(ilon1_atm,ilat1_atm)
  real netlong_atm(ilon1_atm,ilat1_atm)
  real evap_atm(ilon1_atm,ilat1_atm)
  real precip_atm(ilon1_atm,ilat1_atm)
  real dhght_sic(ilon1_sic,ilat1_sic)
  real dfrac_sic(ilon1_sic,ilat1_sic)
  real atmos_lowestlh_atm(ilon1_atm,ilat1_atm)
  real go_solfor(ilat1_ocn,gen_maxnyr)
  real go_tau(2,ilon1_ocn,ilat1_ocn)

  istep_ocn = istep_ocn_bbl
    
  call castSidlDouble2DToSimple2D(tstar_ocn_bbl, tstar_ocn)    
  call castSidlDouble2DToSimple2D(sstar_ocn_bbl, sstar_ocn)
  call castSidlDouble2DToSimple2D(tstar_atm_bbl, tstar_atm)
  call castSidlDouble2DToSimple2D(surf_qstar_atm_bbl, surf_qstar_atm)
  call castSidlDouble2DToSimple2D(hght_sic_bbl, hght_sic)
  call castSidlDouble2DToSimple2D(frac_sic_bbl, frac_sic)
  call castSidlDouble2DToSimple2D(temp_sic_bbl, temp_sic)
  call castSidlDouble2DToSimple2D(albd_sic_bbl, albd_sic)
  call castSidlDouble2DToSimple2D(ocean_stressx2_ocn_bbl, ocean_stressx2_ocn)
  call castSidlDouble2DToSimple2D(ocean_stressy2_ocn_bbl, ocean_stressy2_ocn)
  call castSidlDouble2DToSimple2D(ocean_stressx3_ocn_bbl, ocean_stressx3_ocn)
  call castSidlDouble2DToSimple2D(ocean_stressy3_ocn_bbl, ocean_stressy3_ocn)
  call castSidlDouble2DToSimple2D(co2_atm_bbl, co2_atm)

  call surflux(istep_ocn, tstar_ocn, sstar_ocn, tstar_atm, surf_qstar_atm, hght_sic, &
         frac_sic, temp_sic, albd_sic, ocean_stressx2_ocn, &
         ocean_stressy2_ocn, ocean_stressx3_ocn, ocean_stressy3_ocn, &
         co2_atm, albedo_ocn, latent_ocn, sensible_ocn, netsolar_ocn, &
         netlong_ocn, evap_ocn, precip_ocn, runoff_ocn, surf_latent_atm, &
         surf_sensible_atm, netsolar_atm, netlong_atm, evap_atm, &
         precip_atm, dhght_sic, dfrac_sic, atmos_lowestlh_atm, &
         go_solfor, go_tau)

  call create2dcol(ilon1_ocn, ilat1_ocn, albedo_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, latent_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, sensible_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, netsolar_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, netlong_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, evap_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, precip_ocn_bbl)
  call create2dcol(ilon1_ocn, ilat1_ocn, runoff_ocn_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, surf_latent_atm_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, surf_sensible_atm_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, netsolar_atm_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, netlong_atm_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, evap_atm_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, precip_atm_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, dhght_sic_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, dfrac_sic_bbl)
  call create2dcol(ilon1_atm, ilat1_atm, atmos_lowestlh_atm_bbl)
  call create2dcol(ilat1_ocn, gen_maxnyr, go_solfor_bbl)
  call createCol([0, 0, 0], [1, ilon1_ocn-1, ilat1_ocn-1], go_tau_bbl)

  call castSimple2DtoSidlDouble2D(albedo_ocn, albedo_ocn_bbl)
  call castSimple2DtoSidlDouble2D(latent_ocn, latent_ocn_bbl)
  call castSimple2DtoSidlDouble2D(sensible_ocn, sensible_ocn_bbl)
  call castSimple2DtoSidlDouble2D(netsolar_ocn, netsolar_ocn_bbl)
  call castSimple2DtoSidlDouble2D(netlong_ocn, netlong_ocn_bbl)
  call castSimple2DtoSidlDouble2D(evap_ocn, evap_ocn_bbl)
  call castSimple2DtoSidlDouble2D(precip_ocn, precip_ocn_bbl)
  call castSimple2DtoSidlDouble2D(runoff_ocn, runoff_ocn_bbl)
  call castSimple2DtoSidlDouble2D(surf_latent_atm, surf_latent_atm_bbl)
  call castSimple2DtoSidlDouble2D(surf_sensible_atm, surf_sensible_atm_bbl)
  call castSimple2DtoSidlDouble2D(netsolar_atm, netsolar_atm_bbl)
  call castSimple2DtoSidlDouble2D(netlong_atm, netlong_atm_bbl)
  call castSimple2DtoSidlDouble2D(evap_atm, evap_atm_bbl)
  call castSimple2DtoSidlDouble2D(precip_atm, precip_atm_bbl)
  call castSimple2DtoSidlDouble2D(dhght_sic, dhght_sic_bbl)
  call castSimple2DtoSidlDouble2D(dfrac_sic, dfrac_sic_bbl)
  call castSimple2DtoSidlDouble2D(atmos_lowestlh_atm, atmos_lowestlh_atm_bbl)
  call castSimple2DtoSidlDouble2D(go_solfor, go_solfor_bbl)
  call castSimple3DtoSidlDouble3D(go_tau, go_tau_bbl)

! DO-NOT-DELETE splicer.end(atmosphere.embm.run_surflux)
end subroutine atmosphere_embm_run_surflux_mi


! DO-NOT-DELETE splicer.begin(_miscellaneous_code_end)
! Insert-Code-Here {_miscellaneous_code_end} (extra code)
! DO-NOT-DELETE splicer.end(_miscellaneous_code_end)
