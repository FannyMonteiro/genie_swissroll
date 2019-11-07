! 
! File:          seaice_goldseaice_Impl.F90
! Symbol:        seaice.goldseaice-v0.1.1
! Symbol Type:   class
! Babel Version: 0.11.0
! sidl Created:  20061123 17:49:18 GMT
! Generated:     20061123 17:49:31 GMT
! Description:   Server-side implementation for seaice.goldseaice
! 
! WARNING: Automatically generated; only changes within splicers preserved
! 
! babel-version = 0.11.0
! source-line   = 51
! source-url    = file:/homes/sp1003/genie-babel/genie-seaice/babel/goldseaice.sidl
! 


! 
! Symbol "seaice.goldseaice" (version 0.1.1)
! 


#include "seaice_goldseaice_fAbbrev.h"
#include "sidl_ClassInfo_fAbbrev.h"
#include "seaice_goldseaice_interface_fAbbrev.h"
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

recursive subroutine seaice_goldseaice__ctor_mi(self, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use seaice_goldseaice
  use seaice_goldseaice_impl
  ! DO-NOT-DELETE splicer.begin(seaice.goldseaice._ctor.use)
  ! Insert-Code-Here {seaice.goldseaice._ctor.use} (use statements)
  ! DO-NOT-DELETE splicer.end(seaice.goldseaice._ctor.use)
  implicit none
  type(seaice_goldseaice_t) :: self ! in
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(seaice.goldseaice._ctor)
! Insert-Code-Here {seaice.goldseaice._ctor} (_ctor method)
! DO-NOT-DELETE splicer.end(seaice.goldseaice._ctor)
end subroutine seaice_goldseaice__ctor_mi


! 
! Method:  _dtor[]
! Class destructor called when the class is deleted.
! 

recursive subroutine seaice_goldseaice__dtor_mi(self, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use seaice_goldseaice
  use seaice_goldseaice_impl
  ! DO-NOT-DELETE splicer.begin(seaice.goldseaice._dtor.use)
  ! Insert-Code-Here {seaice.goldseaice._dtor.use} (use statements)
  ! DO-NOT-DELETE splicer.end(seaice.goldseaice._dtor.use)
  implicit none
  type(seaice_goldseaice_t) :: self ! in
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(seaice.goldseaice._dtor)
! Insert-Code-Here {seaice.goldseaice._dtor} (_dtor method)
! DO-NOT-DELETE splicer.end(seaice.goldseaice._dtor)
end subroutine seaice_goldseaice__dtor_mi


! 
! Method:  _load[]
! Static class initializer called exactly once before any user-defined method is dispatched
! 

recursive subroutine seaice_goldseaice__load_mi(exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use seaice_goldseaice
  use seaice_goldseaice_impl
  ! DO-NOT-DELETE splicer.begin(seaice.goldseaice._load.use)
  ! Insert-Code-Here {seaice.goldseaice._load.use} (use statements)
  ! DO-NOT-DELETE splicer.end(seaice.goldseaice._load.use)
  implicit none
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(seaice.goldseaice._load)
! Insert-Code-Here {seaice.goldseaice._load} (_load method)
! DO-NOT-DELETE splicer.end(seaice.goldseaice._load)
end subroutine seaice_goldseaice__load_mi


! 
! Method:  initialise[]
! 

recursive subroutine seaice_goldseaice_initialise_mi(self, alon1_sic_bbl,      &
  alat1_sic_bbl, alon2_sic_bbl, alat2_sic_bbl, alon3_sic_bbl, alat3_sic_bbl,   &
  aboxedge1_lon_sic_bbl, aboxedge1_lat_sic_bbl, aboxedge2_lon_sic_bbl,         &
  aboxedge2_lat_sic_bbl, aboxedge3_lon_sic_bbl, aboxedge3_lat_sic_bbl,         &
  ilandmask1_sic_bbl, ilandmask2_sic_bbl, ilandmask3_sic_bbl,                  &
  koverall_total_bbl, hght_sic_bbl, frac_sic_bbl, temp_sic_bbl, albd_sic_bbl,  &
  test_energy_seaice_bbl, retval, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use seaice_goldseaice
  use sidl_double_array
  use sidl_int_array
  use seaice_goldseaice_impl
  ! DO-NOT-DELETE splicer.begin(seaice.goldseaice.initialise.use)
  ! Insert-Code-Here {seaice.goldseaice.initialise.use} (use statements)
  use genie_control
  use castings
  ! DO-NOT-DELETE splicer.end(seaice.goldseaice.initialise.use)
  implicit none
  type(seaice_goldseaice_t) :: self ! in
  type(sidl_double_1d) :: alon1_sic_bbl ! out
  type(sidl_double_1d) :: alat1_sic_bbl ! out
  type(sidl_double_1d) :: alon2_sic_bbl ! out
  type(sidl_double_1d) :: alat2_sic_bbl ! out
  type(sidl_double_1d) :: alon3_sic_bbl ! out
  type(sidl_double_1d) :: alat3_sic_bbl ! out
  type(sidl_double_1d) :: aboxedge1_lon_sic_bbl ! out
  type(sidl_double_1d) :: aboxedge1_lat_sic_bbl ! out
  type(sidl_double_1d) :: aboxedge2_lon_sic_bbl ! out
  type(sidl_double_1d) :: aboxedge2_lat_sic_bbl ! out
  type(sidl_double_1d) :: aboxedge3_lon_sic_bbl ! out
  type(sidl_double_1d) :: aboxedge3_lat_sic_bbl ! out
  type(sidl_int_2d) :: ilandmask1_sic_bbl ! out
  type(sidl_int_2d) :: ilandmask2_sic_bbl ! out
  type(sidl_int_2d) :: ilandmask3_sic_bbl ! out
  integer (kind=sidl_int) :: koverall_total_bbl ! in
  type(sidl_double_2d) :: hght_sic_bbl ! out
  type(sidl_double_2d) :: frac_sic_bbl ! out
  type(sidl_double_2d) :: temp_sic_bbl ! out
  type(sidl_double_2d) :: albd_sic_bbl ! out
  real (kind=sidl_double) :: test_energy_seaice_bbl ! out
  integer (kind=sidl_int) :: retval ! out
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(seaice.goldseaice.initialise)
! Insert-Code-Here {seaice.goldseaice.initialise} (initialise method)
  
  real, dimension(ilon1_sic) :: alon1_sic
  real, dimension(ilat1_sic) :: alat1_sic
  real, dimension(ilon2_sic) :: alon2_sic
  real, dimension(ilat2_sic) :: alat2_sic
  real, dimension(ilon3_sic) :: alon3_sic
  real, dimension(ilat3_sic) :: alat3_sic
  real, dimension(ilon1_sic+1) :: aboxedge1_lon_sic
  real, dimension(ilat1_sic+1) :: aboxedge1_lat_sic
  real, dimension(ilon2_sic+1) :: aboxedge2_lon_sic
  real, dimension(ilat2_sic+1) :: aboxedge2_lat_sic
  real, dimension(ilon3_sic+1) :: aboxedge3_lon_sic
  real, dimension(ilat3_sic+1) :: aboxedge3_lat_sic
  integer, dimension(ilon1_sic,ilat1_sic) :: ilandmask1_sic
  integer, dimension(ilon2_sic,ilat2_sic) :: ilandmask2_sic
  integer, dimension(ilon3_sic,ilat3_sic) :: ilandmask3_sic
  integer :: koverall_total_int
  
  real, dimension(ilon1_sic,ilat1_sic) :: hght_sic
  real, dimension(ilon1_sic,ilat1_sic) :: frac_sic
  real, dimension(ilon1_sic,ilat1_sic) :: temp_sic
  real, dimension(ilon1_sic,ilat1_sic) :: albd_sic
  real :: test_energy_seaice 

  koverall_total_int = koverall_total_bbl

  call initialise_seaice(alon1_sic, alat1_sic, alon2_sic, &
           alat2_sic, alon3_sic, alat3_sic, aboxedge1_lon_sic, &
	   aboxedge1_lat_sic, aboxedge2_lon_sic, aboxedge2_lat_sic, &
	   aboxedge3_lon_sic, aboxedge3_lat_sic, ilandmask1_sic, &
	   ilandmask2_sic, ilandmask3_sic, koverall_total_int, &
	   hght_sic, frac_sic, temp_sic, albd_sic, test_energy_seaice)

  call create1d(ilon1_sic, alon1_sic_bbl)
  call create1d(ilat1_sic, alat1_sic_bbl)
  call create1d(ilon2_sic, alon2_sic_bbl)
  call create1d(ilat2_sic, alat2_sic_bbl)
  call create1d(ilon3_sic, alon3_sic_bbl)
  call create1d(ilat3_sic, alat3_sic_bbl)
  call create1d(ilon1_sic+1, aboxedge1_lon_sic_bbl)
  call create1d(ilat1_sic+1, aboxedge1_lat_sic_bbl)
  call create1d(ilon2_sic+1, aboxedge2_lon_sic_bbl)
  call create1d(ilat2_sic+1, aboxedge2_lat_sic_bbl)
  call create1d(ilon3_sic+1, aboxedge3_lon_sic_bbl)
  call create1d(ilat3_sic+1, aboxedge3_lat_sic_bbl)

  call create2dcol(ilon1_sic, ilat1_sic, ilandmask1_sic_bbl)
  call create2dcol(ilon2_sic, ilat2_sic, ilandmask2_sic_bbl)
  call create2dcol(ilon3_sic, ilat3_sic, ilandmask3_sic_bbl)

  call create2dcol(ilon1_sic, ilat1_sic, hght_sic_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, frac_sic_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, temp_sic_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, albd_sic_bbl)

  call castSimple1DtoSidlDouble1D(alon1_sic, alon1_sic_bbl)
  call castSimple1DtoSidlDouble1D(alat1_sic, alat1_sic_bbl)
  call castSimple1DtoSidlDouble1D(alon2_sic, alon2_sic_bbl)
  call castSimple1DtoSidlDouble1D(alat2_sic, alat2_sic_bbl)
  call castSimple1DtoSidlDouble1D(alon3_sic, alon3_sic_bbl)
  call castSimple1DtoSidlDouble1D(alat3_sic, alat3_sic_bbl)
  call castSimple1DtoSidlDouble1D(aboxedge1_lon_sic, aboxedge1_lon_sic_bbl)
  call castSimple1DtoSidlDouble1D(aboxedge1_lat_sic, aboxedge1_lat_sic_bbl)
  call castSimple1DtoSidlDouble1D(aboxedge2_lon_sic, aboxedge2_lon_sic_bbl)
  call castSimple1DtoSidlDouble1D(aboxedge2_lat_sic, aboxedge2_lat_sic_bbl)
  call castSimple1DtoSidlDouble1D(aboxedge3_lon_sic, aboxedge3_lon_sic_bbl)
  call castSimple1DtoSidlDouble1D(aboxedge3_lat_sic, aboxedge3_lat_sic_bbl)
  call castSimple2DtoSidlInt2D(ilandmask1_sic, ilandmask1_sic_bbl)
  call castSimple2DtoSidlInt2D(ilandmask2_sic, ilandmask2_sic_bbl)
  call castSimple2DtoSidlInt2D(ilandmask3_sic, ilandmask3_sic_bbl)
  
  call castSimple2DtoSidlDouble2D(hght_sic, hght_sic_bbl)
  call castSimple2DtoSidlDouble2D(frac_sic, frac_sic_bbl)
  call castSimple2DtoSidlDouble2D(temp_sic, temp_sic_bbl)
  call castSimple2DtoSidlDouble2D(albd_sic, albd_sic_bbl)
  
  
  test_energy_seaice_bbl = test_energy_seaice
      
! DO-NOT-DELETE splicer.end(seaice.goldseaice.initialise)
end subroutine seaice_goldseaice_initialise_mi


! 
! Method:  run[]
! 

recursive subroutine seaice_goldseaice_run_mi(self, istep_sic_bbl,             &
  dhght_sic_bbl, dfrac_sic_bbl, ustar_ocn_bbl, vstar_ocn_bbl, hght_sic_bbl,    &
  frac_sic_bbl, temp_sic_bbl, albd_sic_bbl, waterflux_ocn_bbl,                 &
  conductflux_ocn_bbl, test_energy_seaice_bbl, test_water_seaice_bbl,          &
  koverall_bbl, retval, exception)
  use sidl
  use sidl_BaseInterface
  use sidl_RuntimeException
  use seaice_goldseaice
  use sidl_double_array
  use seaice_goldseaice_impl
  ! DO-NOT-DELETE splicer.begin(seaice.goldseaice.run.use)
  ! Insert-Code-Here {seaice.goldseaice.run.use} (use statements)
  use genie_control
  use castings
  ! DO-NOT-DELETE splicer.end(seaice.goldseaice.run.use)
  implicit none
  type(seaice_goldseaice_t) :: self ! in
  integer (kind=sidl_int) :: istep_sic_bbl ! in
  type(sidl_double_2d) :: dhght_sic_bbl ! in
  type(sidl_double_2d) :: dfrac_sic_bbl ! in
  type(sidl_double_2d) :: ustar_ocn_bbl ! in
  type(sidl_double_2d) :: vstar_ocn_bbl ! in
  type(sidl_double_2d) :: hght_sic_bbl ! out
  type(sidl_double_2d) :: frac_sic_bbl ! out
  type(sidl_double_2d) :: temp_sic_bbl ! in
  type(sidl_double_2d) :: albd_sic_bbl ! in
  type(sidl_double_2d) :: waterflux_ocn_bbl ! out
  type(sidl_double_2d) :: conductflux_ocn_bbl ! out
  real (kind=sidl_double) :: test_energy_seaice_bbl ! out
  real (kind=sidl_double) :: test_water_seaice_bbl ! out
  integer (kind=sidl_int) :: koverall_bbl ! in
  integer (kind=sidl_int) :: retval ! out
  type(sidl_BaseInterface_t) :: exception ! out

! DO-NOT-DELETE splicer.begin(seaice.goldseaice.run)
! Insert-Code-Here {seaice.goldseaice.run} (run method)
  integer :: istep_sic
  real, dimension(ilon1_sic,ilat1_sic) :: dhght_sic
  real, dimension(ilon1_sic,ilat1_sic) :: dfrac_sic
  real, dimension(ilon1_ocn,ilat1_ocn) :: ustar_ocn
  real, dimension(ilon1_ocn,ilat1_ocn) :: vstar_ocn
  real, dimension(ilon1_sic,ilat1_sic) :: hght_sic
  real, dimension(ilon1_sic,ilat1_sic) :: frac_sic
  real, dimension(ilon1_sic,ilat1_sic) :: temp_sic
  real, dimension(ilon1_sic,ilat1_sic) :: albd_sic
  real, dimension(ilon1_ocn,ilat1_ocn) :: waterflux_ocn
  real, dimension(ilon1_ocn,ilat1_ocn) :: conductflux_ocn
  real :: test_energy_seaice 
  real :: test_water_seaice 
  integer :: koverall      


  istep_sic = istep_sic_bbl

  call castSidlDouble2DtoSimple2D(dhght_sic_bbl, dhght_sic)
  call castSidlDouble2DtoSimple2D(dfrac_sic_bbl, dfrac_sic)
  call castSidlDouble2DtoSimple2D(ustar_ocn_bbl, ustar_ocn)
  call castSidlDouble2DtoSimple2D(vstar_ocn_bbl, vstar_ocn)
  call castSidlDouble2DtoSimple2D(temp_sic_bbl, temp_sic)
  call castSidlDouble2DtoSimple2D(albd_sic_bbl, albd_sic)

  koverall = koverall_bbl
  
  call gold_seaice(istep_sic, dhght_sic, dfrac_sic, ustar_ocn, &
           vstar_ocn, hght_sic, frac_sic, temp_sic, albd_sic, &
           waterflux_ocn, conductflux_ocn, test_energy_seaice, &
           test_water_seaice, koverall)
 
  call create2dcol(ilon1_sic, ilat1_sic, hght_sic_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, frac_sic_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, temp_sic_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, albd_sic_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, waterflux_ocn_bbl)
  call create2dcol(ilon1_sic, ilat1_sic, conductflux_ocn_bbl)

  call castSimple2DtoSidlDouble2D(hght_sic, hght_sic_bbl)
  call castSimple2DtoSidlDouble2D(frac_sic, frac_sic_bbl)
!  call castSimple2DtoSidlDouble2D(temp_sic, temp_sic_bbl)
!  call castSimple2DtoSidlDouble2D(albd_sic, albd_sic_bbl)
  call castSimple2DtoSidlDouble2D(waterflux_ocn, waterflux_ocn_bbl)
  call castSimple2DtoSidlDouble2D(conductflux_ocn, conductflux_ocn_bbl)
  
  test_energy_seaice_bbl = test_energy_seaice
  test_water_seaice_bbl = test_water_seaice

! DO-NOT-DELETE splicer.end(seaice.goldseaice.run)
end subroutine seaice_goldseaice_run_mi


! DO-NOT-DELETE splicer.begin(_miscellaneous_code_end)
! Insert-Code-Here {_miscellaneous_code_end} (extra code)
! DO-NOT-DELETE splicer.end(_miscellaneous_code_end)
