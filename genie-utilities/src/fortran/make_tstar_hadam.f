      program make_tstar_hadam

c***********************************************************

c     This code can be used to create a temperature field, 
c       identical to one from hadam3.

c     I compile with: 
c     ifort -o make_tstar_hadam_NOCVS make_tstar_hadam.f -lnc1 -lutil1 -lnetcdf -ligcm3 -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-lib/libutil1 -L/home/ggdjl/genie/genie-igcm3/lib -fpp

c     The lines between the asterixes will need to be edited.

c     Please edit the main body of the code also, to tidy it up!

c     Eventually, this piece of code could run with makefile.arc!

c     Before running this code, the hadam results are concatenated using
c     the make_hadam_concat script.


      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
#include '/home/ggdjl/genie/genie-igcm3/src/fortran/param1.cmn'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/param2.cmn'

      include 'names_hadam_mod.inc'
c      include 'names_hadam_pre.inc'

c***********************************************************

c     MONTHS:
      integer nmonths
      parameter(nmonths=12)
      real months(nmonths)

c     IGCM GRID:
      integer ilon1_atm,ilat1_atm
      parameter (ilon1_atm=64,ilat1_atm=32)
      real alon1(ilon1_atm),alat1(ilat1_atm),
     :     aboxedge1_lon(ilon1_atm+1),aboxedge1_lat(ilat1_atm+1)

c     DATA:

      real sst_12_igcm(ilon1_atm,ilat1_atm,nmonths)

      integer work_igcm(ilon1_atm,ilat1_atm)

      real datain_1_igcm(ilon1_atm,ilat1_atm)

c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a

c     LOOPING:
      integer i,j,m,v


c     FOR GRID:
      real ax

      print*,'make_tstar_hadam is starting....'

      do m=1,nmonths
       months(m)=m
      enddo


c     COPIED FROM INITIALISE_ATMOS.F
      ax=360.0/real(mg)
      do i=1,ilon1_atm
         aboxedge1_lon(i)=(i-1.5)*ax
         alon1(i)=(i-1.0)*ax
      end do
      aboxedge1_lon(mg+1)=(mg-0.5)*ax
      call gwtcnr(alat1,jg)
      call gwtbox(aboxedge1_lat,jg)

      print*,'igcm grid set up OK in main'


c     set up work arrays (one everywhere)

      work_igcm(:,:)=1.

      print*,'work arrays initialised'

      call read_bconds(2,12,1,sst_12_igcm,in_sst_filename,
     &      in_sst_lonname,in_sst_latname,in_sst_varname)

      print*,'all data read in'


c     process the data...............

      sst_12_igcm(:,:,:)=sst_12_igcm(:,:,:)-273.16

      print*,'post-processing complete'

      print*,'writing netcdf output'

      ndim=3
      dimname(1,1)='longitude'
      ndims(1)=ilon1_atm
      natts(1)=2
      attdimname(1,1,1)='long_name'
      attdimname(2,1,1)='longitude'
      attdimname(1,2,1)='units'
      attdimname(2,2,1)='degrees east'

      dimname(2,1)='latitude'
      ndims(2)=ilat1_atm
      natts(2)=2
      attdimname(1,1,2)='long_name'
      attdimname(2,1,2)='latitude'
      attdimname(1,2,2)='units'
      attdimname(2,2,2)='degrees north'

      dimname(3,1)='time'
      ndims(3)=nmonths
      natts(3)=2
      attdimname(1,1,3)='long_name'
      attdimname(2,1,3)='Month Number'
      attdimname(1,2,3)='units'
      attdimname(2,2,3)='seconds'


      nvar=1

      v=1
      varname(v,1)='sst'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='SST from AMIP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='oC'


      call ininc(out_sst_filename,
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))

      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)
      call writedim(nc(1),iddim(3,1),months)

      v=1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        sst_12_igcm )
      print*,'written file'
  

      call end_netcdf(1)



      print*,'files written OK'


      end


#include 'read_bconds.f'
