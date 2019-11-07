      program make_bconds

c     THIS CODE MAKES:
c     a file containing some of the boundary conditions necessary to do a 
c       fixed SST simulation with genie.
c     also can be used for slab simulations.
c     It is particularly useful for teh LGM simulations.
c     Make sure that makefile.arc is set to the same as the 
c       compiler precision for this code

c     I Compile with: 
c     ifort -o make_bconds make_bconds.f -lnc1 -lutil1 -lnetcdf -ligcm3 -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-lib/libutil1 -L/home/ggdjl/genie/genie-igcm3/lib 

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
#include "/home/ggdjl/genie/genie-igcm3/src/fortran/param1.cmn"
      include "/home/ggdjl/genie/genie-igcm3/src/fortran/param2.cmn"

c***********************************************************

      include 'names_lgm.inc'
c      include 'names_cnt.inc'

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
      real seaice_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real relief_1_igcm(ilon1_atm,ilat1_atm)
      real landmask_1_igcm(ilon1_atm,ilat1_atm)
      real landice_1_igcm(ilon1_atm,ilat1_atm)
      real vege_1_igcm(ilon1_atm,ilat1_atm)

      integer work_igcm(ilon1_atm,ilat1_atm)

   

c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a

c     LOOPING:
      integer i,j,m,v


c     FOR GRID:
      real ax


      print*,'make_lgm is starting....'

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

      print*,'igcm grid set up OK'


c     set up work arrays (one everywhere)

      work_igcm(:,:)=1.

      print*,'work arrays initialised'

      call read_bconds(2,12,1,sst_12_igcm,in_sst_filename,
     &      in_sst_lonname,in_sst_latname,in_sst_varname)

      call read_bconds(2,12,1,seaice_12_igcm,in_seaice_filename,
     &      in_seaice_lonname,in_seaice_latname,in_seaice_varname)
 
      call read_bconds(1,1,2,relief_1_igcm,in_relief_filename,
     &      in_relief_lonname,in_relief_latname,in_relief_varname)

      call read_bconds(1,1,3,landmask_1_igcm,in_relief_filename,
     &      in_relief_lonname,in_relief_latname,in_relief_varname)

      call read_bconds(1,1,4,vege_1_igcm,in_vege_filename,
     &      in_vege_lonname,in_vege_latname,in_vege_varname)

      call read_bconds(1,1,1,landice_1_igcm,in_landice_filename,
     &      in_landice_lonname,in_landice_latname,in_landice_varname)

      print*,'all data read in'


c     process the data...............

      do j=1,ilat1_atm
      do i=1,ilon1_atm
   

c     Make the orography and land-sea mask
      if (landmask_1_igcm(i,j).ge.0.5) then
         landmask_1_igcm(i,j)=1.0
         relief_1_igcm(i,j)=relief_1_igcm(i,j)
      else
         landmask_1_igcm(i,j)=0.0
         relief_1_igcm(i,j)=0.0
      endif

c     make sure landice is zero if landmask is zero
      landice_1_igcm(i,j)=landice_1_igcm(i,j)*landmask_1_igcm(i,j)

c     make the landice mask binary
      if (landice_1_igcm(i,j).gt.0.3) then
        landice_1_igcm(i,j)=1.0
      else
        landice_1_igcm(i,j)=0.0
      endif

c     replace landice with bare soil
      if (landice_1_igcm(i,j).eq.1) then
        vege_1_igcm(i,j)=24.0
      endif

      if (landmask_1_igcm(i,j).eq.0) then
      vege_1_igcm(i,j)=1.0
      endif

      if ((landmask_1_igcm(i,j).eq.1).and.
     :          (vege_1_igcm(i,j).eq.1.0)) then
      if (vege_1_igcm(i+1,j).ne.1.0) 
     :        vege_1_igcm(i,j)=vege_1_igcm(i+1,j)
      if (vege_1_igcm(i,j+1).ne.1.0) 
     :        vege_1_igcm(i,j)=vege_1_igcm(i,j+1)
      if (vege_1_igcm(i-1,j).ne.1.0) 
     :        vege_1_igcm(i,j)=vege_1_igcm(i-1,j)
      if (vege_1_igcm(i,j-1).ne.1.0) 
     :        vege_1_igcm(i,j)=vege_1_igcm(i,j-1)
      endif 

      if ((landmask_1_igcm(i,j).eq.1).and.
     :          (vege_1_igcm(i,j).eq.1.0)) then
      print*,'Land but no veg at.....',i,j
      stop
      endif


      if ((landice_1_igcm(i,j).eq.1).and.
     :          (vege_1_igcm(i,j).ne.24)) then
      print*,'Ice but not soil at.....',i,j
      stop
      endif


      if (vege_1_igcm(i,j).eq.2) then
      print*,'Ice in vegetation at',i,j
      stop
      endif


      do m=1,nmonths

      sst_12_igcm(i,j,m)=sst_12_igcm(i,j,m)-273.15

      if (seaice_12_igcm(i,j,m).ge.50) then
      sst_12_igcm(i,j,m)=-2.1
      seaice_12_igcm(i,j,m)=1.0
      else
      seaice_12_igcm(i,j,m)=0.0
      endif

      if ( (sst_12_igcm(i,j,m).le.-2.0).and.
     &             (seaice_12_igcm(i,j,m).ne.1.0) ) then
        sst_12_igcm(i,j,m)=-1.9
      endif      


      enddo      


      enddo
      enddo




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


      nvar=6
      v=1

      varname(v,1)='sst'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='tstar from CLIMAP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='degrees C'

      v=2

      varname(v,1)='seaice'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='sea-ice fraction from CLIMAP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=3

      varname(v,1)='orog'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='height'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='metres'

      v=4

      varname(v,1)='lsm'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='land/sea mask'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='1/0'

      v=5

      varname(v,1)='svege'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='vegetation type'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='1 to 26'

      v=6

      varname(v,1)='icefrac'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='ice sheet fraction'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0/1'


      call ininc(out_data_filename,
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
  

      v=2

      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        seaice_12_igcm )
      print*,'written file'


      v=3

      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        relief_1_igcm )
      print*,'written file'

      v=4

      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        landmask_1_igcm )
      print*,'written file'


      v=5

      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        vege_1_igcm )
      print*,'written file'

      v=6

      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        landice_1_igcm )
      print*,'written file'

      call end_netcdf(1)



      print*,'files written OK'



      end


      include 'read_bconds.f'
