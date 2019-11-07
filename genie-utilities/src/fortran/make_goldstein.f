      program make_goldstein

c     THIS CODE MAKES:
c     Peltier land-sea masks, topography, icesheet masks, and vegetation,
c       for 21 snapshots (presently ice-4g)
c     Currently, ice-sheets are replaced with 0kyr vegetation, bare soil, or the 
c       nearest vegetation-type
c     The Caspian sea is an inland lake.  Its size does not change with time.


c     I Compile with: 
c     ifort -o make_goldstein_NOCVS make_goldstein.f -lnc1 -lutil1 -lnetcdf -ligcm3 -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-lib/libutil1 -L/home/ggdjl/genie/genie-igcm3/lib -L/opt/local/intel_fc_81/lib -fpp 

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
#include "/home/ggdjl/genie/genie-igcm3/src/fortran/param1.cmn"
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/param2.cmn'

c***********************************************************

c      include 'names_goldstein.inc'
      include 'names_goldstein_6432.inc'

c***********************************************************

c     TIMES:

c     IGCM GRID:
      integer ilon1_atm,ilat1_atm
      parameter (ilon1_atm=64,ilat1_atm=32)
      real alon1(ilon1_atm),alat1(ilat1_atm),
     :     aboxedge1_lon(ilon1_atm+1),aboxedge1_lat(ilat1_atm+1)

c     DATA:
      real relief_igcm(ilon1_atm,ilat1_atm)
      real landmask_igcm(ilon1_atm,ilat1_atm)
      real landice_igcm(ilon1_atm,ilat1_atm)
      real vege_igcm(ilon1_atm,ilat1_atm)

c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a
c     From old netdata.cmn:
      integer nmaxdims
      parameter(nmaxdims=4)
      integer ndim,nvar,natts(nall),nattsvar(nall),
     :        vdims(nall),vadims(nmaxdims,nall),   
     :        ndims(nall)   
      character   
     :          attdimname(2,nmaxdims,nall)*200,   
     :          attvarname(2,nmaxdims,nall)*200 


c     LOOPING:
      integer i,j,t,v

c     FOR GRID:
      real ax

c     FOR FILE:
      character filename*200

      print*,'make_goldstein is starting....'


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

      print*,'work arrays initialised'

      filename=trim(in_relief_filename)//'_00.nc'
      call read_bconds(1,1,2,relief_igcm,
     &      filename,
     &      in_relief_lonname,in_relief_latname,in_relief_varname)

      filename=trim(in_landmask_filename)//'.nc'
      call read_bconds(1,1,5,landmask_igcm,
     &      filename,
     &      in_landmask_lonname,in_landmask_latname,
     &      in_landmask_varname)

      filename=trim(in_landice_filename)//'_00.nc'
      call read_bconds(1,1,1,landice_igcm,
     &      filename,
     &      in_landice_lonname,in_landice_latname,in_landice_varname)

      call read_bconds(1,1,4,vege_igcm,in_vege_filename,
     &      in_vege_lonname,in_vege_latname,in_vege_varname)

      print*,'all data read in'

c     process the data...............


      do j=1,ilat1_atm
      do i=1,ilon1_atm


c     Make the orography and land-sea mask
      if (landmask_igcm(i,j).ge.0.5) then
         landmask_igcm(i,j)=1.0
      else
         landmask_igcm(i,j)=0.0
         relief_igcm(i,j)=0.0
      endif

c     make orog at least 5m high if we have land 
      if (landmask_igcm(i,j).eq.1) then
        relief_igcm(i,j)=max(5.0,relief_igcm(i,j))
      endif

c     make icesheet south of 65, north of 80
      if (landmask_igcm(i,j).eq.1) then
        if ((alat1(j).lt.-60).or.(alat1(j).gt.80)) then
          landice_igcm(i,j)=1.0
        endif
      endif

c     make the landice mask binary
      if (landice_igcm(i,j).gt.0.3) then
        landice_igcm(i,j)=1.0
      else
        landice_igcm(i,j)=0.0
      endif

c     make se tip of greenland icesheet
      if ((i.eq.59).and.(j.eq.5)) then
        landice_igcm(i,j)=1
      endif

c     make sure landice is zero if landmask is zero
      landice_igcm(i,j)=landice_igcm(i,j)*
     :               landmask_igcm(i,j)

c     replace 0kyr landice with bare soil
      if (landice_igcm(i,j).eq.1) then
        vege_igcm(i,j)=24.0
      endif

c     get rid of any ice in vegetation 
      if (vege_igcm(i,j).eq.2) then
        vege_igcm(i,j)=24.0
      endif

      enddo
      enddo


      do j=1,ilat1_atm
      do i=1,ilon1_atm

c     expand vegetation if needed
      if ((landmask_igcm(i,j).eq.1).and.
     :          (vege_igcm(i,j).eq.1.0)) then

      if (vege_igcm(i+1,j).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i+1,j)
      if (vege_igcm(i,j+1).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i,j+1)
      if (vege_igcm(i-1,j).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i-1,j)
      if (vege_igcm(i,j-1).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i,j-1)

      if (vege_igcm(i,j).eq.1.0) then
      if (vege_igcm(i+1,j+1).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i+1,j+1)
      if (vege_igcm(i+1,j-1).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i+1,j-1)
      if (vege_igcm(i-1,j+1).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i-1,j+1)
      if (vege_igcm(i-1,j-1).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i-1,j-1)
      endif

      if (vege_igcm(i,j).eq.1.0) then
      if (vege_igcm(i+2,j).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i+2,j)
      if (vege_igcm(i,j+2).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i,j+2)
      if (vege_igcm(i-2,j).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i-2,j)
      if (vege_igcm(i,j-2).ne.1.0)
     :        vege_igcm(i,j)=vege_igcm(i,j-2)
      endif

      endif 

      enddo
      enddo



      do j=1,ilat1_atm
      do i=1,ilon1_atm

c     make sure no points with no veg
      if ((landmask_igcm(i,j).eq.1).and.
     :          (vege_igcm(i,j).eq.1.0)) then
      print*,'Land but no veg at.....',i,j
      stop
      endif

c     make sure 0kyr ice has bare soil below
      if ((landice_igcm(i,j).eq.1).and.
     :          (vege_igcm(i,j).ne.24)) then
      print*,'0kyrBP Ice but not soil at.....',i,j
      stop
      endif

c     Make sure no ice in vegetation
      if (vege_igcm(i,j).eq.2) then
      print*,'Ice in vegetation at',i,j
      stop
      endif

c     Have a look at inland lakes 
      if (vege_igcm(i,j).eq.3) then
      print*,'Inland lake at',i,j
c      stop
      endif

      enddo
      enddo

      print*,'post-processing complete'

      print*,'writing netcdf output'

      ndim=2
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

      nvar=4
      v=1
      varname(v,1)='orog'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='height'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='metres'

      v=v+1
      varname(v,1)='lsm'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='land/sea mask'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='1/0'

      v=v+1
      varname(v,1)='icefrac'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='ice sheet fraction'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0/1'

      v=v+1
      varname(v,1)='svege'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='vegetation type'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='1 to 26'

      call ininc(out_data_filename,
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))

      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)

      v=1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        relief_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        landmask_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        landice_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        vege_igcm )
      print*,'written file'

      call end_netcdf(1)


      print*,'files written OK'


      end


#include 'read_bconds.f'
