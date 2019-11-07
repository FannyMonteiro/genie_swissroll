      program make_peltier

c     THIS CODE MAKES:
c     Peltier land-sea masks, topography, icesheet masks, and vegetation,
c       for 21 snapshots (presently ice-4g)
c     Currently, ice-sheets are replaced with 0kyr vegetation, bare soil, or the 
c       nearest vegetation-type
c     The Caspian sea is an inland lake.  It's size does not change with time.


c     I Compile with: 
c     ifort -o make_peltier_NOCVS make_peltier.f -lnc1 -lutil1 -lnetcdf -ligcm3 -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-lib/libutil1 -L/opt/local/intel_fc_81/lib -L/home/ggdjl/genie/genie-igcm3/lib -fpp

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
#include "/home/ggdjl/genie/genie-igcm3/src/fortran/param1.cmn"
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/param2.cmn'

c***********************************************************

      include 'names_pelt.inc'

c***********************************************************

c     TIMES:
      integer ntimes
      parameter(ntimes=22)
      real times(ntimes)
      integer itimes(ntimes)
      character(len=2), dimension(ntimes) :: year

c     IGCM GRID:
      integer ilon1_atm,ilat1_atm
      parameter (ilon1_atm=64,ilat1_atm=32)
      real alon1(ilon1_atm),alat1(ilat1_atm),
     :     aboxedge1_lon(ilon1_atm+1),aboxedge1_lat(ilat1_atm+1)
      integer ilon1hr_atm,ilat1hr_atm
      parameter (ilon1hr_atm=320,ilat1hr_atm=160)
      real alon1hr(ilon1hr_atm),alat1hr(ilat1hr_atm),
     :     aboxedge1hr_lon(ilon1hr_atm+1),
     :     aboxedge1hr_lat(ilat1hr_atm+1)

c     DATA:
      real reliefhr_22_igcm(ilon1hr_atm,ilat1hr_atm,ntimes)
      real relief_22_igcm(ilon1_atm,ilat1_atm,ntimes)
      real landmask_22_igcm(ilon1_atm,ilat1_atm,ntimes)
      real landice_22_igcm(ilon1_atm,ilat1_atm,ntimes)
      real vege_1_igcm(ilon1_atm,ilat1_atm)

      integer work_igcm(ilon1_atm,ilat1_atm)

c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a

c     LOOPING:
      integer i,j,t,v

c     FOR GRID:
      real ax

c     FOR FILE:
      character filename*200

      print*,'make_peltier is starting....'


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

c     Set up times
      do t=1,ntimes
       times(t)=-1*ntimes+t
       itimes(t)=ntimes-t
      write(year(t),'(i2.2)') itimes(t)
      print*,year(t)
      enddo

c     set up work arrays (one everywhere)

      work_igcm(:,:)=1.

      print*,'work arrays initialised'

      do t=1,ntimes

      filename=trim(in_relief_filename)//'_'//year(t)//'.nc'
      call read_bcondshr(1,1,2,relief_22_igcm(:,:,t),
     &      reliefhr_22_igcm(:,:,t),
     &      filename,
     &      in_relief_lonname,in_relief_latname,in_relief_varname)

      filename=trim(in_relief_filename)//'_'//year(t)//'.nc'
      call read_bconds(1,1,3,landmask_22_igcm(:,:,t),
     &      filename,
     &      in_relief_lonname,in_relief_latname,in_relief_varname)

      filename=trim(in_landice_filename)//'_'//year(t)//'.nc'
      call read_bconds(1,1,1,landice_22_igcm(:,:,t),
     &      filename,
     &      in_landice_lonname,in_landice_latname,in_landice_varname)

      enddo

      call read_bconds(1,1,4,vege_1_igcm,in_vege_filename,
     &      in_vege_lonname,in_vege_latname,in_vege_varname)

      print*,'all data read in'

c     process the data...............


      do j=1,ilat1_atm
      do i=1,ilon1_atm

      do t=1,ntimes
   
c     Make the orography and land-sea mask

c     ------------------------
c     OK, make sure that we have the same land-sea mass as at the present.....
c      if (landmask_22_igcm(i,j,t).ge.0.5) then
      if (landmask_22_igcm(i,j,ntimes).ge.0.5) then
c     ------------------------
         landmask_22_igcm(i,j,t)=1.0
      else
         landmask_22_igcm(i,j,t)=0.0
         relief_22_igcm(i,j,t)=0.0
      endif

c     Sort out the Caspian land-sea mask
      landmask_22_igcm(10,9,t)=1
      landmask_22_igcm(10,8,t)=1

c     make the landice mask binary
      if (landice_22_igcm(i,j,t).gt.0.3) then
        landice_22_igcm(i,j,t)=1.0
      else
        landice_22_igcm(i,j,t)=0.0
      endif

c     make sure landice is zero if landmask is zero
      landice_22_igcm(i,j,t)=landice_22_igcm(i,j,t)*
     :               landmask_22_igcm(i,j,t)


      enddo

c     Sort out the Caspian 'vegetation'
      vege_1_igcm(10,8)=24.0
      vege_1_igcm(10,9)=3.0

c     replace 0kyr landice with bare soil
      if (landice_22_igcm(i,j,22).eq.1) then
        vege_1_igcm(i,j)=24.0
      endif

c     get rid of any ice in vegetation 
      if (vege_1_igcm(i,j).eq.2) then
        vege_1_igcm(i,j)=24.0
      endif

      enddo
      enddo


      do j=1,ilat1_atm
      do i=1,ilon1_atm

c     expand vegetation to 21kyrBP values
      if ((landmask_22_igcm(i,j,1).eq.1).and.
     :          (vege_1_igcm(i,j).eq.1.0)) then

      if (vege_1_igcm(i+1,j).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i+1,j)
      if (vege_1_igcm(i,j+1).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i,j+1)
      if (vege_1_igcm(i-1,j).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i-1,j)
      if (vege_1_igcm(i,j-1).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i,j-1)

      if (vege_1_igcm(i,j).eq.1.0) then
      if (vege_1_igcm(i+1,j+1).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i+1,j+1)
      if (vege_1_igcm(i+1,j-1).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i+1,j-1)
      if (vege_1_igcm(i-1,j+1).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i-1,j+1)
      if (vege_1_igcm(i-1,j-1).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i-1,j-1)
      endif

      if (vege_1_igcm(i,j).eq.1.0) then
      if (vege_1_igcm(i+2,j).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i+2,j)
      if (vege_1_igcm(i,j+2).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i,j+2)
      if (vege_1_igcm(i-2,j).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i-2,j)
      if (vege_1_igcm(i,j-2).ne.1.0)
     :        vege_1_igcm(i,j)=vege_1_igcm(i,j-2)
      endif

      endif 

      enddo
      enddo


      do j=1,ilat1_atm
      do i=1,ilon1_atm

      do t=1,ntimes

c     make sure no points with no veg
      if ((landmask_22_igcm(i,j,t).eq.1).and.
     :          (vege_1_igcm(i,j).eq.1.0)) then
      print*,'Land but no veg at.....',i,j,t
      stop
      endif

      enddo

c     make sure 0kyr ice has bare soil below
      if ((landice_22_igcm(i,j,ntimes).eq.1).and.
     :          (vege_1_igcm(i,j).ne.24)) then
      print*,'0kyrBP Ice but not soil at.....',i,j,t
      stop
      endif

c     Make sure no ice in vegetation
      if (vege_1_igcm(i,j).eq.2) then
      print*,'Ice in vegetation at',i,j
      stop
      endif

c     Have a look at inland lakes 
      if (vege_1_igcm(i,j).eq.3) then
      print*,'Inland lake at',i,j
c      stop
      endif

      enddo
      enddo

      print*,'post-processing complete'

      print*,'writing netcdf output'

      ndim=6

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
      ndims(3)=ntimes
      natts(3)=2
      attdimname(1,1,3)='long_name'
      attdimname(2,1,3)='Time number'
      attdimname(1,2,3)='units'
      attdimname(2,2,3)='kyr'

      dimname(4,1)='nrecs'
      ndims(4)=1
      natts(4)=0

      dimname(5,1)='longitude_hr'
      ndims(5)=ilon1hr_atm
      natts(5)=2
      attdimname(1,1,5)='long_name'
      attdimname(2,1,5)='longitude'
      attdimname(1,2,5)='units'
      attdimname(2,2,5)='degrees east'

      dimname(6,1)='latitude_hr'
      ndims(6)=ilat1hr_atm
      natts(6)=2
      attdimname(1,1,6)='long_name'
      attdimname(2,1,6)='latitude'
      attdimname(1,2,6)='units'
      attdimname(2,2,6)='degrees north'

      nvar=6

      v=1
      varname(v,1)='orog'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='height'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='metres'

      v=v+1
      varname(v,1)='lsm'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='land/sea mask'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='1/0'

      v=v+1
      varname(v,1)='icefrac'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
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

      v=v+1
      varname(v,1)='ntimes'
      vdims(v)=1
      vadims(1,v)=loc_dim('nrecs',dimname,nall)
      nattsvar(v)=0

      v=v+1
      varname(v,1)='oroghr'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude_hr',dimname,nall)
      vadims(2,v)=loc_dim('latitude_hr',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='height'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='metres'

      call ininc(out_data_filename,
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))

      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)
      call writedim(nc(1),iddim(3,1),times)
      call writedim(nc(1),iddim(4,1),1.0)

      v=1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        relief_22_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        landmask_22_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        landice_22_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        vege_1_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        real(ntimes) )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        reliefhr_22_igcm )
      print*,'written file'

      call end_netcdf(1)


      print*,'files written OK'


      end

#include 'read_bconds.f'
#include 'read_bcondshr.f'
