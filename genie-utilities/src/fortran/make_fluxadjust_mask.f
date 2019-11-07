      program make_fluxadjust_mask

c***********************************************************

c     This code is to test the pacific - atlantic freshwater flux
c       in a genie run.  It is not guarenteed at present as gives
c       odd results in a long gawater run.  This needs testing (DJL).

c     I compile with: 
c     ifort -o make_fluxadjust_mask_NOCVS make_fluxadjust_mask.f -lnc1 -lutil1 -lnetcdf -ligcm3 -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-lib/libutil1 -L/opt/local/intel_fc_81/lib -L/home/ggdjl/genie/genie-igcm3/lib -fpp

c     The lines between the asterixes will need to be edited.

c     Please edit the main body of the code also, to tidy it up!

c     Eventually, this piece of code could run with makefile.arc!

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
#include '/home/ggdjl/genie/genie-igcm3/src/fortran/param1.cmn'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/param2.cmn'

      include 'names_fluxadjust.inc'

c***********************************************************

c     MONTHS:

c     IGCM GRID:
      integer ilon1_atm,ilat1_atm
      parameter (ilon1_atm=64,ilat1_atm=32)
      real alon1(ilon1_atm),alat1(ilat1_atm),
     :     aboxedge1_lon(ilon1_atm+1),aboxedge1_lat(ilat1_atm+1)

c     DATA:

      real precip_1_igcm(ilon1_atm,ilat1_atm)
      real evap_1_igcm(ilon1_atm,ilat1_atm)
      real runoff_1_igcm(ilon1_atm,ilat1_atm)
      real lsm_1_igcm(ilon1_atm,ilat1_atm)

      real datain_1_igcm(ilon1_atm,ilat1_atm)


c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a

c     LOOPING:
      integer i,j,v

c     FOR GRID:
      real ax

c     For weights:
      real weight_atm(ilon1_atm,ilat1_atm)
      real weightcheck
      real plumin

c     For transport diagnostics
      real at1,at2,at3
      real at_mask(ilon1_atm,ilat1_atm)
      real at_mask_region(ilon1_atm,ilat1_atm)
      real radea
      parameter(radea=6370e3)
      
c     From old netdata.cmn:
      integer nmaxdims
      parameter(nmaxdims=4)
      integer ndim,nvar,natts(nall),nattsvar(nall),
     :        vdims(nall),vadims(nmaxdims,nall),   
     :        ndims(nall)   
      character   
     :          attdimname(2,nmaxdims,nall)*200,   
     :          attvarname(2,nmaxdims,nall)*200 



      print*,'make_fixed_ncep is starting....'

c     COPIED FROM INITIALISE_ATMOS.F
      ax=360.0/real(mg)
      do i=1,ilon1_atm
         aboxedge1_lon(i)=(i-1.5)*ax
         alon1(i)=(i-1.0)*ax
      end do
      aboxedge1_lon(mg+1)=(mg-0.5)*ax
      call gwtcnr(alat1,jg)
      call gwtbox(aboxedge1_lat,jg)

      if (aboxedge1_lat(1).gt.aboxedge1_lat(ilat1_atm)) then
        plumin=1.0
      else
        plumin=-1.0
      endif
      weightcheck=0.0
      do j=1,ilat1_atm
        do i=1,ilon1_atm
          weight_atm(i,j)=plumin*(sin(aboxedge1_lat(j)*2*pi/360.)-
     &              sin(aboxedge1_lat(j+1)*2*pi/360.))*
     &              (aboxedge1_lon(i+1)/360.-
     &              aboxedge1_lon(i)/360.)/2.
          weightcheck=weightcheck+weight_atm(i,j)
        enddo
      enddo
      print*,'Check for weightings = ',weightcheck

      print*,'igcm grid set up OK in main'


c     set up work arrays (one everywhere)

      print*,'work arrays initialised'


      call open_file_nc(trim(in_data_maskname),ncid)
      call get2d_data_nc(ncid,trim(in_mask_varname),ilon1_atm,ilat1_atm,
     &          datain_1_igcm,ifail)
      if (ifail.ne.0) stop 'Problem reading input mask file'
      lsm_1_igcm(:,:)=datain_1_igcm(:,:)
      call close_file_nc(trim(in_data_maskname),ncid)

      print*,'all data read in'

c     process the data...............


      print*,'calculating diagnostics'


c     Now calculate some diagnostics:

      do j=1,ilat1_atm
        do i=1,ilon1_atm
           if ( (alon1(i).gt.290).or.(alon1(i).lt.20) ) 
     :             at_mask(i,j)=-1
           if (alat1(j).lt.-50)
     :             at_mask(i,j)=0
           if (alat1(j).gt.65)
     :             at_mask(i,j)=-1
           if ( (alon1(i).gt.230).and.(alat1(j).gt.50) ) 
     :             at_mask(i,j)=-1
           if ( (alon1(i).gt.265).and.(alat1(j).gt.10) ) 
     :             at_mask(i,j)=-1
           if ( (alon1(i).lt.50).and.(alat1(j).gt.30) ) 
     :             at_mask(i,j)=-1
            if ( (alon1(i).gt.260).and.(alat1(j).gt.20) ) 
     :             at_mask(i,j)=-1
           at_mask(i,j)=at_mask(i,j)*(1-lsm_1_igcm(i,j))
        enddo
      enddo

      do j=1,ilat1_atm
        do i=1,ilon1_atm
          if (at_mask(i,j).eq.0) then
           if (alon1(i).gt.120)
     :             at_mask(i,j)=1
           if (alat1(j).lt.-50)
     :             at_mask(i,j)=0
          endif
           at_mask(i,j)=at_mask(i,j)*(1-lsm_1_igcm(i,j))
        enddo
      enddo



      at_mask_region(:,:)=0.0
      do j=1,ilat1_atm
        do i=1,ilon1_atm
          if (alat1(j).gt.30) then
            at_mask_region(i,j)=at_mask(i,j)*1
          else if (alat1(j).lt.-30) then
            at_mask_region(i,j)=at_mask(i,j)*3
          else
            at_mask_region(i,j)=at_mask(i,j)*2
          endif
        enddo
      enddo

c     12 is for the number of months, 1e12 is for Sverdrups and kg =>m3
      at1=at1*4.0*pi*radea*radea/(1e9)
      at2=at2*4.0*pi*radea*radea/(1e9)
      at3=at3*4.0*pi*radea*radea/(1e9)

      print*,'diagnostics complete'

      print*,'Diagnostics:'
      print*,'Atlantic 1',at1
      print*,'Atlantic 2',at2
      print*,'Atlantic 3',at3

      print*,pi

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

      nvar=2

      v=1
      varname(v,1)='weight'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='weight of gridbox'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='-'

      v=2
      varname(v,1)='water_mask'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='region number'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='-'

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
     :        weight_atm )
      print*,'written file'
  
      v=2
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        at_mask_region )
      print*,'written file'

      call end_netcdf(1)



      print*,'files written OK'


      end


#include 'read_bconds.f'
