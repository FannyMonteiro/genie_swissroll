      program make_fixed_ncep

c***********************************************************

c     This code can be used to create a forcing field for the
c       genie-fixedatmos module.
c     The forcing field s derived from NCEP reanalyses.

c     I compile with: 
c     ifort -o make_fixed_ncep_NOCVS make_fixed_ncep.f -lnc1 -lutil1 -lnetcdf -ligcm3 -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-lib/libutil1 -L/ -L/home/ggdjl/genie/genie-igcm3/lib -fpp

c     The lines between the asterixes will need to be edited.

c     Please edit the main body of the code also, to tidy it up!

c     Eventually, this piece of code could run with makefile.arc!

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
#include '/home/ggdjl/genie/genie-igcm3/src/fortran/param1.cmn'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/param2.cmn'

      include 'names_ncep.inc' 

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

      real precip_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real sensible_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real latent_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real evap_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real netsolar_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real netlong_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real stressx_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real stressy_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real runoff1_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real runoff2_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real runoff3_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real lsm_1_igcm(ilon1_atm,ilat1_atm)
      real index_2_igcm(ilon1_atm,ilat1_atm,2)

      integer work_igcm(ilon1_atm,ilat1_atm)

      real datain_1_igcm(ilon1_atm,ilat1_atm)
      real datain_2_igcm(ilon1_atm,ilat1_atm,2)

c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a

c     LOOPING:
      integer i,j,m,v

c     FOR GRID:
      real ax

c     For weights:
      real weight_atm(ilon1_atm,ilat1_atm)
      real weightcheck
      real plumin

c     For transport diagnostics
      real at1,at2,at3
      real at_mask(ilon1_atm,ilat1_atm)
      real radea
      parameter(radea=6370e3)
      

      print*,'make_fixed_ncep is starting....'

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

      work_igcm(:,:)=1.

      print*,'work arrays initialised'

      call read_bconds(2,12,1,precip_12_igcm,in_precip_filename,
     &      in_precip_lonname,in_precip_latname,in_precip_varname)

       call read_bconds(2,12,1,netsolar_12_igcm,in_netsolar_filename,
     &      in_netsolar_lonname,in_netsolar_latname,in_netsolar_varname)

       call read_bconds(2,12,1,netlong_12_igcm,in_netlong_filename,
     &      in_netlong_lonname,in_netlong_latname,in_netlong_varname)

       call read_bconds(2,12,1,latent_12_igcm,in_latent_filename,
     &      in_latent_lonname,in_latent_latname,in_latent_varname)

       call read_bconds(2,12,1,sensible_12_igcm,in_sensible_filename,
     &      in_sensible_lonname,in_sensible_latname,in_sensible_varname)

       call read_bconds(2,12,1,stressx_12_igcm,in_stressx_filename,
     &      in_stressx_lonname,in_stressx_latname,in_stressx_varname)

       call read_bconds(2,12,1,stressy_12_igcm,in_stressy_filename,
     &      in_stressy_lonname,in_stressy_latname,in_stressy_varname)

       call read_bconds(2,12,6,runoff3_12_igcm,in_runoff_filename,
     &      in_runoff_lonname,in_runoff_latname,in_runoff_varname)

      call open_file_nc(trim(in_data_maskname),ncid)
      call get2d_data_nc(ncid,trim(in_mask_varname),ilon1_atm,ilat1_atm,
     &          datain_1_igcm,ifail)
      if (ifail.ne.0) stop 'Problem reading input mask file'
      lsm_1_igcm(:,:)=datain_1_igcm(:,:)
      call close_file_nc(trim(in_data_maskname),ncid)

      call open_file_nc(trim(in_data_indexname),ncid)
      call get3d_data_nc(ncid,trim(in_index_varname),
     &          ilon1_atm,ilat1_atm,2,
     &          datain_2_igcm,ifail)
      if (ifail.ne.0) stop 'Problem reading input mask file'
      index_2_igcm(:,:,:)=datain_2_igcm(:,:,:)
      call close_file_nc(trim(in_data_indexname),ncid)

      print*,'all data read in'

c     process the data...............

      netsolar_12_igcm(:,:,:)= netsolar_12_igcm(:,:,:)*(-1.0)
      sensible_12_igcm(:,:,:)= sensible_12_igcm(:,:,:)*(-1.0)
      latent_12_igcm(:,:,:)= latent_12_igcm(:,:,:)*(-1.0)
      evap_12_igcm(:,:,:)= latent_12_igcm(:,:,:)/
     :                 (2.5e6)
      netlong_12_igcm(:,:,:)= netlong_12_igcm(:,:,:)*(-1.0)
      stressx_12_igcm(:,:,:)= stressx_12_igcm(:,:,:)*(-1.0)
      stressy_12_igcm(:,:,:)= stressy_12_igcm(:,:,:)*(-1.0)
c     NCEP runoff is in units of kg/ms/6-hours
      runoff3_12_igcm(:,:,:)= runoff3_12_igcm(:,:,:)/21600.0


c     OK, hardest bit: runoff

c     This runoff based on the idea that the runoff is equal to 
c       P-E over land
      do m=1,nmonths
      do j=1,ilat1_atm
        do i=1,ilon1_atm
  
          if (lsm_1_igcm(i,j).eq.1) then
          runoff1_12_igcm(nint(index_2_igcm(i,j,1)),
     :                   nint(index_2_igcm(i,j,2)), m ) = 
     :    runoff1_12_igcm(nint(index_2_igcm(i,j,1)),
     :                   nint(index_2_igcm(i,j,2)), m ) +
     :            (precip_12_igcm(i,j,m)+evap_12_igcm(i,j,m))*
     :        weight_atm(i,j) /
     :        weight_atm(index_2_igcm(i,j,1),
     :                   index_2_igcm(i,j,2) )
          endif
        enddo
      enddo
      enddo

c     This runoff is directly from NCEP
      do m=1,nmonths
      do j=1,ilat1_atm
        do i=1,ilon1_atm
          if (lsm_1_igcm(i,j).eq.1) then
          runoff2_12_igcm(nint(index_2_igcm(i,j,1)),
     :                   nint(index_2_igcm(i,j,2)), m ) = 
     :    runoff2_12_igcm(nint(index_2_igcm(i,j,1)),
     :                   nint(index_2_igcm(i,j,2)), m ) +
     :            runoff3_12_igcm(i,j,m)*
     :        weight_atm(i,j) /
     :        weight_atm(index_2_igcm(i,j,1),
     :                   index_2_igcm(i,j,2) )
          endif
        enddo
      enddo
      enddo


      print*,'post-processing complete'

      print*,'calculating diagnostics'


c     Now calculate some diagnostics:

      do j=1,ilat1_atm
        do i=1,ilon1_atm
           if ( (alon1(i).gt.290).or.(alon1(i).lt.20) ) 
     :             at_mask(i,j)=1
           if (alat1(j).lt.-50)
     :             at_mask(i,j)=0
           if (alat1(j).gt.65)
     :             at_mask(i,j)=1
           if ( (alon1(i).gt.230).and.(alat1(j).gt.50) ) 
     :             at_mask(i,j)=1
           if ( (alon1(i).gt.265).and.(alat1(j).gt.10) ) 
     :             at_mask(i,j)=1
           if ( (alon1(i).lt.50).and.(alat1(j).gt.30) ) 
     :             at_mask(i,j)=1
            if ( (alon1(i).gt.260).and.(alat1(j).gt.20) ) 
     :             at_mask(i,j)=1
           at_mask(i,j)=at_mask(i,j)*(1-lsm_1_igcm(i,j))
        enddo
      enddo

      at1=0.0
      at2=0.0
      at3=0.0
      do m=1,nmonths
      do j=1,ilat1_atm
        do i=1,ilon1_atm
          if (alat1(j).gt.30) then
            at1=at1+weight_atm(i,j)*at_mask(i,j)*
     : (precip_12_igcm(i,j,m)+evap_12_igcm(i,j,m)+
     :  runoff1_12_igcm(i,j,m))
          else if (alat1(j).lt.-30) then
            at3=at3+weight_atm(i,j)*at_mask(i,j)*
     : (precip_12_igcm(i,j,m)+evap_12_igcm(i,j,m)+
     :  runoff1_12_igcm(i,j,m))
          else
            at2=at2+weight_atm(i,j)*at_mask(i,j)*
     : (precip_12_igcm(i,j,m)+evap_12_igcm(i,j,m)+
     :  runoff1_12_igcm(i,j,m))
          endif
        enddo
      enddo
      enddo

c     12 is for the number of months, 1e12 is for Sverdrups and kg =>m3
      at1=at1*4.0*pi*radea*radea/12e9
      at2=at2*4.0*pi*radea*radea/12e9
      at3=at3*4.0*pi*radea*radea/12e9

      print*,'diagnostics complete'


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


      nvar=12

      v=1
      varname(v,1)='prate'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='prate from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='kg/m2/s'

      v=v+1
      varname(v,1)='sensible'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='sensible heat from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='W/m2'

      v=v+1
      varname(v,1)='latent'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='latent heat from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='W/m2'

      v=v+1
      varname(v,1)='evap'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='evaporation (derived) from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='kg/m2/s'

      v=v+1
      varname(v,1)='netsolar'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='net surface solar radiation from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='W/m2'

      v=v+1
      varname(v,1)='netlong'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='net surface longwave radiation from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='W/m2'

      v=v+1
      varname(v,1)='stressx'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='x surface wind stress from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='N/m2'

      v=v+1
      varname(v,1)='stressy'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='y surface wind stress from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='N/m2'

      v=v+1
      varname(v,1)='runoff'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='runoff (derived) from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='kg/m2/s'

      v=v+1
      varname(v,1)='runoff2'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='runoff (direct) from NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='kg/m2/s'

      v=v+1
      varname(v,1)='at_mask'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='Mask for Atlantic'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=v+1
      varname(v,1)='lsm'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='land/sea of NCEP'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'


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
     :        precip_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        sensible_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        latent_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        evap_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        netsolar_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        netlong_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        stressx_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        stressy_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        runoff1_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        runoff2_12_igcm )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        at_mask )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        lsm_1_igcm )
      print*,'written file'


      call end_netcdf(1)


      print*,'files written OK'

      print*,'Diagnostics:'
      print*,'Atlantic 1',at1
      print*,'Atlantic 2',at2
      print*,'Atlantic 3',at3

      end


#include 'read_bconds.f'
