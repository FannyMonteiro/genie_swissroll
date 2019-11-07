      program make_fixed_igcm

c     Compile with: 
c     ifc -o make_fixed_igcm make_fixed_igcm.f -lnc1 -lutil1 -lnetcdf_ifc -ligcm -L/home/ggdjl/genie/genie-lib 

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/netdata.cmn'
      include '/home/ggdjl/genie/genie-igcm3/param1.cmn'
      include '/home/ggdjl/genie/genie-igcm3/param2.cmn'

c***********************************************************
c     To change the spectral orography:
c     (1) Make a utf file of the netcdf orography (use xconv)
c     (2) Edit specutf_lgm.in or specutf_cnt.in in and out names
c     (4) Run specutf.f!
c***********************************************************

      include 'names_igcm.inc' 

c***********************************************************


      integer in_nlon,in_nlat
      parameter(in_nlon=64,in_nlat=32)

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
      real precip_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real sensible_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real latent_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real netsolar_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real netlong_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real stressx_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real stressy_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real runoff_12_igcm(ilon1_atm,ilat1_atm,nmonths)
      real lsm_1_igcm(ilon1_atm,ilat1_atm)

      real datain_12(ilon1_atm,ilat1_atm,1,nmonths)
      real datain_1(ilon1_atm,ilat1_atm)

      integer work_igcm(ilon1_atm,ilat1_atm)


   

c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a

c     LOOPING:
      integer i,j,m,v


c     FOR GRID:
      real ax

      print*,'make_fixed_igcm is starting....'

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

      call open_file_nc(trim(in_data_maskname),ncid)
      call get2d_data_nc(ncid,trim(in_mask_varname),in_nlon,in_nlat,
     &          datain_1,ifail)
      if (ifail.ne.0) stop 'Problem reading input mask file'
      lsm_1_igcm(:,:)=datain_1(:,:)
      call close_file_nc(trim(in_data_maskname),ncid)

      call open_file_nc(trim(in_data_filename),ncid)

      call get4d_data_nc(ncid,trim(in_sst_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      sst_12_igcm(:,:,:)=datain_12(:,:,1,:)

      call get4d_data_nc(ncid,trim(in_precip1_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      precip_12_igcm(:,:,:)=datain_12(:,:,1,:)

c     I have made a change here....the original was in_precip1_varname
      call get4d_data_nc(ncid,trim(in_precip2_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      precip_12_igcm(:,:,:)=precip_12_igcm(:,:,:)+datain_12(:,:,1,:)

      call get4d_data_nc(ncid,trim(in_sensible_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      sensible_12_igcm(:,:,:)=datain_12(:,:,1,:)

      call get4d_data_nc(ncid,trim(in_latent_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      latent_12_igcm(:,:,:)=datain_12(:,:,1,:)

      call get4d_data_nc(ncid,trim(in_netsolar_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      netsolar_12_igcm(:,:,:)=datain_12(:,:,1,:)

      call get4d_data_nc(ncid,trim(in_netlong_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      netlong_12_igcm(:,:,:)=datain_12(:,:,1,:)

      call get4d_data_nc(ncid,trim(in_stressx_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      stressx_12_igcm(:,:,:)=datain_12(:,:,1,:)

      call get4d_data_nc(ncid,trim(in_stressy_varname),in_nlon,in_nlat,
     &     1,nmonths,datain_12,ifail)
      if (ifail.ne.0) stop 'Problem reading input file'
      stressy_12_igcm(:,:,:)=datain_12(:,:,1,:)

c     Haven't done runoff yet.....
      runoff_12_igcm(:,:,:)=0.0

      call close_file_nc(trim(in_data_filename),ncid)



       print*,'all data read in'


c     process the data...............



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


      nvar=10

      v=1
      varname(v,1)='airs'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='airs from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='degrees C'

      v=2
      varname(v,1)='prate'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='prate from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=3
      varname(v,1)='sensible'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='sensible heat from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=4
      varname(v,1)='latent'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='latent heat from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=5
      varname(v,1)='netsolar'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='net surface solar radiation from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=6
      varname(v,1)='netlong'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='net surface longwave radiation from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=7
      varname(v,1)='stressx'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='x surface wind stress from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=8
      varname(v,1)='stressy'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='y surface wind stress from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=9
      varname(v,1)='runoff'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='runoff from IGCM'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      v=10
      varname(v,1)='lsm'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='land/sea of IGCM'
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
     :        sst_12_igcm )
      print*,'written file'
  
      v=2
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        precip_12_igcm )
      print*,'written file'

      v=3
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        sensible_12_igcm )
      print*,'written file'

      v=4
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        latent_12_igcm )
      print*,'written file'

      v=5
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        netsolar_12_igcm )
      print*,'written file'

      v=6
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        netlong_12_igcm )
      print*,'written file'

      v=7
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        stressx_12_igcm )
      print*,'written file'

      v=8
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        stressy_12_igcm )
      print*,'written file'

      v=9
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        runoff_12_igcm )
      print*,'written file'

      v=10
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        lsm_1_igcm )
      print*,'written file'

      call end_netcdf(1)

      print*,'files written OK'

      end

