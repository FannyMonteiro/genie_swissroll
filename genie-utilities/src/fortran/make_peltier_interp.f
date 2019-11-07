      program make_peltier_interp

c***********************************************************

c     This bit of code takes the output from a Hadley-centr run, and 
c       interpolates
c       it onto the peltier grid.
c     Performs an interpolation from 30kyBP to 21kyrBP.

c     I compile with: 
c     ifort -o make_peltier_interp_NOCVS make_peltier_interp.f -lnc1 -lutil1 -ligcm3 -lnetcdf -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-igcm3/lib -L/home/ggdjl/genie/genie-lib/libutil1 -L/opt/local/intel_fc_81/lib -CB

c     The lines between the asterixes will need to be edited.

c     Please edit the main body of the code also, to tidy it up!

c     Eventually, this piece of code could run with makefile.arc!

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'

      include 'names_hadam_peltier_new2.inc'

c***********************************************************

c     PELTIER GRID:
      integer ilon_p,ilat_p
      parameter (ilon_p=360,ilat_p=180)
      real alon_p(ilon_p),alat_p(ilat_p),
     :     elon_p(ilon_p+1),elat_p(ilat_p+1)

      real orog_p(ilon_p,ilat_p)
      real icefrac_p(ilon_p,ilat_p)
      logical work_p(ilon_p,ilat_p)

      real orogwrite_p(ilon_p,ilat_p,nexp)
      real icefracwrite_p(ilon_p,ilat_p,nexp)

      real oroglgm_p(ilon_p,ilat_p)
      real icefraclgm_p(ilon_p,ilat_p)

      real orog_p_array(ilon_p,ilat_p,nexpin)
      real icefrac_p_array(ilon_p,ilat_p,nexpin)

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
      integer i,j,v,n,t,vi,vf,v1,v2

c     for interpolation
      real frac30

      real d18o(nd18o)

      print*,'make_hadam_peltier is starting....'

c     Set up the peltier grid

      do i=1,ilon_p
        alon_p(i)=-179.5+1.0*(i-1)
        elon_p(i)=-180.0+1.0*(i-1)
      enddo
      elon_p(ilon_p+1)=-180.0+1.0*(ilon_p)

      do j=1,ilat_p
        alat_p(j)=89.5-1.0*(j-1)
        elat_p(j)=90.0-1.0*(j-1)
      enddo
      elat_p(ilat_p+1)=-90.0

      print*,'grids set up OK in main'

c     process the data...............


c     OK, read in peltier lgm orog and icefrac data....

      do n=1,nexpin

      write(my_name_2,FMT2) n-1

      call open_file_nc(trim(in_icefraclgm_filename)//
     &     my_name_2//'.nc',ncid)
      call get2d_data_nc(ncid,trim(in_icefraclgm_varname),
     &     ilon_p,ilat_p,
     &     icefraclgm_p,ifail)
      call close_file_nc(trim(in_icefraclgm_filename)//
     &     my_name_2//'.nc',ncid)

      call open_file_nc(trim(in_oroglgm_filename)//
     &     my_name_2//'.nc',ncid)
      call get2d_data_nc(ncid,trim(in_oroglgm_varname),
     &     ilon_p,ilat_p,
     &     oroglgm_p,ifail)
      call close_file_nc(trim(in_oroglgm_filename)//
     &     my_name_2//'.nc',ncid)

         do j=1,ilat_p
          do i=1,ilon_p
             if (oroglgm_p(i,j).lt.0) then
               oroglgm_p(i,j)=0.0
             endif
          enddo
        enddo     

      orog_p_array(:,:,n)=oroglgm_p(:,:)
      icefrac_p_array(:,:,n)=icefraclgm_p(:,:)

      enddo

      print*,'all data read in'

c     Process

      print*,'finished processing'

      icefrac_p(:,:)=nint(icefrac_p(:,:))

      open(1,file=in_d18o_filename,form='formatted')
      read(1,*) d18o
      close(1)

      do n=1,nexp


c     the d18o datafile starts at 30kyrbp, whereas
c       n goes from 0 to 211
  
        t=(nexp-1)/(nexpin-1)

        vi=int(real(n-1)/real(t))+1
        vi=max(vi,1)
        vf=vi+1
        vf=min(vf,nexpin)
        v1=(vi-1)*t+1
        v2=vi*t+1


        frac30=(d18o(nd18o-n+1)-d18o(nd18o-v1+1))/
     &         (d18o(nd18o-v2+1)-d18o(nd18o-v1+1))

        write(6,'(i5,x,i3,x,i2,x,i2,x,i3,x,i3,x,f9.6,x,f9.6)') 
     &        (n-1)*100,n,vi,vf,v1,v2,frac30,d18o(nd18o-n+1)



        do j=1,ilat_p
          do i=1,ilon_p
             orogwrite_p(i,j,n)=
     &           (1.0-frac30)*orog_p_array(i,j,vi) + 
     &               (frac30)*orog_p_array(i,j,vf)
             icefracwrite_p(i,j,n)=
     &           (1.0-frac30)*icefrac_p_array(i,j,vi) + 
     &               (frac30)*icefrac_p_array(i,j,vf)
          enddo
        enddo

      enddo




      ndim=2
      dimname(1,1)='longitude'
      ndims(1)=ilon_p
      natts(1)=2
      attdimname(1,1,1)='long_name'
      attdimname(2,1,1)='longitude'
      attdimname(1,2,1)='units'
      attdimname(2,2,1)='degrees east'

      dimname(2,1)='latitude'
      ndims(2)=ilat_p
      natts(2)=2
      attdimname(1,1,2)='long_name'
      attdimname(2,1,2)='latitude'
      attdimname(1,2,2)='units'
      attdimname(2,2,2)='degrees north'

      do n=1,nexp

      write(my_name,FMT1) ni+n-1

      nvar=1
      v=1
      varname(v,1)='icemask'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='icefrac from Hadley model'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'
      attvarname(2,2,v)='m'
      call ininc(trim(out_icefractime_filename)//my_name//'.nc',
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))

      call writedim(nc(1),iddim(1,1),alon_p)
      call writedim(nc(1),iddim(2,1),alat_p)
      v=1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        icefracwrite_p(:,:,n) )
      print*,'written file'
      call end_netcdf(1)

      nvar=1
      v=1
      varname(v,1)='orog'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='orog from Hadley model'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='m'
      call ininc(trim(out_orogtime_filename)//my_name//'.nc' ,
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))
      call writedim(nc(1),iddim(1,1),alon_p)
      call writedim(nc(1),iddim(2,1),alat_p)
      v=1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        orogwrite_p(:,:,n) )
      print*,'written file'
      call end_netcdf(1)

      enddo







      print*,'netcdf files written OK'

      end


