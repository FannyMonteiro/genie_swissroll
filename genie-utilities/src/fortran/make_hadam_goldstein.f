      program make_hadam_goldstein

c***********************************************************

c     This bit of code takes the output from a Hadley-centr run, and 
c       interpolates
c       it onto the goldstein grid.

c     I compile with: 
c     ifort -o make_hadam_goldstein_NOCVS make_hadam_goldstein.f -lnc1 -lutil1 -ligcm3 -lnetcdf -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-igcm3/lib -L/home/ggdjl/genie/genie-lib/libutil1 -L/opt/local/lib

c     The lines between the asterixes will need to be edited.

c     Please edit the main body of the code also, to tidy it up!

c     Eventually, this piece of code could run with makefile.arc!

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'

c      include 'names_hadam_goldstein_xakxd.inc' 
      include 'names_hadam_goldstein_xatxd.inc' 

c***********************************************************

c     GOLDSTEIN GRID:
      integer ilon_g,ilat_g
      parameter (ilon_g=36,ilat_g=36)
      real alon_g(ilon_g),alat_g(ilat_g),
     :     elon_g(ilon_g+1),elat_g(ilat_g+1)

c     HADLEY GRID:
      integer ilon_h,ilat_h
      parameter (ilon_h=96,ilat_h=73)
      real alon_h(ilon_h),alat_h(ilat_h),
     :     elon_h(ilon_h+1),elat_h(ilat_h+1)

c     DATA:
      real temp_g(ilon_g,ilat_g)
      real lsm_g(ilon_g,ilat_g)
      logical work_g(ilon_g,ilat_g)

      real temp_h(ilon_h,ilat_h)
      real lsm_h(ilon_h,ilat_h)
      logical work_h(ilon_h,ilat_h)


c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a

c     LOOPING:
      integer i,j,v

c     For goldstein grid:
      real pi
      real th0
      real th1
      real s0,s1,ds
      real sv(0:ilat_g)
      real s(0:ilat_g)




      print*,'make_hadam_goldstein is starting....'

c     Set up the hadley grid

      do i=1,ilon_h
        alon_h(i)=3.75*(i-1)
        elon_h(i)=-1.875+3.75*(i-1)
      enddo
      elon_h(ilon_h+1)=-1.875+3.75*(ilon_h)

      do j=1,ilat_h
        alat_h(j)=90-2.5*(j-1)
        elat_h(j)=91.25-2.5*(j-1)
      enddo
      elat_h(1)=90
      elat_h(ilat_h+1)=-90


c     Set up the goldstein grid (copied form initialise_goldstein)

      pi = 4*atan(1.0)
      th0 = - pi/2    
      th1 = pi/2 
      s0 = sin(th0)    
      s1 = sin(th1)     
      ds = (s1-s0)/ilat_g
      do j=0,ilat_g
         sv(j) = s0 + j*ds
         s(j) = sv(j) - 0.5*ds
         if(s(j).lt.-1.0)s(j) = -2.0 - s(j)
      enddo  


      do i=1,ilon_g
         alon_g(i)=360.0*(i-0.5)/real(ilon_g)-260.0
      end do
      do i=1,ilon_g+1
         elon_g(i)=360.0*(i-1.0)/real(ilon_g)-260.0
      end do

      do j=1,ilat_g
         alat_g(j)=asin(s(j))*180.0/pi
      end do
      do j=1,ilat_g+1
         elat_g(j)=asin(sv(j-1))*180.0/pi
      end do

      print*,'grids set up OK in main'



c     process the data...............


c     OK, read in hadley-centre temp and lsm data....

      call open_file_nc(trim(in_temp_filename),ncid)
      call get4d_data_nc(ncid,trim(in_temp_varname),ilon_h,ilat_h,1,1,
     &     temp_h,ifail)
      call close_file_nc(trim(in_temp_filename),ncid)

      call open_file_nc(trim(in_lsm_filename),ncid)
      call get4d_data_nc(ncid,trim(in_lsm_varname),ilon_h,ilat_h,1,1,
     &     lsm_h,ifail)
      call close_file_nc(trim(in_lsm_filename),ncid)

      print*,'all data read in'

c     Process

      temp_h(:,:)=temp_h(:,:)-273.15
      work_h(:,:)=.true.
      work_g(:,:)=.true.


      print*,'finished processing'


c     Interpolate onto goldstein grid

         call awi(ilon_h,ilon_h,elon_h,
     :             ilat_h,elat_h,temp_h,work_h,
     :             ilon_g,ilon_g,elon_g,
     :             ilat_g,elat_g,temp_g,work_g,
     :             ifail)
         if (ifail.ne.0) then
            print*,'Interp_to_ocn:  Error in awi,   itype 1 ',ifail
            stop 1
         end if

      print*,'finished interpolating'



      print*,'writing netcdf output'

      ndim=2
      dimname(1,1)='longitude'
      ndims(1)=ilon_g
      natts(1)=2
      attdimname(1,1,1)='long_name'
      attdimname(2,1,1)='longitude'
      attdimname(1,2,1)='units'
      attdimname(2,2,1)='degrees east'

      dimname(2,1)='latitude'
      ndims(2)=ilat_g
      natts(2)=2
      attdimname(1,1,2)='long_name'
      attdimname(2,1,2)='latitude'
      attdimname(1,2,2)='units'
      attdimname(2,2,2)='degrees north'


      nvar=2

      v=1
      varname(v,1)='temp_1'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='surface temp from Hadley model'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='oC'

      v=v+1
      varname(v,1)='lsm'
      vdims(v)=2
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='Goldstein land-sea mask'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='0-1'

      call ininc(out_data_filename,
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))

      call writedim(nc(1),iddim(1,1),alon_g)
      call writedim(nc(1),iddim(2,1),alat_g)

      v=1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        temp_g )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        lsm_g )
      print*,'written file'


      call end_netcdf(1)


      print*,'netcdf files written OK'

      print*,'writing ascii output'
      open(1,file=out_ascii_filename)
      do j=1,ilat_g
      do i=1,ilon_g
      write(1,*) temp_g(i,j)
      enddo 
      enddo
      close(1)

      print*,'ascii files written OK'

      end


