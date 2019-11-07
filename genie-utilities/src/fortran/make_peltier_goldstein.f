      program make_peltier_goldstein

c***********************************************************

c     This bit of code takes the output from Peltier, and interpolates
c       it onto the goldstein grid.
c     For the 30kyr simulation, you will need to first run
c       make_hadam_peltier
c     For the new 30kyr simulation, you will need to first run
c       make_hadam_peltier and make_peltier_interp
c     

c     I compile with: 
c     ifort -o make_peltier_goldstein_NOCVS make_peltier_goldstein.f -lnc1 -lutil1 -ligcm3 -lnetcdf -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-igcm3/lib -L/home/ggdjl/genie/genie-lib/libutil1 -L/opt/local/intel_fc_81/lib -CB

c     The lines between the asterixes will need to be edited.

c     Please edit the main body of the code also, to tidy it up!

c     Eventually, this piece of code could run with makefile.arc!

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'

c      include 'names_peltier_goldstein_pre.inc' 
c      include 'names_peltier_goldstein_hol.inc' 
c      include 'names_peltier_goldstein_lgm.inc' 
c      include 'names_peltier_goldstein_st3.inc'
c      include 'names_peltier_goldstein_27e.inc'
c      include 'names_peltier_goldstein_24e.inc'
c      include 'names_peltier_goldstein_18e.inc'
c      include 'names_peltier_goldstein_15e.inc'
c      include 'names_peltier_goldstein_12e.inc'
c      include 'names_peltier_goldstein_09e.inc'
c      include 'names_peltier_goldstein_03e.inc'
c      include 'names_peltier_goldstein_21k.inc'  
c      include 'names_peltier_goldstein_30k.inc'  

      include 'names_peltier_goldstein_30k_new.inc'  

c***********************************************************

c     GOLDSTEIN GRID:
      integer ilon_g,ilat_g
      parameter (ilon_g=36,ilat_g=36)
      real alon_g(ilon_g),alat_g(ilat_g),
     :     elon_g(ilon_g+1),elat_g(ilat_g+1)

c     PELTIER GRID:
      integer ilon_p,ilat_p
      parameter (ilon_p=360,ilat_p=180)
      real alon_p(ilon_p),alat_p(ilat_p),
     :     elon_p(ilon_p+1),elat_p(ilat_p+1)

c     DATA:
      real orog_g(ilon_g,ilat_g)
      real orog_write(ilon_g,ilat_g,nexp)
      real icefrac_g(ilon_g,ilat_g)
      real icefrac_write(ilon_g,ilat_g,nexp)
      real lsm_g(ilon_g,ilat_g)
      logical work_g(ilon_g,ilat_g)

      real orog_p(ilon_p,ilat_p)
      real icefrac_p(ilon_p,ilat_p)
      real lsm_p(ilon_p,ilat_p)
      logical work_p(ilon_p,ilat_p)


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
      integer i,j,v,e

c     For goldstein grid:
      real pi
      real th0
      real th1
      real s0,s1,ds
      real sv(0:ilat_g)
      real s(0:ilat_g)

      print*,'make_hadam_goldstein is starting....'

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


c     OK, read in data....

      open(21,file=out_ascii_orogname)
      open(22,file=out_ascii_icefracname)



      do e=1,nexp

      call open_file_nc(trim(in_lsm_filename),ncid)
      call get2d_data_nc(ncid,trim(in_lsm_varname),ilon_g,ilat_g,
     &     lsm_g,ifail)
      call close_file_nc(trim(in_lsm_filename),ncid)

      write(my_name,FMT1) nexp-e

      call open_file_nc(trim(in_orog_filename)//my_name//'.nc',ncid)
      call get2d_data_nc(ncid,trim(in_orog_varname),ilon_p,ilat_p,
     &     orog_p,ifail)
      call close_file_nc(trim(in_orog_filename)//my_name//'.nc',ncid)

      call open_file_nc(trim(in_icefrac_filename)//my_name//'.nc',ncid)
      call get2d_data_nc(ncid,trim(in_icefrac_varname),ilon_p,ilat_p,
     &     icefrac_p,ifail)
      call close_file_nc(trim(in_icefrac_filename)//my_name//'.nc',ncid)

      print*,'all data read in'

c     Prepare

      do j=1,ilat_p
        do i=1,ilon_p
          if (orog_p(i,j).lt.0.0) then
             orog_p(i,j)=0.0
          endif 
        enddo
      enddo

      work_p(:,:)=.true.
      work_g(:,:)=.true.

c     Interpolate onto goldstein grid

         call awi(ilon_p,ilon_p,elon_p,
     :             ilat_p,elat_p,orog_p,work_p,
     :             ilon_g,ilon_g,elon_g,
     :             ilat_g,elat_g,orog_g,work_g,
     :             ifail)
         if (ifail.ne.0) then
            print*,'Interp_to_ocn:  Error in awi,   itype 1 ',ifail
            stop 1
         end if

         call awi(ilon_p,ilon_p,elon_p,
     :             ilat_p,elat_p,icefrac_p,work_p,
     :             ilon_g,ilon_g,elon_g,
     :             ilat_g,elat_g,icefrac_g,work_g,
     :             ifail)
         if (ifail.ne.0) then
            print*,'Interp_to_ocn:  Error in awi,   itype 1 ',ifail
            stop 1
         end if


      print*,'finished interpolating'

c     process...

      do j=1,ilat_g
        do i=1,ilon_g
          if (lsm_g(i,j).eq.-99999.0) then
             lsm_g(i,j)=1.0
          else
             lsm_g(i,j)=0.0
          endif 
        enddo
      enddo

      orog_g(:,:)=orog_g(:,:)*lsm_g(:,:)
      icefrac_g(:,:)=icefrac_g(:,:)*lsm_g(:,:)


c     ****************************************
c     THIS IS THE OLD, INTEGER CODE:

      do j=1,ilat_g
        do i=1,ilon_g
          if (icefrac_g(i,j).gt.0.0) then
          if (icefrac_g(i,j).lt.0.3) then
             icefrac_g(i,j)=0.0
          else
             icefrac_g(i,j)=1.0
          endif 
          endif
        enddo
      enddo

      do j=1,ilat_g
        do i=1,ilon_g
          if (lsm_g(i,j).gt.0.0) then
          if (icefrac_g(i,j).eq.1.0) then
             icefrac_g(i,j)=2.0
          else
             icefrac_g(i,j)=1.0
          endif 
          endif
        enddo
      enddo

c     ****************************************

c     ****************************************
c     THIS IS THE NEW, REAL CODE:

c      do j=1,ilat_g
c        do i=1,ilon_g
c          if (icefrac_g(i,j).gt.1.0) then
c             icefrac_g(i,j)=1.0
c          endif
c          if (icefrac_g(i,j).lt.0.0) then
c             icefrac_g(i,j)=0.0
c          endif
c        enddo
c      enddo

c      do j=1,ilat_g
c        do i=1,ilon_g
c          if (lsm_g(i,j).gt.0.0) then
c          if (icefrac_g(i,j).gt.0.0) then
c             icefrac_g(i,j)=icefrac_g(i,j)+1.0
c          else
c             icefrac_g(i,j)=1.0
c          endif 
c          endif
c        enddo
c      enddo

c     ****************************************

      do j=1,ilat_g
        do i=1,ilon_g
          if (lsm_g(i,j).gt.0.0) then
          if (orog_g(i,j).lt.20) then
             orog_g(i,j)=20.0
          endif 
          endif
        enddo
      enddo


      print*,'finished processing'


      icefrac_write(:,:,e)=icefrac_g(:,:)
      orog_write(:,:,e)=orog_g(:,:)


      print*,'writing ascii output'

      do i=1,ilon_g
      do j=1,ilat_g
      write(21,*) orog_g(i,j)
      enddo 
      enddo

      do i=1,ilon_g
      do j=1,ilat_g
      write(22,*) icefrac_g(i,j)
      enddo 
      enddo


      enddo

      close(21)
      close(22)

      print*,'ascii files written OK'

      print*,'writing netcdf output'

      do i=1,nexp
      expname(i)=nexp-i
      enddo

      ndim=3
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

      dimname(3,1)='time'
      ndims(3)=nexp
      natts(3)=2
      attdimname(1,1,3)='long_name'
      attdimname(2,1,3)='time'
      attdimname(1,2,3)='units'
      attdimname(2,2,3)='kyrBP'

      nvar=2

      v=1
      varname(v,1)='orog'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='orography from peltier'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='m'

      v=v+1
      varname(v,1)='icefrac'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('time',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='fraction of land-ice from peltier'
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
      call writedim(nc(1),iddim(3,1),expname)


      v=1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        orog_write )
      print*,'written file'

      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        icefrac_write )
      print*,'written file'


      call end_netcdf(1)

      print*,'netcdf files written OK'




      end


