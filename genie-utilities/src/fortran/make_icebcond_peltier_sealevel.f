      program make_peltier_goldstein

c***********************************************************

c     This bit of code takes the output from Peltier ICE-4G, 
c       and interpolates it onto the GOLDSTEIN grid, making use of a 
c       sealevel curve.
c     
c
c     I compile with: 
c     ifort -o make_icebcond_peltier_sealevel_NOCVS make_icebcond_peltier_sealevel.f -lnc1 -lutil1 -ligcm3 -lnetcdf -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-igcm3/lib -L/home/ggdjl/genie/genie-lib/libutil1 -L/opt/local/intel_fc_81/lib -CB
c
      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
c
      include 'names_icebcond_peltier_sealevel_thompson_125k.inc'
c
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
      real frac_g(ilon_g,ilat_g)
      real frac_write(ilon_g,ilat_g,nexp)
      real lsm_g(ilon_g,ilat_g)
      logical work_g(ilon_g,ilat_g)

      real orog_p(ilon_p,ilat_p)
      real frac_p(ilon_p,ilat_p)
      real lsm_p(ilon_p,ilat_p)
      logical work_p(ilon_p,ilat_p)

      real,dimension(nexp) :: expname

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
      integer i,j,v,e,t

c     For goldstein grid:
      real pi
      real th0
      real th1
      real s0,s1,ds
      real sv(0:ilat_g)
      real s(0:ilat_g)

c     For options
      logical fix_mask
      parameter(fix_mask=.true.)
      real orog_limit
      parameter(orog_limit=20.0)

c     For sealevel
      integer nsea
      parameter(nsea=251)
      integer ntim
      parameter(ntim=126)
      real sealevel(nsea,2)
      real dummy(2)
      real p_sea(nexp)
      real t_s
      parameter(t_s=2.0)
      real this_sealevel
      integer this_index
      real this_fraction
      real orog_write_new(ilon_g,ilat_g,ntim)
      real frac_write_new(ilon_g,ilat_g,ntim)
      real,dimension(ntim) :: expname_new

      print*,'make_icebcond is starting....'

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



c     *******************************
c     OK, read in data....

c     Read in goldstein land-sea mask
      call open_file_nc(trim(in_lsm_filename),ncid)
      call get2d_data_nc(ncid,trim(in_lsm_varname),ilon_g,ilat_g,
     &     lsm_g,ifail)
      call close_file_nc(trim(in_lsm_filename),ncid)
c     process the goldstein land-sea mask to 1/0...
      do j=1,ilat_g
        do i=1,ilon_g
          if (lsm_g(i,j).eq.-99999.0) then
             lsm_g(i,j)=1.0
          else
             lsm_g(i,j)=0.0
          endif 
        enddo
      enddo

c     open the output files
      open(21,file=out_ascii_orogname)
      open(22,file=out_ascii_icefracname)

      do e=1,nexp
      write(my_name,FMT1) nexp-e

      call open_file_nc(trim(in_orog_filename)//
     &     trim(my_name)//'.nc',ncid)
      call get2d_data_nc(ncid,trim(in_orog_varname),ilon_p,ilat_p,
     &     orog_p,ifail)
      call close_file_nc(trim(in_orog_filename)//
     &     trim(my_name)//'.nc',ncid)

      call open_file_nc(trim(in_icefrac_filename)//
     &     trim(my_name)//'.nc',ncid)
      call get2d_data_nc(ncid,trim(in_icefrac_varname),ilon_p,ilat_p,
     &     frac_p,ifail)
      call close_file_nc(trim(in_icefrac_filename)//
     &     trim(my_name)//'.nc',ncid)

      print*,'all data read in'

c     Prepare

c     OK, we should really use the land-sea mask info, otherwise
c       we will take ocean points into account when getting 
c       orography and ice

      do j=1,ilat_p
        do i=1,ilon_p
          if (orog_p(i,j).le.0.0) then
             orog_p(i,j)=0.0
             work_p(i,j)=.false.
          else
             work_p(i,j)=.true.
          endif 
        enddo
      enddo

      do j=1,ilat_g
        do i=1,ilon_g
          if (lsm_g(i,j).eq.1) then
            work_g(i,j)=.true.
          else
            work_g(i,j)=.false.
          endif
        enddo
      enddo


c     Interpolate onto goldstein grid, using land-sea info
         call awi(ilon_p,ilon_p,elon_p,
     :             ilat_p,elat_p,orog_p,work_p,
     :             ilon_g,ilon_g,elon_g,
     :             ilat_g,elat_g,orog_g,work_g,
     :             ifail)
         if (ifail.ne.0) then
            print*,'Interp_to_ocn:  Error in awi,   itype 1 ',ifail
         end if


      do j=1,ilat_g
        do i=1,ilon_g
          if (lsm_g(i,j).eq.1) then
            work_g(i,j)=.true.
          else
            work_g(i,j)=.false.
          endif
        enddo
      enddo

         call awi(ilon_p,ilon_p,elon_p,
     :             ilat_p,elat_p,frac_p,work_p,
     :             ilon_g,ilon_g,elon_g,
     :             ilat_g,elat_g,frac_g,work_g,
     :             ifail)
         if (ifail.ne.0) then
            print*,'Interp_to_ocn:  Error in awi,   itype 1 ',ifail
         end if

      print*,'finished interpolating'

      if (fix_mask) then

c     Preserve the land-sea mask
      orog_g(:,:)=orog_g(:,:)*lsm_g(:,:)
      frac_g(:,:)=frac_g(:,:)*lsm_g(:,:)

c     Make areas which are = 50 to 1 if icemask and to orog_limit if orog
      do j=1,ilat_g
        do i=1,ilon_g
          if (frac_g(i,j).eq.50) then
             frac_g(i,j)=0.0
          endif
          if (orog_g(i,j).eq.50) then
             orog_g(i,j)=orog_limit
          endif
        enddo
      enddo

c     Make no land less than orog_limit height
      do j=1,ilat_g
        do i=1,ilon_g
          if (orog_g(i,j).lt.orog_limit) then
             orog_g(i,j)=orog_limit
          endif 
        enddo
      enddo

c     Just a sanity-check to ensure ice is between 0 and 1
      do j=1,ilat_g
        do i=1,ilon_g
          if (frac_g(i,j).gt.1.0) then
             frac_g(i,j)=1.0
          endif
          if (frac_g(i,j).lt.0.0) then
             frac_g(i,j)=0.0
          endif
        enddo
      enddo

      endif

c     The ice data should be 0 for sea, 1 for land, and 2 for ice-covered.  
c     Values between 1 and 2 represent fractional ice-coverage
c     Also, give orography a minimum value of 20 m if there is land.
c     This could probably be relaxed.
      do j=1,ilat_g
        do i=1,ilon_g
          if (lsm_g(i,j).gt.0.0) then

          if (frac_g(i,j).gt.0.0) then
             frac_g(i,j)=frac_g(i,j)+1.0
          else
             frac_g(i,j)=1.0
          endif 

          endif
        enddo
      enddo

      print*,'finished processing'

      frac_write(:,:,e)=frac_g(:,:)
      orog_write(:,:,e)=orog_g(:,:)

      print*,'writing ascii output'

      do i=1,ilon_g
      do j=1,ilat_g
      write(21,*) orog_g(i,j)
      enddo 
      enddo

      do i=1,ilon_g
      do j=1,ilat_g
      write(22,*) frac_g(i,j)
      enddo 
      enddo

      enddo

      close(21)
      close(22)

      print*,'Peltier ascii files written OK'

c     OK, now use the sea-level curve to extrapolate this data.....

c     Read in the sea-level data
      open(23,file='/home/ggdjl/genie/genie-utilities/data/input/'//
     c           'thompson_125k_test.dat')
      do j=1,nsea
        read(23,'(2e16.12e3)') sealevel(j,1),sealevel(j,2)
      enddo
      close(23)

c     and open the ascii data 
      open(24,file=out_ascii_orogname_new)
      open(25,file=out_ascii_icefracname_new)

c     This sea-level data is in the opposite direction to the peltier data...
c     First of all, extract the same period as peltier....
      print*,'number of exp    kyrBP     Peltier sealevel'
      do e=1,nexp
        p_sea(e)=sealevel(nexp*t_s+1-e*t_s,2)
        print*,e,nexp-e,p_sea(e)
      enddo

      print*

      print*,'number of exp    kyrBP  new index  New sealevel'
c     Now loop over all time...
      if ((nsea-1)/t_s+1.ne.ntim) then
        print*,'MAJOR PROBLEM!!'
        stop
      endif

      do t=1,ntim
        this_sealevel=sealevel(nsea-(2*(t-1)),2)
c     if sealevel is 0, use pre-industrial value
        if (this_sealevel.eq.0) then
          this_index=nexp-1
          this_fraction=1.0
c     if sealevel is LGM or higher, use LGM value
        elseif (this_sealevel.le.p_sea(1)) then
          this_index=1
          this_fraction=0.0
        else
c     now locate this sealevel within the peltier curve....
          this_index=-999
          this_fraction=0.0
          do e=1,nexp-1
            if (this_sealevel.ge.p_sea(e).and.
     c          this_sealevel.le.p_sea(e+1)) then
              this_index=e
            endif
          enddo
          if (this_index.ne.-999) then
            if (p_sea(this_index+1)-
     c        p_sea(this_index).ne.0) then
                this_fraction=1-(this_sealevel-p_sea(this_index))/
     c            (p_sea(this_index+1)-p_sea(this_index))
            else
              stop
            endif
          else
            stop          
          endif
        endif

        print*,t,ntim-t,this_index,this_sealevel,'=',p_sea(this_index),
     c    '*',this_fraction,'+',p_sea(this_index+1),'*',1-this_fraction

        orog_write_new(:,:,t)=orog_write(:,:,this_index)*this_fraction+
     c    orog_write(:,:,this_index+1)*(1-this_fraction)
        frac_write_new(:,:,t)=frac_write(:,:,this_index)*this_fraction+
     c    frac_write(:,:,this_index+1)*(1-this_fraction)

      do i=1,ilon_g
      do j=1,ilat_g
      write(24,*) orog_write_new(i,j,t)
      enddo 
      enddo

      do i=1,ilon_g
      do j=1,ilat_g
      write(25,*) frac_write_new(i,j,t)
      enddo 
      enddo

      enddo

      close(24)
      close(25)

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
     :        frac_write )
      print*,'written file'

      call end_netcdf(1)

      do i=1,ntim
      expname_new(i)=ntim-i
      enddo

      ndims(3)=ntim
      call ininc(out_data_filename_new,
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))

      call writedim(nc(1),iddim(1,1),alon_g)
      call writedim(nc(1),iddim(2,1),alat_g)
      call writedim(nc(1),iddim(3,1),expname_new)

      v=1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        orog_write_new )
      print*,'written file'
      v=v+1
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        frac_write_new )
      print*,'written file'

      call end_netcdf(1)


      print*,'netcdf files written OK'




      end


