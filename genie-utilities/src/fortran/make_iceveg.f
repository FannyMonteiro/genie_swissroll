      program make_iceveg

c     This bit of code takes in the original IGCM vegetation field
c      (genie-utilities/data/input/vegetation_t21.nc), and generates
c      six output files in genie-utilities/data/output.  Also
c      takes the original land-sea mask: landmask_std_t21.nc

c      1) vegetation_new_t21.nc is similar but includes antarctica.

c      These others are the ones which should be used with Genie: 
c      2) vegetation_soil_t21.nc replaces ice-sheets with bare-soil.
c      3) vegetation_tundra_t21.nc replaces ice-sheets with tundra.

c      When fixedicesheet is used, and the standard land-ice mask is used,
c      these vegetation types will be effectively under the ice and will
c      have no effect.  If a dynamic ice-sheet is used, when the ice melts,
c      the vegetation type will be exposed.
c      
c      Also produces ice-fraction fields, for use in the fixedicesheet:
c      4) icefrac_t21.nc is the control
c      5) icefrac_nogreen_t21.nc has no greenland ice-sheet
c      6) icefrac_noice_t21.nc has no ice-sheets

c     I compile with: 
c     ifort make_iceveg.f -o make_iceveg -lnc1 -lutil1 -ligcm3 -lnetcdf -L/home/ggejs/genie/genie-lib/libnc1 -L/home/ggejs/genie/genie-lib/libutil1 -L/home/ggejs/genie/genie-igcm3/lib -L/opt/local/intel_fc_81/lib -fpp
c     You will have to edit the path-names in the compilation, and in the
c       include statements below.

      implicit none
      include '/home/ggejs/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggejs/genie/genie-igcm3/src/fortran/netdata.cmn'
#include '/home/ggejs/genie/genie-igcm3/src/fortran/param1.cmn'
      include '/home/ggejs/genie/genie-igcm3/src/fortran/param2.cmn'

c     IGCM GRID:
      integer ilon1_atm,ilat1_atm
      parameter (ilon1_atm=64,ilat1_atm=32)
      real alon1(ilon1_atm),alat1(ilat1_atm)


c     DATA:
      real veg_old(ilon1_atm,ilat1_atm)
      real veg_new(ilon1_atm,ilat1_atm)
      real veg_new1(ilon1_atm,ilat1_atm)
      real veg_new2(ilon1_atm,ilat1_atm)
      real ice_new(ilon1_atm,ilat1_atm)
      real ice_new1(ilon1_atm,ilat1_atm)
      real ice_new2(ilon1_atm,ilat1_atm)

      real land(ilon1_atm,ilat1_atm)

c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim

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
      integer i,j

c     FOR GRID:
      real ax
      real nx


c     COPIED FROM INITIALISE_ATMOS.F
      ax=360.0/real(mg)
      do i=1,ilon1_atm
         alon1(i)=(i-1.0)*ax
      end do
      call gwtcnr(alat1,jg)

 
c      read in the orginal vegetation
      call open_file_nc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/input/'//
     &                        'vegetation_t21.nc',ncid)
      call get2d_data_nc(ncid,'svege',ilon1_atm,ilat1_atm,
     &                     veg_old,ifail)
      call close_file_nc('/home/ggdjl/genie_utilities/'//
     &                        'data/input/'//
     &                        'vegetation_t21.nc',ncid)
      print*,'orginal veg read in.....'


c      read in the land-sea mask
      call open_file_nc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/input/'//
     &                        'landmask_std_t21.nc',ncid)
      call get2d_data_nc(ncid,'lsm',ilon1_atm,ilat1_atm,
     &                     land,ifail)
      call close_file_nc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/input/'//
     &                        'landmask_std_t21.nc',ncid)
      print*,'orginal land read in.....'

      do j=1,ilat1_atm
        do i=1,ilon1_atm
           ice_new(i,j)=0.0
           ice_new1(i,j)=0.0
           ice_new2(i,j)=0.0
           veg_new(i,j)=veg_old(i,j)

c          SORT OUT ANTARCTICA.
c          In the original file, Antarctica was missing...
c          Put it in with code (2), which has a snow-free albedo of
c          0.75 and a snow-covered albedo of 0.8
           if ((land(i,j).eq.1).and.(alat1(j).lt.-60)) then
              veg_new(i,j)=2.
           endif
c          SORT OUT GREENLAND (and part of Alaska)
c          In the original they are the same code as inland lakes (3),
c          which is for e.g. the Caspian...
c          we want them to be the same as Antarctica (2) 
            if ((veg_new(i,j).eq.3).and.(alat1(j).gt.60)) then
              veg_new(i,j)=2.
           endif  

           veg_new1(i,j)=veg_new(i,j)
           veg_new2(i,j)=veg_new(i,j)

c          Replace ice-sheet with bare soil (24) or tundra (23).
c          Set the ice-sheet fraction
           if (veg_new(i,j).eq.2) then
             veg_new1(i,j)=24
             veg_new2(i,j)=23
             ice_new(i,j)=1
             ice_new1(i,j)=1
             if ( (alat1(j).gt.40).and.(alon1(i).gt.250) ) then
               ice_new1(i,j)=0
             endif
           endif

           veg_new1(i,j)=veg_new1(i,j)*land(i,j)
           veg_new2(i,j)=veg_new2(i,j)*land(i,j)
           ice_new(i,j)=ice_new(i,j)*land(i,j)
           ice_new1(i,j)=ice_new1(i,j)*land(i,j)
           ice_new2(i,j)=ice_new2(i,j)*land(i,j)

        enddo
      enddo

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

      nvar=1
      varname(1,1)='svege'
      vdims(1)=2
      vadims(1,1)=loc_dim('longitude',dimname,nall)
      vadims(2,1)=loc_dim('latitude',dimname,nall)
      nattsvar(1)=2
      attvarname(1,1,1)='long_name'
      attvarname(2,1,1)='vegetation (including Antarctica)'
      attvarname(1,2,1)='units'
      attvarname(2,2,1)='no units'
      call ininc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/output/'//
     &                    'vegetation_new_t21.nc',
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))
      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)
      call writevar(nc(1),
     :        idvar(loc_dim('svege',varname(1,1),nall),1),
     :        veg_new )
      call end_netcdf(1)

      nvar=1
      varname(1,1)='svege'
      vdims(1)=2
      vadims(1,1)=loc_dim('longitude',dimname,nall)
      vadims(2,1)=loc_dim('latitude',dimname,nall)
      nattsvar(1)=2
      attvarname(1,1,1)='long_name'
      attvarname(2,1,1)='vegetation (including Antarctica), '//
     &             'bare soil in place of icesheet'
      attvarname(1,2,1)='units'
      attvarname(2,2,1)='no units'
      call ininc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/output/'//
     &                    'vegetation_soil_t21.nc',
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))
      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)
      call writevar(nc(1),
     :        idvar(loc_dim('svege',varname(1,1),nall),1),
     :        veg_new1 )
      call end_netcdf(1)

      nvar=1
      varname(1,1)='svege'
      vdims(1)=2
      vadims(1,1)=loc_dim('longitude',dimname,nall)
      vadims(2,1)=loc_dim('latitude',dimname,nall)
      nattsvar(1)=2
      attvarname(1,1,1)='long_name'
      attvarname(2,1,1)='vegetation (including Antarctica), '//
     &             'tundra in place of icesheet'
      attvarname(1,2,1)='units'
      attvarname(2,2,1)='no units'
      call ininc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/output/'//
     &                    'vegetation_tundra_t21.nc',
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))
      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)
      call writevar(nc(1),
     :        idvar(loc_dim('svege',varname(1,1),nall),1),
     :        veg_new2 )
      call end_netcdf(1)

      nvar=1
      varname(1,1)='icefrac'
      vdims(1)=2
      vadims(1,1)=loc_dim('longitude',dimname,nall)
      vadims(2,1)=loc_dim('latitude',dimname,nall)
      nattsvar(1)=2
      attvarname(1,1,1)='long_name'
      attvarname(2,1,1)='ice fraction'
      attvarname(1,2,1)='units'
      attvarname(2,2,1)='0-1'
      call ininc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/output/'//
     &                    'icefrac_t21.nc',
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))
      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)
      call writevar(nc(1),
     :        idvar(loc_dim('icefrac',varname(1,1),nall),1),
     :        ice_new )
      call end_netcdf(1)

      nvar=1
      varname(1,1)='icefrac'
      vdims(1)=2
      vadims(1,1)=loc_dim('longitude',dimname,nall)
      vadims(2,1)=loc_dim('latitude',dimname,nall)
      nattsvar(1)=2
      attvarname(1,1,1)='long_name'
      attvarname(2,1,1)='ice fraction (no Greenland)'
      attvarname(1,2,1)='units'
      attvarname(2,2,1)='0-1'
      call ininc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/output/'//
     &                    'icefrac_nogreen_t21.nc',
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))
      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)
      call writevar(nc(1),
     :        idvar(loc_dim('icefrac',varname(1,1),nall),1),
     :        ice_new1 )
      call end_netcdf(1)

      nvar=1
      varname(1,1)='icefrac'
      vdims(1)=2
      vadims(1,1)=loc_dim('longitude',dimname,nall)
      vadims(2,1)=loc_dim('latitude',dimname,nall)
      nattsvar(1)=2
      attvarname(1,1,1)='long_name'
      attvarname(2,1,1)='ice fraction (no icesheets)'
      attvarname(1,2,1)='units'
      attvarname(2,2,1)='0-1'
      call ininc('/home/ggdjl/genie/genie-utilities/'//
     &                        'data/output/'//
     &                    'icefrac_noice_t21.nc',
     :     nmaxdims,ndim,nvar,
     :     natts,nattsvar,
     :     vdims,vadims,ndims,
     :     dimname(1,1),varname(1,1),
     :     attdimname,attvarname,
     :     nc(1),iddim(1,1),idvar(1,1))
      call writedim(nc(1),iddim(1,1),alon1)
      call writedim(nc(1),iddim(2,1),alat1)
      call writevar(nc(1),
     :        idvar(loc_dim('icefrac',varname(1,1),nall),1),
     :        ice_new2 )
      call end_netcdf(1)


      end
