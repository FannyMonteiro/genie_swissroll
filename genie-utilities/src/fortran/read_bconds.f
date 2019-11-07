c    This routine readsdata from a file, and then 
c    peforms an interpolation onto a new grid.
c    the type of interpolation depends on the value of 'procresstype':
c    1 = Straight interpolation
c    2 = Remove negative values, then interpolate
c    3 = Produce a 0-1 mask, then interpolate
c    4 = no interpolation (input,output grid are identical)
c    5 = make -99999 equal 1, else zero, then interpolate
c          (useful for converting goldstein output)
c    6 = make anything more negative than -9e36 equal zero, then interpolate
c          (useful for ncep runoff)
c    7 = make anything greater than 1e20 equal zero, then interpolate
c          (useful for Hadam data)


      subroutine read_bconds(gettype,nmonths,processtype,
     &      dataout,in_filename,
     &      in_lonname,in_latname,in_varname)

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
#include "/home/ggdjl/genie/genie-igcm3/src/fortran/param1.cmn"
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/param2.cmn'

c     DATA:
      real,allocatable,dimension(:) :: in_lons
      real,allocatable,dimension(:) :: in_lats
      real,allocatable,dimension(:) :: in_lonsedge
      real,allocatable,dimension(:) :: in_latsedge
      real,allocatable,dimension(:,:,:) :: datain_12
      real,allocatable,dimension(:,:) :: datain_1
      integer,allocatable,dimension(:,:) :: work

      integer in_nlon
      integer in_nlat

      integer sign

      integer nmonths
  
      integer gettype,processtype

c     NETCDF + AWI STUFF: 
      integer ncid,ifail,loc_dim
      integer ier,dimid,a

c     LOOPING:
      integer i,j,m,v

c     IGCM GRID:
      integer ilon1_atm,ilat1_atm
      parameter (ilon1_atm=64,ilat1_atm=32)
      real alon1(ilon1_atm),alat1(ilat1_atm),
     :     aboxedge1_lon(ilon1_atm+1),aboxedge1_lat(ilat1_atm+1)

      integer work_igcm(ilon1_atm,ilat1_atm)
      real dataout(ilon1_atm,ilat1_atm,nmonths)
      real dataout_1_igcm(ilon1_atm,ilat1_atm)

c     FOR GRID:
      real ax

      character in_filename*200
      character in_lonname*200
      character in_latname*200
      character in_varname*200


c     COPIED FROM INITIALISE_ATMOS.F
      ax=360.0/real(mg)
      do i=1,ilon1_atm
         aboxedge1_lon(i)=(i-1.5)*ax
         alon1(i)=(i-1.0)*ax
      end do
      aboxedge1_lon(mg+1)=(mg-0.5)*ax
      call gwtcnr(alat1,jg)
      call gwtbox(aboxedge1_lat,jg)

      print*,'igcm grid set up OK in read_bconds'



      call open_file_nc(trim(in_filename),ncid)
      ier=nf_inq_dimid(ncid,trim(in_lonname),dimid)
      ier=nf_inq_dimlen(ncid,dimid,in_nlon)
      ier=nf_inq_dimid(ncid,trim(in_latname),dimid)
      ier=nf_inq_dimlen(ncid,dimid,in_nlat)

      print*,'size:',in_nlon,in_nlat

      allocate(in_lons(in_nlon))      
      allocate(in_lats(in_nlat))
      allocate(in_lonsedge(in_nlon+1))
      allocate(in_latsedge(in_nlat+1))

      allocate(datain_12(in_nlon,in_nlat,nmonths))
      allocate(datain_1(in_nlon,in_nlat))
      allocate(work(in_nlon,in_nlat))
      work(:,:)=1

      print*,'reading data...'

      if (gettype.eq.2) then
      call get3d_data_nc(ncid,trim(in_varname),in_nlon,in_nlat,
     &     nmonths,datain_12,ifail)
      
      else if (gettype.eq.1) then
      call get2d_data_nc(ncid,trim(in_varname),in_nlon,in_nlat,
     &     datain_12,ifail)
      endif

      call get1d_data_nc(ncid,trim(in_lonname),in_nlon,in_lons,ifail)
      call get1d_data_nc(ncid,trim(in_latname),in_nlat,in_lats,ifail)
      call close_file_nc(trim(in_filename),ncid)

      print*,'data read in from file OK'

      do i=2,in_nlon
        in_lonsedge(i)=(in_lons(i-1)+in_lons(i))/2.
      enddo
      in_lonsedge(1)=(in_lons(1)+in_lons(in_nlon)-360)/2.
      in_lonsedge(in_nlon+1)=in_lonsedge(1)+360.

      do j=2,in_nlat
        in_latsedge(j)=(in_lats(j-1)+in_lats(j))/2.
      enddo
      sign=nint((in_lats(1)-in_lats(in_nlat))/
     &        abs(in_lats(1)-in_lats(in_nlat)))
      in_latsedge(1)=90.*sign
      in_latsedge(in_nlat+1)=-90.*sign


      do m=1,nmonths

         if (processtype.eq.2) then
            datain_12(:,:,m)=max(datain_12(:,:,m),0.0)
         endif

         if (processtype.eq.3) then
            do j=1,in_nlat
            do i=1,in_nlon
              if (datain_12(i,j,m).gt.0.0) then 
                datain_12(i,j,m)=1.0
              else
                datain_12(i,j,m)=0.0
              endif
            enddo
            enddo
         endif

         if (processtype.eq.5) then
            do j=1,in_nlat
            do i=1,in_nlon
              if (datain_12(i,j,m).eq.-99999) then 
                datain_12(i,j,m)=1.0
              else
                datain_12(i,j,m)=0.0
              endif
            enddo
            enddo
         endif

         if (processtype.eq.6) then
            do j=1,in_nlat
            do i=1,in_nlon
              if (datain_12(i,j,m).le.-9e36) then 
                datain_12(i,j,m)=0.0
              endif
            enddo
            enddo
         endif

         if (processtype.eq.7) then
            do j=1,in_nlat
            do i=1,in_nlon
              if (datain_12(i,j,m).ge.1e20) then 
                datain_12(i,j,m)=0.0
              endif
            enddo
            enddo
         endif

         datain_1(:,:)=datain_12(:,:,m)
         work(:,:)=1
         work_igcm(:,:)=1


         if (processtype.ne.4) then

         call awi(in_nlon,in_nlon,in_lonsedge,
     :             in_nlat,in_latsedge,datain_1,work,
     :             ilon1_atm,ilon1_atm,aboxedge1_lon,
     :             ilat1_atm,aboxedge1_lat,dataout_1_igcm,work_igcm,
     :             ifail)
         if (ifail.ne.0) then
            print*,'Interp_to_ocn:  Error in awi,   itype 1 ',ifail
            stop 1
         end if
         else
         dataout_1_igcm(:,:)=datain_1(:,:)
         endif

         dataout(:,:,m)=dataout_1_igcm(:,:)

      enddo

      
      print*,'interpolation OK'

      print*,'end of read_bconds'

      deallocate(in_lons)
      deallocate(in_lats)
      deallocate(in_lonsedge)
      deallocate(in_latsedge)
      deallocate(datain_12)
      deallocate(datain_1) 
      deallocate(work) 




      end
