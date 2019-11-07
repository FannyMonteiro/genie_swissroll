      program make_runoff

c***********************************************************

c     This code can be used to create a runoff field of nearest neightbours,
c       and another field which is like Phil's runoff mask but with no endorheic regions.

c     I compile with: 
c     ifort -o make_runoff_NOCVS make_runoff.f -lnc1 -lutil1 -lnetcdf -ligcm3 -L/home/ggdjl/genie/genie-lib/libnc1 -L/home/ggdjl/genie/genie-lib/libutil1 -L/home/ggdjl/genie/genie-igcm3/lib -L/opt/local/intel_fc_81/lib -fpp

c     The lines between the asterixes will need to be edited.

c     Please edit the main body of the code also, to tidy it up!

c     Eventually, this piece of code could run with makefile.arc!

      implicit none
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netcdf.inc'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/netdata.cmn'
#include '/home/ggdjl/genie/genie-igcm3/src/fortran/param1.cmn'
      include '/home/ggdjl/genie/genie-igcm3/src/fortran/param2.cmn'

c     This is the standard land-sea mask:
c      include 'names_runoff.inc' 

c     This is the Peltier land-sea mask:
c      include 'names_runoff_pelt.inc' 

c     This is the goldstein land-sea mask:
c      include 'names_runoff_goldstein.inc'      

c     This is the goldstein 64x32 land-sea mask:
      include 'names_runoff_goldstein6432.inc'

c***********************************************************

c     MONTHS:
      integer nmonths
      parameter(nmonths=2)
      real months(nmonths)

c     IGCM GRID:
      integer ilon1_atm,ilat1_atm
      parameter (ilon1_atm=64,ilat1_atm=32)
      real alon1(ilon1_atm),alat1(ilat1_atm),
     :     aboxedge1_lon(ilon1_atm+1),aboxedge1_lat(ilat1_atm+1)

c     DATA:

      real runoff_1_igcm(ilon1_atm,ilat1_atm,2)
      real runoff_2_igcm(ilon1_atm,ilat1_atm,2)
      real array_igcm(ilon1_atm,ilat1_atm,2)
      real lsm_1_igcm(ilon1_atm,ilat1_atm)

      integer work_igcm(ilon1_atm,ilat1_atm)

      real datain_1_igcm(ilon1_atm,ilat1_atm)

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
      integer i,j,m,v

c     FOR GRID:
      real ax

c     FOR PROCESSING
      integer xlim_1
      parameter(xlim_1=1)

      integer xlim_2
      parameter(xlim_2=5)

      integer ii,jj,itemp,jtemp
      integer x

      print*,'make_runoff is starting....'

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


      call open_file_nc(trim(in_data_maskname),ncid)
      call get2d_data_nc(ncid,trim(in_mask_varname),ilon1_atm,ilat1_atm,
     &          datain_1_igcm,ifail)
      if (ifail.ne.0) stop 'Problem reading input mask file'
      lsm_1_igcm(:,:)=datain_1_igcm(:,:)
      call close_file_nc(trim(in_data_maskname),ncid)

      call open_file_nc(trim(in_data_runoffname),ncid)
      call get3d_data_nc(ncid,trim(in_runoff_varname),
     &          ilon1_atm,ilat1_atm,2,
     &          runoff_1_igcm,ifail)
      if (ifail.ne.0) stop 'Problem reading input runoff file'
      call close_file_nc(trim(in_data_runoffname),ncid)

      print*,'all data read in'


c     process the data...............

      runoff_2_igcm(:,:,:)=0.0

      do j=1,ilat1_atm
        do i=1,ilon1_atm

        array_igcm(i,j,1)=i
        array_igcm(i,j,2)=j

          if (lsm_1_igcm(i,j).eq.1) then


c     Nearest-neighbour mask
            do x=1,xlim_2
              if (runoff_2_igcm(i,j,1).eq.0) then
                do ii=-x,x
                  do jj=-x,x
                    itemp=ii+i
                    if (itemp.gt.ilon1_atm) itemp=itemp-ilon1_atm
                    if (itemp.lt.1) itemp=itemp+ilon1_atm
                    jtemp=jj+j
                    if (jtemp.gt.ilat1_atm) jtemp=ilat1_atm
                    if (jtemp.lt.1) jtemp=1
                    if (lsm_1_igcm(itemp,jtemp).eq.0) then
                      runoff_2_igcm(i,j,1)=itemp
                      runoff_2_igcm(i,j,2)=jtemp
                    endif
                  enddo
                enddo
              endif
            enddo


c     Phil's runoff mask.....
            do x=1,xlim_1
              if (runoff_1_igcm(i,j,1).eq.0) then
                do ii=-x,x
                  do jj=-x,x
                    itemp=ii+i
                    if (itemp.gt.ilon1_atm) itemp=itemp-ilon1_atm
                    if (itemp.lt.1) itemp=itemp+ilon1_atm
                    jtemp=jj+j
                    if (jtemp.gt.ilat1_atm) jtemp=ilat1_atm
                    if (jtemp.lt.1) jtemp=1
                    if (runoff_1_igcm(itemp,jtemp,1).ne.0) then
                      runoff_1_igcm(i,j,1)=runoff_1_igcm(itemp,jtemp,1)
                      runoff_1_igcm(i,j,2)=runoff_1_igcm(itemp,jtemp,2)
                    endif
                  enddo
                enddo
              endif
            enddo

          else

            runoff_1_igcm(i,j,1)=0.0
            runoff_1_igcm(i,j,2)=0.0
            runoff_2_igcm(i,j,1)=0.0
            runoff_2_igcm(i,j,2)=0.0

          endif

        enddo
      enddo

c     check the data

      do j=1,ilat1_atm
        do i=1,ilon1_atm

          if (lsm_1_igcm(i,j).eq.1) then

            if ((runoff_1_igcm(i,j,1).eq.0).or.
     :         (runoff_1_igcm(i,j,2).eq.0)) then
              print*,'problem a) in 1 at ',i,j
c              stop
            endif

            if ((runoff_2_igcm(i,j,1).eq.0).or.
     :         (runoff_2_igcm(i,j,2).eq.0)) then
              print*,'problem a) in 2 at ',i,j
c              stop
            endif

            if (lsm_1_igcm(runoff_1_igcm(i,j,1),
     :                      runoff_1_igcm(i,j,2)).ne.0) then
              print*,'problem b) in 1 at ',i,j
c              stop
            endif

            if (lsm_1_igcm(runoff_2_igcm(i,j,1),
     :                      runoff_2_igcm(i,j,2)).ne.0) then
              print*,'problem b) in 2 at ',i,j
c              stop
            endif

          else

            if ((runoff_1_igcm(i,j,1).ne.0).or.
     :         (runoff_1_igcm(i,j,2).ne.0)) then
              print*,'problem c) in 1 at ',i,j
c              stop
            endif

            if ((runoff_2_igcm(i,j,1).ne.0).or.
     :         (runoff_2_igcm(i,j,2).ne.0)) then
              print*,'problem c) in 2 at ',i,j
c              stop
            endif

          endif



        enddo
      enddo


c     write the data

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

      dimname(3,1)='ij_index'
      ndims(3)=2
      natts(3)=2
      attdimname(1,1,3)='long_name'
      attdimname(2,1,3)='Index Number'
      attdimname(1,2,3)='units'
      attdimname(2,2,3)='no units'

      nvar=2

      v=1
      varname(v,1)='destination_indices'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('ij_index',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='Destination index for runoff'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='no units'

      v=2
      varname(v,1)='simple_array'
      vdims(v)=3
      vadims(1,v)=loc_dim('longitude',dimname,nall)
      vadims(2,v)=loc_dim('latitude',dimname,nall)
      vadims(3,v)=loc_dim('ij_index',dimname,nall)
      nattsvar(v)=2
      attvarname(1,1,v)='long_name'
      attvarname(2,1,v)='Destination index for runoff'
      attvarname(1,2,v)='units'
      attvarname(2,2,v)='no units'


      call ininc(out_data_filename_1,
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
     :        runoff_1_igcm )
      v=2
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        array_igcm )
      print*,'written file'
      call end_netcdf(1)


      call ininc(out_data_filename_2,
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
     :        runoff_2_igcm )
      v=2
      call writevar(nc(1),
     :        idvar(loc_dim(varname(v,1),varname(1,1),nall),1),
     :        array_igcm )
      print*,'written file'
      call end_netcdf(1)







      print*,'files written OK'

      end


