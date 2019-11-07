*
* subroutine tstep.f for program goldstein 
* error in order of co call and rho b.c. corrected 9/2/1
* flux version fully explicit one step second order variable depth
*
c version with isoneutral diffusion 21/5/1
c notes; fe/fn/fa in previous versions could have been scalars
c error in isoneutral coeffs corrected 29/8/2
c upstream weighting included 29/8/2
c variable ds 6/12/4 nre
c 
c AY (12/07/05) : removed surplus input argument to function
c
      subroutine tstepo

#include "ocean.cmn"
      
      real tv, ups(3), ups0, pec(3)

      real fe(maxl), fw(maxl), fn(maxl), fs(maxl,maxi), fa(maxl)
     1   , fb(maxl,maxi,maxj), fwsave(maxl)
      integer i, j, k, l
c      integer istep

c     logical vdt
c     parameter(vdt = .false.)
      parameter(ups0 = 0.)

c KICO Mixed layer scheme needs:
      real mldtsold(maxl,maxi,maxj,maxk), mldrhoold(maxi,maxj,maxk)
      real mldrhonew(maxi,maxj,maxk)
      real mldtstmp(2), mldrhotmp
      integer mldpk(2,maxi,maxj)

c KICO ediff calc needs
      real diffv

      real tec, scc, dzrho, rdzrho, slim, tv1
      real dxrho(4), dxts(maxl,4), dyrho(4), dyts(maxl,4), dzts(maxl)
      integer ina, nnp, knp

      if(diso)then
         scc = ec(2)
         limps = 0
      endif

      if(imld.eq.1)then
c KICO 18/09/07, before calculating fluxes, calculate energy consumed
c or released in mixing surface forcing over top layer. Needed for
c mld calculation, especially if mixed layer is <1 cell thick.
c NOTE: for now, all "energies" in this scheme are calculated in units
c of density x height^2, or (Energy/area)/2*acceleration (due to gravity).
        do j=1,jmax
          do i=1,imax
            mldtstmp(1) = ts(1,i,j,kmax)-ts(1,i,j,kmax+1)
            mldtstmp(2) = ts(2,i,j,kmax)-ts(2,i,j,kmax+1)
            call eos(ec,mldtstmp(1),mldtstmp(2),zro(kmax),ieos,
     1         mldrhotmp)
            mldpelayer1(i,j) = (mldrhotmp-rho(i,j,kmax))
     1                         *z2dzg(kmax,kmax)
          enddo
        enddo
      endif
c KICO initialise diffv
      diffv=diff(2)

c     if(mod(istep,100 ).eq.1)
c     call velc

      dmax = 0

* find timestep 

c     if(vdt)then 
c     endif

* 2nd order explicit step

c lower boundary fluxes

      do 220 j=1,jmax
         do 220 i=1,imax
            do 220 l=1,lmax
               fb(l,i,j) = 0
  220 continue

      do 100 k=1,kmax

c southern boundary fluxes

         j = 1
         do 230 i=1,imax
            do 230 l=1,lmax
               fs(l,i) = 0
  230    continue
         do 100 j=1,jmax
c western boundary fluxes
            i = 1
c           ups(1) = sign(ups0, u(1,imax,j,k))
            pec(1) = u(1,imax,j,k)*dphi/diff(1)
            ups(1) = pec(1) / (2.0 + abs(pec(1)))
            do 210 l=1,lmax
               if (k.ge.max(k1(imax,j),k1(1,j)))then
c western doorway
                  fw(l) = u(1,imax,j,k)*rc(j)*((1.-ups(1))*ts1(l,1,j,k)
     1                  + (1.+ups(1))*ts1(l,imax,j,k))*0.5
                  fw(l) = fw(l) - (ts1(l,1,j,k) - ts1(l,imax,j,k))
     1                       *rc2(j)*diff(1)
               else
                  fw(l) = 0
               endif
               fwsave(l) = fw(l)
  210       continue
            do 100 i=1,imax

c KICO calculate local vertical diffusivity
              if(k.ge.k1(i,j).and.k.lt.kmax)then
c First get vertical density gradient (also needed for isopycnal diff)
                call eosd(ec,ts1(1,i,j,k:k+1),ts1(2,i,j,k:k+1),zw(k),
     1             rdza(k),ieos,dzrho,tec)
                if(dzrho.lt.-1e-12)then
                  rdzrho = 1.0/dzrho
                else
                  rdzrho = -1e12
                endif
                if((iediff.gt.0).AND.(iediff.lt.3))then
c Value of diffv fine for applying diffusivity, but peclet number calc
c is a 1st order approximation if diffusivity is variable.
                   if(ediffpow2i.eq.0)then
                     diffv=ediff0+ediff1(i,j,k)
                   elseif(ediffpow2i.eq.1)then
                     diffv=ediff0+ediff1(i,j,k)*(-rdzrho)
                   elseif(ediffpow2i.eq.2)then
                     diffv=ediff0+ediff1(i,j,k)*sqrt(-rdzrho)
                   else
                     diffv=ediff0+ediff1(i,j,k)*((-rdzrho)**ediffpow2)
                   endif
                   if(diffv.gt.diffmax(k+1)) diffv=diffmax(k+1)
                 endif
c                 if((i.eq.22).OR.(i.eq.10))then
c                   if(j.eq.24)then
c                     print*,i,j,k,'diff(2)',diff(2),'ediff1',
c     1                 ediff1(i,j,k),'diffv',diffv,
c     2                 'diffmax',diffmax(k+1)
c                   endif
c                 endif
               endif

c              ups(1) = sign(ups0, u(1,i,j,k))
               pec(1) = u(1,i,j,k)*dphi/diff(1)
               ups(1) = pec(1) / (2.0 + abs(pec(1)))
c              ups(2) = sign(ups0, u(2,i,j,k))
c rather untidy mask to avoid undefined dsv at jmax nre
               pec(2) = u(2,i,j,k)*dsv(min(j,jmax-1))/diff(1)
               ups(2) = pec(2) / (2.0 + abs(pec(2)))
c              ups(3) = sign(ups0, u(3,i,j,k))
               pec(3) = u(3,i,j,k)*dza(k)/diffv
               ups(3) = pec(3) / (2.0 + abs(pec(3)))
               do l=1,lmax
c flux to east
                  if(i.eq.imax)then
c eastern edge(doorway or wall)
                     fe(l) = fwsave(l)
                  elseif(k.lt.max(k1(i,j),k1(i+1,j)))then
                     fe(l) = 0
                  else
                     fe(l) = u(1,i,j,k)*rc(j)*((1.-ups(1))
     &                   *ts1(l,i+1,j,k) + (1.+ups(1))*ts1(l,i,j,k))*0.5
                     fe(l) = fe(l) - (ts1(l,i+1,j,k) - ts1(l,i,j,k))
     1                       *rc2(j)*diff(1)
                  endif
c flux to north
                  if(k.lt.max(k1(i,j),k1(i,j+1)))then
                     fn(l) = 0
                  else
                     fn(l) = cv(j)*u(2,i,j,k)*((1.-ups(2))
     &                   *ts1(l,i,j+1,k) + (1.+ups(2))*ts1(l,i,j,k))*0.5
                     fn(l) = fn(l) - cv2(j)*(ts1(l,i,j+1,k) -
     1                       ts1(l,i,j,k))*diff(1)
                  endif
c flux above
                  if(k.lt.k1(i,j))then
                     fa(l) = 0
                  elseif(k.eq.kmax)then
                     fa(l) = ts(l,i,j,kmax+1) 
                  else
                     fa(l) = u(3,i,j,k)*((1.-ups(3))*ts1(l,i,j,k+1) +
     1                       (1.+ups(3))*ts1(l,i,j,k))*0.5
                     fa(l) = fa(l) - (ts1(l,i,j,k+1) - ts1(l,i,j,k))
     1                       *rdza(k)*diffv
                  endif
               enddo
               if(diso)then
c isoneutral diffusion
                  if(k.ge.k1(i,j).and.k.lt.kmax)then
                     if(dzrho.lt.-1e-12)then
                        tv1 = 0.0
c tracer loop          
                        do knp=0,1
                           do nnp=0,1
                              ina = 1+nnp + 2*knp
c phi derivatives
                              do l=1,lmax
                                 if(k+knp.ge.k1(i-1+2*nnp,j))then
                                    dxts(l,ina) = (ts1(l,i+nnp,j,k+knp)
     2                                   - ts1(l,i+nnp-1,j,k+knp))
     3                                   *rc(j)*rdphi
                                 else
                                    dxts(l,ina) = 0.
                                 endif
c s-derivatives
                                 if(k+knp.ge.k1(i,j-1+2*nnp))then
                                    dyts(l,ina) = (ts1(l,i,j+nnp,k+knp)
     2                                   - ts1(l,i,j+nnp-1,k+knp))
     &                                   *cv(j-1+nnp)*rdsv(j+nnp-1)
                                 else
                                    dyts(l,ina) = 0.
                                 endif
                              enddo
                              dxrho(ina) = scc*dxts(2,ina)
     1                             -tec*dxts(1,ina)
                              dyrho(ina) = scc*dyts(2,ina)
     1                             -tec*dyts(1,ina)
c calculate diagonal part
                              tv1 = tv1 + dxrho(ina)*dxrho(ina)
     1                             + dyrho(ina)*dyrho(ina)
                           enddo
                        enddo
                        tv1 = 0.25*tv1*rdzrho*rdzrho
c limit flux by factor slim for large slope
                        if(tv1.gt.ssmax(k))then
c                       slim = ssmax(k)/tv1
                           slim = ssmax(k)*ssmax(k)/(tv1*tv1)
c count flux-limited points
                           limps = limps + 1
                        else
                           slim = 1.0
                        endif
                        tv1 = tv1*slim*diff(1)*rdza(k)
c test vertical diffusion number
                        tv = tv1*dt(k)*rdza(k)
                        if(tv.gt.dmax)then
                           dmax = tv
                        endif
                        do l=1,lmax
                           dzts(l) = (ts1(l,i,j,k+1) 
     1                          - ts1(l,i,j,k))*rdza(k)
c add isoneutral vertical flux
                           tv = 0
                           do ina=1,4
                              tv = tv + (2*dzrho*dxts(l,ina)
     1                             - dxrho(ina)*dzts(l))*dxrho(ina) 
     2                             + (2*dzrho*dyts(l,ina)
     3                             - dyrho(ina)*dzts(l))*dyrho(ina)
                           enddo
                           tv = 0.25*slim*diff(1)*tv/(dzrho*dzrho)
                           fa(l) = fa(l) + tv
c                       if(abs(tv).gt.1e0)
c    1                     write(6,'(3i4,8e14.4)')i,j,k,fa(l),tv
                        enddo
                     endif
                  endif
               endif
               do l=1,lmax
                  tv = 0
                  if(k.ge.k1(i,j))then
                     ts(l,i,j,k) = ts1(l,i,j,k) - dt(k)*( - tv +
     1                           (fe(l) - fw(l))*rdphi
     2                           + (fn(l) - fs(l,i))*rds(j)
     3                           + (fa(l) - fb(l,i,j))*rdz(k))
                  endif

                  fw(l) = fe(l)
                  fs(l,i) = fn(l)
                  fb(l,i,j) = fa(l)
               enddo

               call eos(ec,ts(1,i,j,k),ts(2,i,j,k),zro(k),
     1            ieos,rho(i,j,k))
  100 continue

       if(imld.eq.1) then
c KICO 18/09/07, remember old T and S only so we can calculate PE 
c change. Code inefficient at present because only surface mixed points
c are compared (energy released due to eg bottom mixed layers is
c excluded, although surface plumes from coshuffle are accounted for),
c but T and S remembered everywhere.
         do i=1,imax
           do j=1,jmax
             if (k1(i,j).le.kmax)then
             do k=1,kmax
               do l=1,2
                 mldtsold(l,i,j,k)=ts(l,i,j,k)
               enddo
             enddo
             endif
           enddo
         enddo
c KICO 18/09/07 convection scheme itself not changed, but diagnostic
c for mixed layer depth precursor diagnosed in mldpk(1,:,:). If
c coshuffle used, maxmium depth of surface plume is in mldpk(2,:,:)
        call co(ts,mldpk)
c KICO 18/09/07 calculate PE released and multiply it by efficiency
c of recycling of PE. Need to check non-linearities properly
c represented before any commit??!!
        do i=1,imax
          do j=1,jmax
             if (k1(i,j).le.kmax)then
               mldpeconv(i,j)=0
               k=kmax
c               do while(k.ge.mldpk(iconv+1,i,j))
               do while(k.gt.0)
                 call eos(ec,mldtsold(1,i,j,k),mldtsold(2,i,j,k),
     1             zro(k),ieos,mldrhoold(i,j,k))
                 call eos(ec,ts(1,i,j,k),ts(2,i,j,k),
     1             zro(k),ieos,mldrhonew(i,j,k))
                 mldpeconv(i,j)=mldpeconv(i,j)+(mldrhonew(i,j,k)
     1                        -mldrhoold(i,j,k))*z2dzg(k,k)
                 k=k-1
               enddo
               mldpebuoy(i,j)=mldpeconv(i,j)+mldpelayer1(i,j)
               if (mldpebuoy(i,j).gt.0.)then
                 mldpebuoy(i,j)=mldpebuoy(i,j)*mldpebuoycoeff
               endif
c Add wind energy
               mldemix(i,j)=mldpebuoy(i,j)+mldketau(i,j)*mlddec(kmax)
             endif
          enddo
        enddo

c KICO 18/09/07  Kraus Turner scheme added. Static instability
c driven convection  has already been done in co, and will differ from
c standard KT if iconv.eq.1. krausturner uses PE released
c in co and KE from the wind to deepen the mixed layer further.
        do i=1,imax
          do j=1,jmax
             if (k1(i,j).le.kmax)then
               if (mldemix(i,j).gt.0.)then
c Apply krausturner
                 call krausturner(ts(1:lmax,i,j,1:kmax),mldpebuoy(i,j),
     1            mldketau(i,j),mldpk(1,i,j),mld(i,j),mldk(i,j),k1(i,j))
c                print*,'Y',i,j,k1(i,j),mldk(i,j),
c     1             'mld',mld(i,j),'mldpeconv',
c     2             mldpeconv(i,j),'mldpelayer1',mldpelayer1(i,j),
c     3             'mldpebuoy',mldpebuoy(i,j),'mldketau',
c     4             mldketau(i,j),'mldemix',mldemix(i,j),'ts',
c     5             ts(1,i,j,mldk(i,j)),ts(1,i,j,mldk(i,j)+1),
c     6             ts(1,i,j,kmax),ts(1,i,j,kmax+1),
c     7             ts(2,i,j,mldk(i,j)),
c     8             ts(2,i,j,mldk(i,j)+1),ts(2,i,j,kmax),
c     9             ts(2,i,j,kmax+1)
               else
c Not enough energy even to homogenise first layer. The first layer *is*
c still homogeneous for all tracers, but an mld shallower than the
c first cell is output as a diagnostic
                 mldk(i,j)=kmax
                 if(mldpelayer1(i,j).lt.0)then
                  mld(i,j)=zw(kmax-1)*(1-mldemix(i,j)/mldpelayer1(i,j))
                 else
                   mld(i,j)=zw(kmax-1)
                 endif
c                print*,'N',i,j,k1(i,j),mldk(i,j),
c     1             'mld',mld(i,j),'mldpeconv',
c     2             mldpeconv(i,j),'mldpelayer1',mldpelayer1(i,j),
c     3             'mldpebuoy',mldpebuoy(i,j),'mldketau',
c     4             mldketau(i,j),'mldemix',mldemix(i,j),'ts',
c     5             ts(1,i,j,mldk(i,j)),ts(1,i,j,mldk(i,j)+1),
c     6             ts(1,i,j,kmax),ts(1,i,j,kmax+1),
c     7             ts(2,i,j,mldk(i,j)),
c     8             ts(2,i,j,mldk(i,j)+1),ts(2,i,j,kmax),
c     9             ts(2,i,j,kmax+1)
               endif
             endif
           enddo
         enddo

      else
c Not applying mixed layer scheme. Just call convection scheme
        call co(ts,mldpk)
      endif

c KICO, if thermobaricity is on, make sure rho calculation is vertically local
      if(ieos.ne.0)then
         do i=1,imax
           do j=1,jmax
             if (k1(i,j).le.kmax)then
               do k=1,kmax
                 call eos(ec,ts(1,i,j,k),ts(2,i,j,k),zro(k),
     1              ieos,rho(i,j,k))
               enddo
             endif
           enddo
         enddo
      endif

c periodic b.c. for rho (required at wet points)
c isoneutral code also needs ts1 bc.

      do j=1,jmax
         do k=k1(0,j),kmax
            rho(0,j,k) = rho(imax,j,k)
            do l=1,lmax
               ts1(l,0,j,k) = ts(l,imax,j,k)
            enddo
         enddo
         do k=k1(imax+1,j),kmax
            rho(imax+1,j,k) = rho(1,j,k)
            do l=1,lmax
               ts1(l,imax+1,j,k) = ts(l,1,j,k)
            enddo
         enddo
      enddo

      do 10 k=1,kmax
         do 10 j=1,jmax
            do 10 i=1,imax
               do 10 l=1,lmax
                  if(k.ge.k1(i,j)) ts1(l,i,j,k) = ts(l,i,j,k)
   10 continue

      end
