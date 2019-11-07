pro make

; This bit of code makes text files which contain co2
;   and d18o from the vostok ice-core, from 30kyrBP to 
;   pre-industrial.
; It uses the original Vostok data from petit et al, and interpolates
; it to a value for each kyr.  The co2 data is capped at 279ppmv
; The d18o data is first smoothed using a boxcar average, of width 5 points.


close,1
device,decomposed=0
loadct,39
nt=31
nt_hr=301
nt1=301
nj=156
nj1=91
nd1=283
nd2=318
nd3=183
np1=17
np2=80
co2_max=279
junk=strarr(nj,1)
junk1=strarr(nj1,1)

nsnap=11
snapname=strarr(nsnap,3)
snapname(0,*)='st3'
snapname(1,*)='27e'
snapname(2,*)='24e'
snapname(3,*)='lgm'
snapname(4,*)='18e'
snapname(5,*)='15e'
snapname(6,*)='12e'
snapname(7,*)='09e'
snapname(8,*)='hol'
snapname(9,*)='03e'
snapname(10,*)='pre'
snapnum=intarr(nsnap)
snapnum=[30,27,24,21,18,15,12,9,6,3,0]



co2_arr=fltarr(2,nd1)
co2_arr_hr=fltarr(4,nd3)

openr,1,'/home/ggdjl/paleodata/ice_cores/vostok_co2.dat'
readf,1,junk
readf,1,co2_arr
close,1

openr,1,'/home/ggdjl/paleodata/ice_cores/domec_co2.dat'
readf,1,junk1
readf,1,co2_arr_hr
close,1

co2_index=findgen(nt)*1000.0
co2_index_hr=findgen(nt_hr)*100.0

co2_vect=interpol(co2_arr(1,*),co2_arr(0,*),co2_index)
co2_vect_hr=interpol(co2_arr_hr(2,*),co2_arr_hr(1,*),co2_index_hr)

for x=0,nt-1 do begin
  if (co2_index(x) lt co2_arr(0,0)) then begin
    co2_vect(x)=co2_arr(1,0)
  endif
endfor

for x=0,nt_hr-1 do begin
  if (co2_index_hr(x) lt co2_arr_hr(1,0)) then begin
    co2_vect_hr(x)=co2_arr_hr(2,0)
  endif
endfor

co2_vect=co2_vect*(co2_vect le co2_max)+(co2_vect gt co2_max)*co2_max
;co2_vect_hr=co2_vect_hr*(co2_vect_hr le co2_max)+(co2_vect_hr gt co2_max)*co2_max

co2_vect_hrl=interpol(co2_vect,co2_index,co2_index_hr)
co2_vect_hlr=interpol(co2_vect_hr,co2_index_hr,co2_index)
co2_vect_hrr=interpol(co2_vect_hlr,co2_index,co2_index_hr)
co2_vect_new=co2_vect_hr-co2_vect_hrr+co2_vect_hrl

co2_vect_new(210:300)=co2_vect_hrl(210:300)

window,0
plot,co2_arr(0,0:np1),co2_arr(1,0:np1),ystyle=1,yrange=[170,330],thick=3,/nodata
oplot,co2_index,co2_vect,color=50,thick=3
oplot,co2_index_hr,co2_vect_hr,color=100,thick=3
oplot,co2_index_hr,co2_vect_new,color=250,thick=3



openw,1,'/home/ggdjl/genie/genie-utilities/data/output/co2_30k.dat'
printf,1,reverse(co2_vect)/1e6,format='(f8.6)'
close,1
openw,1,'/home/ggdjl/genie/genie-utilities/data/output/co2_21k.dat'
printf,1,reverse(co2_vect(0:21))/1e6,format='(f8.6)'
close,1
openw,1,'/home/ggdjl/genie/genie-utilities/data/output/co2_30k_highres.dat'
printf,1,reverse(co2_vect_new)/1e6,format='(f9.7)'
close,1


for s=0,nsnap-1 do begin
openw,1,'/home/ggdjl/genie/genie-utilities/data/output/co2_'+snapname(s)+'.dat'
printf,1,co2_vect(snapnum(s))/1e6,format='(f8.6)'
close,1
endfor


d18o_arr=fltarr(2,nd2)
openr,1,'/home/ggdjl/paleodata/ice_cores/vostok_o18.dat'
readf,1,junk
readf,1,d18o_arr
close,1
d18o_orig=d18o_arr(1,*)
d18o_arr(1,*)=smooth(d18o_arr(1,*),5,/edge_truncate)
d18o_index=findgen(nt)*1000.0
d18o_index1=findgen(nt1)*100.0
d18o_vect=interpol(d18o_arr(1,*),d18o_arr(0,*),d18o_index)
d18o_vect1=interpol(d18o_orig,d18o_arr(0,*),d18o_index1)
d18o_vect2=interpol(d18o_arr(1,*),d18o_arr(0,*),d18o_index1)

;d18o_vect1(210:300)=d18o_vect2(210:300)

for x=0,nt-1 do begin
  if (d18o_index(x) lt d18o_arr(0,0)) then begin
    d18o_vect(x)=d18o_arr(1,0)
  endif
endfor
for x=0,nt1-1 do begin
  if (d18o_index1(x) lt d18o_arr(0,0)) then begin
    d18o_vect1(x)=d18o_arr(1,0)
  endif
endfor
window,1
plot,d18o_arr(0,0:np2),-1*d18o_arr(1,0:np2),ystyle=1,yrange=[-1.5,0.5],thick=3
oplot,d18o_index,-1*d18o_vect,color=50,thick=3
oplot,d18o_arr(0,0:np2),-1*d18o_orig(0:np2),color=100,thick=3
oplot,d18o_index1,-1*d18o_vect1,color=150,thick=3

openw,1,'/home/ggdjl/genie/genie-utilities/data/input/d18o_30k.dat'
printf,1,reverse(d18o_vect),format='(f9.6)'
close,1

openw,1,'/home/ggdjl/genie/genie-utilities/data/input/d18o_30k_highres.dat'
printf,1,reverse(d18o_vect1),format='(f9.6)'
close,1


stop




end
