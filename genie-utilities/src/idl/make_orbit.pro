
pro make_orbit

close,1

ny=31
nyoff=9
nexp=1
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
expname=strarr(nexp)
inexpname=strarr(nexp,4)
inexpname(0,*)='_old'
;inexpname(1,*)='_new'
expname=strarr(nexp,4)
;expname(0,*)='_old'
;expname(1,*)='_new'
expname(0,*)=''

junk=''

ecc=fltarr(ny)
sob=fltarr(ny)
gam=fltarr(ny)
tau=fltarr(ny)

for z=0,nexp-1 do begin

filename='/home/ggdjl/orbit/orbit'+inexpname(z)+'.dat'
print,filename
openr,1,filename
readf,1,junk
print,junk
for x=0,ny-1 do begin
readf,1,f,f,f,f,f,f,f,f,f,f,d_gam,d_ecc,d_tau,d_sob
ecc(x)=d_ecc
gam(x)=d_gam
tau(x)=d_tau
sob(x)=d_sob
endfor
close,1

;filename='/home/ggdjl/genie/genie-utilities/data/output/orbit_30k_reverse'+strtrim(expname(z))+'.dat'
;print,filename
;openw,1,filename
;for x=0,ny-1 do begin
;  printf,1,ecc(x),format='(f9.5)'
;endfor
;for x=0,ny-1 do begin
;  printf,1,sob(x),format='(f9.5)'
;endfor
;for x=0,ny-1 do begin
;  printf,1,gam(x),format='(f9.5)'
;endfor
;for x=0,ny-1 do begin
;  printf,1,tau(x),format='(f9.5)'
;endfor
;close,1


for s=0,nsnap-1 do begin
filename='/home/ggdjl/genie/genie-utilities/data/output/orbit_'+snapname(s)+strtrim(expname(z))+'.dat'
print,filename
openw,1,filename
  printf,1,ecc(snapnum(s)),format='(f9.5)'
  printf,1,sob(snapnum(s)),format='(f9.5)'
  printf,1,gam(snapnum(s)),format='(f9.5)'
  printf,1,tau(snapnum(s)),format='(f9.5)'
close,1
endfor

ecc=reverse(ecc)
sob=reverse(sob)
gam=reverse(gam)
tau=reverse(tau)

filename='/home/ggdjl/genie/genie-utilities/data/output/orbit_30k'+strtrim(expname(z))+'.dat'
print,filename
openw,1,filename
for x=0,ny-1 do begin
  printf,1,ecc(x),format='(f9.5)'
endfor
for x=0,ny-1 do begin
  printf,1,sob(x),format='(f9.5)'
endfor
for x=0,ny-1 do begin
  printf,1,gam(x),format='(f9.5)'
endfor
for x=0,ny-1 do begin
  printf,1,tau(x),format='(f9.5)'
endfor
close,1

filename='/home/ggdjl/genie/genie-utilities/data/output/orbit_21k'+strtrim(expname(z))+'.dat'
print,filename
openw,1,filename
for x=nyoff,ny-1 do begin
  printf,1,ecc(x),format='(f9.5)'
endfor
for x=nyoff,ny-1 do begin
  printf,1,sob(x),format='(f9.5)'
endfor
for x=nyoff,ny-1 do begin
  printf,1,gam(x),format='(f9.5)'
endfor
for x=nyoff,ny-1 do begin
  printf,1,tau(x),format='(f9.5)'
endfor
close,1


endfor

stop

end
