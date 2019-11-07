
pro plot_orbit_new

loadct,39
close,1

nyr=31
nexp=22
ny=36
nt=100
minplot=-70
maxplot=70
ncontours=28
nseas=12
ny_um=73
nt_um=12

rad=fltarr(ny,nt,nyr)
rad_um=fltarr(ny_um,nt_um,nexp)
xx=intarr(nyr)

seas=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
expname=strarr(nexp)

; FOR REVERSE:
expname=['xamka','xamkv','xamku','xamkt','xamks','xamkr','xamkq','xamkp','xamko','xamkn','xamkm','xamkl','xamkk','xamkj','xamki','xamkh','xamkg','xamkf','xamke','xamkd','xamkc','xamkb']

for x=0,nyr-1 do begin

xx(x)=(nyr-x-1)*100

noughts,xx(x),str,7
filename='/home/ggdjl/genie_kitty/genie-cgoldstein/solfor_'+str+'.dat'
print,filename
openr,1,filename
for j=0,ny-1 do begin
for t=0,nt-1 do begin
readf,1,raddum
rad(j,t,x)=raddum
endfor
endfor

if (x le nexp-1) then begin
for m=0,nt_um-1 do begin
id=ncdf_open('/home/swsvalde/ummodel/data/'+expname(x)+'/'+expname(x)+'a.pdcl'+seas(m)+'.nc')
ncdf_varget,id,'field200',dummy
ncdf_close,id
for j=0,ny_um-1 do begin
rad_um(ny_um-1-j,m,x)=mean(dummy(*,j))
endfor
endfor
endif

close,1





endfor

set_plot,'ps'
!P.FONT=0
device,filename='orbit_plot.eps',/color,/encapsulated
plot,30-findgen(31),rad(1,0,*)-rad(1,0,0),ystyle=1,yrange=[-60,60],xtickname=['30','20','10','0'],xtitle='kyrBP',ytitle='Wm!E-2!N',thick=5,color=10
oplot,30-findgen(31),rad(34,0,*)-rad(34,0,0),thick=5,color=50
oplot,30-findgen(31),rad(1,50,*)-rad(1,50,0),thick=5,color=170
oplot,30-findgen(31),rad(34,50,*)-rad(34,50,0),thick=5,color=230
plots,[0.5,0.6],[0.15,0.15],/normal,thick=5,color=10
plots,[0.5,0.6],[0.18,0.18],/normal,thick=5,color=50
plots,[0.5,0.6],[0.21,0.21],/normal,thick=5,color=170
plots,[0.5,0.6],[0.24,0.24],/normal,thick=5,color=230
xyouts,0.62,0.14,'65oS, Jan',/normal
xyouts,0.62,0.17,'65oN, Jan',/normal
xyouts,0.62,0.20,'65oS, Jul',/normal
xyouts,0.62,0.23,'65oN, Jul',/normal

device,/close



stop
end
