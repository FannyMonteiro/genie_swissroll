
pro plot_orbit

device,decomposed=0
loadct,39
close,1

nyr=30
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
;expname=['xamka','xamkv','xamku','xamkt','xamks','xamkr','xamkq','xamkp','xamko','xamkn','xamkm','xamkl','xamkk','xamkj','xamki','xamkh','xamkg','xamkf','xamke','xamkd','xamkc','xamkb']
; FOR NORMAL:
expname=['xamkv','xamku','xamkt','xamks','xamkr','xamkq','xamkp','xamko','xamkn','xamkm','xamkl','xamkk','xamkj','xamki','xamkh','xamkg','xamkf','xamke','xamkd','xamkc','xamkb']


for x=0,nyr-1 do begin

; FOR REVERSE:
;xx(x)=x*1000
; FOR NORMAL:
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


contours=minplot+(maxplot-minplot)*findgen(ncontours)/(ncontours-1.0)

window,0,xsize=600,ysize=400
contour,transpose(rad(*,*,x)-rad(*,*,0)),/cell_fill,levels=contours,title=strtrim(xx(x),2)+'-'+strtrim(xx(0),2)
contour,/over,transpose(rad(*,*,x)-rad(*,*,0)),levels=contours,/follow

if (x le nexp-1) then begin
window,1,xsize=600,ysize=400
contour,transpose(rad_um(*,*,x)-rad_um(*,*,0)),/cell_fill,levels=contours,title=strtrim(x,2)+' '+expname(x)+'-'+expname(0)
contour,/over,transpose(rad_um(*,*,x)-rad_um(*,*,0)),levels=contours,/follow
endif
wait,5

endfor

stop
end
