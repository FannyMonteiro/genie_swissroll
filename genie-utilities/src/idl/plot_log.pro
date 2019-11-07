pro plot_log

close,1
nt=300
data=fltarr(nt)
st=''

nf=3
njunk=[3,0,4]
junk3=strarr(3)
junk4=strarr(4)
ymin=[200,0.31,180e-6]
ymax=[300,0.36,280e-6]


names=['orog','lice','co2']

for xx=0,nf-1 do begin
openr,1,'/home/ggdjl/genie_ents/genie-cgoldstein/log_30k_'+names(xx)
if (njunk(xx) eq 3) then begin
readf,1,junk3
endif
if (njunk(xx) eq 4) then begin
readf,1,junk4
endif
for x=0,nt-1 do begin
readf,1,st
data(x)=strtrim(strmid(st,55),2)
endfor
close,1
window,xx
plot,data,ystyle=1,yrange=[ymin(xx),ymax(xx)],title=names(xx)

endfor

stop



end
