@readobs.pro
@analyze.pro
@readopac.pro

thick=3.0      ; Postscript line thickness
size=1.4       ; Postscript char size
;
openr,1,'title.txt'
title = ''
readf,1,title,format='(A60)'
close,1
title=strcompress(title)
;
; Read the total flux SED of the model
;
s   = read_spectrum()
o   = readopac(nr=2)
;
; Now plot the SED
; 
!p.thick=thick
!x.thick=thick
!y.thick=thick
!z.thick=thick
!p.charthick=thick
!p.charsize=size
set_plot,'ps'
device,file='spectrum_nfn_8-100.ps',/color,bits_per_pixel=24
plot_spectrum,s,dpc=100,xr=[8,100],/xs,title=title,yr=[2d-9,2d-7],/ys
device,/close
set_plot,'x'
;
set_plot,'ps'
device,file='spectrum_jy_8-100.ps',/color,bits_per_pixel=24
plot_spectrum,s,dpc=100,xr=[8,100],/jy,/xs,title=title,yr=[20,200],/ys
device,/close
set_plot,'x'
;
set_plot,'ps'
device,file='spectrum_nfn_1-1000.ps',/color,bits_per_pixel=24
plot_spectrum,s,dpc=100,xr=[1,1000],title=title,yr=[1d-11,1d-7],/ys
device,/close
set_plot,'x'
;
set_plot,'ps'
device,file='spectrum_jy_1-1000.ps',/color,bits_per_pixel=24
plot_spectrum,s,dpc=100,xr=[1,1000],/jy,title=title,yr=[2,200],/ys
device,/close
set_plot,'x'
;
set_plot,'ps'
device,file='spectrum_nfn_0.1-1000.ps',/color,bits_per_pixel=24
plot_spectrum,s,dpc=100,xr=[0.1,1000],title=title,yr=[1d-11,1d-6],/ys
device,/close
set_plot,'x'
;
set_plot,'ps'
device,file='spectrum_jy_0.1-1000.ps',/color,bits_per_pixel=24
plot_spectrum,s,dpc=100,xr=[0.1,1000],/jy,title=title,yr=[1,1000],/ys
device,/close
set_plot,'x'
;
set_plot,'ps'
device,file='opacity_8-100.ps',/color,bits_per_pixel=24
plotopac,o,xr=[8,100],/xs,title=title
device,/close
set_plot,'x'
;
set_plot,'ps'
device,file='opacity_8-300.ps',/color,bits_per_pixel=24
plotopac,o,xr=[8,300],/xs,title=title
device,/close
set_plot,'x'
;
!p.thick=1
!x.thick=1
!y.thick=1
!z.thick=1
!p.charthick=1
!p.charsize=1


end

