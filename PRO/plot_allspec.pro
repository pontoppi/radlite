@analyze.pro
@natconst.pro

thick=3.0      ; Postscript line thickness
size=1.4       ; Postscript char size

dpc = 140.     ; Distance to the object in parsec
;
@problem_params.pro
;
openr,1,'spectrumincls.dat'
ni=0
readf,1,ni
incl=dblarr(ni)
readf,1,incl
close,1
;
openr,1,'spectrum_1.dat'
nf=0
readf,1,nf
data=dblarr(2,nf)
readf,1,data
close,1
;
freq = transpose(data[0,*])
spec = dblarr(nf,ni)
;
for i=0,ni-1 do begin
   s = read_spectrum(file='spectrum_'+$
         strcompress(string(i+1),/remove_all)+'.dat')
   spec[*,i] = s.spectrum[1:*]
endfor
;
; Set y-range
;
case istar of 
   0: yrange=[1d-12,1d-7]
   1: yrange=[1d-13,1d-8]
   2: yrange=[1d-14,1d-9]
endcase
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
device,file='spectrum_all.ps',/color,bits_per_pixel=24
plot_spectrum,s,dpc=dpc,/nodata,xrange=[0.1,1000],yrange=yrange
for i=0,ni-1 do begin
   s.spectrum=spec[*,i]
   plot_spectrum,s,dpc=dpc,/oplot
endfor
device,/close
set_plot,'x'
!p.thick=1
!x.thick=1
!y.thick=1
!z.thick=1
!p.charthick=1
!p.charsize=1


end

