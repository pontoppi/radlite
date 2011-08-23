
;-----------------------------------------------------------------
;              READ THE RECTANGULAR TELESCOPE IMAGE
;-----------------------------------------------------------------
function readimage,ext
nx=0
ny=0
nf=0
sizepix_x=0.d0
sizepix_y=0.d0
if n_elements(ext) eq 0 then begin
   filename='image.dat'
endif else begin
   filename='image_'+strcompress(string(ext),/remove_all)+'.dat'
endelse
str=findfile(filename)
if(str(0) ne filename) then begin
    print,'Sorry, cannot find ',filename
    print,'Presumably radical exited without succes.'
    print,'See above for possible error messages of radical!'
    stop
endif
openr,1,filename
readf,1,nx,ny,nf
readf,1,sizepix_x,sizepix_y
image=fltarr(nx,ny)
tau=fltarr(nx,ny)
readf,1,image
close,1
if n_elements(ext) eq 0 then begin
   filename='imtau.dat'
endif else begin
   filename='imtau_'+strcompress(string(ext),/remove_all)+'.dat'
endelse
str=findfile(filename)
if(str(0) eq filename) then begin
openr,1,filename
idum1=0
idum2=0
idum3=0
readf,1,idum1,idum2,idum3
if(idum1 ne nx) then stop
if(idum2 ne ny) then stop
if(idum3 ne nf) then stop
dum1=0
dum2=0
readf,1,dum1,dum2
readf,1,tau
close,1
endif
return,{nx:nx,ny:ny,nrfr:nf,sizepix_x:sizepix_x,sizepix_y:sizepix_y,$
        image:image,tau:tau}
end


;-----------------------------------------------------------------
;                   PLOT THE RECTANGULAR IMAGE
;-----------------------------------------------------------------
pro plotimage,a,log=log,pc=pc,au=au,contour=contour,nlevels=nlevels,$
        noimage=noimage,pos=pos,thrcol=thrcol,thres=thres,maxlog=maxlog,$
        saturate=saturate,xticklen=xticklen,yticklen=yticklen,$
        xminor=xminor,yminor=yminor,plottau=plottau,lgrange=lgrange,$
        levels=levels,_EXTRA=_extra
common colcodes,green,yellow,red,blue,black,white
common colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr
;
; Lin or log
;
if keyword_set(log) eq 0 then begin
    imagedata = a.image
endif else begin
    imagedata = a.image
    taudata   = a.tau
    ;;
    ;; Log this image, but make sure to have a limited range
    ;;
    imagedata = alog10(imagedata+1.d-200)
    if not keyword_set(lgrange) then begin
       mx              = max(imagedata)
       mn              = min(imagedata)
       if keyword_set(maxlog) gt 0.d0 then begin
          if mx-mn gt maxlog then begin
             mn = mx - maxlog
          endif
       endif
       imagedata = imagedata - mn 
       i               = where(imagedata lt 0.d0)
       if n_elements(i) gt 1 and i[0] ge 0 then imagedata[i] = 0.d0
    endif else begin
       minl=lgrange[0]
       maxl=lgrange[1]
       i               = where(imagedata lt minl)
       if n_elements(i) gt 1 and i[0] ge 0 then imagedata[i] = minl
       i               = where(imagedata gt maxl)
       if n_elements(i) gt 1 and i[0] ge 0 then imagedata[i] = maxl
       imagedata=imagedata-minl
    endelse
    ;;
    ;; Do the same for the tau
    ;;
    taudata   = alog10(taudata+1.d-200)
    if not keyword_set(lgrange) then begin
       mx              = max(taudata)
       mn              = min(taudata)
       if keyword_set(maxlog) then begin
          if mx-mn gt maxlog then begin
             mn = mx - maxlog
          endif
       endif
       taudata   = taudata - mn
       i               = where(taudata lt 0.d0)
       if n_elements(i) gt 1 and i[0] ge 0 then taudata[i] = 0.d0
    endif else begin
       minl=lgrange[0]
       maxl=lgrange[1]
       i               = where(taudata lt minl)
       if n_elements(i) gt 1 and i[0] ge 0 then taudata[i] = minl
       i               = where(taudata gt maxl)
       if n_elements(i) gt 1 and i[0] ge 0 then taudata[i] = maxl
       taudata=taudata-minl
    endelse
endelse
;
; Interpreting the input
;
distscale=1.d0
distunit='cm'
if keyword_set(pc) then begin
    distscale = 3.24073473938d-19
    distunit  = 'pc'
endif
if keyword_set(au) then begin
    distscale = 6.68449197861d-14
    distunit  = 'AU'
endif
if not keyword_set(nlevels) then nlevels=10
if not keyword_set(xticklen) then xticklen=-0.02
if not keyword_set(yticklen) then yticklen=-0.02
if not keyword_set(xminor) then xminor=-1
if not keyword_set(yminor) then yminor=-1
if not keyword_set(mincol) then mincol=0
if not keyword_set(pos) then begin
    if !d.name eq 'X' then begin
        pos=[70,70,400,400]
    endif else begin
        pos=[0.1,0.1,0.8,0.8]
    endelse
endif
if not keyword_set(plottau) then plottau=0
if plottau eq 0 then begin
    image=imagedata
endif else begin
    image=taudata
endelse
immir=rotate(image,2)
;
x=(findgen(a.nx)/(a.nx-1.d0)-0.5d0)*a.sizepix_x*a.nx*distscale
y=(findgen(a.ny)/(a.ny-1.d0)-0.5d0)*a.sizepix_y*a.ny*distscale
if !d.name eq 'X' then begin
    erase
    if not keyword_set(noimage) then begin
        z=congrid(immir,pos(2),pos(3))
        if keyword_set(saturate) then begin
	    dum = saturate*max(z)
	    z   = (z lt dum)*z + (z ge dum)*dum
	endif
        if keyword_set(thres) then begin
	    dum = thres*max(z)
	    z   = (z gt dum)*z - dum
        endif
        if not keyword_set(lgrange) then begin
           z=z/max(z)
           bimg=bytscl(z,min=mincol,top=!d.table_size)
        endif else begin
           z=(!d.table_size-mincol-1)*z/(maxl-minl)+mincol
           bimg=bytarr(pos(2),pos(3))
           bimg[*,*]=z[*,*]
        endelse
        tv,bimg,pos(0),pos(1)
    endif
    plot,x,y,/noerase,position=[pos(0),pos(1),$
      pos(0)+pos(2),pos(1)+pos(3)],/device,/nodata,$
      xminor=xminor,yminor=yminor,xstyle=1,ystyle=1,$
      xticklen=xticklen,yticklen=yticklen,_EXTRA=_extra
    if keyword_set(contour) then begin
      contour,immir,x,y,nlevels=nlevels,position=[pos(0),pos(1),$
      pos(0)+pos(2),pos(1)+pos(3)],/device,$
      xminor=xminor,yminor=yminor,xstyle=1,ystyle=1,$
      xticklen=xticklen,yticklen=yticklen,/noerase,levels=levels,_EXTRA=_extra
    endif
endif else begin
    if not keyword_set(noimage) then begin
        z = immir
        if keyword_set(saturate) then begin
	    dum = saturate*max(z)
	    z   = (z lt dum)*z + (z ge dum)*dum
	endif
        if keyword_set(thres) then begin
	    dum = thres*max(z)
	    z   = (z gt dum)*z - dum
        endif
        if not keyword_set(lgrange) then begin
           z=z/max(z)
           bimg=bytscl(z,min=0.d0,top=!d.table_size)
        endif else begin
           z=(!d.table_size-0.d0-1)*z/(maxl-minl)+mincol
           bimg=bytarr(n_elements(z[*,0]),n_elements(z[0,*]))
           bimg[*,*]=z[*,*]
        endelse
        tv,bimg,pos(0),pos(1),$
             xsize=pos(2),ysize=pos(3),/norm
    endif
    plot,x,y,/noerase,position=[pos(0),pos(1),$
      pos(0)+pos(2),pos(1)+pos(3)],/norm,/nodata,$
      xminor=xminor,yminor=yminor,xstyle=1,ystyle=1,$
      xticklen=xticklen,yticklen=yticklen,_EXTRA=_extra
    if keyword_set(contour) then begin
      contour,immir,x,y,nlevels=nlevels,position=[pos(0),pos(1),$
      pos(0)+pos(2),pos(1)+pos(3)],/norm,$
      xminor=xminor,yminor=yminor,xstyle=1,ystyle=1,$
      xticklen=xticklen,yticklen=yticklen,/noerase,$
      levels=levels,_EXTRA=_extra
    endif
endelse
end
