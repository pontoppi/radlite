;------------------------------------------------------;
; READ THE DUST OPACITY FILES                          ;
;------------------------------------------------------;
function readopac,nr=nr,freqfile=freqfile,opacfile=opacfile
close,1
IF NOT KEYWORD_SET(opacfile) THEN BEGIN
    if not keyword_set(nr) then nr=1
    filename = 'dustopac_'+strcompress(string(nr),/remove_all)+'.inp'
ENDIF ELSE BEGIN
    filename=opacfile
ENDELSE
;print,"Reading ",filename
openr,1,filename
nf=0
ns=0
readf,1,nf,ns
if nf ge 1 then begin
    cabs = fltarr(nf,ns)
    csca = fltarr(nf,ns)
    dum  = 0.e0
    for k=1,nf do begin
        for is=1,ns do begin
            readf,1,dum
            cabs(k-1,is-1) = dum
        endfor
    endfor
    for k=1,nf do begin
        for is=1,ns do begin
            readf,1,dum
            csca(k-1,is-1) = dum
        endfor
    endfor
    nrtrange=1
    trange=[0.d0,0.d0]
endif else begin
    nf=ns
    ismooth=0
    nrtrange=0
    readf,1,ns,ismooth,nrtrange
    if(ismooth ne 0) then stop   ; No smoothing yet allowed.
    cabs = fltarr(nf,ns,nrtrange)
    csca = fltarr(nf,ns,nrtrange)
    dum  = 0.e0
    trange=fltarr(nrtrange+1)
    for ir=1,nrtrange do begin
        readf,1,a,b
        trange(ir)=b
        for k=1,nf do begin
            for is=1,ns do begin
                readf,1,dum
                cabs(k-1,is-1,ir-1) = dum
            endfor
        endfor
        for k=1,nf do begin
            for is=1,ns do begin
                readf,1,dum
                csca(k-1,is-1,ir-1) = dum
            endfor
        endfor
    endfor
endelse
close,1
IF KEYWORD_SET(freqfile) THEN file=freqfile ELSE file='frequency.inp'
dd=findfile(file,count=count)
if(count le 0) then begin
    print,"Could not find frequency.inp. Taking frequency.dat"
    file='frequency.dat'
    dd=findfile(file,count=count)
    if(count le 0) then begin
       print,"Could not find frequency.dat either"
       stop
    endif
endif
openr,1,file
nnf=0
readf,1,nnf
if nnf ne nf then begin
   print,"ERROR: frequency file has different nr of points as dustopac file"
   stop
endif
freq = fltarr(nf)
wave = fltarr(nf)
for k=1,nf do begin
    readf,1,dum
    freq(k-1) = dum
    wave(k-1) = 2.9979d14 / dum
endfor
close,1
return,{nf:nf,ns:ns,freq:freq,wave:wave,cabs:cabs,csca:csca,nrt:nrtrange,trange:trange}
end



;------------------------------------------------------;
; PLOT THE DUST OPACITIES                              ;
;------------------------------------------------------;
pro plotopac,a,abs=abs,scat=scat,is=is,oplot=oplot,$
             ylin=ylin,irt=irt,mult=mult,xlin=xlin,$
             _extra=_extra
ipl=0
ylog=1
xlog=1
if not keyword_set(irt) then irt=0
if irt ge a.nrt then begin
   print,'irt too large'
   return
endif
xtitle='!4k!X [!4l!Xm]'
ytitle='!4j!X!D!4k!X!N [cm!U2!Ng!U-1!N]'
if keyword_set(xlin) then xlog=0
if keyword_set(ylin) then ylog=0
if not keyword_set(oplot) then oplot=0
if not keyword_set(mult) then mult=1.d0
if keyword_set(abs) then begin
    if n_elements(is) gt 0 then begin
        if is ge 0 then begin
	    if oplot eq 0 then begin
                plot,a.wave,mult*a.cabs(*,is,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,mult*a.cabs(*,is,irt),_extra=_extra
            endelse
        endif else begin
            dum = mult*a.cabs(*,0,irt)
            for is=1,a.ns-1 do begin
                dum = dum + mult*a.cabs(*,is,irt)
            endfor
	    if oplot eq 0 then begin
               plot,a.wave,dum,ylog=ylog,xlog=xlog,_extra=_extra
	    endif else begin
               oplot,a.wave,dum,_extra=_extra
            endelse
        endelse
    endif else begin
        if oplot eq 0 then begin
            plot,a.wave,mult*a.cabs(*,0,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
        endif else begin
            oplot,a.wave,mult*a.cabs(*,0,irt),_extra=_extra
	endelse
        for is=1,a.ns-1 do begin
            oplot,a.wave,mult*a.cabs(*,is,irt)
        endfor
    endelse
    ipl = 1
endif
if keyword_set(scat) then begin
    if n_elements(is) gt 0 then begin
        if is ge 0 then begin
            if oplot eq 0 then begin
                plot,a.wave,mult*a.csca(*,is,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,mult*a.csca(*,is,irt),_extra=_extra
	    endelse
        endif else begin
            dum = mult*a.csca(*,0,irt)
            for is=1,a.ns-1 do begin
                dum = dum + mult*a.csca(*,is,irt)
            endfor
            if oplot eq 0 then begin
                plot,a.wave,dum,ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,dum,_extra=_extra
	    endelse
        endelse
    endif else begin
        if oplot eq 0 then begin
            plot,a.wave,mult*a.csca(*,0,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
        endif else begin
            oplot,a.wave,mult*a.csca(*,0,irt),_extra=_extra
        endelse
        for is=1,a.ns-1 do begin
            oplot,a.wave,mult*a.csca(*,is,irt)
        endfor
    endelse
    ipl = 1
endif
if ipl eq 0 then begin
    dum = mult*a.csca(*,0,irt) + mult*a.cabs(*,0,irt)
    for is=1,a.ns-1 do begin
        dum = dum + mult*a.csca(*,is,irt) + mult*a.cabs(*,is,irt)
    endfor
    if oplot eq 0 then begin
       plot,a.wave,dum,ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
    endif else begin
       oplot,a.wave,dum,_extra=_extra
    endelse
endif
end

;-----------------------------------------------------------------
;                 FIND A FREQUENCY IN A TABLE
;-----------------------------------------------------------------
function find_freq,freq,nu,eps=eps
nf=n_elements(freq)
idx=-100
for inu=0,nf-2 do begin
   if (nu ge freq(inu) and nu le freq(inu+1)) or $
     (nu ge freq(inu+1) and nu le freq(inu)) then begin
      idx = inu
      eps = (nu-freq[inu]) / (freq[inu+1]-freq[inu])
   endif
endfor
if idx lt 0 then begin
   print,"ERROR: nu out of range"
   stop
endif
return,idx
end

;-----------------------------------------------------------------
;                 FIND OPACITY AT CERTAIN WAVELENGTH
;-----------------------------------------------------------------
function findopac,o,lambda=lambda,freq=freq,abs=abs,sca=sca
if n_elements(abs) eq 0 and n_elements(sca) eq 0 then begin
   abs=1
   sca=1
endif else begin
   if n_elements(abs) eq 0 then abs=0
   if n_elements(sca) eq 0 then sca=0
endelse
if keyword_set(lambda) then freq=2.9979e14/lambda
if not keyword_set(freq) then stop
if freq ge max(o.freq) then return,0.d0
if freq le min(o.freq) then return,0.d0
idx=find_freq(o.freq,freq)
eps=(freq-o.freq[idx])/(o.freq[idx+1]-o.freq[idx])
kappa=0.d0
if abs gt 0 then kappa = kappa + $
   (1.d0-eps)*o.cabs(idx)+eps*o.cabs(idx+1)
if sca gt 0 then kappa = kappa + $
   (1.d0-eps)*o.csca(idx)+eps*o.csca(idx+1)
return,kappa
end

