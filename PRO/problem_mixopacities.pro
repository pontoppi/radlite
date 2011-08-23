pro mixopacities,files,fileo,abun
;
nfiles0 = n_elements(files)
nfiles  = nfiles0
for ifile=nfiles0-1,0,-1 do begin
   if files[ifile] eq '' then nfiles=ifiles
endfor
if nfiles lt 2 then return
if n_elements(abun) lt nfiles then stop
ab     = abun
dum    = total(ab)
ab     = ab / dum
;
; Read the first file (which sets the lambda grid)
;
openr,1,files[0]
iformat=0
nf=0
readf,1,iformat
case iformat of
    1: idum=2
    2: idum=3
    3: idum=4
    else: stop
endcase
readf,1,nf
data1=dblarr(idum,nf)
readf,1,data1
close,1
;
; Set the wavelength grid and the main kappa file
;
lambda     = transpose(data1[0,*])
kappa_abs  = ab[0] * transpose(data1[1,*])
if iformat gt 1 then begin
    kappa_sca  = ab[0] * transpose(data1[2,*])
endif else begin
    kappa_sca  = dblarr(nf)
endelse
;
; Read the other files
;
for ifile=1,nfiles-1 do begin
    ;;
    ;; Read opacity file
    ;; 
    openr,1,files[ifile]
    iformat=0
    nf2=0
    readf,1,iformat
    case iformat of
        1: idum=2
        2: idum=3
        3: idum=4
        else: stop
    endcase
    readf,1,nf2
    data2=dblarr(idum,nf2)
    readf,1,data2
    close,1
    ;;
    ;; Interpolate onto the first
    ;;
    kappa_abs = kappa_abs + ab[ifile] * exp(interpol(alog(data2[1,*]),$
                 alog(data2[0,*]),alog(data1[0,*])))
    if iformat gt 1 then begin
        kappa_sca  = kappa_sca + ab[ifile] * exp(interpol(alog(data2[2,*]),$
                 alog(data2[0,*]),alog(data1[0,*])))
    endif
    ;;
endfor
;
; Write the opacity mixture
;
openw,1,fileo
printf,1,3
printf,1,nf
printf,1,' '
for i=0,nf-1 do printf,1,lambda[i],kappa_abs[i],$
     kappa_sca[i],0.d0,format='(4(E13.6,1X))'
close,1
;
end


function readkappa,file
openr,1,file
iformat=0
readf,1,iformat
case iformat of
    1: begin
        nf=0
        readf,1,nf
        data = dblarr(2,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = dblarr(nf)
        g_sca     = dblarr(nf)
    end
    2: begin
        nf=0
        readf,1,nf
        data = dblarr(3,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = transpose(data[2,*])
        g_sca     = dblarr(nf)
    end
    3: begin
        nf=0
        readf,1,nf
        data = dblarr(4,nf)
        readf,1,data
        lambda    = transpose(data[0,*])
        kappa_abs = transpose(data[1,*])
        kappa_sca = transpose(data[2,*])
        g_sca     = transpose(data[3,*])
    end
endcase
close,1	
freq = 2.9979d14/lambda
return,{freq:freq,lambda:lambda,kappa_abs:kappa_abs,kappa_sca:kappa_sca,$
           g_sca:g_sca}
end
