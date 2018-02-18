;========================================================================
;              HERE ALL THE OUTPUT ROUTINES ARE LISTED
;
; These are all the routines for producing the input files for RADMC
; and other programs related to RADMC.
;
; NOTE: This is compatible with the new fast RADMC with diffusion,
;       so that RADICAL-VET is no longer necessary!
;========================================================================

;------------------------------------------------------------------------
;                  WRITE THE GRID AND THE DENSITY
;------------------------------------------------------------------------
pro write_radmc_density,r,theta,rhodust
nr    = n_elements(r)
nt    = n_elements(theta)
nspec = n_elements(rhodust[0,0,*])
openw,1,'radius.inp'
printf,1,nr
printf,1,' '
for i=1,nr do printf,1,r[i-1]
close,1
openw,1,'theta.inp'
printf,1,nt,1
printf,1,' '
for i=1,nt do printf,1,theta[i-1]
close,1
openw,1,'dustdens.inp'
printf,1,nspec,nr,nt,1
printf,1,' '
for ispec=0,nspec-1 do begin
    for ir=1,nr do begin
        for it=1,nt do begin
            printf,1,rhodust[ir-1,it-1,ispec]+1d-90
        endfor
    endfor
    printf,1,' '
endfor
close,1
end

;------------------------------------------------------------------------
;             WRITE THE STELLAR PARAMETERS AND STELLAR SPECTRUM
;
; Writes the input files of RADMC for the stellar parameters
;
; ARGUMENTS:
;   mstar         Mass of the star
;   rstar         Radius of the star
;   tstar         Temperature of the star
;
; KEYWORDS:
;   kurucz        If set, then use Kurucz model (see the files
;                   'interpolate_kurucz' and problem_kurucz.pro)
;   kurdir        Location of the kurucz database
; 
; RETURNS:
;   inupk         Index of the peak of the spectrum
;   lampk         Wavelength (micron) of the peak of the spectrum
;
;------------------------------------------------------------------------
pro write_star,mstar,rstar,tstar,kurucz=kurucz,kurdir=kurdir,$
               inupk=inupk,lampk=lampk
@problem_natconst.pro
;
; Read this (or at least the first) opacity
;
o      = readopac(nr='1')
nf     = n_elements(o.freq)
;
; Create the stellar spectrum (choose between Kurucz or BB)
;
if keyword_set(kurucz) then begin
   lstar = LS * (rstar/RS)^2 * (tstar/TS)^4
   kurucz,tstar,lstar,mstar,/worig,kurdir=kurdir
endif else begin
   openw,1,'starspectrum.inp'
   printf,1,nf
   printf,1,' '
   for inu=0,nf-1 do begin
      printf,1,o.freq[inu],3.14159265359 * (rstar^2/pc^2) $
             * bplanck(o.freq[inu],tstar)
   endfor
   close,1
endelse
;
; Make starinfo.inp
;
openw,1,'starinfo.inp'
printf,1,'1'
printf,1,rstar
printf,1,mstar
printf,1,tstar
close,1
;
; Determine peak of spectrum 
;
openr,1,'starspectrum.inp'
nnf=0
readf,1,nnf
if nnf ne nf then stop
data=dblarr(2,nf)
readf,1,data
close,1
nu     = transpose(data[0,*])
nulnu  = nu * transpose(data[1,*]) * 4 * !dpi * pc^2
;
; Determine the peak of the stellar spectrum
;
inupk  = find_peak_starspec(nu,nulnu)
lampk  = 1d4*cc/nu[inupk]
;
end


;-------------------------------------------------------------------------
;                  WRITE INPUT FILE FOR CHOPDENS
;
; The chopdens program can make sure that the vertical optical depth of
; the disk never exceeds some value. This is done by removing matter from
; the midplane regions which are anyway too optically thick to be observed.
; Note, however, that this is a trick, nothing more. Reasonable values
; are for instance tauchop=1d4,lmbchop=0.55,idxchop=1.d0.
;-------------------------------------------------------------------------
pro write_chopdens,tauchop,lmbchop,idxchop
openw,1,'chopdens.inp'
;
; OLD STYLE:
;
; printf,1,1
; printf,1,tauchop
; printf,1,lmbchop
; printf,1,idxchop
;
; NEW STYLE:
;
printf,1,'taumax  = ',tauchop
printf,1,'lambda  = ',lmbchop
printf,1,'smooth  = ',idxchop
close,1
end



;-------------------------------------------------------------------------
;                       WRITE THE RADMC.INP FILE
;
; This file is not crucial. The system also runs without it. But if you
; want to file-tune the RADMC run, then this file is required. It is 
; written by default by the main problem_models.pro routines.
;-------------------------------------------------------------------------
pro write_radmc,nphot=nphot,iseed=iseed,ifast=ifast,enthres=enthres,$
                cntdump=cntdump,iquant=iquant,ifinstar=ifinstar,$
                npdiff=npdiff,nvstr=nvstr,ivstrt=ivstrt,$
                vserrtol=vserrtol
;
; Include diffusion per default
;
if n_elements(npdiff) eq 0 then npdiff=30
;
; Check something
;
if keyword_set(nvstr) and npdiff eq 0 then begin
    print,'ERROR IN SETUP: If you do vertical structure iteration'
    print,'    then you must also put npdiff gt 0, i.e. you must'
    print,'    use the diffusion module of RADMC to make sure that'
    print,'    the midplane temperatures are smooth.'
    stop
endif
;
; Default
;
if not keyword_set(nvstr) then nvstr=0
if not keyword_set(ivstrt) then ivstrt=1
if not keyword_set(vserrtol) then vserrtol=0.d0   ; Meaning iterate nvstr always
;
; Check if the iquant is specified
;
if n_elements(iquant) eq 0 then iquant=0
;
; Check if nphot is specified
;
if not keyword_set(nphot) then nphot=100000
;
; Check if iseed is specified
;
if not keyword_set(iseed) then iseed=-17933201
;
; One more thing
;
if n_elements(ifinstar) eq 0 then ifinstar=0
;
; Check if ifast is specified. The fast method can be 20x faster, but
; is perhaps not 100% perfectly reliable. Nevertheless it is recommended
; to enable fast by putting ifast=1, since otherwise some disk models will
; never be doable. Fast is 0 by default. 
; (Note: this is still the old fast method; not the one I devised
;        on suggestion of Michiel Min. It is an approximative method
;        by which the thermal calculations are not always done.)
;
if not keyword_set(ifast) then ifast=0
if not keyword_set(enthres) then enthres=1.d-2
;
; By default do not dump a safety dump, but if you specify 
; cntdump, then the safety dump for RADMC will be made every cntdump
; photons
;
if not keyword_set(cntdump) then cntdump=100000000   ; (effectively no dump!)
;
; First make sure to destroy any radmc.inp that is there
;
spawn,'rm -f radmc.inp >& /dev/null'
;
; Now write radmc.inp
;
; OLD STYLE
;
; openw,1,'radmc.inp'
; printf,1,3            ; Iformat of file (2)
; printf,1,nphot        ; Nr of photons
; printf,1,1            ; (For compatibility with other MC code)
; printf,1,iseed        ; Seed for MC random generator
; printf,1,2            ; Method is always '2'
; printf,1,ifast        ; Fast method (Note: this is still the old fast method!)
; printf,1,enthres      ; If fast method, then what is the energy threshold?
; printf,1,cntdump      ; Safety dump every cntdump
; printf,1,0            ; No restart, just a full start
; printf,1,1            ; Decouple temperatures
; printf,1,iquant       ; Include quantum-heating of small grains?
; printf,1,ifinstar     ; Finite-size star?  0= no.
; printf,1,npdiff    ; The diffusion method near the midplane
; printf,1,1d-10        ; Convergence criterion for diffusion method
; printf,1,nvstr        ; Nr of iterations of vertical structure
; printf,1,vserrtol     ; Vertical structure error tolerance
; printf,1,ivstrt       ; Which dust species represents T_gas?
; printf,1,3000         ; Nr of temperatures for emission database
; printf,1,1d-2         ; Minimum temp of database
; printf,1,1d6          ; Maximum temp of database
; close,1
;
; NEW STYLE
;
openw,1,'radmc.inp'
printf,1,'nphot       = ', nphot        
printf,1,'iseed       = ', iseed        
printf,1,'imethod     = ', 2            
printf,1,'ifast       = ', ifast        
printf,1,'enthres     = ', enthres      
printf,1,'cntdump     = ', cntdump      
printf,1,'irestart    = ', 0            
printf,1,'itempdecoup = ', 1            
printf,1,'iquantum    = ', iquant       
printf,1,'istarsurf   = ', ifinstar     
printf,1,'nphotdiff   = ', npdiff       
printf,1,'errtol      = ', 1d-10        
printf,1,'nvstr       = ', nvstr        
printf,1,'vserrtol    = ', vserrtol     
printf,1,'ivstrt      = ', ivstrt       
printf,1,'ntemp       = ', 3000         
printf,1,'temp0       = ', 1d-2         
printf,1,'temp1       = ', 1d6          
close,1
;
; Done...
;
end



;-------------------------------------------------------------------------
;               WRITE THE INPUT FILE FOR THE RAY TRACER
;
; In the past the ray tracer (for spectra and images) was the full 
; radiative transfer code RADICAL. Now the ray-trace module from RADICAL
; is extracted and put into an independent little program called RAYTRACE.
; It can (and should) do exactly the same as RADICAL, but it has some
; new possibilities such as a better gridding of the pixels of the 
; circular images near the inner rim.
;
; In principle RAYTRACE does not need a separate input file, but it can
; have it. Here is a subroutine to write one.
;-------------------------------------------------------------------------
pro write_raytrace,nextra=nextra,nrref=nrref
if n_elements(nextra) eq 0 then nextra=20
if nextra lt 0 then nextra=-nextra
if n_elements(nrref) eq 0 then nrref=10
; BUGFIX 04.04.07: Was rayrace.inp...
;
; OLD STYLE
;
; openw,1,'raytrace.inp'
; printf,1,2        ; Format number
; printf,1,32       ; Nr of pixels in circle
; printf,1,-nextra  ; (abs) Nr of rays inward of inner rim ('-' means adaptive)
; if nrref gt 0 then printf,1,1 else printf,1,0
; if nrref gt 0 then printf,1,nrref    ; Nr of pix refinement INward of inner edge (NEW: Feb 2007)
; printf,1,1        ; Nr of rays per radial grid point of the model
; printf,1,45.      ; Standard inclination (if nothing is given)
; close,1
;
; NEW STYLE
;
openw,1,'raytrace.inp'
printf,1,'nrphiinf    = ',32       
printf,1,'nrrayextra  = ',-nextra  
if nrref gt 0 then printf,1,'imethod     = ',1 else printf,1,'imethod     = ',0
if nrref gt 0 then printf,1,'nrref       = ',nrref
printf,1,'dbdr        = ',1
printf,1,'inclination = ',45.
close,1
end



;-------------------------------------------------------------------------
;     WRITE THE INPUT FILE FOR THE VERTICAL STRUCTURE MODULE OF RADMC
;
; This file is not strictly necessary. If not provided then all the 
; dust species will participate in the vertical structure calculation.
;-------------------------------------------------------------------------
pro write_vstruct,dostruct
nspec=n_elements(dostruct)
openw,1,'vstruct.inp'
printf,1,1        ; Format number
printf,1,nspec    ; Nr of pixels in circle
for ispec=0,nspec-1 do begin
   printf,1,dostruct[ispec]
endfor
close,1
end
