;====================================================================
;
;        IDL ROUTINES FOR READING/ANALYZING RADICAL OUTPUT
;
;                        general part
;
;                    C.P. Dullemond (2000)
;
;====================================================================

;--------------------------------------------------------------------
;              A FUNCTION USED BY THE ROUTINES BELOW
;--------------------------------------------------------------------
function sigmamass,rhodust,r,t,imirt
nr=n_elements(r)
nt=n_elements(t)
ri=dblarr(nr+1)
ti=dblarr(nt+1)
ri[0]=r[0]
ri[nr]=r[nr-1]
ti[0]=0.
if imirt eq 1 then ti[nt]=!dpi/2 else ti[nt]=!dpi
for ir=1,nr-1 do ri[ir]=sqrt(r[ir]*r[ir-1])
for it=1,nt-1 do ti[it]=0.5*(t[it]+t[it-1])
m=0.d0
sigma=dblarr(nr)
for ir=0,nr-1 do begin
   surf = !dpi * ( ri[ir+1]^2 - ri[ir]^2 )
   for it=0,nt-1 do begin
      vol       = (2.*!dpi/3.)*(ri[ir+1]^3-ri[ir]^3)*abs(cos(ti[it])-cos(ti[it+1]))
      m         = m + vol * rhodust[ir,it]
      sigma[ir] = sigma[ir] + vol * rhodust[ir,it] / surf
   endfor
endfor
if imirt eq 1 then begin
   m=m*2
   sigma=sigma*2
endif
return,{sigma:sigma,m:m}
end

;--------------------------------------------------------------------
;              A FUNCTION USED BY THE ROUTINES BELOW
;--------------------------------------------------------------------
function columns,rho,r,theta
sz=size(rho)
nr=sz(1)
nt=sz(2)
if sz[0] gt 2 then begin
    ns = sz[3]
endif else begin
    ns = 1
endelse
colr=fltarr(nr,nt,ns)
colt=fltarr(nr,nt,ns)
for is=0,ns-1 do begin
    for it=0,nt-1 do begin
        colr[0,it,is] = 0.d0
        for ir=1,nr-1 do begin
            colr[ir,it,is] = colr[ir-1,it,is] + $
            0.5*(rho(ir,it,is)+rho(ir-1,it,is))* $
            (r(ir)-r(ir-1))
        endfor
    endfor
    for ir=0,nr-1 do begin
        colt[ir,0,is] = 0.d0
        for it=1,nt-1 do begin
            colt[ir,it,is] = colt[ir,it-1,is] + $
            0.5*(rho(ir,it)+rho(ir,it-1))* $
            (theta(it)-theta(it-1))*r(ir)
        endfor
    endfor
endfor
return,{colr:colr,colt:colt}
end


;====================================================================
;                         READ FUNCTIONS
;====================================================================


;--------------------------------------------------------------------
;                    READ THE GRID INFORMATION
;--------------------------------------------------------------------
function read_grid,dat=dat
nr=0
ntheta=0
nthsmall=0
nfr=0
imirt=0
close,1
;
; Read the radius data
;
file='radius.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file='radius.dat'
    dd=findfile(file,count=count)
    if(count le 0) then stop
endif
openr,1,file
readf,1,nr
r=dblarr(nr)
readf,1,r
close,1
;
; Read the theta data
;
file='theta.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file='theta.dat'
    dd=findfile(file,count=count)
    if(count le 0) then stop
endif
openr,1,file
readf,1,ntheta,imirt
nthsmall=ntheta
theta=dblarr(ntheta)
readf,1,theta
close,1
if imirt eq 1 then begin
   thbk=dblarr(2*ntheta)
   thbk(0:ntheta-1) = theta(*)
   thbk(ntheta:2*ntheta-1) = !pi - rotate(theta(*),2)
   theta = thbk
   ntheta = ntheta*2
endif
;
; Now read the frequencies
;
file='frequency.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file='frequency.dat'
    dd=findfile(file,count=count)
    if(count le 0) then begin
       print,'Could not find frequency.inp/.dat!!'
       nu=[1.0,1.0]
       nnu=0
    endif
endif
filename=file
str=findfile(filename)
if(str(0) eq filename) then begin
    openr,1,filename
    nnu=1
    readf,1,nnu
    nu=dblarr(nnu+1)
    dum=dblarr(nnu)
    readf,1,dum
    nu(1:*) = dum(*)
    close,1
endif else begin
    nu=[0.d0,0.d0]
    nnu=0
endelse
;
return,{nr:nr,ntheta:ntheta,nthsmall:nthsmall,nnu:nnu,r:r,theta:theta,nu:nu}
end


;--------------------------------------------------------------------
;                             READ DENSITY
;--------------------------------------------------------------------
function read_density,dat=dat,ext=ext
grid = read_grid()
;
filebase='density'
if n_elements(ext) gt 0 then begin
   filebase=filebase+'_'+strcompress(string(ext),/remove_all)
endif
file=filebase+'.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file=filebase+'.dat'
    dd=findfile(file,count=count)
    if(count le 0) then stop
endif
openr,1,file
sizer=0
sizet=0
imirt=0
readf,1,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
rho=dblarr(sizer,sizetbig)
for ir=0,sizer-1 do begin
    for it=0,sizet-1 do begin
        readf,1,a
        rho(ir,it) = a
    endfor
    if(imirt eq 1) then begin
        for it=0,sizet-1 do begin
            rho(ir,sizetbig-it-1) = rho(ir,it)
        endfor
    endif
endfor
;
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
;
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,rho:rho}
end



;--------------------------------------------------------------------
;                         READ TEMPERATURE
;--------------------------------------------------------------------
function read_temperature,dat=dat,ext=ext
grid = read_grid()
;
filebase='temperature'
if n_elements(ext) gt 0 then begin
   filebase=filebase+'_'+strcompress(string(ext),/remove_all)
endif
file=filebase+'.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file=filebase+'.dat'
    dd=findfile(file,count=count)
    if(count le 0) then stop
endif
openr,1,file
sizer=0
sizet=0
imirt=0
readf,1,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
temp=dblarr(sizer,sizetbig)
for ir=0,sizer-1 do begin
    for it=0,sizet-1 do begin
        readf,1,a
        temp(ir,it) = a
    endfor
    if(imirt eq 1) then begin
        for it=0,sizet-1 do begin
            temp(ir,sizetbig-it-1) = temp(ir,it)
        endfor
    endif
endfor
;
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
;
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,t:temp}
end



;--------------------------------------------------------------------
;                        READ MAGNETIC FIELD
;--------------------------------------------------------------------
function read_magfield,dat=dat,ext=ext
grid = read_grid()
;
filebase='magfield'
if n_elements(ext) gt 0 then begin
   filebase=filebase+'_'+strcompress(string(ext),/remove_all)
endif
file=filebase+'.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file=filebase+'.dat'
    dd=findfile(file,count=count)
    if(count le 0) then stop
endif
openr,1,file
sizer=0
sizet=0
imirt=0
readf,1,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
b=dblarr(sizer,sizetbig)
for ir=0,sizer-1 do begin
    for it=0,sizet-1 do begin
        readf,1,a
        b(ir,it) = a
    endfor
    if(imirt eq 1) then begin
        for it=0,sizet-1 do begin
            b(ir,sizetbig-it-1) = b(ir,it)
        endfor
    endif
endfor
;
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
;
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,b:b}
end




;--------------------------------------------------------------------
;                        READ MAGNETIC FIELD
;--------------------------------------------------------------------
function read_iontemp,dat=dat,ext=ext
grid = read_grid()
;
filebase='iontemp'
if n_elements(ext) gt 0 then begin
   filebase=filebase+'_'+strcompress(string(ext),/remove_all)
endif
file=filebase+'.inp'
dd=findfile(file,count=count)
if keyword_set(dat) then count=0
if(count le 0) then begin
    file=filebase+'.dat'
    dd=findfile(file,count=count)
    if(count le 0) then stop
endif
openr,1,file
sizer=0
sizet=0
imirt=0
readf,1,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
ti=dblarr(sizer,sizetbig)
for ir=0,sizer-1 do begin
    for it=0,sizet-1 do begin
        readf,1,a
        ti(ir,it) = a
    endfor
    if(imirt eq 1) then begin
        for it=0,sizet-1 do begin
            ti(ir,sizetbig-it-1) = ti(ir,it)
        endfor
    endif
endfor
;
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
;
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,ti:ti}
end




;-----------------------------------------------------------------
; READ THE ITERATED ELECTRON TEMPERATURE
;-----------------------------------------------------------------
function read_electemp,ext=ext,start=start
if not keyword_set(start) then start=0
grid = read_grid()
print,"Reading electemp info"
;
; Check out how many to read
;
if n_elements(ext) gt 0 then begin
    itertrans=0
    filenamefinal = 'electrontemp_'+strcompress(string(ext),/remove_all)+'.dat'
    final=2
endif else begin
    itertrans=0
    openr,1,'electrontemp.info'
    readf,1,itertrans
    print,itertrans
    if itertrans eq -2 then begin
        filenamefinal = 'electrontemp_final.dat'
        readf,1,itertrans
        itertrans = itertrans+1
        final = 1
    endif else begin
        if itertrans eq -3 then begin
            filenamefinal=' '
            readf,1,filenamefinal
            filenamefinal = strtrim(filenamefinal,2)
            readf,1,itertrans
            itertrans = itertrans+1
            final = 2
        endif else begin
            if start eq 0 then begin
                filenamefinal = 'electrontemp_1.dat'
                ;filenamefinal = 'electrontemp_0.dat'
                final = 0
            endif else begin
                filenamefinal = 'electrontemp_'+strcompress(string(start+1),/remove_all)+'.dat'
                final = 0
            endelse
        endelse
    endelse
endelse
close,1
;
;   First read the dimensions and other information
;
openr,1,filenamefinal
isx=0
isy=0
isfr=0
imirt=0
nlevels=0
readf,1,isx,isy,imirt
if(isx ne grid.nr) then begin
    print,'Dimensions in radius not compatible'
endif
if(isy ne grid.nthsmall) then begin
    print,'Dimensions in theta not compatible'
endif
te  = dblarr(isx,(1+imirt)*isy,itertrans)
close,1
;
;
;   Then do the loop over iterations and read each file
;
for iter=0,itertrans-start-1 do begin
    if((final eq 1) and (iter eq itertrans-start-1)) then begin
        filename = 'electrontemp_final.dat'
    endif else begin
        if((final eq 2) and (iter eq itertrans-start-1)) then begin
            filename = filenamefinal
        endif else begin
            filename = 'electrontemp_'+strcompress(string(iter+start+1),/remove_all)+'.dat'
        endelse
    endelse
    openr,1,filename
    readf,1,idum,idum
    dum=0.d0
    for ir=1,isx do begin
        for it=1,isy do begin
	    readf,1,dum
	    te(ir-1,it-1,iter)= dum
            if(imirt eq 1) then begin
                itt = 2*isy + 1 - it
	        te(ir-1,itt-1,iter)= dum
            endif
        endfor
    endfor
    close,1
endfor
close,1
;
; Then read the ion tempereture
;
filename='iontemp.dat'
str=findfile(filename)
if(str(0) eq filename) then begin
    openr,1,filename
    nt=0
    nr=0
    im=0
    readf,1,nr,nt,im
    ti=dblarr(nr,(1+im)*nt)
    for ir=1,nr do begin
        for it=1,nt do begin
	    readf,1,dum
	    ti(ir-1,it-1)= dum
            if(imirt eq 1) then begin
                itt = 2*nt + 1 - it
	        ti(ir-1,itt-1)= dum
            endif
        endfor
    endfor
    close,1
endif else begin
    ti=0.d0
endelse
;
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
;
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,te:te,ti:ti,niter:itertrans-start}
end



;-----------------------------------------------------------------
;                   READ DUST DENSITIES
;-----------------------------------------------------------------
function read_dustdens,ext=ext,dat=dat,col=col
mass = 1

grid = read_grid()
dum=0

nrspec=0
sizer=0
sizet=0
imirt=0

filebase='dustdens'
if n_elements(ext) gt 0 then begin
   filebase=filebase+'_'+strcompress(string(ext),/remove_all)
endif
if not keyword_set(dat) then begin
   filename=filebase+'.inp' 
endif else begin
   filename=filebase+'.dat' 
endelse
dd=findfile(filename,count=count)
if(count eq 0) then begin
   filename=filebase+'.dat'
   dd=findfile(filename,count=count)
   if(count eq 0) then stop
endif
openr,1,filename
readf,1,nrspec,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
dens=dblarr(sizer,sizetbig,nrspec)
idum=0
for ispec=1,nrspec do begin
    for ir=1,sizer do begin
        for it=1,sizet do begin
            dum = 0.d0
            readf,1,dum
            dens(ir-1,it-1,ispec-1) = dum
            if(imirt eq 1) then begin
                itt = 2*sizet + 1 - it
                dens(ir-1,itt-1,ispec-1) = dum
            endif
        endfor
    endfor
endfor
close,1
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
if keyword_set(mass) then begin
   sigmadust=dblarr(grid.nr,nrspec)
   massdust=dblarr(nrspec)
   for ispec=0,nrspec-1 do begin
      qq=sigmamass(dens[*,*,ispec],grid.r,grid.theta,0)
      sigmadust[*,ispec]=qq.sigma
      massdust[ispec]=qq.m
   endfor
endif else begin
   sigmadust=-1
   massdust=-1
endelse
if keyword_set(col) then begin
   qq=columns(dens,grid.r,grid.theta)
   colr=qq.colr
   colt=qq.colt
endif else begin
   colr=-1
   colt=-1
endelse
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        nspec:nrspec,rr:rr,tt:tt,rho:dens,$
        sigmadust:sigmadust,massdust:massdust,colr:colr,colt:colt}
end


;-----------------------------------------------------------------
;                   READ DUST TEMPERATURES 
;-----------------------------------------------------------------
function read_dusttemp,ext=ext,all=all
grid = read_grid()
niter=0
dum=0
nrsizemax=0
openr,1,'dusttemp.info'
readf,1,niter
readf,1,dum
readf,1,nrsizemax
readf,1,lastsave
close,1

iter=0
nrspec=0
sizer=0
sizet=0
imirt=0

if not keyword_set(all) then begin
   if(niter eq -2) then niter='final'
   if n_elements(ext) gt 0 then niter=ext
   filename='dusttemp_'+strcompress(string(niter),/remove_all)+'.dat'
endif else begin
   filename=['']
   for i=1,lastsave do begin
      filename=[filename,'dusttemp_'+strcompress(string(i),/remove_all)+'.dat']
   endfor
   if niter eq -2 then begin
      filename=[filename,'dusttemp_final.dat']
   endif
   filename=filename(1:*)
endelse

nn=n_elements(filename)

for in=1,nn do begin
   openr,1,filename[in-1]
   readf,1,nrspec,sizer,sizet,imirt
   sizetbig=sizet*(1+imirt)
   if in eq 1 then begin
      temp=dblarr(sizer,sizetbig,nrsizemax,nrspec,nn)
      nrgrainsize=intarr(nrspec)
   endif
   idum=0
   for ispec=1,nrspec do begin
      readf,1,idum
      nrgrainsize(ispec-1)=idum
      for isize=1,nrgrainsize(ispec-1) do begin
         for ir=1,sizer do begin
            for it=1,sizet do begin
               dum = 0.d0
               readf,1,dum
               temp(ir-1,it-1,isize-1,ispec-1,in-1) = dum
               if(imirt eq 1) then begin
                  itt = 2*sizet + 1 - it
                  temp(ir-1,itt-1,isize-1,ispec-1,in-1) = dum
               endif
            endfor
         endfor
      endfor
   endfor
   close,1
endfor
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
;
; If the starinfo file is there, then read this
;
dd=findfile('starinfo.inp',count=count)
if(count gt 0) then begin
   openr,1,'starinfo.inp'
   readf,1,iformat
   readf,1,rstar
   readf,1,mstar
   readf,1,tstar
   close,1
   tthin = sqrt(0.5d0*rstar/rr)*tstar   
endif else begin
   tthin=0.d0
endelse
;
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,nspec:nrspec,nsize:nrgrainsize,$
        temp:temp,tthin:tthin}
end



;--------------------------------------------------------------------
;                     READ THE MEAN INTENSITY
;--------------------------------------------------------------------
function read_meanint,ext=ext,all=all
grid = read_grid()
niter=0
dum=0
;openr,1,'intmean.info'
;readf,1,niter
;readf,1,dum
;readf,1,dum
;readf,1,lastsave
;close,1

iter=0
sizer=0
sizet=0
imirt=0
nfr=0

;if not keyword_set(all) then begin
;   if(niter eq -2) then niter='final'
;   if n_elements(ext) gt 0 then niter=ext
;   filename='intmean_'+strcompress(string(niter),/remove_all)+'.dat'
;endif else begin
;   filename=['']
;   for i=1,lastsave do begin
;      filename=[filename,'intmean_'+strcompress(string(i),/remove_all)+'.dat']
;   endfor
;   if niter eq -2 then begin
;      filename=[filename,'intmean_final.dat']
;   endif
;   filename=filename(1:*)
;endelse

;nn=n_elements(filename)

;for in=1,nn do begin
openr,1,'meanint_radmc.dat';        filename[in-1]
readf,1,nfr,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
;if in eq 1 then begin
meanint=dblarr(sizer,sizetbig,nfr)
;endif
for inu=1,nfr do begin
   for it=1,sizet do begin
      for ir=1,sizer do begin
         dum = 0.d0
         readf,1,dum
         meanint(ir-1,it-1,inu-1) = dum
         if(imirt eq 1) then begin
            itt = 2*sizet + 1 - it
            meanint(ir-1,itt-1,inu-1) = dum
         endif
      endfor
   endfor
endfor
close,1
;endfor

rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,meanint:meanint}
end


;--------------------------------------------------------------------
;                     READ THE FULL INTENSITY
;--------------------------------------------------------------------
function read_intensity
grid = read_grid()
sizer=0
sizet=0
imirt=0
imirp=0
idum=0
sizemu=0
sizephi=0

filename='intensity.dat'
openr,1,filename
readf,1,sizer,sizet,sizemu,sizephi,idum,imirt,imirp
sizetbig=sizet*(1+imirt)
sizephibig=sizephi*(1+imirp)
intensity=dblarr(sizephibig,sizemu,sizer,sizetbig)
for ir=1,sizer do begin
  for it=1,sizet do begin
     for imu=1,sizemu do begin
        dum = dblarr(sizephi)
        readf,1,dum
        if(imirp eq 1) then begin
           intensity(sizephi/2:3*sizephi/2-1,imu-1,ir-1,it-1) = dum(*)
           dum = rotate(dum,2)
           intensity(0:sizephi/2-1,imu-1,ir-1,it-1) = dum(sizephi/2:sizephi-1)
           intensity(3*sizephi/2:2*sizephi-1,imu-1,ir-1,it-1) = dum(0:sizephi/2-1)
        endif else begin
           intensity(*,imu-1,ir-1,it-1) = dum(*)
        endelse
        if(imirt eq 1) then begin
            itt = 2*sizet + 1 - it
            dum2 = dblarr(sizephibig)
            dum2(*) = intensity(*,imu-1,ir-1,it-1)
            intensity(*,imu-1,ir-1,itt-1) = rotate(dum2,2)
        endif
     endfor
  endfor
endfor
close,1
;
; Read the mu data
;
file='mu.dat'
dd=findfile(file,count=count)
if(count le 0) then stop
openr,1,file
readf,1,nmu
mu=dblarr(nmu)
readf,1,mu
close,1
;
; Read the phi data
;
file='phi.dat'
dd=findfile(file,count=count)
if(count le 0) then stop
openr,1,file
readf,1,nphi
phi=dblarr(nphi)
readf,1,phi
close,1
;
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nmu:nmu,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),mu:mu,phi:phi,$
        rr:rr,tt:tt,intensity:intensity}
end


;--------------------------------------------------------------------
;                     READ THE SOURCE
;--------------------------------------------------------------------
function read_source,ext=ext,all=all
grid = read_grid()
niter=0
dum=0
openr,1,'source.info'
readf,1,niter
readf,1,dum
readf,1,lastsave
close,1

iter=0
sizer=0
sizet=0
imirt=0
nfr=0

if not keyword_set(all) then begin
   if(niter eq -2) then niter='final'
   if n_elements(ext) gt 0 then niter=ext
   filename='source_'+strcompress(string(niter),/remove_all)+'.dat'
endif else begin
   filename=['']
   for i=1,lastsave do begin
      filename=[filename,'source_'+strcompress(string(i),/remove_all)+'.dat']
   endfor
   if niter eq -2 then begin
      filename=[filename,'source_final.dat']
   endif
   filename=filename(1:*)
endelse

nn=n_elements(filename)

for in=1,nn do begin
   openr,1,filename[in-1]
   readf,1,nfr,sizer,sizet,imirt
   sizetbig=sizet*(1+imirt)
   if in eq 1 then begin
      source=dblarr(4,sizer,sizetbig,nfr,nn)
   endif
   for inu=1,nfr do begin
      for ir=1,sizer do begin
         for it=1,sizet do begin
            dum = dblarr(4)
            readf,1,dum
            source(*,ir-1,it-1,inu-1,in-1) = dum
            if(imirt eq 1) then begin
               itt = 2*sizet + 1 - it
               source(*,ir-1,itt-1,inu-1,in-1) = dum
            endif
         endfor
      endfor
   endfor
   close,1
endfor
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
therm_src   = dblarr(sizer,sizetbig,nfr,nn)
therm_alpha = dblarr(sizer,sizetbig,nfr,nn)
scat_src    = dblarr(sizer,sizetbig,nfr,nn)
scat_alpha  = dblarr(sizer,sizetbig,nfr,nn)
therm_src(*,*,*,*)   = source(0,*,*,*,*)
therm_alpha(*,*,*,*) = source(1,*,*,*,*)
scat_src(*,*,*,*)    = source(2,*,*,*,*)
scat_alpha(*,*,*,*)  = source(3,*,*,*,*)

return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,therm_src:therm_src,therm_alpha:therm_alpha,$
        scat_src:scat_src,scat_alpha:scat_alpha}
end


;--------------------------------------------------------------------
;                READ THE APPROXIMATE OPERATOR
;--------------------------------------------------------------------
function read_ali
grid = read_grid()
niter=0
dum=0

sizer=0
sizet=0
imirt=0
nfr=0

filename='alidiag.dat'
openr,1,filename
readf,1,nfr,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
ali=dblarr(sizer,sizetbig,nfr)
for inu=1,nfr do begin
   for ir=1,sizer do begin
      for it=1,sizet do begin
         dum = 0.d0
         readf,1,dum
         ali(ir-1,it-1,inu-1) = dum
         if(imirt eq 1) then begin
            itt = 2*sizet + 1 - it
            ali(ir-1,itt-1,inu-1) = dum
         endif
      endfor
   endfor
endfor
close,1
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        rr:rr,tt:tt,ali:ali}
end


;--------------------------------------------------------------------
;                       READ THE FLUX FILE
;--------------------------------------------------------------------
function read_flux
grid = read_grid()
;
openr,1,'fluxcons.dat'
nr=0
nf=0
readf,1,nf,nr
lstar=0.d0
lum=dblarr(nr)
lumnu=dblarr(nr,nf)
readf,1,lstar
readf,1,lum
readf,1,lumnu
close,1
return,{nr:grid.nr,nnu:grid.nnu,$
        r:grid.r,nu:grid.nu(1:*),$
        lum:lum,lumnu:lumnu,lstar:lstar}
end




;--------------------------------------------------------------------
;               READ THE OPTICAL DEPTH AT INFINITY
;--------------------------------------------------------------------
function read_tauinf
grid = read_grid()

imirt=0
sizet=0
nfr=0
dum = 0.d0

openr,1,'tauinf.dat'
readf,1,nfr,sizet,imirt
sizetbig=sizet*(1+imirt)
tauinf=dblarr(nfr,sizetbig)

for it=1,sizet do begin
   for inu=1,nfr do begin
      readf,1,dum
      tauinf(inu-1,it-1) = dum
      if(imirt eq 1) then begin
         itt = 2*sizet + 1 - it
         tauinf(inu-1,itt-1) = dum
      endif
   endfor
endfor
close,1
return,{ntheta:grid.ntheta,nnu:grid.nnu,$
        theta:grid.theta,nu:grid.nu(1:*),$
        tauinf:tauinf}
end


;-----------------------------------------------------------------
;                   READ EDDINGTON FLUX
;-----------------------------------------------------------------
function read_eddflux
grid = read_grid()
dum=0

nrspec=0
sizer=0
sizet=0
imirt=0

openr,1,'eddflux.dat'
readf,1,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
hr=dblarr(sizer,sizetbig)
ht=dblarr(sizer,sizetbig)
idum=0
for ir=1,sizer do begin
   for it=1,sizet do begin
      dum = 0.d0
      readf,1,dum
      hr(ir-1,it-1) = dum
      if(imirt eq 1) then begin
         itt = 2*sizet + 1 - it
         hr(ir-1,itt-1) = dum
      endif
   endfor
endfor
if sizet gt 1 then begin
   for ir=1,sizer do begin
      for it=1,sizet do begin
         dum = 0.d0
         readf,1,dum
         ht(ir-1,it-1) = dum
         if(imirt eq 1) then begin
            itt = 2*sizet + 1 - it
            ht(ir-1,itt-1) = -dum
         endif
      endfor
   endfor
endif
close,1
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        nspec:nrspec,rr:rr,tt:tt,hr:hr,ht:ht}
end


;-----------------------------------------------------------------
;                     READ FLUX CORR FACTOR
;-----------------------------------------------------------------
function read_fluxcorr
grid = read_grid()
dum=0

nrspec=0
sizer=0
sizet=0
imirt=0

openr,1,'fluxcorrfact.dat'
readf,1,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
corr=dblarr(sizer,sizetbig)
idum=0
for ir=1,sizer do begin
   for it=1,sizet do begin
      dum = 0.d0
      readf,1,dum
      corr(ir-1,it-1) = dum
      if(imirt eq 1) then begin
         itt = 2*sizet + 1 - it
         corr(ir-1,itt-1) = dum
      endif
   endfor
endfor
close,1
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        nspec:nrspec,rr:rr,tt:tt,corr:corr}
end


;-----------------------------------------------------------------
;                     READ FLUX DIVERGENCE
;-----------------------------------------------------------------
function read_fluxdiv,ext=ext
grid = read_grid()
dum=0

nrspec=0
sizer=0
sizet=0
imirt=0

filebase='fluxdiv'
if n_elements(ext) gt 0 then begin
   filebase=filebase+'_'+strcompress(string(ext),/remove_all)
endif
filename=filebase+'.dat'

openr,1,filename
readf,1,sizer,sizet,imirt
sizetbig=sizet*(1+imirt)
divh=dblarr(sizer,sizetbig)
idum=0
for ir=1,sizer do begin
   for it=1,sizet do begin
      dum = 0.d0
      readf,1,dum
      divh(ir-1,it-1) = dum
      if(imirt eq 1) then begin
         itt = 2*sizet + 1 - it
         divh(ir-1,itt-1) = dum
      endif
   endfor
endfor
close,1
rr=rebin(grid.r,grid.nr,grid.ntheta)
if grid.ntheta gt 1 then $
   tt=transpose(rebin(grid.theta,grid.ntheta,grid.nr)) $ 
   else tt=0
close,1
;
; Now compute the integrals over the flux divergences
;
dum=0.d0
dum1=0.d0
for ir=2,sizer-1 do begin
   for it=2,sizet-1 do begin
      dum = dum + 2*!pi* 0.25 * grid.r(ir-1)^2 * (grid.r(ir)-grid.r(ir-2)) $
                    * abs(cos(grid.theta(it))-cos(grid.theta(it-2))) $ 
                    * divh(ir-1,it-1)
      dum1 = dum1 + 2*!pi* 0.25 * grid.r(ir-1)^2 * (grid.r(ir)-grid.r(ir-2)) $
                    * abs(cos(grid.theta(it))-cos(grid.theta(it-2))) $ 
                    * abs(divh(ir-1,it-1))
   endfor
endfor
if imirt eq 1 then begin
   dum  = dum * 2
   dum1 = dum1 * 2
endif
dum  = dum * 4*!dpi
dum1 = dum1 * 4*!dpi
;
; Read star info
;
openr,1,'starinfo.inp'
readf,1,iformat
readf,1,rstar
readf,1,mstar
readf,1,tstar
close,1
lstar=4*!dpi*rstar^2*5.6703e-5*tstar^4
;
return,{nr:grid.nr,ntheta:grid.ntheta,nnu:grid.nnu,$
        r:grid.r,theta:grid.theta,nu:grid.nu(1:*),$
        nspec:nrspec,rr:rr,tt:tt,divh:divh,extralum:dum,$
        extraabs:dum1,lstar:lstar}
end





;=================================================================
;                     SPECTRUM PLOT ROUTINES     
;=================================================================

;-----------------------------------------------------------------
;                        READ THE SPECTRUM
;-----------------------------------------------------------------
function read_spectrum,ext=ext,file=file,dpc=dpc
fr_units=0
;
; Read possible dust information
;
z=findfile('dustopac.inp',count=count)
dust=count
if dust gt 0 then begin
    print,"Found dustopac.inp, so assuming dust spectrum"
    fr_units = -1
endif

;
; Read the total spectrum
;
if not keyword_set(file) then begin
   if n_elements(ext) eq 0 then begin
      filename='spectrum.dat'
   endif else begin
      filename='spectrum_'+strcompress(string(ext),/remove_all)+'.dat'
   endelse
endif else begin
   filename=file
endelse
openr,1,filename
nrfr=1
readf,1,nrfr
dum=dblarr(2,nrfr)
freq=dblarr(nrfr+1)
spectrum=dblarr(nrfr+1)
readf,1,dum
freq(1:nrfr) = dum(0,0:nrfr-1)
spectrum(1:nrfr) = dum(1,0:nrfr-1)
close,1
;
; If dpc is specified, then the spectrum is not the standard d=1pc
; spectrum, but in fact is a spectrum at a given distance (for instance,
; if you want to use the new option of RADICAL to include an aperture,
; then RADICAL returns the actual spectrum at d=dpc distance). In order
; not to confuse the plot_spectrum() routine one can rescale it back
; to a spectrum at 1pc by specifying dpc here in this read routine. 
;
if keyword_set(dpc) then spectrum=spectrum*dpc^2
;
; Return all...
;
a={spectrum:spectrum,nfr:nrfr,freq:freq,fr_units:fr_units,obs:0}
return,a
end


;-----------------------------------------------------------------
;                   READ AN OBSERVED SPECTRUM
;-----------------------------------------------------------------
function read_obsspec,file,hz=hz,micron=micron,jy=jy,erg=erg,$
      nufnu=nufnu,npt=npt
cc   = 2.9979d10
fr_units = 0
if (not keyword_set(hz)) and (not keyword_set(micron)) then begin
   print,"Must specify either /hz or /micron"
   return,0
endif
if (not keyword_set(jy)) and (not keyword_set(erg)) $
    and (not keyword_set(nufnu)) then begin
   print,"Must specify either /jy or /erg or /nufnu"
   return,0
endif
openr,1,file
if not keyword_set(npt) then readf,1,npt
data=dblarr(2,npt)
readf,1,data
close,1
freq=dblarr(npt+1)
spectrum=dblarr(npt+1)
if keyword_set(hz) then begin
   freq(1:npt) = data(0,0:npt-1)
endif 
if keyword_set(micron) then begin
   freq(1:npt) = cc / (1d-4*data(0,0:npt-1))
   fr_units = -1
endif
if keyword_set(erg) then begin
   spectrum(1:npt) = data(1,0:npt-1)
endif
if keyword_set(jy) then begin
   spectrum(1:npt) = 1d-23 * data(1,0:npt-1)
endif
if keyword_set(nufnu) then begin
   spectrum(1:npt) = data(1,0:npt-1) / freq(1:npt)
endif
a={spectrum:spectrum,nfr:npt,freq:freq,fr_units:fr_units,obs:1}
return,a
end



;-----------------------------------------------------------------
;                       PLOT SPECTRUM 
;
; ARGUMENTS:
;   ev              = 1 --> frequency in electronvolt (default=Hz)
;   kev             = 1 --> frequency in kiloelectronvolt (default=Hz)
;   lnu             = 1 --> L_nu (default nu*L_nu)
;   rpc             = Distance of observer in units of parsec
;                   = 0 (or undef) --> Lum = Total Lum in erg/s
;                     instead of erg / s cm^2
;   lumtot          = 1 --> total luminosity (default=flux at 1 parsec)
;
;-----------------------------------------------------------------
pro plot_spectrum,a,ev=ev,kev=kev,hz=hz,micron=micron,$
        lnu=lnu,fnu=fnu,nulnu=nulnu,nufnu=nufnu,rpc=rpc,$
        dpc=dpc,xlg=xlg,oplot=oplot,jy=jy,lsun=lsun,$
        itheta=itheta,ldat=ldat,ylin=ylin,lum=lum,$
        thick=thick,xthick=xthick,ythick=ythick,$
        charthick=charthick,charsize=charsize,_extra=_extra
if not keyword_set(itheta) then itheta=0
fr_units = a.fr_units
if keyword_set(micron) then begin
    fr_units = -1
endif
if keyword_set(hz) then begin
    fr_units = 0
endif
if keyword_set(ev) then begin
    fr_units = 1
endif
if keyword_set(kev) then begin
    fr_units = 2
endif
case fr_units of
    -1: begin
        xcoord = 2.9979d14 / a.freq(1:*)
        xtitle = '!4k!X [!4l!Xm]'
    end
    0: begin
        xcoord = a.freq(1:*)
        xtitle = '!4m!X [Hz]'
    end
    1: begin
        xcoord = 4.13568842841d-15 * a.freq(1:*)
        xtitle = '!4m!X [eV]'
    end
    2: begin
        xcoord = 4.13568842841d-18 * a.freq(1:*)
        xtitle = '!4m!X [KeV]'
    end
endcase

if keyword_set(dpc) then rpc=dpc
;
; Plot nuFnu or Fnu (same with Lnu)? And what about Fnu vs Lnu?
;
sed=1
ylum=0
if keyword_set(jy) then begin
   sed=0
endif
if keyword_set(fnu) then begin
   ylum=0
   sed=0
endif
if keyword_set(lnu) then begin
   ylum=1
   sed=0
endif
if keyword_set(nufnu) then begin
   ylum=0
   sed=1
endif
if keyword_set(nulnu) then begin
   ylum=1
   sed=1
endif
if keyword_set(lum) then ylum=1 
if keyword_set(jy) then begin
   ylum=0
   ;sed=0
endif
if keyword_set(lsun) then begin
   ylum=1
   sed=1
endif
;
; If the distance is not given, then for an observed spectrum
; we can only plot the flux, while for a computed spectrum we
; can only plot the luminosity
;
if a.obs eq 0 then begin
   if not keyword_set(rpc) and keyword_set(jy) then begin
      print,"Cannot use Jy units of theoretical spectrum if dpc not known."
      return
   endif
   if not keyword_set(rpc) and ylum eq 0 then begin
      ylum=1
      print,"Do not know distance to source. Plotting Luminosity instead."
   endif
endif else begin
   if not keyword_set(rpc) and ylum eq 1 then begin
      ylum=0
      print,"Do not know distance to source. Plotting Flux instead."
   endif
endelse
;
; Which plot to make? Lum or flux?
;
if ylum eq 0 then begin
   ;
   ; Plot spectrum as flux at a certain distance
   ;
   if a.obs eq 0 then begin
      rpc=1.d0*rpc
      distfact = 1.d0/ (rpc^2)
   endif else begin
      distfact = 1.d0
   endelse
   if not keyword_set(jy) then begin
      if sed eq 0 then begin
         lumfact=1.d0
         ytitle='F!I!4m!X!N  [erg cm!E-2!N Hz!E-1!N s!E-1!N]'
      endif else begin
         lumfact=a.freq(1:*)
         ytitle='!4m!XF!I!4m!X!N [erg cm!E-2!N s!E-1!N]'
      endelse
   endif else begin
      if sed eq 0 then begin
         lumfact=1d+23
         ytitle='F!I!4m!X!N [Jy]'
      endif else begin
         lumfact=1d+23*a.freq(1:*)
         ytitle='!4m!XF!I!4m!X!N [JyHz]'
      endelse
   endelse
endif else begin
   ;
   ; Plot spectrum as luminosity
   ;
   if a.obs eq 0 then begin
      distfact = 1.1965280793d38 ; = (1 parsec)^2 = 1.19d38 cm^2
   endif else begin
      rpc=1.d0*rpc
      distfact = rpc^2 * 1.1965280793d38
   endelse
   if sed eq 0 then begin
      lumfact=1.d0
      ytitle='L!I!4m!X!N  [erg Hz!E-1!N s!E-1!N]'
   endif else begin
      if not keyword_set(lsun) then begin
         lumfact=a.freq(1:*)
         ytitle='!4m!XL!I!4m!X!N [erg s!E-1!N]'
      endif else begin
         lumfact=a.freq(1:*)*2.5956986d-34
;         ytitle='!4m!XL!I!4m!X!N [Lsun]'
         ytitle='!4m!XL!I!4m!X!N [L!D!9n!X!N]'
      endelse
   endelse
endelse
if not keyword_set(xlg) then begin
    xlog=1
endif else begin
    xcoord=alog10(xcoord)
    xtitle = 'log '+xtitle
    xlog=0
endelse
if not keyword_set(ylin) then begin
    ylog=1
endif else begin
    ylog=0
endelse
if not keyword_set(oplot) then begin
    plot,xcoord,distfact*lumfact*a.spectrum(1:*,itheta),$
         xtitle=xtitle,ytitle=ytitle,$
         xlog=xlog,ylog=ylog,thick=thick,xthick=xthick,ythick=ythick,$
        charthick=charthick,charsize=charsize,_extra=_extra
endif else begin
    oplot,xcoord,distfact*lumfact*a.spectrum(1:*,itheta),$
         thick=thick,_extra=_extra
endelse
end




;------------------------------------------------------------------------
;                       FUNCTION: PLANCK
;------------------------------------------------------------------------
function bplanck,nnu,TT
nu=nnu*1.d0
T=TT*1.d0
cc=2.9979d10
hh=6.6262d-27
kk=1.3807e-16
n=n_elements(nu)
bpl=dblarr(n)
for i=0,n-1 do begin
  x=hh*nu[i]/(kk*T)
  ;print,x
  if x gt 100.0 then b=0.d0
  if x gt 1.d-3 then b=(2.d0*hh*nu[i]^3/cc^2)/(exp(x)-1.d0) $
    else b=2.0*nu[i]^2*kk*T/cc^2
  bpl[i]=b
endfor
return,bpl
end



;-----------------------------------------------------------------
;               PLOT A BLACKBODY SPECTRUM FOR A STAR
;-----------------------------------------------------------------
pro plot_bb,temp,rsrsun,spectrum=spectrum,_extra=_extra
kk = 1.3807d-16
hh = 6.6262d-27
pi = !dpi
pc = 3.08572d18
fr_units=0
if not keyword_set(rsrsun) then rsrsun=1.0
rstar=6.96d10*rsrsun
nrfr=100
nu1=100*(kk/hh)*temp
nu0=1.d-5*(kk/hh)*temp
freq=nu0*(nu1/nu0)^(dindgen(nrfr)/(nrfr-1.d0))
spectrum=rstar^2*pi*bplanck(freq,temp)/pc^2
freq=[0.d0,freq]
spectrum=[0.d0,spectrum]
fr_units=-1
;
spectrum={spectrum:spectrum,nfr:nrfr,freq:freq,fr_units:fr_units,obs:0}
plot_spectrum,spectrum,_extra=_extra
end


;-----------------------------------------------------------------
;                         PLOT THE GRID
;-----------------------------------------------------------------
pro plotgrid,a,data=data,rmax=rmax,ymax=ymax,ax=ax,az=az,$
             nlevels=nlevels,_extra=_extra
if not keyword_set(rmax) then rmax=max(a.r)
if not keyword_set(az) then az=0
if not keyword_set(ax) then ax=90
if not keyword_set(ymax) then ymax=rmax
if not keyword_set(nlevels) then nlevels=100
nt=n_elements(a.theta)
if max(a.theta) gt !pi/2.d0 then begin
   nt=nt/2
endif
if not keyword_set(data) then begin
   data=1.0+0*a.rr
   surface,data(*,0:nt-1),a.rr(*,0:nt-1)*sin(a.tt(*,0:nt-1)),$
     a.rr(*,0:nt-1)*cos(a.tt(*,0:nt-1)),ax=ax,az=az,$
     xrange=[0,rmax],yrange=[0,ymax],_extra=_extra
endif else begin 
   contour,data(*,0:nt-1),a.rr(*,0:nt-1)*sin(a.tt(*,0:nt-1)),$
     a.rr(*,0:nt-1)*cos(a.tt(*,0:nt-1)),$
     xrange=[0,rmax],yrange=[0,ymax],nlevels=nlevels,_extra=_extra
endelse
end
