pro makemodels
;
; Parameters of the stars
;
  Tstar    = [30000.]
  Marr     = [4.    ]
  Larr     = [3000. ]
;
; Parameters of the frequency grid
;
  lam0     = 1d4
  lam1     = 30.d0
  lam2     = 6.0
  lam3     = 0.01d0
  n01      = 24
  n12      = 22
  n23      = 30
;
;=============================================================
;            DO THE REMAPPING OF THE KURUCZ MODELS
;=============================================================
;
  Msun     = 1.989d33
  Lsun     = 3.86d33
  Rsun     = 6.960d10
  sigma    = 5.671d-5
  grav     = 6.673d-8
  cl       = 2.998d10
  hpl      = 6.626d-27
  bk       = 1.381-16
  pc       = 3.08572d18

  close,/all

  cc       = 2.9979e10
  freq0    = 1e4*cc/lam0
  freq1    = 1e4*cc/lam1
  freq2    = 1e4*cc/lam2
  freq3    = 1e4*cc/lam3
  dlfreq01 = alog(freq1)-alog(freq0)
  dlfreq12 = alog(freq2)-alog(freq1)
  dlfreq23 = alog(freq3)-alog(freq2)
  nfrnew01 = n01
  nfrnew12 = n12
  nfrnew23 = n23
  freqnw01 = exp((findgen(nfrnew01)/(nfrnew01))*dlfreq01+alog(freq0))
  freqnw12 = exp((findgen(nfrnew12)/(nfrnew12))*dlfreq12+alog(freq1))
  freqnw23 = exp((findgen(nfrnew23)/(nfrnew23-1.e0))*dlfreq23+alog(freq2))
  freq     = [freqnw01,freqnw12,freqnw23]
  nf       = n_elements(freq)

  Mstar = Marr*Msun
  Lstar = Larr*Lsun
  Rstar = sqrt(Lstar/4./!DPI/sigma/Tstar^4)
  print,Rstar
  logg = alog10(grav*Mstar/Rstar^2)
  print, logg
  for i=0,n_elements(Mstar) - 1 do begin
      cmd = string(format='(A50,1x,F6.0,1x,F5.1,1x,F5.1,A10)', $
                   "./interpolate_kurucz -q", $
                   Tstar[i],logg[i],0.0, $
                   "> flux.dat")
      spawn,cmd
      openr,1,'flux.dat'
      a=dblarr(2,1221)
      readf,1,a
      close,1
      loglam = a[0,*]
      logflux = a[1,*]
      flux = 10^logflux
      lam = 10^loglam*1e-4
      nu = cl/lam
      s1 = nu
      s2 = flux
      LL = int_tabulated(s1,s2,/sort,/double)
      lnu = flux*Lstar[i]/LL
      s1 = nu
      s2 = lnu
      LL1 = int_tabulated(s1,s2,/double,/sort)
      print,LL,LL1,LL1/Lstar[i]
      lnusmooth = lnu
      imin = 0
      imax=1200
      ;; Smooth the spectrum
;      lnusmooth[imin:imax] = smooth(lnusmooth[imin:imax],11)
      print,nu[0],nu[imax]
      print,lam[0],lam[imax]
      plot,3d14/nu,nu*lnusmooth/Lsun,/xlog,/ylog
      
      lint = interpol(lnusmooth,nu,freq)
      ii=where(lint lt 0.d0)
      lint[ii]=0.d0

      oplot,3d14/freq,freq*lint/Lsun,psy=2
      bnu = freq
      for k=0,n_elements(freq)-1 do bnu(k) =bplanck(freq[k],Tstar[i])
      print,bnu
      oplot,freq,bnu*4*!DPI^2*Rstar[i]^2
      filename='starspectrum_'+strcompress(string(floor(Tstar[i])),/remove_all)+'.inp'
      openw,2,filename
      printf,2,n_elements(freq)
      factor=1.d0/(4*!dpi*pc^2)
      for j = 0,n_elements(lint)-1 do begin
          printf,2,freq[j],lint[j]*factor
      endfor
      close,2
      s1 = nu
      s2 = lnu
      LL1 = int_tabulated(s1,s2,/double,/sort)
      print, "final",LL1/Lstar[i]


  endfor   
  
end
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
  
