pro makemodels
  Tstar = [10000.0,   7400.87,    5477.28,    4053.66,   9520.0,      3000.06 ]
  Marr = [3.16192,   1.68052,   0.893178,   0.474713,      2.4,     0.252304  ]
  Larr = [39.3626,   5.23934,   0.697379,  0.0928243,     47.0,    0.0123553  ]

  freq = [ 2.997899D+11, 3.787589D+11, 4.785287D+11, 6.045801D+11, $
           7.638339D+11, 9.650388D+11, 1.219244D+12, 1.540408D+12, $
           1.946174D+12, 2.458824D+12, 3.106509D+12, 3.924809D+12, $
           4.958651D+12, 6.264832D+12, 7.915065D+12, 1.000001D+13, $
           1.063858D+13, 1.131793D+13, 1.204068D+13, 1.280957D+13, $
           1.362755D+13, 1.449776D+13, 1.542355D+13, 1.640848D+13, $
           1.745628D+13, 1.857099D+13, 1.975688D+13, 2.101850D+13, $
           2.236072D+13, 2.378861D+13, 2.530768D+13, 2.692376D+13, $
           2.864303D+13, 3.047210D+13, 3.241802D+13, 3.448814D+13, $
           3.669045D+13, 3.903340D+13, 4.152596D+13, 4.417777D+13, $
           4.699884D+13, 5.000005D+13, 6.276820D+13, 7.879686D+13, $
           9.891865D+13, 1.241788D+14, 1.558894D+14, 1.956978D+14, $
           2.456717D+14, 3.084072D+14, 3.871629D+14, 4.860299D+14, $
           6.101439D+14, 7.659519D+14, 9.615475D+14, 1.207091D+15, $
           1.515337D+15, 1.902298D+15, 2.388074D+15, 2.997899D+15  ]

  Msun = 1.989d33
  Lsun = 3.86d33
  Rsun = 6.960d10
  sigma = 5.671d-5
  grav = 6.673d-8
  cl = 2.998d10
  hpl = 6.626d-27
  bk = 1.381-16
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
      plot,nu,lnusmooth,/xlog,/ylog
      
      lint = interpol(lnusmooth,nu,freq)
      oplot,freq,lint,psy=2
      bnu = freq
      for k=0,n_elements(freq)-1 do bnu(k) =bplanck(freq[k],Tstar[i])
      print,bnu
      oplot,freq,bnu*4*!DPI^2*Rstar[i]^2
      openw,2,string(format='(%"flux%d.dat")',floor(Tstar[i]))
      for j = 0,n_elements(lint)-1 do begin
          printf,2,freq[j],lint[j]
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
  
