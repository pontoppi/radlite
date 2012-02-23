FUNCTION Pesc, tau
  IF tau LE 0. THEN Pesc = 0.5
  IF tau GT 0. AND tau LE 0.6 THEN Pesc = 0.5+(0.1995*ALOG(tau)-0.2484)*tau-0.04594*tau^2d0
  IF tau GT 0.6 AND tau LE 9.0 THEN Pesc = (1d0-EXP(-1.422*tau))/(3.324*tau+0.2852*tau^2d0)
  IF tau GT 9.0 THEN Pesc = 0.1999/tau * (ALOG(0.4799*tau))^(-0.4195)
  RETURN, Pesc
END
