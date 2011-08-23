FUNCTION Ppump, tau
  IF tau LE 0 THEN Ppump = 1d0
  IF tau GT 0 AND tau LE 0.9 THEN Ppump = 1d0-0.3989*tau+0.09189*tau^2d0-0.01497*tau^3d0
  IF tau GT 0.9 and tau LE 9.0 THEN Ppump = (1d0-EXP(-0.6437*tau))/(0.6295*tau+0.07008*tau^2d0)
  IF tau GT 9.0 THEN Ppump = (0.8204/tau) * (ALOG(0.3367*tau))^(-0.4306)
  RETURN, Ppump
END
