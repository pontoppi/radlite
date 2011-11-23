PRO callback, status, error, bridge, ud
	npop     = bridge->getvar('npop')
        npop_ini = bridge->getvar('npop_ini')
	(*(ud.npop_all))[*,*, ud.i]     = npop
	(*(ud.npop_ini_all))[*,*, ud.i] = npop_ini
END
