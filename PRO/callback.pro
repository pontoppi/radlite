PRO callback, status, error, bridge, ud
	npop     = bridge->getvar('npop')
        ini_npop = bridge->getvar('ini_npop')
	(*(ud.p_npop_all))[*,*, ud.i]     = npop
	(*(ud.p_npop_ini_all))[*,*, ud.i] = ini_npop
END
