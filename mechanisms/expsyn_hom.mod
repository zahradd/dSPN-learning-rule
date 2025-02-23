NEURON {
	POINT_PROCESS ExpSyn_Hom
	RANGE tau, e, i, w0
	NONSPECIFIC_CURRENT i
        POINTER Rs
        RANGE tau_w
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau = 0.1 (ms) <1e-9,1e9>
	e = 0	(mV)
        tau_w = 6 (s)
        w0 = 0.5
}

ASSIGNED {
	v (mV)
	i (nA)
        Rs
}

STATE {
	g (uS)
        w
}

INITIAL {
	g=0
        w = w0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - e)
}

DERIVATIVE state {
	g' = -g/tau
        w' = -Rs*w/tau_w
}

NET_RECEIVE(dummy (uS)) {
	g = g + w
}
