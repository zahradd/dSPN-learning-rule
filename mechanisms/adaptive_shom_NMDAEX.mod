TITLE simple NMDA receptors





NEURON {
	POINT_PROCESS adaptive_shom_NMDAEX
	RANGE g, Alpha, Beta, Erev, gmax, Cdur, iNMDA,mg, Cmax, eta, alpha, Cdur_init, Cdur_factor, treshf
	NONSPECIFIC_CURRENT  iNMDA
	POINTER weight
	USEION ca_nmda READ ca_nmdai WRITE ica_nmda VALENCE 2	
	USEION cal READ cali VALENCE 2

}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	Cmax	= 1	 (mM)           : max transmitter concentration
	Alpha	= 4 (/ms /mM)	: forward (binding) rate (4)
	Beta 	= 0.01   (/ms)   : backward (unbinding) rate
	Erev	= 0	 (mV)		: reversal potential
    mg   = 0      (mM)           : external magnesium concentration
    eta = 0 (/mV)
    alpha = 0 (/mV)
	gmax = 0   (uS)
	Cdur_init = 50
	Cdur_factor =200


}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	iNMDA 		(nA)		: current = g*(v - e)
	g 		(uS)		: conductance
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	synon
    B                       : magnesium block
	ica_nmda        (nA)
	weight
    cali            (mM)
	ca_nmdai        (mM)
	gluextra
	Cdur (ms)		: transmitter duration (rising phase)
}

STATE {Ron Roff }

INITIAL {
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / (Cmax*Alpha + Beta)
	synon = 0
	Cdur	= Cdur_init + weight*Cdur_factor	: transmitter duration (rising phase)
}
BREAKPOINT {
	SOLVE release METHOD cnexp
    B = mgblock(v)
	g = (Ron + Roff)* gmax * B
	iNMDA = g*(v - Erev)
    ica_nmda = 0.175*iNMDA  :(5-10 times more permeable to Ca++ than Na+ or K+, Ascher and Nowak, 1988)
    iNMDA = 0.825*iNMDA
	Cdur = Cdur_init + weight*Cdur_factor

}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
}

FUNCTION mgblock(v(mV)) {
        TABLE
        DEPEND mg
        FROM -140 TO 80 WITH 1000

        : from Jahr & Stevens


	 mgblock = 1 / (1 + mg * eta * exp(-alpha * v) )  :was 0.062, changed to 0.072 to get a better voltage-dependence of NMDA currents, july 2008, kiki

}

: following supports both saturation from single input and
: summation from multiple inputs
: if spike occurs during CDur then new off time is t + CDur
: ie. transmitter concatenates but does not summate
: Note: automatic initialization of all reference args to 0 except first


NET_RECEIVE(dummy, on, nspike, r0, t0 (ms)) {
	: flag is an implicit argument of NET_RECEIVE and  normally 0
        if (flag == 0) { : a spike, so turn on if not already in a Cdur pulse

		nspike = nspike + 1
		if (!on) {
			r0 = r0*exp(-Beta*(t - t0))
			t0 = t
			on = 1
			synon = synon + weight

			state_discontinuity(Ron, Ron + r0)
			state_discontinuity(Roff, Roff - r0)
		}
:		 come again in Cdur with flag = current value of nspike
		net_send(Cdur, nspike)
       }
	if (flag == nspike) { : if this associated with last spike then turn off

		r0 = weight*Rinf + (r0 - weight*Rinf)*exp(-(t - t0)/Rtau)
		t0 = t
		synon = synon - weight
		state_discontinuity(Ron, Ron - r0)
		state_discontinuity(Roff, Roff + r0)
		on = 0
	}
}





