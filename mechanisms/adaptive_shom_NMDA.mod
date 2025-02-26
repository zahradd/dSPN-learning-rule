TITLE simple NMDA receptors
NEURON {
	POINT_PROCESS adaptive_shom_NMDA
	RANGE g, Alpha, Beta, Erev, gmax, Cdur, iNMDA,mg, Cmax, eta, alpha,w0,  treshf,synon,flagx,conc0,mltypeMin,mltype, rate_ltd,rate_ltp,rate_ltd_thrsh,rate_ltp_tresh,tremin,synapticActiveTime
	NONSPECIFIC_CURRENT  iNMDA
	POINTER dopa
	USEION ca_nmda READ ca_nmdai WRITE ica_nmda VALENCE 2	
	USEION cal READ cali VALENCE 2
	USEION ca READ cai VALENCE 2

}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	Cmax	= 1	 (mM)           : max transmitter concentration
	Cdur      =1  (ms)		: transmitter duration (rising phase)
	Alpha	= 4 (/ms /mM)	: forward (binding) rate (4)
	Beta 	= 0.01   (/ms)   : backward (unbinding) rate
	Erev	= 0	 (mV)		: reversal potential
    mg   = 0      (mM)           : external magnesium concentration
    eta = 0  :change in the code
    alpha = 0 (/mV)
	gmax = 0   (uS)
    w0 = 0
    conc0=0 (mM)
    flagx=0 (1)
    mltype=0 (mM)
    mltypeMin=0.00015 (mM)
    rate_ltd=0 (1)
    rate_ltp=0 (1)
    rate_ltd_thrsh=0 (1)
    rate_ltp_tresh=0 (1)
    treshf=0 (mM)    :change in the code
    treshltp=0.02(mM)
    tremin=0 (mM)
    f=0   (ms)
    conc   (mM)



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
	dopa
    cali            (mM)
	ca_nmdai        (mM)
	cai  (mM)
	synapticActiveTime 

	
}

STATE {Ron Roff weight kernel_MidPoint }

INITIAL {
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / (Cmax*Alpha + Beta)
	synon = 0
	weight = w0
    kernel_MidPoint = treshf
    synapticActiveTime = 0

    
   

}

BREAKPOINT {
	SOLVE release METHOD cnexp
    B = mgblock(v)
	g = (Ron + Roff)* gmax*B
	iNMDA = g*(v - Erev)
    ica_nmda = 0.175*iNMDA   :(5-10 times more permeable to Ca++ than Na+ or K+, Ascher and Nowak, 1988)
    iNMDA = 0.825*iNMDA
    weight=weight+deltaW(dopa,weight,ca_nmdai,cali,cai,g,kernel_MidPoint)
    kernel_MidPoint=kernel_MidPoint+MP_kernel(dopa,kernel_MidPoint,conc0,mltype)

   
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
        synapticActiveTime = 1  : set synapse to active
    	: Schedule an event to reset synapticActiveTime in 400 ms
    	net_send(400, -1)
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
	if (flag == -1) {  : Special flag to reset synapticActiveTime
		synapticActiveTime = 0  : reset synapse to inactive
	}
}


FUNCTION deltaW(dopam,we(uS),conc(mM),calic(mM),cai(mM),g_p(uS),kernel_MidPoint(mM))(uS/ms) {
    UNITSOFF
    
    
    if (g_p < 1e-7) {
       conc0=0
       mltype=0
    } else {
    
       if (conc>conc0 ){
          conc0=conc
       }
       if (calic>mltype ){
          mltype=calic
       }
    }
    
	deltaW =((rate_ltd*we*mltype*LTD_kernel_function(mltype)*(dopam-1))+(rate_ltp*LTP_kernel_function(conc0,kernel_MidPoint)*(dopam+1)))*synapticActiveTime*dopam*dopam
	
	UNITSON    
}
FUNCTION MP_kernel(dopam,kernel_MidPoint(mM),conc0,mltype)(uM/ms) {
    UNITSOFF
    MP_kernel =((-rate_ltd_thrsh*MP_kernel_function(conc0,kernel_MidPoint)*(dopam-1))-rate_ltp_tresh*(MP_kernel_function(conc0,kernel_MidPoint)*(dopam+1)))*synapticActiveTime*dopam*dopam
	UNITSON    
}





FUNCTION LTP_kernel_function(calnm(mM),kernel_MidPoint(mM))() {
    UNITSOFF
    LTP_kernel_function= (1-(1 / (1 + exp((-1000*calnm + (1000*kernel_MidPoint))/1))))*(1 / (1 + exp((-1000*calnm+(1000*kernel_MidPoint))/1))) 
    UNITSON
}
FUNCTION MP_kernel_function(calnm(mM),kernel_MidPoint(mM))() {
    UNITSOFF
    MP_kernel_function= (1-(1 / (1 + exp((-1000*calnm + (1000*kernel_MidPoint))/3))))*(1 / (1 + exp((-1000*calnm+(1000*kernel_MidPoint))/3)))
    UNITSON
}

FUNCTION LTD_kernel_function(calnm(mM))() {
    UNITSOFF
    LTD_kernel_function=1 / (1 + exp((-Norm(calnm)+7)/1))
    UNITSON
}


FUNCTION Norm(calnm(mM))() {
    UNITSOFF
    Norm= 100000*calnm
    UNITSON
}
