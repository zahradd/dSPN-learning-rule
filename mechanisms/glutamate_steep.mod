COMMENT
Updated Exp2Syn synapse with Mg-blocked nmda channel.

Defaul values of parameters (time constants etc) set to match synaptic channels in 
striatal medium spiny neurons (Du et al., 2017; Chapman et al., 2003; Ding et al., 2008).

Robert . Lindroos @ ki . se

original comment:
________________
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT






NEURON {
	POINT_PROCESS glutamate_steep
	RANGE tau1_ampa, tau2_ampa, tau1_nmda, tau2_nmda
	RANGE erev_ampa, erev_nmda, g, i
	NONSPECIFIC_CURRENT i
	
	RANGE i_ampa, i_nmda, g_ampa, g_nmda, I, G, mg, q,ratio,mltypeMin
	
    POINTER dopa,expectt
	RANGE eta, w0, wmax,treshf,conc0,flagx,mltype,alpha
	
	USEION ca_nmda READ ca_nmdai WRITE ica_nmda VALENCE 2
	USEION cal READ cali VALENCE 2
}


UNITS {
	(molar) = (1/liter)
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (millimolar)
}


PARAMETER {
	erev_ampa        = 0.0       (mV)
	erev_nmda        = 15.0       (mV)
	
	tau1_ampa   = 1.9       (ms)
    tau2_ampa   = 4.8       (ms)  : tau2 > tau1
    tau1_nmda   = 2.75      (ms)  : old value was 5.63
    tau2_nmda   = 115.5      (ms)  : tau2 > tau1
    tau   = 5      (ms)

    treshf=0.02 (mM)
    treshltp=0.05(mM)
    tremin=0.05 (mM)

    mg          = 0   (mM)  :change in the code
    alpha       = 0    :change in the code
    q           = 2
    w0 = 0.0000188 (uS)
    wmax = 0.0015 (uS)
    ratio       = 1         (1)   : both types give same maximal amplitude of current
    conc0=0 (mM)
    flagx=0 (1)
    mltype=0 (mM)
    eta=0   :change in the code
    intresh=0.0051  (mM)
    mltypeMin=0.0001 (mM)


}


ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor_nmda
	factor_ampa
	i_ampa
	i_nmda
	g_ampa (uS)
	g_nmda (uS)
	block
	I
	G
	dopa
	expectt
	ica_nmda (nA)
	ca_nmdai (mM)
	cali (mM)
}


STATE {
	A (uS)
	B (uS)
	C (uS)
	D (uS)
	weight (uS)
	tresh (mM)

}



INITIAL {
	LOCAL tp
	if (tau1_nmda/tau2_nmda > .9999) {
		tau1_nmda = .9999*tau2_nmda
	}
	if (tau1_ampa/tau2_ampa > .9999) {
		tau1_ampa = .9999*tau2_ampa
	}
	
	: NMDA
	A           = 0
	B           = 0
	tp          = (tau1_nmda*tau2_nmda)/(tau2_nmda - tau1_nmda) * log(tau2_nmda/tau1_nmda)
	factor_nmda = -exp(-tp/tau1_nmda) + exp(-tp/tau2_nmda)
	factor_nmda = 1/factor_nmda
	
	: AMPA
	C           = 0
	D           = 0
	tp          = (tau1_ampa*tau2_ampa)/(tau2_ampa - tau1_ampa) * log(tau2_ampa/tau1_ampa)
	factor_ampa = -exp(-tp/tau1_ampa) + exp(-tp/tau2_ampa)
	factor_ampa = 1/factor_ampa

	weight = w0
	tresh=treshf

}




BREAKPOINT {
	SOLVE state METHOD cnexp
	
	: NMDA
	g_nmda = (B - A)
	block  = MgBlock()
	i_nmda = g_nmda * (v - erev_nmda) * block
	ica_nmda = i_nmda
	
	: AMPA
	g_ampa = (D - C)
	i_ampa = g_ampa * (v - erev_ampa)
	
	: total current
	G = g_ampa + g_nmda
	I = i_ampa
    i = I
}



DERIVATIVE state {
	A' = -A/tau1_nmda*q
	B' = -B/tau2_nmda*q
	C' = -C/tau1_ampa
	D' = -D/tau2_ampa
	weight'= supra(dopa,expectt,weight,ca_nmdai,cali,g_nmda,tresh)
	tresh'=sTresh(dopa,tresh,conc0)

}




NET_RECEIVE(dummy (uS)) {
	A = A + weight*factor_nmda*ratio
	B = B + weight*factor_nmda*ratio
	C = C + weight*factor_ampa
	D = D + weight*factor_ampa
}



FUNCTION MgBlock() {
    
    MgBlock = 1 / (1 + mg * exp(-alpha * v) / eta)
    
}

FUNCTION supra(dopam,expec,we(uS),conc(mM),calic(mM),g_p(uS),tresh(mM))(uS/ms) {
    
    UNITSOFF
    if (expec==1){
       flagx=1
    }
    if(expec==0 && flagx==1){
       conc0=0
       mltype=0
       flagx=0
    }
    if(expec==0 && flagx==0){
       if (conc>conc0){
          conc0=conc
       }
       if (calic>mltype){
          mltype=calic
       }
    }
    if (dopa==2){
       supra=0.0005*funcCal(conc0,tresh)*(wmax-we)*trap(g_p)
    }else{
	   supra =((0.5*we*mltype*funcCalMin(mltype)*(dopam-1))+(0.0005*(wmax-we)*funcCal(conc0,tresh)*(dopam+1)))*trap(g_p)*dopam*dopam
	}
	UNITSON    
}
FUNCTION sTresh(dopam,tresh(mM),conc0)(uM/ms) {
    UNITSOFF
    if (dopa==2){
       sTresh=0
    }else{
	   sTresh =dopam*dopam*((0.006*(tresh-minTres(conc0,tresh))*(dopam-1))+(0.006*(Calter(conc0,tresh)-tresh)*(dopam+1)))*trap(g_nmda)
	}
	UNITSON    
}


FUNCTION trap(g_p(uS))() {
    UNITSOFF
	if (g_p < 1e-7) {
		trap = 0
	} else {
	    trap = 1
        }
    UNITSON
}
FUNCTION Calter(calnm(mM),tresh(mM))(mM)  {
    UNITSOFF
	if (calnm < tresh+treshltp) {
		Calter = tresh
	} else {
	    Calter= calnm-treshltp
        }
    UNITSON
}
FUNCTION minTres(calnm(mM),tresh(mM))(mM)  {
    UNITSOFF
	if (calnm < tresh-tremin) {
		minTres = calnm+tremin
	} else {
	    minTres= tresh
        }
    UNITSON
}

FUNCTION funcCal(calnm(mM),tresh(mM))() {
    UNITSOFF
    funcCal= 4*(1-(1 / (1 + exp((-1000*calnm + (1000*tresh))/5))))*(1 / (1 + exp((-1000*calnm+(1000*tresh))/5))) 
    UNITSON
}
FUNCTION funcCalMin(calnm(mM))() {
    UNITSOFF
    funcCalMin= 1 / (1 + exp((-100000*calnm+(100000*mltypeMin))/1)) 
    UNITSON
}