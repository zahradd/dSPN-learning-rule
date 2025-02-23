TITLE GABA_A synapse with short-term plasticity


NEURON {
    POINT_PROCESS tmGabaA_Plasticity_tresh
    RANGE tau1, tau2, e, i, q
    RANGE tau, tauR, tauF, U, u0,g, mltype,tresh,tresh_min,g_spike,fff,g_active_initial,wmax,mltype_tau,mltype_max,treshf,treshf_min,g_active_tresh, weight_rate
    RANGE failRate, damod, maxMod, w0
    USEION cal READ cali VALENCE 2
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    : q = 2, now included in tau1,tau2 parameters.     
    tau1= 0.25 (ms) : ORIG: 0.5ms
    tau2 = 3.75 (ms)  : ORIG: 7.5ms, tau2 > tau1
    e = -65 (mV)
    tau = 3 (ms)
    tauR = 500 (ms)  : tauR > tau
    tauF = 0 (ms)    : tauF >= 0
    U = 0.1 (1) <0, 1>
    u0 = 0 (1) <0, 1>
    failRate = 0	
    damod = 0
    maxMod = 1
    h = 1000
    w0 = 0  (uS)
    treshf = 0.001(mM)
	treshf_min = 0.00001(mM)
    g_active_initial=0
    wmax= 0.005
    mltype_max=15
    mltype_tau=100
    g_active_tresh=0.1
    weight_rate=0.0055

     

}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
    factor
    x
    cali (mM)
    fff

}

STATE {
    A (uS)
    B (uS)
   	weight_ltp (uS)
   	tresh (mM)
   	tresh_min (mM)
   	mltype
    g_spike
    
}

INITIAL {
    LOCAL tp
    A = 0
    B = 0
    tp = (tau1*tau2)/(tau2-tau1) * log(tau2/tau1)
    factor = -exp(-tp/tau1) + exp(-tp/tau2)
    factor = 1/factor
   	weight_ltp = w0
   	tresh = treshf
   	tresh_min = treshf_min
   	mltype=0 
    g_spike=0

}

BREAKPOINT {
    SOLVE state METHOD cnexp
    g = (B - A) *modulation() * weight_ltp
    i = g*(v - e)
}

DERIVATIVE state {
    A' = -A/tau1
    B' = -B/tau2

   	mltype'=(mltype_max*cali-mltype)/mltype_tau
   	weight_ltp' = supraba(weight_ltp,tresh,tresh_min,mltype,g)
   	tresh' = sTresh(tresh,tresh_min,mltype)
   	tresh_min' = sTresh_min(tresh,tresh_min,mltype)
    g_spike'=-g_spike/100
}

NET_RECEIVE(weight (uS), y, z, u, tsyn (ms)) {
    LOCAL result
    INITIAL {
        y = 0
        z = 0
        u = u0
        tsyn = t
    }
    if ( weight <= 0 ) {
VERBATIM
        return;
ENDVERBATIM
    }
    if( urand() > failRate ) { 
      z = z*exp(-(t-tsyn)/tauR)
      z = z + (y*(exp(-(t-tsyn)/tau) - exp(-(t-tsyn)/tauR)) / (tau/tauR - 1) )
      y = y*exp(-(t-tsyn)/tau)
      x = 1-y-z
      if (tauF > 0) {
          u = u*exp(-(t-tsyn)/tauF)
          u = u + U*(1-u)
      } else {
          u = U
      }
    A = A + weight*factor*x*u / U
    B = B + weight*factor*x*u / U
    y = y + x*u
    g_spike=g_spike + 0.1
    tsyn = t
    }
}

FUNCTION urand() {
    urand = scop_random(1)
}


FUNCTION modulation() {
    : returns modulation factor
    
    modulation = 1 + damod*(maxMod-1)
}
FUNCTION supraba(we(uS),tresh(mM),tresh_min(mM),mltype(mM),ge(uS))(uS/ms) {
    UNITSOFF   
	supraba =weight_rate*FuncCalj(mltype,ge,tresh,tresh_min)*(we)*(wmax-we)
	UNITSON    
}

FUNCTION sTresh(tresh(mM),tresh_min(mM),calnm)(uM/ms) {
    UNITSOFF
    sTresh =0.000000025*(-1/(1+exp(-(h*(calnm)-h*tresh_min)/1))+3/(1+exp(-(h*(calnm)-h*tresh)/1)))
	UNITSON    
}
FUNCTION sTresh_min(tresh(mM),tresh_min(mM),calnm)(uM/ms) {
    UNITSOFF
    sTresh_min =-0.000000008*(-2/(1+exp(-(h*(calnm+0.014)-h*tresh_min)/1))+3/(1+exp(-(h*(calnm+0.014)-h*tresh)/1)))
	UNITSON    
}
FUNCTION FuncCalj(calnm(mM),g(uS),tresh (mM),tresh_min(mM))() {
    UNITSOFF
    if (g_spike <g_active_tresh) {
          fff=-1
    	  FuncCalj =1*(-1/(1+exp(-(h*(calnm)-h*tresh_min)/1))+3/(1+exp(-(h*(calnm)-h*tresh)/1)))
	} else {
          fff=1
	      FuncCalj =-1*(-1/(1+exp(-(h*(calnm)-h*tresh_min)/1))+3/(1+exp(-(h*(calnm)-h*tresh)/1)))
      }
    UNITSON
}
