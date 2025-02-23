
COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 is very small compared to tau1, this is an alphasynapse with time constant tau2.
If tau1/tau2 is very small, this is single exponential decay with time constant tau2.

The factor is evaluated in the initial block 
such that an event of weight 1 generates a
peak conductance of 1.
supraba(dopa,expectt,weight,ca_nmdai,cali,g)
Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS Gaba_p1
	RANGE tau1, tau2, e, i,w0, wmax, checkactiveD,treshf,treshf_min,mltype,conc0,synapticActiveTime,inRefractoryPeriod, tminrate ,    tmaxrate ,    wrate
	NONSPECIFIC_CURRENT i
    POINTER expectt
	RANGE g,flagx,mltype
	USEION ca_nmda READ ca_nmdai
	USEION cal READ cali VALENCE 2
	USEION ca READ cai VALENCE 2

}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1 = 0.1 (ms) <1e-9,1e9>
	tau2 = 20 (ms) <1e-9,1e9>
	tau   = 10       (ms)
	e=0	(mV)
	w0 = 0.0000188 (uS)
	wmax = 0.006 (uS)
	treshf=0.25(mM)
	treshf_min=0.1(mM)
    h=10
    flagx=0 (1)
    mltype=0 (mM)
    checkactiveD=0 (1)
    f=0   (ms)
    conc0=0 (mM)
    tminrate=0.00001
    tmaxrate=0.00005
    wrate=0.055

}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	expectt
	dopa
	ca_nmdai (mM)
	cali (mM)
	cai (mM)
	synapticActiveTime
	inRefractoryPeriod
	last_time_calic_above_threshold


}

STATE {
	A (uS)
	B (uS)
	weight (uS)
	tresh (mM)
	tresh_min (mM)
	

}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	if (tau1/tau2 < 1e-9) {
		tau1 = tau2*1e-9
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	weight = w0
	tresh = treshf
	tresh_min = treshf_min
	synapticActiveTime = 0
	inRefractoryPeriod=0



	

}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
	weight' =supraba(dopa,expectt,weight,tresh,tresh_min,cali+cai,g)
	tresh' = sTresh(expectt,dopa,tresh,tresh_min,mltype)
	tresh_min' = sTresh_min(expectt,dopa,tresh,tresh_min,mltype)


	

}


NET_RECEIVE(dummy) {
  A = A + weight*factor
  B = B + weight*factor
  if (inRefractoryPeriod == 0) {
    if (flag == 0) { 
      synapticActiveTime = 1  
      inRefractoryPeriod = 1
      COMMENT
        Enter refractory period
      ENDCOMMENT
      net_send(400, -1)
    }
  }
  if (flag == -1) {
    synapticActiveTime = 0
    inRefractoryPeriod = 0
    COMMENT
      Exit refractory period
    ENDCOMMENT
  }
}

FUNCTION supraba(dopam,expec,we(uS),tresh(mM),tresh_min(mM),calic(mM),ge(uS))(uS/ms) {
    UNITSOFF
    
      if (calic > 0.0004) {
           last_time_calic_above_threshold = t
           if (calic > mltype) {
               mltype = calic
           }
       } else if (t - last_time_calic_above_threshold >= 300) {
           mltype = 0
       }
       
   
	supraba =wrate*expec*FuncCalj(mltype,dopam,ge,tresh,tresh_min)*we*(wmax-we)
	UNITSON    
}

FUNCTION sTresh(expec,dopam,tresh(mM),tresh_min(mM),calnm)(uM/ms) {
    UNITSOFF
    sTresh =tmaxrate*expec*(-1/(1+exp(-(h*Norm(calnm)-h*tresh_min)/1))+3/(1+exp(-(h*Norm(calnm)-h*tresh)/1)))
	UNITSON    
}
FUNCTION sTresh_min(expec,dopam,tresh(mM),tresh_min(mM),calnm)(uM/ms) {
    UNITSOFF
    sTresh_min =-tminrate*expec*(-2/(1+exp(-(h*Norm(calnm+0.0006)-h*tresh_min)/1))+3/(1+exp(-(h*Norm(calnm+0.0006)-h*tresh)/1)))
	UNITSON    
}
FUNCTION FuncCalj(calnm(mM),dopam,g(uS),tresh (mM),tresh_min(mM))() {
    UNITSOFF
    if (synapticActiveTime==0) {
	  FuncCalj =1*(-1/(1+exp(-(h*Norm(calnm)-h*tresh_min)/1))+3/(1+exp(-(h*Norm(calnm)-h*tresh)/1)))
	  

	} else {
	  FuncCalj =-1*(-1/(1+exp(-(h*Norm(calnm)-h*tresh_min)/1))+3/(1+exp(-(h*Norm(calnm)-h*tresh)/1)))

      }
    UNITSON
}

FUNCTION Norm(calnm(mM))() {
    UNITSOFF
    Norm=250*calnm
    UNITSON
}
