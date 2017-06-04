// ----------------------------------------------------
//
// 		VOCAL SYNTHESIS MODEL FOR AUDIO FILTERING
//					Graham Davis
//				Music 420 - Spring 2017
//					
// ----------------------------------------------------

import("stdfaust.lib");

// ----------------------------------------------------
// 						UI CONTROLS
// ----------------------------------------------------

itfreq = hslider("freq[unit:Hz]", 200, 100, 300, 1) : si.smoo;
lip = hslider("lip", 0.1, 0.1, 0.99, 0.01) : si.smoo;
nostril = hslider("nostril", 0.99, 0.1, 0.99, 0.01) : si.smoo;
velum = hslider("velum", 0.01, 0.0, 0.4, 0.01) : si.smoo;

// Tube width controls
w1 = vslider("w1[style:knob]", 1.2, 1.2, 2.0, 0.05) : si.smoo;
w2 = vslider("w2[style:knob]", 1.6, 1.2, 3.0, 0.05) : si.smoo;
w3 = vslider("w3[style:knob]", 1.8, 2.0, 3.0, 0.05) : si.smoo;
w4 = vslider("w4[style:knob]", 1.2, 1.2, 3.0, 0.05) : si.smoo;

// ----------------------------------------------------
// 						PROCESS METHOD
// ----------------------------------------------------

process = glottis : fullVocalTract <: _,_;

// ----------------------------------------------------
// 					MODEL IMPLEMENTATIONS
// ----------------------------------------------------

fullVocalTract(x) = pm.endChain(pm.chain(pm.inRightWave(x) : lVocalTract : pm.out)) : lips;

// System including both mouth and nasal cavity
upperSystem(x, y, z) = 	((0, y, z) : uVocalTract),
						((0, 0, z) : nasalCavity) : upperSystemOutput;

upperSystemOutput(x1, y1, z1, x2, y2, z2) = 0, y1 + y2, finalOutput(y1, y2);

//finalOutput(vocal, nasal) = (vocal : lips) + (nasal : nostrils);
finalOutput(vocal, nasal) = (vocal : lips);

// Glottis model 
glottis = os.imptrain(itfreq + 0.2*os.osci(3)) : fi.lowpass(1, lCutoff);	

// Lips model
lips = _ : fi.lowpass(2, (1-lip)*lCutoff) : _;			// Lip Transmission 

// Nostrils model
nostrils = _ : fi.lowpass(2, (1-nostril)*nCutoff) : _; // Nostril Transmission 

// Upper vocal tract
uVocalTract = pm.terminations(*(-.3), pm.chain(
	 section(junction(velum, w2, w3), k(w3, w4)) :
	 lsection(k(w3, w4), 0)), *(lipReflection));

// Nasal cavity
nasalCavity = pm.terminations(*(-.3), pm.chain(
	 section(junction(w3, w2, velum), k(velum, 0.8))), *(lipReflection/2));

// Lower vocal tract
lVocalTract = pm.chain(
	 fsection(glottalReflection, k(w1, w2)) : 
	 section(k(w1, w2), k(w2, w3)) :
	 lsection(k(w2, w3), lipReflection))
;

// ----------------------------------------------------
// 				INSTANCE VARIABLES + METHODS		   
// ----------------------------------------------------
lCutoff = 1000;
nCutoff = 1000;

lipReflection = -0.85;
glottalReflection = 0.7;

nMax = 100;									// Max delay
n = (.01556/pm.speedOfSound) * ma.SR;		// Vocal tube length 
tube = pm.openTube(nMax, n);				// Tube wave guide

lsection(lk, rk) = (_, *(1+lk), _) : 
	pm.terminations(*(-lk), pm.chain( tube ), *(rk));

fsection(lk, rk) = (*(1-rk), _, _) : 
	pm.terminations(*(-lk), pm.chain( tube ), *(rk));

section(lk, rk) = (*(1-rk), *(1+lk), _) : 
	pm.terminations(*(-lk), pm.chain( tube ), *(rk));				

area(d) = ((d/2)^2);					// Compute cross sectional area from diameter
k(d1, d2) = (area(d2) - area(d1)) 
			/ (area(d2) + area(d1));	// Calculate reflection coefficient from diameter (cm)


junction(i, j, k) = 0.3;