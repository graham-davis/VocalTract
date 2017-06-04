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

noiseLevel = hslider("noise", 0.05, 0.01, 0.2, 0.01) : si.smoo;
itfreq = hslider("freq[unit:Hz]", 200, 100, 300, 1) : si.smoo;
lip = hslider("lip", 0.5, 0.1, 0.99, 0.01) : si.smoo;
nostril = hslider("nostril", 0.99, 0.1, 0.99, 0.01) : si.smoo;
velum = hslider("velum", 0.01, 0.01, 0.4, 0.01) : si.smoo;

// Tube width controls
w1 = vslider("w1[style:knob]", 10, 2, 36, 1) : si.smoo;
w2 = vslider("w2[style:knob]", 16, 2, 36, 1) : si.smoo;
w3 = vslider("w3[style:knob]", 24, 2, 36, 1) : si.smoo;
w4 = vslider("w4[style:knob]", 14, 2, 36, 1) : si.smoo;


// ----------------------------------------------------
// 						PROCESS METHOD
// ----------------------------------------------------

process = glottis : fullVocalTract <: _,_;

// ----------------------------------------------------
// 					MODEL IMPLEMENTATIONS
// ----------------------------------------------------

fullVocalTract(x) = pm.endChain(pm.chain(pm.inRightWave(x) : lVocalTract : upperSystem : pm.out));

upperSystem(x, y, z) = 	((0, (1-velum)*y, z) : uVocalTract),
						((0, velum*y, z) : nasalCavity) : upperSystemOutput;

upperSystemOutput(x1, y1, z1, x2, y2, z2) = 0, y1 + y2, finalOutput(y1, y2);

finalOutput(vocal, nasal) = (vocal : lips) + (nasal : nostrils);

// Glottis model 
glottis = (0*no.noise) + ((1-noiseLevel)*os.imptrain(itfreq + 0.05*os.osci(3))) : fi.lowpass(1, lCutoff);	

// Lips model
lips = _ : fi.highpass(2, (1-lip)*lCutoff) : _;			// Lip Transmission 

// Nostrils model
nostrils = _ : fi.highpass(2, (1-nostril)*nCutoff) : _; // Nostril Transmission 

// Upper vocal tract
uVocalTract = pm.terminations(0, pm.chain(
	 tube) , rt)
with {
	rt = fi.lowpass(2, (1-lip)*lCutoff);				// Lip reflection
};

// Nasal cavity
nasalCavity = pm.terminations(0, pm.chain(
		tube), 0);

// Lower vocal tract
lVocalTract = pm.lTermination(lt, pm.chain(
	 section(0, k(w1, w2)) : 
	 section(k(w1, w2), k(w2, w3)) : 
	 section(k(w2, w3), k(w3, w4))))
with {
	lt = *(0.7);										// Glottal reflection coefficient
};


// ----------------------------------------------------
// 				INSTANCE VARIABLES + METHODS
// ----------------------------------------------------
nSegments = 6;

lCutoff = 3500;
nCutoff = 1000;

nMax = 100;									// Max delay
n = (.01556/pm.speedOfSound) * ma.SR;		// Vocal tube length 
tube = pm.openTube(nMax, n);				// Tube wave guide
section(lk, rk) = pm.chain(					// full tube section
	pm.terminations(*(-lk), tube, 
	*(rk)), *(1-rk), *(1+lk),_);				

area(d) = ((d/2)^2);					// Compute cross sectional area from diameter
k(d1, d2) = (area(d2) - area(d1)) 
			/ (area(d2) + area(d1));	// Calculate reflection coefficient from diameter (cm)