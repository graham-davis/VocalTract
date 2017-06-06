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

intensity = hslider("intensity", 0.5, 0.0, 1.0, 0.05) : si.smoo;
itfreq = hslider("freq[unit:Hz]", 200, 100, 300, 1) : si.smoo;
voicing = hslider("voicing", 0.7, 0.0, 1.0, 0.01) : si.smoo;
lip = hslider("lip", 0.01, 0.01, 0.8, 0.01) : si.smoo;
velum = hslider("velum", 0.01, 0.0, 0.4, 0.01) : si.smoo;
th = hslider("tongue height", 0.0, 0.0, 0.7, 0.01) : si.smoo;
tc = hslider("tongue center", 8.0, 2.0, 15.0, 1.0) : si.smoo;
lpfc = hslider("cutoff", 0.01, 0.01, 0.95, 0.01) : si.smoo;


w1 = 0.6;
w2 = 0.6;
w3 = 0.6;
w4 = 0.6;
w5 = 1.1;
// TONGUE START
w6 = 1.50 	- tongue(1);
w7 = 1.45 	- tongue(2);
w8 = 1.39 	- tongue(3);
w9 = 1.34 	- tongue(4);
w10 = 1.32 	- tongue(5);
w11 = 1.3 	- tongue(6);
w12 = 1.3 	- tongue(7);
w13 = 1.32	- tongue(8);
w14 = 1.35	- tongue(9);
w15 = 1.39	- tongue(10);
w16 = 1.45	- tongue(11);
w17 = 1.50 	- tongue(12);
w18 = 1.56 	- tongue(13);
w19 = 1.61 	- tongue(14);
w20 = 1.64	- tongue(15);
w21 = 1.65	- tongue(16);
// TONGUE END
w22 = 1.5;
w23 = 1.5;
w24 = 1.5;

// ----------------------------------------------------
// 						PROCESS METHOD
// ----------------------------------------------------

// Glotal input
process = intensity*glottis : fullVocalTract <: _,_;

// VST input
//process = fullVocalTract <: _,_;

// ----------------------------------------------------
// 					MODEL IMPLEMENTATIONS
// ----------------------------------------------------

fullVocalTract(x) = pm.endChain(pm.chain(pm.inRightWave(x) : lVocalTract : lipLowpass : pm.out)) : lips : lpf;

// Glottis model 
glottis = (voicing*voiced) + ((1-voicing) * aspiration) 
with {
	voiced = os.imptrain(freq) : fi.peak_eq(3.5, freq, 50.0) : fi.lowpass(1, voicing*gCutoff + 10);
	aspiration = ((0.05*os.osci(freq*2))+0.1*no.gnoise(20)) : fi.lowpass(1, 10000) : fi.highpass(2, 4000);
	freq = itfreq + 0.05*os.osci(3);
};	

lpf = _ : fi.lowpass(2, (1-lpfc)*lCutoff) : _;

// Lips model
lips = _ : fi.highpass(2, lip*lCutoff) : _;			// Lip Transmission 

lipLowpass = fi.lowpass(2, (1-lip)*lCutoff), _, _;

// Lower vocal tract
lVocalTract = pm.chain(
	 lsection(glottalReflection, k(w1, w2)) : 
	 section(k(w1, w2), k(w2, w3)) :
	 section(k(w2, w3), k(w3, w4)) :
	 section(k(w3, w4), k(w4, w5)) :
	 section(k(w4, w5), k(w5, w6)) :
	 section(k(w5, w6), k(w6, w7)) :
	 section(k(w6, w7), k(w7, w8)) :
	 section(k(w7, w8), k(w8, w9)) :
	 section(k(w8, w9), k(w9, w10)) :
	 section(k(w9, w10), k(w10, w11)) :
	 section(k(w10, w11), k(w11, w12)) :
	 section(k(w11, w12), k(w12, w13)) :
	 section(k(w12, w13), k(w13, w14)) :
	 section(k(w13, w14), k(w14, w15)) :
	 section(k(w14, w15), k(w15, w16)) :
	 section(k(w15, w16), k(w16, w17)) :
 	 section(k(w16, w17), k(w17, w18)) :
	 section(k(w17, w18), k(w18, w19)) :
	 section(k(w18, w19), k(w19, w20)) :
	 section(k(w19, w20), k(w20, w21)) :
	 section(k(w20, w21), k(w21, w22)) :
	 section(k(w21, w22), k(w22, w23)) :
 	 section(k(w22, w23), k(w23, w24)) :
	 rsection(k(w23, w24), lipReflection));

// ----------------------------------------------------
// 				INSTANCE VARIABLES + METHODS		   
// ----------------------------------------------------
gCutoff = 2500;								// Glottal cutoff frequency
lCutoff = 2500;								// Lip cutoff frequency

tractLength = 0.169;						// Vocal tract length in meters
nSections 	= 24;							// Number of sections in model
nTSections 	= 16;							// Number of tongue sections in model
n = pm.l2s(tractLength/nSections);			// Tube length per tract section
tube = pm.openTube(n*2, n);					// Tube wave guide

lipReflection = -0.85;						// Reflection coefficient of lips
glottalReflection = 0.7;					// Reflection coefficient of glottus

// Initial KL junction
lsection(lk, rk) = (*(1-rk), _, _) : 
	pm.terminations(*(-lk), pm.chain( tube ), *(rk));

// Terminating KL junction 
rsection(lk, rk) = (_, *(1+lk), _) : 
	pm.terminations(*(-lk), pm.chain( pm.inRightWave(fricative) : tube ), *(rk))
	with{
		fricative = (lip > 0.6) * (((lip-0.6)/0.3)*0.2*no.gnoise(20)) : fi.bandpass(2, 700, 1600);
	};

// Middle KL junction
section(lk, rk) = (*(1-rk), *(1+lk), _) : 
	pm.terminations(*(-lk), pm.chain( tube ), *(rk));				

// Compute cross sectional area from diameter (disregarding PI)
area(d) = ((d/2)^2);	

// Compute scattering coefficient at junction of two diameters (d1 on left, d2 on right)				
k(d1, d2) = (area(d1) - area(d2)) 
			/ (area(d1) + area(d2));

// Calculate tongue effect on specific tube segment
tongue(index) = ma.fmax(ma.fmin(-trueIndex^2 + th, 1.0), -1.0)
with {
	trueIndex = (1/nTSections)*(index - tc);
};
