function [ase_wl, ase_freq] = ase_freqVPI(f1,fend,deltaF,c,flag)  
if nargin<5
    ase_freq=f1:deltaF:fend+deltaF;
    ase_wl = c./ase_freq;
else 
    ase_wl = [0,0,0];

end