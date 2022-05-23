function [ase_wl, ase_freq] = ase_freqVPI(c,flag)  
if nargin<2
    ase_freq = flip(189.25e12:300e9:198.85e12);
    ase_wl = c./ase_freq;
else 
    ase_wl = [0,0,0];
    ase_freq = [];
end