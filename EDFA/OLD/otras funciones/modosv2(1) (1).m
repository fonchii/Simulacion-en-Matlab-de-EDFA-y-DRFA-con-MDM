function modos = modosv2(radio,n1,n2,lambda)
%%LP modes , Calculate and draw
%Programmed by Ammar Rajab.
%Tishreen University , Latakia , Syria.

% % fiber attributes
if nargin<4
    clc
     radio=25e-6;    %diameter  
     n1=1.46;    %core  
     n2=1.455;   %cladding    
     lambda=1550e-9;
    % fprintf('Modo por defecto: \na=25 μm ; n1=1.46 ; n2 = 1.455 ; λ=1550nm')
end
%a = 2*radio; 
a = radio;
v=(pi*a/lambda)*sqrt(n1^2-n2^2);    %normalized freq

%options = optimset('Display','off');

% % Calculate the modes

m0 = AllZeros(@(x) besselj(0,x), 1, v);
m1 = AllZeros(@(x) besselj(1,x), 1, v);
m2 = AllZeros(@(x) besselj(2,x), 1, v);
m3 = AllZeros(@(x) besselj(3,x), 1, v);
m4 = AllZeros(@(x) besselj(4,x), 1, v);
m5 = AllZeros(@(x) besselj(5,x), 1, v);
m6 = AllZeros(@(x) besselj(6,x), 1, v);
m7 = AllZeros(@(x) besselj(6,x), 1, v);
m8 = AllZeros(@(x) besselj(6,x), 1, v);
m9 = AllZeros(@(x) besselj(6,x), 1, v);
m10 = AllZeros(@(x) besselj(6,x), 1, v);

nmodos = length(m0) + length(m1) + length(m2) + length(m3) + length(m4) + length(m5)...
    + length(m6) + length(m7)+ length(m8)+ length(m9)+ length(m10); 
nmodos_t = length(m0) + 2*length(m1) + 2*length(m2) + 2*length(m3) + 2*length(m4) + 2*length(m5)...
     + 2*length(m6) + 2*length(m7) + 2*length(m8) + 2*length(m9) + 2*length(m10);
 
m = [m0 m1 m2 m3 m4 m5 m6 m7 m8 m9 m10];
msort = sort(m);
while(m(1)<1)
    m = m(2,end)
end
msort_lp = []; msort_LP = [];
% Ordenar modos:
for x=1:length(msort)
    for xx=1:length(m0)
        if m0(xx)==msort(x)
            msort_lp = [msort_lp strcat(' 0',int2str(xx))];
            msort_LP = [msort_LP strcat(' LP0',int2str(xx))];
        end
    end
    for xx=1:length(m1)
        if m1(xx)==msort(x)
            msort_lp = [msort_lp strcat(' 1',int2str(xx))];
            msort_LP = [msort_LP strcat(' LP1',int2str(xx))];
        end
    end
    for xx=1:length(m2)
        if m2(xx)==msort(x)
            msort_lp = [msort_lp strcat(' 2',int2str(xx))];
            msort_LP = [msort_LP strcat(' LP2',int2str(xx))];
        end
    end
    for xx=1:length(m3)
        if m3(xx)==msort(x)
            msort_lp = [msort_lp strcat(' 3',int2str(xx))];
            msort_LP = [msort_LP strcat(' LP3',int2str(xx))];
        end
    end
    for xx=1:length(m4)
        if m4(xx)==msort(x)
            msort_lp = [msort_lp strcat(' 4',int2str(xx))];
            msort_LP = [msort_LP strcat(' LP4',int2str(xx))];
        end
    end
    for xx=1:length(m5)
        if m5(xx)==msort(x)
            msort_lp = [msort_lp strcat(' 5',int2str(xx))];
            msort_LP = [msort_LP strcat(' LP5',int2str(xx))];
        end
    end
    for xx=1:length(m6)
        if m6(xx)==msort(x)
            msort_lp = [msort_lp strcat(' 6',int2str(xx))];
            msort_LP = [msort_LP strcat(' LP6',int2str(xx))];
        end
    end
end

msort_lp = strsplit(msort_lp) ; msort_lp = string(msort_lp(2:end));
msort_LP = strsplit(msort_LP) ; msort_LP = string(msort_LP(2:end));
modos = struct;
modos.freqs = m; modos.freqs_sort = msort; modos.LPsort=msort_lp;
modos.LPPsort=msort_LP;
modos.nmodos=nmodos; modos.nmodos_t=nmodos_t;
modos.m0 = m0; modos.m1 = m1; modos.m2 = m2; modos.m3 = m3;
modos.m4 = m4; modos.m5 = m5; modos.m6 = m6; modos.m7 = m7;
modos.m8 = m8; modos.m9 = m9; modos.m10 = m10;
modos.v = v;