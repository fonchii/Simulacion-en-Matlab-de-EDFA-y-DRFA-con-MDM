function modos = modos(radio,n1,n2,lambda)
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
     fprintf('Modo por defecto: \na=25 μm ; n1=1.46 ; n2 = 1.455 ; λ=1550nm')
end
a = 2*radio; 

v=(pi*a/lambda)*sqrt(n1^2-n2^2);    %normalized freq

% % Calculate the modes

m0=[];m1=[];m2=[];m3=[];m4=[];m5=[];m6=[];

for x=1:10
    m0(end+1)=fzero(@(x) besselj(0,x),[x-1 x]*pi);
    if x*pi>v
        break;
    end
end

for x=1:10
    t=fzero(@(x)besselj(1,x),[x-1 x]*pi);
    if t>0
        m1(end+1)=t;
    end
    if x*pi>v
        break;
    end
end

for x=1:10
    t=fzero(@(x)besselj(2,x),[x-1 x]*pi);
    if t>0
        m2(end+1)=t;
    end
    if x*pi>v
        break;
    end
end

for x=3:10
    t=fzero(@(x)besselj(3,x),[x-1 x]*pi);
    if t>0
        m3(end+1)=t;
    end
    if x*pi>v
        break;
    end
end

for x=3:10
    t=fzero(@(x)besselj(4,x),[x-1 x]*pi);
    if t>0
        m4(end+1)=t;
    end
    if x*pi>v
        break;
    end
end

for x=3:10
    t=fzero(@(x)besselj(5,x),[x-1 x]*pi);
    if t>0
        m5(end+1)=t;
    end
    if x*pi>v
        break;
    end
end

for x=4:10
    t=fzero(@(x)besselj(6,x),[x-1 x]*pi);
    if t>0
        m6(end+1)=t;
    end
    if x*pi>v
        break;
    end
end

nmodos = length(m0) + length(m1) + length(m2); 
nmodos_t = length(m0) + 2*length(m1) + 2*length(m2);
m = [m0 m1 m2 m3 m4 m5 m6];
msort = sort(m);
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
modos.m4 = m4; modos.m5 = m5; modos.m6 = m6;
modos.v = v;


%% GRAFICOS


% %%mesh draw points r & phi
% r=linspace(0,a,400);
% phi=linspace(0,2*pi,400);
% [r,phi]=meshgrid(r,phi);

%%draw modes for m=0

% v = 1;h=3; 
% while v*h<length(m0)
%     v=v+1;end
% 
% f1=figure();f1.WindowState = 'maximized';
% for i=1:length(m0)
%     z=besselj(0,m0(i)*r/a).*cos(0.*phi);
%     x=r.*cos(phi);
%     y=r.*sin(phi);
%     subplot(v,h,i)
%     mesh(x,y,z.^2);view(2);
%     title(strcat('LP 0',int2str(i)));
% end
% 
% %%draw modes for m=1
% for i=1:length(m1)
%     z_a=besselj(1,m1(i)*r/a).*cos(1.*phi);
%     z_b=besselj(1,m1(i)*r/a).*sin(1.*phi);
%     x=r.*cos(phi);
%     y=r.*sin(phi);
%     figure;% WindowState='maximized';
%     subplot 121 ; mesh(x,y,z_a.^2);title(strcat(strcat('LP 1',int2str(i)),'a'));view(2);
%     subplot 122 ; mesh(x,y,z_b.^2);title(strcat(strcat('LP 1',int2str(i)),'b'));view(2);
% end
% 
% %%draw modes for m=2
% 
% for i=1:length(m2)
%     z_a=besselj(2,m2(i)*r/a).*cos(2.*phi);
%     z_b=besselj(2,m2(i)*r/a).*sin(2.*phi);
%     x=r.*cos(phi);
%     y=r.*sin(phi);
%     figure; %WindowState='maximized';
%     subplot 121 ; mesh(x,y,z_a.^2);title(strcat(strcat('LP 2',int2str(i)),'a'));view(2);
%     subplot 122 ; mesh(x,y,z_b.^2);title(strcat(strcat('LP
%     2',int2str(i)),'b'));view(2);
% end
