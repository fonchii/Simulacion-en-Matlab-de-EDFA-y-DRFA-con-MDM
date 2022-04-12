% CROSS SECTION COMPARATIVE

VPI = load('Erbium_VPI.dat');
OS = load('Erbium.dat');

figure(1)
plot(VPI(:,1),VPI(:,2),'DisplayName','VPI');hold on;plot(OS(:,1),OS(:,3),'DisplayName','OptiSystem') ; legend()
title('Emission Cross Section')
figure(2)
plot(VPI(:,1),VPI(:,3),'DisplayName','VPI');hold on;plot(OS(:,1),OS(:,2),'DisplayName','OptiSystem');legend()
title('Absorption Cross Section')
