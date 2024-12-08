clear all; close all;

%Ezmx = waveamp(Ydim,Xdim,isrc,jsrc,dx,dy,lambda);



%  Ezmx = waveamp(100,100,-20,30,50.0E-9,50.0E-9,500.0E-9);
%  Ezmx = Ezmx + waveamp(100,100,-20,60,50.0E-9,50.0E-9,500.0E-9);


Ezmx(1:300,1:300) = 0.0;

lambda = 300.0E-9;

%ipoints = [150:160 170:180];
%ipoints = [131 131] ;
ipoints = [150:170];
%ipoints = [130 170];

for iplane = ipoints 
    
   
    ampsrc=1;
    Ezmx = Ezmx + waveamp(300,300,0,iplane,ampsrc,100.0E-9,100.0E-9,lambda);
end

% for iplane = 180:210
%     Ezmx = Ezmx + waveamp(300,300,-50,iplane,50.0E-9,50.0E-9,500.0E-9);
% end


% plot(squeeze(abs(Ezmx(20,:))));
%
% plot(squeeze(real(Ezmx(20,:))));


figure; hold all;
plot(squeeze(abs(Ezmx(10,:))));

plot(squeeze(abs(Ezmx(50,:))));
%plot(squeeze(abs(Ezmx(100,:))));
%plot(squeeze(abs(Ezmx(150,:))));
%plot(squeeze(abs(Ezmx(200,:))));
plot(squeeze(abs(Ezmx(250,:))));

shading interp;
title('Absolute value Ez along a lines');


figure;
pcolor(abs(Ezmx));
shading interp; axis equal;
title('Absolute value Ez ');


figure;
pcolor(real(Ezmx));
shading interp; axis equal;
title('Real value of Mz');
clim([-150 150])


figure;
surf(real(Ezmx));
shading interp; 
title('Real value of Mz');