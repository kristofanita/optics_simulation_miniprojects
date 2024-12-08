% Gaussian Slit Parameters
w = 10e-6; % slit width in meters (standard deviation of Gaussian)
lambda = 630e-9; % wavelength in meters
d = 200e-6; % distance in meters

% Grid specifications
dx = 0.1e-6; % grid spacing in x-direction
dy = 0.1e-6; % grid spacing in y-direction
Xdim = 300;
Ydim = 300;

% Initialize the field
Ezmx = zeros(Ydim, Xdim);

% Gaussian slit
y_center = Ydim/2;
ampsrc = 1;
for ysrc = 1:Ydim
    amp = ampsrc * exp(-((ysrc - y_center)*dy)^2 / (2*w^2));
    Ezmx = Ezmx + waveamp(Ydim, Xdim, 0, ysrc, amp, dx, dy, lambda);
end

% Visualization of diffraction pattern
figure; hold all;
plot(squeeze(abs(Ezmx(10,:))));
plot(squeeze(abs(Ezmx(50,:))));
plot(squeeze(abs(Ezmx(250,:))));
shading interp;
title('Absolute value Ez along different lines');

figure;
pcolor(abs(Ezmx));
shading interp; axis equal;
title('Absolute value Ez');

figure;
pcolor(real(Ezmx));
shading interp; axis equal;
title('Real value of Ez');
caxis([-150 150]);

figure;
surf(real(Ezmx));
shading interp; 
title('Real value of Ez');
