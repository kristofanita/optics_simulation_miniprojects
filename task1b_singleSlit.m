% Single Slit Parameters
w = 10e-6; % slit width in meters
lambda = 630e-9; % wavelength in meters
d = 200e-6; % distance in meters

% Grid specifications
dx = 0.1e-6; % grid spacing in x-direction
dy = 0.1e-6; % grid spacing in y-direction
Xdim = 300;
Ydim = 300;

% Initialize the field
Ezmx = zeros(Ydim, Xdim);

% Position of the slit along y-axis
y_slit_start = round((Ydim/2) - (w/(2*dy)));
y_slit_end = round((Ydim/2) + (w/(2*dy)));

% Add wave amplitudes for each point on the slit
ampsrc = 1;
for ysrc = y_slit_start:y_slit_end
    Ezmx = Ezmx + waveamp(Ydim, Xdim, 0, ysrc, ampsrc, dx, dy, lambda);
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
