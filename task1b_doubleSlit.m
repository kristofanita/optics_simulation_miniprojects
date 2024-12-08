% Double Slit Parameters
w = 10e-6; % slit width in meters
distance_between_slits = 20e-6; % distance between slits in meters
lambda = 630e-9; % wavelength in meters
d = 200e-6; % distance in meters

% Grid specifications
dx = 0.1e-6; % grid spacing in x-direction
dy = 0.1e-6; % grid spacing in y-direction
Xdim = 300;
Ydim = 300;

% Initialize the field
Ezmx = zeros(Ydim, Xdim);

% Position of the slits along y-axis
y_slit_start1 = round((Ydim/2) - (distance_between_slits/(2*dy)) - (w/(2*dy)));
y_slit_end1 = round(y_slit_start1 + w/dy);

y_slit_start2 = round((Ydim/2) + (distance_between_slits/(2*dy)) - (w/(2*dy)));
y_slit_end2 = round(y_slit_start2 + w/dy);

% Add wave amplitudes for each point on the slits
ampsrc = 1;
for ysrc = y_slit_start1:y_slit_end1
    Ezmx = Ezmx + waveamp(Ydim, Xdim, 0, ysrc, ampsrc, dx, dy, lambda);
end

for ysrc = y_slit_start2:y_slit_end2
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
