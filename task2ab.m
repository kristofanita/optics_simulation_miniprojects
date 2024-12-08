% Parameters for the diffraction grating
lambda = 500e-9; % wavelength in meters
d = 5e-6; % distance between slits in meters (grating period)
number_of_slits = 5; % number of slits
distance_to_screen = 2000e-6; % distance to the screen in meters

% Grid specifications
dx = 0.1e-6; % grid spacing in x-direction
dy = 0.1e-6; % grid spacing in y-direction
Xdim = 800;
Ydim = 800;

% Initialize the field
Ezmx = zeros(Ydim, Xdim);

% Position of the slits along y-axis
ampsrc = 1;
for n = 0:(number_of_slits-1)
    slit_position = round((Ydim/2) + n*d/dy);
    Ezmx = Ezmx + waveamp(Ydim, Xdim, 0, slit_position, ampsrc, dx, dy, lambda);
end
% A)
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



% Theoretical diffraction pattern calculation
theta = linspace(-pi/2, pi/2, Xdim);
d_theta = distance_to_screen / (dx * Xdim); % angular resolution
I_theoretical = (sin(number_of_slits * pi * d * sin(theta) / lambda) ./ ...
                 sin(pi * d * sin(theta) / lambda)).^2;

% Plot comparison
figure;
subplot(2,1,1);
plot(theta, I_theoretical);
title('Theoretical Diffraction Pattern of Grating');
subplot(2,1,2);
plot(abs(Ezmx(round(Ydim/2), :)));
title('Simulated Diffraction Pattern of Grating');


%%
%B)
% Parameters for incoherent light simulation
coherence_length = 5e-6; % coherence length in meters

% Initialize the field
Ezmx = zeros(Ydim, Xdim);

% Random phase shifts
for n = 0:(number_of_slits-1)
    slit_position = round((Ydim/2) + n*d/dy);
    random_phase = exp(1i * 2 * pi * rand);
    if mod(n, 2) == 0
        coherence_length = 5e-6;
    else
        coherence_length = 2e-6;
    end
    Ezmx = Ezmx + waveamp(Ydim, Xdim, 0, slit_position, ampsrc * random_phase, dx, dy, lambda);
end

% Visualization of diffraction pattern with incoherent light
figure; hold all;
plot(squeeze(abs(Ezmx(10,:))));
plot(squeeze(abs(Ezmx(50,:))));
plot(squeeze(abs(Ezmx(250,:))));
shading interp;
title('Absolute value Ez along different lines with incoherent light');

figure;
pcolor(abs(Ezmx));
shading interp; axis equal;
title('Absolute value Ez with incoherent light');

figure;
pcolor(real(Ezmx));
shading interp; axis equal;
title('Real value of Ez with incoherent light');
caxis([-150 150]);

figure;
surf(real(Ezmx));
shading interp; 
title('Real value of Ez with incoherent light');
