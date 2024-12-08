% Parameters
lambda = 500e-9; % wavelength in meters
k = 2 * pi / lambda; % wavenumber
Xdim = 1000;
Ydim = 1000;
dx = 1e-6; % 1 micron spacing
dy = 1e-6; % 1 micron spacing

% Define the positions of the point sources
source1_position = [Ydim/2, Xdim/3];
source2_position = [Ydim/2, 2*Xdim/3];

% Calculate the wave propagation from each source
Ez_source1 = waveamp(Ydim, Xdim, source1_position(1), source1_position(2), 1, dx, dy, lambda);
Ez_source2 = waveamp(Ydim, Xdim, source2_position(1), source2_position(2), 1, dx, dy, lambda);

% Generate the interference pattern (hologram)
Ez_interference = Ez_source1 + Ez_source2;

% Record the interference pattern along a line (e.g., the central row)
hologram = abs(Ez_interference(Ydim/2, :)).^2;

% Visualization
figure;
plot(hologram);
title('Recorded Hologram (Interference Pattern)');
xlabel('Position along the line');
ylabel('Intensity');


%%%%%%%%%%%%%%%%%%%%%%
%B)
% Parameters for reconstruction
illumination_position = source1_position; % Use source1 to illuminate the hologram

% Initialize the field with the illuminating point source
Ez_illumination = waveamp(Ydim, Xdim, illumination_position(1), illumination_position(2), 1, dx, dy, lambda);

% Apply the hologram (recorded interference pattern) as a modulation
Ez_reconstructed = Ez_illumination;
Ez_reconstructed(Ydim/2, :) = Ez_reconstructed(Ydim/2, :) .* hologram;

% Visualization
figure;
imagesc(abs(Ez_reconstructed));
colorbar;
title('Reconstructed Wavefront');
xlabel('x (microns)');
ylabel('y (microns)');