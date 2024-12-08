% Parameters
lambda = 500e-9; % wavelength in meters
k = 2 * pi / lambda; % wavenumber
f = 0.01; % focal length of the lens in meters (10 mm)
Xdim = 1000;
Ydim = 1000;
dx = 1e-6; % 1 micron spacing
dy = 1e-6; % 1 micron spacing

% Generate the phase shift introduced by the lens
[X, Y] = meshgrid((1:Xdim) * dx, (1:Ydim) * dy);
lens_center = [Xdim / 2, Ydim / 2];
x_shift = (X - lens_center(1)) * dx;
y_shift = (Y - lens_center(2)) * dy;
phase_shift = exp(1i * k * (x_shift.^2 + y_shift.^2) / (2 * f));

% Parallel light beams represented as plane waves
Ezmx = ones(Ydim, Xdim);

% Apply the lens phase shift
Ezmx = Ezmx .* phase_shift;

% Propagate the wave to the focal point
Ezmx_focused = waveamp(Ydim, Xdim, Ydim/2, Xdim/2, 1, dx, dy, lambda);

% Visualization
figure;
imagesc(abs(Ezmx_focused));
colorbar;
title('Focal Point Intensity');
xlabel('x (microns)');
ylabel('y (microns)');

%%%%%%%%%%%%%%%%
%B)
% Parameters
source_distance = 0.02; % distance of the source from the lens in meters
image_distance = 0.02; % distance of the image from the lens in meters

% Place the source
source_position = [Ydim/2, Xdim/2];

% Calculate the wave propagation to the lens
Ezmx_to_lens = waveamp(Ydim, Xdim, source_position(1), source_position(2), 1, dx, dy, lambda);

% Apply the lens phase shift
Ezmx_lens = Ezmx_to_lens .* phase_shift;

% Propagate the wave beyond the lens
Ezmx_image = waveamp(Ydim, Xdim, Ydim/2, Xdim/2, 1, dx, dy, lambda);

% Visualization
figure;
imagesc(abs(Ezmx_image));
colorbar;
title('Image of Point Source');
xlabel('x (microns)');
ylabel('y (microns)');


%%%%%%%%%%%%%%%%%%
%C)
% Parameters for two point sources
source_distance_1 = [Ydim/2, Xdim/2 - 10]; % first point source
source_distance_2 = [Ydim/2, Xdim/2 + 10]; % second point source

% Calculate the wave propagation to the lens for each source
Ezmx_to_lens_1 = waveamp(Ydim, Xdim, source_distance_1(1), source_distance_1(2), 1, dx, dy, lambda);
Ezmx_to_lens_2 = waveamp(Ydim, Xdim, source_distance_2(1), source_distance_2(2), 1, dx, dy, lambda);

% Sum the fields
Ezmx_to_lens_total = Ezmx_to_lens_1 + Ezmx_to_lens_2;

% Apply the lens phase shift
Ezmx_lens = Ezmx_to_lens_total .* phase_shift;

% Propagate the wave beyond the lens
Ezmx_image = waveamp(Ydim, Xdim, Ydim/2, Xdim/2, 1, dx, dy, lambda);

% Visualization
figure;
imagesc(abs(Ezmx_image));
colorbar;
title('Image of Two Point Sources');
xlabel('x (microns)');
ylabel('y (microns)');


%%%%%%%%%%%%%%%%%%%
%D)
% Parameters for incoherent light
wavelengths = linspace(450e-9, 650e-9, 10); % range of wavelengths
num_wavelengths = length(wavelengths);
phases = 2 * pi * rand(1, num_wavelengths);

% Initialize the total field
Ezmx_total = zeros(Ydim, Xdim);

for i = 1:num_wavelengths
    lambda = wavelengths(i);
    k = 2 * pi / lambda;
    
    % Generate the phase shift for the current wavelength
    phase_shift = exp(1i * k * (x_shift.^2 + y_shift.^2) / (2 * f));
    
    % Calculate the wave propagation to the lens
    Ezmx_to_lens = waveamp(Ydim, Xdim, source_position(1), source_position(2), 1, dx, dy, lambda);
    
    % Apply the lens phase shift
    Ezmx_lens = Ezmx_to_lens .* phase_shift * exp(1i * phases(i));
    
    % Propagate the wave beyond the lens
    Ezmx_image = waveamp(Ydim, Xdim, Ydim/2, Xdim/2, 1, dx, dy, lambda);
    
    % Sum the fields
    Ezmx_total = Ezmx_total + abs(Ezmx_image).^2;
end

% Visualization
figure;
imagesc(Ezmx_total);
colorbar;
title('Image with Incoherent Light');
xlabel('x (microns)');
ylabel('y (microns)');





















