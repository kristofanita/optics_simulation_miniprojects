% Constants
lambda = 630e-9; % Wavelength in meters
D = 0.01; % Size of the simulation grid in meters
source_separation = 0.002; % Separation between the two point sources (in meters)

% Geometry for light propagation
grid_size = 500; % Size of the simulation grid (number of points)
x = linspace(-D/2, D/2, grid_size); % x-axis for the grid
y = linspace(-D/2, D/2, grid_size); % y-axis for the grid
[X, Y] = meshgrid(x, y);

% Generate spherical waves from two point sources (recording the hologram)
source1_x = grid_size / 2 + round(source_separation / (2 * (x(2) - x(1))));
source2_x = grid_size / 2 - round(source_separation / (2 * (x(2) - x(1))));
source_y = grid_size / 4;

Ez1 = waveamp(grid_size, grid_size, source_y, source1_x, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez2 = waveamp(grid_size, grid_size, source_y, source2_x, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez_total = Ez1 + Ez2;

% Calculate the intensity pattern (interference pattern)
intensity = abs(Ez_total).^2;

% Record the interference pattern along a line (this will be the hologram)
hologram_line = intensity(grid_size / 2, :);

% Plot the hologram line
figure;
plot(x, hologram_line, 'r-', 'LineWidth', 2);
title('Hologram (interference pattern along a line)');
xlabel('x (m)');
ylabel('Intensity');
grid on;

% Use the hologram as a grating to illuminate and generate the image

% Create a new wavefront by using the hologram as a source
Ez_hologram = zeros(grid_size, grid_size);
Ez_hologram(grid_size / 2, :) = hologram_line .* exp(1i * 2 * pi * X(grid_size / 2, :) / lambda);

% Generate the image using the hologram
% Use one of the point sources to illuminate the hologram
source3_y = grid_size / 4;
Ez3 = waveamp(grid_size, grid_size, source3_y, grid_size / 2, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez_reconstructed = Ez_hologram .* exp(1i * 2 * pi * sqrt((X - X(grid_size / 2, :)).^2 + (Y - Y(source3_y, :)).^2) / lambda);

% Calculate the intensity pattern of the reconstructed image
intensity_reconstructed = abs(Ez_reconstructed).^2;

% Plot the intensity pattern of the reconstructed image
figure;
imagesc(x, y, intensity_reconstructed);
title('Reconstructed image from the hologram');
xlabel('x (m)');
ylabel('y (m)');
colorbar;
hold on;

% Plot the point source
plot(0, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Point Source');
legend;
hold off;
