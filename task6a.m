% Constants
lambda = 630e-9; % Wavelength in meters
D = 0.01; % Size of the simulation grid in meters
source_separation = 0.002; % Separation between the two point sources (in meters)

% Geometry for light propagation
grid_size = 500; % Size of the simulation grid (number of points)
x = linspace(-D/2, D/2, grid_size); % x-axis for the grid
y = linspace(-D/2, D/2, grid_size); % y-axis for the grid
[X, Y] = meshgrid(x, y);

% Generate spherical waves from two point sources
source1_x = grid_size / 2 + round(source_separation / (2 * (x(2) - x(1))));
source2_x = grid_size / 2 - round(source_separation / (2 * (x(2) - x(1))));
source_y = grid_size / 4;

Ez1 = waveamp(grid_size, grid_size, source_y, source1_x, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez2 = waveamp(grid_size, grid_size, source_y, source2_x, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez_total = Ez1 + Ez2;

% Calculate the intensity pattern (interference pattern)
intensity = abs(Ez_total).^2;

% Plot the intensity pattern (interference pattern)
figure;
imagesc(x, y, intensity);
title('Interference pattern from two point sources');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

% Record the interference pattern along a line (this will be the hologram)
hologram_line = intensity(grid_size / 2, :);

% Plot the hologram line on the interference pattern
figure;
plot(x, hologram_line, 'r-', 'LineWidth', 2);
title('Hologram (interference pattern along a line)');
xlabel('x (m)');
ylabel('Intensity');
grid on;
