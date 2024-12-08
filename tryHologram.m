% Constants
lambda = 630e-9; % Wavelength in meters
D = 0.0025; % Size of the simulation grid in meters
source_separation = 0.0025; % Separation between the two point sources (in meters)

% Geometry for light propagation
grid_size = 500; % Size of the simulation grid (number of points)
x = linspace(-D, D, grid_size); % x-axis for the grid
y = linspace(0, 0.02, grid_size); % y-axis for the grid
[X, Y] = meshgrid(x, y);

% Positions of the point sources
source_y = grid_size - round(grid_size / 5); % y = 0.018
source1_x = grid_size / 2 + round(source_separation / (2 * (x(2) - x(1))));
source2_x = grid_size / 2 - round(source_separation / (2 * (x(2) - x(1))));

% Generate spherical waves from two point sources
Ez1 = waveamp(grid_size, grid_size, source_y, source1_x, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez2 = waveamp(grid_size, grid_size, source_y, source2_x, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez_total = Ez1 + Ez2;

% Calculate the intensity pattern (interference pattern)
intensity = abs(Ez_total).^2;

% Extract the interference pattern until y = 0
y_zero_idx = round(grid_size * 0.5); % Index corresponding to y = 0
intensity_below_zero = intensity(1:y_zero_idx, :);
y_below_zero = y(1:y_zero_idx);

% Save the interference values along y = 0 to holoGrid
holoGrid = intensity(y_zero_idx, :);

% Plot the intensity pattern (interference pattern) until y = 0
figure;
imagesc(x, y_below_zero, intensity_below_zero);
title('Interference pattern from two point sources (until y = 0)');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

% Plot the point sources on the intensity pattern
hold on;
plot(x(source1_x), 0.01, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Source 1');
plot(x(source2_x), 0.01, 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Source 2');
plot(x, zeros(1, grid_size), 'k-', 'LineWidth', 2, 'DisplayName', 'y = 0 Line');
hold off;
legend;
grid on;

% Plot the saved interference values along y = 0 (holoGrid)
figure;
plot(x, holoGrid, 'k-', 'LineWidth', 2);
title('Interference values along y = 0 (holoGrid)');
xlabel('x (m)');
ylabel('Intensity');
grid on;
