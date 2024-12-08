% RECORD HOLOGRAM
% Constants
lambda = 630e-9; % Wavelength in meters
D = 0.0025; % Size of the simulation grid in meters
source_separation = 0.0025; % Separation between the two point sources (in meters)
threshold = 50; % Threshold for determining low intensity 

% Geometry for light propagation
num_points = 500; % Size of the simulation grid (number of points)
x = linspace(-D, D, num_points); % x-axis for the grid
y = linspace(0, 0.02, num_points); % y-axis for the grid
[X, Y] = meshgrid(x, y);

% Positions of the point sources
source_y = num_points - round(num_points / 5); % y = 0.018
source1_x = num_points / 2 + round(source_separation / (2 * (x(2) - x(1))));
source2_x = num_points / 2 - round(source_separation / (2 * (x(2) - x(1))));

% Generate spherical waves from two point sources
Ez1 = waveamp(num_points, num_points, source_y, source1_x, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez2 = waveamp(num_points, num_points, source_y, source2_x, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez_total = Ez1 + Ez2;

% Calculate the intensity pattern (interference pattern)
intensity = abs(Ez_total).^2;

% Extract the interference pattern at y = 0
y_zero_idx = round(num_points * 0.5); % Index corresponding to y = 0
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
%plot(x, zeros(1, grid_size), 'k-', 'LineWidth', 2, 'DisplayName', 'y = 0 Line');

high_intensity = holoGrid >= threshold;
low_intensity = holoGrid < threshold;
edges = [];
gaps = [];
% Plot the high intensity sections as continuous lines and calculate gaps
for i = 1:length(x) - 1
    if high_intensity(i) && high_intensity(i + 1)
        plot(x(i:i+1), [0 0], 'k.', 'LineWidth', 2);
    % else
    %    gaps = [gaps, x(i + 1) - x(i)];
    %    plot(x(i:i+1), [save_y save_y], 'r--', 'LineWidth', 2);
    end
    if high_intensity(i) && low_intensity(i-1) && high_intensity(i+1)
        plot(x(i), 0, "k|", 'LineWidth', 2)
        edges = [edges, x(i)];
    elseif high_intensity(i) && low_intensity(i+1) && high_intensity(i-1)
        plot(x(i), 0, "k|", 'LineWidth', 2)
        edges = [edges, x(i)];
    end
end
gaps = diff(edges);
hold off;
legend;
grid on;

% Plot the grating positions as a dashed line with gaps
figure;
hold on;
high_intensity = holoGrid >= threshold;
low_intensity = holoGrid < threshold;

% Plot the high intensity sections as continuous lines
for i = 1:length(x) - 1
    if high_intensity(i) && high_intensity(i + 1)
        plot(x(i:i+1), [1 1], 'k-', 'LineWidth', 2);
    end
end

title('Grating based on low intensity values');
xlabel('x (m)');
ylabel('y (arbitrary units)');
ylim([0.95 1.05]);
grid on;
hold off;


%%
%%%%%%%%%%%%%%% KEPALKOTAS %%%%%%%%%%%%%%%%
% Constants
lambda = 630e-9; % Wavelength in meters
d = 0.02; % Distance to detector in meters
h1 = 0.01; % Distance from source to medium boundary in meters
h2 = 0.01; % Distance from medium boundary to detector in meters

% Grid definition
num_points = 500; % Size of the simulation grid (number of points)
x = linspace(0, d, num_points); % x-axis for the grid
y = linspace(-h1, h2, num_points); % y-axis for the grid
[X, Y] = meshgrid(x, y);
distance_per_index = d / num_points;

% Position of the point source
source_x = 0;
source_y = h1;

low_intensity_indices = find(holoGrid < threshold);
gap_sizes = [];
gap_positions = [];

i = 1;
while i <= length(low_intensity_indices)
    start_index = low_intensity_indices(i);
    
    j = i;
    while j < length(low_intensity_indices) && low_intensity_indices(j + 1) == low_intensity_indices(j) + 1
        j = j + 1;
    end
    
    end_index = low_intensity_indices(j);
    
    gap_size = (end_index - start_index + 1) * distance_per_index;
    gap_sizes = [gap_sizes, gap_size];
    
    gap_center = (start_index + end_index) / 2 * distance_per_index;
    gap_positions = [gap_positions, gap_center];

    i = j + 1;
end

% Calculate the path lengths for light propagating from source to each point on the grid
s1 = sqrt((X - source_x).^2 + (Y - source_y).^2);

% Generate the wavefronts emanating from the point source
Ez = exp(1i * 2 * pi * s1 / lambda);

% Initialize the total field
Ez_total = zeros(size(Ez));

% Plot the light paths from the source to the grating
figure;
hold on;

% Plot light paths from source to grating
for i = 1:length(gap_positions)
    plot([source_x gap_positions(i)], [source_y 0], 'r-'); % Red lines for paths to grating
end

% Plot the grating positions
for i = 1:length(gap_positions)
    plot(gap_positions(i), 0, 'k*', 'LineWidth', 2); % Grating slits
end

% Using Huygens principle to calculate wavefronts from each gap
for i = 1:length(gap_positions)-1
    gap_center = (gap_positions(i) + gap_positions(i + 1)) / 2;
    gap_idx = find(abs(x - gap_center) == min(abs(x - gap_center)));
    Ez_gap = waveamp(num_points, num_points, y_zero_idx, gap_idx(1), 1, distance_per_index, y(2) - y(1), lambda);
                    %(Ydim,       Xdim,       isrc,       jsrc,  ampsrc,dx,               dy,           lambda)
    Ez_total = Ez_total + Ez_gap;
end

% Calculate the intensity pattern from the secondary wavelets
intensity_total = abs(Ez_total).^2;

% Extract the interference pattern at y = 0
y_zero_idx = round(num_points * 0.5); % Index corresponding to y = 0
intensity_above_zero = intensity(1:y_zero_idx, :);
y_above_zero = y(1:y_zero_idx);

% Plot the wavefronts from the secondary wavelets using contour plot
imagesc(x, y_above_zero, intensity_above_zero);
title('Light paths with grating');
xlabel('x (m)');
ylabel('y (m)');
colorbar;
plot(source_x, source_y, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Source');

legend;
grid on;
hold off;

figure;
hold on;
imagesc(x, y, intensity_total);
title('Hologram total EZ intenzity')
xlabel('x (m)');
ylabel('y (m)');
colorbar;

figure;
plot(x, holoGrid, 'k-', 'LineWidth', 2);
title('Interference values along y = 0 (holoGrid)');
xlabel('x (m)');
ylabel('Intensity');
grid on;