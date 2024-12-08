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

figure;
hold on;
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

% Plot the wavefronts from the secondary wavelets using imagesc
pcolor(x, y, intensity_total); % Imagesc plot to show intensity pattern
title('Light paths with grating');
xlabel('x (m)');
ylabel('y (m)');
colorbar;
axis xy; % Ensure the y-axis is oriented correctly

hold on;
plot(source_x, source_y, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Source');
plot(d, -h2, 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Detector');

legend;
grid on;
hold off;


%%

