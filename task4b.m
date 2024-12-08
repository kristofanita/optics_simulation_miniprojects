% Constants
lambda = 630e-9; % Wavelength in meters
n_lens = 1.5; % Refractive index of the lens
focal_length = 0.02; % Focal length of the lens in meters

% Geometry
x = linspace(-0.01, 0.01, 1000); % x-axis for the lens
y = linspace(-0.01, 0.01, 1000); % y-axis for the lens
[X, Y] = meshgrid(x, y);

% Calculate the phase shift across the lens surface
delta = (2 * pi / lambda) * (n_lens - 1) * sqrt(X.^2 + Y.^2 + focal_length^2);
tau = exp(-1i * delta); % Transmittance of the lens

% Geometry for light propagation
source_distance = 0.01; % Distance from source to lens (in meters)
focal_distance = 0.02; % Distance from lens to focal point (in meters)

% Generate spherical wave from point source using Huygens-Fresnel principle
Ez_before_lens = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2, size(X, 2) / 2, 1, x(2) - x(1), y(2) - y(1), lambda);

% Apply lens effect (phase shift)
Ez_after_lens = Ez_before_lens .* tau;

% Propagation from the lens to the focal point
Ez_focused = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2, size(X, 2) / 2, Ez_after_lens, x(2) - x(1), y(2) - y(1), lambda);

% Plot the light paths
figure;
hold on;

% Plot the lens as a line
plot([0, 0], [-0.01, 0.01], 'k-', 'LineWidth', 2, 'DisplayName', 'Lens');

% Plot the spherical waves before the lens (from a single point source)
source_x = -source_distance;
source_y = 0;
for i = 1:20:length(x)
    plot([source_x, 0], [source_y, y(i)], 'r-', 'HandleVisibility', 'off');
end

% Plot the focused rays after the lens
for i = 1:20:length(x)
    plot([0, focal_distance], [y(i), 0], 'b-', 'HandleVisibility', 'off');
end

% Add labels and title
title('Light paths through a lens from a point source');
xlabel('x (m)');
ylabel('y (m)');
legend;
grid on;
axis equal;
hold off;

% Display the resultant field
disp('Resultant field:');
disp(Ez_focused);
