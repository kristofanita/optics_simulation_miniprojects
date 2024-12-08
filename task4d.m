% Constants
lambda_center = 550e-9; % Central wavelength in meters (green light)
delta_lambda = 50e-9; % Range of wavelengths (500 nm to 600 nm)
num_wavelengths = 10; % Number of different wavelengths
focal_length = 0.02; % Focal length of the lens in meters

% Geometry for light propagation
source_distance = 0.01; % Distance from source to lens (in meters)
focal_distance = 0.02; % Distance from lens to focal point (in meters)

% Lens diameter
D = 100e-6; % Diameter of the lens in meters

% Generate different wavelengths
lambdas = linspace(lambda_center - delta_lambda / 2, lambda_center + delta_lambda / 2, num_wavelengths);

% Initialize the intensity pattern
total_intensity = zeros(1000, 1000);

% Create a figure for light paths
figure;
hold on;

% Plot the lens as a line
plot([0, 0], [-0.01, 0.01], 'k-', 'LineWidth', 2, 'DisplayName', 'Lens');

% Plot the spherical waves before the lens (from two point sources)
source_x = -source_distance;
source_y1 = 0.0005;
source_y2 = -0.0005;

% Initialize figure for intensity patterns
figure;
tiledlayout(1, num_wavelengths);

for lambda = lambdas
    % Calculate numerical aperture
    theta_max = atan(D / (2 * focal_length));
    NA = sin(theta_max); % Numerical aperture

    % Phase shift introduced by the lens
    x = linspace(-D/2, D/2, 1000); % x-axis for the lens
    y = linspace(-D/2, D/2, 1000); % y-axis for the lens
    [X, Y] = meshgrid(x, y);
    delta = (2 * pi / lambda) * (1 - 1) * sqrt(X.^2 + Y.^2 + focal_length^2);
    tau = exp(-1i * delta); % Transmittance of the lens

    % Generate spherical waves from two close point sources with random phase
    phase1 = rand() * 2 * pi;
    phase2 = rand() * 2 * pi;
    Ez_before_lens1 = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2 + 5, size(X, 2) / 2, exp(1i * phase1), x(2) - x(1), y(2) - y(1), lambda);
    Ez_before_lens2 = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2 - 5, size(X, 2) / 2, exp(1i * phase2), x(2) - x(1), y(2) - y(1), lambda);
    Ez_before_lens = Ez_before_lens1 + Ez_before_lens2;

    % Apply lens effect (phase shift)
    Ez_after_lens = Ez_before_lens .* tau;

    % Propagation from the lens to the focal point
    Ez_focused = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2, size(X, 2) / 2, Ez_after_lens, x(2) - x(1), y(2) - y(1), lambda);

    % Calculate the intensity pattern at the focal point
    intensity = abs(Ez_focused).^2;
    
    % Accumulate the intensity pattern
    total_intensity = total_intensity + intensity;

    % Plot the light paths
    figure(1);
    for i = 1:100:length(x)
        plot([source_x, 0], [source_y1, y(i)], 'r-', 'HandleVisibility', 'off');
        plot([source_x, 0], [source_y2, y(i)], 'r-', 'HandleVisibility', 'off');
        plot([0, focal_distance], [y(i), 0], 'b-', 'HandleVisibility', 'off');
    end

    % Plot the intensity pattern for each wavelength
    figure(2);
    nexttile;
    imagesc(x, y, intensity);
    wavelength = num2str(round(lambda*1e9));
    title([wavelength ' nm']);
    xlabel('x (m)');
    ylabel('y (m)');
    
end
colorbar;
% Normalize the intensity pattern
total_intensity = total_intensity / num_wavelengths;

% Plot the total intensity pattern
figure;
imagesc(x, y, total_intensity);
title('Total intensity pattern at the focal point for white light');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

% Add labels and title to the light paths figure
figure(1);
title('Light paths through a lens from two point sources');
xlabel('x (m)');
ylabel('y (m)');
legend;
grid on;
axis equal;
hold off;
