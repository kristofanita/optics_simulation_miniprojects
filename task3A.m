% Constants
lambda = 630e-9; % Wavelength in meters
k1 = 2 * pi / lambda; % Wavenumber in the first medium
n1 = 1; % Refractive index of the first medium
n2 = 1.5; % Refractive index of the second medium
k2 = n2 * k1; % Wavenumber in the second medium

% Geometry
h1 = 0.01; % Distance from source to medium boundary (in meters)
h2 = 0.01; % Distance from medium boundary to detector (in meters)
d = 0.02; % Total distance (in meters)
num_segments = 200; % Number of 1mm segments
x = linspace(0, d, num_segments); % x-axis divided into segments

% Initialize array for amplitudes
amplitudes = zeros(size(x));

% Calculate amplitudes for each path
for i = 1:length(x)
    % Distance in the first medium
    s1 = sqrt(x(i)^2 + h1^2);
    % Distance in the second medium
    s2 = sqrt((d - x(i))^2 + h2^2);
    % Calculate the complex amplitude
    E = exp(1i * k1 * s1) * exp(1i * k2 * s2);
    amplitudes(i) = E;
end

% Sum the amplitudes to get the resultant field
resultant_field = sum(amplitudes);

% Plot the geometry and light paths
figure;
hold on;

% Plot the source
plot(0, h1, 'ro', 'MarkerSize', 8, 'DisplayName', 'Source');
text(0, h1 + 0.001, 'Source', 'HorizontalAlignment', 'right');

% Plot the detector
plot(d, -h2, 'bo', 'MarkerSize', 8, 'DisplayName', 'Detector');
text(d, -h2 - 0.001, 'Detector', 'HorizontalAlignment', 'left');

% Plot the medium boundary
plot([0 d], [0 0], 'k--', 'DisplayName', 'Medium Boundary');
text(d / 2, 0.001, 'Medium Boundary', 'HorizontalAlignment', 'center');

% Plot the light paths for different x
for i = 1:length(x)
    % Points in first medium
    x1 = linspace(0, x(i), 100);
    y1 = h1 * (1 - x1 / x(i));
    plot(x1, y1, 'r-', 'HandleVisibility', 'off');

    % Points in second medium
    x2 = linspace(x(i), d, 100);
    y2 = -h2 * (x2 - x(i)) / (d - x(i));
    plot(x2, y2, 'b-', 'HandleVisibility', 'off');
end

% Plot a straight path as a reference
plot([0 d], [h1 -h2], 'g-', 'LineWidth', 0.5, 'DisplayName', 'Straight path of light');
theta_straight = atan((h1 + h2) / d);
eqn_straight = sprintf('y = %.2f * x + %.2f', tan(theta_straight), h1);
%text(d/2, (h1 - h2)/2, ['Straight Path: ' eqn_straight], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

title('Light paths with n1 = n2 = 1');
xlabel('x (m)');
ylabel('y (m)');
legend;
grid on;
axis equal;
hold off;

% Display the resultant field
disp('Resultant field:');
disp(resultant_field);


%%

% Define the observation points on the detector screen
obs_points = linspace(-0.01, 0.01, 500); % Observation points on the screen

% Initialize intensity array
intensity = zeros(size(obs_points));

% Calculate the Huygens-Fresnel principle
for i = 1:length(obs_points)
    % Initialize the resultant field for each observation point
    resultant_field = 0;
    for j = 1:length(x)
        % Calculate distances
        s1 = sqrt(x(j)^2 + h1^2);
        s2 = sqrt((d - x(j))^2 + (obs_points(i) + h2)^2);
        % Calculate the complex amplitude
        E = exp(1i * k1 * s1) * exp(1i * k2 * s2);
        resultant_field = resultant_field + E;
    end
    % Calculate the intensity
    intensity(i) = abs(resultant_field)^2;
end

% Plot the intensity pattern
figure;
plot(obs_points, intensity);
title('Interference Pattern on Detector Screen');
xlabel('Position on Detector Screen (m)');
ylabel('Intensity');
grid on;
