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

% Calculate the angles for the light paths using Snell's law
%incident_angle = pi / 4; % Example incident angle of 45 degrees
incident_angle = pi / 3;
theta1 = incident_angle;
theta2 = asin(n1 / n2 * sin(theta1));

% All paths
theta1_all = atan(h1 ./ x);
theta2_all = asin(n1 / n2 * sin(theta1_all));

% Calculate points for the incident and refracted rays
boundary_x = h1 / tan(theta1); % x-coordinate at the medium boundary
incident_x = linspace(0, boundary_x, 100);
incident_y = h1 - tan(theta1) * incident_x;

refracted_x = linspace(boundary_x, d, 100);
refracted_y = tan(theta2) * (refracted_x - boundary_x);

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
    y2 = tan(theta2_all(i)) * (x2 - x(i));
    plot(x2, -y2, 'b-', 'HandleVisibility', 'off'); % Note the change in y direction
end

% Plot the incident and refracted rays
plot(incident_x, incident_y, 'g-', 'LineWidth', 2, 'DisplayName', 'Incident Ray');
plot(refracted_x, -refracted_y, 'g-', 'LineWidth', 2, 'DisplayName', 'Refracted Ray');

title("Light paths with refraction according to Snell's law (n1 < n2)");
xlabel('x (m)');
ylabel('y (m)');
legend;
grid on;
axis equal;
hold off;

% Display the resultant field
disp('Resultant field:');
disp(resultant_field);
