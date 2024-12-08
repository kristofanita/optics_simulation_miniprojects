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

% Generate spherical waves from two close point sources
Ez_before_lens1 = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2 + 10, size(X, 2) / 2, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez_before_lens2 = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2 - 10, size(X, 2) / 2, 1, x(2) - x(1), y(2) - y(1), lambda);
Ez_before_lens = Ez_before_lens1 + Ez_before_lens2;

% Apply lens effect (phase shift)
Ez_after_lens = Ez_before_lens .* tau;

% Propagation from the lens to the focal point
Ez_focused = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2, size(X, 2) / 2, Ez_after_lens, x(2) - x(1), y(2) - y(1), lambda);

% Plot the light paths
figure;
hold on;

% Plot the lens as a line
plot([0, 0], [-0.01, 0.01], 'k-', 'LineWidth', 2, 'DisplayName', 'Lens');

% Plot the spherical waves before the lens (from two point sources)
source_x = -source_distance;
source_y1 = 0.0005;
source_y2 = -0.0005;
for i = 1:20:length(x)
    plot([source_x, 0], [source_y1, y(i)], 'r-', 'HandleVisibility', 'off');
    plot([source_x, 0], [source_y2, y(i)], 'm-', 'HandleVisibility', 'off');
end

% Plot the focused rays after the lens
for i = 1:20:length(x)
    plot([0, focal_distance], [y(i), 0], 'b-', 'HandleVisibility', 'off');
end

% Add labels and title
title('Light paths through a lens from two point sources');
xlabel('x (m)');
ylabel('y (m)');
legend;
grid on;
axis equal;
hold off;

% Plot the intensity pattern at the focal point
intensity = abs(Ez_focused).^2;
figure;
imagesc(x, y, intensity);
title('Intensity pattern at the focal point');
xlabel('x (m)');
ylabel('y (m)');
colorbar;
display("A lencse véges felbontásának megértése érdekében a diffrakció hatását kell megvizsgálnunk. A diffrakció hatására a lencse által fókuszált kép nem tökéletes pont lesz, hanem egy Airy-korong, ami a közeli pontforrások megkülönböztetésének korlátját adja.");


%%

% Constants
lambda = 630e-9; % Wavelength in meters
focal_length = 0.02; % Focal length of the lens in meters

% Geometry for light propagation
source_distance = 0.01; % Distance from source to lens (in meters)
focal_distance = 0.02; % Distance from lens to focal point (in meters)

% Lens diameters for experiments
lens_diameters = [10e-6, 100e-6, 500e-6]; % Different diameters in meters
num_lenses = length(lens_diameters);

% Create a figure for intensity patterns
figure;
tiledlayout(1, num_lenses);

for k = 1:num_lenses
    D = lens_diameters(k);
    theta_max = atan(D / (2 * focal_length));
    NA = sin(theta_max); % Numerical aperture

    % Phase shift introduced by the lens
    x = linspace(-D/2, D/2, 1000); % x-axis for the lens
    y = linspace(-D/2, D/2, 1000); % y-axis for the lens
    [X, Y] = meshgrid(x, y);
    delta = (2 * pi / lambda) * (1 - 1) * sqrt(X.^2 + Y.^2 + focal_length^2);
    tau = exp(-1i * delta); % Transmittance of the lens

    % Generate spherical waves from two close point sources
    Ez_before_lens1 = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2 + 5, size(X, 2) / 2, 1, x(2) - x(1), y(2) - y(1), lambda);
    Ez_before_lens2 = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2 - 5, size(X, 2) / 2, 1, x(2) - x(1), y(2) - y(1), lambda);
    Ez_before_lens = Ez_before_lens1 + Ez_before_lens2;

    % Apply lens effect (phase shift)
    Ez_after_lens = Ez_before_lens .* tau;

    % Propagation from the lens to the focal point
    Ez_focused = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2, size(X, 2) / 2, Ez_after_lens, x(2) - x(1), y(2) - y(1), lambda);

    % Calculate the intensity pattern at the focal point
    intensity = abs(Ez_focused).^2;

    % Plot the intensity pattern
    nexttile;
    imagesc(x, y, intensity);
    title(['Lens Diameter: ' num2str(D*1e6) ' \mum']);
    xlabel('x (m)');
    ylabel('y (m)');
    colorbar;
end

display("két pontszerű forrás akkor különböztethető meg, ha a köztük lévő távolság legalább:")
display("d_min=lambda/(2*NA)")
display("ahol lambda a fény hullámhossza, és 𝑁𝐴 a numerikus apertúra, amely a következőképpen definiálható:")
display("NA=n*sin(theta)")

%%

% Constants
lambda = 630e-9; % Wavelength in meters
focal_length = 0.02; % Focal length of the lens in meters

% Geometry for light propagation
source_distance = 0.01; % Distance from source to lens (in meters)
focal_distance = 0.02; % Distance from lens to focal point (in meters)

% Lens diameters for experiments
lens_diameters = [10e-6, 100e-6, 500e-6]; % Different diameters in meters
num_lenses = length(lens_diameters);

% Create a figure for intensity patterns
figure;
tiledlayout(1, num_lenses);

for k = 1:num_lenses
    D = lens_diameters(k);
    theta_max = atan(D / (2 * focal_length));
    NA = sin(theta_max); % Numerical aperture

    % Phase shift introduced by the lens
    x = linspace(-D/2, D/2, 1000); % x-axis for the lens
    y = linspace(-D/2, D/2, 1000); % y-axis for the lens
    [X, Y] = meshgrid(x, y);
    delta = (2 * pi / lambda) * (1 - 1) * sqrt(X.^2 + Y.^2 + focal_length^2);
    tau = exp(-1i * delta); % Transmittance of the lens

    % Generate spherical waves from two close point sources
    Ez_before_lens1 = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2 + 5, size(X, 2) / 2, 1, x(2) - x(1), y(2) - y(1), lambda);
    Ez_before_lens2 = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2 - 5, size(X, 2) / 2, 1, x(2) - x(1), y(2) - y(1), lambda);
    Ez_before_lens = Ez_before_lens1 + Ez_before_lens2;

    % Apply lens effect (phase shift)
    Ez_after_lens = Ez_before_lens .* tau;

    % Propagation from the lens to the focal point
    Ez_focused = waveamp(size(Y, 1), size(X, 2), size(Y, 1) / 2, size(X, 2) / 2, Ez_after_lens, x(2) - x(1), y(2) - y(1), lambda);

    % Calculate the intensity pattern at the focal point
    intensity = abs(Ez_focused).^2;

    % Plot the intensity pattern
    nexttile;
    imagesc(x, y, intensity);
    title(['Lens Diameter: ' num2str(D*1e6) ' \mum']);
    xlabel('x (m)');
    ylabel('y (m)');
    colorbar;
    
    % Calculate Abbe resolution limit
    abbe_limit = lambda / (2 * NA);
    disp(['Abbe Limit for lens diameter ' num2str(D*1e6) ' µm: ', num2str(abbe_limit), ' meters']);
end

% Add an overall title
sgtitle('Intensity patterns at the focal point for different lens diameters');



