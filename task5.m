% Parameters
c = 3e8;                % Speed of light in vacuum (m/s)
dx = 1e-8;              % Spatial step (10 nm)
dt = dx / (2 * c);      % Time step (Courant condition)
Nx = 2000;              % Number of spatial points
Nt = 3000;              % Number of time steps
lambda = 500e-9;        % Wavelength (500 nm)
f = c / lambda;         % Frequency
omega = 2 * pi * f;     % Angular frequency
t0 = 30;                % Center of the Gaussian pulse
spread = 10;            % Width of the Gaussian pulse

% Material properties
n_air = 1.0;
n_glass = 1.5;
d = 20e-6;              % Thickness of the glass plate (20 microns)

% Refractive index profile
n = ones(1, Nx);
glass_start = round(Nx/3);
glass_end = glass_start + round(d/dx);
n(glass_start:glass_end) = n_glass;

% Field arrays
Ez = zeros(1, Nx);
Hy = zeros(1, Nx);

% Simulation loop
for t = 1:Nt
    % Update magnetic field
    for i = 1:Nx-1
        Hy(i) = Hy(i) + (dt / (mu0 * dx)) * (Ez(i+1) - Ez(i));
    end
    
    % Update electric field
    for i = 2:Nx
        Ez(i) = Ez(i) + (dt / (epsilon0 * n(i)^2 * dx)) * (Hy(i) - Hy(i-1));
    end
    
    % Source (Gaussian pulse)
    Ez(50) = Ez(50) + exp(-0.5 * ((t - t0) / spread)^2) * sin(omega * t * dt);
    
    % Absorbing Boundary Conditions (PML)
    Ez(1) = 0;
    Ez(Nx) = 0;
    
    % Plot fields
    if mod(t, 50) == 0
        plot(Ez);
        ylim([-1 1]);
        title(['Time step: ', num2str(t)]);
        drawnow;
    end
end



%%%%%%%%%%%%%%%%
%B)
% Measure fields before and after the glass plate
incident_field = Ez(1:glass_start-1);
transmitted_field = Ez(glass_end+1:end);

% Calculate power
incident_power = sum(incident_field.^2);
transmitted_power = sum(transmitted_field.^2);
reflection_power = sum((Ez(1:glass_start-1) - incident_field).^2);

% Transmission and reflection coefficients
T = transmitted_power / incident_power;
R = reflection_power / incident_power;

disp(['Transmission Coefficient: ', num2str(T)]);
disp(['Reflection Coefficient: ', num2str(R)]);


%%%%%%%%%%%%%%%%%
%C)
% Adjust thickness for maximum transparency
optimal_thickness = lambda / (2 * n_glass);
disp(['Optimal thickness for maximum transparency: ', num2str(optimal_thickness), ' meters']);