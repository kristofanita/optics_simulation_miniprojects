% Parameters
lambda = 630e-9;
d = 200e-6;
Xdim = 1000;
Ydim = 1000;
dx = 1e-6;
dy = 1e-6;

% Define aperture functions
rect_params.width = 10; % 10 microns
double_slit_params.width = 2; % 2 microns
double_slit_params.separation = 20; % 20 microns
gaussian_params.sigma = 5; % 5 microns

aperture_rect = create_aperture(Xdim, 'rectangular', rect_params);
aperture_double_slit = create_aperture(Xdim, 'double_slit', double_slit_params);
aperture_gaussian = create_aperture(Xdim, 'gaussian', gaussian_params);

% Compute diffraction patterns using waveamp.m
Ez_rect = zeros(Ydim, Xdim);
Ez_double_slit = zeros(Ydim, Xdim);
Ez_gaussian = zeros(Ydim, Xdim);

for x = 1:Xdim
    if aperture_rect(x) > 0
        Ez_rect = Ez_rect + waveamp(Ydim, Xdim, Ydim/2, x, aperture_rect(x), dx, dy, lambda);
    end
    if aperture_double_slit(x) > 0
        Ez_double_slit = Ez_double_slit + waveamp(Ydim, Xdim, Ydim/2, x, aperture_double_slit(x), dx, dy, lambda);
    end
    if aperture_gaussian(x) > 0
        Ez_gaussian = Ez_gaussian + waveamp(Ydim, Xdim, Ydim/2, x, aperture_gaussian(x), dx, dy, lambda);
    end
end

% Compute theoretical Fourier Transforms
f_x = linspace(-1/(2*dx), 1/(2*dx), Xdim);
rect_FT = abs(fftshift(fft(aperture_rect)));
double_slit_FT = abs(fftshift(fft(aperture_double_slit)));
gaussian_FT = abs(fftshift(fft(aperture_gaussian)));

% Visualization
figure;
subplot(3,2,1); plot(aperture_rect); title('Rectangular Aperture');
subplot(3,2,2); imagesc(abs(Ez_rect)); title('Rectangular Aperture Diffraction Pattern');

subplot(3,2,3); plot(aperture_double_slit); title('Double Slit Aperture');
subplot(3,2,4); imagesc(abs(Ez_double_slit)); title('Double Slit Diffraction Pattern');

subplot(3,2,5); plot(aperture_gaussian); title('Gaussian Aperture');
subplot(3,2,6); imagesc(abs(Ez_gaussian)); title('Gaussian Diffraction Pattern');

figure;
subplot(3,1,1); plot(f_x, rect_FT); title('Rectangular Aperture Fourier Transform');
subplot(3,1,2); plot(f_x, double_slit_FT); title('Double Slit Fourier Transform');
subplot(3,1,3); plot(f_x, gaussian_FT); title('Gaussian Fourier Transform');
