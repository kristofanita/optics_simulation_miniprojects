% Theoretical Fourier Transform of Slit Apertures

% Single Slit
slit = zeros(1, Xdim);
slit(y_slit_start:y_slit_end) = 1;
F_slit = fftshift(fft(slit));

% Double Slit
slit = zeros(1, Xdim);
slit(y_slit_start:y_slit_end) = 1;
slit(y_slit_start2:y_slit_end2) = 1;
F_double_slit = fftshift(fft(slit));

% Gaussian Slit
slit = exp(-(((1:Xdim) - Xdim/2) * dy).^2 / (2*w^2));
F_gaussian_slit = fftshift(fft(slit));

% Comparison of Theoretical and Simulated Patterns

% Single Slit
figure;
subplot(2,1,1);
plot(abs(F_slit));
title('Theoretical Fourier Transform of Single Slit');
subplot(2,1,2);
plot(abs(Ezmx(round(Ydim/2), :)));
title('Simulated Diffraction Pattern of Single Slit');

% Double Slit
figure;
subplot(2,1,1);
plot(abs(F_double_slit));
title('Theoretical Fourier Transform of Double Slit');
subplot(2,1,2);
plot(abs(Ezmx(round(Ydim/2), :)));
title('Simulated Diffraction Pattern of Double Slit');

% Gaussian Slit
figure;
subplot(2,1,1);
plot(abs(F_gaussian_slit));
title('Theoretical Fourier Transform of Gaussian Slit');
subplot(2,1,2);
plot(abs(Ezmx(round(Ydim/2), :)));
title('Simulated Diffraction Pattern of Gaussian Slit');