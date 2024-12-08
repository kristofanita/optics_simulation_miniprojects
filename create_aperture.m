function aperture = create_aperture(Xdim, aperture_type, params)
    aperture = zeros(1, Xdim);
    switch aperture_type
        case 'rectangular'
            width = params.width;
            center = Xdim / 2;
            aperture(center - width/2 : center + width/2) = 1;
        case 'double_slit'
            width = params.width;
            separation = params.separation;
            center = Xdim / 2;
            aperture(center - width/2 - separation : center + width/2 - separation) = 1;
            aperture(center - width/2 + separation : center + width/2 + separation) = 1;
        case 'gaussian'
            sigma = params.sigma;
            x = linspace(-Xdim/2, Xdim/2, Xdim);
            aperture = exp(-x.^2 / (2 * sigma^2));
        otherwise
            error('Unknown aperture type');
    end
end
