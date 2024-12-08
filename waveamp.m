%This returns a matrix with the wave amplitudes
function [Ezmx] = waveamp(Ydim,Xdim,isrc,jsrc,ampsrc,dx,dy,lambda)

    [X,Y] = meshgrid((1:Xdim)*dx, (1:Ydim)*dy);
    
    r = sqrt( (isrc*dy-Y).^2 + (jsrc*dx-X).^2 );
    
        k = 2*pi/lambda;
        
        ampsrc; 
        Ezmx = ampsrc*1./sqrt(r).*exp(1i*k*r);
        
       % Ezmx(isrc,jsrc) = 1.0/dx;
    
end
