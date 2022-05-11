function [F, u] = propFFT(F, x, lambda, f)
    % PROPFFT FFT propagation
    % An FFT propagation step, corresponding to projection from one focal
    % point of a lens to the focal point on the other side of the lens.
    
    Nx = length(x);
    dx = x(2) - x(1);
    du = lambda*f/(Nx*dx);
    u = (-Nx/2:Nx/2-1)*du;
    
    F = ifftshift(fft2(fftshift(F)));
    F = F*exp(1i*2*pi*f/lambda)/1i/lambda/f*dx*dx;
   
end