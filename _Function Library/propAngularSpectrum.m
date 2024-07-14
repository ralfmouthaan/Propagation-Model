% Think this works, haven't properly checked it.

function F = propAngularSpectrum(F, x, lambda, z)

    Nx = length(x);
    dx = x(2) - x(1);
    du = 1/dx;
    u = -Nx/2 + 1/2 : Nx/2 - 1/2;
    u = u*du/length(u);
    
    filter = exp(1i*2*pi*z.*sqrt(1/lambda^2 - u.^2 - u.'.^2));
    
    F = fftshift(fft2(fftshift(F)));
    F = F.*filter;
    F = fftshift(ifft2(fftshift(F)));    
   
end
