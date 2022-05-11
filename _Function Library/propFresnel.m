function [F, u] = propFresnel(F, x, lambda, z)
    % PROPFRESNEL Fresnel propagation
    % Author : RPM
    % Implements Eq. 4-17 of Goodman
    % See also Tan lecture notes
    % See also Voelz book
    % Inputs:
    %   - F : fields in source plane
    %   - x : coordinates in source plane
    %   - lambda : wavelength
    %   - z : propagation distance
    % Outputs:
    %   - F: fields in observation plane
    %   - u: coordinates in observation plane
    
    % Calculate Fresnel number
    Lx = max(x) - min(x);
    Fnumber = Lx^2/lambda/z;
    if Fnumber > 1
        warning(['Not operating in Fresnel region, F = ' num2str(Fnumber)])
    end
    
    % Wavevector
    k = 2*pi/lambda;

    % Coord calculations
    Nx = length(x);
    dx = (max(x) - min(x))/(Nx-1);
    [x_mesh, y_mesh] = meshgrid(x, x);
    du = lambda*z/(Nx*dx);
    u = (-Nx/2:Nx/2-1)*du;
    [u_mesh, v_mesh] = meshgrid(u, u);

    % Fresnel calculations
    F = F.*exp(1i*k/(2*z)*(x_mesh.^2 + y_mesh.^2)); 
    F = ifftshift(fft2(fftshift(F)));
    F = F.*exp(1i*k/(2*z)*(u_mesh.^2 + v_mesh.^2));
    F = F * exp(1i*k*z)/(1i*lambda*z)*dx*dx;
    
end