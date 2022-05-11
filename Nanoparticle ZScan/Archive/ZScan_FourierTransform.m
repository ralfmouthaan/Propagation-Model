% Ralf Mouthaan
% University of Cambridge
% June 2020
% 
% Modelling of z-scan system.
% Uses Fourier transforms for propagation.
% This approach tries to use Fourier transforms for all propagations, but I
% think I'm getting into trouble with the fact that I'm not actually in the
% Fourier domain.

clc; clear variables; close all;

%% User-defined parameters

n0 = 1.45; % n = n0 + n2*|E|^2
n2 = 5e-6;
lambda = 532e-9; % Green DPSS
w0 = 5e-3; % width of beam before lens
f = 75e-3; % focal length of lens
d = 10; % Distance from sample to diode. Far enough to be far-field?
L = 300e-6;

%% Derived parameters

k0 = 2*pi/lambda;
wR = lambda*f/(w0*pi);
zR= k0*wR^2/2; % diffraction length of beam
x = linspace(-10*w0, 10*w0, 2000);
[x_mesh, y_mesh] = meshgrid(x, x.');
r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
arrz = linspace(-4*zR,4*zR,100);

%% Calculations

figure;
DiodePower = zeros(size(arrz));
for ii = 1:length(arrz)
    
    z = arrz(ii);
    fprintf('z = %f\n', z);
    
    % Beam before lens
    x = linspace(-10*w0, 10*w0, 2000);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
    F = exp(-r_mesh.^2/w0^2);
    
    imagesc(x, x.', abs(F));
    colorbar
    
    sum(sum(abs(F).^2))
    
    F = propLens(F, x, lambda, f);
    [F, x] = propFresnel(F, x, lambda, f);
    
    sum(sum(abs(F).^2))
    
    imagesc(x, x.', abs(F));
    
    F = F.*exp(1i*k0*n2*L*abs(F).^2);
    
    imagesc(angle(F));
    
    [F, x] = propFresnel(F, x, lambda, d);
    
    imagesc(x, x.', abs(F));
    
    r_mesh = sqrt(x.^2 + x.'.^2);
    DiodePower(ii) = sum(sum(abs(F(r_mesh < 0.1)).^2));
    
end

figure;
plot(arrz, DiodePower);

%% Helper functions
function F = propLens(F, x, lambda, f)
    % PROPLENS Propagation through a lens
    % Author : RPM
    % Based on Eq. 5-10 of Goodman
    % Inputs:
    %   - F : fields in source plane
    %   - x : coordinates in source plane
    %   - lambda : wavelength
    %   - f : lens focal length
    % Outputs:
    %   - F = fields in observation plane

    k = 2*pi/lambda;
    [x_mesh, y_mesh] = meshgrid(x, x);
    F = F .* exp(-1i*k/2/f*(x_mesh.^2 + y_mesh.^2));
    
end
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
