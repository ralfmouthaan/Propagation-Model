% Ralf Mouthaan
% University of Cambridge
% May 2021
% 
% Modelling of z-scan system.
% Uses Fourier transforms for propagation.
% This approach tries to use Fourier transforms for all propagations.

clc; clear variables; close all;

%% User-defined parameters

n0 = 1.45; % Linear refractive index
n2 = 1e-3; % Non-linear refractive index
lambda = 532e-9; % Green 
w0 = 8.04e-6; % beam waist at focus.
z2 = 10; % Distance from sample to diode. THIS HAS BEEN ADJUSTED TO BE FRESNEL REGION
ra = 0.005; % Aperture size. % SMALL SEEMS TO WORK BEST
L = 30e-6; % Path length through cuvette

c = 4200; %J/kg/C (Water)
rho = 997; %kg/m3 (Water)
k = 0.6089; %W/m/C (Water)
t = 0.0001;

%% Derived parameters

k0 = 2*pi/lambda;
z0= k0*w0^2/2; % diffraction length of beam
arrz = linspace(-10*z0,10*z0,50);
D = k/rho/c; % Diffusion constant

%% Calculations without sample

x = linspace(-50*w0, 50*w0, 2000);
r = sqrt(x.^2 + x.'.^2);
F = exp(-r.^2/w0^2);
[F, x] = propFresnel(F, x, lambda, z2);
r = sqrt(x.^2 + x.'.^2);
S = sum(sum(r(r<ra).*abs(F(r<ra)).^2)) / sum(sum(r.*abs(F).^2));
T0 = sum(sum(abs(F(r < ra)).^2));

%% Calculations with sample

T = zeros(size(arrz));
for ii = 1:length(arrz)
    
    % z coordinate
    z = arrz(ii);
    d = z2 - z; % END UP WITH A SUPERIMPOSED SLOPE IF I USE D = Z2 - Z
    fprintf('z = %f\n', z);
    
    % Transverse coordinates
    x = linspace(-50*w0, 50*w0, 2000);
    [x_mesh, y_mesh] = meshgrid(x, x.');
    r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
    
    % Some derived parameters
    Rz = z*(1+z0^2/z^2);
    wz = w0*sqrt(1+z^2/z0^2);
    A = w0/wz*rho*c*1000;
    sigma = wz/2; % Ignore phase component - heating doesn't care.
    
    if z==0
        F = w0/wz*exp(-r_mesh.^2/wz^2);
    else
        F = w0/wz*exp(-r_mesh.^2/wz^2 - 1i*k0*r_mesh.^2/2/Rz);
    end
    
    % Beam profile without sample
    F0 = F; 
    
    % Heating profile
    F_Heating = -A*sigma^2/2/k * ( expint(r_mesh.^2/2/sigma^2) - expint(r_mesh.^2/(2*sigma^2 + 4*D*t)));
    F_Heating(r_mesh == 0) = -A*sigma^2/2/k * log(2*sigma^2/(2*sigma^2 + 4*D*t));
    
    % Phase profile with sample
    dphase = n2*F_Heating*L/lambda*2*pi;
    F = F.*exp(-1i*dphase);
    
    imagesc(dphase);
    colorbar;
    
    % Propagation
    [F, ~] = propFresnel(F, x, lambda, d);
    [F0, x] = propFresnel(F0, x, lambda, d);
    
    % Sum over aperture
    r_mesh = sqrt(x.^2 + x.'.^2);
    T(ii) = sum(sum(abs(F(r_mesh < ra)).^2));
    T0(ii) = sum(sum(abs(F0(r_mesh < ra)).^2));
    
end

xlabel('radius (\mum)');
ylabel('RI');

figure;
plot(arrz/z0, T./T0);
xlabel('z/z_0');
ylabel('T/T_0');

DeltaPhi0 = (max(T./T0) - min(T./T0))/0.406/(1-S)^0.25;
n2_calculated = DeltaPhi0/L/k0;
fprintf('n2 assigned = %e\n', n2);
fprintf('Delta Phi0 calculated = %e\n', DeltaPhi0);
fprintf('n2 calculated = %e\n', n2_calculated);

%% Helper functions
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
