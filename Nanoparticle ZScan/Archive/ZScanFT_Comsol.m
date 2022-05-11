% Ralf Mouthaan
% University of Cambridge
% June 2020
% 
% Modelling of z-scan system.
% Uses Fourier transforms for propagation.
% This approach tries to use Fourier transforms for all propagations.

clc; clear variables; close all;

%% User-defined parameters

n0 = 1.45; % Linear refractive index
n2 = 1e-5; % Non-linear refractive index
lambda = 532e-9; % Green 
w0 = 8.04e-6; % beam waist at focus.
z2 = 10; % Distance from sample to diode. THIS HAS BEEN ADJUSTED TO BE FRESNEL REGION
ra = 0.01; % Aperture size. % SMALL SEEMS TO WORK BEST
L = 300e-6; % Path length through cuvette

%% Derived parameters

k0 = 2*pi/lambda;
z0= k0*w0^2/2; % diffraction length of beam
arrz = linspace(-10*z0,10*z0,50);

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
    
    % Field before sample
    if z==0
        F = w0/wz*exp(-r_mesh.^2/wz^2);
    else
        F = w0/wz*exp(-r_mesh.^2/wz^2 - 1i*k0*r_mesh.^2/2/Rz);
    end
    
    % Phase profile due to sample
    RI = getComsolRI(r_mesh, z);
    dphase = RI*L/lambda*2*pi;
    F = F.*exp(-1i*dphase);
    
    % Propagation
    [F, x] = propFresnel(F, x, lambda, d);
    
    % Sum over aperture
    r_mesh = sqrt(x.^2 + x.'.^2);
    T(ii) = sum(sum(abs(F(r_mesh < ra)).^2));
    
end

xlabel('radius (\mum)');
ylabel('RI');

figure;
plot(arrz/z0, T/T0);
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
function RI = getComsolRI(rmesh, zinterp)

    % Omid has obtained the resultant refractive indices using Comsol.
    % These are saved to a file. This function loads these, and performs
    % some (a lot) of interpolation to get the refractive indices in the
    % right coordinate system.
    
    %% Load data from file
    RI = dlmread('zscan_ri.txt');
    RI = RI(2:end,2:end);
    z = (0:size(RI, 1)-1)*1e-3;
    r = (0:size(RI,2)-1);
    r = r-max(r)/2;
    r = r*0.5e-6;
    
    %% 1D interpolation
    
    rmax = max(r);
    dr = 0.1e-6;
    rmax = round(rmax/dr)*dr;
    rnew = (0:rmax/dr)*dr;

    for i = 1:size(RI, 1)
        RInew(i,:) = interp1(r, RI(i,:), rnew);
    end

    r = rnew;
    RI = RInew;
    clear RInew rnew
    
    %% 2D interpolation

    rinds = rmesh/dr;
    rinds = round(rinds);
    rinds(rinds > rmax/dr) = rmax/dr;

    for i = 1:size(RI, 1)

        RIslice = RI(i, rinds);
        RIslice = reshape(RIslice, size(rmesh, 1), size(rmesh, 2));

        RInew(:,:,i) = RIslice;

    end

    RI = RInew;
    
    %% 3D interpolation
    
    x = 1:size(rmesh, 1);
    RI = interp3(x, x.', z, RI, x, x.', abs(zinterp));

end
