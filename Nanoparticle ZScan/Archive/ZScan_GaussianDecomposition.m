% Ralf Mouthaan
% University of Cambridge
% June 2020
%
% Modelling of z-scan system.
% Uses Gaussian decomposition for propagation.
% Assumes a Gaussian beam throughout.
% This is the approach originally used by Shak-Bahae.
%
% Reference:
% Shak-Bahae et al. "Z-Scan: A Simple And Sensitive Technique For 
% Nonlinear Refraction Measurements", Proc. SPIE, Nonlinear Optical 
% Properties of Materials, 1989.

clc; clear variables; close all;

%% User-defined parameters

% n = n0 + n2*|E|^2

n0 = 1.45; % Linear refractive index
n2 = 0.5e-6; % Non-linear refractive index
lambda = 532e-9; % Green DPSS
w0 = 8.04e-6; % width of beam before lens
z0 = 511e-6;
z2 = 0.5; % Distance from sample to diode. FAR FIELD CONDITION?
ra = 0.001; % Aperture size
L = 300e-6; % Path length through cuvette

%% Derived parameters

k0 = 2*pi/lambda;
rsample = linspace(0, 4*w0, 1000);
rreplay = linspace(0, 5e-2, 1000); % SHOULD BE ABLE TO CALCULATE AN APPROPRIATE RANGE
arrz = linspace(-5*z0,5*z0,100);

%% Error Checks

% SOME ERROR CHECKS


%% Calculations

% With Sample
figure;
T  = zeros(size(arrz));
for ii = 1:length(arrz)
    
    z = arrz(ii);
    Rz = z*(1+z0^2/z^2);
    wz = w0*sqrt(1+z^2/z0^2);
    d = z2 - z;
    g=1+d/Rz;
    
    % Field before sample
    E = w0/wz*exp(-rsample.^2/wz^2 - 1i*k0*rsample.^2/2/Rz);
    
    % Exponential phase profile imposed by sample (like Sheik-Bahae)
    % Can remove the 1+z^2/z0^2 term, as we have already compensated for
    % the fact we are not at the focal point in our calculation of E (see
    % the line before this one).
    dphase0 = k0*n2*L*abs(E(1)).^2; %/(1+z^2/z0^2); 
   %dphase = dphase0 * exp(-2*rsample.^2/wz^2);
    dphase = k0*n2*L*abs(E).^2; 
    
%     % Gaussian decomposition approximation to phase profile imposed by sample
%     % Really they mean a Taylor series decomposition
%     dphaseapprox = zeros(size(rsample));
%     for m = 0:100
%         dphaseapprox = dphaseapprox + (1i*dphase0).^m/factorial(m) .* exp(-2*m*rsample.^2/wz^2);
%     end
%     
%     figure;
%     plot(rsample, imag(exp(1i*dphase)), 'b-');
%     hold on
%     plot(rsample, imag(dphaseapprox), 'rx');

    % Field with phase profile imposed
    E = E.*exp(1i*dphase);
    
    % Field profile at diode without sample phase profile imposed
    % CAN BE MOVED OUT OF FOR LOOP, AS I DON'T THINK THIS CHANGES
    wmo = wz;
    dm = k0*wmo^2/2;
    wm = wmo*sqrt(g^2 + d^2/dm^2);
    Rm = d/(1-g/(g^2 + d^2/dm^2));
    thetam = atan2(d/dm, g); % Does this pick up the right one?
    EReplay_NoSample =  (g^2 + d^2/dm^2)^(-1/2) .* exp( - rreplay.^2/wm^2 - 1i*k0*rreplay.^2/2/Rm + 1i*thetam);
    EReplay_NoSample = EReplay_NoSample*E(1);
    sum(EReplay_NoSample)
    S = sum(rreplay(rreplay<ra).*abs(EReplay_NoSample(rreplay<ra)).^2) / ...
        sum(rreplay.*abs(EReplay_NoSample).^2);
    
    % Field profile at diode, with sample phase profile imposed
    % Calculated by using Gaussian decomposition of field after sample,
    % propagating each Gaussian to replay field plane, and recombining there.
    EReplay_Sample = zeros(size(rreplay));
    for m = 0:10
        
        wmo = wz/sqrt(2*m+1);
        dm = k0*wmo^2/2;
        wm = wmo*sqrt(g^2 + d^2/dm^2);
        Rm = d/(1-g/(g^2 + d^2/dm^2));
        thetam = atan2(d/dm, g); % Does this pick up the right one?
        
        EReplay_Sample = EReplay_Sample + (1i*dphase0).^m/factorial(m) * wmo/wm .* exp( - rreplay.^2/wm^2 - 1i*k0*rreplay.^2/2/Rm + 1i*thetam);
        
    end   
    EReplay_Sample = EReplay_Sample*E(1);
    
    plot(rreplay, abs(EReplay_Sample));
    hold on
    
    T(ii) = sum(rreplay(rreplay<ra).*abs(EReplay_Sample(rreplay<ra)).^2);
    T(ii) = T(ii) / sum(rreplay(rreplay<ra).*abs(EReplay_NoSample(rreplay<ra)).^2);
    
end

figure;
plot(arrz/z0, T);
xlabel('z/z_0');
ylabel('T/T_0');

fprintf('Tmax - Tmin = %f\n', max(T) - min(T));
DeltaPhi0 = 0.406*(1-S)^0.25*(max(T) - min(T));
n2 = DeltaPhi0/L/k0;
fprintf('n2 = %e\n', n2);
