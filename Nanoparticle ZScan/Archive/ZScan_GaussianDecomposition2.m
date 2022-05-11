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
n2 = 1e-5; % Non-linear refractive index
lambda = 532e-9; % Green 
w0 = 8.04e-6; % beam waist at focus.
z2 = 0.5; % Distance from sample to diode. THIS IS NOT FAR-FIELD
ra = 0.0005; % Aperture size. % SMALL SEEMS TO WORK BEST
L = 300e-6; % Path length through cuvette

%% Derived parameters

k0 = 2*pi/lambda;
z0 = k0*w0^2/2; %511e-6; % THIS ANALYSIS DOES NOT WORK IF THIS RELATIONSHIP DOES NOT HOLD
rsample = linspace(0, 4*w0, 1000);
rreplay = linspace(0, 5e-2, 1000); % SHOULD BE ABLE TO CALCULATE AN APPROPRIATE RANGE
arrz =  linspace(-5*z0,5*z0,100);

%% Calculations

% With Sample
figure;
T  = zeros(size(arrz));
T0 = zeros(size(arrz));
for ii = 1:length(arrz)
    
    z = arrz(ii);
    Rz = z*(1+z0^2/z^2);
    wz = w0*sqrt(1+z^2/z0^2);
    d = z2 - z;
    g=1+d/Rz;
    
    % Field before sample at plane z
    if abs(z) == 0
        E = w0/wz*exp(-rsample.^2/wz^2);
    else
        E = w0/wz*exp(-rsample.^2/wz^2 - 1i*k0*rsample.^2/2/Rz);
    end
    
    % Exponential phase profile imposed by sample
    % Can remove the 1+z^2/z0^2 term present in Sheik-Bahae, as E is not E
    % at the focus. Instead, it is E at the plane z.
    dphase = k0*n2*L*abs(E).^2;
    
    % Field profile at diode without sample phase profile imposed
    % CAN BE MOVED OUT OF FOR LOOP, AS I DON'T THINK THIS CHANGES
    wmo = wz;
    dm = k0*wmo^2/2;
    wm = wmo*sqrt(g^2 + d^2/dm^2);
    Rm = d/(1-g/(g^2 + d^2/dm^2));
    thetam = atan2(d/dm, g); % Does this pick up the right one?
    EReplay_NoSample =  (g^2 + d^2/dm^2)^(-1/2) .* exp( - rreplay.^2/wm^2 - 1i*k0*rreplay.^2/2/Rm + 1i*thetam);
    EReplay_NoSample = EReplay_NoSample*E(1);
    S = sum(rreplay(rreplay<ra).*abs(EReplay_NoSample(rreplay<ra)).^2) / ...
        sum(rreplay.*abs(EReplay_NoSample).^2);
    T0(ii) = sum(rreplay(rreplay<ra).*abs(EReplay_NoSample(rreplay<ra)).^2);
    
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
        
        EReplay_Sample = EReplay_Sample + (1i*dphase(1)).^m/factorial(m) * wmo/wm .* exp( - rreplay.^2/wm^2 - 1i*k0*rreplay.^2/2/Rm + 1i*thetam);
        
    end
    EReplay_Sample = EReplay_Sample*E(1);
    
    plot(rreplay, abs(EReplay_Sample));
    hold on
    
    T(ii) = sum(rreplay(rreplay<ra).*abs(EReplay_Sample(rreplay<ra)).^2);
    
end

figure;

figure;
plot(arrz/z0, T);
hold on
plot(arrz/z0, T0);
xlabel('z/z_0');
ylabel('T/T_0');

fprintf('Tmax - Tmin = %f\n', max(T./T0/S) - min(T./T0/S));
DeltaPhi0 = (max(T./T0) - min(T./T0))/0.406/(1-S)^0.25;
fprintf('Calculated Delta Phi 0 %e\n', DeltaPhi0);
n2_calculated = DeltaPhi0/L/k0;
fprintf('n2 set = %e\n', n2);
fprintf('n2 calculated = %e\n', n2_calculated);
