% Ralf Mouthaan
% University of Cambridge
% May 2021
%
% Checking my heat equation maths

clc; clear variables; close all;

%% User-defined

c = 4200; %J/kg/C (Water)
rho = 997; %kg/m3 (Water)
k = 0.6089; %W/m/C (Water)
D = k/rho/c; % Water, dunno about the units.
sigma = 1e-6;

%% Derived values

A = rho*c;
r = linspace(-20e-3,20e-3, 1000);
arrt = linspace(0, 3, 10);

%% Heat Source

f = A/rho/c * exp(-r.^2/2/sigma^2);

figure;
plot(r*1e3, f, 'linewidth', 3);
xlabel('mm');
ylabel('\DeltaT (C)');

arrlegend{1} = 'Heat Source';

%% Temperature curves

hold on
for ii = 1:length(arrt)
    
    t = arrt(ii);
    phi = -A*sigma^2/2/k * ( expint(r.^2/2/sigma^2) - expint(r.^2/(2*sigma^2 + 4*D*t)));
    plot(r*1e3, phi);
    
    arrlegend{ii+1} = sprintf('t = %0.2fs', t);

end

legend(arrlegend);