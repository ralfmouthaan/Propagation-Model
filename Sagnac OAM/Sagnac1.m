% Ralf Mouthaan
% Adelaide University
% November 2026
%
% Debugging Katelyn's sagnac OAM system.
% Katelyn can't seem to get proper OAM beams using her misaligned Sagnac
% system. This script is using an ASM propagation model to investigate what
% the effect of misaligning the mirrors should do in theory.

clc; clear variables; close all;
addpath('../_Function Library')

%%

lambda = 633e-9;
Nx = 1024;
w0 = 1e-3;
d = 10e-2; % Side length of square sagnac

x = linspace(-5e-3, 5e-3, Nx);
[x_mesh, y_mesh] = meshgrid(x, x.');
r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);

M1_tiltx = 0; M1_tilty = 50;
M2_tiltx = 2.5; M2_tilty = 0;
M3_tiltx = -M1_tiltx; M3_tilty = -M1_tilty;
Defocus = 0;

% Original beam
F0 = exp(-(r_mesh.^2)/w0^2).*exp(1i*r_mesh.^2*Defocus*500);

% Beam going round clockwise
F1 = F0;
F1 = propAngularSpectrum(F1, x, lambda, d);
F1 = F1.*exp(1i*x_mesh*M1_tiltx*500);
F1 = F1.*exp(1i*y_mesh*M1_tilty*500);
F1 = propAngularSpectrum(F1, x, lambda, d);
F1 = F1.*exp(1i*x_mesh*M2_tiltx*500);
F1 = F1.*exp(1i*y_mesh*M2_tilty*500);
F1 = propAngularSpectrum(F1, x, lambda, d);
F1 = F1.*exp(1i*x_mesh*M3_tiltx*500);
F1 = F1.*exp(1i*y_mesh*M3_tilty*500);
F1 = propAngularSpectrum(F1, x, lambda, d);

% Beam going round counter-clockwise
F2 = F0;
F2 = propAngularSpectrum(F2, x, lambda, d);
F2 = F2.*exp(1i*x_mesh*-M3_tiltx*500);
F2 = F2.*exp(1i*y_mesh*M3_tilty*500);
F2 = propAngularSpectrum(F2, x, lambda, d);
F2 = F2.*exp(1i*x_mesh*-M2_tiltx*500);
F2 = F2.*exp(1i*y_mesh*M2_tilty*500);
F2 = propAngularSpectrum(F2, x, lambda, d);
F2 = F2.*exp(1i*x_mesh*-M1_tiltx*500);
F2 = F2.*exp(1i*y_mesh*M1_tilty*500);
F2 = propAngularSpectrum(F2, x, lambda, d);

Fsum = F1 - F2; % Difference because Sagnac BS

%%

figure; 
subplot(1,2,1);
imagesc(x*1e3, x*1e3, abs(F1).^2);
xlabel('mm'); ylabel('mm');
axis square;
title('Intensity');
subplot(1,2,2);
imagesc(x*1e3, x*1e3, angle(F1));
xlabel('mm'); ylabel('mm');
axis square;
title('Phase');

figure; 
subplot(1,2,1);
imagesc(x*1e3, x*1e3, abs(F2).^2);
xlabel('mm'); ylabel('mm');
axis square;
title('Intensity');
subplot(1,2,2);
imagesc(x*1e3, x*1e3, angle(F2));
xlabel('mm'); ylabel('mm');
axis square;
title('Phase');

figure; 
subplot(1,2,1);
imagesc(x*1e3, x*1e3, abs(Fsum).^2);
xlabel('mm'); ylabel('mm');
axis square;
title('Intensity');
subplot(1,2,2);
imagesc(x*1e3, x*1e3, angle(Fsum));
xlabel('mm'); ylabel('mm');
axis square;
title('Phase');

