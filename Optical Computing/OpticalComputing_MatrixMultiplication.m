% Ralf Mouthaan
% University of Cambridge
% January 2021
%
% Model of a typical matrix multiplication system. This script models Fresnel propagation and
% propagation through a lens separately.
%
% The system being modelled is:
% SMF --> Lens 1 --> SLM 1 --> Lens 2 --> SLM 2 --> Lens 3 --> Camera
% Two focal lengths between each
%
% Random matrices are used as examples

clc; clear variables; close all;
addpath('../Function Library');

%% User-entered parameters

lambda = 633e-9;
Nx = 1000;
f = 1;

% Set SLM1
M1 = rand(10);
SLM1 = zeros(Nx);
SLM1(Nx/2 - size(M1,1)/2+1:Nx/2 + size(M1,1)/2, ...
     Nx/2 - size(M1,2)/2+1:Nx/2 + size(M1,2)/2) = M1;

% Set SLM2
SLM2 = zeros(Nx);
M2 = rand(10);
SLM2(Nx/2 - size(M2,1)/2:Nx/2 + size(M2,1)/2-1, ...
     Nx/2 - size(M2,2)/2:Nx/2 + size(M2,2)/2-1) = M2; % We have some kind of off-by-one error that I don't understand but am correcting here.

%% Calculations

x = linspace(-1000e-5, 1000e-5, Nx);
SMF = CreateSMF(x);
F = SMF.F; % This is the beam out of the SMF

[F, x] = propFresnel(F, x, lambda, f);
F = propLens(F, x, lambda, f);
[F,x] = propFresnel(F, x, lambda, f); % Illumnation of SLM1

F = F.*SLM1; % Diffraction field of SLM1

[F, x] = propFresnel(F, x, lambda, 2*f);
F = propLens(F, x, lambda, f);
[F,x] = propFresnel(F, x, lambda, 2*f); % Illumination of SLM2

F = F.*flipud(fliplr(SLM2)); % Image is flipped, so flip matrix by which we're multiplying.
                             % Diffraction field of SLM2

[F, x] = propFresnel(F, x, lambda, 2*f);
F = propLens(F, x, lambda, f);
[F,x] = propFresnel(F, x, lambda, 2*f); % Incident on camera

%% Plot Reults

figure; imagesc(M1.*M2); axis square; title('Expected Result');

M1M2 = F(Nx/2 - size(M1,1)/2+1:Nx/2 + size(M1,1)/2, ...
     Nx/2 - size(M1,2)/2+1:Nx/2 + size(M1,2)/2);
 
 figure; imagesc(abs(M1M2)); axis square; title('Observed Result');