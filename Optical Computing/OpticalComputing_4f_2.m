% Ralf Mouthaan
% University of Cambridge
% January 2021
%
% Model of a typical 4f system. This script models Fresnel propagation and
% propagation through a lens separately.
%
% The system being modelled is:
% SMF --> Lens 1 --> SLM 1 --> Lens 2 --> SLM 2 --> Lens 3 --> Camera
% One focal length between each.
%
% This is a correlator, so as an example, a low-pass filter of the mandrill
% is implemented.

clc; clear variables; close all;
addpath('../Function Library');

%% User-entered parameters

lambda = 633e-9;
Nx = 1000;
f = 1;

% Set SLM1
load('mandrill', 'X');
SLM1 = zeros(Nx);
SLM1(Nx/2 - size(X,1)/2+1:Nx/2 + size(X,1)/2, ...
     Nx/2 - size(X,2)/2+1:Nx/2 + size(X,2)/2) = X;

% Set SLM2
SLM2 = zeros(Nx);
SLM2(Nx/2-round(Nx/25):Nx/2+round(Nx/25),Nx/2-round(Nx/25):Nx/2+round(Nx/25)) = 1;

%% Calculations

x = linspace(-1000e-5, 1000e-5, Nx);
SMF = CreateSMF(x);
F = SMF.F; % This is the beam out of the SMF

[F, x] = propFresnel(F, x, lambda, f);
F = propLens(F, x, lambda, f);
[F,x] = propFresnel(F, x, lambda, f); % Illumnation of SLM1

F = F.*SLM1; % Diffraction field of SLM1

[F, x] = propFresnel(F, x, lambda, f);
F = propLens(F, x, lambda, f);
[F,x] = propFresnel(F, x, lambda, f); % Illumination of SLM2

F = F.*SLM2; % Diffraction field of SLM2

[F, x] = propFresnel(F, x, lambda, f);
F = propLens(F, x, lambda, f);
[F,x] = propFresnel(F, x, lambda, f); % Incident on camera

%% Plot Reults

% Plot SLM 1
figure; 
imagesc(SLM1);
axis square;
title('SLM1');
xticks(''); yticks('');
colormap gray;

% Plot SLM2
figure;
imagesc(SLM2);
axis square;
title('SLM2');
xticks(''); yticks('');
colormap gray

% Plot camera
F = flipud(F); % Seems to be flipped. Not sure if this is a bug or real.
figure; 
imagesc(x*1e3, x.'*1e3, abs(F));
axis square;
xlabel('mm'); ylabel('mm');
title('Camera');
colormap gray;