% Ralf Mouthaan
% University of Cambridge
% 
% Script to plot aberrated point spread functions

clc; clear variables; close all;

addpath('Function Library');

f = 500e-3;
z1 = 1;
z2 = f*z1/(z1-f);
lambda = 633e-9;
Nx = 5000;

x = linspace(-0.5e-3,0.5e-3, Nx);
r_mesh = sqrt(x.^2 + x.'.^2);
%F = exp(-r_mesh.^2/1e-10);
F = SpadeTarget(Nx);

figure;
imagesc(x, x.', abs(F));
axis square;
title('Object');

[F,x] = propFresnel(F,x,lambda,z1);
F = propLens(F,x,lambda,f);

[x_mesh, y_mesh] = meshgrid(x, x.');
r_mesh = sqrt(x_mesh.^2 + y_mesh.^2);
theta_mesh = atan2(y_mesh, x_mesh);
Z = Zernike(2,2,r_mesh,theta_mesh);
F = F.*exp(1i*Z.F*10000);
Z = Zernike(2,4,r_mesh,theta_mesh);
F = F.*exp(1i*Z.F*10000);
Z = Zernike(0,4,r_mesh,theta_mesh);
F = F.*exp(1i*Z.F*5000);
 F(r_mesh > 0.01) = 0;

figure;
imagesc(x, x.', abs(F));
axis square;
title('Lens Plane');

[F,x] = propFresnel(F,x,lambda,z2);
F = rot90(F, 2);

figure;
imagesc(x, x.', abs(F));
axis square;
set(gca, 'FontSize', 14);
xticks('');
yticks('');

function [Target] = DotTarget(Nx)
    
    Target = zeros(Nx, Nx);
    Target(Nx/2:Nx/2+1,Nx/2:Nx/2+1) = 8;

end
function [Target] = SmileyTarget(Nx)

    Target = imread('SmileyFace.jpg');
    Target = rgb2gray(Target);
    Target = double(Target);
    
    x = linspace(0, 1, size(Target, 1));
    xi = linspace(0, 1, Nx);    
    
    Target = interp2(x, x.', Target, xi, xi.');
    Target = imbinarize(Target, 100);
    Target = ~Target;
    Target = Target*8;

end
function [Target] = CrossTarget(Nx)

    Target = zeros(Nx, Nx);
    
    Target(Nx/10:Nx-Nx/10,Nx/2-Nx/20:Nx/2+Nx/20) = 8;
    Target(Nx/2-Nx/20:Nx/2+Nx/20,Nx/10:Nx-Nx/10) = 8;

end
function [Target] = HeartTarget(Nx)

    xi = linspace(0, 1, Nx);

    Target = imread('CardSuits.png');
    Target = rgb2gray(Target);
    Target = Target(1:size(Target,1)/2,size(Target,2)/2+1:size(Target,2));
    Target = double(Target);
    
    x = linspace(0, 1, size(Target, 1));
    
    Target = interp2(x, x.', Target, xi, xi.');
    Target = imbinarize(Target, 100);
    Target = ~Target;
    Target = Target*100;

end
function [Target] = SpadeTarget(Nx)

    xi = linspace(0, 1, Nx);

    Target = imread('CardSuits.png');
    Target = rgb2gray(Target);
    Target = Target(1:size(Target,1)/2,1:size(Target,1)/2);
    Target = double(Target);
    
    x = linspace(0, 1, size(Target, 1));
    
    Target = interp2(x, x.', Target, linspace(0, 0.8, Nx), xi.');
    Target = imbinarize(Target, 100);
    Target = ~Target;
    Target = Target*100;

end
function [Target] = DiamondTarget(Nx)

    xi = linspace(0, 1, Nx);

    Target = imread('CardSuits.png');
    Target = rgb2gray(Target);
    Target = Target(size(Target,2)/2+1:size(Target,2),1:size(Target,1)/2);
    Target = double(Target);
    
    x = linspace(0, 1, size(Target, 1));
    
    Target = interp2(x, x.', Target, linspace(0, 0.8, Nx), xi');
    Target = imbinarize(Target, 100);
    Target = ~Target;
    Target = Target*100;

end
function [Target] = ClubTarget(Nx)

    xi = linspace(0, 1, Nx);

    Target = imread('CardSuits.png');
    Target = rgb2gray(Target);
    Target = Target(size(Target,2)/2+1:size(Target,2),size(Target,2)/2+1:size(Target,2));
    Target = double(Target);
    
    x = linspace(0, 1, size(Target, 1));
    
    Target = interp2(x, x.', Target, xi, xi.');
    Target = imbinarize(Target, 100);
    Target = ~Target;
    Target = Target*100;

end