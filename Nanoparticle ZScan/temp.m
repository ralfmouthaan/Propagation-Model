
clc; clear variables; close all;

%% Load data
RI = dlmread('zscan_ri.txt');
RI = RI(2:end,2:end);
z = (0:size(RI, 1)-1)*1e-3;
r = (0:size(RI,2)-1);
r = r-max(r)/2;
r = r*0.5e-6;

%% Plot of raw data

figure;
for i = 1:size(RI, 1)
    plot(r*1e6, RI(i,:));
    hold on
end

xlabel('radius (\mum)')
ylabel('RI');

%% 1D interpolation

rmax = max(r);
dr = 0.1e-6;
rnew = (0:rmax/dr)*dr;

for i = 1:size(RI, 1)
    RInew(i,:) = interp1(r, RI(i,:), rnew);
end

r = rnew;
RI = RInew;
clear RInew rnew

for i = 1:size(RI, 1)
    plot(r*1e6, RI(i,:), 'x');
    hold on
end

xlabel('radius (\mum)')
ylabel('RI');

%% 2D interpolation

x = linspace(-rmax, rmax, 1000);
rinds = sqrt(x.^2 + x.'.^2);
rinds = rinds/dr;
rinds = round(rinds);
rinds(rinds > rmax/dr) = rmax/dr;

for i = 1:size(RI, 1)

    RIslice = RI(i, rinds);
    RIslice = reshape(RIslice, 1000, 1000);
    
    RInew(:,:,i) = RIslice;
    
    figure;
    imagesc(RIslice);
    axis square;
    colorbar
    title(z(i));
    
end

RI = RInew;

%% 3D interpolation

figure;

for zinterp = 0:0.5e-3:4e-3
    RIslice = interp3(x, x.', z, RI, x, x.', zinterp);
    surf(x*1e6, x.'*1e6, RIslice, 'edgecolor', 'none');
    zlim([1.3262 1.3274]);
    title([num2str(zinterp*1000) 'mm']);
    xlabel('\mum');
    ylabel('\mum');
    zlabel('RI');
    colormap jet
    pause(1)
end








