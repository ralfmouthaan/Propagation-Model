% Ralf Mouthaan
% Adelaide University
% November 2026
%
% Debugging Katelyn's sagnac OAM system.
% Katelyn can't seem to get proper OAM beams using her misaligned Sagnac
% system. This script is using an ASM propagation model to investigate what
% the effect of misaligning the mirrors should do in theory.
%
% This version optimises M1_tilty and M2_tiltx to give the most uniform
% intensity around the ring of the output LG beam.

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

M1_tiltx = 0; M1_tilty = 2;     % M1_tilty is optimised
M2_tiltx = 0.1; M2_tilty = 0;   % M2_tiltx is optimised
Defocus = 0;

% Original beam
F0 = exp(-(r_mesh.^2)/w0^2).*exp(1i*r_mesh.^2*Defocus*500);

%% Optimisation

% Cost = RMS non-uniformity of intensity around the ring (std/mean)
CostFn = @(p) RingCost(SimulateSagnac(F0, x, x_mesh, y_mesh, lambda, d, ...
    [M1_tiltx p(1) p(2) M2_tilty]), x, x_mesh, y_mesh);

p_init = [M1_tilty M2_tiltx];
options = optimset('Display', 'iter', 'TolX', 1e-4, 'TolFun', 1e-4, 'MaxFunEvals', 300);
[p_opt, cost_opt] = fminsearch(CostFn, p_init, options);

fprintf('\nInitial: M1_tilty = %.4f, M2_tiltx = %.4f, cost = %.4f\n', ...
    p_init(1), p_init(2), CostFn(p_init));
fprintf('Optimal: M1_tilty = %.4f, M2_tiltx = %.4f, cost = %.4f\n\n', ...
    p_opt(1), p_opt(2), cost_opt);

M1_tilty = p_opt(1);
M2_tiltx = p_opt(2);

%% Simulate optimised system

Fsum = SimulateSagnac(F0, x, x_mesh, y_mesh, lambda, d, ...
    [M1_tiltx M1_tilty M2_tiltx M2_tilty]);

%% Ring uniformity analysis

R = AnalyseRing(Fsum, x, x_mesh, y_mesh);

fprintf('Output power (relative to input): %.4f\n', sum(abs(Fsum(:)).^2)/sum(abs(F0(:)).^2));
fprintf('Beam centre: (%.3f, %.3f) mm\n', R.x0*1e3, R.y0*1e3);
fprintf('Ring radius: %.3f mm\n', R.r_ring*1e3);
fprintf('Min/max intensity ratio: %.4f\n', R.Uniformity_minmax);
fprintf('Intensity visibility (max-min)/(max+min): %.4f\n', R.Visibility);
fprintf('RMS non-uniformity (std/mean): %.4f\n', R.RMS_nonuniformity);
fprintf('Topological charge: %.3f\n', R.l_charge);
fprintf('RMS phase deviation from linear ramp: %.4f rad\n', std(R.phase_dev));

%%

figure;
subplot(1,2,1);
imagesc(x*1e3, x*1e3, abs(Fsum).^2);
xlabel('mm'); ylabel('mm');
axis square;
title('Intensity');
hold on;
plot(R.x_ring*1e3, R.y_ring*1e3, 'w--', 'LineWidth', 1);
plot(R.x0*1e3, R.y0*1e3, 'w+');
hold off;
subplot(1,2,2);
imagesc(x*1e3, x*1e3, angle(Fsum));
xlabel('mm'); ylabel('mm');
axis square;
title('Phase');
hold on;
plot(R.x_ring*1e3, R.y_ring*1e3, 'w--', 'LineWidth', 1);
plot(R.x0*1e3, R.y0*1e3, 'w+');
hold off;

figure;
subplot(3,1,1);
plot(R.r_centres*1e3, R.I_radial);
xline(R.r_ring*1e3, 'r--');
xlabel('Radius (mm)'); ylabel('Mean intensity');
title('Azimuthally averaged radial profile');
subplot(3,1,2);
plot(R.theta*180/pi, R.I_ring_norm);
xlim([0 360]); ylim([0 max(1.5, max(R.I_ring_norm)*1.1)]);
xlabel('\theta (deg)'); ylabel('I / <I>');
title(sprintf('Intensity around ring (min/max = %.3f, std/mean = %.3f)', ...
    R.Uniformity_minmax, R.RMS_nonuniformity));
subplot(3,1,3);
plot(R.theta*180/pi, R.phase_dev);
xlim([0 360]);
xlabel('\theta (deg)'); ylabel('Phase deviation (rad)');
title(sprintf('Phase around ring minus linear ramp (l = %.2f)', R.l_charge));

%% Functions

function Fsum = SimulateSagnac(F0, x, x_mesh, y_mesh, lambda, d, tilts)

    M1_tiltx = tilts(1); M1_tilty = tilts(2);
    M2_tiltx = tilts(3); M2_tilty = tilts(4);
    M3_tiltx = -M1_tiltx; M3_tilty = -M1_tilty;

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

end

function cost = RingCost(Fsum, x, x_mesh, y_mesh)

    % Penalise solutions with no vortex, more than one vortex in the beam
    % (at large tilts the arms stop overlapping and form a fringe pattern
    % with a row of vortices), or not a single-charge vortex
    R = AnalyseRing(Fsum, x, x_mesh, y_mesh);
    if isempty(R) || R.N_vortex ~= 1 || abs(round(R.l_charge)) ~= 1
        cost = 10;
    else
        cost = R.RMS_nonuniformity;
    end

end

function R = AnalyseRing(Fsum, x, x_mesh, y_mesh)

    I = abs(Fsum).^2;

    % Intensity centroid, used as a first guess of the beam position (tilts
    % can shift the beam off-axis)
    xc = sum(I(:).*x_mesh(:))/sum(I(:));
    yc = sum(I(:).*y_mesh(:))/sum(I(:));

    % Beam centre = phase singularity nearest the centroid. The centroid
    % itself is biased towards the brighter side of a non-uniform ring.
    % Singularities are found from the phase circulation around each 2x2
    % pixel plaquette. Only singularities within the beam are counted (the
    % phase is just numerical noise where there is no light).
    ph = angle(Fsum);
    wrap = @(a) angle(exp(1i*a));
    circ = wrap(ph(1:end-1,2:end) - ph(1:end-1,1:end-1)) + ...
           wrap(ph(2:end,2:end)   - ph(1:end-1,2:end)) + ...
           wrap(ph(2:end,1:end-1) - ph(2:end,2:end)) + ...
           wrap(ph(1:end-1,1:end-1) - ph(2:end,1:end-1));
    r_c = sqrt((x_mesh - xc).^2 + (y_mesh - yc).^2);
    [r_sorted, idx_sort] = sort(r_c(:));
    r_99 = r_sorted(find(cumsum(I(idx_sort)) >= 0.99*sum(I(:)), 1));
    in_beam = r_c(1:end-1,1:end-1) < r_99;
    [row_v, col_v] = find(abs(round(circ/(2*pi))) > 0 & in_beam);
    if isempty(row_v)
        R = [];
        return;
    end
    R.N_vortex = length(row_v);
    x_v = (x(col_v) + x(col_v+1))/2;
    y_v = (x(row_v) + x(row_v+1))/2;
    [~, idx_v] = min((x_v - xc).^2 + (y_v - yc).^2);
    R.x0 = x_v(idx_v);
    R.y0 = y_v(idx_v);

    % Ring radius = peak of azimuthally averaged radial intensity profile
    dx = x(2) - x(1);
    r_rel = sqrt((x_mesh - R.x0).^2 + (y_mesh - R.y0).^2);
    r_bins = 0:dx:max(x)-max(abs([R.x0 R.y0]));
    r_idx = discretize(r_rel(:), r_bins);
    valid = ~isnan(r_idx);
    R.I_radial = accumarray(r_idx(valid), I(valid), [length(r_bins)-1 1], @mean);
    R.r_centres = r_bins(1:end-1) + dx/2;
    [~, idx_max] = max(R.I_radial);
    R.r_ring = R.r_centres(idx_max);

    % Sample intensity and phase around the ring
    N_theta = 720;
    R.theta = linspace(0, 2*pi, N_theta+1); R.theta = R.theta(1:end-1);
    R.x_ring = R.x0 + R.r_ring*cos(R.theta);
    R.y_ring = R.y0 + R.r_ring*sin(R.theta);
    I_ring = interp2(x_mesh, y_mesh, I, R.x_ring, R.y_ring, 'linear');
    F_ring = interp2(x_mesh, y_mesh, Fsum, R.x_ring, R.y_ring, 'linear');
    phi_ring = unwrap(angle(F_ring));

    % Uniformity metrics
    R.I_ring_norm = I_ring/mean(I_ring);
    R.Uniformity_minmax = min(I_ring)/max(I_ring);                         % 1 = perfectly uniform
    R.Visibility = (max(I_ring) - min(I_ring))/(max(I_ring) + min(I_ring)); % 0 = perfectly uniform
    R.RMS_nonuniformity = std(I_ring)/mean(I_ring);                         % 0 = perfectly uniform

    % Topological charge from phase winding, and deviation from linear
    % phase ramp
    R.l_charge = (phi_ring(end) - phi_ring(1) + angle(F_ring(1)/F_ring(end)))/(2*pi);
    p = polyfit(R.theta, phi_ring, 1);
    R.phase_dev = phi_ring - polyval(p, R.theta);

end
