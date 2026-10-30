function postprocess_cylinder_Re100()
%POSTPROCESS_CYLINDER_RE100  post-processing of main_cylinder_Re100 results
%
%   Reads  results/AerodynamicForces_Re100.txt  and  results/flowfield_Re100.mat
%   and produces:
%     - C_D(t), C_L(t) on [0,15] and on the paper window [14.5,15]
%     - the 26 uniformly spaced samples of C_D, C_L in [14.5,15] (paper's data)
%     - periodicity check: lift zero up-crossings -> period, frequency,
%       Strouhal number, cycle-to-cycle spread of the peaks, FFT, phase plot
%     - comparison with the DFG 2D-2 benchmark (Schaefer & Turek 1996) and
%       the FeatFlow reference signal shipped with redbKIT
%     - flow fields at t = 15 (|v|, u_x, p, vorticity) and vorticity at four
%       phases of one shedding period (wake / Karman vortex street)
%   All figures (PNG) and a text summary are written to ./results .
%   Works in MATLAB and GNU Octave.

if ~exist('results', 'dir'), mkdir('results'); end
close all;

Re    = 100;
D     = 0.1;    Umean = 1.0;    rho = 1.0;
xc    = 0.2;    yc    = 0.2;    r   = 0.05;
t_win = [14.5 15.0];    % time window used in the paper
t_dev = [10.0 15.0];    % fully-developed window for the periodicity analysis

%% ========================================================================
%  1) Aerodynamic coefficients
%  ========================================================================
fid = fopen(fullfile('results', 'AerodynamicForces_Re100.txt'), 'r');
C   = textscan(fid, '%f %f %f %f', 'HeaderLines', 1);
fclose(fid);
t  = C{1};  CD = C{2};  CL = C{3};
% (forces are already scaled by 2/(rho*Umean^2*D) inside redbKIT)
dt = median(diff(t));

% ---- 26 samples on [14.5, 15] (N_T = 26 of the paper) -------------------
t26  = linspace(t_win(1), t_win(2), 26)';
CD26 = interp1(t, CD, t26);
CL26 = interp1(t, CL, t26);
fid  = fopen(fullfile('results', 'CD_CL_Re100_t14.5-15_26pts.csv'), 'w');
fprintf(fid, 't,C_D,C_L\n');
fprintf(fid, '%.4f,%.8f,%.8f\n', [t26 CD26 CL26]');
fclose(fid);

% ---- periodicity analysis on the fully developed regime -----------------
[tc, Tk, stat] = cycle_analysis(t, CD, CL, t_dev);
f_shed = 1/mean(Tk);
St     = f_shed*D/Umean;

% FFT (Hann window, zero padding)
[fL, PL] = spectrum(t, CL, t_dev);
[fD, PD] = spectrum(t, CD - mean(CD(t >= t_dev(1))), t_dev);
[~, iL]  = max(PL);   f_fft_L = fL(iL);
[~, iD]  = max(PD);   f_fft_D = fD(iD);

% ---- FeatFlow reference signal (shipped with redbKIT) -------------------
R      = load(fullfile('reference', 'bdforces_q2_lv6_dt4'));
t_ref  = R(:,2);  CD_ref = R(:,4);  CL_ref = R(:,5);
[tc_ref, Tk_ref, stat_ref] = cycle_analysis(t_ref, CD_ref, CL_ref, [t_ref(1) t_ref(end)]);
St_ref = D/Umean/mean(Tk_ref);

% ---- self-similarity: last cycle vs previous cycles ---------------------
% max |C(t) - C(t+T)| over the developed regime (T = mean period)
T_mean = mean(Tk);
tt     = (t_dev(1):dt/2:t_dev(2) - T_mean)';
dCL    = max(abs(interp1(t, CL, tt, 'spline') - interp1(t, CL, tt + T_mean, 'spline')));
dCD    = max(abs(interp1(t, CD, tt, 'spline') - interp1(t, CD, tt + T_mean, 'spline')));

%% ========================================================================
%  2) Summary
%  ========================================================================
S = {};
S{end+1} = '=================================================================';
S{end+1} = ' Flow past a cylinder, Re = 100  (redbKIT, P2-P1, BDF2)';
S{end+1} = '=================================================================';
S{end+1} = sprintf(' nu = %.4g,  Umean = %g,  D = %g,  dt = %g,  T = %g', ...
                   D*Umean/Re, Umean, D, dt, t(end));
if exist(fullfile('results', 'flowfield_Re100.mat'), 'file')
    tmp = load(fullfile('results', 'flowfield_Re100.mat'));
    FIELD = tmp.FIELD;
    S{end+1} = sprintf(' mesh: %d triangles, %d P2 nodes, %d free dofs; wall-clock %.1f min', ...
                       FIELD.numElem, FIELD.numNodes, FIELD.numDofs, FIELD.t_cpu/60);
else
    FIELD = [];
end
S{end+1} = '';
S{end+1} = sprintf(' Periodicity analysis on t in [%g, %g]: %d full lift cycles', ...
                   t_dev(1), t_dev(2), numel(Tk));
S{end+1} = sprintf('   period T          = %.5f  (cycle-to-cycle std %.2e)', mean(Tk), std(Tk));
S{end+1} = sprintf('   frequency f       = %.5f  (FFT: f_L = %.4f, f_D = %.4f = %.2f f_L)', ...
                   f_shed, f_fft_L, f_fft_D, f_fft_D/f_fft_L);
S{end+1} = sprintf('   max|C_L(t)-C_L(t+T)| = %.2e,  max|C_D(t)-C_D(t+T)| = %.2e', dCL, dCD);
S{end+1} = '';
S{end+1} = '                        present      FeatFlow ref     DFG 2D-2 range';
S{end+1} = sprintf('   C_D max             %8.4f      %8.4f        [3.2200, 3.2400]', stat.CDmax, stat_ref.CDmax);
S{end+1} = sprintf('   C_D min             %8.4f      %8.4f', stat.CDmin, stat_ref.CDmin);
S{end+1} = sprintf('   C_D mean            %8.4f      %8.4f', stat.CDmean, stat_ref.CDmean);
S{end+1} = sprintf('   C_L max             %8.4f      %8.4f        [0.9900, 1.0100]', stat.CLmax, stat_ref.CLmax);
S{end+1} = sprintf('   C_L min             %8.4f      %8.4f', stat.CLmin, stat_ref.CLmin);
S{end+1} = sprintf('   Strouhal St=fD/U    %8.4f      %8.4f        [0.2950, 0.3050]', St, St_ref);
S{end+1} = sprintf('   (peak spreads over cycles: C_D max %.1e, C_L max %.1e)', ...
                   stat.CDmax_std, stat.CLmax_std);
S{end+1} = '';
S{end+1} = ' 26 samples on [14.5,15] -> results/CD_CL_Re100_t14.5-15_26pts.csv';
S{end+1} = sprintf('   C_D in [%.4f, %.4f],  C_L in [%.4f, %.4f]', min(CD26), max(CD26), min(CL26), max(CL26));
S{end+1} = '=================================================================';
fid = fopen(fullfile('results', 'summary_Re100.txt'), 'w');
fprintf(fid, '%s\n', S{:});
fclose(fid);
fprintf('%s\n', S{:});

%% ========================================================================
%  3) Figures: coefficients
%  ========================================================================
h = figure('visible', 'off', 'position', [100 100 1000 650]);
subplot(2,1,1);
plot(t, CD, 'b-', 'linewidth', 1.2); hold on;
yl = [0 3.6]; fill_window(t_win, yl);
ylim(yl); xlim([0 t(end)]); grid on; xlabel('t [s]'); ylabel('C_D');
title('Re = 100: drag coefficient (shaded: t \in [14.5, 15])');
subplot(2,1,2);
plot(t, CL, 'r-', 'linewidth', 1.2); hold on;
yl = [-1.2 1.2]; fill_window(t_win, yl);
ylim(yl); xlim([0 t(end)]); grid on; xlabel('t [s]'); ylabel('C_L');
title('Re = 100: lift coefficient');
save_fig(h, 'fig1_CD_CL_history_0-15');

h = figure('visible', 'off', 'position', [100 100 1000 650]);
m = t >= t_win(1) - 1e-9 & t <= t_win(2) + 1e-9;
subplot(2,1,1);
plot(t(m), CD(m), 'b-', 'linewidth', 1.5); hold on;
plot(t26, CD26, 'ko', 'markersize', 5, 'markerfacecolor', 'k');
grid on; xlabel('t [s]'); ylabel('C_D'); xlim(t_win);
legend('C_D (\Delta t = 0.01)', '26 uniform samples', 'location', 'eastoutside');
title('Fully developed regime t \in [14.5, 15]');
subplot(2,1,2);
plot(t(m), CL(m), 'r-', 'linewidth', 1.5); hold on;
plot(t26, CL26, 'ko', 'markersize', 5, 'markerfacecolor', 'k');
grid on; xlabel('t [s]'); ylabel('C_L'); xlim(t_win);
legend('C_L (\Delta t = 0.01)', '26 uniform samples', 'location', 'eastoutside');
save_fig(h, 'fig2_CD_CL_window_14.5-15');

% phase-averaged cycles vs FeatFlow reference
h = figure('visible', 'off', 'position', [100 100 1000 650]);
subplot(2,1,1); hold on;
for k = 1:numel(tc)-1
    mk = t >= tc(k) & t <= tc(k+1);
    p1 = plot((t(mk)-tc(k))/(tc(k+1)-tc(k)), CL(mk), 'r-');
end
for k = 1:numel(tc_ref)-1
    mk = t_ref >= tc_ref(k) & t_ref <= tc_ref(k+1);
    p2 = plot((t_ref(mk)-tc_ref(k))/(tc_ref(k+1)-tc_ref(k)), CL_ref(mk), 'k--', 'linewidth', 1.5);
end
grid on; xlabel('phase (t - t_k)/T'); ylabel('C_L'); box on;
legend([p1 p2], sprintf('present, %d cycles in [10,15]', numel(tc)-1), 'FeatFlow reference', 'location', 'northeast');
title('Periodicity: all shedding cycles overlaid in phase');
subplot(2,1,2); hold on;
for k = 1:numel(tc)-1
    mk = t >= tc(k) & t <= tc(k+1);
    p1 = plot((t(mk)-tc(k))/(tc(k+1)-tc(k)), CD(mk), 'b-');
end
for k = 1:numel(tc_ref)-1
    mk = t_ref >= tc_ref(k) & t_ref <= tc_ref(k+1);
    p2 = plot((t_ref(mk)-tc_ref(k))/(tc_ref(k+1)-tc_ref(k)), CD_ref(mk), 'k--', 'linewidth', 1.5);
end
grid on; xlabel('phase (t - t_k)/T'); ylabel('C_D'); box on;
legend([p1 p2], 'present', 'FeatFlow reference', 'location', 'northeast');
save_fig(h, 'fig3_cycles_overlay_vs_reference');

% spectrum + limit cycle
h = figure('visible', 'off', 'position', [100 100 1000 420]);
subplot(1,2,1);
semilogy(fL, PL/max(PL), 'r-', fD, PD/max(PD), 'b-', 'linewidth', 1.2); hold on;
grid on; xlim([0 15]); ylim([1e-6 2]);
xlabel('f [Hz]'); ylabel('normalised PSD');
legend(sprintf('C_L, f = %.3f (St = %.3f)', f_fft_L, f_fft_L*D/Umean), ...
       sprintf('C_D, f = %.3f', f_fft_D), 'location', 'northeast');
title('Spectrum on t \in [10, 15]');
subplot(1,2,2);
md = t >= t_dev(1);
plot(CL(md), CD(md), 'k-', 'linewidth', 1); grid on;
xlabel('C_L'); ylabel('C_D'); title('Limit cycle (C_L, C_D), t \in [10, 15]');
save_fig(h, 'fig4_spectrum_limit_cycle');

%% ========================================================================
%  4) Figures: flow fields
%  ========================================================================
if ~isempty(FIELD) && ~isempty(FIELD.SNAP.t)
    V   = FIELD.vertices;  E = FIELD.elements;  SN = FIELD.SNAP;
    ns  = numel(SN.t);

    % mesh
    h = figure('visible', 'off', 'position', [100 100 1100 330]);
    triplot(E', V(1,:), V(2,:), 'k'); axis equal; axis([0 2.2 0 0.41]);
    title(sprintf('P1 mesh (%d triangles); P2-P1 Taylor-Hood elements', size(E,2)));
    xlabel('x'); ylabel('y');
    save_fig(h, 'fig5_mesh');
    h = figure('visible', 'off', 'position', [100 100 700 600]);
    triplot(E', V(1,:), V(2,:), 'k'); axis equal; axis([0.05 0.45 0.02 0.38]);
    title('mesh around the cylinder'); xlabel('x'); ylabel('y');
    save_fig(h, 'fig5_mesh_zoom');

    W   = vorticity_P1(V, E, SN.ux, SN.uy);
    k15 = ns;   % last snapshot (t = 15)

    plot_field(V, E, sqrt(SN.ux(:,k15).^2 + SN.uy(:,k15).^2), jet(256), [0 2.2], ...
               sprintf('velocity magnitude |v|, t = %.2f', SN.t(k15)), 'fig6_velocity_magnitude_t15');
    plot_field(V, E, SN.ux(:,k15), jet(256), [-0.6 2.2], ...
               sprintf('streamwise velocity u_x, t = %.2f', SN.t(k15)), 'fig6_velocity_ux_t15');
    plot_field(V, E, SN.uy(:,k15), bwr(256), [-1 1], ...
               sprintf('transverse velocity u_y, t = %.2f', SN.t(k15)), 'fig6_velocity_uy_t15');
    pp = SN.p(:,k15);
    plot_field(V, E, pp, jet(256), [min(pp) max(pp)], ...
               sprintf('pressure p, t = %.2f', SN.t(k15)), 'fig6_pressure_t15');
    plot_field(V, E, W(:,k15), bwr(256), [-30 30], ...
               sprintf('vorticity \\omega_z, t = %.2f (colour range clipped to \\pm30)', SN.t(k15)), ...
               'fig6_vorticity_t15');

    % vorticity at 4 phases of one period, starting at the first lift
    % up-crossing after 14.5
    t_start = tc(find(tc >= SN.t(1), 1));
    if isempty(t_start), t_start = SN.t(1); end
    h = figure('visible', 'off', 'position', [100 100 1100 1000]);
    colormap(bwr(256));
    for q = 0:3
        tq = t_start + q*T_mean/4;
        [~, kq] = min(abs(SN.t - tq));
        subplot(4,1,q+1);
        draw_patch(V, E, W(:,kq), [-30 30]);
        title(sprintf('vorticity, t = %.2f  (phase \\approx %d/4 T,  C_L = %+.3f)', ...
              SN.t(kq), q, interp1(t, CL, SN.t(kq))));
    end
    save_fig(h, 'fig7_vorticity_4phases');

    % sampled near-wake transverse velocity (wake oscillation)
    xp = [0.5 1.0 1.5];  yp = 0.2;
    h  = figure('visible', 'off', 'position', [100 100 1000 380]); hold on;
    cols = {'r', 'b', 'k'};  lg = cell(1, numel(xp));
    for j = 1:numel(xp)
        [~, iv] = min((V(1,:)-xp(j)).^2 + (V(2,:)-yp).^2);
        plot(SN.t, SN.uy(iv,:), [cols{j} '-o'], 'markersize', 4);
        lg{j} = sprintf('u_y at (%.2f, %.2f)', V(1,iv), V(2,iv));
    end
    grid on; box on; xlabel('t [s]'); ylabel('u_y'); legend(lg, 'location', 'eastoutside');
    title('wake probes (transverse velocity), t \in [14.5, 15]');
    save_fig(h, 'fig8_wake_probes');
end
fprintf('\nPost-processing done. Figures and summary written to ./results\n');
end


%% ========================================================================
%  Local functions
%  ========================================================================
function [tc, Tk, s] = cycle_analysis(t, CD, CL, win)
% zero (mean) up-crossings of C_L -> cycles; peak statistics per cycle
m   = t >= win(1) & t <= win(2);
tt  = t(m);  cl = CL(m);  cd = CD(m);
cl0 = cl - mean(cl);
i   = find(cl0(1:end-1) < 0 & cl0(2:end) >= 0);
tc  = tt(i) - cl0(i).*(tt(i+1)-tt(i))./(cl0(i+1)-cl0(i));
Tk  = diff(tc);
nc  = numel(Tk);
CLmx = zeros(nc,1); CLmn = CLmx; CDmx = CLmx; CDmn = CLmx; CDav = CLmx;
for k = 1:nc
    mk = tt >= tc(k) & tt <= tc(k+1);
    CLmx(k) = peak3(cl(mk),  1);  CLmn(k) = -peak3(-cl(mk), 1);
    CDmx(k) = peak3(cd(mk),  1);  CDmn(k) = -peak3(-cd(mk), 1);
    CDav(k) = trapz(tt(mk), cd(mk)) / (tt(find(mk,1,'last')) - tt(find(mk,1)));
end
s.CLmax = mean(CLmx); s.CLmin = mean(CLmn);
s.CDmax = mean(CDmx); s.CDmin = mean(CDmn); s.CDmean = mean(CDav);
s.CLmax_std = std(CLmx); s.CDmax_std = std(CDmx);
end

function v = peak3(y, ~)
% maximum refined by a parabola through the 3 points around the discrete max
[v, j] = max(y);
if j > 1 && j < numel(y)
    a = y(j-1); b = y(j); c = y(j+1);
    den = a - 2*b + c;
    if den ~= 0
        v = b - (a - c)^2/(8*den);
    end
end
end

function [f, P] = spectrum(t, y, win)
m  = t >= win(1) & t <= win(2);
y  = y(m);  y = y(:) - mean(y);
n  = numel(y);  dt = median(diff(t(m)));
w  = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
N  = 2^nextpow2(16*n);
Y  = fft(y.*w, N);
P  = abs(Y(1:N/2+1)).^2;
f  = (0:N/2)'/(N*dt);
end

function W = vorticity_P1(V, E, UX, UY)
% omega = dv/dx - du/dy, P1 gradient per element, area-weighted nodal average
T  = E';
x  = V(1,:)';  y = V(2,:)';
x1 = x(T(:,1)); x2 = x(T(:,2)); x3 = x(T(:,3));
y1 = y(T(:,1)); y2 = y(T(:,2)); y3 = y(T(:,3));
A2 = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1);
bx = [y2-y3, y3-y1, y1-y2] ./ [A2 A2 A2];
cy = [x3-x2, x1-x3, x2-x1] ./ [A2 A2 A2];
wA = abs(A2);
den = accumarray(T(:), repmat(wA, 3, 1), [numel(x) 1]);
W  = zeros(numel(x), size(UX,2));
for k = 1:size(UX,2)
    ux = UX(:,k);  uy = UY(:,k);
    we = sum(bx.*uy(T), 2) - sum(cy.*ux(T), 2);
    W(:,k) = accumarray(T(:), repmat(we.*wA, 3, 1), [numel(x) 1]) ./ den;
end
end

function draw_patch(V, E, c, cl)
patch('Faces', E', 'Vertices', V', 'FaceVertexCData', c(:), ...
      'FaceColor', 'interp', 'EdgeColor', 'none');
hold on;
th = linspace(0, 2*pi, 100);
fill(0.2 + 0.05*cos(th), 0.2 + 0.05*sin(th), [0.6 0.6 0.6], 'edgecolor', 'k');
axis equal; axis([0 2.2 0 0.41]); caxis(cl); colorbar; box on;
xlabel('x'); ylabel('y');
end

function plot_field(V, E, c, cmap, cl, ttl, name)
h = figure('visible', 'off', 'position', [100 100 1100 330]);
colormap(cmap);
draw_patch(V, E, c, cl);
title(ttl);
save_fig(h, name);
end

function cmap = bwr(n)
% blue-white-red diverging colormap
x  = linspace(-1, 1, n)';
r  = min(1, 1 + x);  g = 1 - abs(x);  b = min(1, 1 - x);
cmap = [r g b];
end

function fill_window(tw, yl)
patch([tw(1) tw(2) tw(2) tw(1)], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.4], ...
      'edgecolor', 'none', 'facealpha', 0.4);
end

function save_fig(h, name)
print(h, fullfile('results', [name '.png']), '-dpng', '-r150');
close(h);
end
