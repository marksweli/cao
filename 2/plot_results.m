function plot_results()
%PLOT_RESULTS  figures for the multi-fidelity dataset (run build_datasets first)
%
%   results/fig_CD_CL_LF_vs_HF_Re100.png   LF vs HF C_D, C_L at Re = 100 on
%                                         [14.5, 15] (26 dataset samples)
%   results/fig_history_LF_vs_HF_Re100.png full history t in [0, 15]
%   results/fig_LF_HF_maps.png            lift / drag over (t, Re), LF and HF
%                                         (as Fig. 7 of the paper)
%   results/fig_LF_vs_HF_allRe.png        C_D, C_L on [14.5, 15] for all HF Re
%   results/fig_summary_vs_Re.png         mean C_D, max C_L, St versus Re
%   results/summary.csv                   same numbers in a table

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
cfgL = read_yaml(fullfile(root, 'config', 'LF.yaml'));
cfgH = read_yaml(fullfile(root, 'config', 'HF.yaml'));
DL = load(fullfile(root, [cfgL.output.dataset '.mat']));
DH = load(fullfile(root, [cfgH.output.dataset '.mat']));
out = fullfile(root, 'results');
cL = [0.85 0.33 0.10];  cH = [0 0.45 0.74];

% ---- 1) LF vs HF at Re = 100 on the sampling window ----------------------
Re0 = 100;
RL = load_case(root, cfgL, Re0);  RH = load_case(root, cfgH, Re0);
if ~isempty(RL) && ~isempty(RH)
    h = figure('visible', 'off', 'position', [100 100 1000 700]);
    names = {'CD', 'CL'};  lab = {'C_D (drag)', 'C_L (lift)'};
    for k = 1:2
        subplot(2, 1, k); hold on; box on; grid on;
        w = RL.t >= 14.5 - 1e-9;  plot(RL.t(w), RL.(names{k})(w), '-', 'color', cL, 'linewidth', 1.2);
        w = RH.t >= 14.5 - 1e-9;  plot(RH.t(w), RH.(names{k})(w), '-', 'color', cH, 'linewidth', 1.2);
        plot(RL.t_s, RL.([names{k} '_s']), 'o', 'color', cL, 'markersize', 5);
        plot(RH.t_s, RH.([names{k} '_s']), 's', 'color', cH, 'markersize', 5);
        xlim([14.5 15]); ylabel(lab{k}); xlabel('t');
        legend({sprintf('LF (coarse mesh, \\Deltat = %g)', RL.dt), sprintf('HF (fine mesh, \\Deltat = %g)', RH.dt), ...
                'LF dataset samples (26)', 'HF dataset samples (26)'}, 'location', 'eastoutside');
        title(sprintf('%s at Re = %d, t \\in [14.5, 15]', lab{k}, Re0));
    end
    save_fig(h, fullfile(out, 'fig_CD_CL_LF_vs_HF_Re100.png'));

    h = figure('visible', 'off', 'position', [100 100 1100 700]);
    for k = 1:2
        subplot(2, 1, k); hold on; box on; grid on;
        plot(RL.t, RL.(names{k}), '-', 'color', cL, 'linewidth', 1.0);
        plot(RH.t, RH.(names{k}), '-', 'color', cH, 'linewidth', 1.0);
        yl = ylim;  patch([14.5 15 15 14.5], [yl(1) yl(1) yl(2) yl(2)], [0.6 0.6 0.6], ...
                          'facealpha', 0.25, 'edgecolor', 'none');
        xlim([0 15]); ylabel(lab{k}); xlabel('t');
        legend({'LF', 'HF', 'dataset window [14.5, 15]'}, 'location', 'northwest');
        title(sprintf('%s at Re = %d, full history t \\in [0, 15]', lab{k}, Re0));
    end
    save_fig(h, fullfile(out, 'fig_history_LF_vs_HF_Re100.png'));
end

% ---- 2) maps over (t, Re), as Fig. 7 of the paper ------------------------
h = figure('visible', 'off', 'position', [100 100 1200 800]);
D = {DL, DH};  F = {'Fy', 'Fx'};  Fn = {'Lift', 'Drag'};
for r = 1:2
    lims = [min(min(DL.(F{r})(:)), min(DH.(F{r})(:))), max(max(DL.(F{r})(:)), max(DH.(F{r})(:)))];
    for c = 1:2
        subplot(2, 2, 2*(r-1)+c);
        imagesc(D{c}.t, D{c}.Re, D{c}.(F{r}));  set(gca, 'ydir', 'normal');
        caxis(lims); colormap(jet(256)); colorbar;
        xlabel('t'); ylabel('Re');
        title(sprintf('%s %s force F_%s (%d Re x %d t)', D{c}.fidelity, Fn{r}, ...
                      char('y' - (r == 2)), numel(D{c}.Re), numel(D{c}.t)));
    end
end
save_fig(h, fullfile(out, 'fig_LF_HF_maps.png'));

% ---- 3) LF vs HF for all HF Reynolds numbers -----------------------------
h = figure('visible', 'off', 'position', [100 100 1400 800]);
cm = jet(numel(DH.Re));
names = {'CD', 'CL'};
for k = 1:2
    for c = 1:2
        subplot(2, 2, 2*(k-1)+c); hold on; box on; grid on;
        Dc = D{c};
        idx = find(ismember(round(Dc.Re), round(DH.Re)));
        for j = 1:numel(idx)
            col = cm(find(round(DH.Re) == round(Dc.Re(idx(j))), 1), :);
            plot(Dc.t, Dc.(names{k})(idx(j), :), '.-', 'color', col);
        end
        xlabel('t'); ylabel(names{k}); xlim([14.5 15]);
        title(sprintf('%s: %s, Re = %s', Dc.fidelity, names{k}, mat2str(round(Dc.Re(idx)'))));
    end
end
save_fig(h, fullfile(out, 'fig_LF_vs_HF_allRe.png'));

% ---- 4) integral quantities versus Re ------------------------------------
cfgs = {cfgL, cfgH};
T = cell(1, 2);
for c = 1:2
    Res = round(linspace(cfgs{c}.reynolds.min, cfgs{c}.reynolds.max, cfgs{c}.reynolds.n)*1e6)/1e6;
    T{c} = zeros(0, 6);
    for i = 1:numel(Res)
        R = load_case(root, cfgs{c}, Res(i));
        if isempty(R), continue; end
        w = R.t >= 10 - 1e-9;  t = R.t(w);  cd = R.CD(w);  cl = R.CL(w);
        St = strouhal(t, cl, 2*cfgs{c}.physics.U_max/3, cfgs{c}.physics.D);
        T{c}(end+1, :) = [Res(i), mean(cd), max(cd), max(cl), max(cl) - min(cl), St];
    end
end
fid = fopen(fullfile(out, 'summary.csv'), 'w');
fprintf(fid, 'fidelity,Re,CD_mean,CD_max,CL_max,CL_peak_to_peak,St\n');
for c = 1:2
    for i = 1:size(T{c}, 1)
        fprintf(fid, '%s,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f\n', cfgs{c}.name, T{c}(i, :));
    end
end
fclose(fid);
h = figure('visible', 'off', 'position', [100 100 1300 400]);
ttl = {'mean C_D (t \in [10,15])', 'max C_L', 'Strouhal number St = f D / U_{mean}'};
col = [2 4 6];
for k = 1:3
    subplot(1, 3, k); hold on; box on; grid on;
    plot(T{1}(:,1), T{1}(:,col(k)), 'o-', 'color', cL);
    plot(T{2}(:,1), T{2}(:,col(k)), 's-', 'color', cH);
    xlabel('Re'); title(ttl{k}); legend({'LF', 'HF'}, 'location', 'northeast');
end
save_fig(h, fullfile(out, 'fig_summary_vs_Re.png'));
fprintf('figures written to %s\n', out);
end

% -------------------------------------------------------------------------
function R = load_case(root, cfg, Re)
f = fullfile(root, cfg.output.dir, sprintf('Re%06.2f', Re), 'result.mat');
R = [];
if exist(f, 'file'), S = load(f); R = S.R; end
end

function St = strouhal(t, cl, U, D)
% frequency from the mean period between upward zero crossings of C_L - mean
s = cl - mean(cl);
i = find(s(1:end-1) < 0 & s(2:end) >= 0);
if numel(i) < 3, St = NaN; return; end
tc = t(i) - s(i).*(t(i+1) - t(i))./(s(i+1) - s(i));
St = D/(U*mean(diff(tc)));
end

function save_fig(h, f)
print(h, f, '-dpng', '-r150');
close(h);
end
