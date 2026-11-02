function build_datasets(cfg_files)
%BUILD_DATASETS  collect the solved cases into the LF and HF datasets
%
%   BUILD_DATASETS() reads results/<LF|HF>/Re*/result.mat and writes, for
%   each fidelity, the dataset given by output.dataset in the YAML file:
%     results/dataset_LF.mat / .csv   (19 Re x 26 time samples)
%     results/dataset_HF.mat / .csv   (10 Re x 26 time samples)
%
%   .mat contents (struct D):
%     Re  [nRe x 1]      Reynolds numbers
%     t   [1 x n_t]      sampling times, linspace(14.5, 15, 26)
%     CD, CL [nRe x n_t] drag / lift coefficients
%                        C = 2 F / (rho U_mean^2 D)
%     Fx, Fy [nRe x n_t] drag / lift forces on the cylinder, F = C rho U_mean^2 D / 2
%     dt, mesh, nodes, triangles, dofs, description
%   .csv: one row per (Re, t) pair, columns Re,t,CD,CL,Fx,Fy
%
%   Missing cases are reported; the datasets then contain only the
%   available Reynolds numbers.

if nargin < 1
    cfg_files = {fullfile('config', 'LF.yaml'), fullfile('config', 'HF.yaml')};
end
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

for f = 1:numel(cfg_files)
    cfg = read_yaml(fullfile(root, cfg_files{f}));
    Res = round(linspace(cfg.reynolds.min, cfg.reynolds.max, cfg.reynolds.n)*1e6)/1e6;
    ts  = linspace(cfg.sampling.t_start, cfg.sampling.t_end, cfg.sampling.n_t);
    D = struct('fidelity', cfg.name, 'Re', [], 't', ts, 'CD', [], 'CL', [], 'Fx', [], 'Fy', []);
    missing = [];
    for i = 1:numel(Res)
        rf = fullfile(root, cfg.output.dir, sprintf('Re%06.2f', Res(i)), 'result.mat');
        if ~exist(rf, 'file'), missing(end+1) = Res(i); continue; end %#ok<AGROW>
        S = load(rf);  R = S.R;
        if max(abs(R.t_s(:)' - ts)) > 1e-10
            error('sampling times of %s differ from the configuration', rf);
        end
        D.Re(end+1, 1) = R.Re;
        D.CD(end+1, :) = R.CD_s(:)';
        D.CL(end+1, :) = R.CL_s(:)';
        D.Fx(end+1, :) = R.Fx_s(:)';
        D.Fy(end+1, :) = R.Fy_s(:)';
        D.dt = R.dt;  D.mesh = R.mesh_file;
        D.nodes = R.numVertices;  D.triangles = R.numElem;  D.dofs = R.numDofs;
    end
    D.description = sprintf(['%s dataset: 2D flow past a cylinder (DFG benchmark), redbKIT P2-P1 FEM, BDF2, ' ...
        'dt = %g, mesh %s; each Re simulated on [0, %g]; forces on the cylinder sampled at %d uniform ' ...
        'times in [%g, %g]; CD/CL = 2 Fx/Fy / (rho U_mean^2 D), rho = %g, U_mean = %g, D = %g'], ...
        cfg.name, cfg.time.dt, cfg.mesh.file, cfg.time.tf, cfg.sampling.n_t, cfg.sampling.t_start, ...
        cfg.sampling.t_end, cfg.physics.rho, 2*cfg.physics.U_max/3, cfg.physics.D);
    if ~isempty(missing)
        warning('%s: %d cases missing (Re = %s)', cfg.name, numel(missing), mat2str(missing));
    end

    base = fullfile(root, cfg.output.dataset);
    save([base '.mat'], '-struct', 'D', '-v7');
    fid = fopen([base '.csv'], 'w');
    fprintf(fid, 'Re,t,CD,CL,Fx,Fy\n');
    for i = 1:numel(D.Re)
        for j = 1:numel(ts)
            fprintf(fid, '%.2f,%.4f,%.10e,%.10e,%.10e,%.10e\n', D.Re(i), ts(j), ...
                    D.CD(i,j), D.CL(i,j), D.Fx(i,j), D.Fy(i,j));
        end
    end
    fclose(fid);
    fprintf('%s: %d Re x %d t -> %s.mat / .csv\n', cfg.name, numel(D.Re), numel(ts), cfg.output.dataset);
end
end
