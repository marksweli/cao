function R = solve_case(cfg, Re, root, mesh)
%SOLVE_CASE  run (or resume) one Reynolds number of a given fidelity
%
%   R = SOLVE_CASE(CFG, RE, ROOT, MESH) solves the cylinder problem with the
%   settings of the YAML struct CFG at Reynolds number RE. ROOT is folder
%   "2" and MESH = struct('vertices', 'boundaries', 'elements').
%
%   Files in  ROOT/<cfg.output.dir>/Re<RE>/ :
%     forces.txt      time history t, C_D, C_L (written while running)
%     checkpoint.mat  solver state, rewritten every cfg.checkpoint.every_steps
%                     steps; removed when the case is complete
%     sessions.txt    one line per (re)start: start time, restart time,
%                     wall-clock time
%     result.mat      final result (struct R); its presence marks the case
%                     as DONE, so finished cases are never recomputed
%
%   R fields: Re, nu, dt, t, CD, CL, Fx, Fy (forces on the cylinder),
%   t_s/CD_s/CL_s/Fx_s/Fy_s (n_t uniform samples on [t_start, t_end]),
%   newton_its, wall_time, cfg_name, mesh info.

case_dir = case_folder(cfg, Re, root);
res_file = fullfile(case_dir, 'result.mat');
if exist(res_file, 'file')
    S = load(res_file);
    R = S.R;
    fprintf('[%s] %s Re = %6.2f: already done, skipped\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), cfg.name, Re);
    return;
end
if ~exist(case_dir, 'dir'), mkdir(case_dir); end

global CYL_CASE
c = struct();
c.Re   = Re;
c.H    = cfg.physics.H;   c.U_max = cfg.physics.U_max;
c.D    = cfg.physics.D;   c.rho   = cfg.physics.rho;
c.t0   = cfg.time.t0;     c.tf    = cfg.time.tf;    c.dt = cfg.time.dt;
c.BDF_order    = cfg.time.BDF_order;
c.nonlinearity = cfg.time.nonlinearity;
c.newton_tol   = cfg.newton.tol;
c.newton_maxit = cfg.newton.maxit;
c.forces_file  = fullfile(case_dir, 'forces.txt');
c.ckpt_file    = fullfile(case_dir, 'checkpoint.mat');
c.ckpt_every   = cfg.checkpoint.every_steps;
c.verbose      = cfg.output.verbose;
c.save_snapshots = cfg.output.save_snapshots;
c.snap_t0      = cfg.sampling.t_start;
CYL_CASE = c;

resumed = exist(c.ckpt_file, 'file') == 2;
fprintf('\n[%s] ===== %s Re = %6.2f (nu = %.4e, dt = %g): %s =====\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'), cfg.name, Re, 2*c.U_max/3*c.D/Re, c.dt, ...
        ternary(resumed, 'RESUMING from checkpoint', 'starting from t = 0'));

t_start = tic;
fem = cfg.fem;
[~, MESH, DATA, OUT] = NSt_Solver_ckpt(2, mesh.elements, mesh.vertices, mesh.boundaries, ...
                                       fem, 'datafile_cylinder');
wall = toc(t_start);

fid = fopen(fullfile(case_dir, 'sessions.txt'), 'a');
fprintf(fid, '%s  restarted_from_t = %s  wall_clock_s = %.1f\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        num2str(OUT.restarted_from), wall);
fclose(fid);

% ---- assemble result ----------------------------------------------------
Umean = 2*c.U_max/3;
q     = 0.5*c.rho*Umean^2*c.D;          % F = C * q
R = struct();
R.cfg_name   = cfg.name;
R.Re         = Re;
R.nu         = DATA.dynamic_viscosity/DATA.density;
R.dt         = c.dt;
R.t          = OUT.t;
R.CD         = OUT.Fx;
R.CL         = OUT.Fy;
R.Fx         = OUT.Fx*q;
R.Fy         = OUT.Fy*q;
R.newton_its = OUT.newton_its;
[R.t_s, R.CD_s, R.CL_s] = sample_window(R.t, R.CD, R.CL, cfg.sampling);
R.Fx_s       = R.CD_s*q;
R.Fy_s       = R.CL_s*q;
R.mesh_file  = cfg.mesh.file;
R.numVertices = MESH.numVertices;
R.numElem    = MESH.numElem;
R.numDofs    = length(MESH.internal_dof);
R.wall_time_last_session = wall;
R.sessions   = fileread(fullfile(case_dir, 'sessions.txt'));
if cfg.output.save_snapshots
    R.SNAP = OUT.SNAP;
end
if any(R.newton_its(2:end) >= c.newton_maxit)
    warning('solve_case: Newton reached maxit at %d time steps', sum(R.newton_its(2:end) >= c.newton_maxit));
end
if abs(R.t(end) - c.tf) > 1e-8
    error('solve_case: simulation stopped at t = %g instead of tf = %g', R.t(end), c.tf);
end

tmp = fullfile(case_dir, 'result_tmp.mat');
save(tmp, 'R', '-v7');
movefile(tmp, res_file, 'f');
if exist(c.ckpt_file, 'file'), delete(c.ckpt_file); end
fprintf('[%s] %s Re = %6.2f DONE in %.1f min (this session)\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        cfg.name, Re, wall/60);
end

function d = case_folder(cfg, Re, root)
d = fullfile(root, cfg.output.dir, sprintf('Re%06.2f', Re));
end

function [ts, CDs, CLs] = sample_window(t, CD, CL, smp)
ts  = linspace(smp.t_start, smp.t_end, smp.n_t)';
CDs = zeros(size(ts));  CLs = CDs;
for i = 1:numel(ts)
    [d, j] = min(abs(t - ts(i)));
    if d < 1e-8                         % sample coincides with a time step
        CDs(i) = CD(j);  CLs(i) = CL(j);
    else
        CDs(i) = interp1(t, CD, ts(i), 'spline');
        CLs(i) = interp1(t, CL, ts(i), 'spline');
    end
end
end

function r = ternary(c, a, b)
if c, r = a; else, r = b; end
end
