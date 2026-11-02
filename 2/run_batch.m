function run_batch(cfg_files, Re_list)
%RUN_BATCH  batch solver for the multi-fidelity dataset (resumable)
%
%   RUN_BATCH('config/LF.yaml')  -> 19 Re = 70:5:160,  coarse mesh, dt = 0.02
%   RUN_BATCH('config/HF.yaml')  -> 10 Re = 70:10:160, fine mesh,   dt = 0.01
%   RUN_BATCH({'config/HF.yaml','config/LF.yaml'})  -> all cases, HF first
%   RUN_BATCH(CFG_FILE, RE_LIST) solves only the given Reynolds numbers.
%
%   Every case is simulated on t in [0, 15] and stored in
%   results/<LF|HF>/Re<Re>/result.mat (see src/solve_case.m).
%
%   Checkpoint / resume ("断点续算"):
%     * cases with result.mat are skipped (never recomputed);
%     * an interrupted case restarts from its last checkpoint.mat
%       (solver state incl. BDF history, saved every checkpoint.every_steps
%       steps); restarted runs are bit-identical to uninterrupted ones;
%     * so after a crash / hang / Ctrl-C / reboot simply run the same command
%       again (run_all.sh does this automatically through a watchdog).
%   Parallel workers: several MATLAB/Octave processes may run RUN_BATCH at
%   the same time; a case being solved is locked (lock.txt with the process
%   id), the other workers pick the next free case. Locks of dead processes
%   or with no progress for 20 minutes are considered stale.
%   Each case is wrapped in try/catch; results/<LF|HF>/status.csv lists the
%   state of all cases.

if ischar(cfg_files), cfg_files = {cfg_files}; end
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
setup_paths();

t_all = tic;
n_new = 0;
for f = 1:numel(cfg_files)
    cfg = read_yaml(fullfile(root, cfg_files{f}));
    if nargin < 2 || isempty(Re_list)
        Res = linspace(cfg.reynolds.min, cfg.reynolds.max, cfg.reynolds.n);
    else
        Res = Re_list;
    end
    Res = round(Res*1e6)/1e6;
    out_dir = fullfile(root, cfg.output.dir);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    [V, B, E] = load_msh_2D(fullfile(root, cfg.mesh.file));
    mesh = struct('vertices', V, 'boundaries', B, 'elements', E);
    fprintf('\n%s: mesh %s (%d nodes, %d triangles), dt = %g, T = %g, %d Reynolds numbers\n', ...
            cfg.name, cfg.mesh.file, size(V,2), size(E,2), cfg.time.dt, cfg.time.tf, numel(Res));

    for i = 1:numel(Res)
        cdir = fullfile(out_dir, sprintf('Re%06.2f', Res(i)));
        if exist(fullfile(cdir, 'result.mat'), 'file')
            fprintf('%s Re = %6.2f: done\n', cfg.name, Res(i));
            continue;
        end
        if ~claim_case(cdir)
            fprintf('%s Re = %6.2f: being solved by another worker, skipped\n', cfg.name, Res(i));
            continue;
        end
        fprintf('\n---- %s case %d / %d ----\n', cfg.name, i, numel(Res));
        try
            solve_case(cfg, Res(i), root, mesh);
            n_new = n_new + 1;
        catch err
            fprintf(2, '%s Re = %g FAILED: %s\n', cfg.name, Res(i), err.message);
            fid = fopen(fullfile(cdir, 'error.txt'), 'a');
            fprintf(fid, '%s  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), err.message);
            fclose(fid);
        end
        release_case(cdir);
        write_status(out_dir, Res, cfg.name);
    end
    write_status(out_dir, Res, cfg.name);
end
fprintf('\nrun_batch finished in %.1f min (%d cases solved in this session)\n', toc(t_all)/60, n_new);
end

% -------------------------------------------------------------------------
function ok = claim_case(cdir)
if ~exist(cdir, 'dir'), mkdir(cdir); end
lock = fullfile(cdir, 'lock.txt');
me = my_pid();
if exist(lock, 'file')
    pid = str2double(strtrim(fileread(lock)));
    if pid ~= me && ~lock_is_stale(cdir, pid)
        ok = false;
        return;
    end
end
fid = fopen(lock, 'w');  fprintf(fid, '%d\n', me);  fclose(fid);
pause(0.5);                                   % guard against simultaneous claims
ok = str2double(strtrim(fileread(lock))) == me;
end

function release_case(cdir)
lock = fullfile(cdir, 'lock.txt');
if exist(lock, 'file'), delete(lock); end
end

function stale = lock_is_stale(cdir, pid)
stale = false;
if isunix() && ~isnan(pid)
    stale = system(sprintf('kill -0 %d 2>/dev/null', pid)) ~= 0;   % process is dead
end
if ~stale
    files = {'forces.txt', 'lock.txt'};
    last = 0;
    for k = 1:numel(files)
        d = dir(fullfile(cdir, files{k}));
        if ~isempty(d), last = max(last, d(1).datenum); end
    end
    stale = (now - last)*24*60 > 20;           % no progress for 20 minutes
end
end

function pid = my_pid()
if exist('OCTAVE_VERSION', 'builtin')
    pid = getpid();
else
    pid = feature('getpid');
end
end

function write_status(out_dir, Res, name)
fid = fopen(fullfile(out_dir, 'status.csv'), 'w');
fprintf(fid, 'fidelity,Re,status,t_reached\n');
for i = 1:numel(Res)
    cdir = fullfile(out_dir, sprintf('Re%06.2f', Res(i)));
    t_last = NaN;
    if exist(fullfile(cdir, 'result.mat'), 'file')
        s = 'done';
    elseif exist(fullfile(cdir, 'lock.txt'), 'file')
        s = 'running';
    elseif exist(fullfile(cdir, 'error.txt'), 'file')
        s = 'failed';
    elseif exist(fullfile(cdir, 'checkpoint.mat'), 'file')
        s = 'interrupted';
    else
        s = 'pending';
    end
    ff = fullfile(cdir, 'forces.txt');
    if exist(ff, 'file')
        try
            A = dlmread(ff, '', 1, 0);
            t_last = A(end, 1);
        catch
        end
    end
    fprintf(fid, '%s,%.2f,%s,%.4f\n', name, Res(i), s, t_last);
end
fclose(fid);
end
