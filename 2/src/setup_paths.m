function setup_paths()
%SETUP_PATHS  add redbKIT and the local sources to the MATLAB/Octave path
%
%   redbKIT is searched in the environment variable REDBKIT_PATH, then in
%   ../../redbKIT relative to this folder (i.e. next to the repository).
%   redbKIT must have been compiled once (>> make in MATLAB, or
%   setup_redbkit_octave.sh for GNU Octave).

here = fileparts(fileparts(mfilename('fullpath')));    % folder "2"
rk   = getenv('REDBKIT_PATH');
if isempty(rk)
    rk = fullfile(here, '..', '..', 'redbKIT');
end
if ~exist(fullfile(rk, 'FEM_library'), 'dir')
    error('redbKIT not found in "%s": set the environment variable REDBKIT_PATH', rk);
end
addpath(genpath(fullfile(rk, 'FEM_library')));
addpath(genpath(fullfile(rk, 'RB_library')));
mumps = fullfile(rk, 'FEM_library', 'LinearSolver', 'Mumps', 'Libraries');
if exist(mumps, 'dir'), rmpath(genpath(mumps)); end

addpath(fullfile(here, 'src'));
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(here, 'src', 'octave_compat'));
end
end
