%MAIN_CYLINDER_RE100  2D flow past a cylinder at Re = 100 with redbKIT
%
%   Single-case baseline (HF setting) of Section 5 of
%   P. Conti, M. Guo, A. Manzoni, J.S. Hesthaven, "Multi-fidelity surrogate
%   modeling using long short-term memory networks", CMAME 2023.
%
%   - unsteady incompressible Navier-Stokes, rho = 1, nu = 0.1/Re = 1e-3
%   - Omega = (0,2.2)x(0,0.41) \ B_0.05(0.2,0.2)
%   - parabolic inflow with U(t) = 0.75(1-cos(pi t)) (t<1), 1.5 (t>=1)
%   - no-slip on walls/cylinder, zero traction on the outlet, v(x,0) = 0
%   - P2-P1 Taylor-Hood FEM, BDF2, dt = 0.01 (HF time step), T = 15
%
%   Usage (MATLAB or GNU Octave >= 6):
%     1) download redbKIT v2.2 (https://github.com/redbKIT/redbKIT) and
%        compile it once ( >> make  in the redbKIT root folder );
%     2) set REDBKIT_PATH below (or the environment variable REDBKIT_PATH);
%     3) cd to this folder and run  >> main_cylinder_Re100
%   Outputs are written to ./results ; post-processing is done by
%   postprocess_cylinder_Re100.m (can also be run stand-alone afterwards).

clear; clc; close all;

%% ------------------------------------------------------------------------
%  User settings
%  ------------------------------------------------------------------------
REDBKIT_PATH = getenv('REDBKIT_PATH');
if isempty(REDBKIT_PATH)
    REDBKIT_PATH = fullfile(pwd, '..', '..', 'redbKIT');   % <-- edit if needed
end

global CYL_CFG
CYL_CFG = struct();
CYL_CFG.Re           = 100;          % Reynolds number  (nu = 0.1/Re)
CYL_CFG.dt           = 0.01;         % HF time step of the paper (LF: 0.02)
CYL_CFG.tf           = 15;           % final time of the paper
CYL_CFG.BDF_order    = 2;
CYL_CFG.nonlinearity = 'implicit';   % Newton at each time step
CYL_CFG.snap_t0      = 14.5;         % store flow-field snapshots for t >= 14.5
CYL_CFG.snap_every   = 2;            % every 2 steps -> 26 snapshots in [14.5,15]
CYL_CFG.forces_file  = fullfile('results', 'AerodynamicForces_Re100.txt');

mesh_file = fullfile('mesh', 'obstacleDFG');   % gmsh P1 mesh (flags 1..4)
fem       = {'P2', 'P1'};                      % Taylor-Hood elements
dim       = 2;

%% ------------------------------------------------------------------------
%  Paths
%  ------------------------------------------------------------------------
if ~exist(REDBKIT_PATH, 'dir')
    error('redbKIT not found in "%s". Set REDBKIT_PATH.', REDBKIT_PATH);
end
addpath(genpath(fullfile(REDBKIT_PATH, 'FEM_library')));
addpath(genpath(fullfile(REDBKIT_PATH, 'RB_library')));
mumps_dir = fullfile(REDBKIT_PATH, 'FEM_library', 'LinearSolver', 'Mumps', 'Libraries');
if exist(mumps_dir, 'dir'), rmpath(genpath(mumps_dir)); end
addpath(pwd);   % make sure the local NSt_Solver_snap/datafile are used
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(pwd, 'octave_compat'));   % small shims, GNU Octave only
end

if ~exist('results', 'dir'), mkdir('results'); end
diary(fullfile('results', 'run_log.txt'));

%% ------------------------------------------------------------------------
%  Mesh + solve
%  ------------------------------------------------------------------------
% equivalent to redbKIT's msh_to_Mmesh(mesh_file, dim), but MATLAB/Octave
% portable (redbKIT's read_msh is distributed only as MATLAB P-code)
[vertices, boundaries, elements] = load_msh_2D(mesh_file);

t_cpu = tic;
[U, MESH, DATA, SNAP] = NSt_Solver_snap(dim, elements, vertices, boundaries, ...
                                        fem, 'datafile_cylinder_Re100');
t_cpu = toc(t_cpu);
fprintf('\n\nTotal wall-clock time: %.1f s (%.2f h)\n', t_cpu, t_cpu/3600);

%% ------------------------------------------------------------------------
%  Save flow field (P1 vertices) and run info
%  ------------------------------------------------------------------------
FIELD.vertices  = MESH.vertices(1:2, :);
FIELD.elements  = MESH.elements(1:3, :);
FIELD.SNAP      = SNAP;
FIELD.cfg       = CYL_CFG;
FIELD.nu        = DATA.dynamic_viscosity / DATA.density;
FIELD.numDofs   = length(MESH.internal_dof);
FIELD.numNodes  = MESH.numNodes;
FIELD.numElem   = MESH.numElem;
FIELD.t_cpu     = t_cpu;
save(fullfile('results', 'flowfield_Re100.mat'), 'FIELD', '-v7');

diary off;

%% ------------------------------------------------------------------------
%  Post-processing (figures, C_D / C_L statistics, Strouhal number)
%  ------------------------------------------------------------------------
postprocess_cylinder_Re100();
