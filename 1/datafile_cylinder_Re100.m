% DATAFILE  2D flow around a cylinder, Re = 100 (redbKIT NSt_Solver)
%
% Setting of Section 5 of: P. Conti, M. Guo, A. Manzoni, J.S. Hesthaven,
% "Multi-fidelity surrogate modeling using long short-term memory networks",
% CMAME 2023 (see ../yuanzhuraoliu_LSTM.md):
%   Omega = (0,2.2) x (0,0.41) \ B_r(0.2,0.2),  r = 0.05
%   rho = 1,  sigma = -p I + 2 nu eps(v),  Re = U_mean * 2r / nu = 0.1 / nu
%   Gamma_D1 (walls y=0, y=0.41 and cylinder) : v = 0
%   Gamma_D2 (inlet x=0)                       : v = h(x,t)
%       h = ( 4 U(t) y (0.41-y) / 0.41^2 , 0 ),
%       U(t) = 0.75 (1 - cos(pi t))  for t < 1,  U(t) = 1.5 for t >= 1
%   Gamma_N  (outlet x=2.2)                    : sigma n = 0
%   v(x,0) = 0
%
% Mesh boundary flags (mesh/obstacleDFG.geo):
%   1 = inlet x=0, 2 = outlet x=2.2, 3 = top/bottom walls, 4 = cylinder
%
% Optional overrides can be passed through the global struct CYL_CFG
% (fields: Re, dt, tf, BDF_order, nonlinearity, snap_t0, snap_every,
% forces_file). This file is eval'ed by redbKIT's CFD_read_DataFile.

global CYL_CFG
cfg = struct('Re', 100, 'dt', 0.01, 'tf', 15, 'BDF_order', 2, ...
             'nonlinearity', 'implicit', 'snap_t0', 14.5, 'snap_every', 2, ...
             'forces_file', 'results/AerodynamicForces_Re100.txt');
if ~isempty(CYL_CFG)
    fn = fieldnames(CYL_CFG);
    for i = 1:numel(fn)
        cfg.(fn{i}) = CYL_CFG.(fn{i});
    end
end

H      = 0.41;            % channel height
Umax   = 1.5;             % max inflow velocity (t >= 1)
Umean  = 2*Umax/3;        % = 1
D      = 2*0.05;          % cylinder diameter (characteristic length)
rho    = 1.0;
nu     = Umean*D/cfg.Re;  % kinematic viscosity -> 1e-3 for Re = 100

%% PHYSICAL DATA
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.*x.*y);

% Inflow ramp U(t)/Umax : 0.5(1-cos(pi t)) for t<1, 1 for t>=1
data.bcDir_t  = @(t)( 0.5*(1-cos(pi*t)).*(t<1) + 1*(t>=1) );

% Dirichlet data: parabolic profile on the inlet (x == 0), zero elsewhere
data.bcDir{1} = @(x, y, t, param)( data.bcDir_t(t) * 4*Umax*y.*(H-y)./H^2 .* (x==0) + 0.*x.*y );
data.bcDir{2} = @(x, y, t, param)( 0.*x.*y );

% Homogeneous Neumann (do-nothing / zero traction) on the outlet
data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

% Initial condition v(x,0) = 0
data.u0{1}    = @(x, y, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, t, param)(0.*x.*y);

% Boundary flags
data.flag_dirichlet{1} = [1 3 4];
data.flag_neumann{1}   = [2];
data.flag_dirichlet{2} = [1 3 4];
data.flag_neumann{2}   = [2];

% Model parameters (redbKIT uses the dynamic viscosity mu = rho*nu)
data.dynamic_viscosity = rho*nu;
data.density           = rho;

% Nonlinear solver (Newton, fully implicit)
data.NonLinearSolver.tol   = 1e-8;
data.NonLinearSolver.maxit = 30;

% Linear solver
data.LinearSolver.type     = 'backslash';
data.mumps_reordering      = 7;
data.Preconditioner.type   = 'None';
data.Preconditioner.num_subdomains = 0;   % unused; avoids gcp() (Parallel Toolbox / Octave)

% Time discretization: BDF
data.time.BDF_order    = cfg.BDF_order;
data.time.t0           = 0;
data.time.dt           = cfg.dt;
data.time.tf           = cfg.tf;
data.time.nonlinearity = cfg.nonlinearity;

%% Output
% C_D = 2 F_x / (rho Umean^2 D),  C_L = 2 F_y / (rho Umean^2 D)
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(rho*Umean^2*D);
data.Output.DragLift.flag            = 4;
data.Output.DragLift.filename        = cfg.forces_file;

data.Output.Snapshots.t_start        = cfg.snap_t0;
data.Output.Snapshots.every          = cfg.snap_every;
