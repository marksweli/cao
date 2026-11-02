% DATAFILE  2D flow around a cylinder for the LF / HF batch (redbKIT)
%
% Setting of Section 5 of: P. Conti, M. Guo, A. Manzoni, J.S. Hesthaven,
% "Multi-fidelity surrogate modeling using long short-term memory networks",
% CMAME 2023 (see ../../yuanzhuraoliu_LSTM.md); identical to ../../1 except
% that all parameters come from the global struct CYL_CASE set by
% solve_case.m from the YAML configuration:
%   Omega = (0,2.2) x (0,0.41) \ B_0.05(0.2,0.2),  rho = 1,  nu = U_mean D / Re
%   inlet  (flag 1): v = (4 U(t) y (H-y)/H^2, 0),
%                    U(t) = U_max (1-cos(pi t))/2 for t<1, U_max for t>=1
%   outlet (flag 2): sigma n = 0
%   walls (flag 3) and cylinder (flag 4): v = 0;  v(x,0) = 0
% This file is eval'ed by redbKIT's CFD_read_DataFile.

global CYL_CASE
c = CYL_CASE;

H      = c.H;
Umax   = c.U_max;
Umean  = 2*Umax/3;
D      = c.D;
rho    = c.rho;
nu     = Umean*D/c.Re;

%% PHYSICAL DATA
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.*x.*y);

data.bcDir_t  = @(t)( 0.5*(1-cos(pi*t)).*(t<1) + 1*(t>=1) );
data.bcDir{1} = @(x, y, t, param)( data.bcDir_t(t) * 4*Umax*y.*(H-y)./H^2 .* (x==0) + 0.*x.*y );
data.bcDir{2} = @(x, y, t, param)( 0.*x.*y );

data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

data.u0{1}    = @(x, y, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, t, param)(0.*x.*y);

data.flag_dirichlet{1} = [1 3 4];
data.flag_neumann{1}   = [2];
data.flag_dirichlet{2} = [1 3 4];
data.flag_neumann{2}   = [2];

data.dynamic_viscosity = rho*nu;
data.density           = rho;

%% Solvers
data.NonLinearSolver.tol   = c.newton_tol;
data.NonLinearSolver.maxit = c.newton_maxit;

data.LinearSolver.type     = 'backslash';
data.mumps_reordering      = 7;
data.Preconditioner.type   = 'None';
data.Preconditioner.num_subdomains = 0;   % unused; avoids gcp()

%% Time discretization (BDF)
data.time.BDF_order    = c.BDF_order;
data.time.t0           = c.t0;
data.time.dt           = c.dt;
data.time.tf           = c.tf;
data.time.nonlinearity = c.nonlinearity;

%% Output
% C_D = 2 F_x / (rho Umean^2 D),  C_L = 2 F_y / (rho Umean^2 D)
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(rho*Umean^2*D);
data.Output.DragLift.flag            = 4;
data.Output.DragLift.filename        = c.forces_file;

data.Output.Checkpoint.file          = c.ckpt_file;
data.Output.Checkpoint.every         = c.ckpt_every;
data.Output.Verbose                  = c.verbose;
if c.save_snapshots
    data.Output.Snapshots.t_start    = c.snap_t0;
    data.Output.Snapshots.every      = 1;
end
