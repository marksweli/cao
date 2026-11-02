function generate_meshes(cfg_files)
%GENERATE_MESHES  build the LF (coarse) and HF (fine) Gmsh meshes
%
%   GENERATE_MESHES() reads config/LF.yaml and config/HF.yaml and runs
%     gmsh <geo> -2 -format msh22 -setnumber s1 .. -setnumber s3 .. \
%          -setnumber Mesh.Algorithm .. -o <msh>
%   for each of them (Gmsh >= 4 must be on the system path). It then prints
%   the mesh statistics and writes results/mesh_info.txt and the mesh
%   figures results/fig_mesh_LF_HF.png.
%
%   GENERATE_MESHES(CFG_FILES) uses the given cell array of YAML files.

if nargin < 1
    cfg_files = {fullfile('config', 'LF.yaml'), fullfile('config', 'HF.yaml')};
end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));

[st, ver] = system('gmsh --version 2>&1');
if st ~= 0
    error('gmsh not found on the system path');
end
ver = strtrim(ver);
if ~exist(fullfile(here, 'results'), 'dir'), mkdir(fullfile(here, 'results')); end

info = {};
info{end+1} = sprintf('Gmsh version %s', ver);
M = cell(1, numel(cfg_files));
for i = 1:numel(cfg_files)
    cfg = read_yaml(fullfile(here, cfg_files{i}));
    geo = fullfile(here, cfg.mesh.geo);
    msh = fullfile(here, cfg.mesh.file);
    cmd = sprintf(['gmsh "%s" -2 -format msh22 -setnumber s1 %.10g -setnumber s3 %.10g ' ...
                   '-setnumber Mesh.Algorithm %d -o "%s" -v 2'], ...
                  geo, cfg.mesh.s1, cfg.mesh.s3, cfg.mesh.algorithm, msh);
    fprintf('%s\n', cmd);
    [st, out] = system(cmd);
    if st ~= 0
        error('gmsh failed:\n%s', out);
    end
    [V, B, E] = load_msh_2D(msh);
    hmin = min_edge(V, E);  hmax = max_edge(V, E);
    info{end+1} = sprintf(['%s: %s  (s1 = %g, s3 = %g)\n    %d nodes, %d triangles, %d boundary edges ' ...
                           '-> %d Gmsh elements;  edge length in [%.4f, %.4f]'], ...
                          cfg.name, cfg.mesh.file, cfg.mesh.s1, cfg.mesh.s3, size(V,2), size(E,2), ...
                          size(B,2), size(E,2) + size(B,2), hmin, hmax); %#ok<AGROW>
    M{i} = struct('name', cfg.name, 'V', V, 'E', E, 'nb', size(B,2));
end
info{end+1} = 'Paper (Fig. 6): LF 7798 elements / 3899 nodes, HF 16478 elements / 8239 nodes';
fid = fopen(fullfile(here, 'results', 'mesh_info.txt'), 'w');
fprintf(fid, '%s\n', info{:});
fclose(fid);
fprintf('%s\n', info{:});

% figure similar to Fig. 6 of the paper
h = figure('visible', 'off', 'position', [100 100 1100 700]);
for i = 1:numel(M)
    subplot(numel(M), 1, i);
    triplot(M{i}.E(1:3,:)', M{i}.V(1,:), M{i}.V(2,:), 'color', [0.35 0.35 0.35]);
    axis equal; axis([0 2.2 0 0.41]); xlabel('x'); ylabel('y');
    title(sprintf('%s mesh: %d triangles + %d boundary edges = %d elements, %d nodes', M{i}.name, ...
          size(M{i}.E,2), M{i}.nb, size(M{i}.E,2) + M{i}.nb, size(M{i}.V,2)));
end
print(h, fullfile(here, 'results', 'fig_mesh_LF_HF.png'), '-dpng', '-r150');
close(h);
end

function h = min_edge(V, E)
h = min(edge_lengths(V, E));
end

function h = max_edge(V, E)
h = max(edge_lengths(V, E));
end

function l = edge_lengths(V, E)
T = E(1:3, :);
l = [sqrt(sum((V(:,T(1,:)) - V(:,T(2,:))).^2, 1)), ...
     sqrt(sum((V(:,T(2,:)) - V(:,T(3,:))).^2, 1)), ...
     sqrt(sum((V(:,T(3,:)) - V(:,T(1,:))).^2, 1))];
end
