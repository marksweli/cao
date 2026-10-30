function [vertices, boundaries, elements] = load_msh_2D(filename)
%LOAD_MSH_2D  read a 2D gmsh ASCII v2 mesh in redbKIT format
%
%   [VERTICES, BOUNDARIES, ELEMENTS] = LOAD_MSH_2D(FILENAME) returns the
%   same arrays as redbKIT's msh_to_Mmesh(FILENAME, 2), but is written in
%   plain MATLAB/Octave code (redbKIT ships read_msh only as MATLAB P-code,
%   which GNU Octave cannot execute).
%     vertices   : 2 x Nv  node coordinates
%     elements   : 4 x Ne  [node1; node2; node3; physical tag]
%     boundaries : 7 x Nb  [node1; node2; 0; 0; physical tag; 1; 0]

if numel(filename) < 4 || ~strcmp(filename(end-3:end), '.msh')
    filename = [filename, '.msh'];
end
fprintf('\nLoading mesh from file %s ... ', filename);
t0  = tic;
fid = fopen(filename, 'r');
if fid < 0, error('Cannot open %s', filename); end
cleaner = onCleanup(@() fclose(fid));

line = '';
while ~strcmp(strtrim(line), '$Nodes')
    line = fgetl(fid);
    if ~ischar(line), error('$Nodes section not found'); end
end
nN    = str2double(fgetl(fid));
nodes = fscanf(fid, '%f', [4, nN]);
vertices = zeros(2, nN);
vertices(:, nodes(1,:)) = nodes(2:3, :);

while ~strcmp(strtrim(line), '$Elements')
    line = fgetl(fid);
end
nE   = str2double(fgetl(fid));
tri  = zeros(4, nE);  nt = 0;
seg  = zeros(3, nE);  ns = 0;
for k = 1:nE
    v     = sscanf(fgetl(fid), '%d')';
    etype = v(2);  ntags = v(3);
    phys  = v(4);
    conn  = v(4+ntags:end);
    switch etype
        case 1                         % 2-node line
            ns = ns + 1;  seg(:, ns) = [conn(1:2)'; phys];
        case 2                         % 3-node triangle
            nt = nt + 1;  tri(:, nt) = [conn(1:3)'; phys];
    end
end
elements = tri(:, 1:nt);
seg      = seg(:, 1:ns);

boundaries        = zeros(7, ns);
boundaries(1:2,:) = seg(1:2,:);
boundaries(5,:)   = seg(3,:);
boundaries(6,:)   = 1;

fprintf('done in %3.2f s (%d vertices, %d triangles, %d boundary edges)\n\n', ...
        toc(t0), nN, nt, ns);
end
