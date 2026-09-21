function [vertices, faces] = read_vtk_surface(vtk_file)
%READ_VTK_SURFACE Parse an ASCII VTK POLYDATA triangular surface mesh.
%   [VERTICES, FACES] = READ_VTK_SURFACE(VTK_FILE) returns VERTICES as an
%   Nx3 matrix of point coordinates and FACES as an Mx3 matrix of
%   1-based triangle vertex indices.
%
%   Original implementation (Neurodynamics_Project), independent of the
%   Pang et al. 2023 / Vohryzek et al. 2025 toolboxes.

fid = fopen(vtk_file, 'r');
if fid == -1
    error('read_vtk_surface:fileNotFound', 'Could not open %s', vtk_file);
end

vertices = [];
faces = [];

line = fgetl(fid);
while ischar(line)
    trimmed = strtrim(line);
    if startsWith(trimmed, 'POINTS')
        tok = strsplit(trimmed);
        n_points = str2double(tok{2});
        vals = fscanf(fid, '%f', n_points * 3);
        vertices = reshape(vals, 3, n_points)';
    elseif startsWith(trimmed, 'POLYGONS')
        tok = strsplit(trimmed);
        n_polys = str2double(tok{2});
        raw = fscanf(fid, '%f', n_polys * 4); % assumes triangles: [3 i j k]
        raw = reshape(raw, 4, n_polys)';
        faces = raw(:, 2:4) + 1; % 0-based -> 1-based
    end
    line = fgetl(fid);
end
fclose(fid);

if isempty(vertices)
    error('read_vtk_surface:noPoints', 'No POINTS section found in %s', vtk_file);
end
end
