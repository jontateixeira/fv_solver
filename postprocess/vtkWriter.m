function vtkWriter(basename,varargin)
%Write vtkfile format of geometry
%
% SYNOPSIS
%       vtkWriter(filename, 'pn1', pv1, ...)
%
% PARAMETERS
%
%   filename - Name of vtk file
%
%   'pn'/pv  - List of property names/property values.
%              This list will be passed to block node data or cell data of
%              vtk file (i.e. if 'pn' = 'X5' and pv = cheinit [from
%              initcond.chemInitCond_] into vtk file will have the node
%              variable called X5 with values cheinit, where
%              numel(cheinit) = G.nodes.num )
%
% RETURNS
%   vtk file format with geometry data
%
% USAGE
%   1) Just write mesh
%     >> vtkWriter('fivespot','nodes',coords,'cells',cells);
%   2) write pressure field
%     >> vtkWriter('Press','nodes',coords,'cells',cells,'timestep',time,'cycle',step,'pressure',p);
%
% SEE ALSO:
%   fopen, fclose, fprintf.
%
% NOTES
%   to see mesh, you need to install VisIt/paraview program or similar 
%   program that can be read vtk files data.
%

%%  check input

assert(ischar(basename), ...
    'First argument must be a filename.')

assert(mod(numel(varargin), 2) == 0,...
    'Last arguments must be a list of property names/property values')

opt = struct('timestep',0,'cycle',0, 'verbose', false);

vararg = zeros(numel(varargin)/2,2); j = 0;
has_coords = false; has_elem = false; varin = {};
for i=1:2:numel(varargin)
    assert(ischar(varargin{i}),...
        'Last arguments must be a list of property names/property values');
    assert(isnumeric(varargin{i+1}),...
        'Last arguments must be a list of property names/property values');
    switch lower(varargin{i})
        case 'timestep'
            assert(numel(varargin{i+1})==1,...
                'timestep property must be a scalar array 1-by-1');
            opt.timestep = varargin{i+1};
        case 'cycle'
            assert(numel(varargin{i+1})==1,...
                'cycle property must be a scalar array 1-by-1');
            opt.cycle = varargin{i+1};
        case 'coord'
            [nn, ndim] = size(varargin{i+1});
            assert(nn > 1,...
                'coord property must be a scalar array n-by-d. n and d > 1');
            assert(ndim > 1,...
                'coord property must be a scalar array n-by-d. n and d > 1');
            coord = varargin{i+1}; has_coords = true;
        case 'elem'
            [nel, nc] = size(varargin{i+1});
            assert(nn > 1,...
                'elem property must be a scalar array e-by-c. e and c > 1');
            assert(ndim > 1,...
                'elem property must be a scalar array e-by-c. e and c > 1');
            elem = varargin{i+1}; has_elem = true;
        otherwise
            if j==0
                varin = {varargin{i}, varargin{i+1}};
            else
                varin = {varin{:}, varargin{i}, varargin{i+1}}; end
            j = j + 1;
            [vararg(j,1), vararg(j,2)] = size(varargin{i+1});
    end
end
vararg = vararg(1:j,:);
% check mesh info provided
assert(has_coords,'must info coordinate mesh info (coord).');
assert(has_elem,'must info connectivity mesh info (elem).');
nnel = nc - 1; ndim = ndim - double(sum(coord(:,end))==0);

%==================================< vtk file >==
arqvtk=strcat(basename,sprintf('%05d',opt.cycle),'.vtk');
fprintf('--- creating %s file\n', arqvtk);


%% ---------------------------------------------------------- vtk file headline

fvtk = fopen(arqvtk,'w');
fprintf(fvtk,'# vtk DataFile Version 2.0\n');
fprintf(fvtk,'problem created by MMGT - MATLAB Mesh Generator Toolkit\n');
fprintf(fvtk,'ASCII\n');
fprintf(fvtk,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fvtk,'FIELD FieldData 2\n');
fprintf(fvtk,'TIME 1 1 float\n');
fprintf(fvtk,'%f\n',opt.timestep);
fprintf(fvtk,'CYCLE 1 1 int\n');
fprintf(fvtk,'%d\n',opt.cycle);

%% -------------------------------------------------- node coordinates
% fprintf('writing coordinates.....');
fprintf(fvtk,'POINTS %i float\n',nn);
% coordinate
fprintf(fvtk,'%22.15e %22.15e %22.15e\r\n',coord');
% fprintf('Done!!!\n');


%% ---------------------------------------------------- Connectivity
% fprintf('writing connectivities..');
cell_type = ones(nel,1);
fprintf(fvtk,'CELLS  %i   %i \n',nel, nel*(nnel+1));

switch (nnel)
    case 3                           % 3-nodes triangle
        D = [nnel*ones(nel,1), elem(:,1:end-1)-1];
        fprintf(fvtk,'%10i %10i %10i %10i\r\n',D');
        cell_type = 5*cell_type;
    case 4                       % 4-nodes element (tetra. or qrad.)
        D = [nnel*ones(nel,1), elem(:,1:end-1)-1];
        fprintf(fvtk,'%10i %10i %10i %10i %10i\r\n',D');
        cell_type = 9*cell_type;
        if (ndim > 2)
            cell_type(:) = 10;
        end
    case 8                       % 8-nodes hexahedron
        D = [nnel*ones(nel,1), elem(:,1:end-1)-1];
        fprintf(fvtk,'%10i %10i %10i %10i %10i %10i %10i %10i %10i\r\n', D');
        cell_type(:) = 12;
end
% fprintf('Done!!!\n');


%% ---------------------------------------------------- Cells type ID
% fprintf('writing vtk_cell_type...');
fprintf(fvtk,'CELL_TYPES  %i \n',nel);
fprintf(fvtk,'%5i\r\n',cell_type);
% fprintf('Done!!!\n');

%% ----------------------------------------------------------- Point Data flux bcond
tot = 0;
for ivar = 1:size(vararg,1)
    if (vararg(ivar,1) ~= nn), continue; end
    tot = tot + 1;
end
tot = (tot)*nn;

if tot>0
    fprintf(fvtk,'POINT_DATA  %i \n',nn);
%     fprintf('writing point data......');

%% ----------------------------------------------------------- Point Data varargin
    count = 0;
    for ivar = 1:size(vararg,1)

        if (vararg(ivar,1) ~= nn), continue; end

        if vararg(ivar,2)>1
            fprintf(fvtk,'VECTORS %s float\n',varin{2*ivar-1});
            if ndim > 2
                for i=1:nn
                    fprintf(fvtk,'%f ',varin{ivar*2}(i,:));
                    fprintf(fvtk,'\n');
                end
            else
                for i=1:nn
                    fprintf(fvtk,'%f ',[varin{ivar*2}(i,:), 0]);
                    fprintf(fvtk,'\n');
                end
            end
        else
            fprintf(fvtk,'SCALARS %s float\n',varin{2*ivar-1});
            fprintf(fvtk,'LOOKUP_TABLE default\n');
            fprintf(fvtk,'%f\n',varin{ivar*2});
        end
        count = count + 1;
    end
%     fprintf('Done!!!\n');
end

%% ----------------------------------------------------------- Cells Data Materials
tot = 0;
for ivar = 1:size(vararg,1)
    if (vararg(ivar,1) ~= nel), continue; end
    tot = tot + 1;
end
tot = (tot)*nel;

if tot>0

    fprintf(fvtk,'CELL_DATA  %i \n',nel);
%     fprintf('writing cell data.......');

    %% ----------------------------------------------------------- Cells Data
    count = 0;
    for ivar = 1:size(vararg,1)

        if (vararg(ivar,1) ~= nel), continue; end

        if vararg(ivar,2)>1
            fprintf(fvtk,'VECTORS %s float\n',varin{2*ivar-1});
            if ndim > 2
                for i=1:nel
                    fprintf(fvtk,'%f ',varin{2*ivar}(i,:));
                    fprintf(fvtk,'\n');
                end
            else
                for i=1:nel
                    fprintf(fvtk,'%f ',[varin{ivar*2}(i,:), 0]);
                    fprintf(fvtk,'\n');
                end
            end
        else
            fprintf(fvtk,'SCALARS %s float\n',varin{2*ivar-1});
            fprintf(fvtk,'LOOKUP_TABLE default\n');
            fprintf(fvtk,'%f\n',varin{ivar*2});
        end
        count = count + 1;
    end
end
% fprintf('Done!!!\n');
fclose(fvtk);

end
