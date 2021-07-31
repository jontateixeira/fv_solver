function vtpWriter(basename, passo, xyz, tof, nts, varargin)
%Write vtp file format of geometry
%
% SYNOPSIS
%       vtpWriter(filename, S, T)
%       vtpWriter(filename, S, T, 'pn1', pv1, ...)
%
% PARAMETERS
%
%  filename  - Name of vtk file
%
%  S         - Cell array of individual streamlines suitable for calls like
%              streamline(pollock(...)) and streamtube(pollock(...)).
%
%  T         - Cell array of individual streamlines Time-of-flight of coordinate.
%
%  X         - Cell array of streamline segment saturation.
%
%   'pn'/pv  - List of property names/property values.  OPTIONAL.
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
%   1)>> vtpWriter('fivespot',S, T);
%   2)>> vtpWriter('fault',S, T,'X5',cheinit);
%
% SEE ALSO:
%   fopen, fclose, fprintf.
%
% NOTES
%   to see mesh you need install paraview program or similar program that
%   can be read vtk files data.
%

%%  check input

assert(ischar(basename), ...
    'First argument must be a filename.')
 
assert(iscell(xyz), ...
    'Second argument must be a cell array.')
 
assert(iscell(tof), ...
    'Third argument must be a cell array.')

assert(mod(numel(varargin), 2) == 0,...
    'Last arguments must be a list of property names/property values')

% check inputs
oset = cumsum(cellfun(@(x) size(x,1), xyz)); 
dims = max(cellfun(@(x) size(x,2), xyz));
nsl  = numel(tof);

t  = vertcat(tof{:});
idx  = (t ~= inf);
t  = t(idx);
pos = vertcat(xyz{:}); pos = pos(idx,:);
npt  = size(pos,1);
if dims<3
   pos = [pos, zeros(npt,1)];
end

ts = cell2mat( cellfun(@times, nts, ...
                  cellfun(@ones, ...
                     cellfun(@size, tof, 'UniformOutput', false), ...
                  'UniformOutput', false), ...
               'UniformOutput', false));

vararg = zeros(numel(varargin)/2,2);
for i=1:2:numel(varargin)
    assert(ischar(varargin{i}),...
        'Last arguments must be a list of property names/property values');
    [vararg(i,1), vararg(i,2)] = size(varargin{i+1});
    assert(any(vararg(i,:)==nsl), 'Invalid property size() = %u',vararg(i,1))  
end


%==================================< vtk file >==
arqvtk=strcat(basename,sprintf('%04d',passo),'.vtp');
fprintf('** creating %s file\n', arqvtk);

% vtp file headline
fvtp = fopen(arqvtk,'w');
fprintf(fvtp,'<?xml version="1.0"?>\n<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n   <PolyData>\n');
fprintf(fvtp,'      <Piece NumberOfPoints="%u" NumberOfVerts="0" NumberOfLines="%u" NumberOfStrips="0" NumberOfPolys="0">',npt,nsl);

%% streamline properties
fprintf(fvtp,'         <PointData>\n');

% time-of-flight 
fprintf(fvtp,'            <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n','TIME-OF-FLIGHT');
fprintf(fvtp,'               ');
for i=1:npt
    fprintf(fvtp,'%e ',t(i));
    if (mod(i,10)==0)
        fprintf(fvtp,'\n');
        fprintf(fvtp,'               ');
    end
end
if (mod(i,5)~=0), fprintf(fvtp,'\n'); end
fprintf(fvtp,'            </DataArray>\n');

% time step numbers 
fprintf(fvtp,'            <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n','TIME-STEPS');
fprintf(fvtp,'               ');
for i=1:npt
    fprintf(fvtp,'%e ',ts(i));
    if (mod(i,10)==0)
        fprintf(fvtp,'\n');
        fprintf(fvtp,'               ');
    end
end
if (mod(i,5)~=0), fprintf(fvtp,'\n'); end
fprintf(fvtp,'            </DataArray>\n');



for i=1:2:numel(varargin)
   fprintf(fvtp,'            <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="ascii">\n',upper(varargin{i}));
   fprintf(fvtp,'               ');
   prop = varargin{i+1};
   prop = vertcat(prop{:}); prop = prop(idx);
   prop(isnan(prop)) = eps;
   for j=1:npt
       fprintf(fvtp,'%e ',prop(j));
       if (mod(j,10)==0)
           fprintf(fvtp,'\n');
           fprintf(fvtp,'               ');
       end
   end
   if (mod(j,10)~=0), fprintf(fvtp,'\n'); end
   fprintf(fvtp,'            </DataArray>\n');
end

fprintf(fvtp,'         </PointData>\n');
fprintf(fvtp,'         <CellData>\n         </CellData>\n');

%% streamline coordinates and segments
fprintf(fvtp,'         <Points>\n');
fprintf(fvtp,'            <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">\n');
fprintf(fvtp,'               ');
for i = 1:npt
   fprintf(fvtp,'%e %e %e ',pos(i,1),pos(i,2),pos(i,3));
   if (mod(i,3)==0)
       fprintf(fvtp,'\n');
       fprintf(fvtp,'               ');
   end
end
if (mod(i,3)~=0), fprintf(fvtp,'\n'); end
fprintf(fvtp,'            </DataArray>\n         </Points>\n');

fprintf(fvtp,'         <Verts>\n            <DataArray type="Int64" Name="connectivity" format="ascii">\n');
fprintf(fvtp,'            </DataArray>\n            <DataArray type="Int64" Name="offsets" format="ascii">\n');
fprintf(fvtp,'            </DataArray>\n         </Verts>\n');

fprintf(fvtp,'         <Lines>\n            <DataArray type="Int64" Name="connectivity" format="ascii">\n');
fprintf(fvtp,'               ');
for i=1:npt
   fprintf(fvtp,'%u ',i-1);
   if (mod(i,10)==0)
       fprintf(fvtp,'\n');
       fprintf(fvtp,'               ');
   end
end
if (mod(i,10)~=0), fprintf(fvtp,'\n'); end
fprintf(fvtp,'            </DataArray>\n            <DataArray type="Int64" Name="offsets" format="ascii">\n');
fprintf(fvtp,'               ');
ii = 0;
for i=1:nsl
   fprintf(fvtp,'%u ',oset(i));
   ii = ii+1;
   if (mod(ii,9)==0)
       fprintf(fvtp,'\n');
       fprintf(fvtp,'               ');
   end
end
if (mod(ii,9)~=0), fprintf(fvtp,'\n'); end
fprintf(fvtp,'            </DataArray>\n         </Lines>\n');
fprintf(fvtp,'         <Strips>\n            <DataArray type="Int64" Name="connectivity" format="ascii">\n');
fprintf(fvtp,'            </DataArray>\n            <DataArray type="Int64" Name="offsets" format="ascii">\n');
fprintf(fvtp,'            </DataArray>\n         </Strips>\n         <Polys>\n');
fprintf(fvtp,'            <DataArray type="Int64" Name="connectivity" format="ascii">\n');
fprintf(fvtp,'            </DataArray>\n            <DataArray type="Int64" Name="offsets" format="ascii">\n');
fprintf(fvtp,'            </DataArray>\n         </Polys>\n      </Piece>\n   </PolyData>\n</VTKFile>\n');
fclose(fvtp);

end
