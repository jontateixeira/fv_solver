function varargout = pollock_jac(pos, substeps, maxsteps, reverse)
%Trace streamlines in logically Cartesian grid using Pollock approximation.
%
% SYNOPSIS:
%   [S,T,C] = pollock_jac(G, state)
%   [S,T,C] = pollock_jac(G, state, startpos)
%   [S,T,C] = pollock_jac(G, state, 'pn', pv, ...)
%   [S,T,C] = pollock_jac(G, state, startpos, 'pn', pv, ...)
%
% PARAMETERS:
%
%   G         - Cartesian or logically Cartesian grid.
%
%   state     - State structure with field 'flux'.
%
% OPTIONAL PARAMETER
%
%   positions - Matrix of size (N, 1) or (N, d+1), where d is the dimension
%               of the grid, used to indicate where the streamlines should
%               start.
%
%               If the size is (N, 1), positions contains the cell indices
%               in which streamlines should start. Each streamline is
%               started in the the local coordinate (0.5, 0.5, ...). To be
%               precise, this is the mean of the corner points, not the
%               centroid of the cell.
%
%               If the size is (N, d+1), the first column contains cell
%               indices, and the d next columns contain the local
%               coordinates at which to start streamlines.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   substeps  - Number of substeps in each cell, to improve visual quality.
%               Default 5.
%
%   maxsteps  - Maximal number of points in a streamline.
%               Default 1000.
%
%   reverse   - Reverse velocity field before tracing.
%               Default false.
%
% RETURNS:
%
%  S      - Cell array of individual streamlines suitable for calls like
%           streamline(pollock(...)) and streamtube(pollock(...)).
%
%  T      - Time-of-flight of coordinate.
%
%  C      - Cell number of streamline segment, i.e, line segment between
%           two streamline coordinates.
%
% EXAMPLE:
%
%   S = pollock_jac(G, x);
%   % pad with nan's
%   S = reshape([S, repmat({[nan, nan]}, [numel(S),1])]',[], 1);
%   S = vertcat(S{:});
%   plot(S(:,1), S(:,2), 'r-');
%
%   streamline(pollock_jac(G, x));
% SEE ALSO:
%

   global elemfluxedge
    
   positions   = pos;
   if reverse
      elemfluxedge = -elemfluxedge;
   end
     
     [varargout{1:nargout}] = trace(positions, substeps, maxsteps);

end


% ========================================================================
function varargout = trace(pos, substeps, maxsteps)
   

   global elemfluxedge coord neigh prd elem

   d              = size(coord(:,1:2), 2);
   numStreamlines = size(pos,1);
   assert(size(pos, 2) == d+1);

   magic  = 1000;
   XYZ    = nan(numStreamlines, d, magic);
   T      = nan(numStreamlines, magic);
   C      = nan(numStreamlines, magic);
   active = true(numStreamlines, 1);

   % Store crossing coordinates of active streamlines
   [XYZ(active,:,1)] = globalCoordinate(pos(active,1), pos(active, 2:end));
   T(active, 1) = zeros(sum(active), 1);
   C(active, 1) = pos(active,1);
   
   i = 2;
   while any(active)
      % Realloc
      if i+substeps+1 > size(XYZ, 3)
         magic = max(magic, substeps+1);
         XYZ   = cat(3, XYZ, nan(numStreamlines, d, magic));
         T     = cat(2, T,   nan(numStreamlines, magic));
         C     = cat(2, C,   nan(numStreamlines, magic));
      end
      current_cell = pos(active,1);
      
      
      % Take another pollock step
      [pos(active, :), t, xyz] = step(pos(active,:), elemfluxedge, neigh, substeps);
      
      % Store crossing coordinates and, optionally, coordinates along curve
      % trajectory in cell of active streamlines
      for k=1:substeps
         [XYZ(active, :, i+k-1)] = globalCoordinate(current_cell, xyz(:,:,k));
      end
      
      % Cells by NaN - "TESTE"
   
%       idx = and(~isnan(pos(:,2)),~isnan(pos(:,3)));
%       pos = pos(idx,:);  
%       active = active(idx);
%       idx = and(~isnan(xyz(:,1)),~isnan(xyz(:,2)));
%       xyz = xyz(idx,:);
       % Cells by NaN - "TESTE"
       pos(isnan(pos))=0;      

      % Mapeamento para espaço de referência
       pos(active,:) = mapeamento_piola(pos(active,1), XYZ(active, 1, i+k-1), ...
                                      XYZ(active, 2, i+k-1), elem, coord);
%      pos(abs(pos(:,2))<1e-3,2) = 0; pos(abs(pos(:,end))<1e-3,end) = 0;
%      pos(pos(:,2)>.999999,2) = 1; pos(pos(:,end)>0.999999,end) = 1;

     
      T(active, i-1+(1:substeps)) = repmat(t/substeps, [1, substeps]);
%       C(active, i-1+(1:opt.substeps)) = repmat(pos(active, 1), [1, opt.substeps]);
      C(active, i-1+(1:substeps)) = repmat(current_cell, [1, substeps]);

   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   % modified JCT (ismember function
      % Update active flag
      active(active)    =  and(pos(active,1) ~= current_cell, ...
                            ~ismember(current_cell, prd)) ;
   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      i = i+substeps;
      if i > maxsteps, break;end
   end

   %% Pack coordinates in list with streamlines separated by NaN.
   pl = reshape(permute(XYZ, [3,1,2]), [], d);

   i = ~isnan(pl(:,1));
   j = i|[true;i(1:end-1)];
   pl = pl(j,:);

   % Pack streamline coordinates in a cell array suitable for use with
   % Matlab streamline, i.e., as in 'streamline(pollock(G, resSol));'
   flag = isnan(pl(:,1));
   ix = find(flag);
   dd  = diff([0;ix])-1;
   varargout{1} = mat2cell(pl(~flag,:), dd, d);
   if nargout > 1
      T = reshape(T', [], 1);
      T = T(j);
      varargout{2} = mat2cell(T(~flag), dd, 1);
   end
   if nargout > 2
      C = reshape(C', [], 1);
      C = C(j);
      varargout{3} = mat2cell(C(~flag), dd, 1);
   end
end


% ========================================================================
function xyz = globalCoordinate(c, pl)

   global coord elem
   
% Compute global coordinate corresponding to local coorinate p in cells c
% pl  - local positions == [xi,eta]
% c  - cells

   if numel(c)==1, pl = reshape(pl, 1, []); end
   % Compute node weight for quadrilateral or hexahedron
   d = size(coord(:,1:2), 2);
   w = ones(size(pl,1), 2^d);
   for i=1:d
      mask        = logical(bitget((0:2^d-1)', i));
      w(:, mask)  = w(:, mask).* repmat(  pl(:,i), [1, sum( mask)]);
      w(:,~mask)  = w(:,~mask).* repmat(1-pl(:,i), [1, sum(~mask)]);
   end

   % Compute weighted average of corner points
   xyz = zeros(size(pl,1), d);
   for i=1:d
      xi       = coord(:,i);
      xyz(:,i) = sum( w.*reshape(xi(elem(c, [1 2 4 3]))', 2^d, [])', 2);  
   end

end


%% ========================================================================
function [pos, tof, xyz] = step(pos,elemfluxedge, neigh, nsubsteps)
% Update pos array by computing new local coordinate and new cell.
% In addition, compute curve within cell.

   global elemarea
   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   % modified JCT   
   % infinite time number ~ 1000 centuries!
   oo = 2.9e12;
   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   f = elemfluxedge(pos(:,1),:);
   n = neigh(pos(:,1),:);

   
   dims = size(pos, 2)-1;
   T    = nan(size(pos,1),dims);
   for i=1:dims
      T(:,i) = computeTime(pos(:,1+i), f(:,2*i-1:2*i));
   end
   [tof_pseudo, dir] = min(T, [], 2);
   tof = computeTimeOfFlight(pos (:,1), tof_pseudo, pos (:,2:end), f);
   
   tof = abs(tof);
   tof_pseudo = abs(tof_pseudo);
   
   % IF streamline is a source cell (prd or inj cell) pseudoTOF -> oo
   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   % modified JCT
   % bool = tof_pseudo > oo;
   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   % original version
   bool = isinf(tof_pseudo);
   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   if any(bool)
       src_cell  = pos(bool,1);
       flux_src  = f(bool,:).*repmat([-1,1,-1,1],[sum(bool),1]);
       tof_src   = elemarea(src_cell)./sum(abs(flux_src),2);
       tof(bool) = tof_src;
   end

   xyz = zeros(size(pos,1), dims, nsubsteps);
   d   = zeros(size(pos, 1), 1);
   for s=1:nsubsteps
      for i=1:dims
         t = tof_pseudo*s/nsubsteps;
         [xyz(:,i,s), d(:,i)] = computePosition(pos(:,1+i), f(:,2*i-1:2*i), t);
      end
   end

   pos (:,2:end) = xyz(:,:,s);
%    pos(bool,2:3) = 0.5;


   % Find direction to look up neighbor cell
   k  = 2*(dir-1)+d(sub2ind([numel(dir), 3], (1:numel(dir))', dir));
   t  = sub2ind(size(n), (1:numel(k))', k);  
     
   % Update cell number if NOT at boundary.
   % IF at boundary, mark dir with NaN to avoid changing local coordinate
   % below.
   ind         = n(t)==0;
   pos(~ind,1) = n(t(~ind));
%    dir (ind)   = nan;
% 
%    % Change local coordinate when moving to new cell
%    k = sub2ind(size(d), (1:size(dir,1))', dir);
%    k = k(~isnan(k));
%    pos(numel(dir) + k ) = 2-d(k);
%    
    pos = abs(pos);
  
end


% ========================================================================
function t = computeTime(xi, v)
% Compute time needed to reach xi=0 or xi=1 given velocities v=[v1,v2] at
% xi=0 and xi=1.  The formula is
%
%   t = xi/ui  or t = (1-xi)/ui,    if v1 = v2 = ui, and
%
%   t = 1/(v2-v1)*log(ue/ui),       otherwise
%
% where ui=v2*xi+v1*(1-xi) is the velocity at xi, and ue=v2 if ui>0 or
% ue=v1 if ui<0.
   tolerance = 100*eps;

   ui         = v(:,1) + xi.*diff(v, 1, 2);%(:,2)-v(:,1));
   ue         = v(:,    2);
   ue (ui<0)  = v(ui<0, 1);
   arg        = ue./ui;
   t          = inf(size(xi));

   % Linear velocity
   ind        = abs(diff(v, 1, 2)) > tolerance*abs(v(:,1));
   t(ind,:)   = 1./diff(v(ind,:), 1, 2).*log(arg(ind,:));

   % Constant velocity
   ds         = -xi;
   ds(ui > 0) = 1-xi(ui>0);
   t(~ind)    = ds(~ind)./ui(~ind);

   % nan happens for ui=ui=0
   t(arg<0 | isnan(arg))   = inf;
   
end


% ========================================================================
function t = computeTimeOfFlight(c, tof_pseudo, u0, f)
% Compute time needed to reach borders given velocities v=[v1,v2]. The 
% formula is
% 
%       tof = int_0^T (jac(eta,ksi))

 global coord elem 

ns   = size(u0, 1);     % number of streamlines
dims = size(u0, 2);     % dimension of grid

% integratoin weights and points
xi = [-.57735,.57735];
wi = [1, 1];

% Get corner points
d = size(coord(:,1:2), 2); % dimension of grid
xc = zeros(size(u0,1), size(elem(:, 1:4), 2), d);  % cell corner point
for i=1:d
    xn        = coord(:,i);
    xc(:,:,i) = xn(elem(c, 1:4)); 
end

t = zeros(ns,1);     % time of flight
x = u0;              % streamline position0
jat = tof_pseudo./2; % timestep jacobian
for i=1:numel(xi)
    time = (1+xi(i))*jat;
    x0 = x;    
    % integration of ode
    for d=1:dims
        [x(:,d), ~] = computePosition(x0(:,d), f(:,2*d-1:2*d), time);
    end
    j  = wi(i)*computeJac(x, xc);
    t = t + j;
end
t = jat.*t;
end


% ========================================================================
function [x, i] = computePosition(xi, v, t)
% Compute position at time t given start point xi and velocities v=[v1,v2].
%
%   x = xi + v*t,    if v is constant or
%
%   x = xi + (ui*exp((v2-v1)*t) - ui)/(v2-v1), otherwise
%
   tolerance = 100*eps;

   du        = diff(v, 1, 2);
   ui        = v(:,1) + xi.*du;
   i         = 1 + ~(ui<0);
   x         = inf(size(xi));

   ind       = abs(du) > tolerance*abs(v(:,1));

   % linear velocity
   x(ind)    = xi(ind) + ( ui(ind).*exp(du(ind).*t(ind)) - ui(ind))./du(ind);

   % Constant velocity
   x(~ind,:) = xi(~ind,:) + v(~ind, 1).*t(~ind, :);
   x(~ind & t==inf) = xi(~ind & t==inf);
   
   % negative coordinate get to 1e-6
   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   % modified JCT
   % x(x<0) = tolerance;
   % coordinate great than 1 get to 1
   % x(x>1+tolerance) = 1-tolerance;
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% original version
   x(x<0) = 1e-6;
   % ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
end


% =========================================================================
function j = computeJac(x, xc)
% jacobian on cell
%
% J1 = x1*y2 - y1*x2 + x2*y4 - y2*x4 + y1*x4 - x1*y4;
% J2 = y1*x2 - x1*y2 + x1*y3 - y1*x3 + y2*x4 - x2*y4 + x3*y4 - y3*x4
% J3 = y1*x3 - x1*y3 + x2*y3 - y2*x3 + x1*y4 - y1*x4 + y2*x4 - x2*y4
% 
% jac = J1 + J2 * y + J3 * x
% 
    x1 = xc(:,1,1); y1 = xc(:,1,2);
    x2 = xc(:,2,1); y2 = xc(:,2,2);
    x3 = xc(:,3,1); y3 = xc(:,3,2);
    x4 = xc(:,4,1); y4 = xc(:,4,2);
    
    J1 = x1.*y2 - y1.*x2 + x2.*y4 - y2.*x4 + y1.*x4 - x1.*y4;
    J2 = y1.*x2 - x1.*y2 + x1.*y3 - y1.*x3 + y2.*x4 - x2.*y4 + x3.*y4 - y3.*x4;
    J3 = y1.*x3 - x1.*y3 + x2.*y3 - y2.*x3 + x1.*y4 - y1.*x4 + y2.*x4 - x2.*y4;
    j = J1 + J2.*x(:,2) + J3.*x(:,1);
end


