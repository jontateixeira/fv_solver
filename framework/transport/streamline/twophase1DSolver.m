function [Res, dtof] = twophase1DSolver(tf, qelem, nw, no)
%Residual of single point upwind solver for two-phase flow along streamline.
%
% SYNOPSIS:
%   F       = twophase1DSolver(G, state, rock, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophase1DSolver returns function handles for the residual
%   matrix for the implicit upwind-mobility weighted discretization of
%
%      s_t + f(s)_tau = 0
%
%   where f is the fractional flow function,
%
%                  mw(s)
%        f(s) = -------------
%               mw(s) + mo(s)
%
%   mi = kr_i/mu_i is the phase mobility of phase i, mu_i is the phase 
%   viscosity .  The source term is neglected.
%
%   Using a first-order upstream mobility-weighted discretization in space
%   and a backward Euler discretization in time, the residual of the
%   multiples 1d equations that must be solved to move the solution
%   s from time=0 to time=tf, is obtained by calling F(s,dt),
%   defined as
%      F(s, dt) = dt/dtof穂H(s)],
%
%      H(s) = [f_i - f_j]  (cell on regular streamline)
%
%   where f_i = mw_i/(mw_i+mo_i), mo_i and mw_i are phase upwind
%   mobilities at cell i.
%
% REQUIRED PARAMETERS:
%   resSol  - Reservoir solution structure containing valid
%             saturation resSol.s with one value for each cell in
%             the grid.
%
%   G       - Grid data structure discretizing the reservoir model.
%
%   rock    - Struct with fields perm and poro.  The permeability field
%             ('perm') is only referenced when solving problems involving
%             effects of gravity and/or capillary pressure and need not be
%             specified otherwise.
%
%   trans    - two point flux transmissibilities. If given this will be
%              used not rock
%
%   fluid   - Data structure describing the fluids in the problem.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   gravity   - The current gravity in vector form. Defaults to gravity().
%
%
% RETURNS:
%   F       - Residual
%
% EXAMPLE:
%
%
%
% SEE ALSO:
%   implicitTransport, explicitTransport.

% TODO:
%   - implement gravity effects for pressure boundary and wells
%   - multipliers for gravity flux.  Handle this in 'getFlux'

% NOTE we solve
%              __
%        s_t + \/路 [D路f(s) + (G + P)路f(s)路mo(s)] = f(s) Q
%
% for constant vector and scalar fields D (Darcy flux), G (gravity flux),
% P (capillary flux), and Q (source term).
                               
   % Compute the constant tof on regular streamline
   nsl = cellfun(@(x) numel(x), tf, 'UniformOutput', false);
   dtof= cellfun(@(t,n) t(end)/(2*(n-1)), tf, nsl, 'UniformOutput', false);   
   
   % Return function handles
   Res = @(xr, dtf, dt, c) Residual(xr, dtf, dt, c, qelem, nw, no);
%    Res = @(xr, dtf, dt, cells) Residual(xr, dtf, dt, cells, fluid, q);
end


%--------------------------------------------------------------------------

function Res = Residual (resSol, dtof, dt, c, qelem, nw, no)
% F = Residual (resSol, dtof, dt, c, fluid, q)
%  DISCRETISATION
%
%  Compute  F given by a residual function Fk for each grid cell
%
%   Fk(s, s0, dt) = s - s0 + dt/pv路[H(s) - max(q,0) - min(q,0)路f],
%
%   H(s) = sum_i [f_i(dflux_i + mo_i*(gflux_i)+pcflux)]  (faces i of cell)
%
%  where f_i = mw_i/(mw_i+mo_i), mo_i and mw_i are phase upwind
%  mobilities at face i and dflux, gflux, and pcflux are Darcy flux,
%  gravity flux (=face_normal*Kg(rho1-rho2)), and capillary flux,
%  respectively. The other flux function f is the fractional flow function
%  evaluated with cell mobilities.
%
%  Advancing an implicit upwind mobility weighted scheme one time step
%  amounts to solving F(s, s0, dt) = 0, given s0 and dt.
%
%  The upwind mobility weighted flux f_i is simply computed by evaluating
%  the water mobility mw in the upwind cell wrt water flux
%
%    wflux = fw(dflux + mo*gflux).
%
%  Similarily, we evaluate mo in the upwind cell wrt oil flux
%
%    oflux = fo(dflux - mw*gflux).
%
%  The upwind cell for each phase at each internal grid face is found
%  using findPhaseUpwindCells, and computeConstData.  The gravity flux and
%  capillary flux are computed in getFlux.

   global visc satlimit
   
%    compi=[(qsat');1-(qsat')]; %satura玢o nos po鏾s injetor e produtor - CHECAR
   compi = [1 0; 0 1];
   
   s0 = 1.0;
     
   % Compute cell mobilities
   sat = [s0; resSol];     
   
   %Initialize mobility
   m = zeros(length(sat),2);
   % Compute cell mobilities
   Krw1 = ((sat - satlimit(1))/(1 - satlimit(1) - satlimit(2))).^nw;
   Kro1 = ((1 - sat - satlimit(1))/(1 - satlimit(1) - satlimit(2))).^no;
   mobw = Krw1/visc(1);  % water mobility
   mobo = Kro1/visc(2);  % oil mobility
   
   m = [mobw, mobo]; 
 
%    s0 = 1.0;
%    if any(c(1)==q.cell)
%       s0 = q.compi(c(1)==q.cell,1);
%    end

%    % Compute cell mobilities
%    sat = [s0; resSol];
%    mu  = fluid.properties(resSol);
%    kr  = fluid.relperm(sat);
%    m   = bsxfun(@rdivide, kr, mu);
%    clear mu sat kr
    
   f  = m(:,1)./sum(m,2);
   df = diff(f);

   Res   = dt/dtof*df; % (eq. 4.53) Jonathan tese
   

end


%--------------------------------------------------------------------------


