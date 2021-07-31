function [S_old, state] = streamlineMPFAD(p, influx, bflux, S_old, numsline, substeps, maxsteps, tf)
%Explicit 1D transport solver for two-phase flow along streamlines.
%
% DESCRIPTION:
%   Explicit discretization of the transport equation
%
%      s_t + f_tau = 0
%
%   where f is the fractional flow function,
%
%                  mw(s)
%        f(s) = -------------
%               mw(s) + mo(s)
%
%   mi = kr_i/mu_i is the phase mobility of phase i, mu_i is the phase 
%   viscosity.
%
%   We use a first-order upstream mobility-weighted discretization in space
%   and a backward Euler discretization in time. The transport equation is
%   solved on the time interval [0,tf] by calling twophaseJacobian to build
%   a function computing the residual of the discrete system in addition to
%   a function taking care of the update of the solution during the
%   time loop.
%
%
% RETURNS:
%   S_old     - Reservoir solution with updated saturation.
%

    global wells q elem elemfluxedge faces inedge bedge coord courant prd

    opt = struct(...
        'numsline' , numsline   , ...
        'max_dt'   , inf        , ...
        'dt_factor', 0.8        , ...
        'satwarn'  , sqrt(eps)  , ...
        'substeps' , substeps   , ...
        'maxsteps' , maxsteps   );
%         'dt'       , tf);

	influx_cp = influx;
    bflux_cp  = bflux;
      
    % Fluxos por face em cada elemento de maneira ordenada (w,e,s,n)
    
    % Para as faces de contorno 
    elemfluxedge = zeros(size(faces)); 
    edges = [influx; -bflux];
    elemfluxedge = edges(faces);
    
    % Para as faces internas
    for i=1:size(faces,1)
        for j=1:4
            if faces(i,j)<=size(inedge,1)
                if j==1 || j==3
                    if inedge(faces(i,j),3)==i %&& influx(faces(i,j))>0
                        elemfluxedge(i,j)=-influx(faces(i,j));
                    end
%                     if inedge(faces(i,j),4)==i && influx(faces(i,j))<0
%                         elemfluxedge(i,j)=-influx(faces(i,j));
%                     end    
                elseif j==2 || j==4
                    if inedge(faces(i,j),4)==i %&& influx(faces(i,j))<0
                        elemfluxedge(i,j)=-influx(faces(i,j));
                    end
%                     if inedge(faces(i,j),3)==i &&  influx(faces(i,j))<0
%                         elemfluxedge(i,j)=-influx(faces(i,j));
%                     end
                end
            end
        end           
    end
    
   % Compute wells    
   qelem = wells(:,1); qflux = q(qelem); qsat = wells(:,4); prd = qelem(qflux<0)';
   S_old(qelem(1)) =  1;  % saturacao do injetor (inicialmente é 1)
   
   %% Regular tracing
   % streamline launching
   pos  = launching(opt, qelem, qflux, qsat);
   dims = 2;
   vst  = zeros(size(elem,1),1);
   cxyz = []; %#ok<*NASGU>
   ctof = [];
   ccel = [];
   csat = [];
   cnts = [];
   c    = [];
   dtf  = [];
   u    = [];
   is_reverse_trace = false;
  % while any(vst==0)
      % streamline tracing
      [S,T,C] = tracing(pos, substeps, maxsteps, qelem, qflux, any(vst>0));
      
      if is_reverse_trace
          n = numel(C); slseminj = false(n,1);
          % remover linha que nao alcancaram o injetor
          for i = 1:n
              slseminj(i) = qelem(1)~=C{i}(1);
              vst(C{i})= vst(C{i})+ slseminj(i);
          end
          % deixar apenas as sl q alcancaram o injetor
          C = C(~slseminj);
          S = S(~slseminj);
          T = T(~slseminj);
      end
      
      % mapping
      U = cellfun(@(t, c) map_regular(S_old, t, c), T, C, 'UniformOutput', false);
      
      % ------- EQUAÇÃO DA SATURAÇÃO --------%
      [Res, dtof] = twophase1DSolver(T, qelem, qsat);
%     [F, dtof] = twophase1DSolver(T, fluid, q);

      % ---------- Time step estimate from state ---------------
      getdt = @(s, t) estimate_dt(s, t, courant, opt.max_dt);
%     getdt = @(s, t) estimate_dt(s, t, fluid, opt.dt_factor, opt.max_dt);
       
     
      s  = U;
      t  = zeros(numel(s),1);
      n  = zeros(numel(s),1);
      
      while any(t < tf)
         % min(tf-t,getdt(s))
         dt      = min([tf-t, cellfun(getdt, s, dtof)], [], 2);
         % s = s - (t < tf)*F(s,dt,c)
         s       = cellfun(@minus, s, ...
                     cellfun(@times, mat2cell(t < tf,ones(numel(dt),1)), ...
                        cellfun(Res, s, dtof, mat2cell(dt,ones(numel(dt),1)), C, ...
                        'UniformOutput', false),...
                      'UniformOutput', false),...
                   'UniformOutput', false);
         t       = t + dt;
         n(t<tf) = n(t<tf) + 1;
         % Correct possible bad saturations
         s       = cellfun(@(x) correct_saturations(x, opt.satwarn), s, ...
                            'UniformOutput', false);
      end

      
      % mapping to simulation grid
      [sm, tm, cm, i, X] = mapping(s, T, dtof, C);
    %  vst = vst + i; 
%      % if any(vst==0)
%          ci = (1:size(elem,1))'; ci = ci(vst==0);
%          pos = [ci, ci, 0.5*ones(numel(ci),dims)];
%          disp('** missed cells found:')
%          fprintf('%5d\n', ci);
%          is_reverse_trace = true;
%      % end
      c    = [  c; cm];
      u    = [  u; sm];
      dtf  = [dtf; tm];
      
      cxyz = [cxyz; S];
      ctof = [ctof; T];
      ccel = [ccel; C];
      csat = [csat; X];
      cnts = [cnts; mat2cell(n+1,ones(numel(s),1))];
      
  %  end
   
   %% map to simulation grid
   % Compute weithed average of time of flight in each cell that has been
   % crossed by one or more streamlines,
   %
   %    tof_i = sum(t_ij*ct_ij; j) / sum(t_ij; j)
   %
   % where t_ij is the time-of-flight diff through cell i of streamline j

   sat = accumarray(c, u.*dtf, [size(elem,1), 1], [], inf)./...
      accumarray(c, dtf,     [size(elem,1), 1], [], 1);
  
   % celulas onde nao havera slines (dominio ignorado) saturacao nao se
   % alterara
   if  max(elem(:,end))==2001
       sat(elem(:,end)==2001) = S_old(elem(:,end)==2001);
   elseif any(sat==inf)
       sat(sat==inf) = S_old(sat==inf);
   end
   
   % Cells by NaN - "TESTE"
   sat(isnan(sat))=0;
   
   sat = [sat, 1 - sat];
   S_old = sat(:,1);      % saturação de água em cada célula 
   S_old(qelem(1)) =  1;  % saturacao do injetor (inicialmente é 1)
   
   %% last operations
%    if opt.onlygrav,
%       state.flux = flux;
%    end

%    if any(any(isnan(S_old))),
%       error('streamlineTransport: Transport step failed')
%    end

   
   state.xyz  = cxyz;
   state.tof  = ctof;
   state.cell = ccel;
   state.ssat = csat;
   state.time = cnts;
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function s = correct_saturations(s, satwarn)
   % Display error if s > 1+satwarn
   % or s < 0 - satwarn
   i = find(s > 1 + satwarn);
   if ~isempty(i)
      disp('Saturation exceeds 1 in cells:')
      fprintf('%5d %g\n', [i, s(i)] .');
      warning('explicitTransportStreamline: Error larger than satwarn, reduce timestep size')
   end

   i = find(s(:,1) < -satwarn);
   if ~isempty(i)
      disp('Saturation less than 0 in cells:')
      fprintf('%5d %g\n', [i, s(i)] .');
      warning('explicitTransportStreamline: Error larger than satwarn, reduce timestep size')
   end
   % Correct numerical errors
   s(s(:) > 1, 1) = 1 - eps;
   s(s(:) < 0, 1) = 0;
end


%--------------------------------------------------------------------------

function pos = launching(opt, qelem, qflux, qsat)

   global neigh elemfluxedge elem coord

   % getting injection terms
   inj  = qflux>0;
   cell_src = qelem(inj);
      
   % grid dimension 
   griddim = 2;

   % splitting streamline per injetor cells
   cells0 = abs(elemfluxedge(cell_src,:)) > sqrt(eps);  % direction to launch (w,e,s,n)
   nsline = repmat(round(opt.numsline*qflux(inj)./sum(qflux(inj))...
   ./sum(cells0,2)),[1, size(neigh,2)]);
   nsline = nsline.*cells0; % number of streamline per face to launch
   icells = repmat(cell_src,[1,size(nsline,2)]);
   dir    = repmat([1,1,2,2],[size(nsline,1),1]);
   xzero  = repmat([0, 1],[numel(cell_src), griddim]);
   neighbor = neigh(cell_src,:); 

   icells = icells(nsline>0);
   dir    = dir(nsline>0);
   neighbor = neighbor(nsline>0);
   xzero  = xzero(nsline>0);
   nsline = nsline(nsline>0);

   pos = zeros(sum(nsline(:)), griddim + 2); j = 0;
   for i=1:numel(nsline)
      pos(j+1:j+nsline(i),1) = icells(i);
      pos(j+1:j+nsline(i),2) = neighbor(i);
      pos(j+1:j+nsline(i),2+dir(i)) = xzero(i);
      pos(j+1:j+nsline(i),2+3-dir(i)) = linspace(0.001, 0.999, nsline(i));
      j = j + nsline(i);
   end
   % compute global coord (outside points - injector cells)
   xyz = globalCoordinate(pos(:,1), pos(:,3:end));
   % compute local coords for next cells to trace
   pos(:,2:end) = mapeamento_piola(pos(:,2), xyz(:,1), xyz(:,2), elem, coord);
      
end


%--------------------------------------------------------------------------

function [S,T,C] = tracing(pos, substeps, maxsteps, qelem, qflux, reverse)

    global elemarea pormap centelem elem

   % trace streamlines
   [S,T,C] = pollock_jac(pos(:,2:end), substeps, maxsteps, reverse);
%    [S,T,C] = pollock_numeric(G, s,  pos(:,2:end), 'substeps', opt.substeps, ...
%                             'maxsteps', opt.maxsteps, 'reverse', reverse);

   pormap = pormap(1)*ones(size(centelem,1),1);
   porosity = pormap(elem(:,5)); % porosity
   % getting injection terms
   injector = qflux>0;
   cell_src = qelem(injector);
   tof_src  = elemarea(cell_src).*porosity(cell_src)./sum(abs(qflux(injector)),2);

   % source cells
   cells0 = pos(:,1); 
   tof0   = zeros(size(cells0));
   
   for i=1:numel(cell_src)
      tof0(cells0 == cell_src(i)) = tof_src(i);
   end

   % getting production terms
   producer  = qflux<0;
   cell_sink = qelem(producer);
   tof_sink  = elemarea(cell_sink).*porosity(cell_sink)./sum(abs(qflux(producer)),2);

   % completing streamlines on source and sink positions
   poro = porosity;
   if (~reverse)
      for i=1:numel(S)
         x    = S{i};
         t    = T{i};
         c    = C{i};
         t    = poro(c).*t;
         if isinf(t(end))
            tof  = tof_sink(c(end) == cell_sink);
            S{i} = [centelem(cells0(i),1:2); 
                    x(1:end-substeps,:);
                    centelem(c(end),1:2)];
            T{i} = cumsum([0.0; tof0(i); t(2:end-substeps); tof]);
            C{i} = [cells0(i); cells0(i); c(2:end-substeps+1)];
         else
            S{i} = [centelem(cells0(i),1:2); x];
            T{i} = cumsum([0.0; tof0(i); t(2:end)]);
            C{i} = [cells0(i); cells0(i); c(2:end)];
         end
      end
   else
      for i=1:numel(S)
         x    = S{i};
         t    = T{i};
         c    = C{i};
         t    = poro(c).*t;
         if isinf(t(end))
            S{i} = [centelem(cells0(i),1:2);
                    x(end-substeps:-1:1,:)];
            T{i} = cumsum([0.0; tof0(i); t(end-substeps:-1:2)]);
            C{i} = [cells0(i); c(end-substeps+1:-1:1)];
         else
            S{i} = x(end:-1:1,:);
            T{i} = cumsum([0.0; t(end:-1:2)]);
            C{i} = c(end:-1:1);
         end
      end
   end
end


%--------------------------------------------------------------------------

function var = map_regular(sat,t, c)
   
   n  = numel(t);
   x  = sat(c);

   dt = t(end)/(2*(n-1));
   xr = zeros(2*(n-1),1);

   for j=1:2*(n-1)
      lb = (j-1)*dt;
      ub =     j*dt;
      val=0.0;
      for k=2:n
         lb2 = t(k-1); if (lb2 > ub), break; end
         ub2 = t(k);   if (ub2 < lb), continue; end
         low=lb2;      if (lb >= lb2), low=lb; end
         upper=ub2;    if (ub <= ub2), upper=ub; end
         val = val + x(k)*(upper-low);
      end
      xr(j) = val/dt;
   end
   var = xr;
end


%--------------------------------------------------------------------------

function dt = estimate_dt(s, tof, courant, max_dt)
       
   global visc no nw

%    mu = visc; 
%    %% Compute cell mobility and its derivative
%    sat       = s;

   %% Find simple approximation to the maximal wave speed from advective
   % term in reservoir based on derivative of flux on face.
%    df = @(mob, dmob) (mob(:,2).*dmob(:,1) - mob(:,1).*dmob(:,2))./(sum(mob, 2).^2);
%    d  = df(mob, dmob);
   d = calcdfwdS(s,nw,no); 
   m  = max(abs(d), [], 2);


   wavespeed  = max(m);
   
   dt = min([max_dt, courant*tof/wavespeed]);
%     dt = min([max_dt, ncourant*tof/wavespeed]);

end


%--------------------------------------------------------------------------

function [x, dt, c, i, ss] = mapping(s, tf, dtf, c)

   global elem 

   ss = cellfun(@map_streamline, s, tf, dtf, 'UniformOutput', false);
   
   % time of flight diff along streamlines
   dt = cellfun(@diff, tf, 'uniformoutput', false);
   dt = cellfun(@(x) [0; x], dt, 'uniformoutput', false);

   % rearrange to plain double arrays
   dt  = vertcat(dt{:});
   i= (dt ~= inf);
   dt  = dt(i);
   tf = vertcat(tf{:}); tf = tf(i);
   x = vertcat(ss{:});  x  = x(i);
   c  = vertcat(c{:});  c  = c(i);
   
   i = accumarray(c,1,[size(elem,1),1]);
   i(elem(:,end)==2)=1; % celulas q estao fora do calculo da saturacao (nao serao tracadas slines!)

end


%--------------------------------------------------------------------------

function var = map_streamline(s, t, dt)
   n  = numel(t);
   x  = s;
   xs = zeros(n,1);
   xs(1) = 1; % pq definimos o valor da saturacao no poco inj =1 desde sempre! 
              % ver linha 76 
   for i=2:n
      lb = t(i-1);
      ub = t(i);
      val=0.0;
      for j=1:2*(n-1)
         lb2  = (j-1)*dt; if (lb2 > ub), break; end
         ub2  =    j *dt; if (ub2 < lb), continue; end
         low  =      lb2; if (lb >= lb2), low = lb; end
         upper=      ub2; if (ub <= ub2), upper=ub; end
         val = val + x(j)*(upper-low);
      end
      xs(i) = val/(ub-lb);
      if abs(ub-lb)<sqrt(eps)
         xs(i) = xs(i-1);
      end
   end
   var = xs;
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
