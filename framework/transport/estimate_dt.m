function [d_t, v] = estimate_dt(S_old,influx,bflux,CFL,nw,no,inedge,bedge, ...
            d_t, v0, tsalgo, step)
%
global  porousarea q normals

v = 0;
[dfdS1] = calcdfdS1(S_old,nw,no);

% internal faces
dfdS = 0.5*(dfdS1(inedge(:,3))+dfdS1(inedge(:,4)));
poro = 0.5*(porousarea(inedge(:,3)) + porousarea(inedge(:,4)));
d_t0 = (poro'*CFL)./(abs(dfdS.*influx)+abs(q(inedge(:,3))+q(inedge(:,4)))/2);
d_t0 = min(d_t0);

% boundary faces
dfdS = dfdS1(bedge(:,3)); dfdS(bflux < 1e-16) = 1;
poro = porousarea(bedge(:,3));
d_t0 = min([d_t0; (poro'*CFL)./(abs(dfdS.*bflux)+abs(q(bedge(:,3))))]);
if step == 1, d_t = d_t0; end

if tsalgo.dvtol > 0
   %compute the correction factor based on variation of the total velocity 
   % field over time step.
   v     = [bflux; influx]./vecnorm(normals,2,2);
   dv    =  (v - v0);
   dvel  = sqrt(sum(abs(dv).^2));
   mtime = tsalgo.dvtol/dvel;
   mtime = min( max(mtime, min(tsalgo.mtime)), max(tsalgo.mtime));
   d_t   = max([mtime*d_t, d_t0]);
   disp(['velocity variation: ',num2str(dvel)]);
   disp(mtime);
   if abs(d_t-d_t0)<1e-9
       disp(['IMPES timestep: ',num2str(d_t0)]); end
else
   % others methods:
   %      constant timestep size or (mtime x cfl) timestep   or                common cfl timestep
   d_t = d_t*(tsalgo.dvtol == -1) + d_t0*tsalgo.mtime(1)*(tsalgo.dvtol < -1) + d_t0*(tsalgo.dvtol == 0);
end
end

