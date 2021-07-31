clear; clc; format long; close all;
%-------------------------------------------------------------------------%
% Internal metric units:                                                  %
% Pressure [GPa], velocity [m/s], source/sink [m^2/s], viscosity [cP]     %                                                                %
%-------------------------------------------------------------------------%

global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge kmap ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime nflagface auxface gamma ...
    wells porousarea source O1 P1 T1 ve21 ve11 theta21 theta11 eps_ knownb ...
    neta1 esuel F q fract formethod timedimkey phflw rowposit regionfract ...
    region faces neigh state jointnodes weightDMP problemData

%-------------------------------------------------------------------------%
% input file containing all setting of the current study:
inputFile = 'workspace\fracture2_problem.inp';
%-------------------------------------------------------------------------%

[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,...
normals,esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
satlimit,pormap,bcflag,courant,totaltime,kmap,wells,porousarea,ds,tconv,...
fract,elem_post,regionfract,coord_post,formethod,timedimkey,phflw,knownb,...
rowposit,region,numsline,substeps,maxsteps,jointnodes] = preprocessor(inputFile);

%-------------------------------------------------------------------------%
% ---< Finite Volume Method Setup >---
% read problem data
problemData = INI('File', inputFile); problemData.read();
res_folder = char(problemData.problem.outputPath);

%-------------------------------------------------------------------------%
eps_ = sqrt(eps); state = []; VPI(1) = 0; elem(:,end) = region;
cumulateoil(1) = 0; oilrecovery(1) = 1; q = 0; watercut(1,1) = 0; CFL = courant;
bflux = zeros(size(bedge,1),1); influx = zeros(size(inedge,1),1);
% relative permeability data (exponential factors)
nw = problemData.fluids.nw;
no = problemData.fluids.no;

% sources terms
[ source ] = sourceterm( elem );

[n_v, n_d] = size(coord);
[n_c, n_f, n_mat] = deal(numel(elem(:,1)), numel(inedge(:,1))+numel(bedge(:,1)), ...
                         max(region));

% -----------------------------
% time step selection algorithm
v0     = zeros(size(normals,1),1);
tsalgo = problemData.timestep;

% save MATLAB Command Window output screen
diary([res_folder,filesep,mfilename,'.out'])

disp('')
fprintf('--- number of nodes: %i\n--- number of cells: %i\n', n_v, n_c);
fprintf('--- number of faces: %i\n--- number of materials: %i\n', n_f, n_mat);
if tsalgo.dvtol>0
    fprintf(['--- max timestep size: %e\n','--- vel tolerance:     %e\n',...
             '--- multipliers:  [ %.4f, %.4f]\n\n'], tsalgo.dtmax, ...
             tsalgo.dvtol, min(tsalgo.mtime), max(tsalgo.mtime));
end 

% -------------------------
% Flag's assigns
nflag = calflag;

% ========================================================================%
% GEOMETRICAL PARAMETERS
% ========================================================================%

% NLFVPP/MPFAH/NLFVDMP/MPFAQL Setups
if strcmp(formethod,'NLFVPP')==1 || strcmp(formethod,'MPFAQL')==1 ...
   || strcmp(formethod,'NLFVDMP')==1 || strcmp(formethod,'MPFAH')==1
    [csi, p, tol, nit, nflagface, auxface, gamma] = preNLFV(kmap,nflag); 
    [weightDMP] = weightnlfvDMP(kmap);
else
    csi=0;tol=0;p=0;nit=0;weightDMP=0;nflagface=0;auxface=0;gamma=0;
end

%-------------------------------------------------------------------------%
% compute interpolation parameters
[ O1, P1, T1, ve21, ve11, theta21, theta11, neta1 ] = geomparamLPEW2( coord, esurn2 );

%-------------------------------------------------------------------------%
% compute pressure related parameters
[Hesq, Kde, Kn, Kt, Ded] = Kde_Ded_Kt_Kn(kmap);

% ------------------------------------------------------------------------%
% compute saturation related parameters
[N,F,V,~,esuel,~,~,~,S_old,S_cont] = presaturation(wells,nflag);

% ------------------------------------------------------------------------%
% Preprocessor de Streamline
if ds==3
    disp('--- Streamline geometrical setups'); tic;
    [faces, neigh] = teste_unified_preprocess;
    toc; disp('--- Done.')
end


% =========================================================================
% ---< Variables initialization >---
t(1) = 0;
time = 0;                   % total simulation time process
timeref = 0;                % simulation time reference (VPI or truth time)
step = 0;                   % iteration step
VPI(1) = 0;                 % dimensionless simulation time (VPI value)
cumulateoil(1) = 0;         % cummulative oil production
oilrecovery(1) = 1;         % oil recovery
watercut(1,1)  = 0;         % water cut
countime = 0;               % a dimentional simulation time [s]

% ---- modified by JCT ----
finaltime = totaltime(2);   % simulation time
dt_ref = 0;                 % time step reference (VPI or truth time)
% -------------------------

%-------------------------------------------------------------------------%
% Initializaton of variables (initial condition)
[ S_old, step, VPI, oilrecovery, cumulateoil, watercut, dt_ref, countime, time ] = ...
  condinicial( S_old, step, VPI, oilrecovery, cumulateoil, watercut, dt_ref, countime, time );
   
vpi_old = VPI(size(VPI,2)); v = 100*vpi_old; time2 = ceil(v)/100; passo = time2*100;


%-------------------------------------------------------------------------%
% IMPES / SEQ
% initial timestep size: initial     or   constant timestep     [unity conform input param.] 
d_t = eps + finaltime/tsalgo.mtime(1)*(tsalgo.dvtol == -1);
while timeref < finaltime
    fprintf ('- Step %u ------------------------------------\n',step);
    tic;
    
    % ---------------------------------------------------------------------
    % -----< Implicit Pressure >-----
    step = step + 1;
    
    % compute mobilities
    [mobility] = mobilityface(S_old,nw,no,S_cont);
    
    % compute fractional flow
    [f_elem] = fractionalflow(S_old,nw,no);

    % Solve pressure
    [w, s, p, influx, bflux, q] = solverpressure(kmap, mobility, wells,weightDMP,...
           S_old, nflag, V, csi, N, nw, no, Hesq, Kde, nit, tol, Kn, Kt, Ded, p);

    
    % ---------------------------------------------------------------------
    % -----< Explicit Saturation >-----
    if phflw~=0
        % estimate time step size [we use the d_t*factor to convert in s]
        [d_t, v0] = estimate_dt(S_old,influx,bflux,CFL,nw,no,inedge,bedge, ...
                                d_t, v0, tsalgo, step);
        % limit time step size (for sanity)
        % we use "d_t/factor" to go back for input time unit!!
        d_t = min([tsalgo.dtmax, d_t]);
        % finaltime, timeref and d_t are in input time unit!!!! 
        d_t = min([finaltime - timeref, d_t]);
        
        if step==1
            dt_ref = d_t;
            if strcmp(formethod,'NLFVDMP')==1||strcmp(formethod,'NLFVPP')==1
                d_t = 1e-12;
            end
        end

        % Solve Saturation
        [S_old, state] = calcsaturation ( f_elem,S_cont,satlimit,visc,S_old,...
                    influx,bflux,d_t,wells,q,nw,no,elem,bedge,inedge,...
                    elemarea,pormap,ds,region );
                
        % Prodution report
        [VPI,oilrecovery,cumulateoil,watercut,countime]=reportproduction(countime,...
            VPI,wells,f_elem,step,oilrecovery,cumulateoil,watercut,q,d_t,porousarea,...
            S_old, bedge, bflux);

        vpi_old = VPI(step+1);
        countime_old = countime(step+1);
        
        if strcmp(timedimkey,'vpi')==1
            timeref = vpi_old;
        else
            timeref = countime_old;
        end        
    else
       timeref = finaltime; d_t = 0;
    end

    fprintf('p  :\t\t[%f, %f]\n', max(p)*1e4,min(p)*1e4)
    fprintf('sat:\t\t[%f, %f]\n-- PVI:\t\t%f;\n-- simtime:\t\t%g;\t-- dtime:\t%g\n', ...
        max(S_old),min(S_old),vpi_old, timeref, d_t)
    
    if max(S_old)>1.0+eps_ || min(S_old)<-eps_
       supmax=find(S_old>1.0000001);
       S_old(supmax)
       elem(supmax,:)
       submin=find(S_old<0);
       S_old(submin)
       elem(submin,:)
       pause
    end
    
    %% Write outputs
    if or(vpi_old >= time2, vpi_old > problemData.problem.maxVPI)
        
        if max(S_old)<=1.0+eps_ && min(S_old)>=-eps_
            save([res_folder, filesep, 'Saturation'],'S_old');
            save([res_folder, filesep, 'VPI'],'VPI');
            save([res_folder, filesep, 'Countime'],'countime');
            save([res_folder, filesep, 'TimeStep'],'step');
            save([res_folder, filesep, 'Watercut'],'watercut');
            save([res_folder, filesep, 'OilRecovery'],'oilrecovery');
            save([res_folder, filesep, 'CumulateOil'],'cumulateoil');
            save([res_folder, filesep, 'TIME'],'time');
            
            CopiaDeSeguranca(S_old,VPI,countime,step,watercut,oilrecovery,cumulateoil,time);
        end

        if step==1
            dt_ref=d_t; save([res_folder,filesep,'DT'],'dt_ref');
        end
        
        postprocessor(passo, vpi_old, coord_post, elem_post, p, S_old)
        if ds==3
            vtpWriter([res_folder,filesep,'PressOut'],step,state.xyz,state.tof, ...
                      state.time,'saturation',state.ssat,'CFL',state.cfl);
        end
        time2 = time2 + 0.025;
        passo = passo + 1;
    end

    time = time+toc;
    fprintf('-- runtime: %g\n', time)
    fprintf ('----------------------------------------------\n');
    if (vpi_old > problemData.problem.maxVPI)
        break; end
end
diary off
    
%% saving full production data
header = '#       sim_time             PVI       water_cut    oil_recovery     cumulateoil\n';
fid = fopen([res_folder,filesep,'report_production.out'],'w');
fprintf(fid,header);
D = [countime', VPI', watercut', oilrecovery', cumulateoil'];
% D = [VPI', watercut', oilrecovery', cumulateoil'];
save([res_folder, filesep, 'report_production.out'],'D','-ascii','-append');
fclose(fid);

subplot(2,1,1)
plot(VPI,watercut, VPI, oilrecovery)
xlabel('VPI'); ylabel('WATERCUT/OIL-RECOVERY');
subplot(2,1,2)
plot(VPI,cumulateoil)
xlabel('VPI'); ylabel('CUMULATIVE OIL');
