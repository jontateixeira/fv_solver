function [w, s, p, influx, bflux, q]=solverpressure(kmap, mobility,wells,weightDMP,...
            S_old, nflag, V, parameter, N, nw, no, Hesq, Kde,...
            nit, tol, Kn, Kt, Ded, p)

global formethod elem inedge bedge knownb
    
% calculo dos pesos que correspondem ao LPEW2
[ w, s ] = Pre_LPEW_2(kmap, mobility, V, S_old, nw, no, N);

if strcmp(formethod,'MPFAD')==1
% Se for o m�todo Linear -------------------------------------------------%
% Montagem da matriz global
[ M, I ] = globalmatrix(w, s, Kde, Ded, Kn, Kt, nflag, Hesq, wells, mobility);

% calculo das press�es
p = M\I;

% c�lculo das vaz�es 
[influx,bflux] = flowratecalc(p,w,s,Kde,Ded,Kn,Kt,Hesq,nflag,mobility);

% balan�o de fluxo
q = totalflow(influx,bflux,elem,inedge,bedge);
%-------------------------------------------------------------------------%

elseif strcmp(formethod,'NLFV')==1
% Se for o m�todo N�o-Linear ---------------------------------------------%
[pinterp]=pressureinterp(p,nflag,w,s);

% Montagem da matriz global
[ M, I ] = assemblematrix_NLFV(pinterp, parameter, wells, mobility);

% calculo das press�es
[ p, influx, bflux, q ] = iterpicard(M,I,nit,tol,parameter,w,s,p,nflag,wells,mobility);
%-------------------------------------------------------------------------%

elseif strcmp(formethod,'MPFAQL')==1
% Montagem da matriz global ----------------------------------------------%
[ M, I ] = assemblematrix_MPFAQL(parameter,w,s,nflag,weightDMP,wells,mobility);

% calculo das press�es
p = M\I;

% Press�o nodal
[pinterp]=pressureinterp(p,nflag,w,s);

% c�lculo das vaz�es 
[ influx, bflux, q ] = flowratelfvLPEW(parameter,weightDMP,mobility,pinterp,p);
%-------------------------------------------------------------------------%


elseif strcmp(formethod,'MPFAO')==1
% Se for o m�todo TPS Linear ---------------------------------------------%
knownboundlength = getknownboundlength(knownb);

[transmvecleft,transmvecright,knownvecleft,knownvecright,...
      storeinv,Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
      bodyterm] = transmTPS(kmap,knownboundlength);

[p,influx,bflux,q] = solvepressureMPFAO(transmvecleft,...
      transmvecright,knownvecleft,knownvecright,storeinv,Bleft,Bright,...
      wells,mapinv,maptransm,mapknownvec,pointedge,mobility,bodyterm);
%-------------------------------------------------------------------------%

end

end
