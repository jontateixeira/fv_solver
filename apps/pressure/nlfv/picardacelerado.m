function [ p, influx, bflux, flowresult ] = picardacelerado( M,I,nitpicard,tolpicard,...
                                            parameter,w,s,p,nflag,wells,mobility )

global bedge nflagface weightDMP auxface gamma formethod
      
% calculo do residuo Inicial
R0=norm(M*p-I);

%% inicializando dados para iteração Picard
step=0;
er=1;

while tolpicard <= er && step < nitpicard
    %% atualiza iterações
    step = step+1;
    
    % Acelerador de Anderson
    [ p_new ] = picfernacelerador(p,parameter,w,s,nflagface,nflag,gamma,weightDMP,auxface,wells,mobility,formethod);
    
    r = p_new(:)<0; x = find(r==1);
    if max(x)>0
        p_new(x) = 0;
    end
   
    %% Interpolação das pressões na arestas (faces)
    [pinterp_new] = pressureinterp(p,nflagface,nflag,w,s,parameter,weightDMP,mobility);
                                      
    %% Calculo da matriz global
    if strcmp(formethod,'NLFVPP')==1 
        [ M_new, I_new ] = assemblematrix_NLFV(pinterp_new,parameter,wells,mobility);
    elseif strcmp(formethod,'NLFVDMP')==1 
        [ M_new, I_new ] = assemblematrixDMPSY(p,pinterp_new,gamma,parameter,weightDMP,auxface,wells,mobility);
    end
   
    %% Calculo do residuo
    R = norm(M_new*p_new - I_new);
         
    if (R0 ~= 0.0)
        er = abs(R/R0);
    else
        er = 0.0; %exact
    end
       
    %% atualizar
    M = M_new;
    I = I_new;
    
end

p = M\I;

[pinterp] = pressureinterp(p,nflagface,nflag,w,s,parameter,weightDMP,mobility);

[flowrate,flowresult] = flowrateNLFV(p,pinterp,parameter,mobility);

bflux = flowrate(1:size(bedge,1));
influx = flowrate(1+size(bedge,1):size(flowrate,1));

end