function [p,influx,bflux,flowresult]=iterpicard(M,I,...
          nitpicard,tolpicard,parameter,w,s,p,nflag,wells,mobility)

global bedge
      
% calculo do residuo Inicial
R0=norm(M*p-I);

%% inicializando dados para itera��o Picard
step=0;
er=1;

while tolpicard <= er && step < nitpicard
    %% atualiza itera��es
    step = step+1;
       
    p_new = M\I;  % invers�o sem pivotamento
    
    %% Interpola��o das press�es na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflag,w,s);
                                      
    %% Calculo da matriz global
    [ M_new , I_new ] = assemblematrix_NLFV(pinterp_new,parameter,wells,mobility);
   
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

[pinterp] = pressureinterp(p_new,nflag,w,s);

[flowrate,flowresult] = flowrateNLFV(p,pinterp,parameter,mobility);

bflux = flowrate(1:size(bedge,1));
influx = flowrate(1+size(bedge,1):size(flowrate,1));

end