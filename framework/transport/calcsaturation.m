function [ S_old, state ] = calcsaturation ( f_elem,S_cont,satlimit,visc,S_old,...
                     influx,bflux,d_t,wells,q,nw,no,elem,bedge,inedge,...
                     elemarea,pormap,ds,region )
%
state = [];
if ds==2
    
    [S_old]= firstorderstandard(S_old,influx,bflux,q,f_elem,d_t,wells,...
                                S_cont,nw,no,elem,inedge,bedge,pormap,...
                                elemarea,satlimit,visc,region);             % IMPES
      
elseif ds==1
    
    [S_old] = firstorderstdseqimp(S_old,influx,bflux,d_t,wells,q,nw,no,...
                                 elem,bedge,inedge,elemarea,pormap,region); % SEQUENCIAL IMPLÍCITO
      
elseif ds==3
    
    [S_old, state] = streamlineMPFAD_3(S_old,influx,bflux,d_t,nw,no); % STEAMLINES
    
end

            
end
