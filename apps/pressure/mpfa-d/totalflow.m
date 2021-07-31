%função que calcula o somatório das vazões em cada elemento

function q = totalflow(influx,bflux,elem,inedge,bedge)

%-----------------somatório das vazões em cada elemento--------------------
q=zeros(size(elem,1),1);

% faces internas
for iface_int = 1:size(inedge,1)
    
    lef = inedge(iface_int,3);
    
    rel = inedge(iface_int,4);
    
    q(lef) = q(lef)+influx(iface_int);
    
    q(rel) = q(rel)-influx(iface_int);
    
end

%faces da fronteira
for iface_cont=1:size(bedge,1)
    
    lef=bedge(iface_cont,3);
    
    q(lef)=q(lef)+bflux(iface_cont);
    
end

end

















