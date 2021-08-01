function nflag = calflag
    % determinar o flag do nó interior e fronteira de Neumann
    global  coord bedge bcflag wells elem

    nflag=5000*ones(size(coord,1),2);
    bcflag(bcflag(:,1)<200,2)=bcflag(bcflag(:,1)<200,2)*1e-4;
    bcflag(bcflag(:,1)>200,2)=bcflag(bcflag(:,1)>200,2)/86400; % m³/d -> m³/s

    for b=1:size(bedge,1)
        f1=find(bcflag(:,1)==bedge(b,4));
        nflag(bedge(b,1),1)=bcflag(f1,1);
        nflag(bedge(b,1),2)=bcflag(f1,2);
        f2=find(bcflag(:,1)==bedge(b,5));
        nflag(bedge(b,2),1)=bcflag(f2,1);
        nflag(bedge(b,2),2)=bcflag(f2,2);
    end

    if wells(1,1)==0
        wells = double.empty(0,0);
    end
end