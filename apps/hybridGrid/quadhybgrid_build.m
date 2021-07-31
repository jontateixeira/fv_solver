function [ coord, elem, fract, kmap, region, regionfract, flagcorr ] = quadhybgrid_build( coord, elem, ...
                       centelem, fract, kfrat, kmap, spclnodes, region, regionfract, flaglim )
%
elemn=zeros(size(elem,1),size(elem,2));
nc=size(coord,1);
g=1;

t=1; d=1; bnf=0; ebnf=0;
fratnodes=zeros(1,1);
for i=1:size(fract,1)
    for j=1:2
        outrono=fract(i,j);
        if isempty(find(fratnodes==outrono))~=0
            fratnodes(t,1)=outrono;
            t=t+1;
        end
    end
    for j=1:4
        if fract(i,j)~=0 && fract(i,j)<=size(flaglim,1)
            bnf(d)=fract(i,j);
            ebnf(d)=i;
            d=d+1;
        end
    end
    nnf(i,1)=sum(fract(:,1)==fract(i,1))+sum(fract(:,2)==fract(i,2))...
             +sum(fract(:,1)==fract(i,2))+sum(fract(:,2)==fract(i,1));
    cf(i,:)=0.5*(coord(fract(i,1),:)+coord(fract(i,2),:));
end

fractn=fract; confcent=0; v=0; vc=0; tc=1; ipo2=1;

% Percorre os nós que pertencem às fraturas ------------------------------%
for i=1:size(fratnodes,1)
    
    no = fratnodes(i);
    
    clear e v vd vc ce n f c d ff l fm
    
    % Verifica quais são os elementos que compartilham o nó em questão ---%
    esurn = [find(elem(:,1)==no); find(elem(:,2)==no); find(elem(:,3)==no); find(elem(:,4)==no)]';
    
    % Verifica os centróides desses elementos ----------------------------%
    centeresurn = centelem(esurn,:);
       
    % Verifica a que fraturas pertencem esses nós ------------------------%    
    fractsurn = [find(fract(:,1)==no); find(fract(:,2)==no)]'; % Fratura
    outrono = [fract(find(fract(:,1)==no),2); fract(find(fract(:,2)==no),1)]'; % Outro nó da Fratura
    fractreg = regionfract(fractsurn); fractreg = fractreg - 2000;
    
    % Vetor Fratura ------------------------------------------------------%
    for j = 1:size(outrono,2)
        v(j,:) = coord(outrono(j),:) - coord(no,:);
    end
    % Vetor nó - centróide dos elementos ---------------------------------%
    for j = 1:size(esurn,2)
        vc(j,:) = centeresurn(j,:) - coord(no,:);
    end
    
    %---------------------------------------------------------------------%
    % Construção da malha híbrida ----------------------------------------%
    %---------------------------------------------------------------------%
    
    em = zeros(size(v,1),size(esurn,2)); % Matriz com os elementos que estão 
                                     % entre uma fratura e outra, incicando
                                     % que nós novos eles receberão.
    fm = zeros(size(v,1),2); % Matriz com as fraturas que receberão cada novo
                             % nó.
    
    % Determina os vetores que vão deslocar o ponto original -------------%
    if size(v,1)>1 || any(bnf==fratnodes(i))
        [ vd, em ] = deslocvec( fractreg, esurn, v, vc, any(bnf==fratnodes(i)) );
        % Verifica em que fraturas ficarão os nós novos.
        [ fm ] = quadrverify( v, vd, fractsurn, any(bnf==fratnodes(i)));
        if size(vd,1)>2
            for j=2:size(fm,1)
                for k=1:size(fm,2)
                    for h=1:size(fractn,2)
                        if fractn(fm(j,k),h)==no && isempty(find(fm(1,:)==fm(j,k)))==1
                            fractn(fm(j,k),h)=0;
                        end 
                    end
                end
            end
        elseif size(vd,1)==2
            po2(ipo2,1) = fratnodes(i);
            po2(ipo2,2) = size(coord,1)+g;
            ipo2 = ipo2 + 1;
        end
        %-----------------------------------------------------------------%
        % O ponto com o nome original é aquele que foi deslocado pelo
        % primeiro vetor da matriz.
        o=coord(no,:);
        coord(no,:)=o+vd(1,:); % Ponto com o nome original.           
        for j=2:size(vd,1)
            coordn(g,:)=o+vd(j,:); % Ponto com nome novo.
            % Ajusta, nos elementos correspondentes, o nome novo do ponto.
            for k=1:size(em,2)
               if em(j,k)~=0
                  if elemn(em(j,k),1)==0
                      elemn(em(j,k),:)=elem(em(j,k),:);
                  end
                  f=find(elem(em(j,k),:)==no);
                  elemn(em(j,k),f)=nc+g;
               end
            end
            % Ajusta, nas fraturas correspondentes, o nome novo do ponto.
            if size(vd,1)==2
                ff1=find(fract(:,1)==no);
                ff2=find(fract(:,2)==no);
                if fractn(ff1,3)==0
                    fractn(ff1,3)=nc+g;
                else
                    fractn(ff1,4)=nc+g;
                end
                if fractn(ff2,3)==0
                    fractn(ff2,3)=nc+g;
                else
                    fractn(ff2,4)=nc+g;
                end
            elseif size(vd,1)>2
                for k=1:size(fm,2)
                    for h=1:size(fractn,2)
                        if fractn(fm(j,k),h)==0
                             fractn(fm(j,k),h)=nc+g;
                             break
                        end
                    end
                end
                if j==2
                   nwfr=zeros(1,size(fract,2));
                   nwfr(1,1)=no;
                   nwfr(1,2)=nc+g;
                elseif nwfr(1,j)==0
                   nwfr(1,j)=nc+g;
                end
            end
            g=g+1;
        end
        if size(vd,1)>2 
            noz=sum(nwfr~=0);
            noo=min(nwfr(1:noz));
            np=find(spclnodes(:,3)==noo);
            nf=spclnodes(np,2);
            nwfr(5)=nf;
            fractn(size(fractn,1)+1,:)=nwfr;
            confcent(tc,1)=size(fractn,1);
            tc=tc+1;
        end
    end
end

coord(size(coord,1)+1:size(coord,1)+size(coordn,1),:)=coordn;
for i=1:size(elemn,1)
   if elemn(i,1)~=0
      elem(i,:)=elemn(i,:); 
   end
end

fract=fractn;
fract(:,5)=0;
sk=max(region);
regionfract=regionfract-2000+sk;
for i=1:size(kfrat,1)
    kmap(sk+i,:)=[sk+i kfrat(i,1)*kmap(kfrat(i,2),2) kfrat(i,1)*kmap(kfrat(i,2),3)...
                       kfrat(i,1)*kmap(kfrat(i,2),4) kfrat(i,1)*kmap(kfrat(i,2),5)];
end

if confcent(1)~=0
    for i=1:size(confcent,1)
        rfn=zeros(size(fract,1),1);
        q=sum(fract(confcent(i),:)~=0);
        x1=coord(fract(confcent(i),1),:); 
        x2=coord(fract(confcent(i),2),:);
        x3=coord(fract(confcent(i),3),:);
        if q==4
           x4=coord(fract(confcent(i),4),:);
        else
           x4=[0 0 0];
        end
        fe=(x1+x2+x3+x4)/q;
        coord(size(coord,1)+1,:)=fe;
        sfn1=sum(fract(confcent(i),1:4)~=0);
        for j=1:size(fract,1)
            sfn2=sum(fract(j,1:4)~=0);
            for y=1:sfn1
                if isempty(find(fract(j,:)==fract(confcent(i),y)))==0 && rfn(j)==0
                   fract(j,sfn2+1)=size(coord,1);
                   rfn(j)=1;
                end
            end
        end
        fract(confcent(i),:)=zeros(1,5);
    end
end

for i=1:size(fract,1)
    for j=1:5
        if fract(i,j)~=0
            x(j,:)=coord(fract(i,j),:); 
        else
            x(j,:)=[0 0 0];
        end
    end
    q=sum(fract(i,:)~=0);
    if q~=0
        fe=(x(1,:)+x(2,:)+x(3,:)+x(4,:)+x(5,:))/q;
        t=zeros(1,q);
        v(1,:)=x(1,:)-fe;
        for j=2:q
            v(j,:)=x(j,:)-fe;
            uu(j,:)=cross(v(1,:),v(j,:));
            u(j)=uu(j,3);
            t(j)=acos(dot(v(1,:),v(j,:))/(norm(v(1,:))*norm(v(j,:))));
            if u(j)<0,t(j)=(2*pi)-t(j);end
        end
        tq=0;fract2(i,1)=fract(i,1);t1=t;
        for j=1:q-1
            mf=max(t1(1:q));pf=find(t1(1:q)==mf);
            fract2(i,q-tq)=fract(i,pf);t1(pf)=0;tq=tq+1;
        end    
    end
end

clear fract
fract = fract2;

numc = size(coord,1);
fractn = zeros(size(fract,1),6);
fractn(1:size(fract,1),1:size(fract,2)) = fract;
for i=1:size(po2,1)
   coord(size(coord,1)+1,:) = 0.5*(coord(po2(i,1),:)+coord(po2(i,2),:));
   for j=1:size(fractn,1)
      if isempty(intersect(fractn(j,:),po2(i,:)))==0 && all(po2(i,:)==intersect(fractn(j,:),po2(i,:)))
         fractn(j,sum(fractn(j,:)~=0)+1)=size(coord,1);
         flagcorr(i,:) = [po2(i,1) size(coord,1)];
      end
   end
end

for ifractn = 1:(size(fractn,1))
    %"elemnode" receives three or four nodes which constitute the element
    fractnode = fractn(ifractn,1:sum(fractn(ifractn,:) ~= 0));
    %Calculate the centroid
    centfrat(ifractn,:) = mean(coord(fractnode,:));
end  %End of FOR (each element)

[ fractn ] = reordnodelem( fractn, centfrat, coord );

if isempty(find(fractn(:,6)~=0))==0
    fractn2 = fractn;
    regionfract2 = regionfract;
    lin = find(fractn(:,6)~=0);
    for i=1:size(lin,1)
        col = find(fractn(lin(i),:)>numc);
        arr = [fractn(lin(i),:) fractn(lin(i),:)];
        fractn2(lin(i),:) = zeros(1,6);
        fractn2(lin(i),1:4) = arr(col(1):col(1)+3);
        fractn2(size(fractn2,1)+1,1:4) = arr(col(1)+3:col(1)+6); regionfract2(size(fractn2,1))=regionfract(lin(i));
    end
    clear fract regionfract
    fract = fractn2;
    regionfract = regionfract2;
else
    fract = fractn(:,1:5);
end

end
