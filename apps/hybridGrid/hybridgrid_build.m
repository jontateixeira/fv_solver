function [ coord, elem, fract, kmap, region, regionfract ] = hybridgrid_build( coord, elem, ...
           fract, afrat, kfrat, kmap, spclnodes, region, regionfract )
%

elemn=zeros(size(elem,1),size(elem,2));
nc=size(coord,1);
g=1;

t=1;
fratnodes=zeros(1,1);
for i=1:size(fract,1)
    for j=1:2
        n=fract(i,j);
        if isempty(find(fratnodes==n))~=0
            fratnodes(t,1)=n;
            t=t+1;
        end
    end
    nnf(i,1)=sum(fract(:,1)==fract(i,1))+sum(fract(:,2)==fract(i,2))...
             +sum(fract(:,1)==fract(i,2))+sum(fract(:,2)==fract(i,1));
    cf(i,:)=0.5*(coord(fract(i,1),:)+coord(fract(i,2),:));
end

fractn=fract;
ino=zeros(size(fract,1),1);
e=0;v=0;vd=0;vc=0;ce=0;n=0;f=0;c=0;d=0;ff=0;l=0;fm=0;

% Percorre os n�s que pertencem �s fraturas ------------------------------%
for i=1:size(fratnodes,1)
    
    no=fratnodes(i);
    
    clear e v vd vc ce n f c d ff l fm
    
    % Verifica quais s�o os elementos que compartilham o n� em quest�o ---%
    e1=find(elem(:,1)==no);
    e2=find(elem(:,2)==no);
    e3=find(elem(:,3)==no);
    e4=find(elem(:,4)==no);
    t=1;
    if isempty(e1)==0
        for j=1:size(e1,1)
            e(t)=e1(j);
            t=t+1;
        end
    end
    if isempty(e2)==0
        for j=1:size(e2,1)
            e(t)=e2(j);
            t=t+1;
        end
    end
    if isempty(e3)==0
        for j=1:size(e3,1)
            e(t)=e3(j);
            t=t+1;
        end
    end
    if isempty(e4)==0
        for j=1:size(e4,1)
            e(t)=e4(j);
            t=t+1;
        end
    end
    % Verifica os centr�ides desses elementos ----------------------------%
    for j=1:size(e,2)
       x1=coord(elem(e(j),1),:); 
       x2=coord(elem(e(j),2),:); 
       x3=coord(elem(e(j),3),:); 
       if elem(e,4)~=0
           x4=coord(elem(e(j),4),:); 
           q=4;
       else
           x4=[0 0 0];
           q=3;
       end
       ce(j,:)=(x1+x2+x3+x4)/q;
    end
    % Verifica a que fraturas pertencem esses n�s ------------------------%    
    c1=find(fract(:,1)==no); % Fratura
    c2=find(fract(:,2)==no); % Fratura
    n1=fract(c1,2); % Outro n� da Fratura
    n2=fract(c2,1); % Outro n� da Fratura
    t=1;
    if isempty(n1)==0
        for j=1:size(n1,1)
            n(t)=n1(j);
            l(t,1)=c1(j);
            t=t+1;
        end
    end
    if isempty(n2)==0
        for j=1:size(n2,1)
            n(t)=n2(j);
            l(t,1)=c2(j);
            t=t+1;
        end
    end
    ff=l;
    l=regionfract(l);
    l=l-2000;
    % Vetor Fratura ------------------------------------------------------%
    for j=1:size(n,2)
        v(j,:)=coord(n(j),:)-coord(no,:);
    end
    % Vetor n�-centr�ide das fraturas ------------------------------------%
    for j=1:size(e,2)
        vc(j,:)=ce(j,:)-coord(no,:);
    end
    
    %---------------------------------------------------------------------%
    % Constru��o da malha h�brida ----------------------------------------%
    %---------------------------------------------------------------------%
    
    em=zeros(size(v,1),size(e,2)); % Matriz com os elementos que est�o 
                                   % entre uma fratura e outra, incicando
                                   % que n�s novos eles receber�o.
    fm=zeros(size(v,1),2); % Matriz com as fraturas que receber�o cada novo
                           % n�.
    
    % Determina os vetores que v�o deslocar o ponto original -------------%
    if size(v,1)>1
        for y=1:size(v,1)
            d(y)=0;
            c(y)=0;
            for j=1:size(v,1)
                if y~=j
                    cl=cross(v(y,:),v(j,:));
                    c(j)=cl(3);
                    d(j)=acos(dot(v(y,:),v(j,:))/(norm(v(j,:))*norm(v(y,:))));
                    if c(j)<0
                        d(j)=(2*pi)-d(j);
                    end
                end
            end
            if y==size(v,1)
               cl=cross(v(y,:),v(1,:));
               c(y)=cl(3);
               d(y)=acos(dot(v(y,:),v(1,:))/(norm(v(1,:))*norm(v(y,:))));
               if c(y)<0
                  d(y)=(2*pi)-d(y);
               end
            end
            theta=7;
            for j=1:size(d,2)
                if d(j)~=0 && d(j)<theta
                    theta=d(j);
                end
            end
            s=1;
            for j=1:size(e,2)
               e1=cross(v(y,:),vc(j,:));
               ec=e1(3);
               ed=acos(dot(v(y,:),vc(j,:))/(norm(vc(j,:))*norm(v(y,:))));
               if ec<0
                  ed=(2*pi)-ed; 
               end
               if ed<theta
                   em(y,s)=e(j);
                   s=s+1;
               end
            end
            theta=real(theta);
            R=[cos(theta/2) -sin(theta/2) 0;sin(theta/2) cos(theta/2) 0;0 0 1];
            vd(y,:)=0.5*((R*v(y,:)')/norm(v(y,:)))*abs(afrat(l(y))/sin(theta/2));
            % 'vd' � a matriz dos vetores que deslocam o ponto original.
        end
        % Verifica em que fraturas ficar�o os n�s novos.
        [ fm ] = quadrverify( v, vd, ff );
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
        end
        %-----------------------------------------------------------------%
        % O ponto com o nome original � aquele que foi deslocado pelo
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
            regionfract(size(regionfract,1)+1)=nf;
            fractn(size(fractn,1)+1,:)=nwfr;
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
sk=max(region);
regionfract=regionfract-2000+sk;
for i=1:size(kfrat,1)
    kmap(sk+i,:)=[sk+i kfrat(i,1)*kmap(kfrat(i,2),2) kfrat(i,1)*kmap(kfrat(i,2),3)...
                       kfrat(i,1)*kmap(kfrat(i,2),4) kfrat(i,1)*kmap(kfrat(i,2),5)];
end

for i=1:size(fract,1)
    x1=coord(fract(i,1),:); 
    x2=coord(fract(i,2),:);
    x3=coord(fract(i,3),:);
    if fract(i,4)~=0
        x4=coord(fract(i,4),:);
        q=4;
    else
        x4=[0 0 0];
        q=3;
    end
    fe=(x1+x2+x3+x4)/q;
    v1=x1-fe;
    v2=x2-fe;
    v3=x3-fe;
    uu2=cross(v1,v2);
    uu3=cross(v1,v3);
    u2=uu2(3);
    u3=uu3(3);
    t2=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
    t3=acos(dot(v1,v3)/(norm(v1)*norm(v3)));
    if u2<0,t2=(2*pi)-t2;end
    if u3<0,t3=(2*pi)-t3;end
    if t3<t2
       nv=fract(i,2);
       tt=t2;
       fract(i,2)=fract(i,3);
       t2=t3;
       fract(i,3)=nv;
       t3=tt;
    end
    if fract(i,4)~=0
        v4=x4-fe;
        uu4=cross(v1,v4);
        u4=uu4(3);
        t4=acos(dot(v1,v4)/(norm(v1)*norm(v4)));
        if u4<0,t4=(2*pi)-t4;end
        if t4<t3
            nv=fract(i,3);
            tt=t3;
            fract(i,3)=fract(i,4);
            t3=t4;
            fract(i,4)=nv;
            t4=tt;
            if t3<t2
                nv=fract(i,2);
                tt=t2;
                fract(i,2)=fract(i,3);
                t2=t3;
                fract(i,3)=nv;
                t3=tt;
            end
        end
    end
end

end

