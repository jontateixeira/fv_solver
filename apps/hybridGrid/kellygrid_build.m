function [ coord, elem, fract, kmap, region, regionfract ] = kellygrid_build( coord, elem, ...
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
e=0;v=0;vd=0;vc=0;ce=0;n=0;f=0;c=0;d=0;ff=0;l=0;fm=0;z=1;

% Percorre os nós que pertencem às fraturas ------------------------------%
for i=1:size(fratnodes,1)
    
    no=fratnodes(i);
    
    clear e v vd vc ce n f c d ff l fm
    
    % Verifica quais são os elementos que compartilham o nó em questão ---%
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
    % Verifica os centróides desses elementos ----------------------------%
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
    % Verifica a que fraturas pertencem esses nós ------------------------%    
    c1=find(fract(:,1)==no); % Fratura
    c2=find(fract(:,2)==no); % Fratura
    n1=fract(c1,2); % Outro nó da Fratura
    n2=fract(c2,1); % Outro nó da Fratura
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
    % Vetor nó-centróide das fraturas ------------------------------------%
    for j=1:size(e,2)
        vc(j,:)=ce(j,:)-coord(no,:);
    end
    
    %---------------------------------------------------------------------%
    % Construção da malha híbrida ----------------------------------------%
    %---------------------------------------------------------------------%
    
    em=zeros(size(v,1),size(e,2)); % Matriz com os elementos que estão 
                                   % entre uma fratura e outra, incicando
                                   % que nós novos eles receberão.
    fm=zeros(size(v,1),2); % Matriz com as fraturas que receberão cada novo
                           % nó.
    
    % Determina os vetores que vão deslocar o ponto original -------------%
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
            vd(y,:)=0.5*((R*v(y,:)')/norm(v(y,:)))*abs(max(afrat(l))/sin(theta/2));
            % 'vd' é a matriz dos vetores que deslocam o ponto original.
        end
        % Verifica em que fraturas ficarão os nós novos.
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
        % O ponto com o nome original é aquele que foi deslocado pelo
        % primeiro vetor da matriz.
        o=coord(no,:);
        coordk(z,:)=o;
        correspk(z,1)=no;
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
            % SÓ MUDA ESSA PARTE DAS FRATURAS!!! -------------------------%
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
            %-------------------------------------------------------------%
            g=g+1;
        end
        if size(vd,1)>2 
            noz=sum(nwfr~=0);
            noo=min(nwfr(1:noz));
            np=find(spclnodes(:,3)==noo);
            nf=spclnodes(np,2);
            nwfr(5)=nf;
            fractn(size(fractn,1)+1,:)=nwfr;
            regionfract(size(fractn,1),:)=nwfr(5);
        end
        z=z+1;
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

% MODIFICA O FRACTN DAS FRATURAS -----------------------------------------%
sc=size(coord,1);
coord(size(coord,1)+1:size(coord,1)+size(coordk,1),:)=coordk;
fractmod = zeros(size(fract,1),6);
fractmod(1:size(fract,1),1:size(fract,2)-1) = fractn(1:size(fract,1),1:size(fract,2)-1);
fractmod2 = fractmod;
for i=1:size(fract,1)
   for j=1:2
      if isempty(find(correspk==fract(i,j)))==0 && fractmod(i,5)==0
          fractmod(i,5)=sc+find(correspk==fract(i,j));
      elseif isempty(find(correspk==fract(i,j)))==0
          fractmod(i,6)=sc+find(correspk==fract(i,j));
      end
   end
end
for i=1:size(fractmod,1)
   if fractmod(i,6)~=0
       a=fractmod(i,:)';
       c=sum(coord(a,:))/6;
       vn0=coord(fractmod(i,1),:)-c;
       vn0=vn0/norm(vn0);
       for j=1:6
           vn=coord(fractmod(i,j),:)-c;
           vn=vn/norm(vn);
           an(j)=acos(dot(vn0,vn));
           cn=cross(vn0,vn);
           c0(j)=cn(3);
           if c0(j)<0
              an(j)=(2*pi)-an(j);
           end
       end
       ami=min(an(5),an(6));
       ama=max(an(5),an(6));
       t1=1;t2=1;
       for j=1:6
           if an(j)>ami && an(j)<ama
              fratnova1(1,t1)=fractmod(i,j);
              t1=t1+1;
           else
              fratnova2(1,t2)=fractmod(i,j);
              t2=t2+1;
           end
       end
       fratnova1(3:4)=fractmod(i,5:6);
       fractmod2(i,:)=zeros(1,6);
       fractmod2(size(fractmod2,1)+1,:)=zeros(1,6);
       fractmod2(i,1:4)=fratnova1;
       fractmod2(size(fractmod2,1),1:4)=fratnova2;
       fractmod2(size(fractmod2,1),6)=fract(i,5);
       regionfract(size(fractmod2,1))=regionfract(i);
   elseif fractmod(i,4)==0
       fractmod2(i,4)=fractmod(i,5);
   end
end
fractn1(:,1:4)=fractmod2(:,1:4);
fractn1(:,5)=fractmod2(:,6);
fractn1(1:size(fract,1),5)=fract(:,5);
fractn1(:,5)=0;
%-------------------------------------------------------------------------%

% regionfract(size(fract,1)+1:2*size(fract,1)-2)=regionfract(2:size(fract,1)-1);
fract=fractn1;

for i=1:size(kfrat,1)
    kmap(sk+i,:)=[sk+i kfrat(i,1)*kmap(kfrat(i,2),2) kfrat(i,1)*kmap(kfrat(i,2),3)...
                       kfrat(i,1)*kmap(kfrat(i,2),4) kfrat(i,1)*kmap(kfrat(i,2),5)];
end

for i=1:size(fract,1)
    x(1,:)=coord(fract(i,1),:); 
    x(2,:)=coord(fract(i,2),:);
    x(3,:)=coord(fract(i,3),:);
    if fract(i,4)~=0
        q=4;
        x(4,:)=coord(fract(i,4),:);
    else
        q=3;
        x(4,:)=[0 0 0];
    end   
    fe=sum(x)/q;
    vr(1,:)=x(1,:)-fe;
    vr(1,:)=vr(1,:)/norm(vr(1,:));
    ar(1)=0;
    for j=2:q
        vr(j,:)=x(j,:)-fe;
        vr(j,:)=vr(j,:)/norm(vr(j,:));
        cr1=cross(vr(1,:),vr(j,:));
        cr(j)=cr1(3);
        ar(j)=acos(dot(vr(1,:),vr(j,:)));
        if cr(j)<0
            ar(j)=(2*pi)-ar(j);
        end
    end
    for j=1:4
        ap=find(ar==min(ar));
        fr1(j)=fract(i,ap);
        ar(ap)=7;
    end 
    fract(i,1:4)=fr1;
    
    x(1,:)=coord(fract(i,1),:); 
    x(2,:)=coord(fract(i,2),:);
    x(3,:)=coord(fract(i,3),:);
    if fract(i,4)~=0
        x(4,:)=coord(fract(i,4),:);
    end   
    
    if fract(i,4)~=0
        v12=x(2,:)-x(1,:);v14=x(4,:)-x(1,:);v12=v12/norm(v12);v14=v14/norm(v14);td(1)=acos(dot(v12,v14));
        v23=x(3,:)-x(2,:);v21=x(1,:)-x(2,:);v23=v23/norm(v23);v21=v21/norm(v21);td(2)=acos(dot(v21,v23));
        v34=x(4,:)-x(3,:);v32=x(2,:)-x(3,:);v34=v34/norm(v34);v32=v32/norm(v32);td(3)=acos(dot(v34,v32));
        v41=x(1,:)-x(4,:);v43=x(3,:)-x(4,:);v41=v41/norm(v41);v43=v43/norm(v43);td(4)=acos(dot(v43,v41));
        for j=1:4
            if abs(td(j)-pi)<1e-6
                vd=x(j,:)-fe;
                vd=vd*1e-3;
                coord(fract(i,j),:)=coord(fract(i,j),:)+vd;
            end
        end
    else
        v12=x(2,:)-x(1,:);v13=x(3,:)-x(1,:);v12=v12/norm(v12);v13=v13/norm(v13);td(1)=acos(dot(v12,v13));
        v23=x(3,:)-x(2,:);v21=x(1,:)-x(2,:);v23=v23/norm(v23);v21=v21/norm(v21);td(2)=acos(dot(v21,v23));
        v31=x(1,:)-x(3,:);v32=x(2,:)-x(3,:);v31=v31/norm(v31);v32=v32/norm(v32);td(3)=acos(dot(v31,v32));
        for j=1:3
            if abs(td(j)-pi)<1e-6
                vd=x(j,:)-fe;
                vd=vd*1e-3;
                coord(fract(i,j),:)=coord(fract(i,j),:)+vd;
            end
        end
    end   
    
end

end

