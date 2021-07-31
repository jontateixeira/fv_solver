function [ elem, bedgeaux, nbe, fract, coord ] = vecdeslnodebound( coord, elem,...
                                                      fract, afrat, bedgeaux, nbe )
%
coordn=coord;
bedgeauxn=bedgeaux;

N=size(elem,1)-size(fract,1);
nd=max(elem(1:N,5));
t=1;
fratnodes=zeros(1,1);
for i=1:size(fract,1)
    for j=1:4
        n=fract(i,j);
        if n~=0 && isempty(find(fratnodes==n))~=0
            fratnodes(t,1)=n;
            t=t+1;
        end
    end
end

y=1;
for i=1:size(fratnodes,1)
    p=find(bedgeaux(:,1)==fratnodes(i));
    if isempty(p)==0
        if p==1
            a=size(bedgeaux,1)+1;
            g=a-1;
        else
            a=p-1;
            g=a;
        end
        n1=bedgeaux(p,1);
        n2=bedgeaux(p,2);
        n0=bedgeaux(g,1);
        
        f1=find(fract(:,1)==n1);
        f2=find(fract(:,2)==n1);
        f3=find(fract(:,3)==n1);
        if isempty(f1)==0
            C=f1;
            L=fract(C,5)-nd;
        elseif isempty(f2)==0
            C=f2;
            L=fract(C,5)-nd;
        elseif isempty(f3)==0
            C=f3;
            L=fract(C,5)-nd;
        end
        
        v=coord(n2,:)-coord(n1,:);
        vd=0.5*(v/norm(v))*abs(afrat(L));
        
        coordn(n1,:)=coord(n1,:)+vd;
        coordn(size(coord,1)+y,:)=coord(n1,:)-vd;
        
        bedgeauxn(p+1:size(bedgeaux,1)+1,:)=bedgeaux(p:size(bedgeaux,1),:);
        bedgeauxn(p,:)=[size(coord,1)+y n1 0 0 bedgeaux(p,5)];
        bedgeauxn(a,2)=size(coord,1)+y;
        
        h=N+C;
        
        fract(C,4)=size(coord,1)+y;
        
        coord=coordn;
        bedgeaux=bedgeauxn;
        
        x1=coord(fract(C,1),:); 
        x2=coord(fract(C,2),:);
        x3=coord(fract(C,3),:);
        x4=coord(fract(C,4),:);
        fe=(x1+x2+x3+x4)/4;
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
            nv=fract(C,2);
            tt=t2;
            fract(C,2)=fract(C,3);
            t2=t3;
            fract(C,3)=nv;
            t3=tt;
        end
        if fract(C,4)~=0
            v4=x4-fe;
            uu4=cross(v1,v4);
            u4=uu4(3);
            t4=acos(dot(v1,v4)/(norm(v1)*norm(v4)));
            if u4<0,t4=(2*pi)-t4;end
            if t4<t3
                nv=fract(C,3);
                tt=t3;
                fract(C,3)=fract(C,4);
                t3=t4;
                fract(C,4)=nv;
                t4=tt;
                if t3<t2
                    nv=fract(C,2);
                    tt=t2;
                    fract(C,2)=fract(C,3);
                    t2=t3;
                    fract(C,3)=nv;
                    t3=tt;
                end
            end
        end
        
        elem(h,:)=fract(C,:);
        
        for k=1:size(elem,1)
            if elem(k,4)==0
                q=3;
            else
                q=4;
            end                
            if (elem(k,1)==n1 && elem(k,q)==n0)||...
               (elem(k,2)==n1 && elem(k,1)==n0)||...
               (elem(k,3)==n1 && elem(k,2)==n0)||...
               (elem(k,q)==n1 && elem(k,q-1)==n0)
               d=find(elem(k,:)==n1);
               elem(k,d)=size(coord,1);
            end
        end 
        
    end
end

nbe=size(bedgeaux,1);

end

