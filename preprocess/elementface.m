function [F,V,N]=elementface(nflag,inedge)
global bedge elem coord esurn1 esurn2 nsurn1 nsurn2
%

si=size(inedge,1);
F=zeros(size(elem,1),6);

ny=0;
for y=1:size(esurn2,1)-1
    n=esurn2(y+1)-esurn2(y);
    if n>ny
        ny=n;
    end
end

N=zeros(size(coord,1),ny);
V=zeros(2,ny,size(coord,1));

for i=1:size(elem,1)
   t=1;
   for j=1:size(inedge,1)
       if (inedge(j,3)==i)||(inedge(j,4)==i)
           F(i,t)=j;
           t=t+1;
       end
   end
   for j=1:size(bedge,1)
      if bedge(j,3)==i
          F(i,t)=j+si;
          t=t+1;
      end
   end
end 

for y=1:size(coord,1)
    if (isempty(find(bedge(:,1)==y))==0)||(isempty(find(bedge(:,2)==y))==0)
        isbedge=1; % se for de contorno.
    else
        isbedge=0; % se for interno.
    end
    ve=esurn1(esurn2(y)+1:esurn2(y+1));
    if isbedge==0
        for g=1:size(ve,1)-1
            for k=1:size(F,2)
                a=F(ve(g),k);
                f1=find(F(ve(g+1),:)==a);
                f=F(ve(g+1),f1);
                if isempty(f1)==0
                    break
                end
            end
            v(g)=f;
        end
        for k=1:size(F,2)
            a=F(ve(size(ve,1)),k);
            f1=find(F(ve(1),:)==a);
            f=F(ve(1),f1);
            if isempty(f1)==0
                break
            end
        end
        v(size(ve,1))=f;
    else
        for g=1:size(ve,1)-1
            for k=1:size(F,2)
                a=F(ve(g),k);
                f1=find(F(ve(g+1),:)==a);
                f=F(ve(g+1),f1);
                if isempty(f1)==0
                    break
                end
            end
            v(g+1)=f;
        end
        f1=find(F(ve(1),:)>size(inedge,1));
        if size(f1,2)>1
            eb=ve(1);
            pe=find(elem(eb,:)==y);
            cf1=elem(eb,pe:4);
            cf2=elem(eb,1:pe-1);
            cf=horzcat(cf1,cf2);
            bf1=bedge(F(ve(1),f1(1))-size(inedge,1),1:2);
            bf2=bedge(F(ve(1),f1(2))-size(inedge,1),1:2);
            if bf1(1)==y
                bff1=bf1(2);
            else
                bff1=bf1(1);
            end
            if bf2(1)==y
                bff2=bf2(2);
            else
                bff2=bf2(1);
            end
            ff1=intersect(cf,bff1);
            ff2=intersect(cf,bff2);
            if isempty(ff1)==0
                f2(1)=find(cf==ff1);
            else
                f2(1)=0;
            end
            if isempty(ff2)==0
                f2(2)=find(cf==ff2);
            else
                f2(2)=0;
            end
            if f2(1)>f2(2)
                f1=f1(2);
            else
                f1=f1(1);
            end 
        end
        v(1)=F(ve(1),f1);
        f1=find(F(ve(size(ve,1)),:)>size(inedge,1));
        if size(f1,2)>1
            eb=ve(size(ve,1));
            pe=find(elem(eb,:)==y);
            cf1=elem(eb,pe:4);
            cf2=elem(eb,1:pe-1);
            cf=horzcat(cf1,cf2);
            bf1=bedge(F(ve(size(ve,1)),f1(1))-size(inedge,1),1:2);
            bf2=bedge(F(ve(size(ve,1)),f1(2))-size(inedge,1),1:2);
            if bf1(1)==y
                bff1=bf1(2);
            else
                bff1=bf1(1);
            end
            if bf2(1)==y
                bff2=bf2(2);
            else
                bff2=bf2(1);
            end
            ff1=intersect(cf,bff1);
            ff2=intersect(cf,bff2);
            if isempty(ff1)==0
                f2(1)=find(cf==ff1);
            else
                f2(1)=0;
            end
            if isempty(ff2)==0
                f2(2)=find(cf==ff2);
            else
                f2(2)=0;
            end
            if f2(1)<f2(2)
                f1=f1(2);
            else
                f1=f1(1);
            end 
        end
        v(size(ve,1)+1)=F(ve(size(ve,1)),f1);
    end
    n=size(v,2);
    N(y,1:n)=v;
    if isbedge==1 % Se for de contorno
        V(1,1:n-1,y)=v(1:n-1);
        V(2,1:n-1,y)=v(2:n);
    else
        V(1,1:n,y)=v;
        V(1,n,y)=v(1);
        V(2,1:n-1,y)=v(2:n);
        V(2,n,y)=v(n);
    end   
    clear v;
end

end