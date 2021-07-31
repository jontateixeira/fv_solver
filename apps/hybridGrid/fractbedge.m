function [ bedge, inedge, elem, coord, nnode, nbe, centelem ] = fractbedge( bedge,...
                           inedge, elem, flaglim, coord, fract, elem_post, centelem, region )
%
global afrat

t=1; bnf=0; ebnf=0;
for i=1:size(fract,1)
    for j=1:4
        if fract(i,j)~=0 && fract(i,j)<size(flaglim,1)
            bnf(t)=fract(i,j);
            ebnf(t)=i;
            t=t+1;
        end
    end
end
nesf=size(elem_post,1);

t=size(coord,1);
r=ones(size(coord,1),1);
for i=1:size(bedge,1)
    reg = region(bedge(i,3));
    if isempty(find(bnf==bedge(i,1)))==0
        vd=coord(bedge(i,2),:)-coord(bedge(i,1),:);
        vd=0.5*afrat(reg)*vd/norm(vd);
        coord(bedge(i,1),:)=coord(bedge(i,1),:)+vd;
        r(bedge(i,1))=2;
    elseif isempty(find(bnf==bedge(i,2)))==0
        vd=coord(bedge(i,1),:)-coord(bedge(i,2),:);
        vd=r(bedge(i,2))*0.5*afrat(reg)*vd/norm(vd);
        coord(t+1,:)=coord(bedge(i,2),:)+vd;
        pbe=find(elem(bedge(i,3),:)==bedge(i,2));
        elem(bedge(i,3),pbe)=t+1;
        piee=find(inedge(:,3)==bedge(i,3));
        pied=find(inedge(:,4)==bedge(i,3));
        for k=1:size(piee,1)
            if inedge(piee(k),1)==bedge(i,2), inedge(piee(k),1)=t+1; end;
            if inedge(piee(k),2)==bedge(i,2), inedge(piee(k),2)=t+1; end;
            if inedge(piee(k),1)>inedge(piee(k),2),
               di=inedge(piee(k),2);
               inedge(piee(k),2)=inedge(piee(k),1);inedge(piee(k),1)=di;
               de=inedge(piee(k),3);inedge(piee(k),3)=inedge(piee(k),4);
               inedge(piee(k),4)=de;
            end
        end
        for k=1:size(pied,1)
            if inedge(pied(k),1)==bedge(i,2), inedge(pied(k),1)=t+1; end;
            if inedge(pied(k),2)==bedge(i,2), inedge(pied(k),2)=t+1; end;
            if inedge(pied(k),1)>inedge(pied(k),2)
               di=inedge(pied(k),2);
               inedge(pied(k),2)=inedge(pied(k),1);inedge(pied(k),1)=di;
               de=inedge(pied(k),3);inedge(pied(k),3)=inedge(pied(k),4);
               inedge(pied(k),4)=de;
            end
        end
        fbf=find(bnf==bedge(i,2)); 
        elem(ebnf(fbf)+nesf,4)=t+1;
        bedge(size(bedge,1)+1,:)=[t+1 bedge(i,2) ebnf(fbf)+nesf bedge(i,5) bedge(i,5)];
        bedge(i,2)=t+1;
        t=t+1;
    end
end

nnode = size(coord,1);
nbe = size(bedge,1);

for i=nesf:size(elem,1)
   v1=coord(elem(i,1),:)-centelem(i,:);
   if elem(i,4)==0
       noe=elem(i,1:3)';
       C=coord(noe,:);
       centelem(i,:)=(C(1,:)+C(2,:)+C(3,:))/3;
       n=3;
   else
       noe=elem(i,1:4)';
       C=coord(noe,:);
       centelem(i,:)=(C(1,:)+C(2,:)+C(3,:)+C(4,:))/4;
       n=4;
   end
   for j=1:n
       v2=coord(elem(i,j),:)-centelem(i,:);
       bv=cross(v1,v2);
       b=bv(3);
       theta(j)=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
       if b<0
          theta(j)=(2*pi)-theta(j); 
       end
   end 
   ord=zeros(size(theta,2),1);
   theta1=theta;
   for k=1:size(theta,2)
       pos=find(theta1==min(theta1));
       e1(k)=elem(i,pos);
       theta1(pos)=2*pi;
   end
   elem(i,1:n)=e1;
   clear ord theta e1
end

end

