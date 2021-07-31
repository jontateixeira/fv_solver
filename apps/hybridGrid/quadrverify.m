function [ fm ] = quadrverify( v, vd, ff, veri )
%

for i=1:size(vd,1)
   for j=1:size(v,1)
       v1=vd(i,:)/norm(vd(i,:));
       v2=v(j,:)/norm(v(j,:));
       cv3=cross(v1,v2);
       cv(j)=cv3(3);
       dv(j)=acos(dot(v1,v2)/(norm(v1)*norm(v2)));
       if cv(j)<0
           dv(j)=2*pi-dv(j);
       end
   end
   mip=find(dv==min(dv));
   map=find(dv==max(dv));
   if veri, mip=1; map=1; end
   fm(i,:)=[ff(map) ff(mip)];
end

end

