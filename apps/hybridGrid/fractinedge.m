function [ fedge ] = fractinedge( inedge, fract )
%
fedge=zeros(size(inedge,1),1);

for i=1:size(fract,1)
   if fract(i,1)>fract(i,2)
       no1=fract(i,2);
       no2=fract(i,1);
   else
      no1=fract(i,1);
      no2=fract(i,2);
   end
   c1=inedge(:,1);
   f1=find(c1==no1);
   e1=inedge(f1,1:2);
   c2=e1(:,2);
   f2=find(c2==no2);
   fedge(f1(f2))=1;       
end

end

