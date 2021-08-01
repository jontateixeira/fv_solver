function [inedge,bedge] = reordbedge(inedge,bedge,coord,centelem)

for i=1:size(bedge,1)
   v1 = coord(bedge(i,2),:)-coord(bedge(i,1),:);
   v2 = centelem(bedge(i,3),:)-coord(bedge(i,1),:);
   bv=cross(v1,v2);
   b=bv(3);
   if b<0
      vc = [bedge(i,2) bedge(i,1)];
      bedge(i,1:2) = vc;
   end
end

for i=1:size(inedge,1)
    if inedge(i,1)>inedge(i,2)
        c=inedge(i,1);
        inedge(i,1)=inedge(i,2);
        inedge(i,2)=c;
        d=inedge(i,3);
        inedge(i,3)=inedge(i,4);
        inedge(i,4)=d;
    end
end

end