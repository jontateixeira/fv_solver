function [bedge,inedge,knownb] = getinedgebedge(coord,elem,bedgeaux,flaglim,flagcorr)
%
knownb = 0;
faces = zeros(1,4);
bedge = zeros(1,5);
auxboundflag = zeros(size(coord,1),1);

for i=1:size(bedgeaux,1)
    for j=1:2
        auxboundflag(bedgeaux(i,j),1)=bedgeaux(i,5);
    end
end
for i=1:size(flaglim,1)
    auxboundflag(i,1)=flaglim(i);
end
if flagcorr(1)~=0
    for i=1:size(flagcorr,1)
        auxboundflag(flagcorr(i,2),1)=auxboundflag(flagcorr(i,1));
    end
end

t=1;
for i=1:size(elem,1)
    n=sum(elem(i,:)~=0);
    for j=1:n
        if j==n, j1=1; else j1=j+1; end
        f=[elem(i,j) elem(i,j1) i];  
        if isempty(intersect(find(faces(:,1)==f(2)),find(faces(:,2)==f(1))))==0
            c=intersect(find(faces(:,1)==f(2)),find(faces(:,2)==f(1)));   
            faces(c,4)=f(3);
        else
            faces(t,1:3)=f;
            t=t+1;
        end
    end
end

ti=1;tb=1;
for i=1:size(faces,1)
    if faces(i,4)==0
        bedge(tb,1:4)=faces(i,:);
        tb=tb+1;
    else
        inedge(ti,:)=faces(i,:);
        ti=ti+1;
    end
end

for i=1:size(bedge,1)
   if isempty(intersect(find(bedgeaux(:,1)==bedge(i,2)),find(bedgeaux(:,2)==bedge(i,1))))==0 || ... 
       isempty(intersect(find(bedgeaux(:,2)==bedge(i,2)),find(bedgeaux(:,1)==bedge(i,1))))==0
       c=[intersect(find(bedgeaux(:,1)==bedge(i,2)),find(bedgeaux(:,2)==bedge(i,1)));...
          intersect(find(bedgeaux(:,2)==bedge(i,2)),find(bedgeaux(:,1)==bedge(i,1)))];  
       bedge(i,4:5)=[bedgeaux(c,5) bedgeaux(c,5)];
   end 
end

for i=1:size(auxboundflag,1)
    c1=find(bedge(:,1)==i);
    if auxboundflag(i)~=0
        bedge(c1,4)=auxboundflag(i);
    end
end

end

