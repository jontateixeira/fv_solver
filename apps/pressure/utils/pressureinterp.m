function [pressurinterp]=pressureinterp(p,nflag,w,s)

global coord esurn1 esurn2

%% interpolação das pressões nos nós
for no=1:size(coord,1)
    
    nec1=esurn2(no+1)-esurn2(no);
    p1=0;
    
    if nflag(no,1) > 200
        if nflag(no,1) > 201 && nflag(no,1) < 300
            for j=1:nec1
                element1=esurn1(esurn2(no)+j);
                p1=p1+w(esurn2(no)+j)*p(element1);
            end
            p1=p1+s(no,1);
        else
            for j=1:nec1
                element1=esurn1(esurn2(no)+j);
                p1=p1+w(esurn2(no)+j)*p(element1);
            end
        end

    else
        p1=nflag(no,2);
    end

    pressurinterp(no,1)=p1;

end

end
