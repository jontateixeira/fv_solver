function [ vd, em ] = deslocvec( fr, e, v, vc, veri )
%
global afrat

if veri, v(2,:)=-v(1,:); end

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
    vd(y,:)=0.5*((R*v(y,:)')/norm(v(y,:)))*abs(max(afrat(fr))/sin(theta/2));
    % 'vd' é a matriz dos vetores que deslocam o ponto original.
 end

end

