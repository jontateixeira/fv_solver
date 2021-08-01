function[elem,kmap]=adequapermeab(benchmark,kmap,elem)
global centelem coord pormap

switch benchmark
    
       case 'nikitin'        
                            
         for i=1:size(centelem,1)
                    if centelem(i,2)<=(0.51-centelem(i,1))
                        elem(i,5)=i;

                        K(i,1:5)=[i 505 495 495 505];
                    elseif (0.51-centelem(i,1))<centelem(i,2) && centelem(i,2)<(0.99-centelem(i,1))
                        elem(i,5)=i;

                        K(i,1:5)= [i 1000 0 0 10];
                    elseif (0.99-centelem(i,1))<=centelem(i,2) && centelem(i,2)<(1.49-centelem(i,1))
                        elem(i,5)=i;

                        K(i,1:5)=[i 10 0 0 1000];
                    elseif (1.49-centelem(i,1))<=centelem(i,2)
                        elem(i,5)=i;

                        K(i,1:5)=[i 505 495 495 505];
                    end
          end

        
    case 'nikitin1'
                
        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);
        pormap = pormap(1)*ones(size(centelem,1),1);
        %Swept all elements:
        for i = 1:size(centelem,1)
            %Get the vertices:
            vertices = elem(i,1:4);
            %Get only the non zero values
            vertices = vertices(logical(vertices ~= 0));
            %Get the "y" coordinate of each vertex
            ycoordvtx = coord(vertices,2);
            %Calculate the maximum and minimum "y" coordinate
            ycoordmin = min(ycoordvtx);
            ycoordmax = max(ycoordvtx);
            %Define "x" and "y" (centroid of element)
            x = centelem(i,1);
            y = centelem(i,2);
            %Evaluate the position of "x" and "y"
            y_reg1 = -x + 0.5;  %-x + 0.75
            y_reg2 = -x + 1;
            y_reg3 = -x + 1.5;  %-x + 1.25
            %The element is in region 1 or 4
            if (x <= 0.5 && y <= y_reg1) || (x <= 0.5 && ...
                    y_reg1 >= ycoordmin && y_reg1 <= ycoordmax) ...
                    || (x > 0.5 && y >= y_reg3) || (x > 0.5 && ...
                    y_reg3 >= ycoordmin && y_reg3 <= ycoordmax)
                %Definition of permeability components
                k = [505 495; 495 505];
                %The element is in region 2
            elseif (x <= 0.5 && y > y_reg1 && y < y_reg2) || ...
                    (y_reg2 >= ycoordmin && y_reg2 <= ycoordmax) || ...
                    (x > 0.5 && y < y_reg2)
                %Definition of permeability components
                k = [1000 0; 0 10];
                %The element is in region 3
            else
                %Definition of permeability components
                k = [10 0; 0 1000];
            end  %End of IF
            %Build "kmap"
            K(i,:) = [i k(1,1) k(1,2) k(2,1) k(2,2)];
            elem(i,5)=i;
        end  %End of FOR
        kmap=K;
    case 'lamine'
        
        K(1,1:5)=[1 5.5 4.5 4.5 5.5];
        elem(:,5)=1;
    case 'durlofsky'
        
        K(1,1:5)=[1 1 0 0 1];
        elem(:,5)=1;
    case 'shuec1'
        for i = 1:size(centelem,1)
            epsilon= rand(1,1);
            s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
            kmap(i,1:5)=[i s 0 0 s];
            K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
        end
        
    case 'chueh' % teste case FERNANDO
        v=randperm(length(centelem));
        for i = 1:size(centelem,1)

                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end

                K(i,1:5)=[i min(max(Sum,0.01),4) 0 0 min(max(Sum,0.01),4)];
                normKmap(i)=norm([K(i,2:3); K(i,4:5)]);
        end

case 'shuec2'
    %======================================================================
    %Example 31.2: Two-Phase Flow case. Adapted from Chueh et al., 2010 
    %(Very High Permeability Distribution)
    
        %Define number of the randomic values
        N = 40;
        %Define Randomic paramiter
        randcoord = getrandist;
        %Initialize "qsi" and "randcoord"
        qsi = zeros(N,1);

        %Initialize "kmap"
        kmap = zeros(size(centelem,1),5);

        %Swept all elements
        for i = 1:size(centelem,1)
            %Swept all randomic values in order define "randcoord"
            for irand = 1:N
                %Define "randcoord"
                qsi(irand) = exp(-(norm(centelem(i,1:2) - ...
                    randcoord(irand,:))/0.05)^2);
            end  %End of FOR
            
            %Definition of permeability constant
            kconst = min(max(sum(qsi),0.05),0.71);
            
            %Build "kmap"
            kmap(i,:) = [i kconst 0 0 kconst];
            K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
        end  %End of FOR
        
    case 'spe'
        
            %kmap = load('-ascii','SPE10_35.dat');
            n_elem_x = sqrt(size(centelem,1));
            load spe10; K_=KU(:,1:n_elem_x,1:n_elem_x,layer); 
            perm_ = squeeze(K_(1,:,:)); perm = perm_(:);
            for i = 1:length(perm)
                k_map(i,:) = [i perm(i)*[1 0 0 1]];
            end
            kmap = k_map;
            %K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
%==========================================================================        
        
%         v=randperm(length(centelem));
%         savefile = 'perm_chuec2_16x16_v2.mat'; save(savefile,'v');
%         %load('-mat','perm_chuec2_16x16_v2.mat');
%         for i = 1:size(centelem,1)
%             %if centelem(i,1)<0.5
%                 Sum=0;
%                 for m=1:40
%                     xl=centelem(v(m),:);
%                     Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
%                 end
%                 
%                 kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
%                 K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            %else
                
                %epsilon= rand(1,1);
                %s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                %kmap(i,1:5)=[i s 0 0 s];
                %K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            %end
        %end
    case 'shuec3'
        v=randperm(length(centelem));
        for i = 1:size(centelem,1)
            if centelem(i,1)<0.5 && centelem(i,2)<0.5
                epsilon= rand(1,1);
                s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                kmap(i,1:5)=[i s 0 0 s];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            elseif centelem(i,1)>0.5 && centelem(i,2)<0.5
                
                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end
                
                kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            elseif centelem(i,1)<0.5 && centelem(i,2)>0.5
                
                Sum=0;
                for m=1:40
                    xl=centelem(v(m),:);
                    Sum=Sum + exp(-(norm(centelem(i,:)-xl)/0.05)^2);
                end
                
                kmap(i,1:5)=[i min(max(Sum,0.05),0.95) 0 0 min(max(Sum,0.05),0.95)];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            else
                
                epsilon= rand(1,1);
                s= 0.1+ 50*(1+sin(10*(centelem(i,1)+centelem(i,2))))*epsilon;
                kmap(i,1:5)=[i s 0 0 s];
                K(i)=norm([kmap(i,2:3); kmap(i,4:5)]);
                
            end
            
        end
    
end
kmap=K;
% normk = calcnormk(kmap);
pormap = pormap(1)*ones(size(centelem,1),1);
end
function [randcoord] = getrandist

randcoord = ...
   [0.060180000000000   0.039620000000000;
   0.337450000000000   0.025730000000000;
   0.295320000000000   0.096820000000000;
   0.513380000000000   0.039440000000000;
   0.982290000000000   0.005730000000000;
   0.014670000000000   0.226090000000000;
   0.138620000000000   0.239620000000000;
   0.005170000000000   0.468580000000000;
   0.138440000000000   0.506530000000000;
   0.016050000000000   0.509130000000000;
   0.279250000000000   0.426090000000000;
   0.279430000000000   0.303350000000000;
   0.573750000000000   0.155830000000000;
   0.477270000000000   0.354460000000000;
   0.466900000000000   0.331950000000000;
   0.562870000000000   0.238620000000000;
   0.581900000000000   0.218400000000000;
   0.091960000000000   0.862580000000000;
   0.887500000000000   0.722770000000000;
   0.498670000000000   0.678290000000000;
   0.670250000000000   0.926130000000000;
   0.512560000000000   0.928210000000000;
   0.810060000000000   0.927240000000000;
   0.971930000000000   0.936480000000000;
   0.996670000000000   0.968610000000000;
   0.451010000000000   0.818270000000000;
   0.724800000000000   0.540370000000000;
   0.733800000000000   0.586140000000000;
   0.756040000000000   0.513060000000000;
   0.950220000000000   0.405030000000000;
   0.892670000000000   0.426740000000000;
   0.922850000000000   0.387290000000000;
   0.928420000000000   0.441850000000000;
   0.887500000000000   0.246440000000000;
   0.889500000000000   0.201670000000000;
   0.339800000000000   0.652870000000000;
   0.330120000000000   0.729340000000000;
   0.273070000000000   0.722770000000000;
   0.657540000000000   0.478110000000000;
   0.609880000000000   0.582960000000000];
end