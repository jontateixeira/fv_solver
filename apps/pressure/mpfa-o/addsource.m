function [M,mvector] = addsource(M,mvector,wells)
%Define global parameters:
global elemarea

%Case there is a source term in any element
if (size(wells,2) > 1) 
    %Catch the row which report to injector well (verify in saturation 
    %flag column)
    injecrow = find(wells(:,3) ~= 0);
    %Catch the row which report to injector well (verify in saturation 
    %flag column)
    producrow = find(wells(:,3) == 0);
    
    %----------------------------------------------------------------------
    %Apply FLOWRATE or PRESSURE in the INJECTOR well
    
    %1. Applying FLOWRATE:
    if any(wells(injecrow,5) == 0)
        %Define the value of flowrate
        flowratevalue = wells(injecrow,6);
        %Sum all element volumes inside well
        sumelemwell = sum(elemarea(wells(injecrow,1)));
        %Swept all eleemnts with source term
        for inj = 1:length(injecrow)
            %Apply the source to "mvector" (algebric system)
            mvector(wells(injecrow(inj),1)) = ...
                mvector(wells(injecrow(inj),1)) + ...
                (flowratevalue(inj)*elemarea(wells(injecrow(inj),1))/...
                sumelemwell);
        end  %End of FOR (flowrate in injector well)
    %2. Applying PRESSURE:
    elseif any(wells(injecrow,5) > 400 & wells(injecrow,5) < 501)
        %Define the value of pressure
        presvalue = wells(injecrow,6);
        %Swept all eleemnts with source term
        for inj = 1:length(injecrow)
            %Add to independent vector ("mvector") the terms of matrix "M"
            %associated to knouwn pressure
            mvector = ...
                mvector - M(:,wells(injecrow(inj),1))*presvalue(inj);  
            %Attribute the pressure value to position in the independent 
            %vector
            mvector(wells(injecrow(inj),1)) = presvalue(inj);  
            %Null both row and column of global matrix:
            %Null column
            M(:,wells(injecrow(inj),1)) = 0;
            %Null row
            M(wells(injecrow(inj),1),:) = 0;
            %Put 1 in the position (i,i)
            M(wells(injecrow(inj),1),wells(injecrow(inj),1)) = 1;
        end  %End of FOR (injector)
    end  %End of IF (injector)

    %----------------------------------------------------------------------
    %Apply FLOWRATE or PRESSURE in the PRODUCER well
    
    %1. Applying FLOWRATE:
    if any(wells(producrow,5) == 0)
        %Define the value of flowrate
        flowratevalue = wells(producrow,6);
        %Sum all element volumes inside well
        sumelemwell = sum(elemarea(wells(producrow,1)));
        %Swept all elements with source term
        for iprod = 1:length(producrow)
            %Apply the source to "mvector" (algebric system)
            mvector(wells(producrow(iprod),1)) = ...
                mvector(wells(producrow(iprod),1)) + ...
                (flowratevalue(iprod)*...
                elemarea(wells(producrow(iprod),1))/sumelemwell);
        end  %End of FOR (flowrate in producer well)
    %2. Applying PRESSURE:
    elseif any(wells(producrow,5) > 500 & wells(producrow,5) < 601)
        %Define the value of pressure
        presvalue = wells(producrow,6);
        %Swept all eleemnts with source term
        for iprod = 1:length(producrow)
            %Add to independent vector ("mvector") the terms of matrix "M"
            %associated to knouwn pressure
            mvector = mvector - ...
                M(:,wells(producrow(iprod),1))*presvalue(iprod);  
            %Attribute the pressure value to position in the independent 
            %vector
            mvector(wells(producrow(iprod),1)) = presvalue(iprod);  
            %Null both row and column of global matrix:
            %Null column
            M(:,wells(producrow(iprod),1)) = 0;
            %Null row
            M(wells(producrow(iprod),1),:) = 0;
            %Put 1 in the position (i,i)
            M(wells(producrow(iprod),1),wells(producrow(iprod),1)) = 1;
        end  %End of FOR (producer)
    end  %End of IF (producer)
    
end  %End of IF
