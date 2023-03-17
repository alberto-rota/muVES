function [indexes,flip] = findDIR(vessels)
%findDIR gives indexes of the vessel with at least 1 DIR
indexes = [];
flip = [];
for v=1:size(vessels,1)
    if contains(vessels(v,1),'DIR') 
        indexes = [indexes v];
        flip = [flip 0];
    elseif  contains(vessels(v,2),'DIR')
        indexes = [indexes v];
        flip = [flip 1];
    end     
end

end

