function [IDJunction,flip] = findJunctions(vessels,IDVessel,flipIN)
    %findJunctions return id of the downstream vessels
    %   Detailed explanation goes here
    

    %starting vessel retrieval
    starting_vessel = vessels(IDVessel,:);
    
    % Finding id in the junction
    IDJunction = [];
    flip =[];
    points = starting_vessel{3};
    if flipIN
        end_point = points(1,:);
    else 
        end_point = points(end,:);
    end
    
    
    for v=1:size(vessels,1)
        if v~= IDVessel
            vessel = vessels(v,:);
            points_v = vessel{3};
            start_v = points_v(1,:);
            end_v = points_v(end,:);
            if isequal(end_point,start_v) 
                IDJunction = [IDJunction v];
                flip = [flip 0];
            elseif isequal(end_point,end_v)
                IDJunction = [IDJunction v];
                flip = [flip 1];
            end
        end
    end
end