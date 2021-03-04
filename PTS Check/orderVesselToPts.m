function orderedVessel = orderVesselToPts(points,flip)
%orderVesselToPts order vessels for pts (start, end, points). If flip
%enabled, it reverse the vessel
%   Detailed explanation goes here

if flip
    orderedVessel(1,:) = points(end,:);
    orderedVessel(2,:) = points(1,:);
    for i=2:size(points,1)-1
        orderedVessel(i+1,:)=  points(end-i+1,:);
    end
else
    orderedVessel(1,:) = points(1,:);
    orderedVessel(2,:) = points(end,:);
    for i=2:size(points,1)-1
        orderedVessel(i+1,:)=  points(i,:);
    end
end
end

