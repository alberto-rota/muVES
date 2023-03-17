function string = addVesselToString(stringIN,vessel,flip,id)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
string = strcat(stringIN,'BEGIN_ARC\n');
if flip
    string = strcat(string,vessel{2});
    string = strcat(string,'\n');
    string = strcat(string,vessel{1});
    string = strcat(string,'\n');
else
    string = strcat(string,vessel{1});
    string = strcat(string,'\n');
    string = strcat(string,vessel{2});
    string = strcat(string,'\n');
end

points = vessel{3};
orderedVessel = orderVesselToPts(points,flip);
string = strcat(string,sprintf('%d\t%0.6f\t%0.6f\t%0.6f\tstart',id,orderedVessel(1,1),orderedVessel(1,2),orderedVessel(1,3)));
string = strcat(string,'\n');
string = strcat(string,sprintf('%d\t%0.6f\t%0.6f\t%0.6f\tend',id,orderedVessel(2,1),orderedVessel(2,2),orderedVessel(2,3)));
string = strcat(string,'\n');
for i=3:size(orderedVessel,1)
    string = strcat(string,sprintf('%d\t%0.6f\t%0.6f\t%0.6f\tpoint',id,orderedVessel(i,1),orderedVessel(i,2),orderedVessel(i,3)));
    string = strcat(string,'\n');
end
string = strcat(string,'END_ARC\n');

end

