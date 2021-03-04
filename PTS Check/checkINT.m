function good = checkINT(vessels,point)

    count = 0;
   for v=1:size(vessels,1)
            vessel = vessels(v,:);
            points_v = vessel{3};
            start_v = points_v(1,:);
            end_v = points_v(end,:);
            if isequal(point,start_v) && contains(vessel{1},'BC INT')
                count = count +1;
                %disp('found')
            elseif isequal(point,end_v) && contains(vessel{2},'BC INT')
                count = count +1;
                %disp('found')
            end
   end
    
   good = (count > 1);
   
   if good == 2
       str = sprintf('Trivial junction at id %d',v);
       disp(str)
   end

end