function vessels = plotPts(filename)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    
    %importing
    pts_file = importPtsFile(filename);
    
    %checking
    if(size(pts_file,1)==0) 
        error('Empty pts file');
    end
    
    %reading vessels
    vessels = {}; %BC_in BC_out coord
    line = pts_file{1,1};
    line_id = 2;
    vessel_id = 1;
    while ~(strcmp(line,'END_LIST'))
        line = pts_file{line_id,1};
        if (strcmp(line,'BEGIN_ARC'))
            line_id = line_id +1; 
            vessels{vessel_id,1} = strcat(pts_file{line_id,1},pts_file{line_id,2});
            line_id = line_id +1;
            vessels{vessel_id,2} = strcat(pts_file{line_id,1},pts_file{line_id,2});
            coord = [];
            line_id = line_id +1;
            coord(1,:) = [pts_file{line_id,2} pts_file{line_id,3} pts_file{line_id,4}];
            end_id = line_id+1;
            line_id = line_id +2;
            while ~(strcmp(line,'END_ARC'))
                coord(end+1,:) = [pts_file{line_id,2} pts_file{line_id,3} pts_file{line_id,4}];
                line_id = line_id +1;
                line = pts_file{line_id,1};
            end % vessel cycle
            coord(end+1,:) = [pts_file{end_id,2} pts_file{end_id,3} pts_file{end_id,4}];
            vessels{vessel_id,3} = coord;
            vessel_id = vessel_id + 1;
        else
            line_id = line_id +1;   
        end % if
        
    end % while main cycle
    
    %% plotting
    figure
    hold on
    for i=1:size(vessels,1)
        coord = vessels{i,3};
        plot3(coord(:,1),coord(:,2),coord(:,3));
        mean_point_id = floor(1+size(coord,1)/2);
        %text(coord(mean_point_id,1),coord(mean_point_id,2),coord(mean_point_id,3),sprintf('Vessel %d',i));
        %in
        plot3(coord(1,1),coord(1,2),coord(1,3),'ro','MarkerFaceColor','r');
        text(coord(1,1),coord(1,2),coord(1,3),vessels{i,1});
        %out
        plot3(coord(end,1),coord(end,2),coord(end,3),'bo','MarkerFaceColor','b');
        text(coord(end,1),coord(end,2),coord(end,3),vessels{i,2});
    end
    hold off
    
end