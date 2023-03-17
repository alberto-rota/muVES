function [good,problemsID] = checkPts(filename)
    %check the goodness of the network
    %   filename: name of .pts to be used

     %importing
    pts_file = importPtsFile(filename);
    
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
    
    problemsID = [];
    for v=1:size(vessels,1)
        if contains(vessels{v,1},'BC INT')
            points_v = vessels{v,3};
            start_v = points_v(1,:);
            %disp('Searching for INT')
            if ~checkINT(vessels,start_v)
                problemsID = [problemsID v];
                str = sprintf('Problem found in inlet of vessel %d',v);
                disp(str)
            end
        end
        if contains(vessels{v,2},'BC INT')
            points_v = vessels{v,3};
            end_v = points_v(end,:);
            %disp('Searching for INT')
            if ~checkINT(vessels,end_v)
                problemsID = [problemsID v];
                str = sprintf('Problem found in inlet of vessel %d',v);
                disp(str)
            end
        end
        
    end

    good = isempty(problemsID);

end