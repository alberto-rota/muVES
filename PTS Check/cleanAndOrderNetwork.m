function [Network,cleanedID] = cleanAndOrderNetwork(filename, radius)
    %cleanAndOrderNetwork starting from a .pts skeleton, this function
    %order the network as requeswted by the 3d1dcode after a cleaning based
    %on BC (no MIX - MIX condition)
    %   filename: name of .pts to be used
    %   radius: file name of the radius
    
    %% Reading .pts
    
    %checking
    switch nargin
        case 2
            checkRadius = true;
        case 1 
            radius ='empty';
            checkRadius = false;
        case 0 
            error('.pts file not selected')
    end
    
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
    
    %% Ordering and printing pts
    
    %string to be printed
    toBePrinted = "BEGIN_LIST\n";
    checkedVessels = [];
    
    % find a DIR-INT to start
    [DIR_id,flip] = findDIR(vessels);
        
    % start the while for unconnected network
    while length(checkedVessels) < size(vessels,1)
         
        % retrieving the first ID DIR not checked
        current_end =[];
        for i=1:length(DIR_id)
            if ~ismember(DIR_id(i),checkedVessels)
                current_end = DIR_id(i);
                current_flip = flip(i);
                break;
            end
        end
        % check for no DIR left
        if isempty(current_end)
            w_msg = sprintf('Unconnected region inside are left (processed %d out of %d)', length(checkedVessels), size(vessels,1));
            warning(w_msg)
            break;
        else
            %writing first vessel
            checkedVessels = [checkedVessels current_end];
            toBePrinted = addVesselToString(toBePrinted,vessels(current_end,:),current_flip,length(checkedVessels));
        end
        
        %cycle over a connected network
        while ~isempty(current_end)
        
            % - find junctions (of one end)
            [downstreamV,flip_V] = findJunctions(vessels,current_end(1),current_flip(1));
            
            for d=1:length(downstreamV)
                if ~ismember(downstreamV(d),checkedVessels)
                    % - update current ends
                    % - add current end of the network
                    current_end = [current_end downstreamV(d)];
                    current_flip = [current_flip flip_V(d)];
                    % - add added vessels to list
                    checkedVessels = [checkedVessels downstreamV(d)];
                    % - print vessels (if not printed)
                    toBePrinted = addVesselToString(toBePrinted,vessels(downstreamV(d),:),flip_V(d),length(checkedVessels));   
                end
            end
            
            % - update current ends
            %delete old
            current_end = current_end(2:end);
            current_flip = current_flip(2:end);
        end % while for connected network
    
    % if vessels are left in the initial list, restart while
    end %while

    toBePrinted = strcat(toBePrinted,'END_LIST');
    
    % Output
    Network = vessels;
    cleanedID = checkedVessels;
    
    %Write to .pts
    cleanedname= strcat(filename,'_cleaned.pts');
    fileID = fopen(cleanedname,'wt');
    fprintf(fileID,toBePrinted);
    fclose(fileID);
    
    %Handling radius
    if checkRadius
        radius_vector = importPtsRadius(radius);
        radius_String = 'BEGIN_LIST\n';
        for v=1:length(checkedVessels)
            radius_String = strcat(radius_String,sprintf('%s',radius_vector{v}));
            radius_String = strcat(radius_String,'\n');
        end
        radius_String = strcat(radius_String,'END_LIST');

        cleanedname= strcat(radius,'_cleaned.pts');
        fileID = fopen(cleanedname,'wt');
        fprintf(fileID,radius_String);
        fclose(fileID);
    end
end