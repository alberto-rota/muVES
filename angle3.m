function ang = angle3(v1, v2)
%ANGLE3(v1, v2) calculates the angle between the 3D vectors 'v1' and 'v2',
%   of dimension 3x1 or 1x3 . The angle is returned in radians.
% 
%    Example
%     v1 = [1 0 0];
%     v2 = [0 1 0];
%     a = angle3(v1,v2)
% 
%     a =  1.5708

    if all(size(v1) > 1) || all(size(v2) > 1) 
        error('The inputs are not both 3x1 or 1x3 vectors');
    end
    if isrow(v1)
        v1 = v1';
    end
    if isrow(v2)
        v2 = v2';
    end
    ang = atan2(norm(cross(v1,v2)),dot(v1,v2));
end