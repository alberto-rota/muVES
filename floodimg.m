function img = floodimg(img, basepixel)
%FLOODIMAGE highlights the region of white pixels to which 'basepixel'
%	belongs. 'basepixel' is a 1x2 or 2x1 vector containing the coordinates of
%   the basepoint.
%   'img' should be a logical matrix of black and white pixels, and it will
%   be converted to this format if this is not true (with the imbinarize
%   function).
%   FLOODIMAGE returns the image in a 'double' format with a value of 2 on
%   the pixels in the found region
%
%   Example: 
%     img =
% 
%          1     1     1     0     0     0
%          1     1     1     0     0     0
%          1     1     1     0     0     0
%          0     0     0     0     1     1
%          1     1     0     0     1     1
%          1     1     0     0     1     1
% 
%     region = floodimg(img, [2;2])
%     region =
% 
%          2     2     2     0     0     0
%          2     2     2     0     0     0
%          2     2     2     0     0     0
%          0     0     0     0     1     1
%          1     1     0     0     1     1
%          1     1     0     0     1     1

    dimlims = [[1;1] size(img)'];
    if isrow(basepixel)
        basepixel = basepixel';
    end
    if any(basepixel<dimlims(:,1)) |  any(basepixel>dimlims(:,2))
        return;
    end
    if img(basepixel(1), basepixel(2)) == 0 | img(basepixel(1), basepixel(2)) == 2
        return;
    else
        img(basepixel(1), basepixel(2)) = 2;
        img = floodimg(img,[basepixel(1)+1, basepixel(2)]);
        img = floodimg(img,[basepixel(1)-1, basepixel(2)]);
        img = floodimg(img,[basepixel(1), basepixel(2)+1]);
        img = floodimg(img,[basepixel(1), basepixel(2)-1]);
        img = floodimg(img,[basepixel(1)+1, basepixel(2)+1]);
        img = floodimg(img,[basepixel(1)-1, basepixel(2)-1]);
        img = floodimg(img,[basepixel(1)-1, basepixel(2)+1]);
        img = floodimg(img,[basepixel(1)+1, basepixel(2)-1]);
        
    end
end