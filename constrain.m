function ac = constrain(a,range)
% CONSTRAIN(a, range) Holds the value of 'a' inside the interval defined by
%   'range'. If 'a' is outside the range specified, its value will be set
%   either to the maximim or minimum value in 'range'
%
%    Example
%     a = 10;
%     ac = constrain(a, [4 8])
%     ac =  8
% 
%     ac = constrain(a, [4 12])
%     ac =  10
if isnumeric(a) && isnumeric(range) && isscalar(a)
    if (a>max(range))
        ac = max(range);
    elseif(a<min(range))
        ac = min(range);
    else 
        ac = a;
    end
else
    error('Invalid input');
end

