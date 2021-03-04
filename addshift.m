function B = addshift(A, p, dim)
%ADDSHIFT(A, p, dim) shifts matrix 'A' of 'p' places in the direction specified
%   by 'dim'. If 'A' is a NxM matrix and:
%   - dim == 1: Produces a N'xM matrix where N' = N+p;
%   - dim == 2: Produces a NxM' matrix where M' = M+p;
%   For p>0, the matrix is shifted right (or down), adding zeros to the left
%   (or up). For p<0, the matrix is shifted left (or up), adding zeros to the
%   right (or down);
%
%   Example: A = [2 3 4; 5 6 7]
%   B = ADDSHIFT(A,3,1);   
%   B =
%       0     0     0
%       0     0     0
%       0     0     0
%       2     3     4
%       5     6     7

[r, c] = size(A);
if dim==1   
    z = zeros(abs(p),c);
    if p>0
        B = [z; A];
    else
        B = [A; z];
    end
elseif dim==2   
    z = zeros(r,abs(p));
    if p>0
        B = [z A];
    else
        B = [A z];
    end
else
    error('Specify a valid dimension');
end
end
