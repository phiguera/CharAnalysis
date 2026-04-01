function s = charYLabel(transform)
% charYLabel   Return a CHAR y-axis label string for the given transform.
%
%   s = charYLabel(transform)
%
%   transform == 0  ->  'CHAR (pieces cm^-^2 yr^-^1)'
%   transform == 1  ->  'log CHAR (pieces cm^-^2 yr^-^1)'
%   transform == 2  ->  'ln CHAR (pieces cm^-^2 yr^-^1)'

switch transform
    case 1
        s = 'log CHAR (pieces cm^-^2 yr^-^1)';
    case 2
        s = 'ln CHAR (pieces cm^-^2 yr^-^1)';
    otherwise
        s = 'CHAR (pieces cm^-^2 yr^-^1)';
end

end
