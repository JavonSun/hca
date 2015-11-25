% map subjects in s1 to those in s2, if present return the corresponding
% index in s2, otherwise NAN

function map = match(s1, s2)

[c ia ib] = intersect(s1, s2); 
map = zeros(size(s1, 1), 1); 
map(:) = nan;
map(ia) = ib;

end