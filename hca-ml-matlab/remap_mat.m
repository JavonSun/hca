% update kinship matrix with new index (index of subjects in another dataset),
% individuals with no new index will be removed from the returned family
% and kinship matrix

% Input:
%   idx_map -- a vector that contains new index of individuals
%   f1 -- cell array, N_f * 1, each cell contains members of one family,
%         individuals are represented by old index
%   m1 -- cell array, N_f * k, k matrixes of each family in f1

% Output:
%   f2 -- families in which individuals are indexed with new index
%   m2 -- maped matrixes

% Javon, Jan. 26, 2013

function [f2, m2] = remap_mat(idx_map, f1, m1)

assert(length(f1) == size(m1, 1));

ief = []; % index of families became empty when shrinking
f2 = cell(length(f1), 1);
m2 = m1;

for i = 1:length(f1)
    fm = idx_map(f1{i}); % remap
    f2{i} = fm;

    % check if all present
    iNa = isnan(fm); % index of members do not present in the data within the family
    if (isempty(find(~iNa, 1)))
        ief = [ief; i]; % empty family
        continue;
    end
    if (isempty(find(iNa, 1)))
       % all subjects in this family present in the data
       continue;
    end
    % shrink 
    fm(iNa) = [];
    f2{i} = fm;
    for j = 1:size(m1, 2)
       m = m1{i, j};
       m(iNa, :) = [];
       m(:, iNa) = [];
       m2{i, j} = m;
    end
end

if ~isempty(ief)
    f2(ief) = []; % remove empty families
    m2(ief, :) = [];
end

end