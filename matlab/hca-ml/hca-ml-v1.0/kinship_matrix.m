% compute the kinship matrix based on pedigree

% input:
%   pedigree - each row represents one person, three column: column 1 -
%              individual id, column 2 - father id, column 3- monther id
%   fams - a cell array of family, each cell contains all the members of one familly, and 
%          members need to be organized in parent -> progeny order, which
%          can be obtained by calling the function get_families()

% output:
%   kmat - a cell array of kinship matrixes, one per family. The matrix is ordered 
%       in the same as in input fams.

% Javon Sun, 2/16/2011

function kmat = kinship_matrix(pedigree, fams)
    kmat = cell(length(fams), 1);
    for f = 1 : length(fams)
        % compute kinship matrix family by family
        kmat{f} = zeros(length(fams{f}));
        kmat{f}(:) = nan;
        for i = 1 : length(fams{f})           
            ii = fams{f}(i);
            if isnan(pedigree(ii, 2))
                % father missing
                fi = nan;
            else
                fi = find(pedigree(:, 1) == pedigree(fams{f}(i), 2)); % index of father in pedigree marix
                fii = find(fams{f} == fi);  % index of father in the family
            end;
            if isnan(pedigree(ii, 3))
                % mother missing
                mi = nan;
            else
                mi = find(pedigree(:, 1) == pedigree(fams{f}(i), 3)); % index of mother in pedigree marix
                mii = find(fams{f} == mi);  % index of mother in the family
            end;
            
            % kinship of i with himself or herself
            if isnan(fi) || isnan(mi)
                % father or mother missing
                kmat{f}(i, i) = 1 / 2;
            else
                % display(sprintf('f = %d, i = %d', f, i));
                kmat{f}(i, i) = (1 / 2) * (1 + kmat{f}(fii, mii));
            end;
            
            % kinship of i with all foreruners in the family if any
            if i > 1
                for j = i - 1 : -1 : 1
                    if isnan(fi)
                        % father missing
                        f_j = 0;
                    else
                        f_j = kmat{f}(j, fii);
                    end;
                    if isnan(mi)
                        % mother missing
                        m_j = 0;
                    else
                        m_j = kmat{f}(j, mii);
                    end;
                    kmat{f}(i, j) = (1 / 2) * (f_j + m_j);
                    kmat{f}(j, i) = kmat{f}(i, j);
                end;
            end;
        end;
    end;