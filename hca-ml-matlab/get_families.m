% parse a pedigree matrix, extract all members of each family and organize
% them in parent -> progeny order

% input:
%   pedigree - each row represents one person, three column: column 1 -
%              individual id, column 2 - father id, column 3- monther id

% Javon Sun, 2/16/2011
function fams = get_families(pedigree)
    ind_ids = pedigree(:, 1);
    fam_mem = cell(length(ind_ids), 1);    
    % parents = zeros(length(ind_ids), 0); % 1 means the individual is some one's parent, otherwise not
    par_prog = zeros(length(ind_ids), 1); % 1 means the individual is some one's parent, 0 means progeny
    par_fam = zeros(length(ind_ids), 1); % family number of each individual in parents vector
    % progeny = ones(length(ind_ids), 0);  % 1 means the individual is some one's progeny, otherwise not
    
    % initialize, identify those missing both parents as parent
    par_prog(isnan(pedigree(:, 2)) & isnan(pedigree(:, 3))) = 1;
    par_fam(:) = 1 : length(ind_ids);
    par_fam(par_prog == 0) = 0;
    for i = (find(par_prog ==  1))'
        fam_mem{i} = i;
    end;
    
    nProg_unknown_par = length(find(par_prog == 0)); % number of progenies with unknown parents
    while nProg_unknown_par > 0
        prog_unknown_par = find(par_prog == 0);
        for i = prog_unknown_par'
            % scan each progeny with unknown parents, and try to unite
            % families with same progeny
            if (~isnan(pedigree(i, 2)) && par_prog(ind_ids == pedigree(i, 2)) == 0) ...
                    || (~isnan(pedigree(i, 3)) && par_prog(ind_ids == pedigree(i, 3)) == 0)
                % father or mother still unknown
                continue;
            end;
            if ~isnan(pedigree(i, 2)) && ~isnan(pedigree(i, 3))
                dad_fam = par_fam(find(ind_ids == pedigree(i, 2), 1));
                mom_fam = par_fam(find(ind_ids == pedigree(i, 3), 1));
                if dad_fam ~= mom_fam
                    % combine father's family and mother's family using the
                    % smaller family id
                    com_fam = min(dad_fam, mom_fam);
                    if com_fam == dad_fam
                        % mother's family merges into father's family
                        fam_mem{com_fam} = [fam_mem{com_fam} fam_mem{mom_fam}];
                        par_fam(fam_mem{mom_fam}) = com_fam;
                        fam_mem{mom_fam} = [];
                    else
                        % father's family merges into mother's family
                        fam_mem{com_fam} = [fam_mem{com_fam} fam_mem{dad_fam}];
                        par_fam(fam_mem{dad_fam}) = com_fam;
                        fam_mem{dad_fam} = [];
                    end;
                else
                    com_fam = dad_fam;
                end;
            elseif ~isnan(pedigree(i, 2))
                % mother missing
                com_fam = par_fam(find(ind_ids == pedigree(i, 2), 1));
            else
                % father missing
                com_fam = par_fam(find(ind_ids == pedigree(i, 3), 1));
            end;
            % add i to the family
            fam_mem{com_fam} = [fam_mem{com_fam} i];
            par_prog(i) = 1;
            par_fam(i) = com_fam;
        end;
        
        if nProg_unknown_par == length(find(par_prog == 0))
            warning('Incomplete pedigree file: some parents may not be in the file as an individual');
            break;
        end;
        nProg_unknown_par = length(find(par_prog == 0));
    end; 
    
    % prepare the output
    fams = cell(0, 1);
    cnt = 0;
    for i = 1 : length(ind_ids)
        if ~isempty(fam_mem{i})
            cnt = cnt + 1;
            fams{cnt} = fam_mem{i};
        end;
    end;