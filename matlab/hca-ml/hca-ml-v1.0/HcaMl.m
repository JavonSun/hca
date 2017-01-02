% search for linear combination of a set of features by minimizing negative log likilihood 

% min_{\beta}'(sum_i(H_i' * {/phi}_i)^{-1} * H_i){\beta} + \lambda(sum_i(u_i + t_i))
% s.t. {\beta}'H'H{\beta} - N = 0
%      u, t >= 0
%      uXm' * w - um = 0
%      uXf' * w - uf = 0
% H = [X, -X, 1, [0/-1]_m, [0/-1]_f]
% \beta = [u; t; um; uf]
% w = u - v
% uXm = mean(X_m, 1)'  sample male mean
% uXf = mean(X_f, 1)' sample female mean
% X = x - z(z'z)^(-1)z'x

% Input:
%   x -- N * d, data
%   z -- N * p, covariants, set to [] if no covariates
%   gender -- N * 1, 1 means male, 2 means female
%   fams -- N_f * 1, each cell contains one family,
%           subjects are represented by the corresponding index in x
%   kinmat - N_f * 1, kinship matrix, family members alignment in each
%            matrix is consistent with that in corresponding family in fams
%   lambda -- a scalar

% Output:
%   t -- the quantitative score vector
%   w = u - t, the feature weight vector
%   r -- residual after removing from t the effect from covariates, if
%       there is no covariates, r is identical to t
%   v -- coefficient of covariates
%   funval -- the value of the objective function with returned linear model

% Javon, Sep. 3, 2013

function [t, w, r, v, funval] = HcaMl(x, z, gender, fams, kinmat, lambda)

N = size(x, 1);
d = size(x, 2);

if ~isempty(z)
    X = x - (z / (z' * z)) * z' * x;
else
    X = x;
end

bisex = false;
if length(unique(gender)) == 2
    bisex = true;
end

if bisex
    idc_male = zeros(N, 1); % indicator vector of male
    idc_female = zeros(N, 1); % indicator vector of female
    idc_male(gender == 1) = 1;
    idc_female(gender == 2) = 1;

    H = [X -X -idc_male -idc_female];
else
    H = [X -X (-1 * ones(N, 1))];
end

% formulate \sum_i(H_i'{\phi_i}^{-1}H_i)
mat1 = zeros(size(H, 2));
for i = 1:length(fams)
    mat1 = mat1 + H(fams{i}, :)' / (2 * kinmat{i}) * H(fams{i}, :);
end

mat2 = H' * H; % H'H

if bisex
    [r w um uf funval] = biSex(X, gender, N, d, mat1, mat2, lambda);
else
    [r w u funval] = uniSex(X, N, d, mat1, mat2, lambda);
end

if ~isempty(z)
    t = x * w;
    v = ((z' * z) \ z') * t;
else
    t = r;
    v = [];
end
end

function [t w um uf funval] = biSex(x, gender, N, d, mat1, mat2, lambda)

uxm = sum(x(gender == 1, :), 1) / length(find(gender == 1)); % sample mean of male
uxf = sum(x(gender == 2, :), 1) / length(find(gender == 2)); % sample mean of female
Aeq = [uxm -uxm -1 0; uxf -uxf 0 -1];
beq = [0 0]';

% options = optimset('GradObj','on', 'Hessian','user-supplied', 'GradConstr','on', 'MaxFunEvals', 1e04, 'MaxIter', 1e04);
options = optimset('GradObj','on', 'Algorithm', 'active-set', 'GradConstr','on', 'MaxFunEvals', 1e04, 'MaxIter', 1e04);
% Matlab 2014 complains: Hessian option set to 'user-supplied' but no Hessian function provided in options HessFcn
% nor in HessMult. So removed 'Hessian','user-supplied'. 3/25/2015

% initialization
u0 = ones(d, 1);
v0 = zeros(d, 1);
t0 = x * (u0 - v0);
um0 = mean(t0(gender == 1));
uf0 = mean(t0(gender == 2));
beta0 = [u0; v0; um0; uf0];
beta0 = beta0 * sqrt(N / (beta0' * mat2 * beta0)); % make beta0 be a feasible point

beta = fmincon(@objective, beta0, [], [], Aeq, beq, [zeros(2 * d, 1); -Inf; -Inf], [], @constrain, options);
w = beta(1:d) - beta(d+1:2*d);
um = beta(end - 1);
uf = beta(end);
t = x * w;
funval = objective(beta);

    function [o g h] = objective(p)
        o = p' * mat1 * p + lambda * sum(p(1:end-2));
        g = 2 * mat1 * p + lambda * [ones(2 * d, 1); 0; 0];
        h = mat1;
    end

    function [c ceq gc gceq] = constrain(p)
        c = [];
        ceq = p' * mat2 * p - N;
        gc = [];
        gceq = 2 * mat2 * p;
    end

end

function [t w u funval] = uniSex(x, N, d, mat1, mat2, lambda)

ux = sum(x, 1) / N; % sample mean 
Aeq = [ux -ux -1];
beq = 0;

% options = optimset('GradObj','on', 'Hessian','user-supplied', 'GradConstr','on', 'MaxFunEvals', 1e04, 'MaxIter', 1e04);
options = optimset('GradObj','on', 'GradConstr','on', 'MaxFunEvals', 1e04, 'MaxIter', 1e04);
% Matlab 2014 complains: Hessian option set to 'user-supplied' but no Hessian function provided in options HessFcn
% nor in HessMult. So removed 'Hessian','user-supplied'. 3/25/2015

% initialization
u0 = ones(d, 1);
v0 = zeros(d, 1);
t0 = x * (u0 - v0);
ut0 = mean(t0);
beta0 = [u0; v0; ut0];
beta0 = beta0 * sqrt(N / (beta0' * mat2 * beta0)); % make beta0 be a feasible point

beta = fmincon(@objective, beta0, [], [], Aeq, beq, [zeros(2 * d, 1); -Inf], [], @constrain, options);
w = beta(1:d) - beta(d+1:2*d);
u = beta(end);
t = x * w;
funval = objective(beta);

    function [o g h] = objective(p)
        o = p' * mat1 * p + lambda * sum(p(1:end-1));
        g = 2 * mat1 * p + lambda * [ones(2 * d, 1); 0];
        h = mat1;
    end

    function [c ceq gc gceq] = constrain(p)
        c = [];
        ceq = p' * mat2 * p - N;
        gc = [];
        gceq = 2 * mat2 * p;
    end

end









