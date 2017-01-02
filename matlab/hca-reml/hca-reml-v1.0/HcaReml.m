function [t, w, funval] = HcaReml(X, C, G, lambda)
% search linear combination of variables in X that leads to a quantitative
% trait with high heritability estimates using correlation matrix G

% Input: 
%     X - N * d matrix , data on d phenotypic features for N subjectss
%     C - N * c matrix, data on c covariates
%     G - N * N matrix, genetic correlation matrx
%     lambda - the hyper-parameter in the model

% Output:
%     t - a vector with size of N, values for the N subjects on the learned
%         trait
%     w - a vector with size of d, coefficients of learned linear
%         linear combination
%     funval - value fo the objective function

% Javon
% 4/29/2015

N = size(X, 1); % number of subjects
d = size(X, 2); % number of phenotypic features

%  calculate P
if isempty(C)
    C = ones(N, 1);
else
    C = [ones(N, 1) C];
end
tmp = C' / G * C;
tmp = C / tmp * C';
J = eye(N) - tmp / G;
P = G \ J;

% calculate Q 
Q = J' * J;

% compose H
H = [X, -X];

% calculate the two intermedia matrix
mat1 = H' * P * H / N; % normalize loss, so per sample
mat2 = H' * Q * H / N;

options = optimset('Algorithm', 'active-set', 'GradObj','on', 'Hessian','user-supplied', 'GradConstr','on', 'MaxFunEvals', 1e08, 'MaxIter', 1e04);

% initialization
u0 = ones(d, 1);
% u0 = rand(d, 1); % initialize with random number
v0 = zeros(d, 1);
gamma0 = [u0; v0];
gamma0 = gamma0 / sqrt((gamma0' * mat2 * gamma0)); % make gamma0 be a feasible point

% gamma = fmincon(@objective, gamma0, [], [], [], [], zeros(2 * d, 1), [], @constrain, options);
gamma = fmincon(@objective, gamma0, [], [], [], [], zeros(2 * d, 1), [], @constrain, options);
w = gamma(1:d) - gamma(d+1:2*d);
t = X * w;

[funval] = objective(gamma);

% inner function, calculate value, gradient and hessian of the objective
% function
function [o, g, h] = objective(gamma)
    o = gamma' * mat1 * gamma + lambda * sum(gamma) / d;
    g = 2 * mat1 * gamma + lambda * ones(2 * d, 1) / d;
    h = 2 * mat1;
end

% inner function, calculate value, gradient of constraints
function [c, ceq, gc, gceq] = constrain(gamma)
    c = [];
    ceq = gamma' * mat2 * gamma - 1;
    gc = [];
    gceq = 2 * mat2 * gamma;
end

end

