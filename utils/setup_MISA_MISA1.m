function data1 = setup_MISA_MISA1(M,X,W0,dist_type)

nComp = size(W0{M(1)},1);

S = cell(1,max(M));   % Cell array: each cell contains a matrix K x C(m).

% Each k-th row has 0's and 1's to indicate what sources go within the k-th subspace in dataset m
for mm=M
%    S{mm} = sparse(1:nComp,1:nComp,1,nComp,nComp,nComp);
    S{mm} = [eye(nComp-1) [zeros(1, nComp-2) 1]'];
              
end

% params
K=size(S{M(1)},1);   % Number of subspaces

if strcmpi(dist_type,'gaussian')
    % Set Kotz parameters to multivariate Gaussian with variance = 3.2899
    eta = ones(K,1);
    beta = ones(K,1);
    lambda = (1.5/(pi^2))*ones(K,1);
elseif strcmpi(dist_type,'kotz')
    % Set Kotz parameters to multivariate Kotz with variance = 3.2899
    eta = ones(K,1);
    beta = 0.546172370649572*ones(K,1);
    lambda = 0.896584116042220*ones(K,1);
end

% Use relative gradient
gradtype = 'relative';

% Enable scale control
sc = 1;

% Turn off preprocessing (still removes the mean of the data)
preX = false;

ut = utils;
w0 = ut.stackW(W0(M)); % vectorize unmixing matrix for compatibility with Matlab's optimization toolbox

% % Use normalized MSE
 REtype = 'NMSE';

% % Use the transpose of W as the reconstruction approach
 REapproach = 'WT';%'WT'; % 'PINV' for pseudoinverse or W

% % Tolerance level (0 means least error possible)
 RElambda = 5e-5;

% % Other parameters required by the @MISAKRE API but not used
 REref = {};
 REreftype = 'linearreg';
 REreflambda = {.9};
 rC = {[],[]};

%% Initialize MISA object

data1 = MISAK(w0, M, S, X, ...
                beta, eta, lambda, ...
                gradtype, sc, preX);
%                 REtype, REapproach, RElambda, ...
%                 REref, REreftype, REreflambda, rC);
 
% [~, old_M, ~, ~, ~, old_ref] = data1.update(data1.S,1:8,beta,lambda,eta,data1.reference);
