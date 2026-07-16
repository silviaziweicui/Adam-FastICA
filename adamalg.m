%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Moment Estimation algorithm to iteratively optimize a set of weights (MU
% filter) to maximize the sparseness of the source (MU pulse train)
% Input: 
%   w = initial weigths
%   X = whitened signal
%   B = separation matrix of MU filters
%   maxiter = maximal number of iteration before convergence
%   contrastfunc = contrast function


% Output:
%   w = weigths (MU filter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w,k] = adamalg(w, X, B , maxiter,  contrastfunc)

k = 1;
% delta(k) = 1;
delta = ones(1,maxiter);
TOL = 0.0001;  % tolerance between two iterations 0.0001
BBT = B * B';  % Precompute B*B' for orthogonalization
b1 = 0.9; 
b2 = 0.99;
LR = 0.001;  

switch contrastfunc
    case 'skew'
        gp = @(x) 2*x;  
        g = @(x) x.^2;
    case 'kurtosis'
        gp = @(x) 3*x.^2;
        g = @(x) x.^3;
    case 'logcosh'
        gp = @(x) tanh(x);
        g = @(x) log(cosh(x));
    case 'sparse_kurtosis'  
        gp = @(x) 4 * abs(x).^3 .* sign(x); 
        g = @(x) x.^4;    
    case 'exp_neg' 
        gp = @(x) -x .* exp(-x.^2 / 2);
        g = @(x) -exp(-x.^2 / 2);
end

% Initialize Adam parameters
if k ==1
    mt = zeros(size(w));
    vt = zeros(size(w));
end

while delta(k) > TOL && k < maxiter
    % Update weights
    wlast = w; % Save last weights

    % Contrast function
    wTX = w' * X;                 
    term1 = mean(X .* g(wTX), 2); 
    beta = mean(wTX .* g(wTX));   
    gt = term1 - beta * w;  
    gt = gt / (norm(gt) + 1e-10);  

    mt = b1 * mt + (1 - b1) * gt;  
    vt = b2 * vt + (1 - b2) * (gt.^2);
    
    mtt = mt / (1 - b1^k);
    vtt = vt / (1 - b2^k);      

    % Update weights
    dc = mtt./(sqrt(vtt) + 1e-8);  
    w = w - LR * dc;

    % 3b: Orthogonalization
    w = w - BBT * w;

    % 3c: Normalization
    w = w / norm(w);

    % Update convergence criteria
    k = k + 1;
    delta(k) = abs(w' * wlast - 1);    
end