function  [Z] =  DWNNM_ADMM( Y,W_local, NSig, Par )
% This routine solves the following weighted nuclear norm optimization problem with column weights,
%
% min_{X, Z} ||W(Y-X)||_F^2 + ||Z||_w,*  s.t.  X = Z
%
% Inputs:
%        Y      -- 3p^2 x M dimensional noisy matrix, D is the data dimension, and N is the number of image patches.
%         W      --3p^2 x M  dimensional  of local weights
%        NSig -- 3p^2 x 1 dimensional vector of noise weights
%        Par   -- structure of parameters
% Output:
%        Z      -- 3p^2 x M dimensional denoised matrix

% tol = 1e-8;
if ~isfield(Par, 'maxIter')
    Par.maxIter = 10;
end
if ~isfield(Par, 'rho')
    Par.rho = 1;
end
if ~isfield(Par, 'mu')
    Par.mu = 1;
end
if ~isfield(Par, 'display')
    Par.display = true;
end
%add matrix size
[matrix_row,matrix_col]=size(Y);
% Initializing optimization variables
% Intialize the weight matrix W
NSig=repmat(NSig,1,matrix_col);%resize to the same size with W_local
W_n = 1 ./ (NSig+eps);
% Initializing optimization variables
X = zeros(size(Y));
Z = zeros(size(Y));
A = zeros(size(Y));
%% Start main loop
iter = 0;
PatNum       = size(Y,2);
TempC  = Par.Constant * sqrt(PatNum);
while iter < Par.maxIter
    iter = iter + 1;
    
    % update X, fix Z and A
    % min_{X} ||W [W_n* (Y - X)]||_F^2 + 0.5 * rho * ||X - Z + 1/rho * A||_F^2   %W represent hadamard product
    %transform the ordinary matrix to diag matrix
    Y_diag=diag(Y(:));
    Z_diag=diag(Z(:));
    A_diag=diag(A(:));
    W_n_diag=diag(W_n(:));
    W_diag=diag(W_local(:));
    W_all=W_diag*W_n_diag;   %final weight matrix
    I_diag=eye(matrix_col*matrix_row,matrix_col*matrix_row);%unit matrix
    X_diag = (W_all.^2 + 0.5 * Par.rho*I_diag) \ (W_all.^2 * Y_diag + 0.5 * Par.rho * Z_diag - 0.5 * A_diag); 
    %transform the diag matrix to original size
    X=reshape(diag(X_diag),matrix_row,matrix_col);  


    % update Z, fix X and A
    % min_{Z} ||Z||_*,w + 0.5 * rho * ||Z - (X + 1/rho * A)||_F^2
    Temp = X + A/Par.rho;
    [U, SigmaTemp, V] =   svd(full(Temp), 'econ');
    [SigmaZ, svp] = ClosedWNNM(diag(SigmaTemp), 2/Par.rho*TempC, eps);
    Z =  U(:, 1:svp) * diag(SigmaZ) * V(:, 1:svp)';
    %     % check the convergence conditions
    %     stopC = max(max(abs(X - Z)));
    %     if Par.display && (iter==1 || mod(iter,10)==0 || stopC<tol)
    %         disp(['iter ' num2str(iter) ',mu=' num2str(Par.mu,'%2.1e') ...
    %             ',rank=' num2str(rank(Z,1e-4*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
    %     end
    %     if stopC < tol
    %         break;
    %     else
    % update the multiplier A, fix Z and X
    A = A + Par.rho * (X - Z);
    Par.rho = min(1e4, Par.mu * Par.rho);
    %     end
end
return;
