function [ Y_hat, W_hat] = DWNNM_Estimation( NL_mat, Sigma_arr, CurPat,WeiPat, Par )

Y_hat = zeros(size(CurPat));
W_hat    = zeros(size(CurPat));
for  i      =  1 : length(Par.SelfIndex) % For each keypatch group
    Y    =   CurPat(:, NL_mat(1:Par.nlsp,i)); % Non-local similar patches to the key patch
    W    =   WeiPat(:, NL_mat(1:Par.nlsp,i));  %add Non-local weight patches to the key patch
    mY  =   repmat(mean( Y, 2 ),1,Par.nlsp);
    Y    =   Y-mY;
    X 	=   DWNNM_ADMM( Y,W, Sigma_arr(:, Par.SelfIndex(i)), Par); % DWNNM Estimation   Sigma_arr(:, Par.SelfIndex(i))  one clomn of Sigma
    Y_hat(:,NL_mat(1:Par.nlsp,i))  = Y_hat(:,NL_mat(1:Par.nlsp,i))+X+mY;
    W_hat(:,NL_mat(1:Par.nlsp,i))     = W_hat(:,NL_mat(1:Par.nlsp,i))+ones(Par.ps2, Par.nlsp);
end
end

