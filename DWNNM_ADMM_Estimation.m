function [ Y_hat, W_hat] = MCWNNM_ADMM2_Estimation( NL_mat, Sigma_arr, CurPat,WeiPat, Par )

Y_hat = zeros(size(CurPat));
W_hat    = zeros(size(CurPat));
for  i      =  1 : length(Par.SelfIndex) % For each keypatch group
    Y    =   CurPat(:, NL_mat(1:Par.nlsp,i)); % Non-local similar patches to the keypatch
    W    =   WeiPat(:, NL_mat(1:Par.nlsp,i));  %add 对于关键块的所有近邻块，按照索引取出对应的权重块
    mY  =   repmat(mean( Y, 2 ),1,Par.nlsp);
    Y    =   Y-mY;
    X 	=   DWNNM_ADMM( Y,W, Sigma_arr(:, Par.SelfIndex(i)), Par); % DWNNM Estimation    取一列噪声强度Sigma_arr(:, Par.SelfIndex(i))
    Y_hat(:,NL_mat(1:Par.nlsp,i))  = Y_hat(:,NL_mat(1:Par.nlsp,i))+X+mY;
    W_hat(:,NL_mat(1:Par.nlsp,i))     = W_hat(:,NL_mat(1:Par.nlsp,i))+ones(Par.ps2ch, Par.nlsp);
end
end

