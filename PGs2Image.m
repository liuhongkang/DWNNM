function      im_out = PGs2Image(X, W, par)
% Reconstruction
im_out = zeros(par.h, par.w);
im_wei = zeros(par.h, par.w);
r = 1:1:par.maxr;
c = 1:1:par.maxc;
k = 0;

for i = 1:1:par.ps
    for j = 1:1:par.ps
        k = k+1;
        im_out(r-1+i, c-1+j)  =  im_out(r-1+i, c-1+j) + reshape( X(k,:)', [par.maxr par.maxc] );
        im_wei(r-1+i, c-1+j)  =  im_wei(r-1+i, c-1+j) + reshape( W(k,:)', [par.maxr par.maxc] );
    end
end
im_out  =  im_out ./ im_wei;