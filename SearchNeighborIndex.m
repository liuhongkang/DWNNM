function  par  =  SearchNeighborIndex(par)
% This Function Precompute the all the patch indexes in the Searching window
% -NeighborIndex is the array of neighbor patch indexes for each keypatch
% -NumIndex is array of the effective neighbor patch numbers for each keypatch
% -SelfIndex is the index of keypatches in the total patch index array
par.maxr = par.h - par.ps + 1;
par.maxc = par.w - par.ps + 1;%%?????¾å???????? step ä¸? 1 ?¨é?¨å????
r          =  1:par.step:par.maxr;
par.r          =  [r r(end) + 1:par.maxr];
c          =  1:par.step:par.maxc;
par.c          =  [c c(end) + 1:par.maxc];%% ?»å??è¿?ç¨?ä¸?ä»????¨å?³é???¾å????ï¼???ä¸??????????????????¾å????
par.lenr = length(par.r);
par.lenc = length(par.c);
par.ps2 = par.ps^2;  %?¾å????ä¸???ç´?ä¸???
% Total number of patches in the test image
par.maxrc = par.maxr * par.maxc;%??????ï¼?ä¸?????è¦?ä¹?ä»? 3
% Total number of seed patches being processed   ?³é???¾å????ä¸???
par.lenrc = par.lenr * par.lenc;
% index of each patch in image
par.Index     =   (1:par.maxrc);
par.Index    =   reshape(par.Index, par.maxr, par.maxc);  
par.NeighborIndex    =   int32(zeros(4 * par.win^2, par.lenrc));   
par.NumIndex        =   int32(zeros(1, par.lenrc));
par.SelfIndex   =   int32(zeros(1, par.lenrc));

for  i  =  1 : par.lenr
    for  j  =  1 : par.lenc
        row = par.r(i);
        col = par.c(j);
        off = (col-1) * par.maxr + row;   %?³é?????¨æ?????¾å????ä¸???ç´¢å?
        off1 = (j-1) * par.lenr + i;    %?³é????ç´¢å?
        
        % the range indexes of the window for searching the similar patches
        % ?³é?????¨æ?????¾å????ä¸?????ç´¢è????
        rmin    =   max( row - par.win, 1 );
        rmax    =   min( row + par.win, par.maxr );
        cmin    =   max( col - par.win, 1 );
        cmax    =   min( col + par.win, par.maxc );
        
        idx     =   par.Index(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        
        par.NumIndex(off1)  =  length(idx);
        par.NeighborIndex(1:par.NumIndex(off1), off1)  =  idx;
        par.SelfIndex(off1) = off;
    end
end
