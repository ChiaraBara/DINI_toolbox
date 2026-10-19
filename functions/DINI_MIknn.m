function IXY = DINI_MIknn(B,iy,ix,k)
    
    metric='chebychev';
    
    M_Y = B(:,iy);
    M_X = B(:,ix);
     
    M_YX = [M_Y M_X];
    
    N = size(M_Y,1); m = size(M_X,2);
      
    %%% neighbor search in the lower density space 
    [~, distances] =  knnsearch(M_YX,M_YX,'K',k+1,'Distance',metric);
    dd = distances(:,end);
    
    %%% range search in the subspace of lower dimension - M_Y
    [~, distance_Y] =  knnsearch(M_Y,M_Y,'K',N,'Distance',metric);
    count_Y = distance_Y(:,2:end) < dd;
    count_Y = max(k-1, sum(count_Y,2));

    %%% range search in the subspace of lower dimension - M_X
    [~, distance_X] =  knnsearch(M_X,M_X,'K',N,'Distance',metric);
    count_X = distance_X(:,2:end) < dd;
    count_X = max(k-1, sum(count_X,2));

    dd2=2*dd;
    dd2(dd2==0)=[]; % do not accept distance=0
    
    %%% OUTPUT
    IXY = psi(k)+psi(N)-(1/N)*(sum(psi(count_Y +1) + psi(count_X +1)));
    
end