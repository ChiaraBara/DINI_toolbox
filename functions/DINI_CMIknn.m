%% Computation of conditional mutual information measures
%%% INPUT
% B: data matrix (in column)
% iy: index of target variable
% ix: index of source variables
% iz: index of conditining variables
% k: number of neighbors
%%% OUTPUT
% IX: mutual information 
% IX_Z: conditional mutual information

function out = DINI_CMIknn(B,iy,ix,iz,k)

    metric='chebychev';

    M_Y = B(:,iy);
    M_X = B(:,ix);
    M_Z = B(:,iz);

    M_YZX = [M_Y M_Z M_X];
    M_YZ = [M_Y M_Z];
    M_YX = [M_Y M_X];
    M_ZX = [M_Z M_X];

    N = size(M_Y,1); 

    %%% neighbor search in the lower density space 
    [~, distances] =  knnsearch(M_YZX,M_YZX,'K',k+1,'Distance',metric);
    dd = distances(:,end);
    
    %%%%%%%%%%%%% mYZ, number of samples of [Y Z]
    [~, distance_YZ] =  knnsearch(M_YZ,M_YZ,'K',N,'Distance',metric);
    count_YZ = distance_YZ(:,2:end) < dd;
    count_YZ = max(k-1, sum(count_YZ,2));

    %%%%%%%%%%%%% mYX, number of samples of [Y X] 
    [~, distance_YX] =  knnsearch(M_YX,M_YX,'K',N,'Distance',metric);
    count_YX = distance_YX(:,2:end) < dd;
    count_YX = max(k-1, sum(count_YX,2));

    %%%%%%%%%%%%% mZX, number of samples of [Z X] 
    [~, distance_ZX] =  knnsearch(M_ZX,M_ZX,'K',N,'Distance',metric);
    count_ZX = distance_ZX(:,2:end) < dd;
    count_ZX = max(k-1, sum(count_ZX,2));

    %%%%%%%%%%%%% mY, number of samples of Y 
    [~, distance_Y] =  knnsearch(M_Y,M_Y,'K',N,'Distance',metric);
    count_Y = distance_Y(:,2:end) < dd;
    count_Y = max(k-1, sum(count_Y,2));

    %%%%%%%%%%%%% mX, number of samples of X 
    [~, distance_X] =  knnsearch(M_X,M_X,'K',N,'Distance',metric);
    count_X = distance_X(:,2:end) < dd;
    count_X = max(k-1, sum(count_X,2));

    %%%%%%%%%%%%% mZ, number of samples of Z 
    [~, distance_Z] =  knnsearch(M_Z,M_Z,'K',N,'Distance',metric);
    count_Z = distance_Z(:,2:end) < dd;
    count_Z = max(k-1, sum(count_Z,2));
    

    IX = psi(N) + (1/N)*sum(psi(count_YX+1) - psi(count_X+1) - psi(count_Y+1));
    IX_Z = psi(k) + (1/N)*sum(psi(count_Z+1) - psi(count_YZ+1) - psi(count_ZX+1));

    %%% OUTPUT
    out.IX = IX;
    out.IX_Z = IX_Z;

end
   
