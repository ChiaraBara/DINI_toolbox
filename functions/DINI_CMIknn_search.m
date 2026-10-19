%% Computation of conditional mutual information measures when conditioning 
%% on the already selected variables and when adding a new variable
%%% INPUT
% B: data matrix (in column)
% iy: index of target variable
% ix: index of source variables
% iz: index of already selected conditioning variables
% iv: index of tested conditioning variable
% k: number of neighbors
%%% OUTPUT
% IX_Z: conditional mutual information when considering already selected
% variables
% IX_Zv: conditional mutual information when considering already selected
% and tested variables


function out = DINI_CMIknn_search(B,iy,ix,iz,iv,k)
    
    metric='chebychev';
    
    M_Y = B(:,iy);
    M_X = B(:,ix);
    M_Z = B(:,iz);
    M_V = B(:,iv);
     
    M_YZVX = [M_Y M_Z M_V M_X];
    M_YZX = [M_Y M_Z M_X];
    M_YZV = [M_Y M_Z M_V];
    M_ZVX = [M_Z M_V M_X];
    M_YZ = [M_Y M_Z];
    M_ZX = [M_Z M_X];
    M_ZV = [M_Z M_V];
    
    N = size(M_Y,1); 
    

    %%% neighbor search in the lower density space 
    [~, distances] =  knnsearch(M_YZVX,M_YZVX,'K',k+1,'Distance',metric);
    dd = distances(:,end);
    
    %%%%%%%%%%%%% mYZX, number of samples of [Y Z X] 
    [~, distance_YZX] =  knnsearch(M_YZX,M_YZX,'K',N,'Distance',metric);
    count_YZX = distance_YZX(:,2:end) < dd;
    count_YZX = max(k-1, sum(count_YZX,2));

    %%%%%%%%%%%%% mYZV, number of samples of [Y Z V] 
    [~, distance_YZV] =  knnsearch(M_YZV,M_YZV,'K',N,'Distance',metric);
    count_YZV = distance_YZV(:,2:end) < dd;
    count_YZV = max(k-1, sum(count_YZV,2));

    %%%%%%%%%%%%% mZVX, number of samples of [Z V X] 
    [~, distance_ZVX] =  knnsearch(M_ZVX,M_ZVX,'K',N,'Distance',metric);
    count_ZVX = distance_ZVX(:,2:end) < dd;
    count_ZVX = max(k-1, sum(count_ZVX,2));

    %%%%%%%%%%%%% mYZ, number of samples of [Y Z] 
    [~, distance_YZ] =  knnsearch(M_YZ,M_YZ,'K',N,'Distance',metric);
    count_YZ = distance_YZ(:,2:end) < dd;
    count_YZ = max(k-1, sum(count_YZ,2));
    
    %%%%%%%%%%%%% mZX, number of samples of [Z X] 
    [~, distance_ZX] =  knnsearch(M_ZX,M_ZX,'K',N,'Distance',metric);
    count_ZX = distance_ZX(:,2:end) < dd;
    count_ZX = max(k-1, sum(count_ZX,2));

    %%%%%%%%%%%%% mZV, number of samples of [Z V] 
    [~, distance_ZV] =  knnsearch(M_ZV,M_ZV,'K',N,'Distance',metric);
    count_ZV = distance_ZV(:,2:end) < dd;
    count_ZV = max(k-1, sum(count_ZV,2));

    if isempty(M_Z)

        count_YX = count_YZX;
        count_Y = count_YZ;

        %%%%%%%%%%%%% mX, number of samples of X for 
        [~, distance_X] =  knnsearch(M_X,M_X,'K',N,'Distance',metric);
        count_X = distance_X(:,2:end) < dd;
        count_X = max(k-1, sum(count_X,2));

    else

        %%%%%%%%%%%%% mZ, number of samples of Z for 
        [~, distance_Z] =  knnsearch(M_Z,M_Z,'K',N,'Distance',metric);
        count_Z = distance_Z(:,2:end) < dd;
        count_Z = max(k-1, sum(count_Z,2));

    end

    IX_ZV = psi(k) + (1/N)*sum(psi(count_ZV+1) - psi(count_ZVX+1) - psi(count_YZV+1));
    if isempty(M_Z)
        IX_Z = psi(N) + (1/N)*sum(psi(count_YX+1) - psi(count_Y+1) - psi(count_X+1));
    else
        IX_Z = (1/N)*sum(psi(count_YZX+1) + psi(count_Z+1) - psi(count_ZX+1) - psi(count_YZ+1));
    end

    %%% OUTPUT
    out.IX_ZV = IX_ZV;
    out.IX_Z = IX_Z;
    
end