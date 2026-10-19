%% Generating surrogate data with random temporal shift

function xs=DINI_surrcircshift(x, minshift)
      maxshift = numel(x)-minshift;
      lagshift = fix(rand*(maxshift-minshift+1)+minshift);% casual shift between minshift e maxshift 
      xs=circshift(x,lagshift);%ith shift of the values
end
