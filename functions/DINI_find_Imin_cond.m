%% Selection of subgroup of conditioning variables minimizing mutual information
%% and indication on MI and CMI statistical significance
%%% INPUT
% B: data matrix (in column)
% iy: index of target variable
% ix: index of source variables
% k: number of neighbors
% surrtype: surrogate data approach (choices: circshift, shuff)
% nusurr: number of surrogate
% alpha: significance threshold for surrogate
%%% OUTPUT
% iz_min_tmp: upper triangulat matrix reporting the nodes selected to
% minimize the weight of the correspondent link that thus provide
% synergistic information
% MIcheck: equal to 1 if MI is significant
% CMIcheck: equal to 1 if CMI is significant


function [iz_min, MIcheck, CMIcheck] = DINI_find_Imin_cond(B,iy,ix,k,surrtype,nsurr,alpha)

    minshift = 20;

    M_Y = B(:,iy); M_X = B(:,ix);
    iz = setdiff([1:size(B,2)],[iy,ix]);
    M_Z = B(:,iz);

    %%% searching for Zmin
    % first we have to compute I(Y;X) and check its significance: if it is
    % not significant the code end
    MI = MIknn(B,iy,ix,k);
    if surrtype == 'circshift'
        parfor isurr = 1:nsurr
            B_surr = B; B_surr(:,ix) = DINI_surrcircshift(B_surr(:,ix),minshift);
            MI_surr(isurr,1) = MIknn(B_surr,iy,ix,k);
        end
    elseif surrtype == 'shuff'
        parfor isurr = 1:nsurr
            B_surr = B; B_surr(:,ix) = DINI_surrshuf(B_surr(:,ix));
            MI_surr(isurr,1) = MIknn(B_surr,iy,ix,k);
        end
    end
    th_MI = prctile(MI_surr,(1-alpha)*100,1);
    MI_sign = MI>th_MI;
    
    iz_min = [];
    % we continue with the search of Zmin
    if MI_sign == 1  
        
        MIcheck = 1;
        CMIcheck = [];
        IZmin = MI;   
        min_check = 1;
        
        while min_check == 1
            iv_all = setdiff(iz,iz_min);
            if isempty(iv_all)
                min_check = 0;
                continue;
            else
                for i = 1:length(iv_all)
                    out = DINI_CMIknn_search(B,iy,ix,iz_min,iv_all(i),k);
                    CMI_ZV(i) = out.IX_ZV;
                    CMI_Z(i) = out.IX_Z;
                end
                [~,ind_iv_sel] = min(CMI_ZV);
                iv_sel = iv_all(ind_iv_sel);

                % test the significance of the CMI_ZV, indeed if it is not
                % significant the search stop and we consider CMI_Z as
                % minimum, contrarily, we proceed with the evaluation of
                % the significance of the difference

                diffCMI = CMI_Z(ind_iv_sel)-CMI_ZV(ind_iv_sel); 

                if diffCMI > 0        
                    if surrtype == 'circshift'
                        parfor isurr = 1:nsurr
                            B_surr = B; B_surr(:,iv_sel) = DINI_surrcircshift(B_surr(:,iv_sel),minshift);
                            out = CMIknn_search(B_surr,iy,ix,iz_min,iv_sel,k);
                            CMI_ZV_surr(isurr,1) = out.IX_ZV;
                            CMI_Z_surr(isurr,1) = out.IX_Z;
                        end
                    elseif surrtype == 'shuff'
                        parfor isurr = 1:nsurr
                            B_surr = B; B_surr(:,iv_sel) = DINI_surrshuf(B_surr(:,iv_sel));
                            out = CMIknn_search(B_surr,iy,ix,iz_min,iv_sel,k);
                            CMI_ZV_surr(isurr,1) = out.IX_ZV;
                            CMI_Z_surr(isurr,1) = out.IX_Z;
                        end
                    end
                    diffCMI_surr = CMI_Z_surr-CMI_ZV_surr;
                    th_diffCMI = prctile(diffCMI_surr,(1-alpha)*100,1);
                    if diffCMI>th_diffCMI
                        sign_diff = 1;
                    else
                        sign_diff = 0;
                    end

                    if sign_diff == 1

                        %%% if the difference is significant, I test for
                        %%% the significance of the CMImin
                        if surrtype == 'circshift'
                            parfor isurr = 1:nsurr
                                B_surr = B; B_surr(:,ix) = DINI_surrcircshift(B_surr(:,ix),minshift);
                                out_tmp = CMIknn(B_surr,iy,ix,[iz_min,iv_sel],k);
                                CMI_ZV_surr_tmp(isurr,1) = out_tmp.IX_Z;
                            end
                        elseif surrtype == 'shuff'
                            parfor isurr = 1:nsurr
                                B_surr = B; B_surr(:,ix) = DINI_surrshuf(B_surr(:,ix));
                                out_tmp = CMIknn(B_surr,iy,ix,[iz_min,iv_sel],k);
                                CMI_ZV_surr_tmp(isurr,1) = out_tmp.IX_Z;
                            end
                        end
                        th_CMI_tmp = prctile(CMI_ZV_surr_tmp,(1-alpha)*100,1);

                        if CMI_ZV(ind_iv_sel) > th_CMI_tmp
                            iz_min = [iz_min iv_sel];
                            IZmin = [IZmin CMI_ZV(ind_iv_sel)];
                        else
                            min_check = 0;
                            CMIcheck = 0;
                        end

                        clear CMI_ZV CMI_Z CMI_ZV_surr CMI_Z_surr
                    else
                        min_check = 0;
                    end
                else
                    min_check = 0;
                end

            end
        end
    
    else
        MIcheck = 0;
        CMIcheck = [];
    end


end

