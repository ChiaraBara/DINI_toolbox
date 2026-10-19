%% Inference network algorithm
%%% INPUT
% B: data matrix (in column)
% k: number of neighbors
% surrtype: surrogate data approach (choices: circshift, shuff)
% nusurr: number of surrogate
% alpha: significance threshold for surrogate
%%% OUTPUT
% W: upper triangular matrix reporting the weights (range 0-1) of the
% direct link between each pair of nodes
% iz_min_tmp: upper triangulat matrix reporting the nodes selected to
% minimize the weight of the correspondent link that thus provide
% synergistic information

function [W, iz_min_tmp] = DINI_network_inference(B,k,surrtype,nsurr,alpha)

    Imin = zeros(size(B,2));
    I = zeros(size(B,2));
    W = zeros(size(B,2));
    
    for iy = 1:size(B,2)-1
    
        disp(['target ...',num2str(iy),'... of ... ',num2str(size(B,2))]);
    
        for ix = iy+1:size(B,2)
        
            disp(['source ...',num2str(ix),'... of ... ',num2str(size(B,2))]);
        
            [iz_min_tmp{iy,ix}, MI_check, CMI_check] = DINI_find_Imin_cond(B,iy,ix,k,surrtype,nsurr,alpha);
    
            if MI_check == 1
        
                B_tmp = [B(:,iy), B(:,ix),B(:,iz_min_tmp{iy,ix})]; 
                idx_x = 3:size(B_tmp,2); 
    
                if CMI_check == 0

                    W(iy,ix) = 0;

                else

                    if isempty(idx_x) ~= 1
                
                        out_min = DINI_CMIknn(B_tmp,1,2,idx_x,k);
                        Imin(iy,ix) = out_min.IX_Z;
                        I(iy,ix) = out_min.IX;
                
                        W(iy,ix) = Imin(iy,ix)/I(iy,ix);
        
                    else
                        
                        W(iy,ix) = 1;
        
                    end 
                    
                end
            
            else
    
                W(iy,ix) = 0 ;
    
            end
            
        end
    
    end

end