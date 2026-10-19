clear all; close all; clc;

addpath([pwd,'\functions\']);

%%% simulation parameters
c=[0 0.25 0.5 0.75 1];
M=5;
N = 300; %time series length
nsim = 100;

%%% estimation paramters
k = 10;
surrtype = 'circshift';  %% choice: 'surrshift','shuff'
nsurr = 100;
alpha = 0.05;

for isim = 1:nsim

    disp(['sim ...',num2str(isim),'... of ...',num2str(nsim)]);

    B0 = rand(N,M);
    
    for ic=1:numel(c)
    
        B = B0;
    
        B(:,2) = B(:,2) + c(ic)*B(:,1);
        B(:,3) = B(:,3) + (1-c(ic))*B(:,5);
        B(:,4) = B(:,4) + c(ic)*B(:,3) + 0.5*B(:,5);
        
        B = zscore(B);
    
        [W{isim,ic}, iz_min_tmp{isim,ic}] = DINI_network_inference(B,k,surrtype,nsurr,alpha);
    
    end

    clear B0 B

end

%% plot 

th = 0;

for ic = 1:length(c)
    for isim = 1:nsim
        W_tmp{1,ic}(:,:,isim) = W{isim,ic};    
    end
    Wmean{1,ic} = mean(W_tmp{1,ic},3);
end

n_node = size(Wmean{1},2);

% define nodes and edge connections 
nodeedge = nchoosek(1:n_node, 2);
s = nodeedge(:,1)';
t = nodeedge(:,2)';

f1 = figure; 
for ic = 1:length(c)

    W_tmp = Wmean{1,ic}; mask = tril(true(n_node)); 
    W_tmp(mask) = NaN;
    W_tmp = reshape(W_tmp',1,[]);
    W_tmp = W_tmp(~isnan(W_tmp));
    
    % create the weighted undirected graph
    G = graph(s, t, W_tmp);
    G_reduced = rmedge(G, find(G.Edges.Weight == 0));
    
    % plot the graph
    subplot(1,5,ic)
    p = plot(G_reduced,'Layout','circle','EdgeCData', G_reduced.Edges.Weight, 'LineWidth', 1.5);
    clim([0 1]);
    colormap(turbo); 
    p.NodeColor = 'k';
    p.MarkerSize = 2;
    axis square; axis off;
    title(['c = ',num2str(c(ic))]);

    if ic == length(c)
        colorbar;
    end

    clear G G_reduced
end
f1.WindowState = 'maximized';

f2 = figure; 
for ic = 1:length(c)

    W_cond = zeros(n_node);
    mask = tril(true(n_node)); 
    W_cond(mask) = NaN;
    W_cond = reshape(W_cond',1,[]);
   
    for isim = 1:nsim
        W_tmp = W{isim,ic}; 
        pos_sign = find(W_tmp'>th);
        W_cond(pos_sign) = W_cond(pos_sign) + 1;
    end

    W_cond = W_cond(~isnan(W_cond));
        
    % create the weighted undirected graph
    G = graph(s, t, W_cond);
    G_unweighted = rmedge(G, find(G.Edges.Weight <= 10));

    scaledWidths = 4 * G_unweighted.Edges.Weight / max(G_unweighted.Edges.Weight);

    % plot the graph
    subplot(1,5,ic)
    p = plot(G_unweighted,'Layout','circle','LineWidth', scaledWidths, 'EdgeLabel', G_unweighted.Edges.Weight);
    clim([0 nsim]);
    p.EdgeColor = 'k';
    p.NodeColor = 'k';
    p.MarkerSize = 2;
    axis square; axis off;
    title(['c = ',num2str(c(ic))]);

    clear G G_unweighted

end
f2.WindowState = 'maximized';



