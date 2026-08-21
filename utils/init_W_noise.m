function W0 = init_W_noise(P,W0,Y,seed)

M = 1:length(W0);
comps = size(Y{1},1);

rng(seed)

figure('color',[1 1 1],'visible','off')
for mm = M
    % if mod(mm-1,8) == 0
        tt = eye(size(W0{mm},1));
        W0{mm} = W0{mm}+((tt)*(P*1e-2)*(rand(size(W0{mm}))-.5));
    % else
        % W0{mm} = W0{mm-1};
    % end
    tmp = W0{mm} * Y{mm};
    W0{mm} = repmat(1.813799364234218 ./ std(tmp,[],2),1,size(W0{mm},2)) .* W0{mm};
    tmp = W0{mm} * Y{mm};
    subplot(2,max(M),mm)
    histogram(std(tmp,[],2),20)
    subplot(2,max(M),mm+length(M))
    Rs = corr(tmp');
    ix_off = logical(triu(ones(length(Rs)),1));
    histogram(Rs(ix_off),20)
    % histogram(diag(corr(S_hat(target,:)',tmp')),20)
end
saveas(gcf,sprintf('C%03d_S0_hist_P%03d.png',comps,P))

% Visualize cross-modal correlations:
figure('color',[1 1 1],'visible','off')
subplot(2,3,1)
imagesc(abs(corr((W0{1} * Y{1})',(W0{2} * Y{2})')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,2)
imagesc(abs(corr((W0{2} * Y{2})',(W0{3} * Y{3})')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,3)
imagesc(abs(corr((W0{1} * Y{1})',(W0{3} * Y{3})')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,4)
imagesc(abs(corr((W0{1} * Y{1})',(W0{1} * Y{1})')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,5)
imagesc(abs(corr((W0{2} * Y{2})',(W0{2} * Y{2})')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,6)
imagesc(abs(corr((W0{3} * Y{3})',(W0{3} * Y{3})')),[0 .5]), colormap(hot); colorbar
axis equal tight
saveas(gcf,sprintf('C%03d_S0_corr_P%03d.png',comps,P))

% Visualize initial weight distributions:
figure('color',[1 1 1],'visible','off')
for mm = M
    h = histogram(W0{mm}(:));
    h.Normalization = 'probability';
    h.BinWidth = 1e-3;
    hold on
end
saveas(gcf,sprintf('C%03d_W0hist_P%03d.png',comps,P))
