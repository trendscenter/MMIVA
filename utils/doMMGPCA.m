function [whtM, H] = doMMGPCA(X, comps, rec_type)

M = 1:length(X);
N = size(X{1},2);

V = cellfun(@(x) size(x,1), X);

outs = cell(size(X));

cvx = zeros(N);

for mm = M
    cvx_ = cov(X{mm});
    cvx = cvx + cvx_./(length(M)*trace(cvx_)/N);
end

% [V,ll] = eig(cvx);
% ll = diag(ll);
% ll = ll(end:-1:1);

% Subject-level PCA reduction...
%[H,lambda] = eigs(cvx,comps,'largestabs','Tolerance',1e-10,'MaxIterations',1000,'SubspaceDimension',size(cvx,1));
[H,lambda] = eigs(cvx,comps);
% lambda = diag(lambda);
% max(abs(lambda - ll(1:length(lambda))))

% figure
% semilogy(ll(1:(end-1)),'r')
% hold on
% semilogy(lambda, '.b')
% saveas(gcf,'test.png')

A = cellfun(@(x) sqrt(N./(length(M)*sum(x(:).^2)))*(x*H), X, 'Un', 0);
norm_A = cellfun(@(a) sum(a.^2), A, 'Un', 0)';
norm_A = sqrt(sum(cell2mat(norm_A)));
A = cellfun(@(a) a./repmat(norm_A,size(a,1),1),A,'Un',0);
% norm_A = cellfun(@(a) sum(a.^2), A, 'Un', 0)';
% norm_A = sqrt(sum(cell2mat(norm_A)))

if strcmpi(rec_type, 'WT')
    whtM = cellfun(@(a) a', A, 'Un', 0);
elseif strcmpi(rec_type, 'PINV')
    whtM = cellfun(@(a) pinv(a), A, 'Un', 0);
elseif strcmpi(rec_type, 'REG')
    % To Do...
end

% whtM = cellfun(@(x,w) sqrt(N./(length(M)*sum(x(:).^2))) .* w, X, whtM, 'Un', 0);
whtM = cellfun(@(x,w) sqrt(N-1) * sqrt(N./(length(M)*sum(x(:).^2))) * repmat(1./(norm_A'),1,size(w,2)) .* w, X, whtM, 'Un', 0);
% tmp = cellfun(@(w,x) w * x, whtM, X, 'Un', 0);
% whtM = cellfun(@(w,t) repmat(1.813799364234218 ./ std(t,[],2),1,size(w,2)) .* w, whtM, tmp, 'Un', 0);

H = sqrt(N-1) * H';


%     if V(mm) > N
%         cvx = Dm*Dm'; % Scaling this by N-1 will still require adjustment in Wpca1 below bc VV always unit norm
%     else
%         cvx = Dm'*Dm; % cov(Dm); % Multiplying by sqrt(N-1) below instead, for better numerical stability
%     end
%     out = sort(abs(eig(cvx)));
%     out = out(end:-1:1);
%     outs{mm} = out.*(N/sum(out));
%     % out = cumsum(out);
%     % outs{mm} = out.*(N/out(end));
    
% end

% figure 
% semilogy(outs{1},'r')
% hold on
% semilogy(outs{2},'b')
% semilogy(outs{3},'g')
% % set(gca,'xlim',[0 1000])
% saveas(gcf,'test.png')
