function [Dm, varargout] = remove_mean(Dm)
% [Dm] = remove_mean(Dm)
%     Remove mean per column.
%     Assumes stats/pandas representation:
%         rows are observations, cols are data dimensionality (or variables).
% 
%     Parameters
%     ----------
%     Dm : matrix
%         The dataset to be processed.
%
%     Returns
%     -------
%     Dm : matrix
%         The data matrix after removing the mean from its columns.
%     old_mean_col : row vector
%         The mean of each column of Dm BEFORE mean removal.
% 
%     Notes
%     -----
%     Assumes bad data rows/cols were already removed using detect_bad_data.m


% Remove mean per column
old_mean_col = mean(Dm);
Dm = bsxfun(@minus, Dm, old_mean_col);

if nargout > 1
    varargout{1} = old_mean_col;
end