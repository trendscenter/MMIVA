function [Dm, varargout] = do_variance_normalization(Dm)
% [Dm] = do_variance_normalization(Dm)
%     Do variance normalization of rows.
%     Assumes stats/pandas representation:
%         rows are observations, cols are data dimensionality (or variables).
% 
%     Parameters
%     ----------
%     Dm : matrix
%         The dataset to be normalized.
%
%     Returns
%     -------
%     Dm : matrix
%         The data matrix after applying variance normalization to its rows.
%     row_info : struct
%         Information about data rows. Contains the following fields:
%           - old_mean : the mean of the rows BEFORE varaince normalization
%           - old_std : the std_dev of the rows BEFORE varaince normalization
% 
%     Notes
%     -----
%     Assumes bad data rows/cols were already removed using detect_bad_data.m


% remove mean and divide by std for each row
old_mean_row = mean(Dm,2);
old_std_row = std(Dm,[],2);
Dm = bsxfun(@minus, Dm, old_mean_row);
Dm = bsxfun(@rdivide, Dm, old_std_row);

row_info.old_mean = old_mean_row;
row_info.old_std = old_std_row;

if nargout > 1
    varargout{1} = row_info;
end