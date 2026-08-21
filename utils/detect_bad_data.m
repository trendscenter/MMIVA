function [exclusion_info, varargout] = detect_bad_data(Dm, varargin)
% [out] = detect_bad_data(Dm)
%     Detect bad rows/cols in matrix Dm.
%     Detection is based on
%       - the mode (most common) value along each row/col
%       - the standard deviation of each row/col
%     Also flags dataset for exclusion based on
%       - the total number of elements/cells excluded
%       - the expected number of rows/cols
% 
%     Assumes stats/pandas representation:
%         rows are observations, cols are data dimensionality (or variables).
%     Does not perform data mean removal prior to detection.
%     Does not perform data scaling prior to detection.
% 
%     Parameters
%     ----------
%     Dm : matrix
%         The dataset to be inspected.
%     mode_row_tsh : float in [0, 1] range (optional, default is 0.5)
%         The threshold for row exclusion.
%     mode_col_tsh : float in [0, 1] range (optional, default is 0.5)
%         The threshold for column exclusion.
%     total_tsh : float in [0, 1] range (optional) (optional, default is 0.7)
%         The threshold for full dataset exclusion.
%     min_std_row_tsh : float (optional, default is 10*sqrt(eps))
%         The threshold for minimal stddev of a row.
%     min_std_col_tsh : float (optional, default is sqrt(eps))
%         The threshold for minimal stddev of a col.
%     expected_num_rows : integer (optional, ignored if not specified)
%         The expected number of rows.
%     expected_num_cols : integer (optional, ignored if not specified)
%         The expected number of cols.
% 
%     Returns
%     -------
%     exclusion_info : struct
%         Information about detected bad data. Contains the following fields:
%           - bad_rows : logical column vector of rows to be excluded
%           - bad_cols : logical row vector of cols to be excluded
%           - exclude_dataset : logical indicating the entire dataset was flagged for excluson
%           - mismatch_row_dim : logical indicating mismatch between the expected and the actual number of rows
%           - mismatch_col_dim : logical indicating mismatch between the expected and the actual number of cols
%     Dm_out : matrix (optional)
%         The data matrix after applying row and column exclusions.
% 
%     Notes
%     -----
%     

if isempty(varargin)
    mode_row_tsh = .5;      % threshold for row exclusion
    mode_col_tsh = .5;      % threshold for column exclusion
    total_tsh = .7;         % threshold for full dataset exclusion
    min_std_row_tsh  = 10*sqrt(eps);    % threshold for minimal stddev of a row
    min_std_col_tsh  = sqrt(eps);       % threshold for minimal stddev of a column
    expected_num_rows = []; % expected number of rows
    expected_num_cols = []; % expected number of cols
else
    if length(varargin) >= 1
        mode_row_tsh = varargin{1};
        if isempty(mode_row_tsh)
            mode_row_tsh = .5; % Back to default
        else
            assert(isfloat(mode_row_tsh),' detect_bad_data : mode_row_tsh should be numeric ')
            assert(~isnan(mode_row_tsh),' detect_bad_data : mode_row_tsh is NaN ')
            assert(~isinf(mode_row_tsh),' detect_bad_data : mode_row_tsh is Inf ')
            if length(mode_row_tsh) > 1, mode_row_tsh = mode_row_tsh(1); end
            assert(mode_row_tsh >= 0,' detect_bad_data : mode_row_tsh is negative ')
            assert(mode_row_tsh <= 1,' detect_bad_data : mode_row_tsh is larger than 1 ')
        end
    end
    if length(varargin) >= 2
        mode_col_tsh = varargin{2};
        if isempty(mode_col_tsh)
            mode_col_tsh = .5; % Back to default
        else
            assert(isfloat(mode_col_tsh),' detect_bad_data : mode_col_tsh should be numeric ')
            assert(~isnan(mode_col_tsh),' detect_bad_data : mode_col_tsh is NaN ')
            assert(~isinf(mode_col_tsh),' detect_bad_data : mode_col_tsh is Inf ')
            if length(mode_col_tsh) > 1, mode_col_tsh = mode_col_tsh(1); end
            assert(mode_col_tsh >= 0,' detect_bad_data : mode_col_tsh is negative ')
            assert(mode_col_tsh <= 1,' detect_bad_data : mode_col_tsh is larger than 1 ')
        end
    end
    if length(varargin) >= 3
        total_tsh = varargin{3};
        if isempty(total_tsh)
            total_tsh = .7; % Back to default
        else
            assert(isfloat(total_tsh),' detect_bad_data : total_tsh should be numeric ')
            assert(~isnan(total_tsh),' detect_bad_data : total_tsh is NaN ')
            assert(~isinf(total_tsh),' detect_bad_data : total_tsh is Inf ')
            if length(total_tsh) > 1, total_tsh = total_tsh(1); end
            assert(total_tsh >= 0,' detect_bad_data : total_tsh is negative ')
            assert(total_tsh <= 1,' detect_bad_data : total_tsh is larger than 1 ')
        end
    end
    if length(varargin) >= 4
        min_std_row_tsh = varargin{4};
        if isempty(min_std_row_tsh)
            min_std_row_tsh = 10*sqrt(eps); % Back to default
        else
            assert(isfloat(min_std_row_tsh),' detect_bad_data : min_std_row_tsh should be numeric ')
            assert(~isnan(min_std_row_tsh),' detect_bad_data : min_std_row_tsh is NaN ')
            assert(~isinf(min_std_row_tsh),' detect_bad_data : min_std_row_tsh is Inf ')
            if length(min_std_row_tsh) > 1, min_std_row_tsh = min_std_row_tsh(1); end
            assert(min_std_row_tsh >= 0,' detect_bad_data : min_std_row_tsh is negative ')
        end
    end
    if length(varargin) >= 5
        min_std_col_tsh = varargin{5};
        if isempty(min_std_col_tsh)
            min_std_col_tsh = sqrt(eps); % Back to default
        else
            assert(isfloat(min_std_col_tsh),' detect_bad_data : min_std_col_tsh should be numeric ')
            assert(~isnan(min_std_col_tsh),' detect_bad_data : min_std_col_tsh is NaN ')
            assert(~isinf(min_std_col_tsh),' detect_bad_data : min_std_col_tsh is Inf ')
            if length(min_std_col_tsh) > 1, min_std_col_tsh = min_std_col_tsh(1); end
            assert(min_std_col_tsh >= 0,' detect_bad_data : min_std_col_tsh is negative ')
        end
    end
    if length(varargin) >= 6
        expected_num_rows = varargin{6};
        if isempty(expected_num_rows)
            expected_num_rows = []; % Back to default
        else
            assert(isnumeric(expected_num_rows),' detect_bad_data : expected_num_rows should be numeric ')
            assert(~isnan(expected_num_rows),' detect_bad_data : expected_num_rows is NaN ')
            assert(~isinf(expected_num_rows),' detect_bad_data : expected_num_rows is Inf ')
            if length(expected_num_rows) > 1, expected_num_rows = expected_num_rows(1); end
            assert(expected_num_rows >= 0,' detect_bad_data : expected_num_rows is negative ')
            assert(mod(expected_num_rows,1),' detect_bad_data : expected_num_rows is not integer ')
        end
    end
    if length(varargin) >= 7
        expected_num_cols = varargin{7};
        if isempty(expected_num_cols)
            expected_num_cols = []; % Back to default
        else
            assert(isnumeric(expected_num_cols),' detect_bad_data : expected_num_cols should be numeric ')
            assert(~isnan(expected_num_cols),' detect_bad_data : expected_num_cols is NaN ')
            assert(~isinf(expected_num_cols),' detect_bad_data : expected_num_cols is Inf ')
            if length(expected_num_cols) > 1, expected_num_cols = expected_num_cols(1); end
            assert(expected_num_cols >= 0,' detect_bad_data : expected_num_cols is negative ')
            assert(mod(expected_num_cols,1),' detect_bad_data : expected_num_cols is not integer ')
        end
    end
end


[N, V] = size(Dm);  % Number of observations & Dimensionality of the data
bad_rows = false(N,1);  % Initialize
bad_cols = false(1,V);  % Initialize

% Detect trouble rows with too many identical values over columns:
mode_row = mode(Dm,2);
bad_mode_row = sum(Dm == repmat(mode_row,1,V), 2) > floor(mode_row_tsh*V);

% Detect trouble column with too many identical values over rows:
mode_col = mode(Dm);
bad_mode_col = sum(Dm == repmat(mode_col,N,1)) > floor(mode_col_tsh*N);

bad_rows = logical(bad_mode_row);
bad_cols = logical(bad_mode_col);

% Detect near-zero-valued standard deviations:
% Each row
std_row = std(Dm(~bad_rows,~bad_cols),[],2);
bad_std_row = std_row < min_std_row_tsh;

% Each column
std_col = std(Dm(~bad_rows,~bad_cols));
bad_std_col = std_col < min_std_col_tsh;

bad_rows(~bad_rows) = logical(bad_std_row);
bad_cols(~bad_cols) = logical(bad_std_col);

% Detect dataset exclusion by:
if sum(~bad_rows)*sum(~bad_cols) < total_tsh*N*V
    exclude_dataset = true;
else
    exclude_dataset = false;
end

if ~isempty(expected_num_rows)
    if sum(bad_rows) ~= expected_num_rows
        mismatch_row_dim = true;
        exclude_dataset = true;
    else
        mismatch_row_dim = false;
    end
else
    mismatch_row_dim = [];
end

if ~isempty(expected_num_cols)
    if sum(bad_cols) ~= expected_num_cols
        mismatch_col_dim = true;
        exclude_dataset = true;
    else
        mismatch_col_dim = false;
    end
else
    mismatch_col_dim = [];
end

% Gather exclusion info:
exclusion_info.bad_rows = bad_rows;
exclusion_info.bad_cols = bad_cols;
exclusion_info.exclude_dataset = exclude_dataset;
exclusion_info.mismatch_row_dim = mismatch_row_dim;
exclusion_info.mismatch_col_dim = mismatch_col_dim;

% Apply masks to data:
if nargout > 1
    varargout{1} = Dm(~bad_rows,~bad_cols);
end
