%% GRACEEPOCHS
% Obtains observation time windows for GRACE/GRACE-FO Level-2 products,
% supplementing missing solution months and the mission gap with aligned
% calendar/synthetic windows.
%
% Real solution windows remain unchanged. For missing months, outer bounds
% meet adjacent actual solution windows, while internal boundaries align with
% calendar month starts. Epoch midpoints (mean(epochs, 2)) match the
% timestamps returned by GRACE2PLMT / GRACE2PLMT_NEW. All real windows are
% kept, including multiple solutions within the same month.
%
% Syntax
%   epochs = GRACEEPOCHS
%   epochs = GRACEEPOCHS(product)
%   epochs = GRACEEPOCHS(product, timelim)
%   [epochs, synthetic] = GRACEEPOCHS(__)
%
% Input arguments
%   product - Product specification cell array {centre, release, degree}
%       - centre: Data centre ('CSR', 'GFZ', or 'JPL'). Default is 'CSR'.
%       - release: Solution release level ('RL06' or 6). Default is 'RL06'.
%       - degree: Maximum spherical harmonic degree (60 or 96). Default is 60.
%       The default cell array is {'CSR', 'RL06', 60}.
%       Data type: CELL
%       Dimension: [1 x 3]
%   timelim - Time range of interest
%       A two-element vector [startTime, endTime] specifying the time range
%       to restrict output epochs. If empty, the full time span of available
%       data is returned.
%       The default is [] (all available epochs).
%       Data type: DATETIME | DOUBLE (datenum) | []
%       Dimension: [1 x 2] | [2 x 1] | []
%
% Output arguments
%   epochs - Start and end timestamps of each solution window
%       Start time is in the first column and end time is in the second
%       column. Midpoints (mean(epochs, 2)) correspond to the epoch dates.
%       Data type: DATETIME
%       Dimension: [N x 2]
%   synthetic - Logical flag indicating synthetic (gap-filled) windows
%       true for synthetic calendar windows inserted during missing months
%       or the GRACE/GRACE-FO mission gap; false for actual GRACE/GRACE-FO
%       solution windows.
%       Data type: LOGICAL
%       Dimension: [N x 1]
%
% See also
%   FILLGRACEEPOCHS, AVERAGEEPOCHS, SSH2LONLATT, STERIC2LONLATT,
%   GRACE2PLMT, PARSEGRACESOURCEFILE
%
% Created by
%   2026/09/25, En-Chi Lee (williameclee@arizona.edu)

function [epochs, synthetic] = graceepochs(product, timelim)

    arguments (Input)
        product (1, 3) cell = {'CSR', 'RL06', 60}
        timelim {mustBeTimeRange} = []
    end

    arguments (Output)
        epochs (:, 2) datetime
        synthetic (:, 1) logical
    end

    centre = char(product{1});
    release = product{2};

    if isnumeric(release)
        release = sprintf('RL%02d', release);
    end

    if ~isempty(getenv('ORIGINALGRACEDATA'))
        folder = fullfile(getenv('ORIGINALGRACEDATA'), release, centre);
    elseif ~isempty(getenv('GRACEDATA'))
        folder = fullfile(getenv('GRACEDATA'), 'raw', release, centre);
    else
        folder = fullfile(getenv('IFILES'), 'GRACE', 'raw', release, centre);
    end

    if ~strcmpi(release, 'RL06') || ~ismember(product{3}, [60, 96])
        error('ULMO:Epochs:UnsupportedProduct', ...
        'Automatic epoch discovery supports RL06 degree 60 or 96 GSM products.');
    end

    patterns = {'GSM-2_*_BA01_06*', 'GSM-2_*_BB01_06*'};
    files = dir(fullfile(folder, patterns{1 + (product{3} == 96)}));

    if isempty(files)
        error('ULMO:Epochs:MissingReference', 'No reference GRACE files found in %s.', folder);
    end

    epochs = NaT(numel(files), 2);

    for k = 1:numel(files)
        [~, ~, ~, epochs(k, :)] = parsegracesourcefile(fullfile(files(k).folder, files(k).name));
    end

    [epochs, synthetic] = fillgraceepochs(epochs, timelim);
end
