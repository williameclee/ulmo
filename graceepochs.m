%% GRACEEPOCHS - Obtains observation windows for GRACE and GRACE-FO products.
% Reads GSM or CSR mascon solution windows,
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
%   [epochs, synthetic, equatorRadius] = GRACEEPOCHS(__)
%
% Input arguments
%   product - Product specification cell array {centre, release, degree}
%       - centre: 'CSR', 'GFZ', 'JPL' or 'CSR mascon'. Default is 'CSR'.
%       - release: 'RL05' or 'RL06', or numeric 5 or 6. Default is 'RL06'.
%       - degree: RL06 GSM degree 60 or 96. RL05 uses CSR degree 60 and
%           GFZ/JPL degree 90; CSR mascons use RL06 degree 720. These fixed
%           product degrees override the supplied value, as in GRACE2PLMT_NEW.
%           Default degree: 60.
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
%   equatorRadius - Common source reference equatorial radius.
%       Header parameters are read and checked only when this output is
%       requested. All source files for the product must agree within 1e-6 m,
%       including files outside timelim. Synthetic windows have no separate
%       reference radius.
%       Unit: m.
%       Size: scalar.
%
% Notes
%   Select epochs(~synthetic, :) to retain only actual solution windows.
%
% See also
%   FILLGRACEEPOCHS, AVERAGEEPOCHS, SSH2LONLATT, STERIC2LONLATT,
%   GRACE2PLMT, PARSEGRACESOURCEFILE
%
% Created by
%   2026/09/25, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

function [epochs, synthetic, equatorRadius] = graceepochs(product, timelim)

    arguments (Input)
        product (1, 3) cell = {'CSR', 'RL06', 60}
        timelim {mustBeTimeRange} = []
    end

    arguments (Output)
        epochs (:, 2) datetime
        synthetic (:, 1) logical
        equatorRadius (1, 1) double
    end

    centre = validatestring(product{1}, {'CSR', 'GFZ', 'JPL', 'CSR mascon'});
    release = product{2};

    if isnumeric(release)
        validateattributes(release, {'numeric'}, {'scalar', 'integer', 'finite'});
        release = sprintf('RL%02d', release);
    end
    release = upper(string(release));
    validateattributes(product{3}, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
    assert(isscalar(release) && ismember(release, ["RL05", "RL06"]), ...
        'ULMO:Epochs:UnsupportedProduct', 'Epoch discovery supports RL05 and RL06 products.');

    if strcmp(centre, 'CSR mascon')
        assert(release == "RL06", 'ULMO:Epochs:UnsupportedProduct', 'CSR mascons require RL06.');
        pattern = 'GSU-2_*_B---_06*';
    elseif release == "RL05"
        switch centre
            case 'CSR', pattern = 'GSM-2_2*-2*_*_UTCSR_0060_0005*';
            case 'GFZ', pattern = 'GSM-2_2*-2*_*_EIGEN_G---_005a';
            case 'JPL', pattern = 'GSM-2_2*-2*_*_JPLEM_0000_0005*';
        end
    else
        assert(ismember(product{3}, [60, 96]), 'ULMO:Epochs:UnsupportedProduct', ...
            'RL06 GSM epoch discovery supports degree 60 or 96.');
        patterns = {'GSM-2_*_BA01_06*', 'GSM-2_*_BB01_06*'};
        pattern = patterns{1 + (product{3} == 96)};
    end

    if ~isempty(getenv('ORIGINALGRACEDATA'))
        folder = fullfile(getenv('ORIGINALGRACEDATA'), release, centre);
    elseif ~isempty(getenv('GRACEDATA'))
        folder = fullfile(getenv('GRACEDATA'), 'raw', release, centre);
    else
        folder = fullfile(getenv('IFILES'), 'GRACE', 'raw', release, centre);
    end

    files = dir(fullfile(folder, pattern));

    if isempty(files)
        error('ULMO:Epochs:MissingReference', 'No reference GRACE files found in %s.', folder);
    end

    epochs = NaT(numel(files), 2);
    equatorRadius = NaN;
    radii = zeros(numel(files), 1);

    for iFile = 1:numel(files)
        path = fullfile(files(iFile).folder, files(iFile).name);
        if nargout >= 3
            [~, ~, ~, epochs(iFile, :), ~, radii(iFile)] = parsegracesourcefile(path);
        else
            [~, ~, ~, epochs(iFile, :)] = parsegracesourcefile(path);
        end
    end
    if nargout >= 3
        assert(all(abs(radii - radii(1)) < 1e-6), 'ULMO:Epochs:Reference', ...
            'GRACE source files must use a common reference radius.');
        equatorRadius = radii(1);
    end

    [epochs, synthetic] = fillgraceepochs(epochs, timelim);
end
