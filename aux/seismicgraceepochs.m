%% SEISMICGRACEEPOCHS - Reads actual GRACE windows without filling missing months.
%
% Syntax
%   [epochs, equatorRadius] = seismicgraceepochs(product)
%
% Input arguments
%   product - Normalised GRACE {centre, release, degree} cell.
%       Source discovery follows GRACE2PLMT_NEW and ORIGINALGRACEDATA,
%       GRACEDATA/raw, then IFILES/GRACE/raw.
%
% Output arguments
%   epochs - Actual solution windows, sorted by midpoint. Size: (nDates, 2).
%   equatorRadius - Common source reference radius in metres.
%
% See also
%   grace2plmt_new, parsegracesourcefile, graceepochs
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

function [epochs, equatorRadius] = seismicgraceepochs(product)
    if ~isempty(getenv('ORIGINALGRACEDATA'))
        root = getenv('ORIGINALGRACEDATA');
    elseif ~isempty(getenv('GRACEDATA'))
        root = fullfile(getenv('GRACEDATA'), 'raw');
    else
        root = fullfile(getenv('IFILES'), 'GRACE', 'raw');
    end
    folder = fullfile(root, product{2}, product{1});
    if strcmp(product{1}, 'CSR mascon')
        pattern = 'GSU-2_*_B---_06*';
    elseif strcmp(product{2}, 'RL06')
        if product{3} == 60
            pattern = 'GSM-2_*_BA01_06*';
        else
            pattern = 'GSM-2_*_BB01_06*';
        end
    else
        switch product{1}
            case 'CSR', pattern = 'GSM-2_2*-2*_*_UTCSR_0060_0005*';
            case 'GFZ', pattern = 'GSM-2_2*-2*_*_EIGEN_G---_005a';
            case 'JPL', pattern = 'GSM-2_2*-2*_*_JPLEM_0000_0005*';
        end
    end
    files = dir(fullfile(folder, pattern));
    assert(~isempty(files), 'ULMO:seismicgraceepochs:MissingData', ...
        'No GRACE source files found in %s.', folder);
    epochs = NaT(numel(files), 2);
    radii = zeros(numel(files), 1);
    for iFile = 1:numel(files)
        [~, ~, ~, epochs(iFile, :), ~, radii(iFile)] = ...
            parsegracesourcefile(fullfile(files(iFile).folder, files(iFile).name));
    end
    [midpoints, order] = sort(mean(epochs, 2));
    epochs = epochs(order, :);
    assert(all(diff(midpoints) > seconds(0)) && all(epochs(:, 2) > epochs(:, 1)), ...
        'ULMO:seismicgraceepochs:InvalidWindows', 'GRACE windows must have unique midpoints and positive duration.');
    assert(all(abs(radii - radii(1)) < 1e-6), 'ULMO:seismicgraceepochs:Reference', ...
        'GRACE source files must use a common reference radius.');
    equatorRadius = radii(1);
end
