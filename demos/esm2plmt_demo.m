%% ESM2PLMT_DEMO - Displays ESM components for one GRACE epoch.
% Compares native and truncated H, I and O surface-density fields.
% All calculations stay in memory; no products or figures are saved.
%
% Syntax
%   esm2plmt_demo(graceFile)
%   esm2plmt_demo(graceFile, "Loutput", 45)
%
% Input arguments
%   graceFile - GRACE file readable by PARSEGRACESOURCEFILE.
%   Loutput (name-value) - Truncated output degree, no greater than the source.
%       Default degree: 60.
%   ColorLimits (name-value) - Positive symmetric H, I, O display limits.
%       Unit: mm equivalent water height.
%       Default limits: [500, 500, 100].
%
% See also
%   ESM2PLMT, PLOT_ESM_EPOCH, PARSEGRACESOURCEFILE
%
% Created by
%   2026/10/31, En-Chi Lee (williameclee@gmail.com)

function esm2plmt_demo(graceFile, options)

    arguments (Input)
        graceFile {mustBeTextScalar, mustBeFile}
        options.Loutput (1, 1) double {mustBeInteger, mustBePositive} = 60
        options.ColorLimits (1, 3) double {mustBePositive, mustBeFinite} = [500, 500, 100]
    end

    [~, ~, graceDate, intervals] = parsegracesourcefile(graceFile);
    fields = struct;
    nRows = (options.Loutput + 1) * (options.Loutput + 2) / 2;

    for comp = ["H", "I", "O"]
        [raw, ~, rawDates] = esm2plmt(comp, intervals, '', ...
            'Unit', 'SD', 'TimeFormat', 'datetime', 'SaveData', false);
        assert(options.Loutput <= raw(1, end, 1), ...
        'Demo output degree must not exceed the source degree.');
        % Average the exact reporting window, which may differ from a calendar month.
        [coeffs, dates] = averageepochs(rawDates, ...
            permute(raw(:, :, 3:4), [2, 3, 1]), intervals, hours(6), ...
            'CheckCompleteness', true, 'SamplingOrigin', datetime(2000, 1, 1));
        assert(all(isfinite(coeffs), 'all'), 'ULMO:esm2plmt:IncompleteSampling', ...
        'The prepared data do not cover the complete GRACE interval.');
        assert(isequal(dates, graceDate), 'ESM and GRACE epoch dates differ.');
        native = reshape(raw(1, :, :), [], 4);
        native(:, 3:4) = coeffs;
        fields.(comp).nativeSD = native;
        fields.(comp).truncatedSD = native(1:nRows, :);
    end

    plot_esm_epoch(fields, graceDate, 'ColorLimits', options.ColorLimits);
end
