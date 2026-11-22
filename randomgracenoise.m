%% RANDOMGRACENOISE - Draws monthly GRACE Stokes coefficient noise from GFZ formal covariance factors.
%
% Syntax
%   [noise, dates] = randomgracenoise({'GFZ', 'RL06', 96})
%   [noise, dates] = randomgracenoise(product, timelim, Loutput, ...
%       "NumRealisations", 100, "Stream", stream)
%
% Input arguments
%   product - Native product specifier {'GFZ', 'RL06', 96}.
%   timelim (optional) - Inclusive bounds on the stored reference epochs.
%       Accepts datetime or serial date numbers, including infinite bounds.
%       Default range: [] (full record).
%   Loutput (optional) - Maximum output degree, from 2 to 96.
%       Uses the smallest processed file covering the request in GFZ/NEQ
%       under GRACEDATA, or under IFILES/GRACE when GRACEDATA is unset.
%       Default degree: 60.
%   NumRealisations (name-value) - Number of independent noise series.
%       Default count: 1.
%   Stream (name-value) - Scalar RandStream used for all draws.
%       Default stream: current MATLAB global stream.
%
% Output arguments
%   noise - Dimensionless Stokes perturbations in [degree, order, C, S]
%       format. Size: (lm, 4, time, realisation). Degrees 0 and 1 and
%       zonal sine slots have zero perturbation. Each month is independent.
%   dates - Selected GFZ reference epochs as unzoned UTC datetimes.
%       Size: (time, 1). Gaps in the record are preserved.
%
% See also
%   PROCESSGFZNORMALEQUATIONS, READSHNCTIME
%
% Created by
%   2026/11/10
function [noise, dates] = randomgracenoise(product, timelim, Loutput, options)
    arguments (Input)
        product (1, 3) cell
        timelim {mustBeTimeRange} = []
        Loutput (1, 1) double {mustBeMember(Loutput, 2:96)} = 60
        options.NumRealisations (1, 1) double {mustBeInteger, mustBePositive, mustBeFinite} = 1
        options.Stream (1, 1) RandStream = RandStream.getGlobalStream()
    end

    assert(isequal(product{3}, 96) && strcmpi(product{1}, 'GFZ') && strcmpi(product{2}, 'RL06'), ...
        'ULMO:randomgracenoise:UnsupportedProduct', 'Expected the native product {''GFZ'', ''RL06'', 96}.');
    [path, storedDegree] = resolveFactorFile(Loutput);
    dates = readshnctime(path);
    selected = true(size(dates));
    if ~isempty(timelim)
        if isnumeric(timelim)
            timelim = datetime(timelim, 'ConvertFrom', 'datenum');
        elseif ~isempty(timelim.TimeZone)
            timelim.TimeZone = 'UTC';
            timelim.TimeZone = '';
        end
        selected = dates >= timelim(1) & dates <= timelim(2);
    end
    epochs = find(selected);
    dates = dates(selected);
    [~, ~, ~, lmcosi] = addmon(Loutput);
    nEns = options.NumRealisations;
    noise = repmat(lmcosi, 1, 1, numel(epochs), nEns);
    if isempty(epochs)
        return
    end

    n = (Loutput + 1)^2 - 4;
    degree = double(ncread(path, 'degree', 1, n));
    order = double(ncread(path, 'order', 1, n));
    cs = double(ncread(path, 'cs', 1, n));
    nLm = size(lmcosi, 1);
    lmRows = degree .* (degree + 1) / 2 + order + 1;
    slots = lmRows + cs * nLm;

    % Each packed column runs to the last stored row. Copy only its first
    % n-col+1 entries to recover the leading lower factor of the marginal.
    nStored = (storedDegree + 1)^2 - 4;
    columnLengths = nStored - (0:n-1);
    starts = 1 + [0, cumsum(columnLengths(1:end-1))];
    lower = zeros(n);
    coeffs = zeros(2*nLm, nEns);
    for iEpoch = 1:numel(epochs)
        packed = ncread(path, 'noiseFactor', [1 epochs(iEpoch)], [starts(end) 1]);
        for col = 1:n
            lower(col:n, col) = packed(starts(col):starts(col)+n-col);
        end
        coeffs(slots, :) = lower * randn(options.Stream, n, nEns);
        noise(:, 3:4, iEpoch, :) = reshape(coeffs, nLm, 2, 1, nEns);
    end
end

function [path, degree] = resolveFactorFile(Loutput)
    folder = gfzneqfolder();
    listing = dir(fullfile(folder, 'GFZ_NEQ_L*.nc'));
    degrees = 2:96;
    names = compose("GFZ_NEQ_L%d.nc", degrees);
    available = degrees(ismember(names, {listing(~[listing.isdir]).name}));
    degree = min(available(available >= Loutput));
    if isempty(degree)
        error('ULMO:randomgracenoise:MissingFactors', ...
            'Requested degree %d in %s. Available degrees: %s.', Loutput, folder, mat2str(available));
    end
    path = fullfile(folder, sprintf('GFZ_NEQ_L%d.nc', degree));
end
