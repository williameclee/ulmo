%% PLOTSEALEVELTSERIES
% Plots a quick and dirty time series of global mean sea level from gridded data
%
% Notes
%	This is a dirty helper function not meant for rigorous scientific analysis.
%
% See also
%	STERIC2LONLATT, SSH2LONLATT
%
% Last modified by
%	2025/11/03, williameclee@arizona.edu (@williameclee)

function msl = plotsealeveltseries(dates, sl, lon, lat, product, unit, figtitle, options)

    arguments (Input)
        dates datetime {mustBeVector}
        sl (:, :, :) {mustBeNumeric}
        lon {mustBeReal, mustBeVector}
        lat {mustBeReal, mustBeVector}
        product (1, :) char = ''
        unit (1, :) char = ''
        figtitle (1, :) char = ''
        options.Visible (1, 1) logical = true
    end

    if size(sl, 1) == length(lon)
        % ndgrid -> meshgrid
        sl = permute(sl, [2, 1, 3]);
    end

    domain = GeoDomain('oceans', "Buffer", 0.5, "DefaultParams", true);
    domainMask = domainmask(domain, lon, lat, "Format", 'meshgrid', "BeQuiet", true);
    gridWeights = domainMask .* cosd(lat(:));
    gridWeights(~domainMask) = 0;
    gridWeights = gridWeights ./ sum(gridWeights, 'all', 'omitmissing');
    msl = squeeze(sum(sl .* gridWeights, [1, 2], 'omitmissing'));
    msl = msl - mean(msl, 'omitmissing');

    if ~options.Visible
        return
    end

    figure()

    if ~isempty(figtitle)
        set(gcf, "Name", figtitle, 'NumberTitle', 'off');
    end

    clf

    p = plot(dates, msl);

    if ~isempty(product)
        p.DisplayName = sprintf('Product: %s', char(product));
    end

    if isempty(unit)
        ylabel('Sea level [%s]')
    else
        ylabel(sprintf('Sea level [%s]', unit))
    end

    if ~isempty(figtitle)
        title(figtitle)
    end

    if ~isempty(product)
        legend('Location', 'best')
    end

    xlim([min(dates), max(dates)])
end
