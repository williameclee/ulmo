%% PLOT_ESM_EPOCH - Maps native and truncated ESM components side by side.
% Plots a figure for one epoch from ESM2PLMT_DEMO. Synthesises unmasked
% fields, so spectral spreading outside each source domain remains visible.
%
% Syntax
%   plot_esm_epoch(fields, epoch)
%   plot_esm_epoch(__, "ColorLimits", [500, 500, 100])
%
% Input arguments
%   fields - H, I and O structs with nativeSD and truncatedSD lmcosi matrices.
%       Unit: mm equivalent water height.
%   epoch - Scalar datetime for the epoch midpoint.
%   ColorLimits (name-value) - Positive symmetric H, I, O display limits.
%       Unit: mm equivalent water height.
%       Default limits: [500, 500, 100].
%
% Notes
%   Rows share one colour scale across resolutions. Values beyond each
%   display limit saturate; the data are not clipped. These are epoch mean
%   anomalies, not trends. The plotting grid has 0.5-degree spacing.
%
% See also
%   ESM2PLMT_DEMO, ESM2PLMT
%
% Created by
%   2026/10/31, En-Chi Lee (williameclee@gmail.com)

function plot_esm_epoch(fields, epoch, options)

    arguments (Input)
        fields (1, 1) struct
        epoch (1, 1) datetime
        options.ColorLimits (1, 3) double {mustBePositive, mustBeFinite} = [500, 500, 100]
    end

    % Arrange components by row and native/truncated products by column.
    comps = ["H", "I", "O"];
    names = ["Terrestrial hydrology", "Land ice", "Ocean"];
    variants = ["nativeSD", "truncatedSD"];
    grids = cell(3, 2);
    degrees = zeros(3, 2);

    for row = 1:3

        for column = 1:2
            coeffs = fields.(comps(row)).(variants(column));
            % Synthesize on one global grid without a coastal mask so spectral spreading
            % remains visible when comparing the two degree limits.
            [grids{row, column}, lon, lat] = plm2xyz(coeffs, ...
                0.5, [-180, 90, 180, -90], 'BeQuiet', true);
            assert(all(isfinite(grids{row, column}), 'all'));
            degrees(row, column) = coeffs(end, 1);
        end

    end

    %% Plotting
    coast = load('coastlines');
    figHandle = figure('Visible', 'on');
    sgtitle(figHandle, sprintf('ESAESM2 component mass anomalies (%s)', string(epoch, 'MMM yyyy')));
    % Interpolate a diverging palette centered on zero mass anomaly
    anchors = [0.10, 0.25, 0.55; 0.33, 0.65, 0.83; 0.97, 0.97, 0.94; ...
                   0.90, 0.49, 0.32; 0.58, 0.08, 0.16];
    colours = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, 257));

    for row = 1:3

        for column = 1:2
            ax = subplot(3, 2, 2 * (row - 1) + column, 'Parent', figHandle);
            imagesc(ax, lon(:), lat(:), grids{row, column});
            set(ax, 'YDir', 'normal', 'XTick', -180:60:180, 'YTick', -60:30:60, ...
                'Box', 'on', 'Layer', 'top');
            hold(ax, 'on');
            plot(ax, coast.coastlon, coast.coastlat, 'Color', [0.22, 0.24, 0.26]);
            xlim(ax, [-180, 180]);
            ylim(ax, [-90, 90]);
            pbaspect(ax, [2, 1, 1]);
            % Use the same symmetric limits for both resolutions of each component.
            clim(ax, [-options.ColorLimits(row), options.ColorLimits(row)]);
            colormap(ax, colours);

            if column == 1
                descriptor = "{\itL} = {\itL}_{data} =";
            else
                descriptor = "{\itL} =";
            end

            title(ax, sprintf('%s (%s %d)', names(row), descriptor, degrees(row, column)));

            bar = colorbar(ax, 'eastoutside');
            bar.Ticks = linspace(-options.ColorLimits(row), options.ColorLimits(row), 5);
            bar.Label.String = 'mm EWH';
        end

    end

end
