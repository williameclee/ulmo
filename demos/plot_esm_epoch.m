%% PLOT_ESM_EPOCH - Maps native and truncated ESM components side by side.
% Saves a figure for one epoch from ESM2PLMT_DEMO. Synthesizes unmasked
% fields, so spectral spreading outside each source domain remains visible.
%
% Syntax
%   plot_esm_epoch(products, outputFile)
%   plot_esm_epoch(__, "ColorLimits", [500, 500, 100])
%
% Input arguments
%   products - Product paths returned by ESM2PLMT_DEMO. H, I and O must each
%       provide nativeSD and truncatedSD paths to one-epoch MAT files.
%   outputFile - Destination PNG file in an existing folder.
%   ColorLimits (name-value) - Positive symmetric H, I, O display limits.
%       Unit: mm equivalent water height. Default limits: [500, 500, 100].
%
% Notes
%   Rows share one color scale across resolutions. Values beyond each
%   display limit saturate; the data are not clipped. These are epoch mean
%   anomalies, not trends. The plotting grid has 0.5-degree spacing.
%   Requires the MATLAB coastlines example data and PLM2XYZ dependencies.
%
% See also
%   ESM2PLMT_DEMO, ESM2PLMT, PLM2XYZ
%
% Created by
%   2026/10/29

function plot_esm_epoch(products, outputFile, options)
    arguments (Input)
        products (1, 1) struct
        outputFile (1, 1) string
        options.ColorLimits (1, 3) double {mustBePositive, mustBeFinite} = [500, 500, 100]
    end
    components = ["H", "I", "O"];
    names = ["Terrestrial hydrology", "Land ice", "Ocean"];
    variants = ["nativeSD", "truncatedSD"];
    fields = cell(3, 2);
    degrees = zeros(3, 2);
    referenceDate = [];
    for row = 1:3
        for column = 1:2
            data = load(products.(components(row)).(variants(column)));
            assert(size(data.esmPlmt, 1) == 1 && size(data.esmPlmt, 3) == 4 && ...
                isscalar(data.dates), ...
                'Expected one timefirst SD epoch of the requested component.');
            if isempty(referenceDate)
                referenceDate = data.dates;
            else
                assert(isequal(data.dates, referenceDate), 'ESM epoch dates differ.');
            end
            coefficients = reshape(data.esmPlmt, [], 4);
            [fields{row, column}, lon, lat] = plm2xyz(coefficients, ...
                0.5, [-180, 90, 180, -90], 'BeQuiet', true);
            assert(all(isfinite(fields{row, column}), 'all'));
            degrees(row, column) = coefficients(end, 1);
        end
    end
    coast = load('coastlines');
    figureHandle = figure('Visible', 'off', 'Color', 'w', 'Position', [50, 50, 1500, 1150]);
    cleanup = onCleanup(@() close(figureHandle));
    annotation(figureHandle, 'textbox', [0.02, 0.956, 0.96, 0.04], ...
        'String', 'ESAESM2  |  Component mass anomalies', 'FontSize', 23, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    epoch = datetime(referenceDate, 'ConvertFrom', 'datenum');
    caption = "Epoch midpoint: " + string(epoch, 'dd MMM yyyy HH:mm') + " UTC" + newline + ...
        "Mean mass anomalies, not trends  |  Equivalent water height [mm]";
    annotation(figureHandle, 'textbox', [0.02, 0.904, 0.96, 0.05], 'String', caption, ...
        'FontSize', 13, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    anchors = [0.10, 0.25, 0.55; 0.33, 0.65, 0.83; 0.97, 0.97, 0.94; ...
        0.90, 0.49, 0.32; 0.58, 0.08, 0.16];
    colors = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, 257));
    left = [0.055, 0.55];
    bottom = [0.645, 0.355, 0.065];
    for row = 1:3
        for column = 1:2
            position = [left(column), bottom(row), 0.365, 0.238];
            ax = axes(figureHandle, 'Position', position);
            imagesc(ax, lon(:), lat(:), fields{row, column});
            set(ax, 'YDir', 'normal', 'FontSize', 11, 'TickDir', 'out', ...
                'XTick', -180:60:180, 'YTick', -60:30:60, 'Box', 'on', ...
                'Layer', 'top', 'LineWidth', 0.6);
            hold(ax, 'on');
            plot(ax, coast.coastlon, coast.coastlat, 'Color', [0.22, 0.24, 0.26], 'LineWidth', 0.45);
            xlim(ax, [-180, 180]);
            ylim(ax, [-90, 90]);
            pbaspect(ax, [2, 1, 1]);
            clim(ax, [-options.ColorLimits(row), options.ColorLimits(row)]);
            colormap(ax, colors);
            if column == 1
                descriptor = "native degree ";
            else
                descriptor = "truncated to degree ";
            end
            title(ax, names(row) + "  |  " + descriptor + degrees(row, column), ...
                'FontSize', 14, 'FontWeight', 'normal');
            if row == 3
                xlabel(ax, 'Longitude [°]');
            end
            if column == 1
                ylabel(ax, 'Latitude [°]');
            end
            bar = colorbar(ax, 'eastoutside');
            bar.Ticks = linspace(-options.ColorLimits(row), options.ColorLimits(row), 5);
            bar.Label.String = 'mm';
            bar.FontSize = 10;
            ax.Position = position;
            bar.Position = [left(column) + 0.375, bottom(row), 0.012, 0.238];
        end
    end
    annotation(figureHandle, 'textbox', [0.01, 0.001, 0.98, 0.025], 'String', ...
        ['Same color scale within each row; values beyond the limits saturate.  ', ...
        'No coastal mask or SLE adjustment applied.'], 'FontSize', 12, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    drawnow;
    exportgraphics(figureHandle, outputFile, 'Resolution', 160);
end
