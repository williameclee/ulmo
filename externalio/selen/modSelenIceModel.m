%% MODSELENICEMODEL - Modifies a ice model in SELEN format
% This function modifies the ice model used in a SELEN run by changing the
% volume change and time delay of a specified ice sheet. The modified model
% is saved as a new file in the same format as the original model. This
% function is intended to be used for sensitivity tests of the ice model in
% SELEN.
%
% Syntax
%   modSelenIceModel(model, icesheet, newModel, volumeChangeFrac, timeDelay)
%
% Input arguments
%   model - Name of the SELEN model to modify (without the .pix extension)
%       The model must be located in the DATA folder of the SELEN
%       installation.
%   icesheet (optional) - Name of the ice sheet to modify
%       See ICESHEETPOLY for recognised ice sheet names.
%       The default ice sheet to modify is "fennoscandian".
%   newModel (optional) - Name of the new model to save
%       If empty, the new name is programmatically generated based on the
%       original model name and the modifications made.
%   volumeChangeFrac (optional) - Fractional change in the volume of the
%       specified ice sheet
%       For example, a value of 0.1 means a 10% increase in volume, while a
%       value of -0.2 means a 20% decrease in volume.
%       The default value is 0, meaning no change in volume.
%   timeDelay (optional) - Time delay in years for the specified ice sheet
%       A positive value means the ice sheet changes later than in the
%       original model, while a negative value means it changes earlier.
%       The default value is 0, meaning no time delay.
%
% See also
%   ICESHEETPOLY
%
% References
%   Spada, Giorgio & Melini, Daniele. (2019). SELEN4 (SELEN version 4.0): a
%       Fortran program for solving the gravitationally and topographically
%       self-consistent sea-level equation in glacial isostatic adjustment
%       modeling. Geosci. Model Dev.
%       https://doi.org/10.5194/gmd-12-5055-2019
%
% Author
%   2026/04/07, En-Chi Lee (williameclee@arizona.edu)

function modSelenIceModel(model, icesheet, newModel, volumeChangeFrac, timeDelay)

    arguments (Input)
        model {mustBeTextScalar} = "i6g-R44r-i"
        icesheet {mustBeTextScalar} = "fennoscandian"
        newModel {mustBeTextScalar} = ""
        volumeChangeFrac (1, 1) double {mustBeReal, mustBeFinite} = 0
        timeDelay (1, 1) double {mustBeInteger, mustBeReal, mustBeFinite} = 0
    end

    modelPath = fullfile(getenv("SELEN"), "DATA", sprintf("%s.pix", model));
    icesheet = icesheetCode(icesheet);

    ice = readmatrix(modelPath, "FileType", "text");
    isIcesheet = labelIceSheet(ice(:, 2:3), icesheet);

    iceDiff = ice(:, 5:end) - ice(:, end);

    if volumeChangeFrac ~= 0
        iceDiff(isIcesheet, :) = iceDiff(isIcesheet, :) * (1 + volumeChangeFrac);
    end

    if timeDelay > 0
        iceDiff(isIcesheet, timeDelay + 1:end) = iceDiff(isIcesheet, 1:end - timeDelay);

        iceDiff(isIcesheet, 1:timeDelay) = repmat(iceDiff(isIcesheet, 1), [1, timeDelay]);

    elseif timeDelay < 0
        iceDiff(isIcesheet, 1:end + timeDelay) = iceDiff(isIcesheet, -timeDelay + 1:end);
    end

    iceNew = [ice(:, 1:4), ice(:, end) + iceDiff];

    if ~isempty(newModel)
        newModel = sprintf("%s-%s_V%s_T%s", newModel, upper(icesheet), ...
            replace(replace(sprintf("%+03.0f", volumeChangeFrac * 100), "+", "p"), "-", "n"), ...
            replace(replace(sprintf("%+03d", timeDelay), "+", "p"), "-", "n"));
    end

    newModelPath = fullfile(getenv("SELEN"), "DATA", newModel);

    fid = fopen(newModelPath, "w");
    cleanUp = onCleanup(@() fclose(fid));

    for i = 1:size(iceNew, 1)
        fprintf(fid, "%12d %20f %20f %25f     ", iceNew(i, 1:4));
        fprintf(fid, " %11.0f", iceNew(i, 5:end));
        fprintf(fid, "\n");
    end

end

%% Subfunctions
function code = icesheetCode(name)

    arguments (Input)
        name {mustBeTextScalar}
    end

    arguments (Output)
        code (1, 1) string
    end

    switch lower(name)
        case {"fennoscandian", "fennoscandian ice sheet", "fis"}
            code = "FIS";
        case {"laurentide", "laurentide ice sheet", "lis"}
            code = "LIS";
        case {"greenland", "greenland ice sheet", "gis"}
            code = "GIS";
        case {"patagonian", "patagonian ice sheet", "pis"}
            code = "PIS";
        case {"antarctic", "antarctic ice sheet", "ais"}
            code = "AIS";
        case {"west antarctic", "west antarctic ice sheet", "wais"}
            code = "WAIS";
        case {"east antarctic", "east antarctic ice sheet", "eais"}
            code = "EAIS";
        otherwise
            error("Unrecognised ice sheet name: %s", name)
    end

end

function isIcesheet = labelIceSheet(lonlat, icesheet)

    arguments (Input)
        lonlat (:, 2) {mustBeNumeric}
        icesheet {mustBeTextScalar}
    end

    arguments (Output)
        isIcesheet (:, 1) logical
    end

    isIcesheet = isinterior(icesheetPoly(icesheet), lonlat(:, 1), lonlat(:, 2));

end
