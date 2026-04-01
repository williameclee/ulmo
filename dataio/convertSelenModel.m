%% CONVERTSELENMODEL - Converts Selen outputs to the format used by SLEPIAN
% Converts the geoid and vertical displacement outputs from the SELEN model
% to the format used by SLEPIAN_DELTA. This is not intended to be called
% directly by users.
%
% Syntax
%   convertSelenModel(modelFolder)
%
% Input arguments
%   modelFolder - folder containing the SELEN model outputs
%
% References
%   https://github.com/geodynamics/selen
%
% Author
%	2026/04/01, En-Chi Lee (williameclee@arizona.edu)

function convertSelenModel(modelFolder)

    arguments (Input)
        modelFolder {mustBeFolder}
    end

    dataFolder = fullfile(modelFolder, "FPR");
    % Detect the degree of the model from the model folder name
    [~, modelFolderSelf] = fileparts(modelFolder);
    modelName = regexp(modelFolderSelf, "RUN_([A-Za-z0-9-_]+)", "tokens", "once");
    modelInfo = regexp(modelFolderSelf, "RUN_([A-Za-z0-9]+)-R(\d+)-L(\d+)-I(\d+)", "tokens", "once");

    % Prepare grid for interpolation
    L = str2double(modelInfo{3});
    h = 360 / (2 * L);
    lon = 0:h:360;
    lat = -90:h:90;
    [lonn, latt] = meshgrid(lon, lat);
    %% Geoid and surface density
    geoidPath = fullfile(dataFolder, "gdot.pix");
    geoidXyz = readmatrix(geoidPath, "FileType", "text");
    geoidXyz(:, 3) = geoidXyz(:, 3) / 1e3; % mm -> m
    F = scatteredInterpolant(geoidXyz(:, 1), geoidXyz(:, 2), geoidXyz(:, 3), "nearest", "nearest");
    geoidMesh = F(mod(lonn, 360), latt);
    geoidPlm = xyz2plm_new(flip(geoidMesh), L);

    sdPlm = convertgravity(geoidPlm, "POT", "SD");

    % Save the converted model
    lmcosiM = sdPlm;
    save(fullfile(modelFolder, sprintf("%s_SD.mat", modelName)), "lmcosiM");
    lmcosiM = geoidPlm;
    save(fullfile(modelFolder, sprintf("%s_POT.mat", modelName)), "lmcosiM");

    %% Vertical displacement
    vlmPath = fullfile(dataFolder, "udot.pix");
    vlmXyz = readmatrix(vlmPath, "FileType", "text");
    F = scatteredInterpolant(vlmXyz(:, 1), vlmXyz(:, 2), vlmXyz(:, 3), "nearest", "nearest");
    vlmMesh = F(mod(lonn, 360), latt);
    vlmPlm = xyz2plm_new(flip(vlmMesh), L);

    % Save the converted model
    lmcosiM = vlmPlm;
    save(fullfile(modelFolder, sprintf("%s_VLM.mat", modelName)), "lmcosiM");
end
