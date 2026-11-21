%% STERICCOVERAGE - Identifies ocean and deep-ocean cells in the prepared record.
%
% Syntax
%   masks = STERICCOVERAGE(source, gridNames, available, nDates)
%
% Input arguments
%   source - Structure of steric grids or a matfile object for the native grids.
%   gridNames - Grid variable names ordered by component and depth.
%       Depth columns are full-depth, shallow and deep.
%   available - Logical array identifying the supplied grids in gridNames.
%   nDates - Number of timestamps in the supplied grids.
%
% Output arguments
%   masks - Structure containing oceanMask and, for split products, deepOceanMask.
%       Each is a logical array in the spatial orientation of the supplied grids.
%       A finite sample at any supplied timestamp establishes coverage.
%
% See also
%   fillabsentdeepsteric, steric2lonlatt, steric2gmsl
%
% Created by
%   2026/11/21, En-Chi Lee (williameclee@gmail.com)

function masks = stericcoverage(source, gridNames, available, nDates)
    masks = struct();
    inMemory = isstruct(source);
    if ~inMemory
        sourcePath = source.Properties.Source;
        usePartialReads = H5F.is_hdf5(sourcePath);
    end

    maskNames = ["oceanMask", "deepOceanMask"];
    columns = {1:2, 3};

    for iMask = 1:numel(maskNames)
        fieldNames = gridNames(:, columns{iMask});
        fieldNames = fieldNames(available(:, columns{iMask}));
        covered = false;

        for name = fieldNames(:)'
            if inMemory
                covered = covered | any(isfinite(source.(name)), 3);
            elseif usePartialReads
                % Read bounded blocks rather than retaining whole records for GMSL.
                blockSize = 12;
                for firstDate = 1:blockSize:nDates
                    lastDate = min(firstDate + blockSize - 1, nDates);
                    block = source.(char(name))(:, :, firstDate:lastDate);
                    covered = covered | any(isfinite(block), 3);
                end
            else
                % Older MAT formats cannot read parts of a variable efficiently.
                field = load(sourcePath, char(name));
                covered = covered | any(isfinite(field.(name)), 3);
            end
        end

        if ~isempty(fieldNames)
            masks.(maskNames(iMask)) = covered;
        end
    end
end
