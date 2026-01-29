%% GRACE2SLEPT
% Takes GRACE/GRACE-FO gravimetry data created by GRACE2PLM and projects
% this data into the requested Slepian basis.
%
% Syntax
%   [slept, stdSlept, dates] = GRACE2SLEPT(product, domain, L)
%   [slept, stdSlept, dates] = ...
%       GRACE2SLEPT(product, r, buf, L, phi, theta, omega, Truncation)
%   [slept, stdSlept, dates] = GRACE2SLEPT(__, 'Name', Value)
%   [__, domain, G, eigfun, V, N] = GRACE2SLEPT(__)
%
% Input arguments
%   product - The name of the data product
%       This is a cell array with three parts:
%       - the data center ('CSR', 'GFZ', or 'JPL'),
%       - the release level ('RL04', 'RL05', or 'RL06'),
%       - the dataproduct bandwidth.
%       The default product is {'CSR', 'RL06', 60}.
%       Data type: cell array
%   r - Radius of the concentration region in degrees
%       - A string of the domain of interest. It must be a function that
%           returns the coordinates of the domain.
%       - A cell array of the domain of interest (string) and the buffer
%           value (scalar): {domain, buf}.
%       - An N-by-2 matrix of the domain of interest in [lon, lat] format.
%       - A GeoDomain object.
%   buf - Distance in degrees that the region outline will be enlarged
%              by BUFFERM [default: 0]
%   L - Bandwidth of the window
%       - A scalar value for the highpass bandwidth of the window.
%       - A 2-element vector for the bandpass bandwidth of the window.
%       The default value is 60.
%       Data type: numeric (array)
%   phi, theta, omega - Longitude, colatitude, and anticlockwise azimuthal
%       rotation of the center of the tapers in degrees
%   Truncation - Number of largest eigenfunctions in which to expand
%       - 'N': The Shannon number of eigenfunctions.
%       The default is all of them.
%   Unit - Unit of the output
%       See GRACE2PLMT (GRACE2PLMT_NEW) for the available options. This
%       function does not check the vailidity of the input.
%       The default option is whatever the default is for GRACE2PLMT (as of
%       right now, SD).
%       Data type: char | []
%   LoveNumSource - Source of load Love numbers
%       - 'ISSM': Love numbers from the Ice Sheet System Model (ISSM)
%       - 'Wahr': Love numbers as used in slepian_delta package
%       - 'ALMA3 <model>': Love numbers from the ALMA3 models (currently
%           has to be added manually)
%       The default source is 'ISSM'.
%   Deg1Correction, C20Correction, C30Correction - Logical flag to apply
%       the degree 1, C20, or C30 corrections using the techinical notes
%       See GRACE2PLMT (GRACE2PLMT_NEW) for more details.
%       The default options are whatever the default is for GRACE2PLMT (as
%       of right now, true).
%       Data type: logical | []
%   RecomputeDegree1 - Logical flag to recompute the degree 1 coefficients
%       See GRACE2PLMT (GRACE2PLMT_NEW) for more details.
%       The default option is whatever the default is for GRACE2PLMT (as
%       of right now, no recomputation).
%       Data type: logical | char | cell array | ([])
%   TimeRange - Time range of the output
%       When specified, the output will be truncated to the specified time
%       range.
%       The default time range is [] (all available data).
%       Data type: datetime array | ([numeric] array)
%   ForceNew - Logical flag to force reprocess of the data
%       The default option is false.
%       Data type: logical | ([numeric])
%	SaveData - Logical flag to save the data
%		- true: Save the data to disk.
%		- false: Do not save the data to disk.
%       When recomputing the degree 1 coefficients, no data is saved.
%		The default option is true.
%		Data types: logical | ([numeric])
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is false.
%		Data types: logical | ([numeric])
%
% Output arguments
%   slept - Slepian coefficients of the gravity field
%       Units: m/s^2 (POT) | kg/m^2 (SD)
%       Data type: DOUBLE
%       Dimension: [ndates x truncation] | [truncation x ndates]
%           (depending on the OutputFormat input)
%   stdSlept - Standard deviation of the Slepian coefficients
%       The data type and dimension are the same as for slept.
%   dates - Time stamps of the coefficients
%       The time stamps are the midpoints of the time intervals of the
%       input data.
%       Datatype: DATATIME | DOUBLE
%           (depending on the TimeFormat input)
%       Dimension: [ndates x 1]
%   domain - The domain
%       If there was buffering, this will be a XY array of coordinates,
%       which you can use with SPHAREA to get the Shannon number.
%   G - The unitary matrix of localisation coefficients
%   eigfun - A cell array with cosine/sine coefficients eigenfunctions
%   V - The eigenvalues in this ordering
%   N - The Shannon number
%
% See also
%   GRACE2PLMT (GRACE2PLMT_NEW), PLM2SLEP
%
% Last modified by
%   2026/01/29, williameclee@arizona.edu (@williameclee)
%   2024/08/30, williameclee@arizona.edu (@williameclee)
%   2022/05/18, charig@princeton.edu (@harig00)
%   2012/06/26, fjsimons@alum.mit.edu (@fjsimons)

function varargout = grace2slept_new(varargin)
    %% Initialisation
    % Parse inputs
    [product, domain, L, phi, theta, omega, unit, loveNumSrc, timelim, ...
         truncation, ...
         deg1corr, c20corr, c30corr, redoDeg1, ...
         outputFmt, timeFmt, forceNew, saveData, beQuiet, callChain] = ...
        parseinputs(varargin{:});

    % Figure out if it's low-pass or band-pass
    bp = length(L) == 2;
    maxL = max(L);
    % The spherical harmonic dimension
    if ~bp
        ldim = (L + 1) ^ 2;
    else
        ldim = (L(2) + 1) ^ 2 - L(1) ^ 2;
    end

    % Check if you want the Shannon number of eigenfunctions
    if strcmp(truncation, 'N')
        truncation = ceil((L(end) + 1) ^ 2 * domain.Spharea);
    else
        truncation = conddefval(truncation, ldim);
    end

    % Output file
    outputPath = getoutputfile(domain, L, ...
        product, truncation, unit, bp, ...
        deg1corr, c20corr, c30corr, loveNumSrc);

    %% Loading existing data
    % If this expansion already exists, load it.  Otherwise, or if we force
    % it, make a new one (e.g. if you added extra months to the database).
    vars = {'slept', 'stdSlept', 'dates'};

    if ~forceNew && exist(outputPath, 'file') && all(ismember(vars, who('-file', outputPath)))
        data = load(outputPath, vars{:});

        if ~beQuiet
            fprintf('[ULMO>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">localised GRACE data</a>.\n', ...
                callchaintext(callChain), outputPath, outputPath);
        end

        [slept, stdSlept, dates] = ...
            formatoutput(data.slept, data.stdSlept, data.dates, timelim, outputFmt, timeFmt);

        if nargout <= 4
            varargout = {slept, stdSlept, dates, domain};
            return
        end

        [G, CC, V, N] = ...
            getslepianbasis(domain, L, phi, theta, omega, maxL, truncation, beQuiet);
        varargout = {slept, stdSlept, dates, domain, G, CC, V, N};

        return
    end

    %% Processing GRACE data
    % Use GRACE2PLMT to get the GRACE data
    [plmt, stdPlmt, dates] = ...
        grace2plmt_new(product, "Unit", unit, "LoveNumSource", loveNumSrc, ...
        "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "Deg1Correction", deg1corr, "C20Correction", c20corr, "C30Correction", c30corr, "RecomputeDegree1", redoDeg1, ...
        "BeQuiet", beQuiet, "CallChain", callChain);

    % Limit everything to the window bandwidth
    if size(plmt, 2) > addmup(maxL)
        plmt = plmt(:, 1:addmup(maxL), 1:4);
        stdPlmt = stdPlmt(:, 1:addmup(maxL), 1:4);
    elseif size(plmt, 2) < addmup(maxL)
        [order, degree] = addmon(maxL);
        plmt(:, 1:addmup(maxL), 1) = degree;
        plmt(:, 1:addmup(maxL), 2) = order;
    end

    [G, CC, V, N, ronmW] = ...
        getslepianbasis(domain, L, phi, theta, omega, maxL, truncation, beQuiet);

    nDates = length(dates);
    % GRACE
    plmst = reshape(plmt(:, :, 3:4), [nDates, size(plmt, 2) * 2]);
    plmst = plmst(:, ronmW(1:(maxL + 1) ^ 2)); % nDates x nCoeffs
    slept = (G' * plmst')';
    % GRACE STD
    stdPlmst = reshape(stdPlmt(:, :, 3:4), [nDates, size(stdPlmt, 2) * 2]);
    stdPlmst = stdPlmst(:, ronmW(1:(maxL + 1) ^ 2)); % nDates x nCoeffs
    stdSlept = sqrt(((G .^ 2)' * (stdPlmst .^ 2)')');

    if saveData
        save(outputPath, vars{:}, '-v7.3');

        if ~beQuiet
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">localised GRACE data</a>.\n', ...
                callchaintext(callChain), outputPath, outputPath);
        end

    end

    [slept, stdSlept, dates] = ...
        formatoutput(slept, stdSlept, dates, timelim, outputFmt, timeFmt);

    % Collect output
    varargout = {slept, stdSlept, dates, domain, G, CC, V, N};

end

%% Subfunctions
% Parse inputs
function varargout = parseinputs(varargin)
    dfOpt.Product = {'CSR', 'RL06', 60};
    dfOpt.Domain = 'greenland';
    dfOpt.L = 60;
    dfOpt.CapSpecs = 0;
    dfOpt.Truncation = [];
    dfOpt.Unit = [];
    dfOpt.LoveNumSrc = 'ISSM';
    dfOpt.OutputFormat = 'timefirst';
    dfOpt.TimeFormat = 'datetime';
    dfOpt.ForceNew = false;
    dfOpt.SaveData = true;

    ip = inputParser;
    addOptional(ip, 'Product', dfOpt.Product, ...
        @(x) (iscell(x) && length(x) == 3) || isempty(x));
    addOptional(ip, 'Domain', dfOpt.Domain, ...
        @(x) (ischar(x)) || isstring(x) || iscell(x) || ...
        isa(x, 'GeoDomain') || (isnumeric(x) && (size(x, 2) == 2 || isscalar(x))) || ...
        (isempty(x)));
    addOptional(ip, 'L', dfOpt.L, ...
        @(x) (isnumeric(x) && (length(x) <= 2)) || (isempty(x)));
    addOptional(ip, 'phi', dfOpt.CapSpecs, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(ip, 'theta', dfOpt.CapSpecs, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(ip, 'omega', dfOpt.CapSpecs, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isempty(x)));
    addOptional(ip, 'Truncation', dfOpt.Truncation, ...
        @(x) ((isnumeric(x) && isscalar(x) && x > 0) || ...
        strcmp(x, 'N')) || (isempty(x)));
    addOptional(ip, 'Unit', dfOpt.Unit, ...
        @(x) ischar(x) || isstring(x) || isempty(x));
    addOptional(ip, 'TimeRange', [], ...
        @(x) isempty(x) || ((isdatetime(x) || isnumeric(x))));
    addOptional(ip, 'ForceNew', dfOpt.ForceNew, ...
        @(x) ((isnumeric(x) || islogical(x)) && isscalar(x)) || isempty(x));
    addOptional(ip, 'MoreRegionSpecs', {}, @iscell);
    addParameter(ip, 'LoveNumSource', dfOpt.LoveNumSrc, @(x) ischar(x) || isstring(x));
    addParameter(ip, 'SaveData', [], ...
        @(x) ((isnumeric(x) || islogical(x)) && isscalar(x)) || isempty(x));
    addParameter(ip, 'BeQuiet', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'Deg1Correction', [], ...
        @(x) ((isnumeric(x) || islogical(x)) && isscalar(x)) || isempty(x));
    addParameter(ip, 'C20Correction', [], ...
        @(x) ((isnumeric(x) || islogical(x)) && isscalar(x)) || isempty(x));
    addParameter(ip, 'C30Correction', [], ...
        @(x) ((isnumeric(x) || islogical(x)) && isscalar(x)) || isempty(x));
    addParameter(ip, 'RecomputeDegree1', [], ...
        @(x) islogical(x) || iscell(x) || ischar(x) || isempty(x));
    addParameter(ip, 'OutputFormat', dfOpt.OutputFormat, ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'TimeFormat', dfOpt.TimeFormat, ...
        @(x) ischar(validatestring(x, {'datetime', 'datenum'})));
    addParameter(ip, 'CallChain', {}, @iscell);

    if (length(nargin) >= 3) && ...
            (isa(varargin{1}, {'char', 'string'}) && any(strcmpi(varargin{1}, {'JPL', 'CSR', 'GFZ'}))) && ...
            ((isa(varargin{2}, {'char', 'string'}) && startsWith(varargin{2}, 'RL')) || isnumeric(varargin{2})) && ...
            (isnumeric(varargin{3}) && varargin > 0)
        varargin = [{{varargin{1}, varargin{2}, varargin{3}}}, varargin(4:end)];
    end

    parse(ip, varargin{:});

    product = conddefval(ip.Results.Product, dfOpt.Product);
    domain = conddefval(ip.Results.Domain, dfOpt.Domain);
    L = conddefval(ip.Results.L, dfOpt.L);
    phi = conddefval(ip.Results.phi, dfOpt.CapSpecs);
    theta = conddefval(ip.Results.theta, dfOpt.CapSpecs);
    omega = conddefval(ip.Results.omega, dfOpt.CapSpecs);
    unit = char(conddefval(ip.Results.Unit, dfOpt.Unit));
    loveNumSource = ...
        char(conddefval(ip.Results.LoveNumSource, dfOpt.LoveNumSrc));
    domainSpecs = ip.Results.MoreRegionSpecs;
    J = conddefval(ip.Results.Truncation, dfOpt.Truncation);
    timeRange = ip.Results.TimeRange;
    forceNew = conddefval(logical(ip.Results.ForceNew), dfOpt.ForceNew);
    saveData = conddefval(logical(ip.Results.SaveData), dfOpt.SaveData);
    beQuiet = logical(ip.Results.BeQuiet);
    deg1corr = logical(ip.Results.Deg1Correction);
    c20corr = logical(ip.Results.C20Correction);
    c30corr = logical(ip.Results.C30Correction);
    redoDeg1 = ip.Results.RecomputeDegree1;
    outputFmt = ip.Results.OutputFormat;
    timeFmt = ip.Results.TimeFormat;
    callChain = [ip.Results.CallChain, {mfilename}];

    if isnumeric(product{2})
        product{2} = ['RL0', num2str(product{2})];
    end

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, 'file')
        domain = GeoDomain(domain, domainSpecs{:});
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2}, ...
            domainSpecs{:});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:}, domainSpecs{:});
    end

    disp(redoDeg1)

    if (~isempty(redoDeg1) && ~(islogical(redoDeg1) && ~redoDeg1)) && ...
            (~isempty(ip.Results.SaveData) && ip.Results.SaveData)
        warning(sprintf('%s:CannotSaveData', upper(mfilename)), ...
        'Saving data with recompute degree 1 is not supported')
        saveData = false;
    end

    varargout = ...
        {product, domain, L, ...
         phi, theta, omega, unit, loveNumSource, timeRange, J, ...
         deg1corr, c20corr, c30corr, redoDeg1, ...
         outputFmt, timeFmt, ...
         forceNew, saveData, beQuiet, callChain};
end

function outputPath = getoutputfile(domain, L, ...
        product, truncation, unit, bp, ...
        deg1corr, c20corr, c30corr, loveNumSrc)

    % Folder
    if ~isempty(getenv('GRACE'))
        outputFolder = fullfile(getenv('GRACE'), ...
        'SlepianExpansions');
    else
        outputFolder = fullfile(getenv('IFILES'), ...
            'GRACE', 'SlepianExpansions');
    end

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % File name
    if bp
        Lstr = sprintf('%i-%i', L(1), L(2));
    else
        Lstr = sprintf('%i', L);
    end

    productStr = [product{1}, '_', product{2}, '_', num2str(product{3})];

    if strcmpi(unit, 'SD')
        unit = sprintf('%s_%s', unit, loveNumSrc);
    end

    switch domainType(domain)
        case 'polar'
            r = domain;

            outputFile = sprintf( ...
                'grace2slept-%s-CAP-%i-%s-%i-%s.mat', ...
                productStr, r, Lstr, truncation, unit);

        case {'geodomain', 'lonlat'}

            switch domainType(domain)
                case 'geodomain'
                    domainId = domain.Id;
                case 'lonlat'
                    domainId = hash(domain, 'sha1');
            end

            % The name of the save file
            outputFile = sprintf( ...
                'grace2slept-%s-%s-%s-%i-%s.mat', ...
                productStr, domainId, Lstr, truncation, unit);

    end

    if ~deg1corr
        outputFile = strrep(outputFile, '.mat', '-nDeg1.mat');
    end

    if ~c20corr
        outputFile = strrep(outputFile, '.mat', '-nC20.mat');
    end

    if ~c30corr
        outputFile = strrep(outputFile, '.mat', '-nC30.mat');
    end

    outputPath = fullfile(outputFolder, outputFile);
end

% Get the auxiliary data (GLMALPHA, etc.)
function [G, CC, V, N, ronmW] = ...
        getslepianbasis(domain, L, phi, theta, omega, maxL, truncation, beQuiet)
    [~, ~, ~, lmcosiW, ~, ~, ~, ~, ~, ronmW] = addmon(maxL);

    if phi == 0 && theta == 0 && omega == 0
        [G, V, ~, ~, N] = glmalpha_new( ...
            domain, L, "J", truncation, "BeQuiet", beQuiet);
    else
        % Need to get a complete GLMALPHA but for the rotated basis
        % Definitely, "single-order" has lost its meaning here, but the MTAP
        % will still identify what the order of the unrotated original was
        [G, V, ~, ~, N] ...
            = glmalphapto(domain, L, phi, theta, omega);
        % Since GLMALPHAPTO currently has no option to limit a basis to J,
        % do it here
        G = G(:, 1:truncation);
    end

    % Sort by decreasing eigenvalue
    [V, sortId] = sort(V, 'descend');
    G = G(:, sortId);

    % If you don't do this, the eigenfunctions are ordered in the way
    %   that they correspond to single-orders back when, unrotated, they
    %   belonged to a polar cap, and the eigenvalues are sorted within
    %   these blocks. This is useful for, e.g. SPIE2009_1 a la SDSNEEUW.
    % Collect the eigenvector output into a format that PLM2XYZ knows how to interpret
    CC = cell(1, size(G, 2));

    for j = 1:size(G, 2)
        % Create the blanks
        cosi = lmcosiW(:, 3:4);
        % Stick in the coefficients of the 1st eigentaper
        cosi(ronmW) = G(:, j);
        % Construct the full matrix
        CC{j} = [lmcosiW(:, 1:2), cosi];
    end

end

% Format outputs
function [slept, stdSlept, dates] = ...
        formatoutput(slept, stdSlept, dates, timelim, outputFmt, timeFmt)

    if ~isempty(timelim)
        isValidTime = dates >= timelim(1) & dates <= timelim(2);
        slept = slept(isValidTime, :);
        stdSlept = stdSlept(isValidTime, :);
        dates = dates(isValidTime);
    end

    if strcmp(outputFmt, 'traditional')
        slept = permute(slept, [2, 1]);
        stdSlept = permute(stdSlept, [2, 1]);
    end

    if strcmp(timeFmt, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end

end

% Figute out which type of domain is used
function domainType = domainType(domain)

    if isa(domain, 'GeoDomain')
        domainType = 'geodomain';
    elseif isnumeric(domain) && isscalar(domain)
        domainType = 'polar';
    elseif isnumeric(domain) && size(domain, 2) == 2
        domainType = 'lonlat';
    else
        error('Unknown domain type')
    end

end
