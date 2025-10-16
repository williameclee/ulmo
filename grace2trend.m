%% GRACE2TREND
% Reads GRACE data and computes the mass trend projected onto a set of
% Slepian functions, with optional corrections applied.
%
% Syntax
%   [date, total] = grace2trend(domain)
%   [date, total] = grace2trend(domain, L, timelim, fitwhat)
%   [date, total] = grace2trend(domain, L, timelim, fitwhat, product, gia)
%   [date, total] = grace2trend(__, 'Name', Value)
%   [date, total, totalerror, totalfit, f, ferror] = grace2trend(__)
%   [__, slept, sleptSig, sleptRes] = grace2trend(__)
%
% Input arguments
%   domain - The domain of interest
%       - A string of the domain of interest. It must be a function that
%           returns the coordinates of the domain.
%       - A cell array of the domain of interest (string) and the buffer
%           value (scalar): {domain, buf}.
%       - An N-by-2 matrix of the domain of interest in [lon, lat] format.
%       - A GeoDomain object.
%   L - The bandwith of the spherical harmonic coefficients
%       - A scalar L, denoting the maximum degree of the SH coefficients.
%       - A 1-by-2 vector [Lmin, Lmax], denoting the minimum and maximum
%           degree of the SH coefficients.
%   timelim - The time range of GRACE data
%       A 2-by-1 datetime or timenum array specifying the starting and
%       ending time of the GRACE data used.
%       The default range uses all data available.
%   fitwhat - The functions that you would like to fit to the time series
%       data
%       The format of this array is as follows:
%       [P, T1, T2, T3, ...] where
%       - P is either 0/1/2/3 to fit up to either a mean/linear/quadratic/
%           cubic function (if increasing the power of the function reduces
%           variance enough)
%       - Tn is the period in days of a function (i.e. 365.0)
%       Any # of desired periodic functions [days] (including zero) can be
%       included.
%   product - The name of the data product
%       This is a cell array with three parts:
%       - the data center ('CSR', 'GFZ', or 'JPL'),
%       - the release level ('RL04', 'RL05', or 'RL06'),
%       - the dataproduct bandwidth.
%       The default product is {'CSR', 'RL06', 60}.
%   gia - Whether to apply the GIA correction or which model to use
%       When set to FALSE, no corrections are applied.
%       The default option is the ICE6G_D (VM5a) solution as used in the
%       CSR mascon solution.
%   Deg1Correction, C20Correction, C30Correction - Wether you want to apply
%       the degree 1, C20, or C30 corrections
%   ForceNew - Wether or not you want to force recomputation
%       The default option is false.
%   SaveData - Wether or not you want to save the output data
%       The default option is true.
%   BeQuiet - Wether or not you want to suppress output
%       The default option is false.
%
% Output arguments
%   date - Time stamps of the time series in DATETIME format
%   total - The total mass localised in the domain
%       The unit is in gigaton (10^12 kg).
%   totalerror - The uncertainty (1-sigma) associated with the total time
%       series
%   totalfit - A linear fit to the total
%   f - The parameters for the best-fit line
%       The units are in Gt/yr^p.
%       The second entry is the trend.
%   ferror - The uncertainty associated with the fit parameters
%   slept - The coefficients of each Slepian function
%   sleptSig - The least-squares fitted function for each Slepian
%       coefficient evaluated at those months in the same format
%   sleptRes - Residual time series for each Slepian coefficient
%
% Last modified by
%   2025/10/16, williameclee@arizona.edu (@williameclee)

function varargout = grace2trend(domain, varargin)
    %% Initialisation
    [domain, L, timeRange, fitwhat, product, gia, deg1cor, c20cor, ...
         c30cor, ~, beQuiet, saveData, forceNew, normalise, onlytrend, callChain] = ...
        parseinputs(domain, varargin{:});

    [outputPath, ~, outputExists] = outputfile(domain);
    corrections = {gia, deg1cor, c20cor, c30cor};

    %% Loading precomuted data
    if nargout <= 6 && ~forceNew

        if outputExists
            load(outputPath, 'G2T');

            iExp = experimentnumber( ...
                G2T, domain, L, timeRange, product, corrections);

            if ~isnan(iExp)
                date = G2T(iExp).Date;
                total = G2T(iExp).Total;
                totalerror = G2T(iExp).TotalError;
                totalfit = G2T(iExp).TotalFit;
                fitparams = G2T(iExp).FitParams;
                fitparamerrors = G2T(iExp).FitParamErrors;

                if onlytrend
                    totalfit = totalfit(:, 2);
                    fitparams = fitparams(2);
                    fitparamerrors = fitparamerrors(2);
                end

                if normalise
                    [total, totalerror, totalfit, fitparams, fitparamerrors] = ...
                        mass2weq({total, totalerror, totalfit, fitparams, fitparamerrors}, domain);
                    % total = mass2weq(total, domain);
                    % totalerror = mass2weq(totalerror, domain);
                    % totalfit = mass2weq(totalfit, domain);
                    % fitparams = mass2weq(fitparams, domain);
                    % fitparamerrors = mass2weq(fitparamerrors, domain);
                end

                varargout = ...
                    {date, total, totalerror, totalfit, fitparams, fitparamerrors};

                return
            end

        end

    end

    %% Loading data
    [slept, ~, date, ~, G, ~, V, N] = ...
        grace2slept_new(product, domain, L, ...
        "TimeRange", timeRange, "Unit", 'SD', ...
        "Deg1Correction", deg1cor, "C20Correction", c20cor, ...
        "C30Correction", c30cor, "BeQuiet", beQuiet);

    % Correct for GIA
    if gia
        [~, giaSlept] = ...
            gia2slept(date, gia, domain, L, "BeQuiet", beQuiet);
        slept = slept - giaSlept;
    end

    %% Computing trends
    [sleptSig, sleptRes, ~, ~, total, alphavarall, fitparams, ...
         fitparamerrors, totalfit] = ...
        slept2resid_new(slept, date, fitwhat, "Domain", domain, ...
        "Unit", 'year', "BeQuiet", beQuiet);

    totalerror = sqrt(alphavarall);
    total = total(:)';
    date = date(:)';

    %% Saving and output
    if saveData

        if outputExists
            load(outputPath, 'G2T');
            h = length(G2T) + 1;
        else
            h = 1;
            G2T = struct;
        end

        G2T(h).Domain = domain;
        G2T(h).L = L;
        G2T(h).TimeRange = timeRange;
        G2T(h).Product = product;
        G2T(h).Corrections = corrections;
        G2T(h).Date = date;
        G2T(h).Total = total;
        G2T(h).TotalError = totalerror;
        G2T(h).TotalFit = totalfit;
        G2T(h).FitParams = fitparams;
        G2T(h).FitParamErrors = fitparamerrors;

        save(outputPath, 'G2T');

        if ~beQuiet
            fprintf('[ULMO>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">data</a>.\n', ...
                callchaintext(callChain), outputPath, outputPath);
        end

    end

    if onlytrend
        totalfit = totalfit(:, 2);
        fitparams = fitparams(2);
        fitparamerrors = fitparamerrors(2);
    end

    if normalise
        [total, totalerror, totalfit, fitparams, fitparamerrors] = ...
            mass2weq({total, totalerror, totalfit, fitparams, fitparamerrors}, domain);
        % total = mass2weq(total, domain);
        % totalerror = mass2weq(totalerror, domain);
        % totalfit = mass2weq(totalfit, domain);
        % fitparams = mass2weq(fitparams, domain);
        % fitparamerrors = mass2weq(fitparamerrors, domain);
    end

    varargout = ...
        {date, total, totalerror, totalfit, fitparams, ...
         fitparamerrors, slept, sleptSig, sleptRes, N, G, V};

end

%% Subfunctions
function varargout = parseinputs(varargin)
    fitwhatD = [3, 365.25, 182.625, 161]; % 161: S2 tide
    ip = inputParser;
    addRequired(ip, 'Domain', ...
        @(x) ischar(x) || isstring(x) || iscell(x) || isa(x, 'GeoDomain'));
    addOptional(ip, 'L', 30, @isnumeric);
    addOptional(ip, 'TimeRange', ...
        [datetime(2003, 1, 1), datetime(2022, 12, 31)], ...
        @(x) isdatetime(x) && length(x) == 2);
    addOptional(ip, 'FitWhat', fitwhatD, ...
        @(x) isnumeric(x));
    addOptional(ip, 'DataProduct', {'CSR', 'RL06', 60}, ...
        @(x) iscell(x) && length(x) == 3);
    addParameter(ip, 'GIACorrection', 'ICE-6G_D', ...
        @(x) ischar(x) || istruefalse(x));
    addParameter(ip, 'Deg1Correction', true, @istruefalse);
    addParameter(ip, 'C20Correction', true, @istruefalse);
    addParameter(ip, 'C30Correction', true, @istruefalse);
    addParameter(ip, 'Truncation', [], @isnumeric);
    addParameter(ip, 'BeQuiet', 0.5, @(x) istruefalse(x, true));
    addParameter(ip, 'SaveData', true, @istruefalse);
    addParameter(ip, 'ForceNew', false, @istruefalse);
    addParameter(ip, 'Normalisation', true, @(x) ischar(x) || istruefalse(x));
    addParameter(ip, 'OnlyTrend', false, @(x) ischar(x) || istruefalse(x));
    addParameter(ip, 'CallChain', {mfilename}, @iscell);
    parse(ip, varargin{:});

    domain = ip.Results.Domain;
    L = ip.Results.L;
    timeRange = ip.Results.TimeRange;
    fitwhat = conddefval(ip.Results.FitWhat, fitwhatD);
    DataProduct = ip.Results.DataProduct;
    giaModel = ip.Results.GIACorrection;
    deg1cor = ip.Results.Deg1Correction;
    c20cor = ip.Results.C20Correction;
    c30cor = ip.Results.C30Correction;
    truncation = ip.Results.Truncation;
    beQuiet = uint8((double(ip.Results.BeQuiet) * 2));
    saveData = ip.Results.SaveData;
    forceNew = ip.Results.ForceNew;
    normalise = logical(ip.Results.Normalisation);
    onlytrend = logical(ip.Results.OnlyTrend);
    callChain = [ip.Results.CallChain, {mfilename}];

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain);
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{:});
    end

    if isnumeric(timeRange)
        timeRange = datetime(timeRange, "ConvertFrom", 'datenum');
    end

    if islogical(giaModel) && giaModel
        giaModel = 'Paulson07';
    end

    varargout = ...
        {domain, L, timeRange, fitwhat, DataProduct, giaModel, ...
         deg1cor, c20cor, c30cor, truncation, beQuiet, saveData, forceNew, normalise, onlytrend, callChain};
end

function [outputPath, outputFolder, outputExists] = outputfile(domain)
    outputFolder = fullfile(getenv('IFILES'), upper(mfilename));

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    outputFile = sprintf('%s-%s.mat', mfilename, domain.Domain);
    outputPath = fullfile(outputFolder, outputFile);

    outputExists = exist(outputPath, "file");

end

function iExp = ...
        experimentnumber(G2T, domain, L, timeRange, product, corrections)
    iExp = nan;

    for iiExp = 1:length(G2T)

        if isequal(G2T(iiExp).Domain, domain) && ...
                isequal(G2T(iiExp).L, L) && ...
                isequal(G2T(iiExp).TimeRange, timeRange) && ...
                isequal(G2T(iiExp).Product, product) && ...
                isequal(G2T(iiExp).Corrections, corrections)
            iExp = iiExp;
            break
        end

    end

end
