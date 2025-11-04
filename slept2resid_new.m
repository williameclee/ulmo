%% SLEPT2RESID
% Takes a time series of Slepian coefficients and fits a desired
% combination of functions (e.g. secular, annual, semiannual, etc.) to
% each coefficient in order to separate 'signal' from residual.
% You can choose to fit either a mean, linear, quadratic, or cubic fuction
% as your 'base' function to each Slepian coefficient, by using 'FITWHAT'.
% In these cases of higher functions, they are only used if they pass an
% F-test for variance reduction. For example, a cubic term is only used if
% it is significantly better than a quadratic function.
%
% If you also provide the Slepian functions (CC) and region (domain) for
% this localisation, then we assume that you want the fitting to the
% integrated functions, i.e. if you have surface density, then using the
% integrated Slepian function would give you total mass change in a region.
% In addition, there will be a fitting of the combination (total) of the
% Slepian functions up to the Shannon number for this localisation.
%
% Syntax
%   [sleptSig, sleptRes, ftests] = slept2resid(slept, dates, fitwhat, ...
%       givenerrors, specialterms, CC, domain)
%   [__, extravalues, total, alphavarall, totalparams, ...
%       totalparamerrors, totalfit, eigfunInt, alphavar] = slept2resid(__)
%
% Input arguments
%   sleptCoeffs - The time series of Slepian coefficients
%       This should be a two dimensional matrix (not a cell array), where
%       the first dimension is time, and the second dimension are Slepian
%       coefficients sorted by global eigenvalue.
%   dates - An array of dates corresponding to the SLEPT timeseries
%       It is assumed that DATES matches SLEPT. If DATES is longer than
%       SLEPT then it is assumed that these are extra dates that you would
%       like SLEPTSIG to be evaluated at. This is useful if you want to
%       interpolate to one of the GRACE missing months that is important,
%       such as Jan 2011.
%       The input can be in DATENUM or DATETIME format.
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
%   givenerrors - Given errors, if available
%       In this case a weighted inversion is performed. It should be the
%       same dimensions of SLEPT.
%   specialterms - A cell array such as {2, 'periodic', 1460}
%       At the moment this is pretty specific to our needs, but could be
%       expanded later.
%   CC - The localisation Slepian functions
%       This is can be the cell array obtained from GRACE2SLEPT, or the
%       localisation matrix G from GLMALPHA.
%   domain - The region of interest
%   J - Number of largest eigenfunctions in which to expand
%       The default value is the rounded Shannon number.
%
% Output arguments
%   sleptSig - The least-squares fitted function for each Slepian
%       coefficient evaluated at those months in the same format
%   sleptRes - Residual time series for each Slepian coefficients
%       They are ordered as they were given, presumably by eigenvalue
%       [nmonths x (Lwindow+1)^2]
%   ftests - A matrix, such as [0 1 1] for each Slepian coefficient, on
%       whether the fits you requested passed an F-test for significance.
%   extravalues - These are the values of SLEPTSIG evaluated at your extra
%       dates which were tacked onto DATES
%   total - The time series of the combined mass change from J Slepian
%       functions (i.e. the combined data points)
%       The unit is gigaton (Gt, 10^12 kg).
%   alphavarall - The time averaged variance on each data point (from error
%       propogation from each individual function)
%       These values are the same for every point.
%   totalparams - The parameters of the fits to the total
%       This is a 4-by-j matrix where j are the fits you wanted. Zeros fill
%       out the unused parameters. For example, if you want just a 2nd
%       order polynomial, you will get
%       [intercept intercept; slope slope; 0 quadratic; 0 0]
%       The unit is Gt/day^p.
%   totalparamerrors - The 95% confidence value on this slope
%       This is computed using the critical value for a t-distribution.
%        At the moment, totalparams and totalparamerrors and just the
%       values for a linear fit. Sometime later maybe change this to be
%       potentially a quadratic fit as well.
%       The unit is Gt/yr.
%   totalfit - Datapoints for the best fit line, so you can plot it, and
%       datapoints for the confidence intervals of this linear fit, so you
%       can plot them.
%       To plot these you should use totalfit(:, 2) + totalfit(:, 3) and
%       totalfit(:, 2) - totalfit(:, 3).
%   eigfunInt - This is a vector of the integrals of the Slepian functions,
%       up to J
%       With this we can multiply by the data or SLEPTSIG to get the
%       changes in total mass over time for each function.
%       The unit is Gt.
%   alphavar - The variances on the individual alpha function time series
%       (integrals)
%       This is calculated by making a covariance matrix from SLEPTRES and
%       then doing the matrix multiplication with the eigenfunction
%       integrals.
%
% Last modified by
%   2025/11/03, williameclee@arizona.edu (@williameclee)
%   2012/06/26, charig@princeton.edu (@harig00)

function varargout = slept2resid_new(varargin)
    %% Initialisation
    % Parse inputs
    [slept, date, fitwhat, givenerrors, spTerms, domain, truncation, ...
         unit, doNormalisation, transFlag, beQuiet, callChain] = parseinputs(varargin{:});

    % Format dates

    if length(date) == size(slept, 1) % do nothing
        dateExtra = [];
        moredates = false;
    elseif length(date) > size(slept, 1) % pull off extra dates
        dateExtra = date((size(slept, 1) + 1):end);
        date = date(1:size(slept, 1));
        moredates = true;
    else
        error('Is dates shorter than slept?')
    end

    nMonth = length(date);

    % We will do a scaling/normalisation to improve the solution
    dateMean = mean(date);
    dateStd = std(date);
    dateNml = (date - dateMean) / dateStd; % normalisation
    dateExtraNml = (dateExtra - dateMean) / dateStd;

    % Periodic components
    nOmega = length(fitwhat(2:end));

    if nOmega > 0
        % The (not angular) frequencies being fitted in [1/days]
        freq = 1 ./ fitwhat(2:end);
        % Rescale these to our new xprime
        freq = freq * dateStd;
    else
        freq = [];
    end

    % Slepian coefficients
    nSleps = size(slept, 2);
    % Preallocate the residuals and the evaluated fitted function set
    sleptSig = zeros(size(slept));
    sleptRes = zeros(size(slept));
    extravalues = zeros([length(dateExtra), nSleps]);
    ftests = zeros([nSleps, 3]);
    % Figure out the orders and degrees of this setup
    % BUT orders and degrees have lost their meaning since slept should be
    % ordered by eigenvalue

    %% G matrix assembly
    % We will have the same number of G matrices as order of polynomial
    % fit. These matrices are smallish, so make all 3 regardless of whether
    % you want them all.
    G1 = []; % For line fits
    G2 = []; % For quadratic fits
    G3 = []; % For cubic fits
    % Mean term
    if fitwhat(1) >= 0
        G1 = [G1, ones(size(dateNml'))];
        G2 = [G2, ones(size(dateNml'))];
        G3 = [G3, ones(size(dateNml'))];
    end

    % Secular term
    if fitwhat(1) >= 1
        G1 = [G1, dateNml'];
        G2 = [G2, dateNml'];
        G3 = [G3, dateNml'];
    end

    % Quadratic term
    if fitwhat(1) >= 2
        G2 = [G2, dateNml' .^ 2];
        G3 = [G3, dateNml' .^ 2];
    end

    % Cubic term
    if fitwhat(1) == 3
        G3 = [G3, dateNml' .^ 3];
    end

    % Periodic terms
    if nOmega > 0
        % Angular frequency in radians/(rescaled day) of the periodic terms
        phase = repmat(freq, nMonth, 1) * 2 * pi ... % change into angular
            .* repmat(dateNml', 1, nOmega);
        G1 = [G1, cos(phase), sin(phase)];
        G2 = [G2, cos(phase), sin(phase)];
        G3 = [G3, cos(phase), sin(phase)];

        % Create our specialterms G if we have it. At the moment this is
        % just an additional periodic function, but in the future we could
        % add something else.
        if ~isnan(spTerms{1})
            % Strip off the previous periodic terms here and REPLACE with
            % freqSpec
            spFreq = [freq, dateStd / spTerms{3}];
            spPhase = repmat(spFreq, nMonth, 1) * 2 * pi ... % angular
                .* repmat(dateNml', 1, length(spFreq));
            GSpec1 = [G1(:, 1:end - 2 * nOmega), ...
                          cos(spPhase), sin(spPhase)];
            GSpec2 = [G2(:, 1:end - 2 * nOmega), ...
                          cos(spPhase), sin(spPhase)];
            GSpec3 = [G3(:, 1:end - 2 * nOmega), ...
                          cos(spPhase), sin(spPhase)];
        else
            spFreq = [];
            spPhase = [];
            GSpec1 = [];
            GSpec2 = [];
            GSpec3 = [];
        end

    else
        phase = [];
        spFreq = [];
        spPhase = [];
        GSpec1 = [];
        GSpec2 = [];
        GSpec3 = [];
    end

    %% Compute the fits
    isSpecialterm = 1:nSleps == spTerms{1};
    % Since each Slepian coefficient has different errors, each will have a
    % different weighting matrix.  Thus we loop over the coefficients.
    parfor iSlep = 1:nSleps
        [sleptSig(:, iSlep), sleptRes(:, iSlep), ftests(iSlep, :), ...
             extravalues(:, iSlep)] = slept2resid_fitslept( ...
            slept(:, iSlep), fitwhat, givenerrors(:, iSlep), ...
            isSpecialterm(iSlep), freq, phase, spFreq, spPhase, ...
            dateNml, dateExtraNml, moredates, ...
            G1, G2, G3, GSpec1, GSpec2, GSpec3, nMonth);

    end

    %% Returning requested output
    if transFlag
        sleptSigOut = sleptSig';
        sleptResOut = sleptRes';
        extravaluesOut = extravalues';
    else
        sleptSigOut = sleptSig;
        sleptResOut = sleptRes;
        extravaluesOut = extravalues;
    end

    varargout = {sleptSigOut, sleptResOut, ftests, extravaluesOut};

    % Total combined fitting

    % If we have the parameters for this localization, and we requested the
    % total fit, then let's do that.

    if nargout < 5
        return
    end

    % Get the residual covariance
    [Cab] = slepresid2cov(sleptRes);

    % Calculate the bandwdith for this basis
    L = sqrt(size(slept, 2)) - 1;
    % This should be an integer
    if (floor(L) ~= L)
        error('Something fishy about your L');
    end

    if isa(domain, 'GeoDomain')
        XY = domain.Lonlat;
    else
        % Coordinates or a string, either works
        XY = domain;
    end

    % Calculate the Shannon number for this basis
    truncation = conddefval(truncation, ceil((L + 1) ^ 2 * spharea(XY)));

    % Make the coefficients with reference to some mean
    % If they already are, then this won't matter
    sleptdelta = slept(1:nMonth, :) ...
        - repmat(mean(slept(1:nMonth, :), 1), nMonth, 1);

    % COMBINE

    % We want to take the Slepian functions and combine them to get total mass.
    % For signal, this means integrating the functions and adding them.  For
    % the error, this means using the error propogation equation, where we
    % compute (int)*(covar)*(int)'.  Since the slepcoffs are constants that
    % just come forward, we can do the integration of the eigenfunctions
    % first (and once for each function), then multiply by slepcoffs to
    % get the monthly values.  This is much faster.

    if isa(domain, 'GeoDomain') || ismatrix(domain)
        G = glmalpha_new(domain, L, "BeQuiet", beQuiet, "CallChain", callChain);

        try
            eigfunINT = integratebasis_new( ...
                G, domain, truncation, "BeQuiet", beQuiet, "CallChain", callChain);
        catch
            eigfunINT = integratebasis_new( ...
                G, domain, truncation, "BeQuiet", beQuiet, "ForceNew", true, " CallChain", callChain);
        end

    else
        G = glmalpha_new(domain, L, phi, theta, omega, "BeQuiet", beQuiet, "CallChain", callChain);
        eigfunINT = integratebasis_new( ...
            G, domain, truncation, phi, theta, "BeQuiet", beQuiet, "CallChain", callChain);
    end

    % Since Int should have units of (fn * m^2), need to go from fractional
    % sphere area to real area.  If the fn is surface density, this output is
    % in kilograms.  Then change the units from kg to Gt in METRIC tons
    if doNormalisation
        eigfunINT = eigfunINT * 4 * pi * 6370000 ^ 2/1e12;
    end

    functionintegrals = eigfunINT;

    % Now multiply by the appropriate slepcoffs to get the months
    % This becomes alpha by months
    %functimeseries=repmat(eigfunINT',1,nmonths).*sleptdelta(:,1:N)';
    %functimeseries = sleptdelta(:,1:N)';

    % Here do the total sum of the data
    eigfunINT = eigfunINT(:);

    try
        eigfunINT = eigfunINT(1:truncation);
    catch
        disp(size(eigfunINT))
        disp(truncation)
        error('truncation is too large')
    end

    total = eigfunINT' * sleptdelta(:, 1:truncation)';

    % Get the error
    thevars = diag(Cab(1:truncation, 1:truncation))';
    alphavar = eigfunINT .^ 2 .* thevars;
    % Now the combined error with covariance
    alphavarall = eigfunINT' * Cab(1:truncation, 1:truncation) * eigfunINT;

    % FITTING

    % We have uniform estimated error, which will be different than the polyfit
    % estimated residuals because ours account for sinusoidal signals.  So
    % pass the new error to our function for replacement, so
    % that the fitting confidence intervals reflect that

    [fit, delta, totalparams, paramerrors] = ...
        timeseriesfit([date', total'], alphavarall, 1, 1);

    % Make a matrix for the line, and 95% confidence in the fit
    totalfit = [date', fit, delta];

    % Make the error valid for a year
    totalparamerrors = paramerrors * days(years(1));

    if strcmp(unit, 'year')
        totalparams = totalparams * days(years(1));
    end

    % Collect the expanded output
    varargout = ...
        {sleptSigOut, sleptResOut, ftests, extravaluesOut, ...
         total, alphavarall, totalparams, totalparamerrors, totalfit, ...
         functionintegrals, alphavar};

end

%% Subfunctions
function varargout = parseinputs(varargin)
    sleptD = ...
        'grace2slept("CSR", "greenland", 0.5, 60, [], [], [], [], "SD")';
    dateD = datetime(2004, 1:12, 1);
    fitwhatD = [3, 365.0];
    givenerrorsD = [];
    specialtermsD = {NaN};
    domainD = [];
    ND = [];

    ip = inputParser;
    addOptional(ip, 'slept', sleptD);
    addOptional(ip, 'date', dateD, @(x) isnumeric(x) || isdatetime(x));
    addOptional(ip, 'Fit', fitwhatD);
    addOptional(ip, 'givenerrors', givenerrorsD);
    addOptional(ip, 'specialterms', specialtermsD);
    addOptional(ip, 'Domain', domainD, ...
        @(x) ischar(x) || isstring(x) || iscell(x) ...
        || isa(x, 'GeoDomain') || isnumeric(x) || isempty(x));
    addOptional(ip, 'N', ND);
    addOptional(ip, 'MoreDomainSpecs', {});
    addParameter(ip, 'Unit', 'original', @(x) ischar(x) ...
        && ismember(x, {'original', 'year'}));
    addParameter(ip, 'BeQuiet', false, @(x) islogical(x) || isnumeric(x));
    addParameter(ip, 'Normalisation', true, @istruefalse);
    addParameter(ip, 'CallChain', {}, @iscell);
    parse(ip, varargin{:});
    slept = conddefval(ip.Results.slept, sleptD);
    date = conddefval(ip.Results.date, dateD);
    fitwhat = conddefval(ip.Results.Fit, fitwhatD);
    givenerrors = conddefval(ip.Results.givenerrors, givenerrorsD);
    specialterms = conddefval(ip.Results.specialterms, specialtermsD);
    domain = conddefval(ip.Results.Domain, domainD);
    N = conddefval(ip.Results.N, ND);
    moreDomainSpecs = ip.Results.MoreDomainSpecs;
    unit = ip.Results.Unit;
    doNormalisation = logical(ip.Results.Normalisation);
    beQuiet = logical(ip.Results.BeQuiet);
    callChain = [ip.Results.CallChain, {mfilename}];

    if isdatetime(date)
        date = datenum(date); %#ok<DATNM>
    end

    date = date(:)';

    transFlag = false;

    if size(slept, 1) ~= length(date) && size(slept, 2) == length(date)
        slept = slept';
        transFlag = true;
    end

    if isduration(fitwhat)
        % Convert durations to datetimes
        fitwhat = days(fitwhat);
    end

    % Change the domain to a GeoDomain object if appropriate
    if ischar(domain) || isstring(domain) && exist(domain, "file")
        domain = GeoDomain(domain, moreDomainSpecs{:});
    elseif iscell(domain) && length(domain) == 2
        domain = ...
            GeoDomain(domain{1}, "Buffer", domain{2}, moreDomainSpecs{:});
    elseif iscell(domain) && length(domain) >= 3
        domain = GeoDomain(domain{1}, moreDomainSpecs{:}, domain{2:end});
    end

    if strcmp(slept, sleptD)
        % Evaluate the specified expression
        [slept, ~, date, domain, ~, ~] = eval(slept);
    else
        givenerrors = conddefval(givenerrors, ones(size(slept)));
    end

    varargout = {slept, date, fitwhat, givenerrors, specialterms, domain, N, unit, doNormalisation, transFlag, beQuiet, callChain};
end
