%% CONVERTGRAVITY
% Converts between gravity, potential, and surface mass density spherical
% harmonic coefficients.
%
% Notes
%	The functionalities are being added as needed. not all conversions are
%   supported.
%
% Authored by
%   2025/05/22, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/08/`3, williameclee@arizona.edu (@williameclee)

function output = convertgravity(varargin)
    %% Initialisation
    [input, inputUnit, outputUnit, gravityParam, equatorRadius, ...
         inputFmt, timeDim, L, avgDensity] = parseinputs(varargin{:});

    if strcmpi(inputUnit, outputUnit)
        output = input;
        return
    end

    switch strjoin({inputUnit, outputUnit})
        case strjoin({'GRAV', 'POT'})
            output = grav2pot(input, equatorRadius, inputFmt, timeDim);
        case strjoin({'GRAV', 'SD'})
            output = grav2pot(input, equatorRadius, inputFmt, timeDim);
            output = pot2sd(output, L, avgDensity, inputFmt, timeDim);
        case strjoin({'POT', 'GRAV'})
            output = pot2grav(input, equatorRadius, inputFmt, timeDim);
        case strjoin({'POT', 'SD'})
            output = pot2sd(input, L, avgDensity, inputFmt, timeDim);
        case strjoin({'SD', 'POT'})
            output = sd2pot(input, L, avgDensity, inputFmt, timeDim);
        case strjoin({'SD', 'GRAV'})
            output = sd2pot(input, L, avgDensity, inputFmt, timeDim);
            output = pot2grav(output, equatorRadius, inputFmt, timeDim);
    end

end

%% Subfunctions
function pot = grav2pot(grav, equatorRadius, inputFmt, timeDim)

    switch inputFmt
        case 'lmcosi'

            if ismatrix(grav)
                pot = grav;
                pot(:, 3:4) = grav(:, 3:4) * equatorRadius;
                return
            end

            switch timeDim
                case 'timefirst'
                    pot = grav;
                    pot(:, :, 3:4) = grav(:, :, 3:4) * equatorRadius;
                case 'traditional'
                    pot = grav * equatorRadius;
                    pot(:, 3:4, :) = grav(:, 3:4, :) * equatorRadius;
            end

        case {'cosi', 'L'}
            pot = grav * equatorRadius;
    end

end

function grav = pot2grav(pot, equatorRadius, inputFmt, timeDim)

    switch inputFmt
        case 'lmcosi'

            if ismatrix(pot)
                grav = pot;
                grav(:, 3:4) = pot(:, 3:4) / equatorRadius;
                return
            end

            switch timeDim
                case 'timefirst'
                    grav = pot;
                    grav(:, :, 3:4) = pot(:, :, 3:4) / equatorRadius;
                case 'traditional'
                    grav = pot * equatorRadius;
                    grav(:, 3:4, :) = pot(:, 3:4, :) / equatorRadius;
            end

        case {'cosi', 'L'}
            grav = pot / equatorRadius;
    end

end

function sd = pot2sd(pot, degree, avgDensity, inputFmt, timeDim)

    switch inputFmt
        case 'lmcosi'

            if ismatrix(pot)
                degree = pot(:, 1);
            else

                switch timeDim
                    case 'timefirst'
                        degree = pot(:, :, 1);
                    case 'traditional'
                        degree = pot(:, 1, :);
                end

            end

            lovenum = lovenumber(degree, 'loadinggravitationalpotential');
            % lovenum = lovenums('Wahr', degree);
            % lovenum = reshape(lovenum(:, 2), size(degree));

            convFactor = avgDensity / 3 .* (2 * degree + 1) ./ (1 + lovenum);

            sd = pot;

            if ismatrix(pot)
                sd(:, 3:4) = pot(:, 3:4) .* convFactor;
            else

                switch timeDim
                    case 'timefirst'
                        sd(:, :, 3:4) = pot(:, :, 3:4) .* convFactor;
                    case 'traditional'
                        sd(:, 3:4, :) = pot(:, 3:4, :) .* convFactor;
                end

            end

        case 'cosi'
            error('Not implemented yet for COSI format');
        case 'L'
            lovenum = lovenumber(degree, 'loadinggravitationalpotential');
            sd = avgDensity / 3 * (2 .* degree + 1) ./ (1 + lovenum) .* pot;
    end

end

function pot = sd2pot(sd, degree, avgDensity, inputFmt, timeDim)

    switch inputFmt
        case 'lmcosi'

            if ismatrix(sd)
                degree = sd(:, 1);
            else

                switch timeDim
                    case 'timefirst'
                        degree = sd(:, :, 1);
                    case 'traditional'
                        degree = sd(:, 1, :);
                end

            end

            lovenum = lovenumber(degree, 'loadinggravitationalpotential');
            % lovenum = lovenums('Wahr', degree);
            % lovenum = reshape(lovenum(:, 2), size(degree));

            convFactor = 3 / avgDensity .* (1 + lovenum) ./ (2 * degree + 1);

            pot = sd;

            if ismatrix(sd)
                pot(:, 3:4) = sd(:, 3:4) .* convFactor;
            else

                switch timeDim
                    case 'timefirst'
                        pot(:, :, 3:4) = sd(:, :, 3:4) .* convFactor;
                    case 'traditional'
                        pot(:, 3:4, :) = sd(:, 3:4, :) .* convFactor;
                end

            end

        case 'cosi'
        case 'L'
            lovenum = lovenumber(degree, 'loadinggravitationalpotential');
            % lovenum = lovenums('Wahr', degree);
            % lovenum = lovenum(:, 2);
            convFactor = 3 / avgDensity .* (1 + lovenum) ./ (2 * degree + 1);
            pot = sd .* convFactor;
    end

end

% Parse input arguments
function varargout = parseinputs(varargin)
    ip = inputParser;
    addRequired(ip, 'Plm', @isnumeric);
    addRequired(ip, 'InputUnit', ...
        @(x) ischar(validatestring(x, {'GRAV', 'POT', 'SD'})));
    addRequired(ip, 'OutputUnit', ...
        @(x) ischar(validatestring(x, {'GRAV', 'POT', 'SD'})));
    addOptional(ip, 'GravityParam', fralmanac('GM_EGM96', 'Earth'), ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addOptional(ip, 'EquatorRadius', fralmanac('a_EGM96', 'Earth'), ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'InputFormat', 'lmcosi', ...
        @(x) ischar(validatestring(x, {'lmcosi', 'cosi', 'L'})));
    addParameter(ip, 'TimeDim', 'timefirst', ...
        @(x) ischar(validatestring(x, {'timefirst', 'traditional'})));
    addParameter(ip, 'L', [], ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(ip, 'EarthDensity', 5517, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(ip, varargin{:});

    input = ip.Results.Plm;
    inputUnit = upper(ip.Results.InputUnit);
    outputUnit = upper(ip.Results.OutputUnit);
    gravityParam = ip.Results.GravityParam;
    equatorRadius = ip.Results.EquatorRadius;
    inputFmt = ip.Results.InputFormat;
    timeDim = ip.Results.TimeDim;
    L = round(ip.Results.L);
    avgDensity = ip.Results.EarthDensity;

    if strcmpi(inputFmt, 'L') && isempty(L)
        error(sprintf('%s:InvalidInput:DegreeNotSpecified', mfilename), ...
        'Degree not specified for L format');
    end

    varargout = ...
        {input, inputUnit, outputUnit, gravityParam, equatorRadius, ...
         inputFmt, timeDim, L, avgDensity};

end
