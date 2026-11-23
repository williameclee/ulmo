%% CONVERTGRAVITY - Converts between gravity, potential, and surface mass density spherical harmonic coefficients.
%
% Notes
%	The functionalities are being added as needed. not all conversions are
%   supported.
%
% Authored by
%   2025/05/22, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%   2026/11/23, En-Chi Lee (williameclee@gmail.com)

function output = convertgravity(input, inUnit, outUnit, options)

    arguments (Input)
        input {mustBeNumeric}
        inUnit (1, 1) string {mustBeMember(inUnit, ["GRAV", "POT", "SD"])}
        outUnit (1, 1) string {mustBeMember(outUnit, ["GRAV", "POT", "SD"])}
        options.GravityParam (1, 1) {mustBeNumeric, mustBePositive} = fralmanac('GM_EGM96', 'Earth')
        options.EquatorRadius (1, 1) {mustBeNumeric, mustBePositive} = fralmanac('a_EGM96', 'Earth')
        options.InputFormat (1, 1) string ...
            {mustBeMember(options.InputFormat, ["lmcosi", "cosi", "L"])} = "lmcosi"
        options.TimeDim (1, 1) string ...
            {mustBeMember(options.TimeDim, ["timefirst", "traditional"])} = "timefirst"
        options.L {mustBeNumeric, mustBeScalarOrEmpty, mustBePositive} = []
        options.EarthDensity (1, 1) {mustBeNumeric, mustBePositive} = 5517
        options.LoveNumSource (1, 1) string = "ISSM"
    end

    L = round(options.L);

    if options.InputFormat == "L" && isempty(L)
        error('ULMO:convertgravity:InvalidInput:DegreeNotSpecified', ...
        'Degree not specified for L format');
    end

    if inUnit == outUnit
        output = input;
        return
    end

    switch inUnit + " " + outUnit
        case "GRAV POT"
            output = grav2pot(input, options.EquatorRadius, options.InputFormat, options.TimeDim);
        case "GRAV SD"
            output = grav2pot(input, options.EquatorRadius, options.InputFormat, options.TimeDim);
            output = pot2sd(output, L, options.EarthDensity, options.InputFormat, options.TimeDim, options.LoveNumSource);
        case "POT GRAV"
            output = pot2grav(input, options.EquatorRadius, options.InputFormat, options.TimeDim);
        case "POT SD"
            output = pot2sd(input, L, options.EarthDensity, options.InputFormat, options.TimeDim, options.LoveNumSource);
        case "SD POT"
            output = sd2pot(input, L, options.EarthDensity, options.InputFormat, options.TimeDim, options.LoveNumSource);
        case "SD GRAV"
            output = sd2pot(input, L, options.EarthDensity, options.InputFormat, options.TimeDim, options.LoveNumSource);
            output = pot2grav(output, options.EquatorRadius, options.InputFormat, options.TimeDim);
    end

end

%% Subfunctions
function pot = grav2pot(grav, eqRad, inputFmt, timeDim)

    switch inputFmt
        case 'lmcosi'

            if ismatrix(grav)
                pot = grav;
                pot(:, 3:4) = grav(:, 3:4) * eqRad;
                return
            end

            switch timeDim
                case 'timefirst'
                    pot = grav;
                    pot(:, :, 3:4) = grav(:, :, 3:4) * eqRad;
                case 'traditional'
                    pot = grav * eqRad;
                    pot(:, 3:4, :) = grav(:, 3:4, :) * eqRad;
            end

        case {'cosi', 'L'}
            pot = grav * eqRad;
    end

end

function grav = pot2grav(pot, eqRad, inputFmt, timeDim)

    switch inputFmt
        case 'lmcosi'

            if ismatrix(pot)
                grav = pot;
                grav(:, 3:4) = pot(:, 3:4) / eqRad;
                return
            end

            switch timeDim
                case 'timefirst'
                    grav = pot;
                    grav(:, :, 3:4) = pot(:, :, 3:4) / eqRad;
                case 'traditional'
                    grav = pot * eqRad;
                    grav(:, 3:4, :) = pot(:, 3:4, :) / eqRad;
            end

        case {'cosi', 'L'}
            grav = pot / eqRad;
    end

end

function sd = pot2sd(pot, degree, avgDensity, inputFmt, timeDim, lovenumSrc)

    switch inputFmt
        case 'lmcosi'

            if ismatrix(pot)
                degree = pot(:, 1);
            else

                switch timeDim
                    case 'timefirst'
                        degree = pot(1, :, 1);
                    case 'traditional'
                        degree = pot(:, 1, 1);
                end

            end

            lovenum = lovenumber(degree, 'loadinggravitationalpotential', "Source", lovenumSrc);

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
            lovenum = lovenumber(degree, 'loadinggravitationalpotential', "Source", lovenumSrc);
            sd = avgDensity / 3 * (2 .* degree + 1) ./ (1 + lovenum) .* pot;
    end

end

function pot = sd2pot(sd, degree, avgDensity, inputFmt, timeDim, lovenumSrc)

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

            lovenum = lovenumber(degree, 'loadinggravitationalpotential', "Source", lovenumSrc);

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
            lovenum = lovenumber(degree, 'loadinggravitationalpotential', "Source", lovenumSrc);
            convFactor = 3 / avgDensity .* (1 + lovenum) ./ (2 * degree + 1);
            pot = sd .* convFactor;
    end

end
