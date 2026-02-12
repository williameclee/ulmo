%% LOVENUMBER
% Returns the Love number for a given degree l.
% The love numbers are obtained from the ISSM repository (https://github.com/ISSMteam/ISSM.git).
%
% Input arguments
%   l - Degree of the Love number
%   type - Type of the Love number
%       'loadingverticaldisplacement': Loading Love number for vertical displacements
%       'loadinggravitationalpotential': Loading Love number for gravitational potentials
%       'loadinghorizontaldisplacement': Loading Love number for horizontal displacements
%       'tidalverticaldisplacement': Tidal Love number for vertical displacements
%       'tidalgravitationalpotential': Tidal Love number for gravitational potentials
%       'tidalhorizontaldisplacement': Tidal Love number for horizontal displacements
%   frame - Reference frame
%       'CM': Centre of mass frame.
%       'CF': Centre of figure frame.
%       The default frame is the CF frame.
%   source - Source of the Love numbers
%       'ISSM': Love numbers from the ISSM repository.
%       'Wahr': Love numbers used in the slepian_delta package.
%       'ALMA3 <model>': Love numbers from ALMA3 files for the specified Earth model (e.g. 'ALMA3 PREM').
%           As of right now, the love numbers have to be manually computed using the command line and not automated.
%       The default source is 'ISSM'.
%
% Output arguments
%   ln - Love numbers
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2026/02/12, williameclee@arizona.edu (@williameclee)
%     - Added support for ALMA3 love numbers
%   2025/08/03, williameclee@arizona.edu (@williameclee)

function ln = lovenumber(l, varargin)
    %% Initialisation
    ip = inputParser;
    addRequired(ip, 'l', @isnumeric);
    addOptional(ip, 'type', 'loadinggravitationalpotential', ...
        @(x) ischar(validatestring(lower(x), ...
        {'loadingverticaldisplacement', 'loadinggravitationalpotential', ...
         'loadinghorizontaldisplacement', 'tidalverticaldisplacement', ...
         'tidalgravitationalpotential', 'tidalhorizontaldisplacement', ...
         'LLN geoid', 'LLN VLM', 'TLN geoid', 'TLN VLM', 'LLN HLM', 'TLN HLM'})));
    addOptional(ip, 'frame', 'CF', ...
        @(x) ischar(validatestring(upper(x), {'CM', 'CF'})));
    addParameter(ip, "Source", "ISSM", ...
        @(x) (ischar(x) || isstring(x)) && ...
        (startsWith(upper(x), "ALMA3") || ischar(validatestring(upper(x), {'ISSM', 'WAHR'}))));
    parse(ip, l, varargin{:});
    l = ip.Results.l;
    type = ip.Results.type;
    frame = ip.Results.frame;
    source = ip.Results.Source;

    if max(l) > 10000
        error('Love number must be less than or equal to 10000');
    end

    switch type
        case 'LLN geoid'
            type = 'loadinggravitationalpotential';
        case 'LLN VLM'
            type = 'loadingverticaldisplacement';
        case 'LLN HLM'
            type = 'loadinghorizontaldisplacement';
        case 'TLN geoid'
            type = 'tidalgravitationalpotential';
        case 'TLN VLM'
            type = 'tidalverticaldisplacement';
        case 'TLN HLM'
            type = 'tidalhorizontaldisplacement';
    end

    %% Loading Love numbers
    switch upper(char(source))
        case 'ISSM'
            ln = lovenumber_issm(l, type);
        case 'WAHR'
            % Load default slepian_delta love numbers
            if ~strcmpi(frame, 'CF')
                error('Wahr et al. (1998) love numbers are only available in the CF frame, but %s was requested.', ...
                    frame);
            elseif ~strcmpi(type, 'loadinggravitationalpotential')
                error('Wahr et al. (1998) love numbers are only available for loading gravitational potentials, but %s was requested.', ...
                    type);
            end

            ln = lovenums('Wahr', l);
            ln = ln(:, 2);
        otherwise

            if startsWith(upper(source), "ALMA3")
                ln = lovenumber_alma3(l, type, source);
            else
                error('Unknown love number source: %s', source);
            end

    end

    if strcmpi(frame, 'CF') && any(l == 1)
        l1id = find(l == 1);
        % from Blewitt, 2003, JGR
        if strcmpi(type, 'loadingverticaldisplacement')
            ln(l1id) = -0.269;
        elseif strcmpi(type, 'loadinggravitationalpotential')
            ln(l1id) = 0.021;
        elseif strcmpi(type, 'loadinghorizontaldisplacement')
            ln(l1id) = 0.134;
        end

    end

    ln = reshape(ln, size(l));

end

function ln = lovenumber_issm(l, type)

    if max(l) > 10000
        error('Love number must be less than or equal to 10000');
    end

    ln = load(fullfile(getenv('IFILES'), 'LOVENUMS', 'lovenumbers-ISSM.mat'), type);

    ln = ln.(type);
    ln = ln(l + 1); % l starts at 0
end

function ln = lovenumber_alma3(l, type, model)

    if startsWith(upper(model), "ALMA3")
        % Remove the prefix
        modelParsed = regexprep(modelUpper, "^ALMA3\s*", "");

        if strlength(modelParsed) == 0
            error('Model name must be specified after the ALMA3 prefix (e.g., "ALMA3 ICE").');
        end

        model = modelParsed;
    end

    dataFolder = fullfile(getenv('IFILES'), 'LOVENUMS');

    switch type
        case 'loadinggravitationalpotential'
            dataFile = sprintf("k_%sLLN-*.dat", model);
        case 'loadingverticaldisplacement'
            dataFile = sprintf("h_%sLLN-*.dat", model);
        case 'loadinghorizontaldisplacement'
            dataFile = sprintf("l_%sLLN-*.dat", model);
        case 'tidalgravitationalpotential'
            dataFile = sprintf("k_%sTLN-*.dat", model);
        case 'tidalverticaldisplacement'
            dataFile = sprintf("h_%sTLN-*.dat", model);
        case 'tidalhorizontaldisplacement'
            dataFile = sprintf("l_%sTLN-*.dat", model);
    end

    files = dir(fullfile(dataFolder, dataFile));

    if isempty(files)
        error('No ALMA3 Love number file found for model %s and type %s.', model, type);
    elseif length(files) > 1
        % Choose the last file when sorted alphabetically
        [~, idx] = sort({files.name});
        files = files(idx);
        warning('Multiple ALMA3 Love number files found for model %s and type %s. Using %s.', ...
            model, type, files(end).name);
        files = files(end);
    end

    dataPath = fullfile(dataFolder, files(1).name);

    ln = readmatrix(dataPath, "NumHeaderLines", 5);
    ln = interp1(ln(:, 1), ln(:, 2), l, 'linear', 'extrap');
    % Replace degree 0 & 1 love number with 0
    ln(l == 0) = 0;
    ln(l == 1) = 0;
end
