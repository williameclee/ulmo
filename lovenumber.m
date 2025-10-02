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
%
% Output arguments
%   ln - Love numbers
%
% Authored by
%   2024/11/20, williameclee@arizona.edu (@williameclee)
%
% Last modified by
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
    parse(ip, l, varargin{:});
    l = ip.Results.l;
    type = ip.Results.type;
    frame = ip.Results.frame;

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
    try
        ln = load(fullfile(getenv('IFILES'), 'lovenumbers.mat'), type);
    catch
        error('Check your love number file at %s', fullfile('IFILES', 'lovenumbers.mat'));
    end

    ln = ln.(type);

    if strcmpi(frame, 'CF')
        % from Blewitt, 2003, JGR
        if strcmpi(type, 'loadingverticaldisplacement')
            ln(2) = -0.269;
        elseif strcmpi(type, 'loadinggravitationalpotential')
            ln(2) = 0.021;
        elseif strcmpi(type, 'loadinghorizontaldisplacement')
            ln(2) = 0.134;
        end

    end

    ln = ln(l + 1); % l starts at 0
    ln = reshape(ln, size(l)); 

end
