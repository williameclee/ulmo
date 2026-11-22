%% RESOLVEEARTHQUAKES - Resolves event names and mixed epicentre inputs.
%
% Syntax
%   [epicentres, names, expectedDates] = resolveearthquakes(events)
%
% Input arguments
%   events - Scalar name, string/name-cell vector, numeric (nEvents, 2)
%       [longitude, latitude] matrix, or cell vector mixing scalar names
%       and numeric (1, 2) pairs. Coordinates are degrees east and north.
%
% Output arguments
%   epicentres - Canonical coordinates, size (nEvents, 2).
%       Longitudes use [-180, 180); longitude at either pole is zero.
%   names - Normalised catalogue names or encoded coordinates, (nEvents, 1).
%   expectedDates - Catalogue event days, or NaT for unnamed coordinates.
%
% Notes
%   Catalogue and sources: docs/grace-seismic.md. Matching coordinates
%   within 1e-10 degrees resolve to the same catalogue entry as its name.
%
% See also
%   grace2seismic, num2filename
%
% Created by
%   2026/11/12, En-Chi Lee (williameclee@arizona.edu)

function [epicentres, names, expectedDates] = resolveearthquakes(events)
    catalogue = [95.982, 3.295; 142.373, 38.297; -72.733, -35.909];
    catalogueNames = ["Sumatra"; "Tohoku"; "Chile"];
    catalogueDates = [datetime(2004, 12, 26); datetime(2011, 3, 11); datetime(2010, 2, 27)];
    aliases = {{'sumatra', 'sumatra2004', 'sumatraandaman', 'indianocean'}, ...
        {'tohoku', 'tohoku2011', 'japan'}, {'chile', 'chile2010', 'maule'}};
    if isnumeric(events)
        assert(ismatrix(events) && size(events, 2) == 2, ...
            'ULMO:resolveearthquakes:InvalidInput', 'Epicentres must have two columns.');
        events = num2cell(events, 2);
    elseif ischar(events)
        events = {events};
    elseif isstring(events)
        events = cellstr(events);
    end
    assert(iscell(events) && isvector(events) && ~isempty(events), ...
        'ULMO:resolveearthquakes:InvalidInput', 'Supply a nonempty event vector.');
    epicentres = zeros(numel(events), 2);
    names = strings(numel(events), 1);
    expectedDates = NaT(numel(events), 1);
    for iEvent = 1:numel(events)
        event = events{iEvent};
        if ischar(event) || (isstring(event) && isscalar(event) && ~ismissing(event))
            assert(isrow(event) || isstring(event), ...
                'ULMO:resolveearthquakes:InvalidInput', 'Each name must be scalar text.');
            name = lower(replace(string(event), ["ō", "Ō"], "o"));
            name = char(regexprep(name, '[^a-z0-9]', ''));
            index = find(cellfun(@(list) any(strcmp(name, list)), aliases), 1);
            assert(~isempty(index), 'ULMO:resolveearthquakes:UnknownName', ...
                'Unknown earthquake "%s". Use Sumatra, Tohoku, Chile or [longitude, latitude].', event);
            coordinate = catalogue(index, :);
        else
            assert(isnumeric(event) && isreal(event) && isequal(size(event), [1, 2]) && ...
                all(isfinite(event)) && abs(event(2)) <= 90, ...
                'ULMO:resolveearthquakes:InvalidInput', ...
                'Each coordinate entry must be a finite [longitude, latitude] pair with latitude in [-90, 90].');
            coordinate = double(event);
            if coordinate(1) < -180 || coordinate(1) >= 180
                coordinate(1) = mod(coordinate(1) + 180, 360) - 180;
            end
            if abs(coordinate(2)) == 90, coordinate(1) = 0; end
            index = find(all(abs(catalogue - coordinate) <= 1e-10, 2), 1);
        end
        if isempty(index)
            names(iEvent) = string(num2filename(coordinate(1))) + "_" + num2filename(coordinate(2));
        else
            coordinate = catalogue(index, :);
            names(iEvent) = catalogueNames(index);
            expectedDates(iEvent) = catalogueDates(index);
        end
        epicentres(iEvent, :) = coordinate;
    end
end
