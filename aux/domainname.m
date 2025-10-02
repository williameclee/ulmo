%% DOMAINNAME
% Returns the full name of a domain given its function name.
%
% Syntax
%   name = domainname(domain)
%   name = domainname(domain, fmt)
%
% Input arguments
%   domain - The domain name
%       The domain name should be the name of a function that returns the
%       domain vertices.
%   fmt - The format of the domain name
%       - 'abbrevation' or 'abbr' - Abbrevation
%       - 'short' - Short name
%       - 'long' - Long name
%       The default value is 'short'.
%
% Output arguments
%   name - The full name of the domain
%
% Example
%   >>  domainname('namerica')
%   'N America'
%   >>  domainname('namerica', 'long')
%   'North America'
%
% Last modified by
%   2025/10/01, williameclee@arizona.edu (@williameclee)

function domainName = domainname(domain, format)

    arguments (Input)
        domain {mustBeA(domain, {'char', 'string', 'cell'})}
        format {mustBeA(format, {'char', 'string'}), mustBeMember(format, {'abbrevation', 'abbr', 'short', 'long'})} = 'short'
    end

    arguments (Output)
        domainName (1, :) char
    end

    if iscell(domain)
        domain = domain{1};

        if ~ischar(domain) && ~isstring(domain)
            error ('ULMO:domainname:invalidDomain', ...
            'If domain is a cell array, the first element should be a string.');
        end

    end

    domain = lower(domain);

    if strcmp(format, 'abbr')
        format = 'abbrevation';
    end

    switch domain
        case 'namerica'

            switch format
                case 'abbrevation'
                    domainName = 'NA';
                case 'short'
                    domainName = 'N America';
                case 'long'
                    domainName = 'North America';
            end

        case 'samerica'

            switch format
                case 'abbrevation'
                    domainName = 'SA';
                case 'short'
                    domainName = 'S America';
                case 'long'
                    domainName = 'South America';
            end

        case 'continents'
            domainName = 'All continents';

        case 'oceans'

            switch format
                case 'abbrevation'
                    domainName = 'Global';
                case 'short'
                    domainName = 'Global mean';
                case 'long'
                    domainName = 'Global mean';
            end

        case 'atlantic'

            switch format
                case 'abbrevation'
                    domainName = 'Atl.';
                case 'short'
                    domainName = 'Atlantic';
                case 'long'
                    domainName = 'Atlantic Ocean';
            end

        case 'natlantic'

            switch format
                case 'abbrevation'
                    domainName = 'N. Atl.';
                case 'short'
                    domainName = 'N Atlantic';
                case 'long'
                    domainName = 'North Atlantic Ocean';
            end

        case 'satlantic'

            switch format
                case 'abbrevation'
                    domainName = 'S. Atl.';
                case 'short'
                    domainName = 'S Atlantic';
                case 'long'
                    domainName = 'South Atlantic Ocean';
            end

        case 'pacific'

            switch format
                case 'abbrevation'
                    domainName = 'Pac.';
                case 'short'
                    domainName = 'Pacific';
                case 'long'
                    domainName = 'Pacific Ocean';
            end

        case 'npacific'

            switch format
                case 'abbrevation'
                    domainName = 'N. Pac.';
                case 'short'
                    domainName = 'N Pacific';
                case 'long'
                    domainName = 'North Pacific Ocean';
            end

        case 'spacific'

            switch format
                case 'abbrevation'
                    domainName = 'S. Pac.';
                case 'short'
                    domainName = 'S Pacific';
                case 'long'
                    domainName = 'South Pacific Ocean';
            end

        case 'indian'

            switch format
                case 'short'
                    domainName = 'Indian';
                case 'long'
                    domainName = 'Indian Ocean';
                otherwise
                    domainName = 'Indian';
            end

        otherwise
            domainName = capitalise(domain);

    end

end
