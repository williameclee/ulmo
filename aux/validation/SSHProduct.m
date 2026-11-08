%% SSHPRODUCT - Identify prepared sea surface height products for ssh2lonlatt.
% NASASSH selects the native weekly maps. NASASSHMonthly selects calendar-month means.
%
% See also
%   SSH2LONLATT, PROCESSSSHDATANASA, PROCESSSSHMONTHLYNASA
%
% Created by
%   2025/10/01, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/08, En-Chi Lee (williameclee@gmail.com)

classdef SSHProduct < uint8

    enumeration
        MEaSUREs (0)
        ORAS5 (1)
        DT2024 (2)
        CMEMS (3)
        NASASSH (4)
        NASASSHMonthly (5)
    end

end
