domain = GeoDomain('oceans');

h = 3;
lon = 0:h:360;
lat = -90:h:90;

[mask, ~, ~] = domainmask(domain, lon, lat, ...
    ForceNew = true, SaveData = false);

assert(isequal(size(mask), [length(lat), length(lon)]));
assert(~any(all(mask, 1)));

lonOrigin = 20;
lon = (0:h:360) + lonOrigin;
lat = -90:h:90;

[mask, ~, ~] = domainmask(domain, lon, lat, ...
    ForceNew = true, SaveData = false);

assert(isequal(size(mask), [length(lat), length(lon)]));
assert(~any(all(mask, 1)));
