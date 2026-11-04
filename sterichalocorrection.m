%% STERICHALOCORRECTION - Monthly global halosteric offsets on native support.
% Syntax
%   h = sterichalocorrection(path, indices, hasDeep, applyCorrection)
%
% reads a native steric MAT file. Positive indices select source months;
% zero indices denote missing months and return NaN. Offsets are in mm.
% The full-column mask requires finite, nonzero values in >=80% of all
% source months. Each month uses identical cosine-latitude weights for the
% full and shallow fields; deep is their difference (zero without deep data).
% With applyCorrection=false, applied offsets are zero but estimated offsets
% are retained when available. This removes a monthly mean, not a fitted trend.
function h = sterichalocorrection(path, indices, hasDeep, applyCorrection)
    m = matfile(path);
    c = load(path, 'lat');
    metadata = whos(m);
    fields = {metadata(arrayfun(@(item) all(item.size > 0), metadata)).name};

    if ~applyCorrection && ~ismember('halostericSl', fields)
        h.total = nan(numel(indices), 1);
        h.total(indices > 0) = 0;
        h.shallow = h.total;
        h.deep = h.total;
        h.estimatedTotal = nan(size(h.total));
        h.estimatedShallow = h.estimatedTotal;
        h.estimatedDeep = h.estimatedTotal;
        h.applied = false; h.supportMask = [];
        h.note = 'No halosteric correction requested.';
        return
    end

    a = double(m.halostericSl);
    if size(a, 1) ~= numel(c.lat), a = permute(a, [2, 1, 3]); end
    mask = sum(isfinite(a), 3) >= .8 * size(a, 3) & sum(a ~= 0 & isfinite(a), 3) >= .8 * size(a, 3);
    w = double(mask) .* cosd(double(c.lat(:)));
    h.total = nan(numel(indices), 1); h.shallow = h.total; h.deep = h.total;

    if hasDeep

        if ismember('shallowHalostericSl', fields)
            b = double(m.shallowHalostericSl);
            if size(b, 1) ~= numel(c.lat), b = permute(b, [2, 1, 3]); end
        else
            d = double(m.deepHalostericSl);
            if size(d, 1) ~= numel(c.lat), d = permute(d, [2, 1, 3]); end
            b = a - d;
        end

    end

    for k = find(indices(:)' > 0)
        v = a(:, :, indices(k)); valid = isfinite(v) & mask;

        if hasDeep
            u = b(:, :, indices(k)); valid = valid & isfinite(u);
        else
            u = v;
        end

        den = sum(w(valid)); assert(den > 0);
        h.total(k) = sum(v(valid) .* w(valid)) / den * 1000;
        h.shallow(k) = sum(u(valid) .* w(valid)) / den * 1000;
        h.deep(k) = h.total(k) - h.shallow(k);
    end

    h.estimatedTotal = h.total;
    h.estimatedShallow = h.shallow;
    h.estimatedDeep = h.deep;

    if ~applyCorrection
        h.total(isfinite(h.total)) = 0; h.shallow(isfinite(h.shallow)) = 0; h.deep(isfinite(h.deep)) = 0;
    end

    h.applied = applyCorrection;
    h.supportMask = mask;
    h.note = 'Native global full-column halosteric support, >=80% finite/nonzero months; identical weights for depth components. SIO retains production-loader total default.';
end
