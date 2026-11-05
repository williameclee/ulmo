%% STERIC2GMSL - Read or cache one native global-mean steric time series.
%
% Syntax
%   [gmsl, dates] = STERIC2GMSL(product, Type="halosteric")
%   [gmsl, dates] = STERIC2GMSL(product, Type="thermosteric", Depth="shallow")
%   [gmsl, dates] = STERIC2GMSL(__, SaveData=false, ForceNew=true)
%
% Input arguments
%   product - Prepared product identifier accepted by StericProduct.
%       The native MAT file must already exist at STERICFILEPATH(product).
%   Type (name-value, required) - Selected component.
%       total, thermosteric, halosteric, non-thermosteric (total minus
%       thermosteric), or non-halosteric (total minus halosteric).
%       No correction is applied. There is no default selection.
%   Depth (name-value) - full, shallow, or deep, as supplied by the product.
%       Unsupported component/depth combinations raise an error.
%       Default depth: full.
%   Unit (name-value) - m or mm, for the returned series only.
%       Default unit: m.
%   TimeFormat (name-value) - datetime or datenum, for returned dates only.
%       Default format: datetime.
%   SaveData (name-value) - Whether to append all calculated means to the
%       native MAT file. Existing complete caches can still be read if false.
%       Default flag: true.
%   ForceNew (name-value) - Whether to recalculate all applicable means.
%       With SaveData=false, no file is changed.
%       Default flag: false.
%   BeQuiet (name-value) - Whether to suppress messages.
%       0 and 0.5 print calculation/save messages; 1 suppresses all messages.
%       Default level: 0.5.
%
% Input file content
%   lon - Regular longitude coordinate vector spanning one complete cycle.
%       A duplicated cyclic endpoint is permitted.
%       Size: vector with nLon elements.
%       Unit: degrees east.
%   lat - Strictly monotone latitude coordinate vector in [-90, 90].
%       Spacing may be uniform or nonuniform.
%       Size: vector with nLat elements.
%       Unit: degrees north.
%   dates - Native datetime or datenum vector, strictly increasing.
%       Size: vector with N elements.
%       Unit: days for datenum values.
%   stericSl, thermostericSl, halostericSl - Available native component grids.
%       All grids must share an orientation. Square spatial grids are
%       interpreted as (lat, lon, time), as in the native aggregate format.
%       Missing or empty components are unavailable, not zero observations.
%       Size: (nLat, nLon, N) or (nLon, nLat, N).
%       Unit: m.
%   shallow<Component>Sl, deep<Component>Sl - Supplied depth-layer grids.
%       Component is Steric, Thermosteric, or Halosteric. Each split requires
%       full, shallow and deep grids. Total and halosteric splits must agree
%       when both components are present; missing layers are errors.
%       Size: same as the corresponding full grid.
%       Unit: m.
%
% Output arguments
%   gmsl - Selected component global mean on common observed support.
%       Size: (N, 1).
%       Unit: Unit.
%   dates - Native timestamps in the selected TimeFormat.
%       Size: (N, 1).
%       Unit: days for datenum output.
%
% Output file content
%   <component>Gmsl - Individual cached mean vectors.
%       Components are steric, thermosteric, halosteric, nonThermosteric and
%       nonHalosteric. Supported shallow/deep names prefix and capitalise
%       the component, e.g. shallowHalostericGmsl. Original variables remain.
%       Size: (N, 1).
%       Unit: m.
%   gmslCreationDate - Datetime when newly calculated means were saved.
%       Size: (1, 1).
%       Reusing a cache or setting SaveData=false leaves this value unchanged.
%
% Notes
%   All applicable means are calculated together, then only one is returned.
%   The fixed support reference is full halosteric, otherwise full total,
%   otherwise full thermosteric. Cells need >=80% finite and >=80% finite
%   nonzero reference samples. Each date uses common finite support across
%   all supplied components/layers. Uniform latitude grids use cosine-latitude
%   weights. Nonuniform grids use spherical band areas between arithmetic
%   latitude midpoints: abs(diff(sind(boundaries))). Outer boundaries extend
%   half the adjacent spacing, clipped to [-90, 90]. Spacing differences
%   within 3e-5 degrees are treated as coordinate rounding. A duplicated
%   longitude endpoint is excluded. No supported area gives NaN. Valid zero
%   layers remain in the denominator. Residuals are differences on this
%   common support; total need not equal thermosteric plus halosteric.
%   Supplied layers are averaged directly; full means are their sum and
%   checked against the full grid with tolerance 32*eps(source precision)
%   times max(1, the weighted mean absolute full/shallow/deep values).
%   Source precision is single if any supplied grid is single, else double.
%   Complete caches bypass grid reads. Source/convention changes require
%   explicit ForceNew=true and manual rebuilding of affected downstream
%   caches. This function does not interpolate, remove drift or rebuild grids.
%
% See also
%   steric2lonlatt, stericfilepath, stericfieldschema
%
% Created by
%   2026/11/04, En-Chi Lee (williameclee@gmail.com)
%
% Last modified
%   2026/11/05, OpenAI Codex

function [gmsl, dates] = steric2gmsl(product, options)
    arguments (Input)
        product (1, 1) StericProduct
        options.Type (1, 1) string {mustBeMember(options.Type, ...
            ["total", "thermosteric", "halosteric", "non-thermosteric", "non-halosteric"])}
        options.Depth (1, 1) string {mustBeMember(options.Depth, ["full", "shallow", "deep"])} = "full"
        options.Unit {mustBeTextScalar, mustBeMember(options.Unit, ["m", "mm"])} = "m"
        options.TimeFormat {mustBeTextScalar, mustBeMember(options.TimeFormat, ["datetime", "datenum"])} = "datetime"
        options.SaveData (1, 1) logical = true
        options.ForceNew (1, 1) logical = false
        options.BeQuiet (1, 1) double {mustBeMember(options.BeQuiet, [0, 0.5, 1])} = 0.5
    end

    if ~isfield(options, 'Type')
        error('ULMO:steric2gmsl:TypeRequired', 'Specify Type explicitly.');
    end

    [means, dates, schema] = stericgmsldata(product, Type=options.Type, ...
        Depth=options.Depth, ForceNew=options.ForceNew, ...
        SaveData=options.SaveData, BeQuiet=options.BeQuiet);
    iType = schema.types == options.Type;
    iDepth = schema.depths == options.Depth;
    cacheNames = schema.meanNames;
    gmsl = double(means.(cacheNames(iType, iDepth)));
    if strcmp(options.Unit, 'mm')
        gmsl = 1e3 * gmsl;
    end
    if strcmp(options.TimeFormat, 'datenum')
        dates = datenum(dates); %#ok<DATNM>
    end
end
