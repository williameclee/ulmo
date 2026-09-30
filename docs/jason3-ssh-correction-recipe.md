# Jason-3 correction recipe: source contract and DT2024 limits

Step 1 of [the implementation plan](jason3-ssh-correction-plan.md), reviewed
against primary sources on 2026-09-30. This is a scientific specification;
no correction is enabled by this change.

**Decision:** the downloaded JPL file supports a temporal Jason-3 range
adjustment. It is insufficient by itself to reproduce the ESA global correction
or a calibrated correction for CMEMS all-satellite DT2024. Define the source
processing below, but keep DT2024 integration blocked on the spatial reduction
and mission-transfer evidence listed at the end. Step 1's scientific checkpoint
is therefore still open; a documented unknown is not a validated default.

## Product and source identity

The target is monthly CMEMS
`cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1M-m`, variable `sla`.
The local subset covers January 1993 through April 2025. This is different from
the two-satellite, 0.25-degree C3S product used in the ESA example [2]. Neither
the filename `DT2024.mat` nor the word DT2024 identifies a mission-transfer
algorithm.

The input already downloaded is
`$IFILES/SSH/JASON_3_PD_CORRECTION_20230925.txt`, Brown, Willis and Fournier,
DOI [10.5067/J3L2G-PDCOR][1], version F. The text header is dated 2023-09-25;
the catalog lists publication on 2023-12-12. Use the file's actual timestamps,
not the catalog's open-ended coverage, to determine support. The supplementary
adjustment was not included in Jason-3 GDR-F [1].

Ignore `HDR` header lines and blank lines; require four fields per data row:

| Column | Meaning | Interpretation |
| --- | --- | --- |
| 1 | Jason-3 cycle | Integer identifier, not an elapsed-time axis |
| 2 | Pass | Integer identifier within the cycle |
| 3 | Pass midpoint | `yyyyMMdd'T'HHmmss`; preserve supplied mission times without local-time conversion |
| 4 | Range adjustment | Signed centimeters |

Read-only inspection of the downloaded file gives:

| Property | Observed value |
| --- | --- |
| Data rows | 68,813 |
| First midpoint | 2016-02-12 01:39:15, cycle 0/pass 117 |
| Last midpoint | 2023-08-18 02:55:53, cycle 348/pass 254 |
| Duplicate or non-increasing timestamps | None |
| Value range | -0.044572 to +0.190390 cm |
| Largest gap | 17.915011574 days, 2022-04-07 13:08:35 to 2022-04-25 11:06:12 |

That gap spans the move to the interleaved orbit, not a long series of ordinary
missing passes: DUACS identifies cycle 226 as the final nominal-orbit cycle and
cycle 300 as the first interleaved cycle [5]. The presence of correction records
after April 2022 does not mean Jason-3 remains the reference mission.

## Sign and temporal processing

For a file value `p` in centimeters, return the signed range correction
`c = 0.01*p` in meters. Apply it as `sshCorrected = sshOriginal - c`.
For example, `p = +0.1 cm` lowers SSH by 1 mm, and `p = -0.04 cm` raises
SSH by 0.4 mm. Do not take an absolute value or demean the source silently.
This sign is specified by JPL [1].

The following is a **ULMO temporal approximation**, not ESA's spatially
weighted global series. It can be implemented independently in the helper,
with its limited meaning documented:

1. Validate finite values, integer cycle/pass identifiers and strictly
   increasing unique timestamps before any requested-time filtering. Reject
   malformed rows and conflicting timestamps rather than averaging them away.
2. Within one orbit phase, linearly interpolate the pass-midpoint values in
   elapsed time. JPL permits temporal linear interpolation of missing passes
   because the adjustment changes slowly [1]. Retain bracketing records outside
   the requested interval. Do not interpolate across the nominal/interleaved
   orbit gap for an automatically supported monthly result.
3. Form full calendar-month means of this piecewise-linear series in meters.
   Split intervals at month boundaries and integrate trapezoids exactly:
   `cMonth = sum((cLeft+cRight)/2 .* elapsedSeconds) / monthSeconds`.
   Include all knots inside the month, plus interpolated boundary values.
   This weights elapsed time, not the number of available passes. Leap-year
   February has its actual duration. A boundary needs an exact sample or a
   bracketing pair; do not extrapolate to manufacture support.
4. Label each monthly mean at the same epoch as its native DT2024 monthly SSH
   value. For a standalone series use the true calendar midpoint (halfway
   between consecutive month starts), subject to confirming the loader's
   convention in step 2. When called by the loader, use its actual dates.
5. Only then apply the requested `Interpolation` method to monthly values for
   other output epochs. Preserve the loader's supported methods, date formats
   and meter/millimeter units. The option changes this final resampling, not
   the pass-to-month integration. `QueryDates` must match the output epochs
   exactly; it must not bridge an unsupported month, even with nearest or spline.

Full calendar support alone permits March 2016--March 2022 and May
2022--July 2023 for this temporal approximation. February 2016 and August 2023
are partial source months; April 2022 crosses the orbit gap. Reject them, and
dates outside support, with a clear error identifying the unavailable month.
These intervals are **not a declaration that CMEMS can be corrected throughout
them**. A request containing an unsupported epoch fails as a whole; no silent
truncation, zero fill, endpoint hold, or partial-month renormalization.

The default all-time helper request returns the supported monthly records,
retaining the April 2022 gap; interpolation must treat the two supported
segments separately. A shortened request must give the same monthly values as
the corresponding subset of a full request.

For monthly SSH, point evaluation at the middle of the month, a pass-count
average, and nearest-midpoint downsampling in `interptemporal` are not
substitutes for the integral above. The calendar averaging specification is
a proposed ULMO convention; it is not claimed to reproduce a published CMEMS
correction series.

## What is required for a global L4 adjustment

ESA's SLBC algorithm uses C3S daily grids and computes the Jason-3 adjustment
by interpolating onto altimetry tracks, aggregating into 3-degree longitude by
1-degree latitude cells, and computing a global mean. It then subtracts that
same global value everywhere [2, section 2.1.3, printed page 18]. AVISO describes
area/ocean-coverage weighting and approximately 10-day aggregation for its
global-mean processing [3]. This supports a spatially uniform approximation;
it does not establish a simple average of pass values as equivalent.

The four-column file has no along-track locations, valid-ocean mask, spatial
weights or CMEMS mission-combination weights. Therefore:

- For the temporal approximation above, no additional download is needed.
- For the ESA-style global approximation, obtain a separately documented global
  Jason-3 correction series, or obtain the track observations, quality selection,
  ocean mask and weighting procedure needed to derive it. A corrected GMSL
  series alone is insufficient: differences may include GIA, sea-state-bias
  and other processing changes [3]. No suitable correction-only global series
  has been identified or downloaded in this step.
- Even a global series remains a uniform approximation for CMEMS all-satellite
  L4. It cannot recover the spatially varying effect of rerunning cross-calibration
  and mapping with corrected Jason-3 observations. Results outside Jason-3's
  approximately 66-degree latitude coverage need the same explicit qualification.

Leclercq et al. apply the radiometer adjustment to DT2024 GMSL and C3S Pacific
grids [6, Methods/Data]. That is evidence of scientific use, but those methods
do not supply the missing all-satellite CMEMS transfer recipe.

## Mission transitions and coverage: no automatic transfer yet

Keep two quantities distinct: `cJ3(t)`, the Jason-3 instrumental range adjustment,
and `cDT(t)`, the scalar adjustment adopted for a merged CMEMS record. The source
file defines the former. It does not uniquely define the latter.

| Boundary | Evidence | Required implementation decision |
| --- | --- | --- |
| Before Jason-3 | File begins in February 2016 | Do not infer a CMEMS start date or pre-2016 zero extension from file coverage |
| Entry into the reference record | Intermission biases are estimated during tandem overlap [3] | Establish the target product's overlap, datum and effective transfer, rather than subtracting the first sample to force a zero |
| Sentinel-6 handover | DUACS's 2022-04-06 reference switch is explicitly NRT [4] | Do not hard-code that date for reprocessed DT2024 |
| Interleaved Jason-3 | DUACS reintroduced Jason-3 after its orbit move [5] | Do not assume its contribution to all-satellite L4 vanishes at handover |
| End of text file | Last midpoint is in August 2023 | Do not extend the drift or hold the last value through April 2025 |

Why the offset matters: if the new mission is tied to the old mission through
an overlap mean, changing the old mission by `-cJ3(t)` changes the transferred
datum by an overlap-weighted mean of that correction. Under that simplified
model a successor may inherit a constant offset even after Jason-3 ceases to
be the reference. This is an inference from tandem calibration, not a measured
CMEMS offset. It neither justifies dropping the adjustment to zero nor proves
that holding the value on a switch date is correct. Regional cross-calibration
and continued use of interleaved Jason-3 add further differences.

Consequently no automatic `cJ3 -> cDT` conversion is approved by this recipe.
Neither the raw, initial-value-anchored, nor April-2022-capped September audit
sensitivity scenario is a validated production recipe. Until the transfer is
specified, a future DT2024 correction request must report an unsupported recipe,
rather than return plausible but unverified corrected SSH. Ordinary uncorrected
calls remain available.

## Review and next-step acceptance

The independent helper work may proceed with source parsing, sign/units,
calendar integration and explicit coverage errors. Synthetic tests should cover
both signs; nonconstant sub-month knots; leap February; month-boundary bracketing;
request-window invariance; partial months; the orbit gap; and each final
interpolation method without crossing that gap. Do not require the user's archive
in unit tests.

Before the loader integration checkpoint can pass, record:

1. The chosen global-reduction input and algorithm, or an explicit scientific
   decision to adopt the temporal approximation instead, including its limits.
2. CMEMS DT2024-specific entry/exit dates, overlap weighting and datum transfer,
   including the treatment of interleaved Jason-3 and the supported final date.
3. A reference calculation across both transitions and a measured trend effect;
   do not tune that result to reconcile MEaSUREs and DT2024.

This dependency is the remaining part of step 1. It is not solved by implementing
an interpolation function. No provider messages were sent and no new data
dependencies were downloaded apart from the public algorithm document.

Retain the agreed storage rules: no helper disk cache, no cache versioning,
uncorrected native/interpolated MAT files, and subtraction once in memory only.
Do not add a spatially shared systematic correction uncertainty as independent
cell noise in `sshSigma`.

## Sources

[1]: https://podaac.jpl.nasa.gov/dataset/JASON_3_PD_CORRECTION
[2]: https://climate.esa.int/media/documents/2024-10-31_SLBC_CCI-DT-041-MAG_ATBD_D2-4_V2_0_signed.pdf
[3]: https://www.aviso.altimetry.fr/en/data/products/ocean-indicators-products/mean-sea-level/processing-and-corrections.html
[4]: https://duacs.cls.fr/duacs-system-description/operational-news/duacs-nrt-system-impacting-version-changes/apr-2022-duacs-19-1-0/
[5]: https://duacs.cls.fr/duacs-system-description/operational-news/duacs-nrt-system-impacting-version-changes/may-2022-duacs-19-2-0/
[6]: https://www.nature.com/articles/s43247-025-03149-5

- [1 — JPL dataset catalog and downloaded file header][1].
- [2 — ESA SLBC ATBD, issue 2.0, internal date 2025-03-31][2].
  The filename contains 2024-10-31; the document's internal issue/date identifies
  the version inspected here. Sections 2.1.2--2.1.4 were read in full.
- [3 — AVISO processing, global averaging and intermission biases][3].
- [4 — DUACS NRT Sentinel-6 introduction][4].
- [5 — DUACS Jason-3 orbit change and reintroduction][5].
- [6 — Leclercq et al. (2026), Methods/Data][6].
