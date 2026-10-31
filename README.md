# The Box

**A MATLAB exploration of indoor wireless multipath using room geometry.**

The Box approximates a room as a rectangular enclosure and predicts received signal strength from a direct line-of-sight path plus one reflection from each wall, the floor, and the ceiling. It compares that prediction with measurements collected using Digi XBee Series 2 (S2) 2.4 GHz Zigbee radios.

The project asks whether a simple geometric model can explain how received power changes as a receiver moves through a room. It is an experimental model with explicit assumptions, rather than a calibrated indoor coverage predictor.

## Getting started

Open this repository as the current folder in MATLAB and run:

```matlab
testBoxModel        % Check signal math, reflection geometry, and data loading
TheBox_SIMULATION  % Compare the model with the September measurements
PlotRecorded_Data  % Inspect the densely sampled October measurements
```

The project uses base MATLAB, including `readcell` for the bundled `.xlsx` files. No RF or Statistics toolbox is required. The current implementation has been exercised in MATLAB R2026b; compatibility with older releases has not been established.

`TheBox_SIMULATION` clears workspace variables and closes existing figures before running. It locates the measurement workbook relative to the script; keep the helper functions on MATLAB's path.

## What the simulation shows

Running `TheBox_SIMULATION` creates three figures:

1. **Measured versus predicted RSSI.** A dense model curve is compared with measured mean power at each recorded distance. Shading spans the two antenna-orientation power means.
2. **Absolute prediction error.** One bar per measured distance shows the magnitude of the difference between predicted and measured levels. Predictions are evaluated at the exact measurement locations.
3. **Room and propagation paths.** A 3D view shows the antennas, direct path, and six single-bounce paths at one representative receiver position.

RSSI plots display positive magnitudes labeled **RSSI (-dBm)**. A plotted value of 43 represents -43 dBm, so **lower values mean stronger reception**. Internal calculations use signed dBm. Errors are expressed in dB.

The console prints root-mean-square error (RMSE), mean absolute error (MAE), and signed bias. The workspace retains predictions and errors in `measuredRayDbm`, `measuredMeanDbm`, `errorDb`, `absoluteErrorDb`, and `errorSummary`. A positive bias means the model predicts a stronger signal than measured.

The main comparison does not fit a trend or optimize parameters against the distance sweep. `fitLogDistanceModel.m` remains available as an optional analysis helper, but the simulation does not call it.

## How the model works

For each receiver position, the model:

1. Mirrors the transmitter across each surface to locate a specular reflection point using the image-source method.
2. Calculates the direct and six reflected path lengths.
3. Assigns each path an amplitude based on distance spreading and a complex reflection coefficient.
4. Calculates propagation phase from path length and carrier wavelength.
5. Sums complex amplitudes, then converts the squared magnitude into received power.

The essential relationships are:

$$
h=\sum_{p=0}^{6}\Gamma_p\left(\frac{d_0}{L_p}\right)^{n/2}
 e^{-j2\pi(L_p-L_0)/\lambda},
\qquad
P_{\mathrm{rx,dBm}}=P_{\mathrm{ref,dBm}}+20\log_{10}|h|.
$$

Here $L_p$ is a path length, $L_0$ is the direct-path length, $d_0$ is the reference distance, $\lambda$ is wavelength, and $\Gamma_p$ is a complex field reflection coefficient. The direct path uses $\Gamma_0=1$. The default power-distance exponent is $n=2$, corresponding to free-space spreading along each ray. Subtracting $L_0$ selects a common phase reference without changing predicted power.

See the [educational primer](BACKGROUND.md) for definitions, derivations, and worked examples, or the [code guide](CODE_GUIDE.md) for the MATLAB implementation.

## Measurements and orientation averaging

The bundled experiments contain five packet readings in each of two antenna orientations at every position. Each recorded reading followed a new packet arrival; the repeated samples were not repeated queries of one unchanged packet measurement.

| Workbook | Positions | Distance range | Original unit | Use |
|---|---:|---|---|---|
| `dataCollection_9_10_22.xlsx` | 25 | 0.9144–8.2296 m | Feet | Model comparison with stored room geometry |
| `dataCollection_10_9_22.xlsx` | 37 | 3.9624–4.8768 m | Inches | Dense measurement visualization |

Both workbooks also include a separate 1 m reference row. The loader excludes this row from the spatial sweep.

For XBee S2 `DB` values, hexadecimal bytes encode the magnitude of negative dBm: `21` means -33 dBm and `2B` means -43 dBm. The loader's `negative-magnitude` option implements this conversion. See [Digi's XBee/XBee-PRO ZB user guide](https://www.digi.com/resources/documentation/digidocs/pdfs/90000976.pdf) and the [S2 measurement notes](CODE_GUIDE.md#s2-packet-measurements-and-the-collection-procedure).

The plots average the ten readings in **linear power**, then convert back to dBm. The shaded range spans the independently power-averaged orientations. It shows observed orientation spread, **not a confidence interval**. Averaging two orientations reduces dependence on a single sampled orientation but does not remove antenna directionality or polarization effects. Raw samples and separate-orientation statistics remain available from `loadRecordedData`.

`PlotRecorded_Data` uses this same display convention for October. It does not assume that October used the September geometry or generate an October prediction.

## Configure an experiment

Edit the configuration near the beginning of [`TheBox_SIMULATION.m`](TheBox_SIMULATION.m):

| Setting | Current value | Meaning |
|---|---|---|
| `room_m` | `[12.7127, 14.0462, 3]` | Room dimensions along x, y, z, in meters |
| `txPos_m` | `[2.2479, 7.8486, 1]` | Transmitter coordinates, in meters |
| `rxHeight_m` | `1` | Receiver height, in meters |
| `frequency_Hz` | `2.4e9` | Assumed carrier center frequency |
| `pathExponent` | `2` | Per-ray power-distance exponent |
| `powerReflectance` | Six entries of `0.03` | Illustrative reflected power fractions |
| `reflectionPhase_rad` | Six entries of `pi` | Illustrative reflection phase shifts |

Surface settings are ordered **x=0, x=room length, y=0, y=room width, floor, ceiling**. Power reflectance is converted to field magnitude using its square root; reflection phase is configured separately. The defaults are not measured material properties or angle-dependent Fresnel calculations.

The receiver follows positive x from the transmitter, with the same y-coordinate and configured height. With the default equal antenna heights, the x-offset equals Tx–Rx separation. Changing the route or height also requires updating distance-axis and measurement comparisons appropriately. Antennas must be distinct and strictly inside the room.

The separate September 1 m measurement supplies a provisional reference of -33 dBm. It may contain multipath itself, so it is not an independently isolated direct-path calibration. Actual channel frequency, heights, surface coefficients, and reference scaling should be established before interpreting exact null locations or absolute accuracy.

## Current comparison and limitations

With the bundled September data and current default settings, comparison against mean linear power over both orientations gives:

| Metric | Value |
|---|---:|
| RMSE | 3.66 dB |
| MAE | 2.93 dB |
| Bias, predicted minus measured | -0.15 dB |

These describe this dataset and configuration. They do not establish accuracy in other rooms, and agreement in average level does not establish accurate prediction of spatial fading. The reference scale comes from the separate recorded 1 m reading; surface coefficients remain illustrative.

The model is narrowband, scalar, and limited to single-bounce reflections. It omits antenna radiation patterns, explicit polarization, multiple bounces, furniture, people, diffraction, diffuse scattering, wall transmission, receiver noise, and bandwidth averaging. The historical results do not determine whether channel congestion caused the remaining mismatch.

The tests check mathematical behavior and workbook regression values, rather than establishing experimental predictive accuracy. For quantitative validation, confirm the experimental inputs and evaluate any estimated parameters on held-out measurements.

## Repository guide

| File | Purpose |
|---|---|
| [`BACKGROUND.md`](BACKGROUND.md) | Detailed educational primer on RSSI, interference, reflection, and geometry |
| [`CODE_GUIDE.md`](CODE_GUIDE.md) | Function contracts, data flow, configuration, plots, and measurement notes |
| [`TheBox_SIMULATION.m`](TheBox_SIMULATION.m) | Main model and September comparison |
| [`boxReflectionPaths.m`](boxReflectionPaths.m) | Image-source reflection geometry |
| [`combinePaths.m`](combinePaths.m) | Coherent amplitude sum and power conversion |
| [`loadRecordedData.m`](loadRecordedData.m) | Historical workbook decoding and averaging |
| [`PlotRecorded_Data.m`](PlotRecorded_Data.m) | October measurement plot |
| [`plotBoxModel.m`](plotBoxModel.m), [`plotcube.m`](plotcube.m) | Room and antenna visualization |
| [`fitLogDistanceModel.m`](fitLogDistanceModel.m) | Optional signed-dBm log-distance fit |
| [`testBoxModel.m`](testBoxModel.m) | Numerical invariants and workbook checks |
| [`References/`](References/) | Background reading retained with the project |
