# The Box

A geometric indoor wireless multipath prototype: line of sight plus one specular reflection from each surface of a rectangular room.

Read the [educational primer](BACKGROUND.md) for a step-by-step introduction to RSSI, waves, interference, reflection, and the model's equations. Each equation's symbols and operations are explained with worked examples.

The separate [code guide](CODE_GUIDE.md) walks through the MATLAB files, function inputs and outputs, workbook decoding, configuration, plots, and checks.

## Run

From this repository directory in MATLAB:

```matlab
testBoxModel
TheBox_SIMULATION
PlotRecorded_Data
```

The simulation plots the September measurements alongside a coherent seven-ray prediction. The second script plots the densely sampled October measurements. Both show mean linear power across the two orientations, with their spread shaded. Raw orientation measurements remain available in the loader output. Base MATLAB is sufficient; the original RF and Statistics toolbox dependencies have been removed.

## Assumptions to configure

In `TheBox_SIMULATION.m`, set the actual channel center frequency, room and antenna coordinates, and each surface's power reflectance and reflection phase. Defaults of 2.400 GHz and identical 3% power reflectance with a 180-degree phase reversal are illustrative, not measured material properties. The separately recorded 1 m reading provides a provisional power scale and may itself contain multipath.

`loadRecordedData` explicitly assumes hexadecimal samples encode magnitudes of negative dBm readings. Verify the radio documentation; `signed-byte` is available for two's-complement bytes. All downstream values use signed dBm. **`fitLogDistanceModel` now accepts signed dBm and returns a signed reference level.**

This is a narrowband, scalar, single-bounce model. It does not establish whether historical mismatches were caused by channel congestion. See the primer's [calibration discussion](BACKGROUND.md#9-compare-a-geometric-prediction-with-an-empirical-trend) and the code guide's [experimental assumptions](CODE_GUIDE.md#11-experimental-assumptions-and-next-changes) before comparing absolute power or interference nulls.
