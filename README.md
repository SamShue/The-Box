# The Box

A geometric indoor wireless multipath prototype: line of sight plus one specular reflection from each surface of a rectangular room.

Read [the background and mathematical model](BACKGROUND.md) for RSSI, signal superposition, reflection amplitude and phase, image-source geometry, and experimental limitations.

## Run

From this repository directory in MATLAB:

```matlab
testBoxModel
TheBox_SIMULATION
PlotRecorded_Data
```

The simulation plots the September measurements alongside a coherent seven-ray prediction and an independent fitted trend. The second script plots the densely sampled October measurements. Both preserve the two measurement rotations. Base MATLAB is sufficient; the original RF and Statistics toolbox dependencies have been removed.

## Assumptions to configure

In `TheBox_SIMULATION.m`, set the actual channel center frequency, room and antenna coordinates, and each surface's power reflectance and reflection phase. Defaults of 2.400 GHz and identical 3% power reflectance with a 180-degree phase reversal are illustrative, not measured material properties. The separately recorded 1 m reading provides a provisional power scale and may itself contain multipath.

`loadRecordedData` explicitly assumes hexadecimal samples encode magnitudes of negative dBm readings. Verify the radio documentation; `signed-byte` is available for two's-complement bytes. All downstream values use signed dBm. **`fitLogDistanceModel` now accepts signed dBm and returns a signed reference level.**

This is a narrowband, scalar, single-bounce model. It does not establish whether historical mismatches were caused by channel congestion. See [calibration and limitations](BACKGROUND.md#historical-experiments-and-calibration) before comparing absolute power or interference nulls.
