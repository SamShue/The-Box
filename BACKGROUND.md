# Background: RSSI, multipath, and The Box

## 1. What RSSI measures

**Received Signal Strength Indicator (RSSI)** is a receiver's estimate of received radio signal strength. Depending on the device, it may be a relative index, an encoded register value, or a calibrated power estimate in dBm. Its measurement bandwidth, averaging interval, and treatment of noise and other transmitters depend on the radio. RSSI alone does not identify individual propagation paths or distinguish every cause of a weak link. See the [IEEE RSSI overview](https://technav.ieee.org/topic/received-signal-strength-indicator/).

For a power measurement in dBm:

$$
P_{\mathrm{dBm}}=10\log_{10}\left(\frac{P}{1\,\mathrm{mW}}\right),
\qquad
P_{\mathrm{mW}}=10^{P_{\mathrm{dBm}}/10}.
$$

Thus -40 dBm is stronger than -50 dBm: it represents ten times the power. A difference in dB is a ratio; a value in dBm is an absolute power relative to 1 mW.

The historical spreadsheets contain hexadecimal samples such as `21` and `2B`. Hexadecimal is a number representation, not a power unit: these decode to 33 and 43. The original code treated them as positive magnitudes of negative dBm readings. The updated loader preserves that **unverified assumption** explicitly, converting them to -33 and -43 dBm. It also supports two's-complement bytes through `signed-byte`. Neither interpretation establishes the actual radio calibration; check the device documentation before interpreting predictions as absolute power.

A mean in dB and a mean in linear power answer different questions:

$$
\overline{P}_{\mathrm{dBm,log}}=\frac{1}{N}\sum_i P_{i,\mathrm{dBm}},
\qquad
\overline{P}_{\mathrm{dBm,power}}=10\log_{10}\left(\frac{1}{N}\sum_i10^{P_{i,\mathrm{dBm}}/10}\right).
$$

The loader returns both, and retains the raw samples and separate rotation means. The plots use mean logarithmic readings and label them accordingly.

## 2. How wireless signals combine

A transmitter radiates an electromagnetic wave. At a receiver, its direct arrival and reflected copies superpose as **fields or voltages**, before their combined power is measured. Adding their dBm readings is not a valid substitute.

For a sinusoid at frequency $f$, write one received component as:

$$
v_p(t)=\operatorname{Re}\{a_p e^{j2\pi ft}\},
\qquad a_p=|a_p|e^{j\phi_p}.
$$

The complex number $a_p$ is a phasor: its magnitude describes amplitude and its angle describes phase. With compatible voltage normalization and impedance, the total power is proportional to:

$$
P\propto\left|\sum_p a_p\right|^2.
$$

For two paths, this gives:

$$
P=P_0+P_1+2\sqrt{P_0P_1}\cos(\phi_1-\phi_0).
$$

Equal-amplitude arrivals in phase double the amplitude and quadruple the power, a gain of about 6.02 dB over either arrival alone. Equal-amplitude arrivals 180 degrees apart cancel in the ideal model. Unequal amplitudes leave a residual. With multiple reflections, reflected paths also interfere with one another.

For a path of length $L$, propagation introduces delay and carrier phase:

$$
\tau=\frac{L}{c},\qquad \lambda=\frac{c}{f},\qquad
\phi_{\mathrm{prop}}=-\frac{2\pi L}{\lambda},
$$

using the $e^{j2\pi ft}$ time convention. Only relative phase affects total power. Multiplying every component by the same unit-magnitude phasor changes no predicted power. `combinePaths.m` therefore subtracts the first path length in the phase exponent; its returned channel amplitude is relative to that path's propagation phase.

At 2.400 GHz, wavelength is approximately 12.49 cm. A 1 cm change in relative path length changes relative phase by approximately 29 degrees. This explains why small movements or geometry errors can move interference peaks and nulls.

These equations describe coherent copies of the same narrowband signal. Independent transmissions generally do not maintain a fixed phase relationship, so a busy channel cannot be modeled by adding another permanently coherent ray. For modulated signals with appreciable delay spread across their bandwidth, use a delayed-signal model:

$$
y(t)=\sum_p g_p\,x(t-\tau_p)+\eta(t).
$$

Here $x(t)$ is the transmitted complex envelope, $g_p$ includes carrier phase and path gain, and $\eta(t)$ represents noise and unmodeled interference. The Box currently uses a single-frequency approximation rather than this full waveform model. See [MathWorks' multipath channel formulation](https://www.mathworks.com/help/comm/ref/comm.raytracingchannel-system-object.html).

## 3. What a reflection changes

A smooth planar surface produces a specular reflection: the incident and reflected rays lie in the plane containing the surface normal, with equal angles to that normal. A reflected signal also changes amplitude, phase, and potentially polarization.

In a scalar model, represent one reflection by a complex **field reflection coefficient**:

$$
\Gamma=|\Gamma|e^{j\theta_\Gamma},\qquad
a_{\mathrm{reflected}}=\Gamma\,a_{\mathrm{incident}}.
$$

The magnitude scales field amplitude; the angle adds phase. For a wave reflected back into the same medium, the corresponding power reflectance is $R=|\Gamma|^2$. Consequently, 3% power reflectance means a field coefficient magnitude of $\sqrt{0.03}\approx0.173$, not 0.03. A negative real field coefficient represents a phase reversal of 180 degrees.

The actual coefficient depends on incidence angle, frequency, material electromagnetic properties, and polarization. Fresnel equations describe the two polarization components at an ideal planar interface; a full polarized model must track their bases along the path. Reflection is not universally a 180-degree reversal, nor do all room surfaces have the same coefficient. Surface roughness and objects can also produce diffuse scattering rather than one specular ray. See [MathWorks' reflection and polarization model](https://www.mathworks.com/help/antenna/ug/ray-tracing-for-wireless-communications.html).

For now, The Box accepts one configurable complex scalar coefficient for each surface. The demonstration uses 3% power reflectance and a pi-radian phase shift on every surface. **These are illustrative inputs, not measured or Fresnel-derived material properties.** They replace the old ambiguous `reflectance = 0.03`; that choice does not claim a physical interpretation for the original constant.

## 4. The Box model

### Geometry and image sources

The original idea is to approximate an indoor environment as a rectangular enclosure and predict spatial multipath from its surfaces. A rectangular box is sufficient; the three side lengths need not be equal.

The current model contains seven rays: one unobstructed line of sight and one single-bounce reflection from each of the four walls, floor, and ceiling. Coordinates use meters in an axis-aligned box from $(0,0,0)$ to $(L,W,H)$. The transmitter and receiver must be distinct and strictly inside it.

To find a reflection from a plane through $q$ with unit normal $\hat n$, mirror the transmitter $t$ across the plane:

$$
t'=t-2\hat n\,[\hat n\cdot(t-q)].
$$

The line from the image transmitter to receiver $r$ is:

$$
\ell(u)=t'+u(r-t').
$$

Its intersection with the reflecting plane occurs at:

$$
u=\frac{\hat n\cdot(q-t')}{\hat n\cdot(r-t')},\qquad s=\ell(u).
$$

Accept the reflection only if the intersection lies on the segment and within the finite surface. The unfolded path length equals:

$$
L_p=\|t-s\|+\|s-r\|=\|t'-r\|.
$$

`boxReflectionPaths.m` implements this for all six surfaces. Its output order is LOS, x=0, x=L, y=0, y=W, z=0, z=H. It works for arbitrary interior antenna positions, including unequal heights. The old midpoint projection happens to work on the original symmetric route, but is not valid for general positions.

### Path amplitudes and received power

Define $d_0$ as a reference distance and $P_{\mathrm{ref,dBm}}$ as the power of a single reference path. The scalar channel is:

$$
h=\sum_{p=0}^{6}\Gamma_p\left(\frac{d_0}{L_p}\right)^{n/2}
\exp\left(-j\frac{2\pi L_p}{\lambda}\right),
\qquad \Gamma_0=1.
$$

The amplitude exponent is $n/2$ because power is proportional to squared amplitude. Predicted received power is:

$$
P_{\mathrm{rx,dBm}}=P_{\mathrm{ref,dBm}}+20\log_{10}|h|.
$$

The default is $n=2$, representing free-space power spreading along each ray, with reflection attenuation handled by $\Gamma_p$. An empirical indoor path-loss exponent summarizes many propagation effects; using it on every ray while also modeling reflections risks counting some effects twice. Exact cancellation gives negative infinity dBm in the ideal equations; the simulation adds no artificial noise floor.

`fitLogDistanceModel.m` separately fits the measured trend:

$$
P_{\mathrm{trend,dBm}}(d)=A-10n\log_{10}(d/d_0).
$$

This is linear least squares in $\log_{10}(d/d_0)$ and requires no Statistics toolbox. Its inputs and reference output are now **signed dBm**, an intentional change from the original positive-magnitude convention. The fitted trend is plotted for comparison and does not set the ray amplitudes or their exponent.

### Historical experiments and calibration

| Workbook | Distance unit | Samples per position | Treatment |
|---|---|---|---|
| `dataCollection_9_10_22.xlsx` | Feet | Five per rotation, two rotations | Default simulation comparison |
| `dataCollection_10_9_22.xlsx` | Inches | Five per rotation, two rotations | Dense measurement plot |

Both contain a separate `1 meter` row. The loader detects units from the header, excludes that row from the distance sweep, and reads it as the reference measurement. The September simulation covers the complete measured distance range, 0.9144 to 8.2296 m, with a maximum plotting step of 5 mm. Fine simulation sampling does not increase experimental resolution.

Under the negative-magnitude decoding assumption, September's measured trend gives $A\approx-38.64$ dBm and $n\approx1.257$, whereas its separate 1 m row gives -33 dBm. Its rotation means differ by as much as 11 dB at the same distance. Keep rotations separate: changes in orientation can change antenna gain, polarization coupling, and the strength of individual paths.

The default ray model uses the separately recorded -33 dBm reference as a **provisional scale**, rather than fitting to the same multipath sweep being predicted. That reference may itself contain multipath, so it is not an isolated direct-path calibration. For quantitative validation, obtain an independent reference under controlled conditions or estimate model parameters on training measurements and assess predictions on held-out positions.

The September room is 12.7127 by 14.0462 by 3 m, with Tx at (2.2479, 7.8486, 1) m. Rx moves along positive x at the same y-coordinate and height. Height and carrier settings remain assumptions inherited from the prototype. The October plot does not automatically reuse the September room geometry: confirm that geometry before making a dense model-to-measurement comparison.

### Limits and a useful next experiment

The current model omits multiple-bounce paths, furniture, people, obstructions, transmission through walls, diffraction, diffuse scattering, antenna directionality, explicit polarization, receiver noise, and bandwidth averaging. Single-bounce rays are a starting approximation, not a complete room channel. Even exact arithmetic cannot compensate for unknown channel frequency or inaccurate antenna positions.

For another experiment, record the actual carrier frequency, radio RSSI conversion, antenna orientation and phase-center positions, room dimensions, surface materials, and channel activity. Repeat a dense spatial sweep with fixed orientation and repeat readings at each position. Compare repeatability and separate orientation sweeps before attributing mismatch to congestion. The September one-foot spacing exceeds two wavelengths; the October one-inch spacing is better suited to observing fine spatial fading. Model validation should distinguish an average power trend from accurate locations and depths of interference nulls.

## Implementation map and checks

| File | Purpose |
|---|---|
| `TheBox_SIMULATION.m` | Configure assumptions, sweep September geometry, plot predictions and measurements |
| `boxReflectionPaths.m` | Image-source single-bounce geometry |
| `combinePaths.m` | Coherent scalar amplitude sum and dBm conversion |
| `loadRecordedData.m` | Explicit RSSI decoding, units, raw samples, rotations, reference row |
| `fitLogDistanceModel.m` | Independent empirical trend using base MATLAB |
| `PlotRecorded_Data.m` | Plot October measurements by rotation |
| `testBoxModel.m` | Physical invariants, arbitrary geometry, fitting, and workbook regression checks |

Run the tests and both scripts from the repository directory. The numerical core and plotting use base MATLAB; reading `.xlsx` uses `readcell`. No RF or Statistics toolbox is required.

The repository also includes historical reading in [`References/`](References/), including Andrea Goldsmith's *Wireless Communications* (2005). The external documentation linked above supplies additional background; the equations here describe this implementation's stated scalar approximation.
