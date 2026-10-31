# An educational primer on RSSI and indoor multipath

An indoor radio signal can reach a receiver along several routes. One copy may travel directly, another may reflect from the floor, and others may reflect from walls or the ceiling. These copies can strengthen or weaken one another. Moving an antenna can therefore make the received signal weaker even when it moves closer to the transmitter.

The Box explores this behavior by representing a room as a rectangular enclosure. To understand that model, we first need to understand what a receiver measures, how waves combine, what reflections change, and how geometry determines a wave's phase.

This primer develops those ideas without requiring prior knowledge of radio propagation or complex numbers. The [code guide](CODE_GUIDE.md) explains how the MATLAB files implement them, including the historical measurements and configuration settings.

## Reading the notation

An equation is a compact description of a sequence of operations. Symbols will be introduced where they are used, but these conventions recur:

| Notation | Meaning |
|---|---|
| $x^2$ | Square of $x$: multiply $x$ by itself |
| $\sqrt{x}$ | Square root: the nonnegative number whose square is $x$ |
| $\log_{10}(x)$ | Base-10 logarithm: the exponent needed to make $10$ equal $x$ |
| $\sum_p x_p$ | Add the values of $x_p$ over the indicated path indices |
| $x_p$ | Value of $x$ for path number $p$; the subscript is a label |
| $\lvert z\rvert$ | Magnitude of a scalar or complex number |
| $\lVert\mathbf{x}\rVert$ | Length of a position or direction vector |
| $\mathbf{a}\cdot\mathbf{b}$ | Dot product of two vectors |
| $\propto$ | Is proportional to: equal up to a constant scale factor |
| $\approx$ | Is approximately equal to |

Quantities inside logarithms must be dimensionless ratios. A distance in meters divided by another distance in meters is dimensionless, as is a power divided by a reference power.

## 1. RSSI: what the receiver reports

### Power, amplitude, and the measurement

**Received Signal Strength Indicator**, abbreviated **RSSI**, is a radio receiver's estimate of received signal strength. Some radios report a relative index, while others report a calibrated estimate of power. A raw register value needs the manufacturer's conversion rule before it can be interpreted physically. See the [IEEE overview of RSSI](https://technav.ieee.org/topic/received-signal-strength-indicator/).

Power is the rate at which energy is transferred, measured in watts (W). One milliwatt (mW) is one thousandth of a watt. Received radio powers can be very small, so expressing them on a logarithmic scale is convenient.

Amplitude describes the size of an oscillating field or voltage. It is related to power, but is not itself power. For a sinusoidal voltage across a fixed resistance, doubling the amplitude quadruples the average power. This squared relationship will explain both interference and the factor of 20 in an equation later in this primer.

An RSSI reading also depends on the receiver's measurement bandwidth, averaging time, and design. It does not tell us how much signal arrived along each path. Depending on the measurement, other transmitters and noise can contribute to what is reported.

### From watts to dBm

The unit **dBm** expresses a power relative to 1 mW:

$$
P_{\mathrm{dBm}}=10\log_{10}\left(\frac{P}{1\,\mathrm{mW}}\right).
$$

Every part has a role:

- $P$ is the physical power being measured, with the same units as the reference in the denominator.
- $1\,\mathrm{mW}$ is the fixed reference power. Dividing by it produces a unitless ratio.
- $\log_{10}$ converts a multiplicative ratio into an additive number. For example, $\log_{10}(100)=2$ because $10^2=100$.
- The factor $10$ puts this logarithmic power ratio on the decibel scale.
- $P_{\mathrm{dBm}}$ is the resulting numerical level in dBm.

For 1 mW, the ratio is 1, its logarithm is 0, and the level is 0 dBm. For a weaker signal, the ratio is less than 1 and the logarithm is negative. **Negative dBm does not mean negative physical power.**

| Physical power | Ratio to 1 mW | Level |
|---|---:|---:|
| 1 mW | 1 | 0 dBm |
| 0.1 mW | 0.1 | -10 dBm |
| 0.0001 mW | 0.0001 | -40 dBm |
| 0.00001 mW | 0.00001 | -50 dBm |

Thus -40 dBm is stronger than -50 dBm: it represents ten times the power.

To undo the logarithm, raise 10 to the appropriate exponent:

$$
P=(1\,\mathrm{mW})\,10^{P_{\mathrm{dBm}}/10}.
$$

First divide the dBm reading by 10, then use the result as an exponent, then multiply by 1 mW. For -40 dBm this gives $10^{-4}$ mW, or 0.0001 mW.

### dB and dBm serve different purposes

A power ratio can be expressed in **dB**:

$$
G_{\mathrm{dB}}=10\log_{10}\left(\frac{P_2}{P_1}\right).
$$

$P_1$ is the starting power and $P_2$ is the resulting power. A ratio of 2 is about +3.01 dB; a ratio of 4 is about +6.02 dB. A ratio of 0.5 is about -3.01 dB.

A dBm level specifies absolute power relative to 1 mW. A dB difference specifies a ratio. You can add a gain in dB to a level in dBm to find another level in dBm. You cannot add two dBm levels to find their combined physical power.

### Averaging readings

Suppose two measurements are -40 and -50 dBm. Their arithmetic mean on the logarithmic scale is:

$$
\overline{P}_{\mathrm{log,dBm}}=\frac{1}{N}\sum_{i=1}^{N}P_{i,\mathrm{dBm}}.
$$

Here $N$ is the number of measurements, $i$ identifies each measurement, and the summation means to add all the readings before dividing by $N$. For our two readings, this is $(-40-50)/2=-45$ dBm. It describes the mean logarithmic level.

To find average physical power, convert each reading to mW first:

$$
\overline{P}_{\mathrm{mW}}=\frac{1}{N}\sum_{i=1}^{N}10^{P_{i,\mathrm{dBm}}/10},
\qquad
\overline{P}_{\mathrm{power,dBm}}=10\log_{10}(\overline{P}_{\mathrm{mW}}).
$$

In this notation $\overline{P}_{\mathrm{mW}}$ is the numerical power value in mW. The average is $(0.0001+0.00001)/2=0.000055$ mW, or approximately -42.60 dBm. This differs from -45 dBm because the logarithm is nonlinear. Neither statistic should silently be substituted for the other.

## 2. A radio wave has amplitude and phase

### Describing an oscillation

A simple received voltage can be written as:

$$
v(t)=V\cos(2\pi ft+\phi).
$$

The symbols describe a repeating oscillation:

- $v(t)$ is the instantaneous voltage at time $t$.
- $V$ is the peak voltage amplitude, a nonnegative number.
- $t$ is time in seconds.
- $f$ is frequency in hertz (Hz), meaning cycles per second.
- $2\pi$ converts cycles into radians. One full cycle is $2\pi$ radians, or 360 degrees.
- $\phi$ is the phase offset in radians: where the wave starts within its cycle.
- $\cos$ converts the angle into an oscillating value between -1 and 1.

At 2.4 GHz, $f=2.4\times10^9$ Hz: 2.4 billion cycles occur each second. The voltage changes sign during each cycle, but the average power delivered to a resistance is nonnegative.

### Why power depends on amplitude squared

For an ideal resistance $R_0$, instantaneous power is $v(t)^2/R_0$. The average of $\cos^2$ over a full cycle is $1/2$, so the average power of the sinusoid is:

$$
P=\frac{V^2}{2R_0}.
$$

$V$ is the peak amplitude and $R_0$ is resistance in ohms. If $V$ doubles while $R_0$ stays fixed, $V^2$ grows by a factor of four. The exact voltage-to-power conversion matters when modeling a circuit; for relative channel predictions, its constant factor can be absorbed into a reference power.

### Complex numbers as bookkeeping for phase

Keeping track of many shifted cosines becomes easier with complex numbers. Define $j$ by $j^2=-1$. A complex number has a real part and an imaginary part, which can be drawn as horizontal and vertical coordinates on a plane.

Euler's identity connects this plane to angles:

$$
e^{j\phi}=\cos\phi+j\sin\phi.
$$

$e$ is the base of the natural exponential. The value $e^{j\phi}$ lies on the unit circle at angle $\phi$. Its magnitude is 1, so multiplying by it rotates a complex number without changing its magnitude.

A **phasor** stores an oscillation's amplitude and phase together:

$$
a=Ve^{j\phi},\qquad
v(t)=\operatorname{Re}\{a e^{j2\pi ft}\}.
$$

$a$ is the phasor, $V$ is its magnitude, and $\phi$ is its angle. The operation $\operatorname{Re}$ takes the real part. Substituting Euler's identity gives the same cosine wave introduced above. The imaginary part is a mathematical way to retain the part of the phase information a single cosine value would lose.

For $a=x+jy$, its magnitude is $|a|=\sqrt{x^2+y^2}$. Its squared magnitude is $|a|^2=aa^*$, where the star means complex conjugation: replace $x+jy$ with $x-jy$. This identity is useful when deriving power.

The phasor is not a new physical signal being added to the radio. It is a representation of an ordinary oscillation, with a shared time dependence factored out.

## 3. Multipath: combine amplitudes before calculating power

### Superposition

In the linear propagation approximation, electric fields add at the receiver. Compatible received voltages therefore add too. With $M$ paths numbered $0$ through $M-1$:

$$
a_{\mathrm{total}}=\sum_{p=0}^{M-1}a_p.
$$

$p$ labels a path, and $a_p$ contains that path's amplitude and phase. The summation adds complex numbers, including their signs and angles. Only after that sum do we calculate received power:

$$
P_{\mathrm{total}}\propto|a_{\mathrm{total}}|^2.
$$

This is why neither adding dBm values nor simply adding each path's individual power describes coherent interference.

### Deriving the two-path equation

Consider two arrivals $a_0=V_0e^{j\phi_0}$ and $a_1=V_1e^{j\phi_1}$. Expanding the squared magnitude gives:

$$
\begin{aligned}
|a_0+a_1|^2
&=(a_0+a_1)(a_0^*+a_1^*)\\
&=|a_0|^2+|a_1|^2+a_0a_1^*+a_0^*a_1\\
&=V_0^2+V_1^2+2V_0V_1\cos(\phi_1-\phi_0).
\end{aligned}
$$

The first two terms are the individual squared amplitudes. The last term is the **interference term**: it depends on their phase difference. The two cross-products combine into a cosine because they are complex conjugates.

Using the same voltage-to-power factor for both arrivals turns this into:

$$
P_{\mathrm{total}}=P_0+P_1+2\sqrt{P_0P_1}\cos(\Delta\phi),
\qquad \Delta\phi=\phi_1-\phi_0.
$$

$P_0$ and $P_1$ are the powers each arrival would produce alone, in linear units such as mW. Their geometric mean $\sqrt{P_0P_1}$ sets the size of the interference term. The cosine determines its sign.

For two equal-power arrivals, let each have power $P_s$:

| Phase difference | Cosine | Combined power | Interpretation |
|---|---:|---:|---|
| $0$ radians (0 degrees) | 1 | $4P_s$ | Fully constructive |
| $\pi/2$ radians (90 degrees) | 0 | $2P_s$ | Cross-term is zero |
| $\pi$ radians (180 degrees) | -1 | 0 | Fully destructive in the ideal model |

The constructive case is +6.02 dB relative to either arrival alone. Perfect cancellation requires equal amplitudes as well as opposite phases. If amplitudes differ, the smallest possible squared amplitude is $(V_0-V_1)^2$.

There is no conservation-of-energy contradiction when two arrivals add constructively at a location. Interference redistributes energy across space; other locations can have lower field strength.

With seven paths, every pair can contribute an interference term. Summing all seven complex amplitudes automatically includes those relationships.

## 4. Travel distance creates delay and phase shift

A wave takes time to travel. If a path has length $L$:

$$
\tau=\frac{L}{c},\qquad \lambda=\frac{c}{f}.
$$

$\tau$ is travel time in seconds, $L$ is distance in meters, $c\approx299{,}792{,}458$ m/s is the propagation speed used here, and $\lambda$ is wavelength in meters. Indoor air is approximated as having the vacuum speed of light. The wavelength is the distance traveled during one cycle.

A delayed sinusoid is evaluated at $t-\tau$:

$$
\cos(2\pi f(t-\tau)+\phi)
=\cos(2\pi ft+\phi-2\pi f\tau).
$$

Expanding the parentheses shows the extra phase is negative under our time convention. Substituting $\tau=L/c$ and $\lambda=c/f$ gives:

$$
\phi_{\mathrm{prop}}=-2\pi f\tau=-\frac{2\pi L}{\lambda}=-kL,
\qquad k=\frac{2\pi}{\lambda}.
$$

$k$ is the **wavenumber**, measured in radians per meter. Multiplying it by meters produces a phase in radians. In phasor form, propagation therefore contributes $e^{-jkL}$, a unit-magnitude rotation. Path attenuation is a separate factor.

For two paths, their propagation phase difference is $-k(L_1-L_0)$. Increasing the difference by one wavelength adds a full rotation, returning to the same relative phase. Increasing it by half a wavelength adds 180 degrees.

At 2.400 GHz:

$$
\lambda\approx0.1249\,\mathrm{m},\qquad
|\Delta\phi|=\frac{2\pi(0.01)}{0.1249}\approx0.503\,\mathrm{rad}\approx28.8^\circ.
$$

This calculation describes a 1 cm change in **relative path length**, not necessarily a 1 cm movement of the receiver. Receiver movement changes different path lengths by different amounts.

Only relative phase determines the power sum. Rotating every phasor by the same angle rotates the total phasor but preserves its length. This allows us to choose one arrival as the phase reference.

## 5. Reflections change amplitude, phase, and direction

### Specular reflection and polarization

For a smooth planar surface, a specular reflection behaves geometrically like a mirror. The incoming ray, outgoing ray, and surface normal lie in one plane. The incidence angle equals the reflection angle, with both measured from the normal. A **normal** is a direction perpendicular to the surface.

Real electromagnetic waves also have **polarization**, the orientation and motion of their electric field. Surface response can differ for electric-field components parallel and perpendicular to the plane of incidence. Material properties, angle, frequency, and polarization determine the reflected field. Fresnel equations describe these components for ideal planar interfaces. A scalar model approximates this behavior with one effective coefficient per ray; it does not track the full field orientation. See [MathWorks' reflection and polarization explanation](https://www.mathworks.com/help/antenna/ug/ray-tracing-for-wireless-communications.html).

### The complex reflection coefficient

Represent a reflection by:

$$
\Gamma=\rho e^{j\theta},\qquad a_{\mathrm{after}}=\Gamma a_{\mathrm{before}}.
$$

$\Gamma$ (capital gamma) is the field reflection coefficient. $\rho=|\Gamma|$ is its magnitude, and $\theta$ is its phase shift in radians. The before-and-after amplitudes refer to the reflection interaction; travel to and from the surface must still be accounted for.

If the incident phasor is $Ve^{j\phi}$, multiplication gives:

$$
a_{\mathrm{after}}=\rho V e^{j(\phi+\theta)}.
$$

This explicitly shows both changes: amplitude becomes $\rho V$ and phase becomes $\phi+\theta$. A positive real coefficient introduces no phase change; a negative real coefficient introduces a 180-degree reversal. An arbitrary complex coefficient can introduce another angle. A reflection does not universally reverse the phase.

For a reflected wave returning into the same medium, the power reflectance is:

$$
R=|\Gamma|^2=\rho^2,\qquad \rho=\sqrt{R}.
$$

Here $R$ is a dimensionless power fraction, distinct from the resistance $R_0$ used earlier. If $R=0.03$, then $\rho\approx0.1732$. With a 180-degree phase shift, $\Gamma\approx-0.1732$. Using 0.03 as the field magnitude would instead imply a power fraction of $0.03^2=0.0009$, or 0.09%.

In a path with one reflection, a useful decomposition is:

$$
a_p=(\text{reference amplitude})\,(\text{spreading factor})\,
\underbrace{\rho_p e^{j\theta_p}}_{\text{reflection}}\,
\underbrace{e^{-jkL_p}}_{\text{propagation phase}}.
$$

The total path phase contains both $\theta_p$ and $-kL_p$. Reflection phase cannot be replaced by propagation phase or ignored when predicting nulls.

## 6. Distance reduces the amplitude of each path

An ideal source spreading a fixed total power uniformly in three dimensions distributes that power over a sphere. A sphere's area grows as $4\pi d^2$, so power density decreases as $1/d^2$. With fixed antenna properties, received free-space power follows this inverse-square dependence in the far field.

A more general power-distance relationship is:

$$
P(d)=P_{\mathrm{ref}}\left(\frac{d_0}{d}\right)^n.
$$

$P_{\mathrm{ref}}$ is the linear received power at reference distance $d_0$. $d$ is the new distance. $n$ is the **power path-loss exponent**. With $n=2$, doubling distance leaves one quarter of the power. The reference must apply in the regime where this model is valid; the expression is not a physical description arbitrarily close to an antenna.

Since power is proportional to amplitude squared, taking a square root gives the amplitude ratio:

$$
\frac{V(d)}{V(d_0)}=\left(\frac{d_0}{d}\right)^{n/2}.
$$

This is the origin of $n/2$ in the channel equation. For $n=2$, amplitude decreases as $1/d$, while power decreases as $1/d^2$.

Converting the power law to dBm gives:

$$
P_{\mathrm{dBm}}(d)=P_{\mathrm{ref,dBm}}-10n\log_{10}(d/d_0).
$$

The minus sign comes from exchanging $d_0/d$ for $d/d_0$ inside the logarithm. With $n=2$, doubling distance changes the level by $-20\log_{10}(2)\approx-6.02$ dB.

A measured indoor exponent can summarize reflection, obstruction, and other effects together. It need not equal 2. Using that aggregate exponent on each geometric ray while explicitly adding reflections can count some effects twice. The Box therefore starts with free-space spreading on each ray and treats a fitted measurement trend separately.

## 7. The Box: use room geometry to find seven paths

### Coordinates and distances

Approximate the room by six flat surfaces, with an interior extending from $(0,0,0)$ to $(D_x,D_y,D_z)$. The side lengths need not be equal. Use $D_x,D_y,D_z$ for room dimensions so that $L_p$ can consistently mean a path length.

A position is a three-component vector in meters. Let $\mathbf{t}$ be the transmitter position and $\mathbf{r}$ the receiver position. Their direct distance is:

$$
L_0=\lVert\mathbf{r}-\mathbf{t}\rVert
=\sqrt{(r_x-t_x)^2+(r_y-t_y)^2+(r_z-t_z)^2}.
$$

Subtracting coordinates gives the displacement along each axis. Squaring, adding, and taking the square root is the three-dimensional Pythagorean theorem.

The model keeps one direct path and six paths with exactly one reflection: one from each wall, one from the floor, and one from the ceiling. It assumes these routes are unobstructed.

### Mirror the transmitter across a surface

To locate a physically valid reflection, imagine mirroring the transmitter through the surface. A straight line from that image to the receiver unfolds the reflected path.

Describe the reflecting plane using a point $\mathbf{q}$ on it and a unit normal $\hat{\mathbf{n}}$. A unit normal has length 1. For the wall $x=0$, one choice is $\mathbf{q}=(0,0,0)$ and $\hat{\mathbf{n}}=(1,0,0)$.

For vectors, the dot product is:

$$
\mathbf{a}\cdot\mathbf{b}=a_xb_x+a_yb_y+a_zb_z.
$$

When one vector is a unit direction, the dot product measures the signed component of the other along that direction. Therefore:

$$
\delta=\hat{\mathbf{n}}\cdot(\mathbf{t}-\mathbf{q})
$$

is the transmitter's signed perpendicular distance from the plane. The displacement $\delta\hat{\mathbf{n}}$ takes us from the plane toward the transmitter. Subtract it once to reach the plane and twice to reach the mirror position:

$$
\mathbf{t}'=\mathbf{t}-2\delta\hat{\mathbf{n}}.
$$

For an axis-aligned wall at coordinate $x=w$, this simply replaces $t_x$ with $2w-t_x$, leaving the other coordinates unchanged.

### Find the reflection point

A line from image transmitter $\mathbf{t}'$ toward receiver $\mathbf{r}$ can be written:

$$
\boldsymbol{\ell}(u)=\mathbf{t}'+u(\mathbf{r}-\mathbf{t}').
$$

$u$ is a unitless fraction along the line. At $u=0$ we are at the image, at $u=1$ at the receiver, and at intermediate values between them. A point on the reflecting plane has zero normal displacement from $\mathbf{q}$, so solve:

$$
\hat{\mathbf{n}}\cdot(\boldsymbol{\ell}(u)-\mathbf{q})=0.
$$

Substitute the line expression and isolate $u$:

$$
u=\frac{\hat{\mathbf{n}}\cdot(\mathbf{q}-\mathbf{t}')}
{\hat{\mathbf{n}}\cdot(\mathbf{r}-\mathbf{t}')},
\qquad \mathbf{s}=\boldsymbol{\ell}(u).
$$

The numerator is the normal displacement needed to reach the plane. The denominator is the normal displacement across the complete image-to-receiver line. Their ratio tells us what fraction of that line reaches the surface. $\mathbf{s}$ is the reflection point.

If the denominator is zero, the line is parallel to the plane and this formula cannot produce a unique intersection. A candidate must also have $0<u<1$ and lie within the actual finite wall, rather than elsewhere on its infinite plane.

The reflected length is:

$$
L_p=\lVert\mathbf{s}-\mathbf{t}\rVert+\lVert\mathbf{r}-\mathbf{s}\rVert
=\lVert\mathbf{r}-\mathbf{t}'\rVert.
$$

The mirror construction preserves the first segment's length and unfolds the two segments into a straight line. It also enforces equal incidence and reflection angles.

### A worked reflection example

Let the reflecting wall be $x=0$, transmitter $\mathbf{t}=(2,2,1)$ m, and receiver $\mathbf{r}=(4,5,1)$ m. Assume the wall is large enough to contain the reflection point.

The image is $\mathbf{t}'=(-2,2,1)$ m. The direction from image to receiver is $(6,3,0)$ m. To reach $x=0$ from $x=-2$ requires 2 m of the available 6 m x-displacement, so $u=2/6=1/3$. Thus:

$$
\mathbf{s}=(-2,2,1)+\frac{1}{3}(6,3,0)=(0,3,1)\,\mathrm{m}.
$$

The reflected length is $\sqrt{6^2+3^2}=\sqrt{45}\approx6.708$ m. The direct length is $\sqrt{2^2+3^2}=\sqrt{13}\approx3.606$ m. Their length difference, together with the reflection phase, determines their relative phase at the receiver.

Projecting the transmitter-receiver midpoint onto a wall would not generally produce this point. The image method accounts for unequal distances from the transmitter and receiver to the wall.

## 8. Assemble the channel and predict power

For each path, combine spreading, reflection, and propagation phase into a **dimensionless relative amplitude**:

$$
h_p=\Gamma_p\left(\frac{d_0}{L_p}\right)^{n/2}e^{-jkL_p}.
$$

Here is what each factor contributes:

| Factor | Meaning |
|---|---|
| $p$ | Path index; 0 identifies the direct path |
| $L_p$ | Complete length of path $p$, in meters |
| $d_0$ | Reference path length, in meters |
| $(d_0/L_p)^{n/2}$ | Amplitude change due to distance |
| $\Gamma_p$ | Complex field coefficient for the reflection; 1 for the direct path |
| $k=2\pi/\lambda$ | Propagation phase per meter |
| $e^{-jkL_p}$ | Propagation phase rotation for that path |
| $h_p$ | Path amplitude relative to the reference amplitude |

The seven-path channel is their sum:

$$
h=\sum_{p=0}^{6}h_p
=\sum_{p=0}^{6}\Gamma_p\left(\frac{d_0}{L_p}\right)^{n/2}e^{-jkL_p}.
$$

$h$ is generally complex. Its magnitude tells us the total amplitude relative to a reference single path; its angle gives the total phase relative to our chosen clock. Because power follows amplitude squared:

$$
P_{\mathrm{rx}}=P_{\mathrm{ref}}|h|^2.
$$

$P_{\mathrm{ref}}$ must be a linear power calibrated for the reference single path. Convert to dBm using the logarithm of a product:

$$
\begin{aligned}
P_{\mathrm{rx,dBm}}
&=P_{\mathrm{ref,dBm}}+10\log_{10}(|h|^2)\\
&=P_{\mathrm{ref,dBm}}+20\log_{10}|h|.
\end{aligned}
$$

The factor 20 comes from $\log_{10}(x^2)=2\log_{10}(x)$. It applies to the amplitude ratio; a power ratio would use 10.

For example, if the combined normalized amplitude has magnitude 0.1 and the reference is -33 dBm, the result is $-33+20\log_{10}(0.1)=-53$ dBm. If the magnitude is 2, it is about -26.98 dBm.

A convenient equivalent expression removes the direct path's common phase:

$$
\widetilde{h}=\sum_{p=0}^{6}\Gamma_p\left(\frac{d_0}{L_p}\right)^{n/2}
 e^{-jk(L_p-L_0)}.
$$

Since $\widetilde{h}=e^{jkL_0}h$ and $|e^{jkL_0}|=1$, both give the same received power. The direct term now has zero propagation phase. This is a change of phase reference, not removal of interference.

## 9. Compare a geometric prediction with an empirical trend

A smooth log-distance trend summarizes how measured level varies with distance:

$$
P_{\mathrm{trend,dBm}}(d)=A-10n_{\mathrm{trend}}\log_{10}(d/d_0).
$$

$A$ is the fitted level at $d_0$, because $\log_{10}(1)=0$. $n_{\mathrm{trend}}$ is a fitted power-distance exponent. The subscript distinguishes it from a per-ray spreading exponent. The trend does not specify wall reflection points or predict individual interference nulls.

To see how fitting works, define one explanatory value per measurement:

$$
x_i=-10\log_{10}(d_i/d_0),\qquad y_i=A+n_{\mathrm{trend}}x_i+\epsilon_i.
$$

$d_i$ is measured distance, $y_i$ is measured level in signed dBm, and $\epsilon_i$ is the discrepancy between the line and the measurement. Although distance enters through a logarithm, this expression is linear in the two unknown parameters $A$ and $n_{\mathrm{trend}}$.

Least squares chooses those parameters to minimize:

$$
\sum_{i=1}^{N}\left[y_i-(A+n_{\mathrm{trend}}x_i)\right]^2.
$$

Each bracket is a residual in dB. Squaring prevents positive and negative residuals from canceling; summing measures total squared discrepancy. It is a descriptive fit, not proof that every discrepancy is random noise.

A reference reading taken inside a room can already contain multipath. Using it as an isolated direct-path calibration introduces uncertainty. Likewise, fitting a trend to a multipath sweep and then adding explicit reflections to that fitted curve mixes two descriptions of the same effects. Prefer independently calibrated path scaling, or estimate parameters on training measurements and evaluate on separate held-out measurements.

## 10. Where the approximation ends

The seven-ray calculation is **narrowband**: it evaluates interference at one carrier frequency. A modulated signal contains a range of frequencies. Different paths delay its changing envelope as well as its carrier.

A more general complex-envelope model is:

$$
y(t)=\sum_{p=0}^{M-1}g_p\,x(t-\tau_p)+\eta(t).
$$

$x(t)$ is the transmitted envelope, $x(t-\tau_p)$ is a delayed copy, $g_p$ is its complex gain including carrier phase, $y(t)$ is the received envelope, and $\eta(t)$ represents noise and unmodeled interference. When the envelope changes negligibly over the differences between path delays, the delayed envelopes can be approximated by a common envelope and the scalar phasor sum becomes useful. See [MathWorks' channel formulation](https://www.mathworks.com/help/comm/ref/comm.raytracingchannel-system-object.html).

One useful dimensionless check is $B\Delta\tau$, where $B$ is occupied bandwidth in Hz and $\Delta\tau$ is the difference between relevant path delays in seconds. Values much smaller than 1 support the common-envelope approximation. Across larger bandwidths, phase relationships vary with frequency, and a single-frequency null need not describe the receiver's bandwidth-averaged power.

Independent transmitters do not normally maintain a fixed phase relationship with the desired signal. Their contributions cannot simply be added as permanent coherent reflection paths. Congestion may affect readings and packet reception, but RSSI alone cannot establish it as the cause of a model mismatch.

A real room also contains furniture, people, rough surfaces, obstructions, multiple reflections, diffraction around edges, and transmission through materials. Antennas have direction-dependent gain and polarization. The Box's scalar, single-bounce model omits these effects. Its value is that each retained effect is explicit and inspectable.

A useful experiment fixes antenna orientation, records accurate positions and actual channel frequency, repeats measurements at each position, and samples densely enough to observe spatial fading. Repeating the same route reveals temporal variation; repeating with another orientation reveals orientation-associated changes. A quiet channel helps control one source of variation, while accurate geometry and calibration remain necessary for predicting specific peaks and nulls.

## Further reading

The linked references provide additional context. The repository's [`References/`](References/) directory also contains historical reading, including Andrea Goldsmith's *Wireless Communications* (2005). For the exact implemented assumptions, parameter values, and MATLAB data flow, continue with the [code guide](CODE_GUIDE.md).
