# NACA TN-2474 — Full Paper Transcription

**Title:** Empirical Relation Between Induced Velocity, Thrust, and Rate of Descent of a
Helicopter Rotor as Determined by Wind-Tunnel Tests on Four Model Rotors

**Authors:** Walter Castles, Jr. and Robin B. Gray
**Institution:** Georgia Institute of Technology (under NACA sponsorship)
**Date:** September 1951 (received May 31, 1950)

> Transcribed from page_01.png – page_31.png.  Mathematical equations are rendered in
> plain text; for the original typeset forms see the source images.  Section and equation
> numbers follow the paper exactly.

---

## Summary

*(See abstract.md for the full summary text.)*

---

## Introduction

Contemporary vortex or momentum theory does not yield a useful answer for the mean induced
velocity of a helicopter rotor in vertical descent.  Thus the relation between the mean
induced velocity, thrust, and rate of descent must, at present, be determined
experimentally.

Both wind-tunnel and flight tests, the only readily available experimental methods of
determining the empirical induced-velocity relation for vertical descent, offer certain
difficulties.  It is not easy to maintain the desired zero horizontal component of velocity
on flight tests or to measure the rate of descent with sufficient accuracy.  Also, existing
helicopters tend to suffer a loss of control at the larger rates of descent; consequently,
it is hard to obtain steady-state flight data in this range.  Wind-tunnel model rotor tests,
on the other hand, present certain problems in determining the equivalent free-stream descent
velocities corresponding to the test conditions and also in measuring or deducing the
operating rotor blade angles.

Recent flight-test data have been in serious disagreement with the only previous
wind-tunnel test data, and Glauert's empirical induced-velocity relation, 1/f against 1/F,
based on the results of this test.  This disagreement has cast certain doubts on the
usefulness of wind-tunnel model rotor tests in general.

The present test program was undertaken in an effort to check Glauert's curve of 1/f
against 1/F over the useful range of vertical descent and, if possible, to find the sources
of discrepancy between the model and free-flight results.  It was also desired to evaluate
the principal effects of blade taper and twist on the descent characteristics and to obtain
a sequence of smoke filament and tuft photographs to show the approximate flow patterns.

A reexamination of the test methods and the procedure used to reduce the data for
Glauert's curve of 1/f against 1/F indicated that the probable sources of error were the
use of a wind tunnel having a closed test section and an energy ratio considerably greater
than unity, the determination of the equivalent free-stream descent velocity by means of a
static-pressure wall tap, and the neglect of the dynamic twist of the blades.  The present
tests incorporate a correction for all these errors.

---

## Symbols

| Symbol | Definition |
|--------|-----------|
| a | blade-element lift-curve slope |
| A | rotor-disk area |
| b | number of blades in rotor |
| C_d_e | blade-element profile drag coefficient |
| c_e | equivalent blade chord: c_e = (b/πR²) ∫₀ᴿ c dr |
| σ_e | effective solidity: σ_e = bc_e/(πR) = (b/πR²) ∫₀ᴿ c dr |
| c₁ | blade-element lift coefficient |
| c₀ | extended blade root chord |
| C_T | thrust coefficient: T/(ρπR²(ΩR)²) |
| C_Q | torque coefficient: Q/(ρπR³(ΩR)²) |
| ΔC_Q | increment in torque coefficient over that for zero thrust and zero rate of descent |
| f | thrust coefficient based on descent velocity: T/(ρA V²) |
| F | thrust coefficient based on resultant velocity at rotor |
| G | rotor torque |
| P | radius of blade element |
| R | rotor radius |
| R_e | effective radius: R_e = R − ½ c_tip |
| t | taper factor: (c_tip/c₀ − 1) |
| T | rotor thrust |
| V | descent velocity |
| V_i | induced velocity at rotor, measured with respect to fixed coordinates |
| x | nondimensional radius of blade element (r/R) |
| θ₀, θ₁, θ₂ | coefficients in power series representing θ_e as function of x |
| θ_e | blade-element pitch angle |
| θ₀.₇₅R | extended linear blade root pitch angle |
| θ_tip | extended linear blade tip pitch angle |
| λ₁ | nondimensional induced velocity: V_i / √(T/(2ρA)) |
| λ₂ | nondimensional descent velocity: V / √(T/(2ρA)) |
| ρ | mass density of air |
| σ₁…σ₇ | solidity factors (weighted integrals — see appendix) |
| φ | inflow angle at blade element, positive below rotor plane |
| ω | angular velocity of rotor |

---

## Description of Apparatus

### Wind tunnel

The wind tunnel, in the open jet of which these tests were conducted, is of conventional
closed-return type having a 5-to-1 contraction in the collector.  The open jet is circular,
of 9-foot diameter, and 10 feet long.  There is a flat wake-diffusion screen cover just
below the bottom boundary of the open jet.  The rotors were installed so that the rotor
hub was in the transverse and longitudinal center line of the open jet.

For the present tests an 18- by 18-mesh wake-diffusion screen was installed 4 feet
downstream of the rotor plane.  This screen reduced the tunnel energy ratio to
approximately 0.7.  In addition, a precision electric tachometer was installed on the
wind-tunnel motor and gear driven to read 10 times wind-tunnel propeller shaft speed.
This permitted the wind-tunnel propeller speed to be read to 2½ revolutions per minute.

### Rotor test stand

The rotor test stand is a self-contained unit which mounts in a load member support socket
of the normal wind-tunnel balance as shown in figure 1.  The power is furnished by a
three-phase, 15-horsepower, 1750-revolution-per-minute, wound rotor induction motor.  The
motor shaft is coupled to the pinion in the gear box contained in the 6-inch-pipe tee by
a vertical drive shaft.  The pinion drives, at half motor speed, a ring gear to which the
differential planetary carrier is attached as shown in figure 2.  The reaction from the
miter gear set in the differential carrier is restrained from rotating by the torque-measuring
strain-gauge arm.  The front shaft from the miter gear set drives the separately supported
rotor hub through a ball-bearing slip and universal joint in the center of the rotor.  This
ball-bearing slip and universal joint can transmit no thrust nor steady moments to the hub
other than the driving torque.

The hub is supported and rotates on two ball bearings mounted in a three-arm strain-gauge
spider, the outer ends of the arms of which are ball-bearing mounted in the test stand.
Inside the front of the hub there is a pitch-change motor which screws the outer hub fairing
back and forth and through ball-bearing connections changes the pitch of the rotor blades.
Extending from the front of this hub fairing is a revolution counter driven by the
pitch-change unit.  The blade-angle reading of this counter, which was read during the
tests by means of a stroboscope unit and field glasses, corresponds to a rotor-blade-angle
change of 0.04° to 0.08° over the range of blade angles covered in these tests.

### Gages and instrumentation

The thrust and torque strain gages were read by means of SR-4 bridges arranged in
temperature-compensating circuits.  The thrust could be read to approximately ±0.05
pound.  The torque on the heavier arm used with the 6-foot-diameter rotors could be read
to approximately ±0.004 foot-pound and on the lighter arm used with the 4-foot-diameter
rotor, to approximately ±0.002 foot-pound.

The rotor speed was measured by means of a neon flash lamp actuated by a set of breaker
points on the motor shaft.  This flash lamp illuminated once each rotor revolution a
suitably lined disk driven by a small synchronous motor.  The rotor speed was continuously
maintained constant by means of a three-phase lye barrel rheostat in the motor armature
circuit in such a manner that the line of the lined disk on the synchronous motor shaft
remained approximately stationary.  Assuming that power line frequency was constant, the
error in rotor speed was probably not over 12 revolutions per minute.  In order to be
certain that the desired image harmonic was the one being controlled, the rotor-drive-motor
speed was also read on an electric tachometer.

The tunnel-off fare thrust was not measurable.  The fare torque at 1600 revolutions per
minute, including windage on the hub and blade pockets, was 0.18 foot-pound and constant
after the gears had been run in.  In addition, there was a torque loss of approximately
0.10 foot-pound through the differential gear set.  The torque loss in the ring and pinion
was not reflected in the torque-arm strain-gage reading.

### Model rotor blades

The model rotor blades, each set of which had an effective solidity of 0.05 and NACA 0015
blade airfoil sections, were of the "rigid" type with no initial coning angle.  To keep
torsional deflections to a minimum the blades were designed so that the chordwise locations
of the centers of gravity, elastic axes, and aerodynamic centers of the blade elements
coincided, approximately, and lay on the quarter-chord pitch-change axis.

The large centrifugal loads arising from the full-scale design tip speed and the
above-mentioned torsional considerations necessitated building the blades with a solid
alloy-steel leading edge with the remaining back to approximately the quarter-chord point.

The constant-chord, untwisted blades were constructed with a hollow magnesium
trailing-edge section riveted to the steel leading edge as shown in figure 3.  The twisted
set of blades and the tapered blades were constructed with solid laminated mahogany
trailing-edge sections fastened to the steel leading edge with machine screws.  The blades
were hand-worked to contour, and the final finish obtained, while not aerodynamically smooth,
had no significant imperfections.

The diameter of the hub fairing was 6½ inches.  The blades were of true contour from a
radius of 6.09 inches out, and the blade tips were finished off square.  Rotor static
balance was obtained by ring balance weights held in position around the hub blade sockets
by set screws.

---

## Test Procedure

After the installation of the rotor test stand with hub, but without blades, and the tunnel
propeller and tunnel piezometer reading were calibrated against the tunnel jet velocity as
measured by a standard pitot tube and micromanometer.  Vertical and horizontal jet-centerline
velocity surveys were taken in the plane of rotation.  It was found that the jet velocity
distribution was satisfactory except for a local region directly ahead of the 5-inch-diameter
rotor test-stand support located 18 inches downstream from the plane of rotation.  Within
this local region the velocities were 3½ percent below the average.

The procedure for a typical run was as follows: After sufficient warm-up time for the
lubricating oil pressure to stabilize, tare thrust and torque readings were obtained at the
hub without blades at various rotor speeds and wind-tunnel velocities.  The hub was then
dismounted, a set of blades installed, and the rotor balanced.  The blade angle at the
three-quarter-radius point was set by placing them with a precision inclinometer to within
approximately ±2 minutes of equal angles of about 6° by adjusting the clamps on the
blade-pitch-change arms.  The hub was then reinstalled and the previous calibration of the
thrust strain system was checked by wire, pulley, and weights.

The rotor was then brought up to test speed (1200 or 1600 rpm) and, after another
warm-up period, a reading was taken of the reference blade angle and torque for zero thrust.
The accuracy with which the zero-thrust blade angle could be determined was not too
satisfactory, as explained below.

The blade angle was then set by trial and error at the value giving the desired thrust
coefficient as indicated by the thrust strain-gage setup, and a reading was taken of the
torque and blade angle.

The tunnel fan was then started and the above procedure and readings were repeated for
each successive increment in tunnel velocity that could be obtained from the taps on the
wind-tunnel-motor-armature rheostat.  In addition, a reading of the wind-tunnel propeller
speed was obtained at each of these points.  The run was terminated at that velocity
increment at which the measured torque reached a zero or negative value or, on the small
rotor, at a rate of descent known to be in the windmill brake state.

After a descent run at each desired value of C_T was obtained on a given rotor, the tunnel
was completely blocked by placing a layer of paper over the wake-diffusion screen.  A
hovering test run, reading thrust and torque against blade angle, was then made as a check
on the hovering points obtained in the vertical-descent runs.

At the larger rates of power-on descent the thrust and torque fluctuated in an irregular
manner.  An attempt was made in each such case to read the average values.

A chordwise bending fatigue failure occurred on one of the twisted blades while operating
at C_T ≈ 0.004 at 1600 revolutions per minute and a large rate of power-on descent.  Thus
the hovering check run was not obtained on these blades.

---

## Reduction of Data

As previously mentioned, a certain difficulty was experienced in obtaining the reference
blade angle for zero thrust.  Each of the following available methods appeared likely to
introduce certain errors:

1. Assuming the thrust was zero when the calibrated blade angle at the three-quarter-radius
   point was zero (for untwisted blades)
2. Assuming the calibrated zero point of the thrust strain-gage setup was the same with the
   rotor stationary and rotating
3. Assuming the thrust was zero when the tufts on the wires in the vicinity of the rotor
   were undisturbed (for untwisted blades)

In the first case, inaccuracies in the construction of the blades and subsequent warpage due
to operating stresses were likely to introduce appreciable error.  In the second case, the
accuracy is uncertain because of the impossibility of checking the zero-thrust calibration
point on the rotating rotor with the blades installed.  Although the calibration factor on
the thrust strain-gage spider remained constant during the period of the tests, there were
many small zero shifts.  Because of the low slope of the hovering values of C_T against θ
near the zero-thrust point, the small zero shifts in the thrust could have been translated
into appreciable changes in the zero-thrust blade angle.  In the third case, induced
velocities of the order of ±2 feet per second or less, equivalent to a zero-thrust blade-angle
shift of approximately ±0.3° at the three-quarter-radius point at 1200 revolutions per minute,
could not be detected by the tufts.

In general, method two was assumed to give the correct result unless shown to be obviously
in error by method one or three.

### Determination of operating blade angle

In order to present the final data in useful form, it was necessary to determine the
operating blade angle at the three-quarter-radius point on the rotating blades.  As
previously mentioned, the blades had a symmetrical airfoil section and were designed so
that the blade-element elastic axes and aerodynamic centers were very nearly coincident.
Also, the calculated chordwise deflections of the blades under the applied torques were
very small.  Therefore, the theoretical twist due to the air forces acting on the blades
should have been negligible and it was thus assumed that the only torque acting to twist
the blades was that arising from the inclination of the principal axis of inertia of the
blade sections to the plane of rotation.

The resulting dynamic torque was calculated as a function of the blade angle and rotor
angular velocity for each of the rotors.  The spring constants of each blade were then
measured experimentally (at three stations on the tapered blades), and the spring constant
was then corrected for each blade.  The dynamic twist between the hub and the three-quarter-
radius point was then calculated using the calculated dynamic-torque-loading curve and the
experimentally determined spring-constant curve.  Over the range of blade angles of these
tests the dynamic twist was very nearly a linear function of the blade angle, and the
operating blade angle at the three-quarter-radius point for the various rotors was given
with sufficient accuracy by the expressions:

For the 6-foot-diameter rotor with constant-chord, untwisted blades:
- θ₀.₇₅R = 0.963θ_root at 1600 rpm
- θ₀.₇₅R = 0.978θ_root at 1200 rpm

For the 4-foot-diameter rotor with constant-chord, untwisted blades:
- θ₀.₇₅R = 0.960θ_root at 1600 rpm
- θ₀.₇₅R = 0.965θ_root at 1200 rpm

For the rotor with untwisted tapered blades:
- θ₀.₇₅R = 0.936(θ_root − 7.79) at 1600 rpm
- θ₀.₇₅R = 0.964(θ_root − 7.79) at 1200 rpm

For the rotor with twisted constant-chord blades.  The deflections of the pitch-change
linkage were negligible.

The operating blade angle at the three-quarter-radius point was thus found by subtracting
the blade pitch-counter reading for zero thrust from the blade pitch-counter reading for
the test point in question, reading the equivalent blade root angle from the calibration
curve, and converting this root angle to the value at the three-quarter-radius point by
means of the appropriate equation above.

As a result of the dynamic twist, the actual twist of the rotating blade was slightly
different for each test point.  However, all comparisons have been made on the basis of
the initial static blade twist.

### Determination of equivalent free-stream descent velocity

The equivalent free-stream descent velocity was obtained from the calibration curve of
wind-tunnel jet velocity against wind-tunnel propeller speed, and the aerodynamic velocity
correction for dynamic blade twist was negligible.

In the absence of any applicable theory for, or useful measurements of, the flight range
covered in these tests, it was necessary to make the assumption that the induced velocity
was uniform over the rotor in order to calculate the values of λ₁ against λ₂ or 1/f
against 1/F from the test data.

### Equations used to reduce the data

As a result of this necessary supposition it followed from the assumed geometry that the
inflow angle β was a good approximation to take the inflow angle β as a small angle and
consider all blade elements as untwisted.  Thus the thrust can be written as

    T = ½ ρbc σ₁² r² dr  ... (1)

and the torque, as

    Q = ½ ρbc(c_d_e − a·φ) r² dr  ... (2)

Then, for a linear twist where the blade angle θ at nondimensional radius x = r/R was
given by the expression

    θ = θ₀ + θ₁x  ... (3)

and for an arbitrary plan form denoted by the solidity factors

    σ₁ = (b/πR²) ∫₀ᴿ c dr  ... (4)
    σ₂ = (b/πR²) ∫₀ᴿ cr dr  ... (5)
    σ₃ = (b/πR²) ∫₀ᴿ cr² dr  ... (6)

and so forth, the equation used to calculate λ₁ from the test values of the thrust
coefficient, blade angle, and rate of descent reduced to

    λ₁ = [σ₃(aθ₀σ₃ + aθ₁σ₄ − ΔC_Q/σ₃)] / (σ₃√(σ₃² σ_T))  +  λ₂  ... (7)

where λ₂ is the nondimensional velocity of descent:

    λ₂ = V / √(T/(2ρA))  ... (8)

    λ₁ = V_i / √(T/(2ρA))  ... (9)

and a is the two-dimensional lift-curve slope corrected for the Reynolds number and Mach
number at the three-quarter-radius point.

Similarly, upon writing the variation of the profile drag coefficient ΔC_d_e for the
symmetrical airfoil as

    ΔC_d_e = D₀ + D₁α² + D₂α⁴  ... (10)

where α_r is the blade-element angle of attack, the solution of the torque equation for λ₁
from torque gives:

    λ₁ = [large expression — see equation (11) in source image]  ... (11)

where ΔC_Q is the value of C_Q for the test point minus the value at zero thrust and zero
rate of descent (i.e., minus the value due to minimum profile drag coefficient).

The Reynolds numbers and Mach numbers at the three-quarter-radius points and the
corresponding estimated values of the lift-curve slope are given for each run in tables I
to VIII.  A value of δ₂ = 1.25 was used to reduce the data, as explained in the section
"Analysis and Discussion."

The values of 1/f and 1/F for the comparison were obtained from the conversion formulas

    1/f = λ₂²  ... (12)

    1/F = (λ₁ − λ₂)²  ... (13)

It is to be noted that the radical of equation (11) may go imaginary if, through
experimental errors, the measured torque coefficient ΔC_Q is too large for the measured
extended blade root angle θ₀.  This was the case for those test points listed in the tables
where the value of ΔC_Q is given but the value of λ₁ (torque) is missing.

---

## Results

The results of the force tests on each rotor are presented in the form of graphs of
θ₀.₇₅R and ΔC_Q against V/ΩR for constant values of C_T/σ_e and as graphs of the
equivalent values of λ₁ against λ₂.  The experimental values for the individual test
points are given in tables I to VIII.

Figures 4, 5, 6, and 7 show the values of θ₀.₇₅R against V/ΩR at constant C_T/σ_e for the
6-foot-diameter rotors having constant-chord, untwisted blades; tapered blades; and twisted
blades; and for the 4-foot-diameter rotor with constant-chord, untwisted blades,
respectively.  Figures 8, 9, 10, and 11, respectively, show the variation of ΔC_Q with
V/ΩR at constant C_T/σ_e for the four rotors.  Figures 12, 13, 14, and 15 show the
variation of λ₁ with λ₂ as calculated from the previous test points.

Figures 16 and 17 show the comparison on λ₁ against λ₂ and 1/f against 1/F coordinates of
the experimental values obtained from the data on the 6-foot-diameter rotor having
constant-chord, untwisted blades with the values from Glauert's empirical curve of 1/f
against 1/F from reference 1; the full-scale values of references 2 and 3; and the values
given by the simple momentum theory.

Sketches of the flow patterns deduced from the photographs of the tufts and smoke
streamers are shown for values of the nondimensional descent velocity λ₂ of 0, 0.3, 1.0,
1.35, 1.7, and 2.0, respectively, in figures 18 to 23.  Figure 24 shows one of the original
smoke photographs taken at λ₂ = 0.3.

The comparison of the theoretical and experimental hovering values of C_T against θ₀.₇₅R
with the values from the end points of the vertical-descent tests are given in figures 25
to 27.  Figures 28 to 30 show the similar comparison for the values of ΔC_Q against C_T.

---

## Analysis and Discussion

### Simple momentum considerations

Consider the case of an actuator disk of area A exerting a thrust T and descending at a
velocity V in a perfect fluid.  For this hypothetical case there is no apparent reason why
a normal wake should not exist with a flow pattern of the type shown in figure 31.
Consequently, the simple momentum theory could be used to determine the relation between T,
V, and the induced velocity with respect to fixed coordinates V_i.  Upon applying the usual
momentum and energy relations it is found that

    V_i = V/2 + √((V/2)² + T/(2ρA))  ... (14)

and, since

    λ₂ = V / √(T/(2ρA))

and

    λ₁ = V_i / √(T/(2ρA))

it follows that

    λ₁ = λ₂/2 + √((λ₂/2)² + 1)  ... (15)

The values of λ₁ given by the above equation might reasonably be expected to constitute a
lower limit on the values that can be obtained on an actual rotor for those
vertical-descent conditions where the air flow through the plane of rotation is
predominantly in a downward direction.

### Formation of a "vortex ring" type flow pattern

Consider, as a second approximation, the case of the actuator disk in a slightly viscous
fluid.  For the hovering flight condition, the principal effect of the fluid viscosity on
the wake is to cause the entrapment of air along the periphery of the wake.  Consequently,
the diameter of a given section of the stream tube enclosing the wake increases with time
and distance from the rotor plane.  An analysis of the similar phenomenon associated with
the expansion of a free jet is given in reference 5.

From the standpoint of elementary vortex theory, the actuator-disk wake can be considered
to be composed of a close succession of vortex rings of very small strength.  The effect of
fluid viscosity upon a vortex filament of one of these rings is to cause a continual increase
in core diameter with time.  Consequently, after a certain increase of time, the strength of
the circulation of a filament measured at any given radius from the axis of the vortex will
decrease with time, as explained in reference 6.  As the impulse of each succession of rings
in the wake tends to remain constant, this implies an increase in ring diameter with time or
distance from the rotor and a decrease in velocity of the corresponding point of the wake
and the corresponding velocity of progression of the rings.  Thus, if the actuator disk is
slowly allowed to descend from the hovering condition, the downward velocity of the axes of
the wake vortex rings will, at some distance below the rotor plane, be less than the descent
velocity of the disk.  When the folded-back sheet has passed above the rotor plane the
induced velocities are in such a direction as to cause it to contract and roll up into the
"vortex ring" type flow pattern observed at small rates of descent.

For steady-state descent, the strength of the "bound vortex ring" formed by the rolling up
of the wake vortex sheet cannot increase with time.  Therefore, the vorticity continually
shed from the rotor and entering the "bound vortex ring" at the same rate as it enters.
Turbulent air exchange between the flow in the ring and the surrounding free-stream flow
appears to be the balancing factor.

At the small rates of descent the scale of the turbulence that is, the volume of the
individual masses of air torn from the "bound vortex ring," appears to be small.  As the
steady-state rate of descent is increased, the scale of the individual vortex masses
increases until, at the higher rates of power-on descent, the turbulence becomes severe
enough to cause fluctuations in the rotor forces and moments.

### Determination of equivalent free-stream velocity

Upon investigation it appeared that the method of testing should duplicate the free-stream
flow patterns in the vicinity of the rotor with sufficient accuracy but that the measured
wind-tunnel velocity would be a poor indication of the equivalent free-stream velocity.
For example, at the hovering and points where the model rotor wake was directed back into
the entrance cone of the tunnel, the net tunnel flow corresponding to the free-flight
hovering condition would obviously be some small flow in a reverse direction and to the
zero value of the free-stream rate of descent.  For the hovering condition the correct
tunnel flow would appear to be that which would correspond to obtaining the static-pressure
field about the model rotor if the energy ratio of each stream tube passing through the
rotor and traversing the circuit of the tunnel were unity, that is, if it were the condition
of a rotor in free flight.  The magnitude of this flow is not very amenable to calculation.
However, it may be noted that this is the quantity of air that would be driven through the
wind-tunnel circuit by the free-stream static-pressure field about the model rotor if the
energy ratio of the stream tubes leaving the rotor were unity.

An open-jet, closed-return wind tunnel with an energy ratio of unity is equivalent
aerodynamically to a frictionless open-return tunnel of the type shown in figure 33.
Consequently, by analogy, the open-jet, closed-return tunnel having an energy ratio of
unity, the total head in a stream tube entering the jet is equal to the velocity head at
the plane of the wind-tunnel fan regardless of the changes in tunnel velocity or tunnel
velocity distribution occurring in the jet.  This is true for the hypothetical open-return
tunnel regardless of the changes in tunnel velocity or tunnel velocity distribution occurring
in the jet because there is no jet-type wake.  However, the net tunnel flow is in a reverse
direction for the stream tubes originating in the wake of a hovering or slowly descending
rotor, and consequently the tunnel energy ratio is higher than that for the tunnel as a
whole.  However, the net tunnel flow is in a reverse direction for the stream tubes
originating in the wake and, consequently, the tunnel energy ratio is higher than that
for stream tubes at the higher rates of descent where the flow is in the normal direction.
This tends to compensate for the diffusion of the wake in the tunnel.

As a compromise, the tunnel energy ratio was reduced for the present tests to a value of
approximately 0.7 by the installation of an 18- by 18-mesh screen, as previously noted.

It was impractical to measure directly the pressure rise through the plane of the
wind-tunnel propeller on account of the very small pressure differences involved.
Therefore, a calibration was obtained of the wind-tunnel jet velocity against wind-tunnel
propeller speed for zero model rotor thrust.  Then, making the approximation that the
pressure rise through the plane of the wind-tunnel propeller was unchanged by any change
in wind-tunnel inflow velocity due to model rotor thrust, the equivalent free-stream
velocity at the measured wind-tunnel propeller speed for a given test point could be
obtained from this calibration curve.

On the present tests with the very low tunnel energy ratio the wind-tunnel propeller-blade-
element inflow angles were small enough that the approximation that the pressure rise through
the plane of the wind-tunnel propeller was independent of the model rotor thrust did not
introduce any large errors.

As a check on the hovering data obtained from the end points of the vertical-descent tests,
additional hovering runs were made with the tunnel blocked at the wake-diffusion screen.
This screen was almost at the circuit of the tunnel "below" the rotor.  Thus the virtual
ground plane was at some distance greater than the edge of the tunnel exit cone, five-sixths
of a rotor diameter, "below" the rotor, and the ground effect was very small, probably
measurable in view of the too perfect agreement with the simple independence of
blade-element theory.

The agreement between the values of the hovering blade angles and torque coefficients
obtained from the end points of the vertical-descent tests with those obtained from the
hovering runs with the tunnel blocked was very satisfactory in the low velocity range.

### Discussion of discrepancies between present data and those of reference 7

The test data of reference 7, used to determine the "vortex ring" portion of Glauert's
curve of 1/f against 1/F, were obtained on a 3-bladed, solid-brass-model of a 3½-inch
chord in a 7- by 7-foot square, closed-jet, indraft tunnel.  The model blade angles were
adjustable, but not controllable, and runs were obtained at various combinations of descent
velocity and zero lift.  The reference blade angle for zero lift was obtained from a
static-pressure wall tap located 8 feet ahead of the plane of the test rotor.

The values of 1/F calculated from the test data in reference 7 are all too high in the
vicinity of hovering because of the use of the statically set blade angle without any
correction for dynamic twist.  For example, the dynamic twist at the 1500 revolutions per
minute on the model of reference 7 is of the order of 1½ percent of the set root blade
angle.  If this were to be taken into account, the value of 1/F would be reduced from 2 to
approximately 1.5.  The remaining difference between the residual value of 1.5 and
full-scale flight-test values of the order of 1.1 corresponds to a difference in blade
angle of only 1/2° which may have been partially attributable to an increase in induced
velocity due to the proximity of the closed jet walls.

The tunnel velocity as measured by the static-pressure wall tap was used in reference 7 as
the descent or ascent velocity at the rotor.  As previously noted, the tunnel velocity
contains an increment of the induced velocity of the rotor, it is higher than the equivalent
free-stream velocity in the vertical-ascent range and lower than the equivalent free-stream
velocity in the vertical-descent range.  Thus, it could be reasoned that the values of 1/f
given by Glauert's curve are too high in the vertical-ascent range and too low in the
correct values of 1/F.  The loop in the original data of reference 7 at the hovering point
would appear to be due to the change in sign of the tunnel velocity in order to obtain the
equivalent free-stream velocity.

### Calculation of full-scale blade angle and torque coefficient for given thrust coefficient and rate of descent from experimental curves of λ₁ against λ₂

The customary assumption of the independence of blade elements in calculating the thrust
and torque of a helicopter rotor in the vertical-descent regime from the experimentally
derived values of 1/f against 1/F or λ₁ against λ₂ appears to be of doubtful validity to
the relations were necessarily calculated from the experimental data on the basis of an
assumed uniform normal component of velocity over the rotor disk, and it would seem that
the same assumption should be used for inverse computations.  Second, that part of the
induced flow due to the vortex distribution in the wake will be considerably changed by
the large-scale turbulent mixing of the wake air at the higher rates of power-on descent.
In other words, that part of the induced flow due to the vortex filaments shed from a
blade at a given radius probably do not remain at the proportional wake radius long enough
for the approximation of the independence of blade elements to be applicable.

Thus, making the same assumptions and approximations that were used to calculate the values
of λ₁ against λ₂ from the experimental data, namely, that the induced velocity is uniform,
the blades are everywhere unstalled, the inflow angle φ can be considered a small angle,
and the tip loss can be neglected, it follows for blades of given plan form denoted by the
solidity factors

    σ₁ = (b/πR²) ∫₀ᴿ c dr
    σ₂ = (b/πR³) ∫₀ᴿ cr dr

and so forth and having a linear twist where the blade angle θ at nondimensional radius x
is given by the expression θ = θ₀ + θ₁x, and the torque coefficient is expressed by the
equation

    C_Q = (b/2πR³) × [long expression involving σ factors and λ₁, λ₂]  ... (16)

The value of a to be used in the above equations is that for the approximate blade angle
of attack at the three-quarter-radius point from the appropriate interpolation of the
lift-curve slope of the airfoil at the Mach number and Reynolds number.  Also, the
coefficients D₀, D₁, and D₂ in the equation for the profile drag c_d_e are determined for
the blade airfoil from the aerodynamic drag polar.

At the higher rates of power-on descent, a certain reduction in the value of a obtained
from a given rotor would appear from the data to be necessary.  However, the full-scale
application of the model test results would appear to require the uncorrected model test
values of a, as explained in the discussion.

---

## Concluding Remarks

The mean nondimensional induced velocities calculated from the present tests are
considerably less for hovering and very small rates of descent and considerably larger for
the higher rates of descent than those given by Glauert's curve of 1/f against 1/F (where
F is the thrust coefficient based on the resultant velocity at the rotor).  The major
portion of the disagreement can be accounted for and is due to the fact that previously no
correction was made for dynamic blade twist and the measured tunnel velocity was taken as
the free-stream velocity.

The present data are in good agreement with full-scale flight-test results at the hovering
and autorotation ends of the range, but the peak values of the nondimensional induced
velocity at the large rates of power-on descent are higher than those obtained from the
full-scale flight tests reported by Stewart of references 2, as shown in figure 16.  This
discrepancy at the large rates of power-on descent may have been largely due to the
inability, through loss of control, to maintain the desired flight condition in the
equilibrium flow to be established.

The primary effects of the 3/1 blade taper were to decrease slightly the mean induced
velocity at the small rates of descent and to increase the rate of descent for autorotation
by approximately 3 percent over that for the rotor with constant-chord, untwisted blades
operating at the same thrust coefficient.

Linear twist of 12° increased the "ideal" nondimensional rate of descent for autorotation
by about 10 percent compared with the value for the rotor with the constant-chord, untwisted
blades.  The peak value of the mean nondimensional induced velocity was increased by
approximately 24 percent and it occurred at a nondimensional rate of descent that was about
17 percent higher than for the rotor with constant-chord, untwisted blades.  Also, the
fluctuations in the forces and moments on the rotor with the twisted blades were very much
larger at the higher rates of power-on descent than for the rotors with the tapered or
constant-chord, untwisted blades.  As in the case of the tapered blades, the mean induced
velocity of the rotor with the twisted blades was slightly less, at hovering, than that for
the constant-chord, untwisted blades.

There were no observable fluctuations in forces or moments on any of the rotors in the
autorotation range.

Within the range and accuracy of these tests there were no significant differences in the
curves of nondimensional induced velocity λ₁ against the nondimensional descent velocity
λ₂ due to variations in the thrust coefficient, rotor speed, or rotor diameter.

The present data should be more applicable to full-scale, free-flight calculations than
data from previous model rotor, vertical-descent tests on account of the inclusion of a
correction for the dynamic blade twist and the more exact method used to determine the
equivalent free-stream descent velocity.

*Georgia Institute of Technology*
*Atlanta, Ga., May 31, 1950*

---

## Appendix

### Integrated Thrust Equations for Hovering Rotors with Linearly Tapered and/or Twisted Blades

Assuming independence of blade elements and neglecting the rotation of the slipstream, it
follows from the momentum theory that the thrust dT at radius r can be expressed as

    dT = 4πρV_i² r dr  ... (A1)

Also, upon making the inflow angle φ is a small angle, blade-element analysis gives

    dT = ½ ρb c σ₁² r² dr  ... (A2)

Thus, for rotors having blades with a linear taper where the chord c at nondimensional
radius x can be denoted in terms of the extended blade root chord c₀ and the taper factor t
by the expression

    c = c₀(1 + tx)  ... (A3)

the inflow angle φ is

    φ = −V_i/(Ωr) + √(σ₀c₀(1 + tx)/(8x))  ... (A4)

where σ₀ = bc₀/(πR), the solidity of the extended blade root chord.  However, from the
geometry

    φ = −θ + c₁/a  ... (A5)

or, for rotors having blades with a linear twist where the blade angle θ at nondimensional
radius x can be expressed in terms of the extended blade root pitch angle θ₀ and the twist
θ₁ as

    θ = θ₀ + θ₁x  ... (A6)

it follows that

    c₁/a = θ₀ + θ₁x + V_i/(Ωr) − √(σ₀c₀(1 + tx)/(8x))  ... (A7)

This expression can, for convenience, be factored, giving

    16c₁/(aσ₀²) = θ₀ + θ₁x + [expression in V_i terms] × [radical term]  ... (A8)

where

    θ₀ = 16θ₀/(aσ₀)  ... (A9)
    θ₁ = 16θ₁/(aσ₀)  ... (A10)

Setting up the expression for the thrust coefficient where

    C_T = T/(ρπR²(ΩR)²) = (b/(πR²)) ∫₀ᴿ (½ c c₁(ΩR)²x²) dx  ... (A11)

and substituting the previous value of c₁ given by equation (A7) yields

    32C_T/(a²σ₀²) = ∫₀¹ [θ₀ + θ₁x + (1 + tx)^(1/2) × terms]² × (1 + tx)x² dx  ... (A12)

Integrating the first three terms of equation (A12) and expanding the factor (1 + tx)^(3/2)
in the fourth term by means of the binomial theorem give

    32C_T/(a²σ₀²) = θ₀·I + θ₁·II + III + IV + V + ...  ... (A13)

where

    I   = ∫₀¹ (1 + tx)x² dx = A  ... (A14) — see figure
    II  = ∫₀¹ (1 + tx)x³ dx = B  ... (A15) — see figure
    III = ∫₀¹ (1 + tx)^(3/2) x² dx — see eqs (A16–A18)
    IV  = ...
    V   = ...
    A   = [t(2θ₀ + t) + 2(θ₀ + t) + 1]^(3/2)  ... (A20)

and, for the case of interest where θ₁ is negative,

    B = (1/8θ₁) × [(2θ₀ + 4θ₁ + t)(2θ₀ + 2θ₁ + t + 1)^(1/2) − 2θ₀ − t]
      + [(2θ₀ + t)² − 8θ₁] / [16θ₁√(−2θ₁)]
      × { sin⁻¹[(2θ₀ + 4θ₁ + t)/((2θ₀+t)² − 8θ₁)^(1/2)]
        − sin⁻¹[(2θ₀ + t)/((2θ₀+t)² − 8θ₁)^(1/2)] }  ... (A21)

It is to be noted that the angles in the above equation are in the first or fourth quadrant
depending upon whether the arc sine is positive or negative.

The maximum error introduced in the value of C_T because the sixth and higher terms of the
binomial expansion were dropped is less than 1/2 percent for the extreme case where θ₁ = −0.2
radian and t = −2/3.

For θ₁ = 0 the latter terms of equation (A14) become imaginary.  Thus for tapered but
untwisted blades it is necessary to set θ₁ = 0 before integrating.  Then for tapered but
untwisted blades

    32C_T/(a²σ₀²) = θ₀(1/3 + t/4) + 1/2 + 2/3·t + 1/4·t² + I₁ + II₁  ... (A22)

where

    I₁  = −1/(3t(2θ₀+t)) × (A₁) + 1/(3t(2θ₀+t)) × (I₁) + 1/(t(2θ₀+t)) × (B₁)  ... (A23)
    II₁ = −1/(4(2θ₀+t)) × (A₁) − 5(θ₀+t)/(4(2θ₀+t)) × (I₁) + 1/(4(2θ₀+t)) × (B₁)  ... (A24)
    A₁  = [t(2θ₀+t) + 2(θ₀+t) + 1]^(3/2)  ... (A25)
    B₁  = (1/(2t(2θ₀+t))) × [t(2θ₀+t) + θ₀+t][t(2θ₀+t) + 2(θ₀+t) + 1]^(1/2)
          − (θ₀²)/(2t(2θ₀+t)√(−t(2θ₀+t)))
          × [sin⁻¹((t(2θ₀+t) + θ₀+t)/θ₀) − sin⁻¹(t/(2θ₀+t)/θ₀)]  ... (A26)

For t = 0, the latter terms of equation (A22) become imaginary.  Thus, for constant-chord
but untwisted blades it is necessary to set t = 0 before integrating.  Then for
constant-chord, untwisted blades

    32C_T/(a²σ₀²) = θ₀/3 + (1 − 3θ₀)(1 + 2θ₀)^(3/2) − 1 / (15θ₀²)  ... (A27)

The normal procedure, in using the previous equations which eliminate most of the labor
involved in the customary trial-and-error solution for the radial distribution of c₁, is as
follows for the usual case where it is desired to take tip loss into account:

1. Calculate R_e, the effective radius, where R_e = R − ½ c_tip.
2. Calculate σ₀ = bc₀/(πR_e), C_T' = T/(πρ²R_e⁴), and the factor 32C_T'/(a²σ₀²).
3. Calculate the values of θ₁ = 16θ₁/(aσ₀) and t = c_tip/c₀ − 1.
4. Calculate the values of θ₀ = 16θ₀/(aσ₀) for several assumed values of θ₀ and
   determine the value of θ₀ which will yield the desired thrust coefficient.
5. Calculate the values of θ₀, σ₁, σ₂ from equations (A7), (A11), and the airfoil drag
   polar, giving the desired value of thrust coefficient.
6. Calculate the radial distribution of c_d_e, based on the effective radius, by graphical
   integration where

        C_Q' = ∫₀¹ (½ c c_d_e x³) dx + ∫₀¹ (½ c c_l·tan(α)) x³ dx

   and E' = R/R_e.  The value of c_d_e existing at x = 1 can be assumed to extend to R'.

7. Calculate the value of the torque coefficient

        C_Q = Q/(ρπR⁴(ΩR)²)

If the radial distributions of the blade air loads are desired, they can be calculated in
the usual manner from the results of items 1 through 5 above.

---

## References

1. Glauert, H.: The Analysis of Experimental Results in the Windmill Brake and Vortex Ring
   States of an Airscrew.  R. & M. No. 1026, British A.R.C., 1926.

2. Stewart, W.: Flight Testing of Helicopters.  Jour. R.A.S., vol. 52, no. 449, May 1948,
   pp. 261–292; discussion, pp. 493–101.

3. Gessow, Alfred: Flight Investigation of Effects of Rotor-Blade Twist on Helicopter
   Performance in the High-Speed and Vertical-Autorotative-Descent Conditions.  NACA TN
   1666, 1948.

4. Gustafson, F. B., and Gessow, Alfred: Flight Tests of the Sikorsky HNS-1 (Army YR-4B)
   Helicopter.  II — Hovering and Vertical-Flight Performance with the Original and an
   Alternate Set of Main-Rotor Blades, Including a Comparison with Hovering Performance
   Theory.  NACA MR L5D09a, 1945.

5. Prandtl, L.: The Mechanics of Viscous Fluids.  Spread of Turbulence.  Vol. III of
   Aerodynamic Theory, div. G, sec. 25, W. F. Durand, ed., Julius Springer (Berlin), 1935,
   pp. 172–173.

6. Prandtl, L.: The Mechanics of Viscous Fluids.  Some Examples of Exact Solutions.  Vol.
   III of Aerodynamic Theory, div. G, sec. 10, W. F. Durand, ed., Julius Springer (Berlin),
   1935, pp. 68–69.

7. Lock, C. N. H., Bateman, H., and Townsend, H. C. H.: An Extension of the Vortex Theory
   with Applications to Airscrews of Small Pitch, Including Experimental Results.  R. & M.
   No. 1014, British A.R.C., 1926.
