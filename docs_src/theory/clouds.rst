Clouds
======

Bärli treats clouds as a list of independent **cloud modules**. Every module provides, on the
spectral grid and per layer, the extinction optical depth :math:`\tau`, the single-scattering
albedo :math:`\omega` and the Henyey-Greenstein asymmetry parameter :math:`g`. The model sums
the modules, scales the result with the cloud ramp and hands it to the radiative transfer
together with the gas opacities. Three module types exist: a prescribed table (``fixed``), a
prescribed deck of Mie particles (``mie_slab``) and the Ackerman & Marley (2001) eddy
sedimentation cloud (``ackerman_marley``). The parameter reference is in
:doc:`../user_guide/clouds`.

Combining modules and gas
-------------------------

For every spectral point and layer the enabled modules :math:`k` are combined with the only rule
that is consistent for a mix of absorbing and scattering decks,

.. math::

   \tau = \sum_k \tau_k, \qquad
   \omega = \frac{\sum_k \tau_k\, \omega_k}{\tau}, \qquad
   g = \frac{\sum_k \tau_k\, \omega_k\, g_k}{\sum_k \tau_k\, \omega_k},

i.e. the albedo is extinction-weighted and the asymmetry scattering-weighted. In the radiative
transfer the layer optical depth of gas and cloud are added, the gas Rayleigh scattering being
isotropic:

.. math::

   \tau_j = \tau^{\mathrm{gas}}_{\mathrm{abs}} + \tau^{\mathrm{gas}}_{\mathrm{sca}}
   + \tau_c (1 - \omega_c) + \tau_c\, \omega_c, \qquad
   \omega_j = \frac{\tau^{\mathrm{gas}}_{\mathrm{sca}} + \tau_c \omega_c}{\tau_j}, \qquad
   g_j = \frac{\tau_c\, \omega_c}{\tau^{\mathrm{gas}}_{\mathrm{sca}} + \tau_c\, \omega_c}\; g_c ,

with :math:`\omega_j` clipped at 0.999999. The cloud absorption also enters the absorption
coefficient that weights the local radiative-equilibrium residual of the temperature solvers,
:math:`\kappa^{\mathrm{tot}} = \kappa^{\mathrm{gas}} + \kappa^{\mathrm{cloud}}`: with a cloud,
local equilibrium reads :math:`(\kappa_{\mathrm{gas}} + \kappa_{\mathrm{cloud}})(J - B) = 0`, and
weighting with the gas alone drives the Newton to a wrong root.

Operator splitting
------------------

Cloud modules are updated **once per outer iteration** from the committed atmosphere
(temperature, structure, composition, :math:`K_{zz}`). The temperature corrector's trial
evaluations and Jacobians see the cloud as frozen opacity; the radiative Jacobian never contains
:math:`\partial\tau_c/\partial T`. Convergence of a cloudy model therefore relies on damping the
feedback between the cloud structure and the temperature, which is the job of the relaxation
knobs of the microphysical module and of the model-level cloud ramp.

The cloud ramp
--------------

A thick cloud deck switched on in a deep atmosphere that is still far from equilibrium is
pathological for the Newton (measured: a :math:`\tau = 5` deck at 25 to 60 bar put the bottom
level 1000 K above its neighbour during the mixing-length ramp). The committed cloud optical
depths therefore enter the radiative transfer multiplied by a scale factor:

1. the scale is zero until the corrector's easy start is over;
2. at onset, the scale is the largest value that keeps the total optical depth at 1 micron below
   ``tau_start``;
3. on every outer iteration whose residual is below ``gate`` the scale grows by ``factor``
   (default 2), so the Newton re-converges between doublings and the greenhouse of a thick deck
   builds up gradually;
4. a step that lifts the residual above ``backoff`` has left the Newton's contraction region: the
   step is undone and the growth factor reduced to its square root (measured: a doubling into a
   thick forsterite plus iron deck took the residual from :math:`7 \times 10^{-3}` to 9 and froze
   the run at minimum damping for 90 iterations; a growth factor of 1.2 converged).

A warm start applies the full cloud at once. ``cloud_ramp = "none"`` disables the ramp.

Prescribed clouds
-----------------

``fixed``
~~~~~~~~~

A table of :math:`(\tau, \omega, g)` over pressure levels and wavelengths, read from a file or
handed over at run time, the hook for an external cloud or microphysics code. The table's row at
level :math:`p_k` holds the optical depth of the input layer between :math:`p_k` and the next
level towards lower pressure. The table is mapped **conservatively** onto the model layers: the
cumulative optical depth :math:`C(x) = \sum_{k' < k} \tau_{k'}` at the input levels is
interpolated linearly in :math:`x = -\ln p` onto the model levels and differenced,

.. math::

   \tau^{\mathrm{model}}_j = C(x_{j+1}) - C(x_j),

so the column optical depth is independent of the model grid. :math:`\omega` and :math:`g` are
interpolated in :math:`\ln p` at the layer midpoints and linearly in wavelength, held constant
outside the table's wavelength range; a single-wavelength table is grey.

``mie_slab``
~~~~~~~~~~~~

A deck of Mie particles of one condensate between two pressures, with the column number density
:math:`N_{\mathrm{col}}` distributed uniformly in :math:`\ln p` and a lognormal size distribution
with geometric mean radius :math:`r_g` and geometric standard deviation :math:`\sigma`
(:math:`\sigma = 1` is monodisperse). The optics per unit column follow from the Mie table
(below); with ``optical_depth`` given, the column is
:math:`N_{\mathrm{col}} = \tau_{\mathrm{ref}} / \sigma_{\mathrm{ext}}(\lambda_{\mathrm{ref}})`,
so the deck has exactly the requested optical depth at the spectral point closest to the
reference wavelength. Nothing depends on the atmospheric state but the pressure grid, so the
optics are computed once.

Mie optics
----------

The Mie efficiencies of every condensate are computed with `LX-MIE
<https://github.com/NewStrangeWorlds/LX-MIE>`_ (Kitzmann & Heng 2018) from the complex
refractive index :math:`n + ik` of the compilation shipped in ``data/refractive_indices/``
(Kitzmann & Heng 2018; the source of each data set is given in the file header). Because the
efficiencies of a size *distribution* are smooth in wavelength, they are tabulated on a coarse
grid of resolving power :math:`R = 200` and 120 log-spaced radii between :math:`10^{-4}` and
542 micron (virga's radius grid, with twice its number of radii), with the single-size resonance
ripple averaged over six sub-bins per radius bin, and then interpolated onto the opacity-sampling
grid: linearly in :math:`\ln\lambda`, and :math:`\ln Q` linear in :math:`\ln r`, which is exact
for the Rayleigh power laws and flat in the geometric limit. The cost of this is about 1 % against
a direct Mie call and a factor of about 100 in run time. Outside the data range the refractive
index is held constant.

A population :math:`\{(r_j, N_j)\}` gives the per-layer optics

.. math::

   \tau_{\mathrm{ext}} = \sum_j N_j\, \pi r_j^2\, Q_{\mathrm{ext}}(r_j), \qquad
   \omega = \frac{\sum_j N_j \pi r_j^2 Q_{\mathrm{sca}}}{\tau_{\mathrm{ext}}}, \qquad
   g = \frac{\sum_j N_j \pi r_j^2 Q_{\mathrm{sca}}\, g_j}{\sum_j N_j \pi r_j^2 Q_{\mathrm{sca}}} .

The lognormal distribution is integrated with 41 points equally spaced in :math:`\ln(r/r_g)`
over :math:`\pm 4\ln\sigma`, with trapezoid weights renormalised so that the discrete sum equals
the column exactly.

The Ackerman & Marley cloud
---------------------------

The ``ackerman_marley`` module is the eddy-sedimentation cloud of Ackerman & Marley (2001) for
**one condensate**, in the form of virga's layer solver (Rooney et al. 2022), but fed by the
condensation of the equilibrium chemistry instead of vapour-pressure fits. Several condensates
are several modules whose optics are summed; the particles do not interact.

Condensation source
~~~~~~~~~~~~~~~~~~~

The ``equilibrium`` chemistry with ``condensation = "rainout"`` stores at every level the amount
of each condensate that forms *at that level* from the vapour arriving from below, and depletes
the gas phase above accordingly (:doc:`chemistry`). Above the cloud base the vapour is at
saturation in both the Ackerman & Marley picture and the rainout, so the only quantity the
sedimentation efficiency changes is the suspended condensate. The condensation increment of layer
:math:`j` (between levels :math:`j` and :math:`j+1`) as a mass mixing ratio is

.. math::

   \Delta q_j = \frac{n_{\mathrm{cond}, j+1}\, M}{N_A\, \rho_{j+1}} ,

with the formula weight :math:`M` of the condensate from FastChem's stoichiometry. The module
never touches the composition; the gas-phase depletion is the chemistry's job. Only the rainout
mode has this increment semantics; the ``equilibrium`` condensation mode stores the full local
condensate and would be double counted.

Suspended condensate
~~~~~~~~~~~~~~~~~~~~

Ackerman & Marley balance the upward eddy transport of the total condensable against the
sedimentation of the condensate,

.. math::

   -K_{zz}\, \frac{\partial q_t}{\partial z} - f_{\mathrm{sed}}\, w_\ast\, q_c = 0, \qquad
   w_\ast = \frac{K_{zz}}{L},

whose solution between sources decays as :math:`\exp(-f_{\mathrm{sed}}\, z/L)`. The module uses
the closed-form layer recursion of virga, walking upward from the bottom,

.. math::

   q_{c,j} = q_{c,j-1}\, \exp\!\left(-\frac{f_{\mathrm{sed},j}\, \Delta z_j}{L_j}\right) + \Delta q_j ,

with the mixing length :math:`L_j` equal to the layer-mean scale height (``mixing_length =
"scale_height"``) or the scale height times :math:`\max(0.1, \nabla/\nabla_{\mathrm{ad}})`
(``"lapse_ratio"``). The layer :math:`K_{zz}` is the mean of the bounding levels floored at
``kzz_min``. A tail below :math:`10^{-6}` of the maximum is cut, because it carries no optical
depth but would produce nanometre particles. The sedimentation efficiency is constant
(``fsed_mode = "const"``) or follows the exponential form of Rooney et al. (2022) above the
cloud base,

.. math::

   f_{\mathrm{sed},j} = \varepsilon + (f_{\mathrm{sed}} - \varepsilon)
   \exp\!\left(-\frac{z_j - z_{\mathrm{base}}}{b\, L_j}\right) .

Particle sizes
~~~~~~~~~~~~~~

In every cloudy layer the radius :math:`r_w` whose terminal fall velocity equals the convective
velocity scale :math:`w_\ast` is found by Brent's method, and the local power law
:math:`v_f \propto (r/r_w)^\alpha` is fitted by least squares over the radius grid. With
:math:`\ell = \tfrac{1}{2}\ln^2\sigma` the lognormal moments give (Ackerman & Marley Eqs. 13 to
18)

.. math::

   r_g = r_w\, f_{\mathrm{sed}}^{1/\alpha}\, e^{-(\alpha+6)\ell}, \qquad
   r_{\mathrm{eff}} = r_g\, e^{5\ell}, \qquad
   N_j = \frac{3 M_j}{4\pi \rho_p r_g^3}\, e^{-9\ell}, \qquad
   M_j = \rho_j\, q_{c,j}\, \Delta z_j ,

the geometric mean radius, the effective radius, the column number density and the condensate
mass column of the layer, with the bulk density :math:`\rho_p` of the condensate.

Fall velocity
~~~~~~~~~~~~~

The terminal velocity follows virga's implementation of Ackerman's original. The gas is
characterised by the kinetic-theory mean free path and the Lennard-Jones corrected viscosity
(Ackerman & Marley Eq. B2; Rosner 2000) for a molecule of diameter :math:`d` and well depth
:math:`\varepsilon_k` (defaults for H\ :sub:`2`),

.. math::

   \lambda_{\mathrm{mfp}} = \frac{1}{\sqrt{2}\, n\, \pi d^2}, \qquad
   \eta = \frac{5}{16}\, \frac{\sqrt{\pi k_B T m}}{\pi d^2}\, \frac{1}{1.22\,(T/\varepsilon_k)^{-0.16}} .

Stokes flow with the Cunningham slip correction :math:`\beta = 1 + 1.26\, \mathrm{Kn}` gives

.. math::

   v = \beta\, \frac{2}{9}\, \frac{(\rho_p - \rho)\, g\, r^2}{\eta} ;

if the resulting Reynolds number :math:`\mathrm{Re} = 2 r \rho v/\eta` exceeds 1, the
Pruppacher & Klett (1978) polynomial fit of :math:`\ln \mathrm{Re}` in :math:`\ln(C_d
\mathrm{Re}^2)` is used up to :math:`\mathrm{Re} = 10^3`, and above that the Newton regime with a
constant drag coefficient :math:`C_d = 0.45`.

Damping of the cloud structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cloud is recomputed only while the corrector's residual is below ``gate`` (the first
structure is always taken): a cloud that re-forms on every iteration from a :math:`K_{zz}` that
is itself still ramping keeps the Newton from settling. The committed mass column and radius per
layer are then relaxed in log space towards the new values with weight ``relax``, capped at
``max_change`` dex per update, and committed only if the largest change over the layers that
carry mass exceeds ``tolerance`` dex. A layer that becomes cloudy starts at the fraction ``seed``
of its target; a layer the cloud leaves decays by the cap per update. Both the mass and the
radius are damped because committing the physical radius at once lets :math:`\tau \sim M/r` jump
even when :math:`M` is capped, which was measured to produce a period-two limit cycle between the
deck and the detached convective patches it creates. The optics are recomputed from the committed
structure at every update; the onset of the cloud *opacity* is the cloud ramp's job.

Condensates
~~~~~~~~~~~

The condensate table maps a name to its refractive-index file, bulk density and FastChem
condensate symbol. Condensation comes exclusively from FastChem's thermochemical data; there are
no vapour-pressure curves in Bärli, and the rainout equilibrium decides which species form. In a
solar-composition brown dwarf, for instance, forsterite takes the magnesium before enstatite can
form, so an MgSiO\ :sub:`3` module is legitimately empty and Mg\ :sub:`2`\ SiO\ :sub:`4` is the
silicate to use.

Limitations
-----------

* Homogeneous spherical particles; no mixed grains, no nucleation, no interaction between
  condensates.
* The lognormal distribution is truncated at :math:`\pm 4\ln\sigma`; radii outside the Mie
  table (:math:`10^{-4}` to 542 micron) are clamped.
* The cloud is lagged in the temperature Newton; convergence depends on the damping knobs and
  the ramp.
* Particle sizes scale with :math:`K_{zz}`; a prescribed :math:`K_{zz}` far below the convective
  value under a deck leads to limit cycles with the detached convective patches the deck creates
  (:doc:`../user_guide/convection`).
