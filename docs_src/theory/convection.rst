Convection and vertical mixing
==============================

Convection enters Bärli in two forms: the **mixing-length theory** (``mlt_dry``,
``mlt_moist``), in which the convective flux is a smooth function of the super-adiabaticity
and is solved for together with the radiative flux inside the Newton corrector, and the
classical **convective adjustment** (``dry``, ``moist``), which resets unstable layers to the
adiabat after every temperature step. The eddy diffusion coefficient :math:`K_{zz}` that
parametrises vertical mixing for the chemistry and the clouds is derived from the same
mixing-length flux law.

Stability criterion
-------------------

All schemes test a **link** between the levels :math:`i` (lower) and :math:`i+1` (upper) with
the actual logarithmic temperature gradient

.. math::

   \nabla_i = \frac{\ln T_i - \ln T_{i+1}}{\ln p_i - \ln p_{i+1}}

against the scheme's neutrality gradient averaged over the two levels,

.. math::

   \nabla_{\mathrm{ad},i} = \tfrac{1}{2}\left[\nabla_{\mathrm{ad}}(T_i, p_i, \vec n_i)
   + \nabla_{\mathrm{ad}}(T_{i+1}, p_{i+1}, \vec n_{i+1})\right].

The link is convectively unstable if :math:`\nabla_i > \nabla_{\mathrm{ad},i}`: the
Schwarzschild criterion in :math:`(T, p)` form. The neutrality gradient is the dry adiabatic
gradient for the ``*_dry`` and ``dry`` schemes and the saturated pseudo-adiabat for the
``*_moist`` and ``moist`` schemes (:doc:`atmosphere`). There is no composition term (Ledoux
criterion), no overshoot, and no radiative-loss (efficiency) correction of the convective
elements.

Mixing-length convection
------------------------

The flux law
~~~~~~~~~~~~

In the mixing-length scheme every level remains a free unknown of the temperature solver and
convection enters the energy balance as a flux carried on the links. For the link between the
levels :math:`i` and :math:`i+1`, with the link-averaged density :math:`\rho`, temperature
:math:`T_m`, heat capacity :math:`c_p` and scale height :math:`H_p`, the convective flux is

.. math::
   :label: mlt-flux

   F_{c,i} = C_i\, x_i^{3/2}, \qquad
   x_i = \max\!\left(\nabla_i - \nabla_{\mathrm{ad},i},\, 0\right), \qquad
   C_i = \tfrac{1}{2}\, \rho\, c_p\, T_m\, \frac{\lambda^2}{H_p}\sqrt{\frac{g}{H_p}} .

This is the standard mixing-length result in the limit of efficient convection, written in
the form used by the AGNI model: with the convective velocity and the parcel temperature
excess

.. math::

   w_{\mathrm{conv}} = \lambda \sqrt{\frac{g\, x}{H_p}}, \qquad
   \Delta T = T\, x\, \frac{\lambda}{H_p},

the flux is :math:`F_c = \tfrac{1}{2}\rho\, c_p\, w_{\mathrm{conv}}\, \Delta T`. Up to the
numerical prefactor this is the Böhm-Vitense form with the radiative losses of the elements
neglected.

The **mixing length** is

.. math::

   \lambda = \begin{cases}
     \alpha\, H_p & \text{brown dwarfs and gas planets,}\\[4pt]
     \dfrac{k z}{1 + k z / (\alpha H_p)}, \quad k = 0.4 & \text{terrestrial planets,}
   \end{cases}

where :math:`z` is the altitude of the link above the bottom boundary. The second form is the
Blackadar wall law, :math:`\lambda \to k z` close to the surface and :math:`\lambda \to \alpha
H_p` far from it; it is appropriate for a solid surface but not for the artificial domain cut
of a self-luminous object, which therefore uses :math:`\alpha H_p` throughout. The parameter
``alpha`` defaults to 1. Links above the pressure ceiling ``min_pressure`` carry no convective
flux (:math:`C_i = 0`).

Because :math:`C_i` is large wherever the atmosphere is dense, the law acts as a stiff penalty
that drives :math:`x \to 0^+` in convective regions: the converged temperature profile is
insensitive to ``alpha``, which is exposed mainly for sensitivity tests. The eddy diffusion
coefficient, however, depends on it (see below).

How the flux enters the energy balance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``ratio_ul`` corrector consumes the flux law directly (:doc:`temperature_iteration`):

* the level-collocated convective flux :math:`F_{c}(p_i) = \tfrac{1}{2}(F_{c,i-1} + F_{c,i})`
  (one-sided at the domain ends) is added to the radiative flux in every flux-based residual
  row and in the flux-conservation metric;
* the volumetric convective heating :math:`q_i = -(\dd F_c/\dd z)_i / (4\pi)` enters the
  local radiative-equilibrium rows as a conservative difference of the two adjacent link
  fluxes;
* the Unsöld-Lucy correction anchors the **total** flux :math:`F_{\mathrm{rad}} + F_c` to the
  target flux, so that at convergence
  :math:`F_{\mathrm{rad}}(p) + F_{\mathrm{conv}}(p) = F_\star` at every level.

The thermodynamic prefactors :math:`\rho`, :math:`c_p`, :math:`T_m`, :math:`H_p`,
:math:`\lambda` and :math:`\nabla_{\mathrm{ad}}` are frozen at the committed profile within one
correction (the same operator split as the frozen opacity); only :math:`\nabla(T)` is live, so
the analytic link derivatives

.. math::

   \frac{\partial F_{c,i}}{\partial T_i} = \frac{3}{2}\, \frac{C_i \sqrt{x_i}}{T_i\,\Delta\ln p_i},
   \qquad
   \frac{\partial F_{c,i}}{\partial T_{i+1}} = -\frac{3}{2}\, \frac{C_i \sqrt{x_i}}{T_{i+1}\,\Delta\ln p_i}

are the exact residual gradient. Both vanish as :math:`x \to 0^+`: the flux law is
:math:`C^1` across the radiative-convective boundary, which therefore never has to be located
explicitly. No convective mask, no adiabatic slaving and no boundary controller are needed.

The ``convective`` flags are set diagnostically after each Newton step: every level bounding a
link with :math:`F_c > 0` is flagged, and the iteration log prints their number as ``N_conv``.

The easy-start ramp
~~~~~~~~~~~~~~~~~~~

Started from a profile with large super-adiabaticities, the stiff flux law would throw the
Newton far outside its contraction region. The corrector therefore scales the law with a
homotopy factor :math:`s_f` that starts at :math:`3 \times 10^{-4}` and doubles once the
previous stage has settled (residual and step below a stage tolerance of 0.3), until
:math:`s_f = 1`. While ramping, the reported residual is floored so that convergence cannot be
declared on a weakened convection; the ramp costs about fifteen iterations.

The ramp is skipped (:math:`s_f = 1` from the start) when

* the profile is already on the stable side everywhere, i.e. the largest link
  super-adiabaticity below the pressure ceiling is below 0.02, **and** the object is
  surface-driven (terrestrial) or irradiation-dominated
  (:math:`\Fint < 10^{-2}(\mu_\star S + \Fint)`); the ``adiabat`` initial profile with its 2 %
  stable bias is built for this;
* or :meth:`compute` is called with ``warm_start=True``.

A self-luminous object with :math:`\Fint > 0` always ramps on a cold start (measured: the brown
dwarf test case did not converge in 100 iterations without the ramp and needed 18 with it).

Moist mixing-length convection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mlt_moist`` is the same flux law with the saturated pseudo-adiabat as neutrality gradient.
Because that gradient is evaluated for a saturated parcel regardless of the actual humidity,
the pressure ceiling ``min_pressure`` must exclude the stratosphere (:math:`10^{-2}` bar for
the Earth example); otherwise the whole stratosphere is read as moist-unstable and flattens
onto the pseudo-adiabat.

Convective adjustment
---------------------

The ``dry`` and ``moist`` schemes perform the classical hard adjustment after each temperature
step of the relaxation and flux-divergence solvers. Up to ``max_sweeps`` sweeps walk the links
from the bottom; each sweep stops at the first level above ``min_pressure``. For an unstable
link the pair is reset to the neutrality gradient, :math:`T_{i+1} = r\, T_i` with
:math:`r = (p_{i+1}/p_i)^{\nabla_{\mathrm{ad},i}}`, conserving the enthalpy of the pair with
the mean heat capacities :math:`c_{p,i}`, :math:`c_{p,i+1}`:

.. math::

   T_i^{\mathrm{new}} = \frac{c_{p,i} T_i + c_{p,i+1} T_{i+1}}{c_{p,i} + c_{p,i+1}\, r},
   \qquad T_{i+1}^{\mathrm{new}} = r\, T_i^{\mathrm{new}} .

The bottom level is the anchor set by the flux balance and is never changed; the bottom link
only adjusts the level above it. The pair enthalpy is conserved without pressure weighting of
the two levels. A sweep without any adjustment ends the loop early. Both levels of an adjusted
link are flagged convective; the top of a contiguous run of flagged levels is the
radiative-convective boundary of these schemes.

With the Newton-type correctors the ``dry``/``moist`` schemes are not applied as a post-step
adjustment: the correctors handle convection internally by slaving the convective levels to
the scheme's neutrality gradient (:doc:`temperature_iteration`). The explicit adjustment is
used by the relaxation solvers and the ``helios`` scheme. With ``convection = "none"`` nothing
is adjusted and the corrector runs purely radiatively.

Initial adiabat
---------------

The ``adiabat`` initial profile integrates upward from the surface temperature along the
**active scheme's** neutrality gradient, reduced by 2 %, and floors at the stratosphere
temperature:

.. math::

   T_i = \max\!\left(T_{\mathrm{strat}},\;
   T_{i-1}\,\left(\frac{p_i}{p_{i-1}}\right)^{0.98\, \nabla_{\mathrm{ad}}(T_{i-1}, p_{i-1}, \vec n_{i-1})}\right).

The chemistry is evaluated before (isothermal) and after the integration, because the gradient
depends on the composition and temperature-dependent chemistry such as the Manabe-Wetherald
humidity depends on the result. The 2 % bias puts every link on the stable side of the scheme's
own threshold, which is what allows the mixing-length corrector to skip its easy-start ramp.

.. _theory-kzz:

Vertical mixing: the eddy diffusion coefficient
-----------------------------------------------

The eddy diffusion coefficient :math:`K_{zz}(p)` is an atmosphere property consumed by the
quench chemistry, the Ackerman & Marley clouds and the kinetics coupling. It is selected by the
``kzz`` module spec (:doc:`../user_guide/convection`).

Prescribed profiles
~~~~~~~~~~~~~~~~~~~

* ``constant``: :math:`K_{zz} = \mathrm{value}` at every level.
* ``power_law``: :math:`K_{zz}(p) = \mathrm{value}\,(p/p_{\mathrm{ref}})^{\mathrm{slope}}`
  with :math:`p_{\mathrm{ref}} =` ``pressure`` in bar. A negative slope (mixing growing
  upward) is the GCM-fitted form for the radiative envelope of irradiated planets (e.g.
  Parmentier et al. 2013 for HD 209458b).

Self-consistent mixing-length profile (``mlt``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In convective links the coefficient follows from the same flux law as the convection, with the
unramped super-adiabaticity :math:`x` and the convective velocity
:math:`w_{\mathrm{conv}} = \lambda\sqrt{g x / H_p}`:

* ``scaling = "velocity"`` (default, AGNI's choice):

  .. math::

     K_{zz} = \lambda\, w_{\mathrm{conv}} = \lambda^2 \sqrt{\frac{g\, x}{H_p}} ;

* ``scaling = "flux"`` (Ackerman & Marley 2001; Charnay et al. 2015; the recipe of picaso),
  with :math:`F_{\mathrm{conv}}` from the flux law :eq:`mlt-flux` and the specific gas constant
  :math:`R_{\mathrm{spec}} = k_B/(\bar\mu\, m_u)`:

  .. math::

     K_{zz} = \frac{H_p}{3}\left(\frac{\lambda}{H_p}\right)^{4/3}
     \left(\frac{R_{\mathrm{spec}}\, F_{\mathrm{conv}}}{\rho\, c_p}\right)^{1/3} .

Link values are averaged onto the levels (one-sided at the domain ends). Because the
super-adiabaticity and with it :math:`K_{zz}` go to zero continuously at the
radiative-convective boundary, the value at the top of a zone is tiny and numerically erratic
while the solver iterates. Each zone is therefore **anchored** at the maximum convective
:math:`K_{zz}` within one pressure decade below the zone top, :math:`K_{\mathrm{rcb}}` (the
plain maximum would sit at the domain bottom, since :math:`K_{zz}` keeps growing with the scale
height). Every radiative level above the zone, up to the next zone, receives

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - ``radiative``
     - :math:`K_{zz}(p)`
     - Use
   * - ``constant`` (default)
     - :math:`K_{\mathrm{rcb}}`
     - continues the zone's mixing upward
   * - ``power_law``
     - :math:`K_{\mathrm{rcb}}\,(p/p_{\mathrm{rcb}})^{\mathrm{slope}}`
     - ``slope < 0``: mixing grows upward, irradiated planets (Lee 2024; AGNI uses -0.4);
       ``slope > 0``: decays towards ``min``, self-luminous objects whose radiative zone is
       stirred only by overshoot and gravity waves (Freytag et al. 2010: :math:`10^5` to
       :math:`10^7` cm\ :sup:`2` s\ :sup:`-1`)
   * - ``fixed``
     - ``value``
     - a prescribed value above every zone (retrieved photospheric values of T dwarfs are
       :math:`10^4` to :math:`10^8`, Miles et al. 2020; Mukherjee et al. 2022)

Levels below the deepest zone keep that zone's bottom value, and the whole profile is floored
at ``min`` (:math:`10^4` cm\ :sup:`2` s\ :sup:`-1` by default), which is also the value when
nothing is convective.

The ``mlt`` profile is intended for the mixing-length convection schemes. With the adjustment
schemes it still runs (with :math:`\alpha = 1`, no wall law and no ceiling), but a converged
adjustment leaves almost no super-adiabaticity, so the profile collapses to the floor.

Relaxation of the self-consistent profile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The profile is refreshed once per outer iteration from the committed temperature profile,
whose super-adiabaticity is the stiff quantity the Newton is still solving for; early in a run
the raw value swings by orders of magnitude between iterations, and the loop
:math:`K_{zz} \to` quench composition :math:`\to` opacity :math:`\to T \to K_{zz}` is a Picard
iteration with a loop gain of about 2.5 on the brown dwarf test case, which oscillates
undamped. The committed profile is therefore damped:

1. **Hold.** If the new raw profile finds no convective zone at all, the previous profile is
   kept: a momentarily stable profile in the middle of the iteration must not collapse
   :math:`K_{zz}` to the floor.
2. **Dead band.** With the logarithmic residual :math:`r_i = \ln K_{zz,i}^{\mathrm{raw}} - \ln
   K_{zz,i}^{\mathrm{prev}}`, the previous profile is kept if :math:`\max_i |r_i|` is below
   ``tolerance`` dex (0.05 by default). A change of 0.05 dex moves a quench temperature by about
   0.3 %, far inside the uncertainty of the quench timescales, but every committed change kicks
   the Newton residual up by an order of magnitude; freezing :math:`K_{zz}` lets the Newton
   finish quadratically.
3. **Aitken relaxation** (Küttler & Wall 2008) in log space: with the previous residual
   :math:`r^{\mathrm{prev}}`,

   .. math::

      \omega \leftarrow -\,\omega\, \frac{\sum_i r_i^{\mathrm{prev}} (r_i - r_i^{\mathrm{prev}})}
      {\sum_i (r_i - r_i^{\mathrm{prev}})^2}, \qquad \omega \in [0.05, 1],
      \qquad
      \ln K_{zz,i}^{\mathrm{new}} = \ln K_{zz,i}^{\mathrm{prev}} + \omega\, r_i .

   The user's ``relax`` is the initial weight. Aitken estimates the optimal weight
   :math:`1/(1 + \mathrm{gain})` from the last two residuals; a fixed weight of 0.5 converges
   but only linearly.

Both guards leave the converged fixed point untouched. The relaxation state persists across
:meth:`compute` calls of the same model object.

Operator splitting
~~~~~~~~~~~~~~~~~~

:math:`K_{zz}` is refreshed once per outer iteration from the committed profile and handed to
the chemistry and cloud modules. The trial evaluations inside the Newton corrector do not
update it: like the composition's temperature dependence and the cloud optics, :math:`K_{zz}`
is an operator-split (lagged) quantity, so the corrector never sees :math:`K_{zz}(T)` feedback.
