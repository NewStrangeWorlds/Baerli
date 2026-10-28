Radiative-convective equilibrium and the temperature iteration
==============================================================

This chapter describes the mathematical problem Bärli solves and the numerical schemes the
``solver`` option selects. The practical reference of the solver parameters and of the iteration
log is in :doc:`../user_guide/solvers`.

The problem
-----------

The radiative transfer returns, for the current temperature profile, the wavelength-integrated
net flux :math:`F_i = F_{\uparrow,i} - F_{\downarrow,i}` at every level (upward positive, the
stellar beam folded into the downward flux) and the mean intensity :math:`J_{\nu,i}`. Radiative
equilibrium states that no energy is deposited in any level, i.e. that the net flux is the same
everywhere:

.. math::
   :label: flux-constancy

   F_i = F_\star \quad \text{for all } i, \qquad
   F_\star = \begin{cases}
     \sigma \Teff^4 & \text{brown dwarf}\\
     \sigma \Tint^4 & \text{gas planet}\\
     0 & \text{terrestrial planet.}
   \end{cases}

Two local forms of the same statement are used by the solvers:

* the **heating rate**, :math:`\dd F/\dd p = 0` or equivalently the per-level net heating
  :math:`\sum_\nu w_\nu\, \dd F_\nu/\dd\tau_\nu = 0`, and
* **local radiative equilibrium**, emitted equals absorbed,

  .. math::

     \sum_\nu w_\nu\, \kappa_{\nu,i}\, B_\nu(T_i) = \sum_\nu w_\nu\, \kappa_{\nu,i}\, J_{\nu,i},

  where :math:`w_\nu` are the quadrature weights of the spectral grid and :math:`\kappa` the
  absorption coefficient of gas and clouds.

With convection, :eq:`flux-constancy` becomes :math:`F_{\mathrm{rad}} + F_{\mathrm{conv}} =
F_\star`: the mixing-length schemes add a smooth convective flux to the energy balance, the
adjustment schemes replace the flux condition in unstable layers by the adiabatic gradient
(:doc:`convection`). For a terrestrial planet the bottom level carries the surface energy balance
:math:`F_0 = 0`.

The residuals are made dimensionless with a flux scale chosen per object: :math:`F_\star` for a
brown dwarf, :math:`\mu_\star S + \Fint` for a gas planet and :math:`\max(1, \mu_\star S)` for a
terrestrial planet, :math:`S` being the instellation. For an irradiated planet :math:`\Fint` is
often tiny compared with the stellar flux the streams carry; judging the imbalance relative to
:math:`\Fint` alone would demand a stream accuracy no sampled-opacity solver can deliver, so the
honest scale is the flux the atmosphere actually transports.

Temperature Jacobians
---------------------

Both radiative transfer backends can return, together with the fluxes, the analytic derivatives
of the radiation field with respect to the level temperatures at frozen opacity (the
"Planck-only" response):

.. math::

   \frac{\partial F_i}{\partial T_j}, \qquad
   \frac{\partial}{\partial T_j}\left(\frac{\dd F}{\dd \tau}\right)_i, \qquad
   \sum_\nu w_\nu \kappa_{\nu,i} \frac{\partial J_{\nu,i}}{\partial T_j} .

The derivative with respect to the lower boundary temperature is folded into the column of
level 0. The Newton-type solvers (``ratio_ul``, ``flux_divergence``, ``ptc``) use these; the
relaxation solvers do not need them.

The outer iteration
-------------------

All three object classes run the same outer loop (:doc:`../getting_started/overview`): per
iteration the chemistry (with a refreshed :math:`K_{zz}`), the hydrostatic structure, the gas
opacities, the cloud modules and the cloud ramp are evaluated from the committed profile, then the
radiative transfer, then the corrector's :math:`T` update, then, for the relaxation schemes, the
convective adjustment and Ng acceleration. The loop ends when the corrector's residual falls
below ``convergence_threshold``, when the iteration stalls, or after ``max_iterations``.

Within one corrector step everything but the Planck source is **frozen**: composition, structure,
opacity, cloud optics and :math:`K_{zz}` are operator-split quantities that the inner Newton never
differentiates. The split is exact at the fixed point, since all of them are refreshed from the
committed profile once per outer iteration. The Jacobian is therefore the Planck-only response
above, and a corrector's inner iterations re-evaluate the true residual with the opacity
recomputed only where the scheme needs it.

The ``ratio_ul`` scheme
-----------------------

The default and recommended solver is a Newton method on the **collocated ratio residual with the
full Unsöld-Lucy correction**, derived from the solver of Wogan's clima model. Every level is an
unknown (with mixing-length convection) or every radiative level plus one degree of freedom per
convective zone (with the adjustment schemes).

The residual
~~~~~~~~~~~~

For a radiative level :math:`i` the residual is local radiative equilibrium written as a ratio,
augmented by the convective heating and by the Unsöld-Lucy term:

.. math::
   :label: ratio-residual

   g_i = \frac{\sum_\nu w_\nu \kappa_{\nu,i} J_{\nu,i} + q_i}{\sum_\nu w_\nu \kappa_{\nu,i} B_\nu(T_i)} - 1
   - \frac{\pi\, \delta B_i}{\sigma T_i^4},
   \qquad q_i = -\frac{1}{4\pi}\left(\frac{\dd F_c}{\dd z}\right)_i ,

with :math:`F_c` the mixing-length flux on the links (zero for the adjustment schemes). The
ratio form has an :math:`O(1)` eigenvalue for the grid-scale (Nyquist) mode, so its root carries
no checkerboard; a net-flux residual, being an odd angular moment, has a near-null Nyquist
eigenvalue and converges to a sawtooth. The diagonal :math:`-\langle \partial \ln B/\partial T
\rangle_\kappa` survives :math:`\kappa \to 0`, so transparent levels neither checkerboard nor
limit-cycle.

The **Unsöld-Lucy term** supplies the flux level and gradient the purely local ratio cannot see.
From the Eddington first moment :math:`\dd J/\dd\tau = 3F/(4\pi)` a flux error propagates into a
Planck error; the correction is integrated downward from the top of the atmosphere with the
Planck-mean layer optical depths :math:`\Delta\tau`,

.. math::

   \delta B_{n-1} = \frac{F_{n-1} + F_{c,n-2} - F_\star}{2\pi}, \qquad
   \delta B_i = \delta B_{i+1} + \frac{3}{4\pi}
   \left[\tfrac{1}{2}\big((F_i - F_\star) + (F_{i+1} - F_\star)\big) + F_{c,i}\right] \Delta\tau_{i+1} .

The flux enters as a one-sided cumulative integral, which has only the constant null mode and
suppresses Nyquist content, where a per-level flux difference amplifies it. The coefficients are
the Eddington closure, not tuning parameters. With mixing-length convection the conserved column
quantity is :math:`F_{\mathrm{rad}} + F_c`, hence the convective terms.

The remaining rows are flux balances. The bottom level of a self-luminous object or the surface
of a terrestrial planet carries

.. math::

   g_0 = \frac{F_0 + F_{c,0} - F_\star}{c_{\mathrm{eff},0}}, \qquad
   c_{\mathrm{eff},i} = \frac{\Delta p_i}{g}\, c_{p,i},

the heat capacity of the level's pressure slab; the surface loses :math:`F_{c,0}` as sensible
heat into the convective column. With the adjustment schemes each contiguous convective zone
:math:`[l, u]` is one unknown (its bottom temperature), its interior is slaved to the adiabat,
:math:`T_k = T_l \prod (p_k/p_{k-1})^{\nabla_{\mathrm{ad}}}`, and its row is the energy budget of
the whole zone,

.. math::

   g = \frac{F_u + F_{c,u} - F^{\mathrm{in}}_l}{F_{\mathrm{norm}}}, \qquad
   F^{\mathrm{in}}_l = \begin{cases} F_\star & l = 0 \\ F_{l-1} & \text{else,}\end{cases}

normalised by :math:`F_{\mathrm{norm}} = \max(1, |F_\star|, |F_\downarrow(\mathrm{TOA})|)` rather
than by the zone's heat capacity, which would make a small imbalance invisible to the Newton. The
first radiative level above a slaved zone gets the backward flux difference
:math:`(F_i - F_{i-1})/c_{\mathrm{eff},i}` as its row, because nothing else ties the local
solution to the top of the adiabat.

Band closure
~~~~~~~~~~~~

The sampled spectral points miss the Doppler-narrow line cores. A band-closure correction adds,
per 30 cm\ :sup:`-1` band, the difference between the native line-opacity integral and its
sampled representation, weighted with a frozen band-mean radiation field and an escape
probability :math:`\exp(-\tau_{\mathrm{core}})`, with a cubic taper that switches it off for
:math:`h c \nu / k T \ge 10`. The band-mean field is frozen over the inner solve so that residual
and Jacobian stay a consistent pair.

The inner solve
~~~~~~~~~~~~~~~

The Jacobian of :eq:`ratio-residual` is assembled from the Planck-only backend Jacobians, with
the same one-sided integral applied to the flux Jacobian for the Lucy term and the analytic
mixing-length link derivatives, and rebuilt at every accepted step. The nonlinear system is
solved with the **NLEQ-ERR affine-covariant damped Newton** method of Deuflhard:

1. ordinary Newton correction :math:`J\, \delta x = -g` (dense LU);
2. trial step :math:`x + \lambda\, \delta x` with the damping factor starting at
   :math:`\min(1, 2\lambda_{\mathrm{prev}})`, the **true** residual at the trial point and the
   simplified correction :math:`\overline{\delta x} = -J^{-1} g(x_t)`;
3. the natural monotonicity test :math:`\|\overline{\delta x}\| / \|\delta x\| \le 1 - \lambda/4`
   accepts the step, otherwise :math:`\lambda` is halved down to :math:`\lambda_{\min} = 10^{-4}`;
4. at :math:`\lambda_{\min}` a step is taken only if it does not increase the residual norm, and
   after three consecutive such floor steps the inner solve returns so that the outer loop can
   relinearise with fresh opacities;
5. the inner solve terminates when the scaled step falls below :math:`10^{-7}`, committing that
   last step.

Affine covariance makes the globalisation invariant to the scale gap between the ratio rows
(:math:`\sim 10^{-2}`) and the flux rows (:math:`\sim 10^3`), and the monotonicity test rejects
any overshoot that drives the denominator of the ratio towards its pole. After the inner solve a
final evaluation with the opacity recomputed makes the model state consistent for the driver. No
post-hoc smoothing is applied: a filter after the solve mutates a profile the Newton never saw,
and the two would fight.

Convective zones with the adjustment schemes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``dry`` or ``moist`` convection the convective mask is detected on a smoothed probe copy of
the profile and may only move between settled inner solves, by at most one layer per zone
boundary, with a dead band of ``mask_band`` levels (2 for terrestrial planets, 0 for self-luminous
objects, whose deep radiative band at :math:`\Delta\tau \gg 1` locks a checkerboard in if the
boundary is placed one level off). Near convergence the unsmoothed link above the zone top is
tested with a 2 % super-adiabaticity threshold (picaso's grow-only ratchet), an anti-overshoot
rule retreats from a cold inversion above the zone, and a lockout breaks grow/retract cycles. The
mixing-length schemes replace this machinery entirely.

Convergence metric
~~~~~~~~~~~~~~~~~~

Two numbers are computed at the committed point:

* the **flux-conservation residual**, the largest total-flux imbalance of any row relative to
  :math:`F_{\mathrm{norm}}`, which is the honest error bar; its floor is not zero, since a
  multi-stream sampled-opacity solver leaves a per-layer flux divergence of order :math:`10^{-3}`
  in the deep atmosphere;
* the **own residual**, :math:`\max_i |g_i|`, which the ratio scheme drives to round-off.

The reported residual is :math:`\max(\text{own residual}, \max |\Delta T/T| \text{ of the inner
solve})`, set to 1 after a large mask move and floored at 0.5 while the mixing-length ramp is
running. It is tested against ``convergence_threshold``; the TOA net-flux error of a self-luminous
object is reported separately at convergence because the optically thin skin is excluded from the
gate.

The ``flux_divergence`` scheme
------------------------------

The same class, mask and Jacobian machinery with clima's heat-capacity-weighted backward flux
difference as the residual of every radiative level,

.. math::

   g_i = \frac{F_i - F_{i-1} + F^{\mathrm{lvl}}_{c,i} - F^{\mathrm{lvl}}_{c,i-1}}{c_{\mathrm{eff},i}},

the heating rate :math:`\dd T/\dd t` of the slab, and the zone rows divided by the zone's heat
capacity. No Lucy term is needed, since the flux difference carries level and gradient by
itself. The inner solver is MINPACK's ``hybrj`` (Powell's dogleg trust region with rank-one
updates, through Eigen's port). The root conserves flux exactly per level (:math:`10^{-7}` where
the ratio scheme reaches :math:`6 \times 10^{-5}`) but carries a grid-scale sawtooth, the Nyquist
degeneracy of an odd angular moment. Use it when per-level flux conservation matters more than a
smooth profile. It requires an adjustment convection scheme; with the mixing-length flux law it
was measured not to converge.

The ``ptc`` scheme
------------------

The legacy linearised flux Newton with pseudo-transient continuation, after picaso's formulation.
A grey optical depth from the wavenumber-mean extinction defines a regime weight
:math:`\zeta = \tau/(\tau + 1)`; the optically thin skin (:math:`\zeta < 0.1`) is set by direct
inversion of local radiative equilibrium, a single contiguous deep convective zone grows from the
bottom and is slaved to the adiabat, and the remaining levels are flux unknowns.

* **Self-luminous objects**: an inner Newton iterates the *linearised* flux
  :math:`F_i + \sum_j (\partial F_i / \partial T_j)(T_j - T_j^{\mathrm{cur}}) = F_\star` with a
  backtracking line search and a Laplacian smoothing term on adjacent unknowns that suppresses the
  checkerboard near-null mode and vanishes at the inner fixed point. No radiative transfer solve is
  needed inside the loop.
* **Terrestrial planets**: every radiative level is updated by its local-RE step, the whole
  column is shifted uniformly to drive :math:`F(\mathrm{TOA}) \to 0` (an Unsöld-Lucy level
  correction), and the single flux unknown is the surface, advanced by pseudo-transient
  continuation :math:`(J + C_0/\Delta t)\, s = -(F_0 - F_\star)` with Deuflhard's second-order
  step-size rule. This path converges on the temperature change once a convective zone exists,
  because the troposphere can never reach radiative flux balance.

The driver's ``max_change`` clip applies to this scheme.

The relaxation schemes
----------------------

``time_stepping``
~~~~~~~~~~~~~~~~~

Explicit relaxation on the flux divergence with the local radiative timescale,

.. math::

   \Delta T_i = \gamma\, \frac{p_i}{4\sigma T_i^3}\left(\frac{\dd F}{\dd p}\right)_i, \qquad
   \Delta T_0 = \gamma\, \frac{F_\star - F_0}{4\sigma T_0^3}\ \ (\text{self-luminous bottom anchor}),

with the centered flux divergence in the interior and one-sided differences at the ends. The
bottom anchor follows HELIOS and uses the bottom net flux, not the emergent one: for an optically
thick base the emergent flux is decoupled from the deep temperature. The driver clips the change
to ``max_change``, applies the convective adjustment, the Ng acceleration and, for terrestrial
planets, the explicit surface model (a 50 cm slab ocean stepped by 1000 s) with a Shapiro filter
on the profile. Convergence is judged on :math:`\max |\Delta T / T|`.

``time_stepping_lre``
~~~~~~~~~~~~~~~~~~~~~

The same step blended with an opacity-weighted local radiative equilibrium update,

.. math::

   T_{\mathrm{LRE},i} = T_i \left(\frac{\int \kappa_{\nu,i} J_{\nu,i}\, \dd\nu}
   {\int \kappa_{\nu,i} B_\nu(T_i)\, \dd\nu}\right)^{1/4}, \qquad
   T_i \leftarrow T_i + \mathrm{lre\_fraction}\, (T_{\mathrm{LRE},i} - T_i),

which is well conditioned where the time stepping is slow (the optically thin skin).

``helios``
~~~~~~~~~~

A port of the pseudo-time relaxation of HELIOS (Malik et al. 2017, 2019). Each level gets a
dimensionless imbalance :math:`r_i`, by default the backward flux difference
:math:`(F_{i-1} - F_i)/F_{\mathrm{scale}}` with the bottom anchor :math:`(F_\star -
F_0)/F_{\mathrm{scale}}` (``stencil = "centered"`` and ``"forward"`` and the local heating
``residual = "heating"`` are alternatives), and moves by

.. math::

   \Delta T_i = s_i\, \mathrm{sgn}(r_i)\, |r_i|^{q}, \qquad q = 0.1,

capped at ``max_step``. With :math:`q \approx 0.1` the step is nearly residual-blind: a deep layer
with a tiny imbalance moves by about as many kelvin per iteration as a thin one, which beats the
stiffness between the optically thin skin and the diffusive deep atmosphere. Every
``adapt_interval`` iterations a level that oscillated (moved less than half of what its steps
summed to) has its step :math:`s_i` divided by ``step_shrink``, otherwise multiplied by
``step_grow``. Convergence is judged on the largest flux imbalance of the non-convective levels.
This is not a physical time step. Measured on the gas-planet test case, the backward stencil
converged in about 300 iterations to within 5 K of the Newton solution in the deep atmosphere,
the centered stencil floored with a growing checkerboard and the forward stencil diverged.

Ng acceleration
~~~~~~~~~~~~~~~

The relaxation schemes may be accelerated every ``ng_interval`` iterations with the
extrapolation of Ng (1974) in the form of Olson, Auer & Buchler (1986): from four stored
iterates :math:`T_0 \ldots T_3` with differences :math:`d_k = T_k - T_{k-1}` and second
differences :math:`D_1 = d_2 - d_1`, :math:`D_2 = d_3 - d_2`, the normal equations give
coefficients :math:`c_1, c_2` and the accelerated profile :math:`T^\ast = T_3 - c_2 d_3 - c_1
d_2`. The extrapolation is skipped when the system is nearly singular or would produce
temperatures below 1 K.

Stall detection
---------------

A corrector that proposes no change at all (:math:`\max|\Delta T/T| < 10^{-9}`) for three
iterations while the residual stays above the threshold has hit a floor it cannot pass; a residual
that has not improved by 1 % over the last 25 iterations signals a damped Newton stuck at its
minimum damping. Both end the :meth:`compute` call with ``False``. The window is longer than the
mixing-length ramp, which legitimately holds the residual at its floor for about twelve iterations,
and a changed cloud clears the history.

.. _theory-initial-profiles:

Initial profiles
----------------

``adiabat``
   The adiabat of the active convection scheme from the surface temperature, reduced by 2 % and
   floored at the stratosphere temperature (:doc:`convection`). Clima's start; with the stable
   bias it lets the mixing-length corrector skip its ramp.

``milne``
   The grey Milne-Eddington solution with Rosseland optical depth :math:`\tau_R = \kappa_R\, p /
   g`,

   .. math::

      T^4 = \tfrac{3}{4}\, \Teff^4\, \big(q(\tau_R) + \tau_R\big),

   with the Hopf function :math:`q(\tau)` approximated by a rational fit in :math:`\log_{10}
   \tau` between the limits :math:`q(0) = 0.5774` and :math:`q(\infty) = 0.7104`. Qualitatively
   wrong for an irradiated planet.

``const``
   Isothermal. For an irradiated gas planet an isothermal start at the equilibrium temperature
   puts the deep atmosphere hundreds of kelvin too low, and the first pass is a long
   damped-Newton crawl.

``guillot``
   The analytic irradiated profile of Guillot (2010), the recommended start for irradiated gas
   planets. With :math:`\tau = \kappa_{\mathrm{ir}}\, p / g` and :math:`\gamma =
   \kappa_{\mathrm{vis}}/\kappa_{\mathrm{ir}}`, the isotropic form (Guillot's Eq. 29) with the
   redistribution factor :math:`f`,

   .. math::

      T^4 = \tfrac{3}{4}\Tint^4\left(\tfrac{2}{3} + \tau\right)
      + \tfrac{3}{4} T_{\mathrm{irr}}^4\, f \left[\tfrac{2}{3} + \frac{1}{\gamma\sqrt{3}}
      + \left(\frac{\gamma}{\sqrt{3}} - \frac{1}{\gamma\sqrt{3}}\right) e^{-\gamma\tau\sqrt{3}}\right],

   or the collimated-beam form (Eq. 27) with the zenith-angle cosine :math:`\mu`,

   .. math::

      T^4 = \tfrac{3}{4}\Tint^4\left(\tfrac{2}{3} + \tau\right)
      + \tfrac{3}{4} T_{\mathrm{irr}}^4\, \mu \left[\tfrac{2}{3} + \frac{\mu}{\gamma}
      + \left(\frac{\gamma}{3\mu} - \frac{\mu}{\gamma}\right) e^{-\gamma\tau/\mu}\right],

   floored at 50 K. :math:`f = 1/4` is the global average, :math:`f = 1/2` the dayside average;
   :math:`T_{\mathrm{irr}} = (S/\sigma)^{1/4}` for the instellation :math:`S`.
