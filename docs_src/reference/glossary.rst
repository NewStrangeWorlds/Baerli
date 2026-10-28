Glossary
========

.. glossary::
   :sorted:

   level
      One of the ``nb_grid_points`` pressure points of the grid. Temperatures, number
      densities and fluxes are defined on levels. Level 0 is the bottom of the atmosphere.

   layer
      The interval between two adjacent levels. Optical depths, cloud properties and the
      single-scattering albedo are layer quantities; there are ``nb_grid_points - 1`` layers.

   module spec
      The description of a pluggable component: a type keyword plus named parameters, written
      as ``"type"``, ``("type", {parameters})`` or ``{"type": ..., parameters}``.

   RCE
      Radiative-convective equilibrium: the state in which the net energy flux through every
      layer equals the internal flux (or zero for a planet without internal heat), carried by
      radiation in stable layers and by radiation plus convection in unstable ones.

   RCB
      The radiative-convective boundary, the top of a convective zone.

   corrector
      The temperature-correction scheme selected by the ``solver`` spec: the object that turns
      the current radiation field into a new temperature profile.

   residual
      The measure of the departure from equilibrium that a corrector drives to zero, and the
      convergence metric printed in the iteration log. Each solver defines its own.

   easy start
      The homotopy with which the mixing-length corrector begins a cold run: the convective
      flux law is ramped in over the first iterations. Skipped on a warm start.

   warm start
      A :meth:`compute` call that declares the current profile a converged solution of a
      slightly changed problem, so that ramps are skipped and the full cloud is applied at
      once.

   opacity sampling
      Evaluating the radiative transfer on a subset of the high-resolution wavenumber list of
      line-by-line cross sections, with the exact opacity at each sampled point; as opposed to
      correlated-*k* or band models.

   instellation
      The stellar flux at the planet's orbit, integrated over wavelength and measured
      perpendicular to the beam, in erg cm\ :sup:`-2` s\ :sup:`-1`.

   zenith angle
      In Bärli always the **cosine** :math:`\mu_\star` of the angle between the stellar beam and
      the vertical; the flux entering the atmosphere is :math:`\mu_\star S`.

   fsed
      The sedimentation efficiency of the Ackerman & Marley cloud model: the ratio of the
      mass-weighted sedimentation velocity of the condensate to the convective velocity scale.

   Kzz
      The eddy diffusion coefficient in cm\ :sup:`2` s\ :sup:`-1` that parametrises vertical
      mixing in the quench chemistry, the clouds and the kinetics coupling.

   quench
      The freezing of a chemical abundance at the pressure where the mixing timescale becomes
      shorter than the chemical timescale.

   rainout
      FastChem's condensation mode in which the elements locked in condensates at a level are
      removed from the gas phase of all levels above it.
