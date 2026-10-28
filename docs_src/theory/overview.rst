Overview of the model
=====================

Bärli is a one-dimensional radiative-convective equilibrium model: it finds the temperature
profile :math:`T(p)` of a plane-parallel atmosphere for which the net energy flux through every
level equals the flux the object has to carry, with radiation and convection as the transport
mechanisms and with the chemical composition, the clouds and the vertical mixing evaluated
self-consistently for that profile. The theory chapters document what is actually implemented:

* :doc:`atmosphere`: the pressure grid, the hydrostatic structure and the thermodynamic data
  (heat capacities, adiabatic gradients, water saturation);
* :doc:`radiative_transfer`: the monochromatic multi-stream solvers, their boundary conditions
  and the analytic temperature Jacobians;
* :doc:`opacities`: opacity sampling of line cross sections, collision-induced absorption,
  continua and Rayleigh scattering;
* :doc:`chemistry`: equilibrium chemistry with condensation, prescribed compositions, the
  quench approximation and the coupling to chemical kinetics;
* :doc:`clouds`: prescribed and microphysical (Ackerman & Marley) clouds with Mie optics;
* :doc:`convection`: mixing-length convection, convective adjustment and the eddy diffusion
  coefficient;
* :doc:`temperature_iteration`: the Newton and relaxation schemes that drive the profile to
  equilibrium, and the initial profiles.

The state of a model
--------------------

A model is defined on :math:`N` pressure levels, equidistant in :math:`\log p` from the bottom
(level 0) to the top. Its state is the temperature at every level, the number densities of the
species of Bärli's species table, the mean molecular weight, the altitude and scale height from
hydrostatic equilibrium, the eddy diffusion coefficient, the condensed phase of the equilibrium
chemistry and the structure of the cloud modules. From this state the opacities, the radiation
field and the emergent spectrum follow.

The fixed-point problem
-----------------------

The radiative transfer maps the state to a radiation field, and the radiation field, through the
energy balance, to a temperature correction; the chemistry, the clouds and the mixing map the
temperature back to a new state. Bärli organises this as an **outer iteration** in which all
state-dependent quantities are refreshed once from the committed profile, and an **inner solve**
in which only the temperature changes while composition, opacity, clouds and :math:`K_{zz}` are
frozen (operator splitting). The Newton solvers linearise the frozen-opacity radiation field with
the analytic Jacobians of the radiative transfer; the relaxation solvers step the temperature with
a pseudo-time. The split is exact at the fixed point, so the converged state is a genuine
self-consistent solution, and the freezing is what keeps the inner problem smooth enough for a
damped Newton.

Several feedback loops around the outer iteration have gains above unity and would oscillate if
iterated naively: the eddy diffusion coefficient derived from the convective flux law, the quench
pressures of the disequilibrium chemistry, the structure of a microphysical cloud, and the
composition returned by a coupled kinetics code. Each of these is damped with the same device, an
Aitken dynamic relaxation in log space with a dead band below which a change is not committed.
This is a recurring design pattern in Bärli (:ref:`theory-kzz`, :doc:`chemistry`,
:doc:`clouds`), and the ramps with which the mixing-length flux law and the cloud opacity are
switched on at the start of a cold run belong to the same family of globalisation devices.

Heritage
--------

Many of the schemes are ports or re-implementations of methods from other atmosphere codes,
acknowledged in the source: the collocated residual and the convective-zone handling of Wogan's
clima model, the pseudo-time relaxation of HELIOS (Malik et al. 2017, 2019), the linearised flux
Newton and the convective-boundary ratchet of picaso, the mixing-length formulation and the
velocity-scaled :math:`K_{zz}` of AGNI, the cloud layer solver and fall-velocity routines of virga
(Ackerman & Marley 2001; Rooney et al. 2022), the Guillot (2010) initial profile from BeAR, and the
opacity-sampling machinery that descends from the BeAR retrieval code. The equilibrium chemistry
is FastChem, the radiative transfer DisORT++ and AddingDoublingRT, the Mie optics LX-MIE, all
from the same software family. The bibliography collects the references
(:doc:`../reference/bibliography`).
