Radiative transfer
==================

Bärli computes the radiation field monochromatically at every point of the spectral grid with one
of two plane-parallel multi-stream solvers, DISORT++ or AddingDoublingRT, and integrates the
results over wavenumber. Both solvers return, besides the fluxes and the mean intensity, the
analytic derivatives of the radiation field with respect to the level temperatures that the Newton
solvers need. The user-facing options are in :doc:`../user_guide/radiative_transfer`.

The monochromatic problem
-------------------------

For every wavenumber :math:`\tilde\nu` of the grid the atmosphere is a stack of :math:`N-1`
homogeneous layers, each with an optical depth :math:`\Delta\tau_j`, a single-scattering albedo
:math:`\omega_j` and an azimuthally averaged phase function. The solvers integrate

.. math::

   \mu \frac{\dd I(\tau, \mu)}{\dd\tau} = I(\tau, \mu) - (1 - \omega) B_{\tilde\nu}(T(\tau))
   - \frac{\omega}{2} \int_{-1}^{1} p(\mu, \mu')\, I(\tau, \mu')\, \dd\mu'
   - \frac{\omega}{4\pi}\, S\, p(\mu, -\mu_0)\, e^{-\tau/\mu_0}

with the Planck source linear in optical depth inside each layer and a collimated stellar beam of
flux :math:`S` (per unit area normal to the beam) incident at :math:`\mu_0`. Because the opacities
are sampled, the Planck function is evaluated at the single wavenumber of each point,

.. math::

   B_{\tilde\nu}(T) = \frac{2 h c^2 \tilde\nu^3}{\exp(h c \tilde\nu / k T) - 1}
   \qquad [\mathrm{erg\,s^{-1}\,cm^{-2}\,sr^{-1}\,(cm^{-1})^{-1}}].

Layer properties
~~~~~~~~~~~~~~~~

The gas absorption and scattering coefficients :math:`\kappa^{\mathrm{abs}}_i`,
:math:`\kappa^{\mathrm{sca}}_i` (cm\ :sup:`-1`) are defined on the levels and integrated across
each layer with the trapezoidal rule; the cloud modules provide per-layer optical depths,
albedos and asymmetry parameters directly (:doc:`clouds`):

.. math::

   \tau^{\mathrm{gas}}_{\mathrm{abs}, j} = (z_{j+1} - z_j)\, \frac{\kappa^{\mathrm{abs}}_j + \kappa^{\mathrm{abs}}_{j+1}}{2},
   \qquad
   \Delta\tau_j = \tau^{\mathrm{gas}}_{\mathrm{abs}} + \tau^{\mathrm{gas}}_{\mathrm{sca}}
   + \tau_c (1 - \omega_c) + \tau_c\, \omega_c ,

.. math::

   \omega_j = \frac{\tau^{\mathrm{gas}}_{\mathrm{sca}} + \tau_c\, \omega_c}{\Delta\tau_j}
   \ (\le 0.999999), \qquad
   g_j = \frac{\tau_c\, \omega_c}{\tau^{\mathrm{gas}}_{\mathrm{sca}} + \tau_c\, \omega_c}\, g_c .

A layer with :math:`g_j > 0` uses a Henyey-Greenstein phase function with that asymmetry
parameter (Legendre moments :math:`g^l`); a layer without cloud scattering uses the Rayleigh phase
function. As soon as a layer contains cloud scattering, the whole layer scatters with the single
Henyey-Greenstein function of the scattering-weighted :math:`g`. No delta-M scaling is requested
from either library.

Boundary conditions
~~~~~~~~~~~~~~~~~~~

* **Top**: no incident thermal radiation. For irradiated objects the stellar beam enters with the
  spectral flux :math:`S_{\tilde\nu}` normal to the beam and the direction cosine
  :math:`\mu_0 =` ``zenith_angle``, so the direct flux at the top of the atmosphere is
  :math:`\mu_0 S_{\tilde\nu}`.
* **Bottom, self-luminous objects** (brown dwarfs, gas planets): the diffusion lower boundary of
  a semi-infinite atmosphere. The upwelling intensity at the bottom interface is

  .. math::

     I^\uparrow(\mu) = B(T_{\mathrm{bot}}) + \mu\, \frac{\dd B}{\dd\tau}, \qquad
     \frac{\dd B}{\dd\tau} = \frac{B(T_{\mathrm{bot}}) - B(T_{N-2})}{\Delta\tau_{\mathrm{last}}},

  the gradient of the Planck function across the lowest layer. This is what lets the internal
  flux :math:`\sigma \Tint^4` escape through the bottom and anchors the deep temperature.
* **Bottom, terrestrial planets**: a Lambertian surface with the spectral albedo
  :math:`A_{\tilde\nu}` of the surface module, emitting :math:`(1 - A_{\tilde\nu})
  B(T_{\mathrm{surf}})`. In the Newton solvers the surface is tied to the bottom level
  (:math:`T_{\mathrm{surf}} = T_0`, no independent degree of freedom; the surface derivative is
  folded into the bottom-level column of the Jacobian).

For each spectral point the thermal source is switched on only if :math:`B_{\tilde\nu}(T_{\max})
> 10^{-35}` in cgs units, which avoids underflow in the far ultraviolet of cold atmospheres.

The two backends
----------------

DISORT++
~~~~~~~~

`DisORT++ <https://github.com/NewStrangeWorlds/DisORT>`_ is a C++17 implementation of the
discrete-ordinate method of Stamnes et al. (1988). Bärli uses its flux-only solver with a
compile-time number of streams (4, 6, 8, 10, 12, 14, 16, 32 or 64 in total over both
hemispheres), fixed-size Eigen matrices, a Lambertian surface, thermal emission and the
analytic temperature Jacobians. ``nb_streams`` is the total number of streams; the default is 4.

AddingDoublingRT
~~~~~~~~~~~~~~~~

`AddingDoublingRT <https://github.com/NewStrangeWorlds/AddingDoublingRT>`_ implements the
matrix-operator (adding-doubling) method of Plass, Kattawar & Catchings (1973) with the
normalisation of Hansen (1971), with Gauss-Legendre quadrature per hemisphere and analytic
temperature Jacobians. ``nb_streams`` is the number of quadrature points **per hemisphere**
(at least 2; specialised kernels exist for 2, 4, 8, 16 and 32); the default is 2. It splits the
net flux exactly into a thermal and a stellar part and bands the Jacobian at an effective
absorption depth of 40.

Both backends agree to the accuracy of their angular quadrature; the adding-doubling solver with
two streams per hemisphere is the faster choice for terrestrial and gas-planet models with many
spectral points, DISORT with four streams the reference for brown dwarfs.

Outputs
-------

Per spectral point and level the wrappers return the upward and downward fluxes (the attenuated
direct beam folded into the downward flux), the net flux :math:`F = F_\uparrow - F_\downarrow`
(upward positive) and the mean intensity :math:`J = \frac{1}{4\pi}\oint I\, \dd\Omega` including
the direct beam, converted from the libraries' SI units to erg cm\ :sup:`-2` s\ :sup:`-1` per
cm\ :sup:`-1`. The spectral integrals over the grid are trapezoidal sums with the weights
:math:`w_k`:

.. math::

   F^{\mathrm{tot}}_i = \sum_k w_k\, F_{k,i}, \qquad
   J^{\mathrm{tot}}_i = \sum_k w_k\, J_{k,i}, \qquad
   \left.\frac{\dd F}{\dd p}\right|_i = \frac{F^{\mathrm{tot}}_{i+1} - F^{\mathrm{tot}}_{i-1}}{p_{i+1} - p_{i-1}}

(centered in the interior, one-sided at the ends). The emergent **spectrum** is the upward flux at
the top of the atmosphere converted to per micron, :math:`F_\lambda = F_{\tilde\nu}\,
10^4/\lambda^2` with :math:`\lambda` in micron.

The solvers also return their own flux divergence per level,

.. math::

   \left.\frac{\partial F}{\partial\tau}\right|_i = 4\pi\, (1 - \omega)\, (J_i - B_i),

whose weighted sum over the grid is the ``net_heating``; it is used instead of a
:math:`\kappa (J - B)` reconstruction so that zero net heating is exactly equivalent to a constant
solver net flux.

Temperature Jacobians
---------------------

With ``compute_jacobian`` set the backends return, analytically and at frozen opacity, the
derivatives of the radiation field with respect to every level temperature, which the wrappers
integrate over the grid:

.. math::

   \frac{\partial F^{\mathrm{tot}}_i}{\partial T_j} = \sum_k w_k
   \left(\frac{\partial F_{\uparrow, k, i}}{\partial T_j} - \frac{\partial F_{\downarrow, k, i}}{\partial T_j}\right),
   \qquad
   \frac{\partial}{\partial T_j}\left(\frac{\partial F}{\partial\tau}\right)_i, \qquad
   \sum_k w_k\, \kappa_{k,i}\, \frac{\partial J_{k,i}}{\partial T_j},

the third one weighted with the total absorption coefficient of gas and clouds, because that is
the weight the ratio residual uses, so that :math:`\kappa` cancels in its diagonal. The derivative
with respect to the bottom boundary temperature is folded into the column of level 0. The chain
rule through the Planck function uses :math:`\dd B/\dd T = B\, (x/T)/(1 - e^{-x})` with
:math:`x = h c \tilde\nu / k T`. Composition, structure and opacity are frozen in these
derivatives (:doc:`temperature_iteration`).

Parallelisation
---------------

The spectral loop is parallelised with OpenMP; each thread owns one solver configuration and
workspace, and the per-thread partial sums of the integrated quantities and Jacobians are reduced
in a fixed thread order for run-to-run reproducibility. The radiative transfer saturates all cores,
so Bärli models should be run one at a time.
