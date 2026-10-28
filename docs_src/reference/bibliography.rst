Bibliography
============

Methods and reference codes
---------------------------

* Ackerman, A. S. & Marley, M. S. 2001, ApJ 556, 872: *Precipitating condensation clouds in
  substellar atmospheres* (the eddy-sedimentation cloud model, Eqs. 13 to 18 and Eq. B2).
* Batalha, N. et al., *virga*: the Python implementation of the Ackerman & Marley model whose
  layer solver, fall-velocity routine and radius grid Bärli follows.
* Rooney, C. M., Batalha, N. E., Gao, P. & Marley, M. S. 2022, ApJ 925, 33: closed-form layer
  recursion and the exponential :math:`f_{\mathrm{sed}}` profile.
* Pruppacher, H. R. & Klett, J. D. 1978, *Microphysics of Clouds and Precipitation*: the
  drag-coefficient fit of the fall velocity for :math:`1 < \mathrm{Re} \le 10^3`.
* Rosner, D. E. 2000, *Transport Processes in Chemically Reacting Flow Systems*: Lennard-Jones
  viscosity correction.
* Kitzmann, D. & Heng, K. 2018, MNRAS 475, 94: LX-MIE and the compilation of refractive indices
  of condensates (the sources of the individual data sets are given in the headers of the files
  in ``data/refractive_indices/``).
* Stock, J. W., Kitzmann, D., Patzer, A. B. C. & Sedlmayr, E. 2018, MNRAS 479, 865; Stock, J. W.,
  Kitzmann, D. & Patzer, A. B. C. 2022, MNRAS 517, 4070; Kitzmann, D., Stock, J. W. & Patzer,
  A. B. C. 2024, MNRAS 527, 7263: FastChem, FastChem 2 and FastChem Cond (the rainout
  approximation, Sect. 3.6 of the third paper).
* Asplund, M., Grevesse, N., Sauval, A. J. & Scott, P. 2009, ARA&A 47, 481: solar element
  abundances.
* Zahnle, K. J. & Marley, M. S. 2014, ApJ 797, 41: quench timescales (Eqs. 12, 13, 14, 32, 40,
  43, 44).
* Manabe, S. & Wetherald, R. T. 1967, J. Atmos. Sci. 24, 241: the relative-humidity profile
  (Eq. 2).
* Murphy, D. M. & Koop, T. 2005, QJRMS 131, 1539: saturation vapour pressure of water over ice
  (Eq. 7) and liquid (Eq. 10).
* McBride, B. J., Zehe, M. J. & Gordon, S. 2002, NASA/TP-2002-211556: the NASA polynomials of
  the heat capacities.
* Stamnes, K., Tsay, S.-C., Wiscombe, W. & Jayaweera, K. 1988, Appl. Opt. 27, 2502: DISORT, the
  discrete-ordinate method implemented by DisORT++.
* Plass, G. N., Kattawar, G. W. & Catchings, F. E. 1973, Appl. Opt. 12, 314; Hansen, J. E. 1971,
  J. Atmos. Sci. 28, 1400; Wiscombe, W. J. 1977, J. Atmos. Sci. 34, 1408: the matrix-operator
  (adding-doubling) method and the delta-M scaling implemented by AddingDoublingRT.
* Grimm, S. L. & Heng, K. 2015, ApJ 808, 182; Grimm, S. L. et al. 2021, ApJS 253, 30: HELIOS-k,
  the opacity calculator whose output format Bärli reads.
* Helling, Ch. & Jørgensen, U. G. 1998, A&A 337, 477: the composite-Planck covering distribution
  of spectral points.
* John, T. L. 1988, A&A 193, 189: H\ :sup:`-` bound-free and free-free opacity fits.
* Gruszka, M. & Borysow, A. 1997, Icarus 129, 172; Gruszka, M. & Borysow, A. 1998, Mol. Phys.
  93, 1007; Baranov, Y. I., Lafferty, W. J. & Fraser, G. T. 2004, J. Mol. Spectrosc. 228, 432:
  CO\ :sub:`2`-CO\ :sub:`2` collision-induced absorption (the port by R. Wordsworth from the LMD
  generic GCM).
* Mlawer, E. J. et al. 2012, Phil. Trans. R. Soc. A 370, 2520: the MT_CKD water-vapour
  continuum (version 4.2 is used).
* Lee, H.-W. & Kim, H. I. 2004, MNRAS 347, 802: Rayleigh scattering by atomic hydrogen.
* Sneep, M. & Ubachs, W. 2005, JQSRT 92, 293; Thalman, R., Zarzana, K. J., Tolbert, M. A. &
  Volkamer, R. 2014, JQSRT 147, 171: refractive indices and King factors of the Rayleigh cross
  sections.
* Gueymard, C. A. 2004, Solar Energy 76, 423: the solar spectrum in ``data/stellar_spectra``.
* Guillot, T. 2010, A&A 520, A27: the analytic irradiated temperature profile (Eqs. 27 and 29).
* Deuflhard, P. 2004, *Newton Methods for Nonlinear Problems* (Springer), and Deuflhard, P.
  2002, ZIB report 02-14: the NLEQ-ERR affine-covariant damped Newton method and the
  pseudo-transient continuation step-size rule.
* Ng, K.-C. 1974, J. Chem. Phys. 61, 2680; Olson, G. L., Auer, L. H. & Buchler, J. R. 1986,
  JQSRT 35, 431: Ng acceleration.
* Küttler, U. & Wall, W. A. 2008, Comput. Mech. 43, 61: Aitken dynamic relaxation.
* Moré, J. J., Garbow, B. S. & Hillstrom, K. E. 1980, MINPACK (``hybrj``): the trust-region
  solver of the ``flux_divergence`` scheme, through Eigen's port.
* Unsöld, A. 1955, *Physik der Sternatmosphären*; Lucy, L. B. 1964, SAO Special Report 167: the
  Unsöld-Lucy temperature correction.
* Blackadar, A. K. 1962, J. Geophys. Res. 67, 3095: the mixing-length wall law.
* Wogan, N. F. et al. 2025, *clima* / *photochem*: the AdiabatClimate solver whose residual,
  convective-zone handling and adiabat initialisation the ``ratio_ul`` and ``flux_divergence``
  schemes re-implement.
* Malik, M. et al. 2017, AJ 153, 56; Malik, M. et al. 2019, AJ 157, 170: HELIOS, the
  pseudo-time relaxation of the ``helios`` scheme and the bottom flux anchor.
* Batalha, N. E. et al. 2019, ApJ 878, 70; Mukherjee, S. et al. 2023, ApJ 942, 71: picaso, the
  linearised flux Newton of the ``ptc`` scheme and the grow-only convective-boundary ratchet.
* Nicholls, H. et al. 2024, MNRAS: AGNI, the mixing-length formulation and the velocity-scaled
  :math:`K_{zz}`.
* Kitzmann, D. et al. 2020, ApJ 890, 174: BeAR (Helios-r2), the retrieval code from which the
  opacity machinery and the Guillot profile descend.
* Tsai, S.-M. et al. 2017, ApJS 228, 20; Tsai, S.-M. et al. 2021, ApJ 923, 264: VULCAN, the
  kinetics code re-implemented as neoVULCAN.

Vertical mixing
---------------

* Charnay, B. et al. 2015, ApJ 813, L1: flux-scaled :math:`K_{zz}` from the convective flux.
* Parmentier, V., Showman, A. P. & Lian, Y. 2013, A&A 558, A91: the tracer-transport
  :math:`K_{zz}` power law for HD 209458b.
* Lee, E. K. H. 2024: GCM-fitted :math:`K_{zz}` profiles for irradiated planets.
* Freytag, B. et al. 2010, A&A 513, A19: mixing in the radiative zones of self-luminous objects.
* Miles, B. E. et al. 2020, AJ 160, 63; Mukherjee, S. et al. 2022, ApJ 938, 107: retrieved
  photospheric :math:`K_{zz}` of T dwarfs.
