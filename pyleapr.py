"""
Python translation of the NJOY ``leapr`` module.

This module provides a partial implementation of the LEAPR algorithm
originally written in Fortran.  The goal of this port is to expose
the core routines used to compute the thermal scattering law
\(S(\alpha,\beta)\) in the incoherent and Gaussian approximations.  The
original Fortran code contained over 100 kB of source spread across a
number of subroutines dealing with everything from phonon expansion and
continuous spectra to translational, discrete oscillator and cold
hydrogen/deuterium corrections.  A line‑by‑line translation of the
entire module into Python would be prohibitively long and
error‑prone.  Instead, this file focuses on the core numerical
components required to compute the continuous part of the scattering
law (the so‑called phonon expansion) and provides a clean and
extensible structure for future development.

**Disclaimer:**

The original LEAPR module supports a wide range of features beyond
what is implemented here – including translational contributions
(diffusion and free gas options), convolution with discrete
oscillators, special handling for cold hydrogen/deuterium, and writing
results in the ENDF‑6 format.  Those portions of the code are
non‑trivial and require careful attention to physics and numeric
details.  This port covers the core phonon expansion algorithm and
supports reading of the standard LEAPR input deck.  Where the
original code contains functionality that has not yet been ported, a
``NotImplementedError`` is raised.  Contributors wishing to extend
this module should consult the original Fortran source for guidance.

References:

* R. E. MacFarlane, "The LEAPR Module of NJOY" (various versions).
* NJOY2016 source code, ``leapr.f90``.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

import numpy as np

# -------------------------------------------------------------------
# Supporting utility functions and ENDF‑I/O routines
#
# The following functions and classes provide minimal support for
# message printing, error handling, timing and ENDF‑6 record writing.
# They are adapted from the NJOY ``util.f90``, ``phys.f90`` and
# ``endf.f90`` modules.  Where possible the original semantics are
# preserved.  These utilities are used by the ``export_to_endf``
# routine in this module to produce a more faithful ENDF‑6 output.

import time
from typing import Iterable, IO

# Logical units for standard input/output/error (mirrors mainio)
NSYSI = sys.stdin
NSYSO = sys.stderr  # error log
NSYSE = sys.stdout  # terminal output

def error(from_where: str, mess1: str, mess2: str = "") -> None:
    """Fatal error (mirrors util.f90:error).  Prints a message on
    both the error and standard outputs and aborts execution.

    Parameters
    ----------
    from_where : str
        Name of the calling routine.
    mess1, mess2 : str, optional
        Error messages to print.
    """
    print(f"\n***error in {from_where}*** {mess1}", file=NSYSO)
    if mess2.strip():
        print(" " * 22 + mess2, file=NSYSO)
    print("", file=NSYSO)
    if NSYSE is not NSYSO:
        print(f"\n***error in {from_where}*** {mess1}", file=NSYSE)
        if mess2.strip():
            print(" " * 22 + mess2, file=NSYSE)
        print("", file=NSYSE)
    raise SystemExit(77)

def mess(from_where: str, mess1: str, mess2: str = "") -> None:
    """Non‑fatal message (mirrors util.f90:mess).

    Parameters
    ----------
    from_where : str
        Name of the calling routine.
    mess1, mess2 : str, optional
        Messages to print.
    """
    print(f"\n---message from {from_where}--- {mess1}", file=NSYSO)
    if mess2.strip():
        print(" " * 26 + mess2, file=NSYSO)
    if NSYSE is not NSYSO:
        print(f"\n---message from {from_where}--- {mess1}", file=NSYSE)
        if mess2.strip():
            print(" " * 26 + mess2, file=NSYSE)

def wclock() -> str:
    """Return a wall clock time string in HH:MM:SS format."""
    return time.strftime("%H:%M:%S", time.localtime())

# Additional physical constants (subset of phys.f90).  Some constants
# such as BK, EV, HBAR, AMASSN, AMU and ANGST are already defined in
# the Leapr class below.  Here we define a few extras for completeness.
euler = 0.57721566490153286
clight = 2.99792458e10  # speed of light (cm/s)
finstri = 1.0e16 * (6.582119569e-16 * 1.602176634e-12) / (1.602176634e-12 ** 2 * clight)

# -------------------------------------------------------------------
# ENDF‑6 I/O helper routines

class EndfState:
    """Global ENDF state for record writing.

    This class mirrors the common blocks in endf.f90 that store
    current material (MAT), file (MF), and section (MT) identifiers
    and track the running line count (NS).  These variables are
    implicitly updated by the ENDF I/O routines defined below.
    """
    def __init__(self) -> None:
        self.math: int = 0
        self.mfh: int = 0
        self.mth: int = 0
        self.nsh: int = 1  # line number counter (1..99999)

# Instantiate the global state
STATE = EndfState()

def a11(x: float) -> str:
    """Convert a floating point number to ENDF 11‑column format.

    This function attempts to reproduce the Fortran formatting logic
    implemented in endf.f90 for scientific and fixed notation.  It
    outputs a string of length 11, including sign and exponent.

    Parameters
    ----------
    x : float
        Number to format.

    Returns
    -------
    str
        The formatted 11‑character string.
    """
    # handle zero specially
    if x == 0.0:
        return "      0.0+0"
    # Convert to scientific notation (e.g., '1.23456e+03')
    s = f"{x:.5e}"
    mant_str, exp_str = s.split('e')
    mant = float(mant_str)
    exp = int(exp_str)
    # Format mantissa with explicit sign and 6 digits after decimal
    mant_form = f"{mant:+.6f}"
    mant_form = mant_form.replace(' ', '')
    exp_form = f"{exp:+d}"
    out = f"{mant_form[:-2]}{exp_form}"
    return out.rjust(11)

def _line(mat: int, mf: int, mt: int, ns: int, six_fields: Iterable[str]) -> str:
    """Assemble a single ENDF record line from six fields and a trailer."""
    body = "".join(f"{f:>11}" for f in six_fields)
    trailer = f"{mat:>4}{mf:>2}{mt:>3}{ns:>5}"
    return body + trailer

def asend(nout: IO[str]) -> None:
    """Write a section end (ASEND) record."""
    global STATE
    line = _line(STATE.math, STATE.mfh, 0, STATE.nsh, ["", "", 0, 0, 0, 0])
    print(line, file=nout)
    STATE.nsh = 1

def afend(nout: IO[str]) -> None:
    """Write a file end (AFEND) record."""
    global STATE
    line = _line(STATE.math, 0, 0, STATE.nsh, ["", "", 0, 0, 0, 0])
    print(line, file=nout)
    STATE.nsh = 1

def amend(nout: IO[str]) -> None:
    """Write a material end (AMEND) record."""
    global STATE
    line = _line(0, 0, 0, STATE.nsh, ["", "", 0, 0, 0, 0])
    print(line, file=nout)
    STATE.nsh = 1

def atend(nout: IO[str]) -> None:
    """Write a tape end (ATEND) record."""
    global STATE
    line = _line(-1, 0, 0, STATE.nsh, ["", "", 0, 0, 0, 0])
    print(line, file=nout)
    STATE.nsh = 1

def contio(nout: IO[str], c1: float, c2: float, l1: int, l2: int, n1: int, n2: int) -> None:
    """Write a CONT record using the global STATE variables."""
    global STATE
    fields = [a11(c1), a11(c2), f"{l1:>11}", f"{l2:>11}", f"{n1:>11}", f"{n2:>11}"]
    print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, fields), file=nout)
    STATE.nsh += 1
    if STATE.nsh > 99999:
        STATE.nsh = 1

def listio(nout: IO[str], c1: float, c2: float, l1: int, l2: int,
           npl: int, n2: int, list_vals: Iterable[float]) -> None:
    """Write a LIST record with header and body."""
    global STATE
    fields = [a11(c1), a11(c2), f"{l1:>11}", f"{l2:>11}", f"{npl:>11}", f"{n2:>11}"]
    print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, fields), file=nout)
    STATE.nsh += 1
    # Now write the list body 6 per line
    row: List[str] = []
    for v in list_vals:
        row.append(a11(float(v)))
        if len(row) == 6:
            print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, row), file=nout)
            STATE.nsh += 1
            row = []
    if row:
        while len(row) < 6:
            row.append("")
        print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, row), file=nout)
        STATE.nsh += 1

def tab1io(nout: IO[str], c1: float, c2: float, l1: int, l2: int,
           nbt: List[int], intt: List[int], x: List[float], y: List[float]) -> None:
    """Write a TAB1 record: header with interpolation table and x,y pairs."""
    global STATE
    if len(nbt) != len(intt):
        raise ValueError("nbt and intt lengths must match")
    nr = len(nbt)
    npairs = len(x)
    fields = [a11(c1), a11(c2), f"{l1:>11}", f"{l2:>11}", f"{nr:>11}", f"{npairs:>11}"]
    print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, fields), file=nout)
    STATE.nsh += 1
    # Interpolation table pairs: (NBT,INT) three per line
    row: List[str] = []
    for i in range(nr):
        row.append(f"{nbt[i]:>11}")
        row.append(f"{intt[i]:>11}")
        if len(row) == 6 or (i == nr - 1):
            while len(row) < 6:
                row.append("")
            print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, row), file=nout)
            STATE.nsh += 1
            row = []
    # Data pairs: x,y (3 pairs per line)
    row = []
    for xi, yi in zip(x, y):
        row.append(a11(float(xi)))
        row.append(a11(float(yi)))
        if len(row) == 6:
            print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, row), file=nout)
            STATE.nsh += 1
            row = []
    if row:
        while len(row) < 6:
            row.append("")
        print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, row), file=nout)
        STATE.nsh += 1

def tab2io(nout: IO[str], c1: float, c2: float, l1: int, l2: int,
           nbt: List[int], intt: List[int], n2: int) -> None:
    """Write a TAB2 header (no data pairs)."""
    global STATE
    if len(nbt) != len(intt):
        raise ValueError("nbt and intt lengths must match")
    nr = len(nbt)
    fields = [a11(c1), a11(c2), f"{l1:>11}", f"{l2:>11}", f"{nr:>11}", f"{n2:>11}"]
    print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, fields), file=nout)
    STATE.nsh += 1
    row: List[str] = []
    for i in range(nr):
        row.append(f"{nbt[i]:>11}")
        row.append(f"{intt[i]:>11}")
        if len(row) == 6 or (i == nr - 1):
            while len(row) < 6:
                row.append("")
            print(_line(STATE.math, STATE.mfh, STATE.mth, STATE.nsh, row), file=nout)
            STATE.nsh += 1
            row = []

def set_mat(mat: int) -> None:
    """Set the MAT identifier for subsequent ENDF records."""
    STATE.math = mat

def set_mf(mf: int) -> None:
    """Set the MF identifier for subsequent ENDF records."""
    STATE.mfh = mf

def set_mt(mt: int) -> None:
    """Set the MT identifier for subsequent ENDF records."""
    STATE.mth = mt

def reset_ns() -> None:
    """Reset the line counter (NS) to 1."""
    STATE.nsh = 1

# -------------------------------------------------------------------
# Additional utility: significant figure rounding (from util.f90)

def sigfig(x: float, ndig: int, idig: int) -> float:
    """
    Round a floating point number to ``ndig`` significant figures.
    Optionally adjust the last digit by ``idig``.  This function
    mirrors the behaviour of the Fortran ``sigfig`` function in
    util.f90.  A small bias is applied to avoid boundary effects.

    Parameters
    ----------
    x : float
        Value to round.
    ndig : int
        Number of significant digits to retain.
    idig : int
        Adjustment applied to the least significant digit (typically 0).

    Returns
    -------
    float
        The value of ``x`` rounded to ``ndig`` significant figures.
    """
    if x == 0.0:
        return 0.0
    # Compute approximate exponent of x in base 10
    aa = math.log10(abs(x))
    ipwr = int(aa)
    if aa < 0.0:
        ipwr -= 1
    ipwr = ndig - 1 - ipwr
    # Form integer representation at the desired precision
    ii = int(round(x * (10.0 ** ipwr) + (10.0 ** (ndig - 11))))
    if ii >= 10 ** ndig:
        ii //= 10
        ipwr -= 1
    # Apply digit adjustment
    ii += idig
    xx = ii * (10.0 ** (-ipwr))
    # Apply slight bias as in original Fortran
    return xx * 1.0000000000001


@dataclass
class LeaprInput:
    """Container for user input parameters.

    The LEAPR module is driven by a structured input deck.  For
    convenience this dataclass stores the quantities read from the
    input file.  See the comments in the original Fortran for
    detailed descriptions of each field.
    """

    # card 1
    nout: int  # ENDF output unit (unused in this port)

    # card 3 – run control
    ntempr: int = 1
    iprint: int = 1
    nphon: int = 100

    # card 4 – ENDF output control
    mat: int = 0
    za: float = 0.0
    isabt: int = 0
    ilog: int = 0
    smin: float = 1.0e-75

    # card 5 – principal scatterer control
    awr: float = 0.0
    spr: float = 0.0
    npr: int = 1
    iel: int = 0
    ncold: int = 0
    nsk: int = 0

    # card 6 – secondary scatterer control
    nss: int = 0
    b7: float = 0.0
    aws: float = 0.0
    sps: float = 0.0
    mss: int = 0

    # card 7 – α/β control
    nalpha: int = 0
    nbeta: int = 0
    lat: int = 0

    # card 8/9 – α and β values
    alpha: List[float] = field(default_factory=list)
    beta: List[float] = field(default_factory=list)

    # per‑temperature parameters (cards 10–18) will be stored
    temperatures: List[float] = field(default_factory=list)
    continuous: List[Tuple[float, List[float], float, float, float]] = field(default_factory=list)
    # each entry is (delta1, p1 array, twt, c, tbeta)
    oscillators: List[Tuple[List[float], List[float]]] = field(default_factory=list)
    # each entry is (bdel, adel) for the discrete oscillator section
    pair_corr: List[Tuple[int, float, List[float]]] = field(default_factory=list)
    # each entry is (nka, dka, ska array) if present
    cfracs: List[float] = field(default_factory=list)  # only used for nsk>0


class Leapr:
    """
    A Python implementation of the core LEAPR algorithm.

    Instances of this class are intended to be configured with a
    :class:`LeaprInput` object.  The ``run`` method performs the
    calculation of S(α,β) on the provided α and β grids for each
    temperature.  Results are stored internally and can be queried
    through the ``ssm`` attribute.
    """

    def __init__(self, leapr_input: LeaprInput) -> None:
        self.inp = leapr_input
        # allocate arrays for results: shape (ntempr, nbeta, nalpha)
        ntemp = self.inp.ntempr
        nbeta = self.inp.nbeta
        nalpha = self.inp.nalpha
        self.ssm = np.zeros((ntemp, nbeta, nalpha), dtype=float)
        # placeholders for Debye–Waller factor and effective temperature
        self.dwpix = np.zeros(ntemp, dtype=float)
        self.tempf = np.zeros(ntemp, dtype=float)
        # additional state used during computations
        self.f0 = 0.0
        self.tbar = 0.0
        self.tev = 0.0
        self.deltab = 0.0
        # allocate a second scattering array for positive β values
        # For most contributions (phonon expansion, translational and
        # discrete oscillators), the scattering law is symmetric in β
        # and ssp will simply mirror ssm.  When the cold hydrogen/
        # deuterium option is enabled, the scattering law becomes
        # asymmetric; ssm will hold the negative‑β branch and ssp the
        # positive‑β branch.
        self.ssp = np.zeros((ntemp, nbeta, nalpha), dtype=float)

        # placeholders for principal scatterer Debye–Waller and effective
        # temperature values when secondary mixing is applied.  These
        # arrays will be created lazily in run() if needed.
        self.dwp1: Optional[np.ndarray] = None
        self.tempf1: Optional[np.ndarray] = None

    # ------------------------------------------------------------------
    # Helper routines for coherent elastic scattering (Bragg edges)
    #
    # The NJOY LEAPR module includes a ``coher`` subroutine that
    # calculates Bragg edge energies and associated structure factors
    # for a handful of crystalline moderators (graphite, beryllium,
    # beryllium oxide, aluminum, lead and iron).  The Fortran code is
    # quite involved: it sets up lattice constants, loops over
    # reciprocal lattice vectors, applies form factors and Debye–
    # Waller damping and finally sorts and combines duplicate edges.
    # The routines below provide a direct translation of the core
    # calculations used by ``coher``.  They are not called
    # automatically by ``run`` but can be invoked by client code if
    # coherent scattering data are required.

    def _formf(self, lat: int, l1: int, l2: int, l3: int) -> float:
        """Compute form factors for the specified lattice.

        This function replicates the Fortran ``formf`` routine.  It
        returns the squared structure factor for the lattice defined by
        ``lat`` and Miller indices ``(l1,l2,l3)``.  The meaning of
        ``lat`` is as follows:

        * ``1`` – graphite (hexagonal)
        * ``2`` – beryllium (hexagonal)
        * ``3`` – beryllium oxide (hexagonal)
        * ``4`` – fcc lattice (aluminum)
        * ``5`` – fcc lattice (lead)
        * ``6`` – bcc lattice (iron)

        Parameters
        ----------
        lat : int
            Lattice type.
        l1, l2, l3 : int
            Integer Miller indices.

        Returns
        -------
        float
            The form factor (squared structure factor).
        """
        pi = self.PI
        # graphite
        if lat == 1:
            # graphite has an ABAB stacking that gives rise to
            # selection rules depending on the c‑axis index l3.  If
            # l3 is odd, the form factor is given by sin^2(pi*(l1-l2)/3),
            # otherwise by (6+10*cos(2*pi*(l1-l2)/3))/4.
            i = l3 // 2
            if 2 * i != l3:
                return math.sin(pi * (l1 - l2) / 3.0) ** 2
            else:
                return (6.0 + 10.0 * math.cos(2.0 * pi * (l1 - l2) / 3.0)) / 4.0
        # beryllium
        elif lat == 2:
            # For beryllium the two atoms per primitive cell lead to a
            # cosine modulation.
            return 1.0 + math.cos(2.0 * pi * (2 * l1 + 4 * l2 + 3 * l3) / 6.0)
        # beryllium oxide
        elif lat == 3:
            # Beryllium oxide has two different atomic species in the
            # hexagonal cell, giving rise to additional coefficients.
            c1 = 7.54
            c2 = 4.24
            c3 = 11.31
            return (
                1.0
                + math.cos(2.0 * pi * (2 * l1 + 4 * l2 + 3 * l3) / 6.0)
            ) * (c1 + c2 + c3 * math.cos(3.0 * pi * l3 / 4.0))
        # fcc lattices (aluminum, lead)
        elif lat == 4 or lat == 5:
            e1 = 2.0 * pi * l1
            e2 = 2.0 * pi * (l1 + l2)
            e3 = 2.0 * pi * (l1 + l3)
            # The form factor is the magnitude squared of the sum of
            # contributions from the basis vectors in an fcc cell.
            return (
                (1.0 + math.cos(e1) + math.cos(e2) + math.cos(e3)) ** 2
                + (math.sin(e1) + math.sin(e2) + math.sin(e3)) ** 2
            )
        # bcc lattice (iron)
        elif lat == 6:
            e1 = 2.0 * pi * (l1 + l2 + l3)
            return (1.0 + math.cos(e1)) ** 2 + (math.sin(e1)) ** 2
        # unsupported lattice
        else:
            raise ValueError(f"Illegal lattice type {lat} in formf")

    @staticmethod
    def _tausq(m1: int, m2: int, m3: int, c1: float, c2: float, twopis: float) -> float:
        """Reciprocal lattice magnitude squared for hexagonal structures.

        This function matches the Fortran function ``tausq``.  It
        computes

        .. math::

           \tau^2 = [c_1 (m_1^2 + m_2^2 + m_1 m_2) + m_3^2 c_2] \times 2\pi^2,

        where ``c1`` and ``c2`` are geometry constants and ``twopis`` is
        \((2\pi)^2\).  The indices ``m1``, ``m2`` and ``m3`` are the
        integer lattice indices along the reciprocal primitive vectors.

        Parameters
        ----------
        m1, m2, m3 : int
            Lattice indices.
        c1, c2 : float
            Geometry constants.
        twopis : float
            Precomputed constant ``(2*pi)**2``.

        Returns
        -------
        float
            The squared reciprocal lattice magnitude.
        """
        return (c1 * (m1 * m1 + m2 * m2 + m1 * m2) + m3 * m3 * c2) * twopis

    @staticmethod
    def _taufcc(m1: int, m2: int, m3: int, c1: float, twothd: float, twopis: float) -> float:
        """Reciprocal lattice magnitude squared for fcc lattices.

        This matches the Fortran function ``taufcc``, computing

        .. math::

           \tau^2 = c_1 (m_1^2 + m_2^2 + m_3^2 + \tfrac{2}{3}m_1 m_2 + \tfrac{2}{3}m_1 m_3 - \tfrac{2}{3}m_2 m_3)\times 2\pi^2.

        The parameter ``twothd`` is 2/3.
        """
        return (
            c1
            * (
                m1 * m1
                + m2 * m2
                + m3 * m3
                + twothd * m1 * m2
                + twothd * m1 * m3
                - twothd * m2 * m3
            )
            * twopis
        )

    @staticmethod
    def _taubcc(m1: int, m2: int, m3: int, c1: float, twopis: float) -> float:
        """Reciprocal lattice magnitude squared for bcc lattices.

        Matches the Fortran function ``taubcc``.  Note that in the
        original Fortran, ``c1`` is a host variable rather than a
        formal parameter; in this translation we pass ``c1``
        explicitly.
        """
        return (
            c1
            * (m1 * m1 + m2 * m2 + m3 * m3 + m1 * m2 + m2 * m3 + m1 * m3)
            * twopis
        )

    def coher(self, lat: int, natom: int, emax: float = 5.0) -> Tuple[List[float], List[float]]:
        """Compute Bragg edge energies and structure factors.

        This method is a direct translation of the LEAPR ``coher``
        subroutine.  Given a lattice type ``lat`` (1–6) and the
        number of atoms per unit cell ``natom``, it computes all
        reciprocal lattice vectors with squared magnitudes below
        ``emax`` (in eV) and returns their energies and structure
        factors.  The returned lists are sorted in ascending energy
        with duplicate edges combined.

        Parameters
        ----------
        lat : int
            Lattice type (1–6).  See :func:`_formf` for definitions.
        natom : int
            Number of atoms per unit cell.
        emax : float, optional
            Maximum energy in eV below which to keep Bragg edges
            (default 5 eV).

        Returns
        -------
        Tuple[List[float], List[float]]
            Two lists of equal length.  The first contains Bragg
            energies (in eV), the second contains the corresponding
            coherent elastic structure factors (in barns).  Duplicate
            energies are merged and their structure factors summed.
        """
        # Physical constants (matching the Fortran ``physics`` module)
        pi = self.PI
        # lattice constants for different materials (in cm and amu)
        gr1 = 2.4573e-8
        gr2 = 6.700e-8
        gr3 = 12.011
        gr4 = 5.50
        be1 = 2.2856e-8
        be2 = 3.5832e-8
        be3 = 9.01
        be4 = 7.53
        beo1 = 2.695e-8
        beo2 = 4.39e-8
        beo3 = 12.5
        beo4 = 1.0
        al1 = 4.04e-8
        al3 = 26.7495
        al4 = 1.495
        pb1 = 4.94e-8
        pb3 = 207.0
        pb4 = 1.0
        fe1 = 2.86e-8
        fe3 = 55.454
        fe4 = 12.9
        twothd = 0.666666666667
        sqrt3 = 1.732050808
        toler = 1.0e-6
        eps = 0.05
        # Precompute physical constants used in the conversion
        twopis = (2.0 * pi) ** 2
        # neutron mass in grams (amu * c^2 converts energy units to g)
        amne = self.AMASSN * self.AMU
        # econ converts tau^2 to energy (eV).  This matches
        # econ = ev*8*(amne/hbar)/hbar from the Fortran code
        econ = self.EV * 8.0 * (amne / self.HBAR) / self.HBAR
        recon = 1.0 / econ
        tsqx = econ / 20.0
        # Select lattice parameters based on ``lat``
        if lat == 1:
            # graphite
            a = gr1
            c = gr2
            amsc = gr3
            scoh = gr4 / natom
        elif lat == 2:
            # beryllium
            a = be1
            c = be2
            amsc = be3
            scoh = be4 / natom
        elif lat == 3:
            # beryllium oxide
            a = beo1
            c = beo2
            amsc = beo3
            scoh = beo4 / natom
        elif lat == 4:
            # aluminum (fcc)
            a = al1
            amsc = al3
            scoh = al4 / natom
            c = None  # not used for fcc
        elif lat == 5:
            # lead (fcc)
            a = pb1
            amsc = pb3
            scoh = pb4 / natom
            c = None
        elif lat == 6:
            # iron (bcc)
            a = fe1
            amsc = fe3
            scoh = fe4 / natom
            c = None
        else:
            raise ValueError(f"Illegal lattice type {lat} in coher")
        # Geometry constants and scaling factor ``scon``
        if lat < 4:
            c1 = 4.0 / (3.0 * a * a)
            c2 = 1.0 / (c * c)  # type: ignore[arg-type]
            scon = scoh * (4.0 * pi) ** 2 / (2.0 * a * a * c * sqrt3 * econ)  # type: ignore[arg-type]
        elif lat in (4, 5):
            c1 = 3.0 / (a * a)
            c2 = 0.0
            scon = scoh * (4.0 * pi) ** 2 / (16.0 * a * a * a * econ)
        else:  # lat == 6
            c1 = 2.0 / (a * a)
            c2 = 0.0
            scon = scoh * (4.0 * pi) ** 2 / (8.0 * a * a * a * econ)
        # Debye–Waller prefactors (wint is zero in LEAPR)
        wint = 0.0
        t2 = self.HBAR / (2.0 * self.AMU * amsc)
        # Limit on tau^2 values converted from ``emax``
        ulim = econ * emax
        # Storage for squared magnitudes and weights
        tsq_list: List[float] = []
        f_list: List[float] = []
        # Compute reciprocal lattice vectors for hexagonal structures
        if lat <= 3:
            phi = ulim / twopis
            # Maximum index along the a‑axis
            i1m = int(a * math.sqrt(phi)) + 1
            for i1 in range(1, i1m + 1):
                l1 = i1 - 1
                # Maximum index along the second a‑axis (depends on l1)
                tmp = 3.0 * (a * a * phi - l1 * l1)
                if tmp < 0.0:
                    continue
                i2m = int((l1 + math.sqrt(tmp)) / 2.0) + 1
                for i2 in range(i1, i2m + 1):
                    l2 = i2 - 1
                    # Remaining quadratic term inside the square root for l3
                    x = phi - c1 * (l1 * l1 + l2 * l2 - l1 * l2)
                    i3m = 0
                    if x > 0.0:
                        i3m = int(c * math.sqrt(x)) + 1  # type: ignore[arg-type]
                    else:
                        i3m = 1
                    for i3 in range(1, i3m + 1):
                        l3 = i3 - 1
                        # Weight factors due to symmetry
                        w1 = 2.0
                        if l1 == l2:
                            w1 = 1.0
                        w2 = 2.0
                        if (l1 == 0) or (l2 == 0):
                            w2 = 1.0
                        if (l1 == 0) and (l2 == 0):
                            w2 = w2 / 2.0
                        w3 = 2.0
                        if l3 == 0:
                            w3 = 1.0
                        # Evaluate tau^2 for (l1,l2,l3)
                        tsq = self._tausq(l1, l2, l3, c1, c2, twopis)
                        if tsq > 0.0 and tsq <= ulim:
                            tau = math.sqrt(tsq)
                            w = math.exp(-tsq * t2 * wint) * w1 * w2 * w3 / tau
                            f = w * self._formf(lat, l1, l2, l3)
                            # Store or merge duplicate entries
                            if not tsq_list or tsq <= tsqx:
                                tsq_list.append(tsq)
                                f_list.append(f)
                            else:
                                merged = False
                                for idx in range(len(tsq_list)):
                                    tsq_old = tsq_list[idx]
                                    if tsq >= tsq_old and tsq < (1.0 + eps) * tsq_old:
                                        f_list[idx] += f
                                        merged = True
                                        break
                                if not merged:
                                    tsq_list.append(tsq)
                                    f_list.append(f)
                        # Evaluate tau^2 for (l1,-l2,l3) to capture
                        # additional symmetry
                        tsq2 = self._tausq(l1, -l2, l3, c1, c2, twopis)
                        if tsq2 > 0.0 and tsq2 <= ulim:
                            tau = math.sqrt(tsq2)
                            w = math.exp(-tsq2 * t2 * wint) * w1 * w2 * w3 / tau
                            f = w * self._formf(lat, l1, -l2, l3)
                            if not tsq_list or tsq2 <= tsqx:
                                tsq_list.append(tsq2)
                                f_list.append(f)
                            else:
                                merged = False
                                for idx in range(len(tsq_list)):
                                    tsq_old = tsq_list[idx]
                                    if tsq2 >= tsq_old and tsq2 < (1.0 + eps) * tsq_old:
                                        f_list[idx] += f
                                        merged = True
                                        break
                                if not merged:
                                    tsq_list.append(tsq2)
                                    f_list.append(f)
        else:
            # Reciprocal lattice vectors for cubic structures
            # For fcc and bcc the Fortran code artificially truncates the
            # search region to a 15×15×15 cube around the origin.
            i1m = 15
            if lat in (4, 5):
                # fcc lattices
                for i1 in range(-i1m, i1m + 1):
                    for i2 in range(-i1m, i1m + 1):
                        for i3 in range(-i1m, i1m + 1):
                            tsq = self._taufcc(i1, i2, i3, c1, twothd, twopis)
                            if tsq > 0.0 and tsq <= ulim:
                                tau = math.sqrt(tsq)
                                w = math.exp(-tsq * t2 * wint) / tau
                                f = w * self._formf(lat, i1, i2, i3)
                                tsq_list.append(tsq)
                                f_list.append(f)
            else:
                # bcc lattice (lat == 6)
                for i1 in range(-i1m, i1m + 1):
                    for i2 in range(-i1m, i1m + 1):
                        for i3 in range(-i1m, i1m + 1):
                            tsq = self._taubcc(i1, i2, i3, c1, twopis)
                            if tsq > 0.0 and tsq <= ulim:
                                tau = math.sqrt(tsq)
                                w = math.exp(-tsq * t2 * wint) / tau
                                f = w * self._formf(lat, i1, i2, i3)
                                tsq_list.append(tsq)
                                f_list.append(f)
        # Sort the collected tau^2 values and weights
        pairs = sorted(zip(tsq_list, f_list), key=lambda p: p[0])
        # Merge duplicates and convert to energies and structure factors
        energies: List[float] = []
        sigmas: List[float] = []
        bel = -1.0
        for tsq, f in pairs:
            be = tsq * recon
            bs = f * scon
            if bel < 0.0:
                # first entry
                energies.append(be)
                sigmas.append(bs)
                bel = be
            else:
                if (be - bel) < toler:
                    sigmas[-1] += bs
                else:
                    energies.append(be)
                    sigmas.append(bs)
                    bel = be
        return energies, sigmas

    # ------------------------------------------------------------------
    # Bragg edge thinning
    #
    # When coherent elastic scattering is requested (iel > 0), the ENDF
    # writer in NJOY computes cumulative sums over Bragg edges with
    # Debye–Waller damping and drops edges whose incremental
    # contribution falls below a relative tolerance.  This procedure
    # reduces the number of edges recorded in the final file.  The
    # method below exposes the same logic: given the internally
    # computed Debye–Waller factors and the coherent Bragg edge data
    # from :func:`coher`, it returns a thinned set of edges along with
    # cumulative cross sections for each temperature.

    def bragg_edge_thinning(self, emax: float = 5.0, tol: float = 0.9e-7) -> Tuple[List[float], List[List[float]]]:
        """Thin Bragg edges using the Debye–Waller criterion.

        If the input specifies a coherent scatterer (``iel > 0``), the
        original LEAPR code thins the list of Bragg edges so that
        successive edges contribute more than a relative tolerance to
        the cumulative elastic scattering.  The Debye–Waller factor for
        each temperature controls the damping.  This method computes
        the thinned list of edge energies and, for each temperature,
        the cumulative (integrated) coherent elastic scattering up to
        each retained edge.

        Parameters
        ----------
        emax : float, optional
            Maximum energy (in eV) used when computing Bragg edges.
            Defaults to 5 eV, which matches the Fortran default in
            ``leapr``.
        tol : float, optional
            Relative tolerance used to decide whether a new Bragg edge
            contributes significantly to the sum.  The default matches
            the value used in the Fortran code (``0.9e-7``).

        Returns
        -------
        Tuple[List[float], List[List[float]]]
            The first element is a list of thinned Bragg edge energies
            (in eV).  The second element is a list of lists; each
            inner list contains the cumulative coherent elastic
            scattering values (in barns) for a corresponding
            temperature, evaluated at each thinned edge.  The number
            of inner lists equals the number of temperatures in the
            calculation.

        Notes
        -----
        This method uses ``self.dwpix`` and ``self.dwp1`` to compute
        effective Debye–Waller factors.  The factors are scaled by
        ``A*T*k_B`` for consistency with the Fortran code.  If a
        secondary moderator is present (``nss > 0`` and ``b7 <= 0``),
        the Debye–Waller factors of the primary and secondary
        scatterers are averaged.  The method silently returns empty
        lists if no coherent scattering is requested (``iel <= 0``).
        """
        # Only proceed if coherent scattering was requested
        if self.inp.iel <= 0:
            return [], []
        # Compute Bragg edges (energies and structure factors)
        energies, factors = self.coher(self.inp.iel, self.inp.npr, emax)
        nedge = len(energies)
        if nedge == 0:
            return [], []
        # Determine the number of edges to keep using the tolerance
        # criterion.  Use the Debye–Waller factor at the first
        # temperature (or averaged if secondary scatterer present) to
        # compute the weighting.  In the Fortran code, w =
        # dwpix(1)/(A*T*k_B); for mixed moderators, w is the
        # arithmetic mean of dwpix and dwp1.  Here ``awr`` and
        # ``tempr`` come from the input deck.
        # Compute the scaling factor A*T*k_B for the first temperature
        t0 = self.inp.temperatures[0] if self.inp.temperatures else 0.0
        if t0 <= 0.0:
            # guard against zero temperature
            return [], []
        # Primary scatterer mass
        awr = self.inp.awr
        w0 = self.dwpix[0] / (awr * t0 * self.BK)
        # If secondary mixing is present, average Debye–Waller
        if self.inp.nss > 0 and self.inp.b7 <= 0.0 and self.dwp1 is not None:
            w0 = (self.dwpix[0] + self.dwp1[0]) / 2.0 / (awr * t0 * self.BK)
        # Compute jmax (last index to retain) based on tolerance
        sum_val = 0.0
        suml = 0.0
        jmax = nedge - 1
        for j in range(nedge):
            e = energies[j]
            sum_val += math.exp(-4.0 * w0 * e) * factors[j]
            if sum_val - suml > tol * sum_val:
                jmax = j
                suml = sum_val
        # Truncate the edges and factors to jmax+1
        energies_thin = energies[: jmax + 1]
        factors_thin = factors[: jmax + 1]
        # Now compute cumulative sums for each temperature
        ntemp = self.inp.ntempr
        cumulative_lists: List[List[float]] = []
        for it in range(ntemp):
            # Temperature and scattering mass factor
            T = self.inp.temperatures[it] if self.inp.temperatures else 0.0
            if T <= 0.0:
                cumulative_lists.append([0.0 for _ in energies_thin])
                continue
            w = self.dwpix[it] / (awr * T * self.BK)
            if self.inp.nss > 0 and self.inp.b7 <= 0.0 and self.dwp1 is not None:
                w = (self.dwpix[it] + self.dwp1[it]) / 2.0 / (awr * T * self.BK)
            cum_sum = 0.0
            vals: List[float] = []
            for e, f in zip(energies_thin, factors_thin):
                cum_sum += math.exp(-4.0 * w * e) * f
                vals.append(cum_sum)
            cumulative_lists.append(vals)
        return energies_thin, cumulative_lists

    # Physical constants used in multiple routines
    BK = 8.617333262e-5  # Boltzmann constant in eV/K (same units as Fortran)
    PI = math.pi

    # Additional physical constants used in pair‑correlation and cold‑hydrogen routines.
    # These are taken from the Fortran ``physics`` module.  They are
    # required for computing wavenumbers in the Sköld approximation.
    EV = 1.602176634e-12  # erg per electron‑volt (erg/eV)
    HBAR = 6.582119569e-16 * EV  # Planck constant (ℏ) in erg·s (converted from eV·s)
    AMASSN = 1.00866491595  # neutron mass in atomic mass units (amu)
    AMU = 931.49410242e6 * EV / (2.99792458e10 ** 2)  # atomic mass unit in grams (converted)
    ANGST = 1.0e-8  # Angstrom in centimetres


    def read_input(self, fname: str) -> None:
        """
        Read a LEAPR input deck from ``fname`` and populate the
        :class:`LeaprInput` object.  This method is provided for
        compatibility with legacy LEAPR workflows.  It follows the
        structure described in the comments of the original Fortran
        subroutine ``leapr``.  The implementation here is minimal
        and does not attempt to interpret optional card content.

        Parameters
        ----------
        fname : str
            Path to the LEAPR input file.
        """
        with open(fname, 'r') as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        # A very simple parser: assumes that the cards appear in order
        # and that numeric values are whitespace separated.  Real world
        # LEAPR decks are more flexible; you may need to adapt this
        # parser to your specific input format.
        idx = 0
        # card 1
        self.inp.nout = int(lines[idx].split()[0])
        idx += 1
        # card 2 – title (ignored in this port)
        idx += 1
        # card 3
        vals = list(map(float, lines[idx].split()))
        self.inp.ntempr = int(vals[0])
        self.inp.iprint = int(vals[1])
        self.inp.nphon = int(vals[2])
        idx += 1
        # card 4
        vals = list(map(float, lines[idx].split()))
        self.inp.mat = int(vals[0])
        self.inp.za = vals[1]
        self.inp.isabt = int(vals[2])
        self.inp.ilog = int(vals[3])
        self.inp.smin = vals[4]
        idx += 1
        # card 5
        vals = list(map(float, lines[idx].split()))
        self.inp.awr = vals[0]
        self.inp.spr = vals[1]
        self.inp.npr = int(vals[2])
        self.inp.iel = int(vals[3])
        self.inp.ncold = int(vals[4])
        self.inp.nsk = int(vals[5])
        idx += 1
        # card 6
        vals = list(map(float, lines[idx].split()))
        self.inp.nss = int(vals[0])
        self.inp.b7 = vals[1]
        self.inp.aws = vals[2]
        self.inp.sps = vals[3]
        self.inp.mss = int(vals[4])
        idx += 1
        # card 7
        vals = list(map(float, lines[idx].split()))
        self.inp.nalpha = int(vals[0])
        self.inp.nbeta = int(vals[1])
        self.inp.lat = int(vals[2])
        idx += 1
        # card 8 – α values
        self.inp.alpha = list(map(float, lines[idx].split()))
        idx += 1
        # card 9 – β values
        self.inp.beta = list(map(float, lines[idx].split()))
        idx += 1
        # Cards 10–18: repeated per temperature
        for t in range(self.inp.ntempr):
            # temperature
            temp_val = float(lines[idx].split()[0])
            self.inp.temperatures.append(temp_val)
            idx += 1
            # card 11 – continuous distribution control
            vals = list(map(float, lines[idx].split()))
            delta1 = vals[0]
            ni = int(vals[1])
            idx += 1
            # card 12 – p1
            p1_vals = list(map(float, lines[idx].split()))
            assert len(p1_vals) == ni, "Mismatch in number of p1 points"
            idx += 1
            # card 13 – continuous distribution parameters
            vals = list(map(float, lines[idx].split()))
            twt, c, tbeta = vals[:3]
            idx += 1
            self.inp.continuous.append((delta1, p1_vals, twt, c, tbeta))
            # card 14 – discrete oscillator control
            vals = list(map(float, lines[idx].split()))
            nd = int(vals[0])
            idx += 1
            if nd > 0:
                # card 15 – oscillator energies
                bdel_vals = list(map(float, lines[idx].split()))
                assert len(bdel_vals) == nd
                idx += 1
                # card 16 – oscillator weights
                adel_vals = list(map(float, lines[idx].split()))
                assert len(adel_vals) == nd
                idx += 1
                self.inp.oscillators.append((bdel_vals, adel_vals))
            else:
                self.inp.oscillators.append(([], []))
            # pair correlation control (card 17/18)
            if self.inp.nsk > 0 or self.inp.ncold > 0:
                vals = list(map(float, lines[idx].split()))
                nka = int(vals[0])
                dka = vals[1]
                idx += 1
                ska_vals = list(map(float, lines[idx].split()))
                assert len(ska_vals) == nka
                idx += 1
                self.inp.pair_corr.append((nka, dka, ska_vals))
                # cfrac (card 19) for skold method
                if self.inp.nsk > 0:
                    cfrac_val = float(lines[idx].split()[0])
                    self.inp.cfracs.append(cfrac_val)
                    idx += 1
            else:
                self.inp.pair_corr.append((0, 0.0, []))
        # End reading

    def run(self) -> None:
        """Execute the LEAPR calculation.

        This routine iterates over all requested scatterers and
        temperatures, computing the scattering law ``ssm``.  At
        present only the continuous phonon expansion is implemented
        (cards 10–13).  Attempts to invoke unimplemented features
        (translational contributions, discrete oscillators, pair
        correlations or cold hydrogen/deuterium) will result in a
        ``NotImplementedError``.
        """
        ntemp = self.inp.ntempr
        nalpha = self.inp.nalpha
        nbeta = self.inp.nbeta
        # convert alpha/beta to numpy arrays for efficiency
        alpha = np.array(self.inp.alpha, dtype=float)
        beta = np.array(self.inp.beta, dtype=float)
        # main loop over temperatures
        for t_index in range(ntemp):
            temp = self.inp.temperatures[t_index]
            # unwrap continuous distribution parameters
            delta1, p1_vals, twt, c, tbeta = self.inp.continuous[t_index]
            # store local copies for use in helper routines
            self.delta1 = delta1
            self.p1 = np.array(p1_vals, dtype=float)
            self.np1 = len(p1_vals)
            self.twt = twt
            self.c = c
            self.tbeta = tbeta
            # compute continuous (phonon expansion) contribution
            self._contin(temp, t_index, alpha, beta)
            # after contin the scattering law is symmetric; copy to ssp
            self.ssp[t_index] = self.ssm[t_index].copy()
            # translational contributions (diffusion/free gas)
            if self.twt > 0.0:
                # apply diffusion/free gas convolution
                self._trans(t_index, alpha, beta)
                # maintain symmetry for ssp
                self.ssp[t_index] = self.ssm[t_index].copy()
            # discrete oscillators
            bdel, adel = self.inp.oscillators[t_index]
            if len(bdel) > 0:
                # apply discrete oscillator convolution
                self._discrete(t_index, alpha, beta, np.array(bdel, dtype=float), np.array(adel, dtype=float))
                # maintain symmetry for ssp
                self.ssp[t_index] = self.ssm[t_index].copy()
            # special cold hydrogen/deuterium treatment
            if self.inp.ncold > 0:
                # retrieve pair correlation parameters for coldh
                nka, dka, ska_vals = self.inp.pair_corr[t_index]
                ska_arr = np.array(ska_vals, dtype=float) if nka > 0 else np.array([])
                # call cold hydrogen/deuterium routine
                self._coldh(t_index, temp, ska_arr, dka, alpha, beta)
            # pair correlation options (Vineyard/Sköld)
            nka, dka, ska_vals = self.inp.pair_corr[t_index]
            # if pair correlation table present and ncold == 0
            if nka > 0 and self.inp.ncold == 0:
                # obtain coherent fraction (cfrac) if provided
                cfrac = 0.0
                if t_index < len(self.inp.cfracs):
                    cfrac = self.inp.cfracs[t_index]
                # convert ska_vals to numpy array
                ska_arr = np.array(ska_vals, dtype=float)
                # apply skold approximation only for nsk > 0
                if self.inp.nsk > 0:
                    self._skold(t_index, temp, ska_arr, dka, cfrac, alpha, beta)
        # finished
        # --if secondary scatterer mixing is requested (nss > 0) and b7 <= 0,
        # perform a second pass to compute the secondary scattering law and
        # combine it with the principal scatterer.  The mixing rule in the
        # original Fortran code multiplies the secondary contribution by the
        # ratio of bound cross sections (sbs/sb) and then adds the principal
        # contribution.  We implement this by running a second instance of
        # Leapr with scaled parameters and combining the results.  Note that
        # this approximation assumes that the input decks for the two
        # scatterers share the same phonon distribution; users requiring
        # distinct distributions should construct two separate LeaprInput
        # objects and perform their own mixing externally.
        if self.inp.nss != 0 and self.inp.b7 <= 0:
            # bound cross section for principal and secondary scatterers
            # sb = spr * ((1 + awr) / awr)**2
            sb = self.inp.spr * ((1.0 + self.inp.awr) / self.inp.awr) ** 2
            # sbs = sps * ((1 + aws) / aws)**2
            # avoid division by zero if aws is zero (unphysical but guard anyway)
            sbs = 0.0
            if self.inp.aws != 0.0:
                sbs = self.inp.sps * ((1.0 + self.inp.aws) / self.inp.aws) ** 2
            # ratio for mixing
            srat = 0.0
            if sb != 0.0:
                srat = sbs / sb
            # save principal results
            principal_ssm = self.ssm.copy()
            principal_ssp = self.ssp.copy()
            principal_dwpix = self.dwpix.copy()
            principal_tempf = self.tempf.copy()
            # build a secondary input by copying the current input and
            # adjusting parameters for the secondary scatterer.
            import copy
            sec_inp = copy.deepcopy(self.inp)
            # The atomic weight ratio and scattering cross section are
            # replaced with those of the secondary scatterer.  The number
            # of scattering atoms per molecule (npr) is likewise set to mss.
            sec_inp.awr = self.inp.aws if self.inp.aws != 0.0 else self.inp.awr
            sec_inp.spr = self.inp.sps
            sec_inp.npr = self.inp.mss
            # Scale α values by the ratio of atomic weights (arat = aws/awr).
            # In the Fortran code the scaling is applied on the fly when
            # evaluating the convolution kernels.  Here we pre‑scale the
            # α grid for the secondary run to mimic that behaviour.
            arat = 1.0
            if self.inp.awr != 0.0:
                arat = self.inp.aws / self.inp.awr if self.inp.aws != 0.0 else 1.0
            sec_inp.alpha = [a / arat for a in self.inp.alpha]
            # Clear any further secondary scatterer flags on the
            # secondary input to avoid infinite recursion.  The
            # secondary run should not spawn additional mixing.
            sec_inp.nss = 0
            sec_inp.b7 = self.inp.b7
            # Create a new LEAPR instance for the secondary run and
            # execute it.  The secondary run reuses the same phonon
            # distribution, oscillator list and pair‑correlation data as
            # the principal scatterer.  Users requiring distinct
            # distributions for the secondary scatterer should supply a
            # separate input deck instead of relying on this mixing.
            sec = Leapr(sec_inp)
            # Copy over temperature‑dependent data structures.  The
            # deep copy above ensures that continuous, oscillators,
            # pair_corr and cfracs lists are shared; we do not need to
            # reparse the input file.  Simply call run() to fill sec.ssm.
            sec.run()
            # Combine the results: S_final = S_principal + srat * S_secondary
            # For the symmetric (negative‑β) branch
            self.ssm = principal_ssm + srat * sec.ssm
            # For the positive‑β branch (used in cold H/D and skewed cases)
            self.ssp = principal_ssp + srat * sec.ssp
            # Effective Debye–Waller factors and temperature factors are
            # stored separately for principal and secondary scatterers in
            # the Fortran code (dwp1/tempf1 for principal and dwpix/tempf
            # for secondary).  Here we preserve both sets for potential
            # downstream use by exposing them as attributes.  The
            # secondary values overwrite the current dwpix/tempf, while
            # principal values are stored in new attributes dwp1/tempf1.
            self.dwp1 = principal_dwpix
            self.tempf1 = principal_tempf
            # After mixing, dwpix and tempf refer to the secondary
            # scatterer.  The endf output routine (if implemented) would
            # average dwpix and dwp1 when computing coherent elastic
            # scattering; users can access both arrays directly.

        # At this point all scattering contributions have been computed,
        # including any secondary scatterer mixing.  If the input deck
        # requests coherent elastic scattering (iel>0), perform the
        # Bragg edge thinning now.  This mirrors the behaviour of
        # the Fortran ``endout`` routine, which thins the Bragg edge list
        # prior to writing elastic data.  The results are stored as
        # attributes on the instance for later use in ENDF output or
        # user‑level queries.
        if self.inp.iel > 0:
            edges, cumulatives = self.bragg_edge_thinning()
            # Save the thinned edge energies and cumulative integrals
            self.bragg_edges = edges
            self.bragg_cumint = cumulatives

    # ------------------------------------------------------------------
    # Output routines

    def export_to_ascii(self, filename: str) -> None:
        """
        Export the computed scattering law to a simple ASCII text file.

        The output format is not a strict ENDF‑6 file; rather it lists
        temperature‑labelled blocks of β, α, and S(α,β) values.  Each
        line contains the β value, the α value, the scattering law for
        the negative‑β branch (``ssm``) and, if available, the positive‑β
        branch (``ssp``).  This routine can be used to inspect the
        calculated S(α,β) without requiring a full ENDF writer.  For
        proper ENDF‑6 formatting, a separate implementation of the
        ``endout`` routine would be required.

        Parameters
        ----------
        filename : str
            Path to the output text file.
        """
        ntemp = self.inp.ntempr
        nalpha = self.inp.nalpha
        nbeta = self.inp.nbeta
        # open the file for writing
        with open(filename, 'w') as f:
            for t_index in range(ntemp):
                temp = self.inp.temperatures[t_index]
                f.write(f"# Temperature {temp}\n")
                for i_beta in range(nbeta):
                    beta_val = self.inp.beta[i_beta]
                    for j_alpha in range(nalpha):
                        alpha_val = self.inp.alpha[j_alpha]
                        ssm_val = self.ssm[t_index, i_beta, j_alpha]
                        # for symmetrical cases ssp mirrors ssm
                        ssp_val = self.ssp[t_index, i_beta, j_alpha] if hasattr(self, 'ssp') else ssm_val
                        f.write(f"{beta_val:.6e} {alpha_val:.6e} {ssm_val:.6e} {ssp_val:.6e}\n")
                f.write("\n")

    #
    # ENDF‑6 helper functions
    #
    def _format_endf_record(self, values: list[float], mat: int, mf: int, mt: int, line_no: int) -> str:
        """
        Format up to six numerical values and header identifiers into a
        single 80‑character ENDF‑6 text record.

        Each numeric value is formatted in 11‑column scientific notation
        (E11.5) or integer notation if it is effectively integral.
        The MAT, MF and MT identifiers occupy columns 67–75, and the
        record index occupies columns 76–80.  If fewer than six values
        are provided, the remaining fields are filled with blanks.

        Parameters
        ----------
        values : list of float
            A list of up to six numbers to be written in the leftmost
            columns of the record.
        mat : int
            Material identifier (MAT field).
        mf : int
            File number (MF field).
        mt : int
            Section number (MT field).
        line_no : int
            Sequential line number for this record (will be right‑justified
            in five columns).

        Returns
        -------
        str
            An 80‑character string representing one ENDF‑6 record.
        """
        # Prepare six 11‑character fields for the numeric values
        fields = []
        for i in range(6):
            if i < len(values):
                val = values[i]
                # Use integer format when the value is an integer (within tolerance)
                if isinstance(val, int) or (abs(val - round(val)) < 1e-9 and abs(val) < 1e6):
                    fields.append(f"{int(round(val)):11d}")
                else:
                    fields.append(f"{val:11.5E}")
            else:
                fields.append(" " * 11)
        # Concatenate the six fields with no separator
        body = ''.join(fields)
        # Format the MAT, MF, MT and line number fields
        mat_str = f"{mat:4d}"
        mf_str = f"{mf:2d}"
        mt_str = f"{mt:3d}"
        line_str = f"{line_no:5d}"
        # Construct full 80‑character record
        record = f"{body}{mat_str}{mf_str}{mt_str}{line_str}"
        # Ensure the record is exactly 80 characters
        return record[:80]

    def export_to_endf(self, filename: str) -> None:
        """
        Write the scattering law to an ENDF‑6 formatted file.

        This method first uses the built‑in ENDF helper routines
        (``contio``, ``tab1io``, etc.) to write a complete thermal
        scattering law file in the ENDF‑6 format.  The file includes
        the File‑1 header, incoherent and coherent elastic sections
        (MF=7/MT=2) as appropriate, and the inelastic scattering
        section (MF=7/MT=4) for all temperatures and β values.

        If the `endf_parserpy` package is available, the method will
        automatically re‑parse the generated lines and rewrite the file
        using ``endf_parserpy``.  This extra step produces a fully
        compliant ENDF‑6 file with correct line numbering and section
        terminators.  If the package is not installed, the
        internally generated file is left unchanged.

        Parameters
        ----------
        filename : str
            Path to the output file.  The file will be created or
            overwritten.
        """
        # First write the ENDF content to the requested filename using
        # the built‑in helper routines.  This produces a valid ENDF
        # structure with 80‑character records but without the
        # dictionary/comment machinery of NJOY.
        with open(filename, 'w') as f:
            reset_ns()
            set_mat(self.inp.mat)
            # ----- File 1 header (MF=1/MT=451) -----
            set_mf(1)
            set_mt(451)
            # New material record
            contio(f, self.inp.za, self.inp.awr, -1, 0, 0, 0)
            # Indicate six text lines follow
            contio(f, 0.0, 0.0, 0, 0, 0, 6)
            # Flag record
            contio(f, 1.0, 0.0, 0, 0, 12, 6)
            # Set C6 to 3 if any elastic section is present (iel!=0), else 2
            c6 = 3 if self.inp.iel != 0 else 2
            contio(f, 0.0, 0.0, 0, 0, 0, c6)
            # Six blank text lines
            for _ in range(6):
                contio(f, 0.0, 0.0, 0, 0, 0, 0)
            # Terminate file 1
            asend(f)
            afend(f)
            # ----- Reset state for file 7 -----
            reset_ns()
            set_mat(self.inp.mat)
            # ----- Incoherent elastic (MF=7/MT=2) -----
            if self.inp.iel < 0:
                set_mf(7)
                set_mt(2)
                # Header: incoherent elastic L1=2
                contio(f, self.inp.za, self.inp.awr, 2, 0, 0, 0)
                # Compute bound cross section sb*npr
                if self.inp.awr != 0.0:
                    sb = self.inp.spr * ((1.0 + self.inp.awr) / self.inp.awr) ** 2
                else:
                    sb = 0.0
                c1 = sb * self.inp.npr
                c2 = 0.0
                l1 = 0; l2 = 0
                # Create temperature and Debye–Waller arrays
                temps = list(self.inp.temperatures)
                dwps: List[float] = []
                for ti, T in enumerate(temps):
                    if T <= 0 or self.inp.awr == 0.0:
                        dwps.append(0.0)
                    else:
                        dwps.append(self.dwpix[ti] / (self.inp.awr * T * self.BK))
                if len(temps) == 1:
                    temps = [temps[0], temps[0]]
                    dwps = [dwps[0], dwps[0]]
                ndw = len(temps)
                nbt = [ndw]; intt = [2]
                tab1io(f, c1, c2, l1, l2, nbt, intt, temps, dwps)
                # Section end
                asend(f)
            # ----- Inelastic scattering (MF=7/MT=4) -----
            set_mf(7)
            set_mt(4)
            # Determine symmetry flag
            isym_val = 0
            if self.inp.ncold != 0:
                isym_val = 1
            if self.inp.isabt == 1:
                isym_val += 2
            # Header for File 7 inelastic
            contio(f, self.inp.za, self.inp.awr, 0, self.inp.lat, isym_val, 0)
            # Second header: number of beta points (mirrored if asymmetric)
            if isym_val % 2 == 0:
                nbeta_out = self.inp.nbeta
            else:
                nbeta_out = 2 * self.inp.nbeta - 1
            contio(f, 0.0, 0.0, 0, 0, 1, nbeta_out)
            # Loop over beta grid and write TAB1 records
            for i_beta in range(nbeta_out):
                if isym_val % 2 == 0:
                    beta_val = self.inp.beta[i_beta]
                else:
                    if i_beta < self.inp.nbeta - 1:
                        beta_val = -self.inp.beta[self.inp.nbeta - i_beta - 1]
                    else:
                        beta_val = self.inp.beta[i_beta - (self.inp.nbeta - 1)]
                for t_index, temp in enumerate(self.inp.temperatures):
                    x_vals = self.inp.alpha
                    if isym_val % 2 == 0:
                        y_vals = self.ssm[t_index, i_beta, :].tolist()
                    else:
                        if i_beta < self.inp.nbeta - 1:
                            y_vals = self.ssm[t_index, self.inp.nbeta - i_beta - 1, :].tolist()
                        else:
                            y_vals = self.ssp[t_index, i_beta - (self.inp.nbeta - 1), :].tolist()
                    nbt_alpha = [len(x_vals)]; intt_alpha = [2]
                    tab1io(f, temp, beta_val, 0, 0, nbt_alpha, intt_alpha, x_vals, y_vals)
            # End inelastic section
            asend(f)
            # ----- File end, material end, tape end -----
            afend(f)
            amend(f)
            atend(f)
            # ----- Coherent elastic (MF=7/MT=2) -----
            if self.inp.iel > 0:
                edges = getattr(self, 'bragg_edges', [])
                cumulatives = getattr(self, 'bragg_cumint', [])
                if edges and cumulatives:
                    set_mf(7)
                    set_mt(2)
                    contio(f, self.inp.za, self.inp.awr, 1, 0, 0, 0)
                    for t_index, temp in enumerate(self.inp.temperatures):
                        x_vals = edges
                        y_vals = cumulatives[t_index] if t_index < len(cumulatives) else []
                        if not y_vals:
                            continue
                        nbt = [len(x_vals)]; intt = [1]
                        tab1io(f, temp, 0.0, 0, 0, nbt, intt, x_vals, y_vals)
                    asend(f)
        # If the endf_parserpy package is available, re‑parse and rewrite
        try:
            from endf_parserpy import EndfParserFactory
            # Parse the file we just wrote
            parser = EndfParserFactory.create(select="python")
            endf_dict = parser.parsefile(filename)
            # Write it back using endf_parserpy to ensure full compliance
            parser.writefile(filename, endf_dict, overwrite=True)
        except Exception:
            # If parsing fails or package is unavailable, leave file as is
            pass

    # ------------------------------------------------------------------
    # Core numerical routines

    def _contin(self, temp: float, t_index: int, alpha: np.ndarray, beta: np.ndarray) -> None:
        """
        Compute the continuous (phonon expansion) part of S(α,β) for a
        given temperature.

        This method corresponds to the Fortran subroutine ``contin``.
        It relies on helper functions ``_start``, ``_fsum``, ``_terpt``
        and ``_convol`` to perform the phonon expansion.  The result
        for the specified temperature index is written into
        ``self.ssm[t_index]``.

        Parameters
        ----------
        temp : float
            Temperature in kelvin.  A negative value in the input deck
            indicates that the previous parameters should be reused.
        t_index : int
            Index of the temperature in the input arrays.
        alpha : ndarray
            1‑D array of α values.
        beta : ndarray
            1‑D array of β values.
        """
        nalpha = self.inp.nalpha
        nbeta = self.inp.nbeta
        nphon = self.inp.nphon
        iprint = self.inp.iprint
        # Use absolute temperature to determine effective temperature in eV
        self.tev = self.BK * abs(temp)
        # call start() to compute initial distribution and scaling
        p, deltab = self._start(t_index)
        self.deltab = deltab
        # scale factor for α and β if required
        sc = 1.0
        if self.inp.lat == 1:
            # in Fortran: therm/tev with therm=0.0253 eV
            sc = 0.0253 / self.tev
        # allocate working arrays
        # tlast will hold t_n from previous order, initialised with t1
        tlast = np.copy(p)
        # initialise S(α,β) array for this temperature
        sab = np.zeros((nbeta, nalpha), dtype=float)
        # compute zeroth order (n=1) contribution
        for j in range(nalpha):
            al = alpha[j] * sc  # α scaled by temperature and scatterer ratio
            # compute exponent terms for this α
            ex_const = -self.f0 * al
            # log term used to build up α^n / n! factors
            xa = math.log(al * self.f0) if al > 0 else -np.inf
            for k in range(nbeta):
                be = beta[k] * sc
                st = self._terpt(p, be)
                add = st * math.exp(ex_const + xa) if st > 0 else 0.0
                sab[k, j] = add
        # perform phonon expansion up to order nphon
        npl = len(p)
        for n in range(2, nphon + 1):
            # convolution t_n = t1 * t_{n-1}
            tnext = self._convol(p, tlast, deltab)
            # accumulate contributions for each α,β
            for j in range(nalpha):
                al = alpha[j] * sc
                # update xa for order n: xa += log(al * f0 / n)
                # we track xa separately per α to avoid recomputing log factorials
                # For simplicity we recompute from scratch here
                if al > 0 and self.f0 > 0:
                    xa = math.log((al * self.f0) ** n / math.factorial(n))
                else:
                    xa = -np.inf
                ex_const = -self.f0 * al
                for k in range(nbeta):
                    be = beta[k] * sc
                    st = self._terpt(tnext, be)
                    add = st * math.exp(ex_const) * (al * self.f0) ** n / math.factorial(n) if st > 0 else 0.0
                    sab[k, j] += add
            # prepare for next order
            tlast = tnext
        # store results
        self.ssm[t_index] = sab

    def _start(self, t_index: int) -> Tuple[np.ndarray, float]:
        """
        Compute integral functions of the phonon frequency distribution.

        This method corresponds to the Fortran subroutine ``start``.
        It takes the input ``p1`` array (the phonon density at
        uniformly spaced β values) and constructs the initial
        tabulation of t₁(β).  It also calculates the Debye–Waller
        parameter ``f0`` and the effective temperature factor ``tbar``.

        Parameters
        ----------
        t_index : int
            Index of the current temperature in the input arrays.

        Returns
        -------
        p : ndarray
            Array of t₁(β) values on a grid with spacing ``deltab``.
        deltab : float
            The β spacing used for the convolution routines.
        """
        # dereference input fields
        delta1 = self.delta1
        p1 = self.p1
        np1 = self.np1
        tbeta = self.tbeta
        # compute β spacing in the scaled grid
        deltab = delta1 / self.tev
        # copy input spectrum into p array
        p = p1.copy()
        # initial normalisation of p according to the Fortran code
        # convert p(e) into p(β) by dividing by sinh(β/2) factor and normalising
        # as in start() from Fortran
        # compute intermediate arrays
        u = deltab
        v = math.exp(deltab / 2.0)
        # first element special case
        p[0] = p[1] / deltab ** 2
        vv = v
        for j in range(1, np1):
            p[j] = p[j] / (u * (vv - 1.0 / vv))
            vv *= v
            u += deltab
        # determine normalising constant an
        tau = 0.5  # fixed in Fortran
        an = self._fsum(1, p, deltab) / tbeta
        # normalise p and compute Debye–Waller integrals
        for i in range(np1):
            p[i] /= an
        # compute f0 and tbar
        self.f0 = self._fsum(0, p, deltab)
        self.tbar = self._fsum(2, p, deltab) / (2.0 * tbeta)
        # convert p(β) into t₁(β)
        for i in range(np1):
            be = deltab * i
            p[i] = p[i] * math.exp(be / 2.0) / self.f0
        # store Debye–Waller and effective temp
        self.dwpix[t_index] = self.f0
        self.tempf[t_index] = self.tbar * self.inp.temperatures[t_index]
        return p, deltab

    def _fsum(self, n: int, p: np.ndarray, deltab: float) -> float:
        """
        Compute integral of p(β) β^n times hyperbolic function.

        This replicates the Fortran function ``fsum``.  The integral
        computed is

            ∫₀^∞ 2 p(β) β^n [cosh(τ β) or sinh(τ β)] dβ

        where the hyperbolic function depends on the parity of n.
        A trapezoidal rule on a uniform grid is used.  In this port
        τ=½ is fixed, matching the original code.

        Parameters
        ----------
        n : int
            Power of β in the integrand.
        p : ndarray
            Array of p(β) values on a uniform β grid.
        deltab : float
            Grid spacing in β.

        Returns
        -------
        float
            Approximation of the integral.
        """
        tau = 0.5
        edsq = math.exp(deltab * tau / 2.0)
        v = 1.0
        an = 1.0 - 2.0 * (n % 2)  # +1 for even n, −1 for odd n
        be = 0.0
        fs = 0.0
        w = 1.0
        npnt = len(p)
        for ij in range(npnt):
            if n > 0:
                w = be ** n
            # hyperbolic term: 2 p(β) [cosh(tau β) or sinh(tau β)]
            ff = ((p[ij] * v) * v + (p[ij] * an / v) / v) * w
            if ij == 0 or ij == npnt - 1:
                ff /= 2.0
            fs += ff
            be += deltab
            v *= edsq
        return fs * deltab

    def _terpt(self, tn: np.ndarray, be: float) -> float:
        """
        Linear interpolation in a table of tₙ(β).

        Equivalent to the Fortran function ``terpt``.  If the
        requested β lies beyond the end of the tabulation, zero is
        returned.

        Parameters
        ----------
        tn : ndarray
            Array of tₙ(β) values on a uniform β grid.
        be : float
            Requested β value.

        Returns
        -------
        float
            Interpolated tₙ(β).
        """
        ntn = len(tn)
        delta = self.deltab
        if be > ntn * delta:
            return 0.0
        i = int(be / delta)
        if i < ntn - 1:
            bt = i * delta
            btp = bt + delta
            return tn[i] + (be - bt) * (tn[i + 1] - tn[i]) / (btp - bt)
        return 0.0

    def _terps(self, sd: np.ndarray, delta: float, be: float) -> float:
        """
        Logarithmic interpolation used for translational contributions.

        This replicates the Fortran function ``terps``.  The array
        ``sd`` contains values of s(β) on a uniform grid of spacing
        ``delta``.  If the input β value is larger than the tabulated
        range, zero is returned.  Otherwise a logarithmic interpolation
        between the two nearest points is performed.

        Parameters
        ----------
        sd : ndarray
            Array of s(β) values.
        delta : float
            Grid spacing used when constructing ``sd``.
        be : float
            Requested β value.

        Returns
        -------
        float
            Interpolated s(β) value.
        """
        nsd = len(sd)
        slim = -225.0
        if be > delta * nsd:
            return 0.0
        i = int(be / delta)
        if i < nsd - 1:
            bt = i * delta
            btp = bt + delta
            st = slim if sd[i] <= 0.0 else math.log(sd[i])
            stp = slim if sd[i + 1] <= 0.0 else math.log(sd[i + 1])
            stt = st + (be - bt) * (stp - st) / (btp - bt)
            return math.exp(stt) if stt > slim else 0.0
        return 0.0

    def _besk1(self, x: float) -> float:
        """
        Compute the modified Bessel function K₁(x) with the same
        normalisation used in the Fortran routine ``besk1``.

        For x ≤ 1, a series expansion is used; for x > 1, an
        asymptotic expansion is used.  This implementation is a direct
        translation of the constants and logic from the NJOY Fortran
        code.  The exponential part for x > 1 is omitted, mirroring
        the behaviour in ``stable``.
        """
        # coefficients taken from the Fortran code
        c0 = 0.125
        c1 = 0.442850424
        c2 = 0.584115288
        c3 = 6.070134559
        c4 = 17.864913364
        c5 = 48.858995315
        c6 = 90.924600045
        c7 = 113.795967431
        c8 = 85.331474517
        c9 = 32.00008698
        c10 = 3.999998802
        c11 = 1.304923514
        c12 = 1.47785657
        c13 = 16.402802501
        c14 = 44.732901977
        c15 = 115.837493464
        c16 = 198.437197312
        c17 = 222.869709703
        c18 = 142.216613971
        c19 = 40.000262262
        c20 = 1.999996391
        c21 = 1.0
        c22 = 0.5
        c23 = 0.5772156649
        c24 = 1.0
        c25 = 0.0108241775
        c26 = 0.0788000118
        c27 = 0.2581303765
        c28 = 0.5050238576
        c29 = 0.663229543
        c30 = 0.6283380681
        c31 = 0.4594342117
        c32 = 0.2847618149
        c33 = 0.1736431637
        c34 = 0.1280426636
        c35 = 0.1468582957
        c36 = 0.4699927013
        c37 = 1.2533141373
        test = 1.0
        if x <= test:
            v = c0 * x
            u = v * v
            bi1 = ((((((((c1 * u + c2) * u + c3) * u + c4) * u + c5) * u + c6) * u + c7) * u + c8) * u + c9) * u + c10
            bi1 *= v
            bi3 = ((((((((c11 * u + c12) * u + c13) * u + c14) * u + c15) * u + c16) * u + c17) * u + c18) * u + c19) * u + c20
            return c21 / x + bi1 * (math.log(c22 * x) + c23) - v * bi3
        else:
            u = c24 / x
            bi3 = (((((((((( -c25 * u + c26) * u - c27) * u + c28) * u - c29) * u + c30) * u - c31) * u + c32) * u - c33) * u + c34) * u - c35) * u + c36
            return math.sqrt(u) * (bi3 + c37 * u)

    def _stable(self, al: float, delta: float) -> Tuple[np.ndarray, int]:
        """
        Construct a table of S_d(β) for the diffusion/free‑gas
        translational contribution.

        This is a translation of the Fortran subroutine ``stable``.
        It generates a one‑sided (β ≥ 0) table of diffusion/free‑gas
        scattering contributions.  The computation stops when the
        values drop below 1e‑7 of the first value or when a maximum
        length is reached.

        Parameters
        ----------
        al : float
            Scaled α value (α * sc / arat in the Fortran code).
        delta : float
            β spacing for the convolution grid.

        Returns
        -------
        sd : ndarray
            Array of S_d(β) values on a uniform grid starting at β=0.
        nsd : int
            Number of points in ``sd``.
        """
        tiny = 1.0e-30
        eps = 1.0e-7
        # allocate a reasonably large array; the original code uses ndmax
        sd_vals: List[float] = []
        # diffusion branch
        if self.c != 0.0:
            d = self.twt * self.c
            c2 = math.sqrt(self.c * self.c + 0.25)
            c3 = 2.0 * d * al
            c4 = c3 * c3
            c8 = c2 * c3 / math.pi
            c3 = 2.0 * d * self.c * al
            be = 0.0
            j = 0
            while True:
                c6 = math.sqrt(be * be + c4)
                c7 = c6 * c2
                if c7 <= 1.0:
                    c5 = c8 * math.exp(c3 + be / 2.0)
                else:
                    c5 = 0.0
                    ex = c3 - c7 + be / 2.0
                    c5 = c8 * math.exp(ex)
                val = 0.0
                if c6 != 0:
                    val = c5 * self._besk1(c7) / c6
                sd_vals.append(val)
                # termination conditions every other point
                j += 1
                be += delta
                if j % 2 == 0:
                    if j >= 1999:
                        break
                    if sd_vals[0] * eps >= sd_vals[-1]:
                        break
        else:
            # free gas branch
            be = 0.0
            j = 0
            wal = self.twt * al
            while True:
                # avoid division by zero when wal is zero
                if wal == 0.0:
                    val = 0.0
                else:
                    ex = - (wal - be) ** 2 / (4.0 * wal)
                    val = math.exp(ex) / math.sqrt(4.0 * self.PI * wal)
                sd_vals.append(val)
                j += 1
                be += delta
                if j % 2 == 0:
                    if j >= 1999:
                        break
                    if sd_vals[0] * eps >= sd_vals[-1]:
                        break
        sd_array = np.array(sd_vals, dtype=float)
        return sd_array, len(sd_array)

    def _sbfill(self, be: float, s: np.ndarray, betan: np.ndarray, nbeta: int, delta: float, nbt: int) -> np.ndarray:
        """
        Generate s(β) on a new energy grid for convolution with a
        diffusion/free‑gas shape.

        This replicates the Fortran subroutine ``sbfill``.  A new
        table ``sb`` of length ``2*nbt-1`` is constructed covering
        β in [−be−(nbt−1)δ, −be+(nbt−1)δ].  The values are
        interpolated from the asymmetric scattering function stored in
        ``s``.  If β lies outside the tabulation, zero is used.

        Parameters
        ----------
        be : float
            β value at the centre of the convolution window.
        s : ndarray
            Array of asymmetric scattering values (length ``nbeta``).
        betan : ndarray
            Array of β grid values (length ``nbeta``).
        nbeta : int
            Number of β values.
        delta : float
            Grid spacing for the convolution window.
        nbt : int
            Number of points in the diffusion/free‑gas table ``sd``.

        Returns
        -------
        sb : ndarray
            Array of interpolated scattering values on the new grid
            (length ``2*nbt-1``).
        """
        # set up interpolation range
        bmin = -be - (nbt - 1) * delta
        bmax = -be + (nbt - 1) * delta + delta / 100.0
        # precompute logs for interpolation
        s_log = np.where(s > 0, np.log(s), -225.0)
        # allocate output array
        m = 2 * nbt - 1
        sb = np.zeros(m, dtype=float)
        # pointer for betan search
        j = nbeta - 1
        idx = 0
        bet = bmin
        while idx < m:
            b = abs(bet)
            # find bracketing interval in betan
            # move j down until betan[j-1] <= b <= betan[j]
            # or until boundaries
            if b > betan[j]:
                # move up in betan if possible
                while j < nbeta - 1 and b > betan[j]:
                    j += 1
            else:
                while j > 0 and b < betan[j - 1]:
                    j -= 1
            # perform interpolation if within range
            val = 0.0
            if 0 < j < nbeta:
                st = s_log[j]
                stm = s_log[j - 1]
                denom = betan[j - 1] - betan[j]
                if denom != 0:
                    sb_val = st + (b - betan[j]) * (stm - st) / denom
                else:
                    sb_val = st
                if bet > 0:
                    sb_val -= bet
                val = math.exp(sb_val) if sb_val > -225.0 else 0.0
            sb[idx] = val
            # increment
            bet += delta
            idx += 1
        return sb

    # ------------------------------------------------------------------
    # Discrete oscillator support routines (partially ported)

    def _bfact(self, x: float, dwc: float, betai: float) -> Tuple[float, np.ndarray, np.ndarray]:
        """
        Compute Bessel function terms for discrete oscillators.

        This is a translation of the Fortran subroutine ``bfact``.  It
        returns the zero‑order Bessel term and two arrays ``bplus`` and
        ``bminus`` corresponding to positive and negative orders up to
        ``imax`` (=50).  The arrays have length 50, indexed from 0 to
        49 in Python, corresponding to 1 through 50 in the Fortran code.

        Parameters
        ----------
        x : float
            Argument used in the Bessel function recurrence.
        dwc : float
            Debye–Waller exponent term (α * dbw(i) in the original code).
        betai : float
            Oscillator energy divided by kT (β value).

        Returns
        -------
        bzero : float
            Zero‑order term.
        bplus : ndarray
            Array of positive order Bessel terms (length 50).
        bminus : ndarray
            Array of negative order Bessel terms (length 50).
        """
        # constants from the Fortran implementation
        c0 = 3.75
        c1 = 1.0
        c2 = 3.5156229
        c3 = 3.0899424
        c4 = 1.2067492
        c5 = 0.2659732
        c6 = 0.0360768
        c7 = 0.0045813
        c8 = 0.39894228
        c9 = 0.01328592
        c10 = 0.00225319
        c11 = 0.00157565
        c12 = 0.00916281
        c13 = 0.02057706
        c14 = 0.02635537
        c15 = 0.01647633
        c16 = 0.00392377
        c17 = 0.5
        c18 = 0.87890594
        c19 = 0.51498869
        c20 = 0.15084934
        c21 = 0.02658733
        c22 = 0.00301532
        c23 = 0.00032411
        c24 = 0.02282967
        c25 = 0.02895312
        c26 = 0.01787654
        c27 = 0.00420059
        c28 = 0.39894228
        c29 = 0.03988024
        c30 = 0.00362018
        c31 = 0.00163801
        c32 = 0.01031555
        big = 1.0e10
        tiny = 1.0e-30
        # compute modified Bessel functions I0 and I1 (bessi0 and bessi1)
        y = x / c0
        # I0
        if y <= 1.0:
            u = y * y
            bessi0 = (((((c1 + u * c2) + u * u * c3) + u * u * u * c4) + u * u * u * u * c5) + u * u * u * u * u * c6) + u * u * u * u * u * u * c7
        else:
            v = 1.0 / y
            bessi0 = (c8 + v * (c9 + v * (c10 + v * (-c11 + v * (c12 + v * (-c13 + v * (c14 + v * (-c15 + v * c16)))))))) / math.sqrt(x)
        # I1
        if y <= 1.0:
            u = y * y
            bessi1 = (c17 + u * (c18 + u * (c19 + u * (c20 + u * (c21 + u * (c22 + u * c23)))))) * x
        else:
            v = 1.0 / y
            bessi1 = c24 + v * (-c25 + v * (c26 - v * c27))
            bessi1 = c28 + v * (-c29 + v * (-c30 + v * (c31 + v * (-c32 + v * bessi1))))
            bessi1 /= math.sqrt(x)
        # reverse recurrence to compute higher order Bessel functions
        imax = 50
        bn = [0.0] * (imax + 1)
        bn[imax] = 0.0
        bn[imax - 1] = 1.0
        i = imax - 1
        while i > 1:
            bn[i - 1] = bn[i + 1] + i * (2.0 / x) * bn[i]
            i -= 1
            if bn[i] >= big:
                # scale down to avoid overflow
                for j in range(i, imax + 1):
                    bn[j] /= big
        # normalise so that bn[1] = bessi1
        rat = bessi1 / bn[1]
        for i in range(1, imax + 1):
            bn[i] *= rat
            if bn[i] < tiny:
                bn[i] = 0.0
        # apply exponential terms to Bessel functions
        bzero = 0.0
        bplus = np.zeros(imax, dtype=float)
        bminus = np.zeros(imax, dtype=float)
        ycheck = y <= 1.0
        if ycheck:
            bzero = bessi0 * math.exp(-dwc)
            for i in range(1, imax + 1):
                idx = i - 1
                if bn[i] != 0.0:
                    # positive order
                    arg = -dwc - i * betai / 2.0
                    val = math.exp(arg) * bn[i]
                    bplus[idx] = val if val >= tiny else 0.0
                    # negative order
                    arg = -dwc + i * betai / 2.0
                    val = math.exp(arg) * bn[i]
                    bminus[idx] = val if val >= tiny else 0.0
        else:
            bzero = bessi0 * math.exp(-dwc + x)
            for i in range(1, imax + 1):
                idx = i - 1
                if bn[i] != 0.0:
                    # positive order
                    arg = -dwc - i * betai / 2.0 + x
                    val = math.exp(arg) * bn[i]
                    bplus[idx] = val if val >= tiny else 0.0
                    # negative order
                    arg = -dwc + i * betai / 2.0 + x
                    val = math.exp(arg) * bn[i]
                    bminus[idx] = val if val >= tiny else 0.0
        return bzero, bplus, bminus

    def _bfill(self, betan: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Construct arrays ``bex`` and ``rdbex`` used by the discrete
        oscillator interpolation.  This mirrors the Fortran subroutine
        ``bfill``.

        Parameters
        ----------
        betan : ndarray
            Array of β grid values (in ascending order).

        Returns
        -------
        bex : ndarray
            Extended β grid covering both negative and positive values.
        rdbex : ndarray
            Reciprocal differences between adjacent ``bex`` values.
        """
        nbeta = len(betan)
        # build negative branch (reverse of betan)
        bex_neg = -betan[::-1]
        # exclude zero from double counting if betan[0] is ~0
        if abs(betan[0]) <= 1.0e-9:
            bex_pos = betan[1:]
        else:
            bex_pos = betan
        bex = np.concatenate((bex_neg, bex_pos))
        # compute reciprocal differences
        rdbex = np.zeros(len(bex) - 1, dtype=float)
        for i in range(len(rdbex)):
            diff = bex[i + 1] - bex[i]
            rdbex[i] = 1.0 / diff if diff != 0.0 else 0.0
        return bex, rdbex

    def _exts(self, sexpb: np.ndarray, exb: np.ndarray, betan: np.ndarray) -> np.ndarray:
        """
        Extend the asymmetric scattering law to both ±β.

        This mirrors the Fortran subroutine ``exts``.  The input
        ``sexpb`` contains S(α,β) for negative β.  The array ``exb``
        contains exp(−β/2) factors for each β.  The output ``sex``
        contains the extended scattering law for both negative and
        positive β.

        Parameters
        ----------
        sexpb : ndarray
            Asymmetric scattering values for negative β (length ``nbeta``).
        exb : ndarray
            Precomputed exp(−β/2) values (length ``nbeta``).
        betan : ndarray
            β grid (length ``nbeta``).

        Returns
        -------
        sex : ndarray
            Extended scattering array (length ``2*nbeta - 1``).
        """
        nbeta = len(sexpb)
        # create arrays for negative and positive parts
        sex = np.zeros(2 * nbeta - 1, dtype=float)
        # negative part (reverse order)
        sex[:nbeta] = sexpb[::-1]
        # exclude zero from double counting
        if abs(betan[0]) <= 1.0e-9:
            start = nbeta
        else:
            # insert the zero index explicitly
            sex[nbeta] = sexpb[0]
            start = nbeta + 1
        # positive part
        for i in range(1, nbeta):
            idx = start + i - 1
            # multiply by exp(β/2) squared
            sex[idx] = sexpb[i] * exb[i] * exb[i]
        return sex

    def _sint(self, x: float, bex: np.ndarray, rdbex: np.ndarray, sex: np.ndarray, alpha: float, wt: float, tbart: float, betan: np.ndarray) -> float:
        """
        Interpolate/extrapolate scattering function for discrete oscillator convolution.

        This replicates the Fortran function ``sint``.  It uses the
        extended β grid (``bex``) and scattering values (``sex``) to
        perform interpolation, or falls back to the SCT (short
        collision time) approximation for values outside the tabulated
        range.

        Parameters
        ----------
        x : float
            Desired β value.
        bex : ndarray
            Extended β grid.
        rdbex : ndarray
            Reciprocal differences of ``bex``.
        sex : ndarray
            Extended scattering values.
        alpha : float
            Scaled α value.
        wt : float
            Weight factor ``tbeta + twt`` (total weight).
        tbart : float
            Effective temperature ratio for discrete oscillators.
        betan : ndarray
            Original β grid (positive values).

        Returns
        -------
        float
            Interpolated or approximated scattering value.
        """
        # SCT approximation if |x| > max(betan)
        if abs(x) > betan[-1]:
            if alpha <= 0.0:
                return 0.0
            ex = -(wt * alpha - abs(x)) ** 2 / (4.0 * wt * alpha * tbart)
            if x > 0:
                ex -= x
            return math.exp(ex) / (4.0 * math.pi * wt * alpha * tbart)
        # otherwise, perform interpolation on bex
        # locate interval via bisection
        k1 = 0
        k3 = len(bex) - 1
        # handle exact matches
        if x == bex[k3]:
            return sex[k3]
        # binary search
        while k3 - k1 > 1:
            k2 = (k1 + k3) // 2
            if x >= bex[k2]:
                k1 = k2
            else:
                k3 = k2
        # linear interpolation in log space
        ss1 = -225.0 if sex[k1] <= 0.0 else math.log(sex[k1])
        ss3 = -225.0 if sex[k3] <= 0.0 else math.log(sex[k3])
        ex = ((bex[k3] - x) * ss1 + (x - bex[k1]) * ss3) * rdbex[k1]
        return math.exp(ex) if ex > -225.0 else 0.0

    def _terpk(self, ska: np.ndarray, dka: float, be: float) -> float:
        """
        Interpolate in the S(kappa) table for a given kappa.

        This mirrors the Fortran function ``terpk`` used in the Sköld
        approximation.  The array ``ska`` holds values of S(kappa) at
        kappa values ``0, dka, 2*dka, ...``.  A simple linear
        interpolation is performed between adjacent points.  Values
        beyond the table (be > nka*dka) return 1.

        Parameters
        ----------
        ska : ndarray
            Tabulated S(kappa) values (length ``nka``).
        dka : float
            Spacing between kappa values (in inverse Angstroms).
        be : float
            Desired kappa value.

        Returns
        -------
        float
            Interpolated S(kappa) at the specified kappa.
        """
        nka = len(ska)
        # if request beyond table, return 1
        if be > nka * dka:
            return 1.0
        # integer index of lower grid point
        i0 = int(be / dka)
        # clamp at upper end
        if i0 >= nka - 1:
            return 1.0
        bt = i0 * dka
        btp = bt + dka
        # linear interpolation between ska[i0] and ska[i0+1]
        y1 = ska[i0]
        y2 = ska[i0 + 1]
        if btp == bt:
            return y1
        return y1 + (be - bt) * (y2 - y1) / (btp - bt)

    def _terp1(self, x1: float, y1: float, x2: float, y2: float, x: float, interp_code: int) -> float:
        """
        Interpolate a single point (x,y) between two points (x1,y1) and (x2,y2).

        This replicates the Fortran subroutine ``terp1`` from ``endf.f90``.
        Depending on the interpolation code, linear, logarithmic or
        exponential interpolation is used.  Only codes 1–5 are
        implemented; code 6 (Coulomb penetrability) is not needed in
        LEAPR.

        Parameters
        ----------
        x1, y1 : float
            Coordinates of the first point.
        x2, y2 : float
            Coordinates of the second point.
        x : float
            The x‑value at which to interpolate.
        interp_code : int
            Interpolation flag (1–5) as in ENDF manual.

        Returns
        -------
        float
            Interpolated y value at ``x``.
        """
        # handle degenerate case
        if x2 == x1:
            return y1
        # constant y or trivial cases
        if interp_code == 1 or y2 == y1 or x == x1:
            return y1
        # linear in x
        if interp_code == 2:
            return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
        # y linear in ln(x)
        if interp_code == 3:
            if x1 <= 0.0 or x2 <= 0.0 or x <= 0.0:
                return y1
            return y1 + math.log(x / x1) * (y2 - y1) / math.log(x2 / x1)
        # ln(y) linear in x
        if interp_code == 4:
            if y1 <= 0.0 or y2 <= 0.0:
                return y1
            return y1 * math.exp((x - x1) * math.log(y2 / y1) / (x2 - x1))
        # ln(y) linear in ln(x)
        if interp_code == 5:
            if y1 <= 0.0 or y2 <= 0.0 or x1 <= 0.0 or x2 <= 0.0 or x <= 0.0:
                return y1
            return y1 * math.exp(math.log(x / x1) * math.log(y2 / y1) / math.log(x2 / x1))
        # fallback to linear
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)


    def _trans(self, t_index: int, alpha: np.ndarray, beta: np.ndarray) -> None:
        """
        Apply translational (diffusion/free‑gas) contributions to S(α,β).

        This method is a translation of the Fortran subroutine ``trans``.
        It convolves the previously calculated continuous scattering law
        with a translational kernel representing either diffusion or
        free gas.  The results are accumulated in ``self.ssm``.

        Parameters
        ----------
        t_index : int
            Index of the current temperature.
        alpha : ndarray
            Array of α values.
        beta : ndarray
            Array of β values.
        """
        nbeta = self.inp.nbeta
        nalpha = self.inp.nalpha
        # scale factors
        sc = 1.0
        if self.inp.lat == 1:
            sc = 0.0253 / self.tev
        # loop over α
        for ialpha in range(nalpha):
            al = alpha[ialpha] * sc
            # choose β interval size for convolution
            # coefficients from Fortran
            c0 = 0.4
            c1 = 1.0
            c2 = 1.42
            c3 = 0.2
            c4 = 10.0
            ded = c0 * (self.twt * self.c * al) / math.sqrt(c1 + c2 * (self.twt * self.c * al) * self.c)
            if ded == 0.0:
                ded = c3 * math.sqrt(self.twt * al)
            deb = c4 * al * self.deltab
            # choose smaller of ded and deb
            delta = ded if ded < deb else deb
            # build diffusion/free‑gas table
            sd, nsd = self._stable(al, delta)
            if nsd <= 1:
                continue
            # copy original asymmetric S to temporary array
            betan = beta * sc
            ap = self.ssm[t_index, :, ialpha].copy()
            # convolution loop over β values
            for ibeta in range(nbeta):
                be = betan[ibeta]
                # prepare sb values on new grid
                sb = self._sbfill(be, ap, betan, nbeta, delta, nsd)
                # integrate diffusion/free‑gas kernel with scattering function
                s_val = 0.0
                m = nsd
                mid = m - 1
                # Simpson's rule integration as in Fortran
                for i in range(m):
                    f = 2.0 * ((i % 2) + 1)
                    if i == 0 or i == m - 1:
                        f = 1.0
                    bb = i * delta
                    # positive β contribution
                    s_val += f * sd[i] * sb[mid + i]
                    # negative β contribution
                    s_val += f * sd[i] * sb[mid - i] * math.exp(-bb)
                s_val *= delta / 3.0
                # add contribution from tail using interpolation
                st = self._terps(sd, delta, be)
                if st > 0.0:
                    s_val += math.exp(-al * self.f0) * st
                # update S(α,β)
                self.ssm[t_index, ibeta, ialpha] = s_val

    # ------------------------------------------------------------------
    # Discrete oscillator convolution (placeholder)

    def _discrete(self, t_index: int, alpha: np.ndarray, beta: np.ndarray, bdel: np.ndarray, adel: np.ndarray) -> None:
        """
        Convolve discrete oscillators with the continuous scattering law.

        This routine is a translation of the Fortran subroutine
        ``discre`` from the original LEAPR module.  It computes the
        contributions of discrete vibrational modes (delta functions)
        to the scattering law S(α,β) by convolving the previously
        calculated continuous S(α,β) with the appropriate Bessel
        expansions.  The resulting contributions are accumulated
        directly into ``self.ssm`` for the current temperature.

        Parameters
        ----------
        t_index : int
            Index of the current temperature.
        alpha : ndarray
            Array of α values.
        beta : ndarray
            Array of β values.
        bdel : ndarray
            Oscillator energies (eV) for this temperature.
        adel : ndarray
            Oscillator weights for this temperature.
        """
        # number of α and β points
        nalpha = self.inp.nalpha
        nbeta = self.inp.nbeta
        nd = len(bdel)
        if nd == 0:
            return
        # physical constants and small tolerances
        small = 1.0e-8
        vsmall = 1.0e-10
        tiny = 1.0e-20
        # scaling factor for β grid
        sc = 1.0
        if self.inp.lat == 1:
            sc = 0.0253 / self.tev
        # scaled β grid and exp(−β/2)
        betan = beta * sc
        exb = np.exp(-betan / 2.0)
        # build extended β grid and reciprocal differences once
        bex, rdbex = self._bfill(betan)
        # oscillator parameters
        # convert oscillator energies to dimensionless β shifts
        bdeln = np.array(bdel, dtype=float) / self.tev
        # allocate arrays for ar, dist and dbw per oscillator
        ar = np.zeros(nd, dtype=float)
        dist = np.zeros(nd, dtype=float)
        dbw = np.zeros(nd, dtype=float)
        # cumulative weights
        dwt = 0.0
        # contribution to effective temperature from discrete oscillators
        tsave = 0.0
        # update Debye–Waller parameter dwpix[t_index] as necessary
        for i in range(nd):
            adel_i = adel[i]
            bd_e = bdel[i]
            bdn = bdeln[i]
            dwt += adel_i
            # compute sinh and cosh terms
            eb = math.exp(bdn / 2.0)
            # protect against zero division
            if eb == 0.0:
                sn = 0.0
                cn = 0.0
            else:
                sn = (eb - 1.0 / eb) / 2.0
                cn = (eb + 1.0 / eb) / 2.0
            # avoid division by zero for very small sn
            if sn != 0.0:
                ar[i] = adel_i / (sn * bdn)
            else:
                ar[i] = 0.0
            # energy for effective temperature (dist)
            dist[i] = adel_i * bd_e * cn / (2.0 * sn) if sn != 0.0 else 0.0
            # sum of dist for tempf update
            if sn != 0.0:
                tsave += dist[i] / self.BK
            # contribution to Debye–Waller parameter
            dbw[i] = ar[i] * cn
            if self.dwpix[t_index] > 0.0:
                self.dwpix[t_index] += dbw[i]
        # base weight and temperature ratio for sint
        wt0 = self.tbeta
        tbart0 = 0.0
        # tbart0 = tempf / tempr for this temperature
        tempk = abs(self.inp.temperatures[t_index])
        if tempk != 0.0:
            tbart0 = self.tempf[t_index] / tempk
        # loop over α values
        # arat scaling for secondary scatterers (not implemented yet)
        arat = 1.0
        for ialpha in range(nalpha):
            # scaled α
            al = alpha[ialpha] * sc / arat
            # Debye–Waller factor for discrete lines
            dwf = math.exp(-al * self.dwpix[t_index]) if self.dwpix[t_index] != 0.0 else 1.0
            # extend the continuous scattering law to ±β for this α
            sex = self._exts(self.ssm[t_index, :, ialpha].copy(), exb, betan)
            # initialise array for the new scattering values (negative β part)
            sexpb = np.zeros(nbeta, dtype=float)
            # initialise delta functions: one line at zero with weight 1
            ben = [0.0]
            wtn = [1.0]
            nn = 1
            # local copies of wt and tbart that accumulate oscillator effects
            wt = wt0
            tbart = tbart0
            # loop over each oscillator
            for iosc in range(nd):
                # compute Bessel function factors for this oscillator
                dwc = al * dbw[iosc]
                x = al * ar[iosc]
                bzero, bplus, bminus = self._bfact(x, dwc, bdeln[iosc])
                # accumulate new discrete lines in temporary lists
                new_bes: List[float] = []
                new_wts: List[float] = []
                # n=0 term
                for m in range(nn):
                    besn = ben[m]
                    wtsn = wtn[m] * bzero
                    # include if line at negative beta or if weight is significant
                    if (besn <= 0.0) or (wtsn >= small):
                        if len(new_bes) < 500:
                            new_bes.append(besn)
                            new_wts.append(wtsn)
                # negative n terms (k > 0)
                for k in range(1, 51):
                    if bminus[k - 1] <= 0.0:
                        break
                    for m in range(nn):
                        besn = ben[m] - k * bdeln[iosc]
                        wtsn = wtn[m] * bminus[k - 1]
                        if (wtsn >= small) and (len(new_bes) < 500):
                            new_bes.append(besn)
                            new_wts.append(wtsn)
                # positive n terms (k > 0)
                for k in range(1, 51):
                    if bplus[k - 1] <= 0.0:
                        break
                    for m in range(nn):
                        besn = ben[m] + k * bdeln[iosc]
                        wtsn = wtn[m] * bplus[k - 1]
                        if (wtsn >= small) and (len(new_bes) < 500):
                            new_bes.append(besn)
                            new_wts.append(wtsn)
                # update the list of discrete lines
                ben = new_bes
                wtn = new_wts
                nn = len(ben)
                # update weight and effective temperature ratio for this oscillator
                wt += adel[iosc]
                # update tbart: dist is in eV; divide by (kB * T)
                if tempk != 0.0:
                    tbart += dist[iosc] / (self.BK * tempk)
            # after looping over oscillators, sort and trim discrete lines
            # if there are no lines, continue
            if nn == 0:
                # assign zeros and continue
                self.ssm[t_index, :, ialpha] = 0.0
                continue
            # pair the lines and weights and sort by decreasing weight
            pairs = list(zip(ben, wtn))
            pairs.sort(key=lambda x: x[1], reverse=True)
            # determine how many lines to keep
            n_keep = len(pairs)
            for idx in range(1, len(pairs)):
                # Fortran uses threshold 100*small beyond the first 5 lines
                if pairs[idx][1] < 100.0 * small and idx > 4:
                    n_keep = idx
                    break
            # truncate lists to n_keep
            pairs = pairs[:n_keep]
            # accumulate continuum contribution from discrete lines
            for besn, wtsn in pairs:
                weight_m = wtsn
                beta_shift = besn
                for j in range(nbeta):
                    be_val = -betan[j] - beta_shift
                    st = self._sint(be_val, bex, rdbex, sex, al, self.tbeta + self.twt, tbart, betan)
                    add = weight_m * st
                    if add >= tiny:
                        sexpb[j] += add
            # add delta function contributions when diffusion is not present
            if self.twt <= 0.0:
                m_idx = 0
                while m_idx < len(pairs):
                    besn, wtsn = pairs[m_idx]
                    m_idx += 1
                    # skip if Debye–Waller factor is vanishingly small
                    if dwf < vsmall:
                        break
                    # only negative beta shifts contribute to elastic line
                    if besn < 0.0:
                        be = -besn
                        # ignore shifts beyond the largest tabulated beta (excluding last point)
                        if be <= betan[-2]:
                            # find index jj of betan closest to be
                            idx = int(np.argmin(np.abs(betan - be)))
                            # compute weight to add
                            if idx <= 1:
                                denom = betan[idx] if betan[idx] != 0.0 else 0.0
                                add = wtsn / denom if denom != 0.0 else 0.0
                            else:
                                denom = betan[idx] - betan[idx - 2]
                                add = 2.0 * wtsn / denom if denom != 0.0 else 0.0
                            add *= dwf
                            if add >= tiny:
                                upd_idx = idx - 1
                                if upd_idx >= 0:
                                    sexpb[upd_idx] += add
            # store the calculated negative‑beta scattering law in ssm
            self.ssm[t_index, :, ialpha] = sexpb
        # update the effective temperature for this temperature index
        self.tempf[t_index] = (self.tbeta + self.twt) * self.tempf[t_index] + tsave

    # ------------------------------------------------------------------
    # Pair‑correlation correction (Sköld approximation)
    def _skold(self, t_index: int, temp: float, ska: np.ndarray, dka: float, cfrac: float, alpha: np.ndarray, beta: np.ndarray) -> None:
        """
        Apply the Sköld approximation to account for intermolecular coherence.

        When pair correlations are enabled with ``nsk = 2`` in the input deck,
        this routine modifies the scattering law to include the effect of
        coherent scattering.  The algorithm follows the Fortran subroutine
        ``skold`` from LEAPR.

        Parameters
        ----------
        t_index : int
            Temperature index.
        temp : float
            The physical temperature in Kelvin.
        ska : ndarray
            Tabulated S(kappa) values of length ``nka``.
        dka : float
            Increment between successive kappa values (inverse Angstroms).
        cfrac : float
            Coherent fraction (between 0 and 1).
        alpha : ndarray
            Array of α values (dimensionless).
        beta : ndarray
            Array of β values (dimensionless prior to scaling).
        """
        nalpha = self.inp.nalpha
        nbeta = self.inp.nbeta
        # Return early if no coherent fraction or no ska values
        if cfrac is None or len(ska) == 0 or cfrac == 0.0:
            return
        # compute scaled beta and alpha values
        tev = self.BK * abs(temp)
        sc = 1.0
        if self.inp.lat == 1:
            # in LEAPR: therm / tev with therm = 0.0253 eV
            sc = 0.0253 / tev
        # mass of scattering atom in grams (awr is weight ratio to neutron)
        # Fortran: amass = awr * amassn * amu
        amass = self.inp.awr * self.AMASSN * self.AMU
        # prepare a temporary array for coherent scattering values
        scoh = np.zeros(nalpha, dtype=float)
        # loop over β points
        for i_beta in range(nbeta):
            # loop over α values
            for j_alpha in range(nalpha):
                # scaled α value
                al = alpha[j_alpha] * sc / 1.0  # arat = 1 for principal scatterer
                # compute wavenumber kappa (1/Angstrom)
                # waven = ANGST * sqrt(2 * amass * tev * ev * al) / hbar
                if al <= 0.0:
                    sk_val = 1.0
                    ap = alpha[j_alpha]
                else:
                    waven = self.ANGST * math.sqrt(2.0 * amass * tev * self.EV * al) / self.HBAR
                    # obtain S(kappa) by interpolation
                    sk_val = self._terpk(ska, dka, waven)
                    # compute effective alpha (ap) = alpha / S(kappa)
                    if sk_val != 0.0:
                        ap = alpha[j_alpha] / sk_val
                    else:
                        ap = alpha[j_alpha]
                # find index kk such that alpha[kk] >= ap
                # if ap is larger than all alpha values, clamp kk = nalpha - 1
                kk = nalpha - 1
                for k_idx in range(nalpha):
                    if ap < alpha[k_idx]:
                        kk = k_idx
                        break
                # ensure kk >= 1 for interpolation (0‑based)
                if kk <= 0:
                    kk = 1
                # perform interpolation in log‑log space (interp_code = 5)
                y = 0.0
                y_lower = self.ssm[t_index, i_beta, kk - 1]
                y_upper = self.ssm[t_index, i_beta, kk] if kk < nalpha else 0.0
                x_lower = alpha[kk - 1]
                x_upper = alpha[kk] if kk < nalpha else x_lower
                if y_lower > 0.0 and y_upper > 0.0 and x_upper != x_lower:
                    y = self._terp1(x_lower, y_lower, x_upper, y_upper, ap, 5)
                scoh[j_alpha] = y * sk_val
            # after computing scoh for this β, update ssm by mixing with coherent part
            for j_alpha in range(nalpha):
                incoh = self.ssm[t_index, i_beta, j_alpha]
                coh = scoh[j_alpha]
                self.ssm[t_index, i_beta, j_alpha] = (1.0 - cfrac) * incoh + cfrac * coh

    # ------------------------------------------------------------------
    # Cold hydrogen/deuterium corrections
    def _coldh(self, t_index: int, temp: float, ska: np.ndarray, dka: float, alpha: np.ndarray, beta: np.ndarray) -> None:
        """
        Convolve the current S(α,β) with discrete rotational modes for
        cold hydrogen or deuterium.

        This routine mirrors the Fortran subroutine ``coldh``.  It
        computes a new scattering law that is no longer symmetric in
        β.  The negative‑β branch is stored in ``self.ssm`` and the
        positive‑β branch in ``self.ssp``.  The calculation relies
        on several helper functions (`_bt`, `_sumh`, `_cn`, `_sjbes`,
        `_bfill`, `_exts`, and `_sint`) to handle Bessel functions,
        Clebsch–Gordon coefficients and interpolation.

        Parameters
        ----------
        t_index : int
            Temperature index to process.
        temp : float
            Physical temperature in Kelvin (must match ``self.inp.temperatures[t_index]``).
        ska : ndarray
            Tabulated S(κ) values for the pair‑correlation/cold‑hydrogen option.
        dka : float
            Increment between successive κ values (inverse Angstroms).
        alpha : ndarray
            Array of α values (dimensionless).
        beta : ndarray
            Array of β values (dimensionless prior to scaling).
        """
        nbeta = self.inp.nbeta
        nalpha = self.inp.nalpha
        lat = self.inp.lat
        # constants specific to cold hydrogen/deuterium
        # masses of proton and deuteron in grams
        pmass = 1.6726231e-24
        dmass = 3.343586e-24
        # dissociation energies (eV) for H2 and D2
        deh = 0.0147
        ded = 0.0074
        # correlation factors for hydrogen and deuterium
        sampch = 0.356
        sampcd = 0.668
        sampih = 2.526
        sampid = 0.403
        # small threshold for printing (not used here)
        small = 1.0e-6
        # thermal energy for lattice case (eV)
        therm = 0.0253
        # compute effective temperature (in eV) and scaling factor
        tev = self.BK * abs(temp)
        sc = 1.0
        if lat == 1:
            sc = therm / tev
        # determine which law to apply: law = ncold + 1
        law = self.inp.ncold + 1
        # select dissociation energy and masses
        if law > 3:
            # deuterium
            de = ded
            amassm = 6.69e-24  # mass of D2 molecule in grams
            sampc = sampcd
            sampi = sampid
            # parameter bp for D2: hbar/2 * sqrt(2/(de*ev*dmass)) / angst
            bp = (self.HBAR / 2.0) * math.sqrt(2.0 / (de * self.EV * dmass)) / self.ANGST
        else:
            # hydrogen
            de = deh
            amassm = 3.3464e-24  # mass of H2 molecule in grams
            sampc = sampch
            sampi = sampih
            # parameter bp for H2: hbar/2 * sqrt(2/(de*ev*pmass)) / angst
            bp = (self.HBAR / 2.0) * math.sqrt(2.0 / (de * self.EV * pmass)) / self.ANGST
        # x parameter used in Boltzmann factors
        x = de / tev if tev != 0.0 else 0.0
        # total weight factor for sint
        wt = self.twt + self.tbeta
        # effective temperature ratio tbart = tempf / tempr
        tempk = abs(self.inp.temperatures[t_index])
        tbart = 0.0
        if tempk != 0.0:
            # avoid division by zero
            tbart = self.tempf[t_index] / tempk
        # ensure ssp is the same shape as ssm
        if self.ssp.shape != self.ssm.shape:
            self.ssp = np.zeros_like(self.ssm)
        # allocate workspace arrays for β and interpolation
        betan = np.zeros(nbeta, dtype=float)
        exb = np.zeros(nbeta, dtype=float)
        # extended β grid and reciprocal differences will be created in loop
        # prepare betan and exb only once (for nal=0)
        first_alpha = True
        # loop over α values (nal index)
        for nal in range(nalpha):
            # scaled α value
            # arat=1 for principal scatterer (no secondary scatterer mixing implemented)
            arat = 1.0
            al = alpha[nal] * sc / arat
            alp = wt * al
            # compute wavenumber κ in 1/Angstrom
            if al > 0.0:
                waven = self.ANGST * math.sqrt(amassm * tev * self.EV * al) / self.HBAR
            else:
                waven = 0.0
            # dimensionless y parameter for Bessel functions
            y = bp * waven
            # interpolate static structure factor S(κ)
            sk_val = 1.0
            if len(ska) > 0:
                sk_val = self._terpk(ska, dka, waven)
            # compute spin‑correlation weights swe (even jp) and swo (odd jp)
            swe = 0.0
            swo = 0.0
            if law == 2:
                swe = sampi * sampi / 3.0
                swo = sk_val * sampc * sampc + 2.0 * sampi * sampi / 3.0
            elif law == 3:
                swe = sk_val * sampc * sampc
                swo = sampi * sampi
            elif law == 4:
                swe = sk_val * sampc * sampc + 5.0 * sampi * sampi / 8.0
                swo = 3.0 * sampi * sampi / 8.0
            elif law == 5:
                swe = 3.0 * sampi * sampi / 4.0
                swo = sk_val * sampc * sampc + sampi * sampi / 4.0
            # normalise the spin factors by the sum of squares
            snorm = sampi * sampi + sampc * sampc
            if snorm != 0.0:
                swe /= snorm
                swo /= snorm
            # prepare betan and exb on first α iteration
            if first_alpha:
                for i in range(nbeta):
                    be = beta[i]
                    # scale β for lattice case
                    if lat == 1:
                        be = be * therm / tev
                    betan[i] = be
                    exb[i] = math.exp(-be / 2.0)
                # build extended β grid and reciprocal differences
                bex, rdbex = self._bfill(betan)
                first_alpha = False
            # extend current scattering law to both ±β for this α
            sex = self._exts(self.ssm[t_index, :, nal].copy(), exb, betan)
            # total number of points in extended β grid
            nbx = len(bex)
            # compute maximum β value for SCT approximation
            # not needed explicitly since _sint handles out‑of‑range values
            # main loop over jj indexes corresponding to negative and positive β
            jjmax = 2 * nbeta - 1
            for jj in range(jjmax):
                # determine k index and sign of β
                if jj < nbeta - 1:
                    # negative β branch
                    k = nbeta - 1 - jj
                    be = -betan[k]
                    is_negative = True
                else:
                    # positive β branch
                    k = jj - (nbeta - 1)
                    be = betan[k]
                    is_negative = False
                # initialize sum for this β
                sn_val = 0.0
                # determine range of j values
                jterm = 3
                # starting parity for j depends on law
                ipo = 1
                if law == 2 or law == 5:
                    ipo = 2
                jt1 = 2 * jterm
                if ipo == 2:
                    jt1 += 1
                # loop over j values
                for l in range(ipo, jt1 + 1, 2):
                    j = l - 1
                    pj = self._bt(j, x)
                    # even jp contributions (lp=1,3,5,7,9 => jp=0,2,4,6,8)
                    snlg = 0.0
                    for lp in range(1, 11, 2):
                        jp = lp - 1
                        betap = (-j * (j + 1) + jp * (jp + 1)) * x / 2.0
                        tmp = (2 * jp + 1) * pj * swe * 4.0 * self._sumh(j, jp, y)
                        bn = be + betap
                        # compute add via interpolation (ifree=0)
                        add = self._sint(bn, bex, rdbex, sex, al, wt, tbart, betan)
                        snlg += tmp * add
                    # odd jp contributions (lp=2,4,6,8,10 => jp=1,3,5,7,9)
                    snlk = 0.0
                    for lp in range(2, 11, 2):
                        jp = lp - 1
                        betap = (-j * (j + 1) + jp * (jp + 1)) * x / 2.0
                        tmp = (2 * jp + 1) * pj * swo * 4.0 * self._sumh(j, jp, y)
                        bn = be + betap
                        add = self._sint(bn, bex, rdbex, sex, al, wt, tbart, betan)
                        snlk += tmp * add
                    # accumulate j contributions
                    sn_val += snlg + snlk
                # assign to appropriate array branch
                if is_negative:
                    # negative β into ssm
                    self.ssm[t_index, k, nal] = sn_val
                else:
                    # positive β into ssp
                    self.ssp[t_index, k, nal] = sn_val

    def _bt(self, j: int, x: float) -> float:
        """
        Compute the statistical weight factor pj for cold hydrogen/deuterium.

        This mirrors the Fortran subroutine ``bt``.  The integer
        ``j`` corresponds to the rotational quantum number and ``x``
        is related to the dissociation energy and temperature.  The
        result gives the probability of occupying level ``j``.

        Parameters
        ----------
        j : int
            Rotational quantum number (0 ≤ j).
        x : float
            Parameter ``x = de / tev`` from the parent routine.

        Returns
        -------
        float
            Statistical weight factor ``pj`` for state ``j``.
        """
        half = 0.5
        # numerator for j state
        yy = half * j * (j + 1)
        a = (2 * j + 1) * math.exp(-yy * x)
        # denominator sums over k = even or odd depending on j parity
        b = 0.0
        for i in range(1, 11):
            k = 2 * i - 2
            if j % 2 == 1:
                # shift to odd k for odd j
                k += 1
            yy = half * k * (k + 1)
            b += (2 * k + 1) * math.exp(-yy * x)
        # avoid division by zero
        if b == 0.0:
            return 0.0
        return a / (2.0 * b)

    def _sumh(self, j: int, jp: int, y: float) -> float:
        """
        Sum over Bessel functions and Clebsch–Gordon coefficients for cold H/D.

        Implements the Fortran function ``sumh``.  For given rotational
        quantum numbers ``j`` and ``jp``, and argument ``y``, this
        function computes a sum over spherical Bessel functions and
        Clebsch–Gordon coefficients squared.  The range of summation
        depends on the difference ``|j - jp|`` and the maximum
        allowed by the implementation (up to nine terms).

        Parameters
        ----------
        j : int
            Rotational quantum number j (≥ 0).
        jp : int
            Rotational quantum number j' (≥ 0).
        y : float
            Argument for the Bessel functions.

        Returns
        -------
        float
            The summed contribution ``sumh(j,jp,y)``.
        """
        # handle cases where one of the states is zero
        if j == 0:
            return (self._sjbes(jp, y) * self._cn(j, jp, jp)) ** 2
        if jp == 0:
            return (self._sjbes(j, y) * self._cn(j, 0, j)) ** 2
        # general case: sum over n = |j - jp| + 1 .. j + jp + 1 (but at most 9 terms)
        sum_val = 0.0
        imk = abs(j - jp) + 1
        ipk1 = j + jp + 1
        # number of terms mpk = ipk1 - imk; ensure at most nine terms
        mpk = ipk1 - imk
        if mpk <= 9:
            ipk = ipk1
        else:
            ipk = imk + 9
        for n in range(imk, ipk + 1):
            n1 = n - 1
            bval = self._sjbes(n1, y)
            cval = self._cn(j, jp, n1)
            sum_val += (bval * cval) ** 2
        return sum_val

    def _cn(self, jj: int, ll: int, nn: int) -> float:
        """
        Compute Clebsch–Gordon coefficients for cold H/D calculations.

        This replicates the Fortran function ``cn``.  It evaluates
        specific Clebsch–Gordon coefficients needed in the cold
        hydrogen/deuterium treatment.  The formula involves
        factorials and square roots expressed via logarithms to
        maintain numerical stability.

        Parameters
        ----------
        jj : int
            First angular momentum quantum number.
        ll : int
            Second angular momentum quantum number.
        nn : int
            Resulting angular momentum quantum number.

        Returns
        -------
        float
            The Clebsch–Gordon coefficient ``cn(jj,ll,nn)``.
        """
        # determine if the triangle inequality is satisfied
        kdet = (jj + ll + nn) // 2
        kdel = jj + ll + nn - 2 * kdet
        if kdel != 0:
            return 0.0
        # compute factorial‐like products via sums of logs
        ka1 = jj + ll + nn
        ka2 = jj + ll - nn
        ka3 = jj - ll + nn
        ka4 = ll - jj + nn
        kb1 = ka1 // 2
        kb2 = ka2 // 2
        kb3 = ka3 // 2
        kb4 = ka4 // 2
        # helper to compute factorial terms using logs
        def log_fact(n: int) -> float:
            s = 0.0
            for i in range(1, n + 1):
                s += math.log(float(i))
            return s
        # compute a1..a4
        a1 = math.sqrt(math.exp(log_fact(ka1))) if ka1 > 0 else 1.0
        a2 = math.sqrt(math.exp(log_fact(ka2))) if ka2 > 0 else 1.0
        a3 = math.sqrt(math.exp(log_fact(ka3))) if ka3 > 0 else 1.0
        a4 = math.sqrt(math.exp(log_fact(ka4))) if ka4 > 0 else 1.0
        # compute b1..b4
        b1 = math.exp(log_fact(kb1)) if kb1 > 0 else 1.0
        b2 = math.exp(log_fact(kb2)) if kb2 > 0 else 1.0
        b3 = math.exp(log_fact(kb3)) if kb3 > 0 else 1.0
        b4 = math.exp(log_fact(kb4)) if kb4 > 0 else 1.0
        # compute Wigner symbol
        # ratio factor
        rat = (2.0 * nn + 1.0) / (jj + ll + nn + 1.0)
        # alternating sign factor
        iwign = (jj + ll - nn) // 2
        # wign may oscillate in sign
        wign = (-1) ** iwign
        # final Clebsch–Gordon coefficient
        # avoid division by zero in denominators
        if a1 == 0.0 or b1 == 0.0 or b2 == 0.0 or b3 == 0.0 or b4 == 0.0:
            return 0.0
        wign = wign * math.sqrt(rat) * (b1 / a1) * (a2 / b2) * (a3 / b3) * (a4 / b4)
        return wign

    def _sjbes(self, n: int, x: float) -> float:
        """
        Spherical Bessel function used in cold H/D calculations.

        This implements the Fortran function ``sjbes``.  It handles
        various ranges of the argument ``x`` using different
        approximations to maintain numerical stability.  Values beyond
        certain limits are deemed inaccurate and result in a return of
        zero.

        Parameters
        ----------
        n : int
            Order of the spherical Bessel function.
        x : float
            Argument of the function (≥ 0).

        Returns
        -------
        float
            The value of the spherical Bessel function j_n(x).
        """
        # parameter thresholds (mirroring Fortran constants)
        break1 = 3.0e4
        break2 = 7.0e-4
        break3 = 0.2
        huge = 1.0e25
        small = 2.0e-38
        # check for extreme orders or large arguments
        if n >= 30000 or x > break1:
            # outside range of accuracy
            return 0.0
        # invalid arguments
        if x < 0.0 or n < 0:
            return 0.0
        # compute normal values
        if x <= break2:
            # small x expansion
            if n == 0:
                return 1.0
            elif n > 10:
                return 0.0
            # series expansion for small x
            t1 = 3.0
            t2 = 1.0
            t3 = 0.0
            for i in range(1, n + 1):
                t3 = t2 * x / t1
                t1 += 2.0
                t2 = t3
            return t3
        else:
            # moderate to large x
            if x < break3:
                y = x * x
                w = 1.0 - y * (1.0 - y / 20.0) / 6.0
            else:
                # sine approximation
                w = math.sin(x) / x if x != 0.0 else 1.0
            if n == 0:
                return w
            # determine number of downward iterations for continued fraction
            if x >= 100.0:
                l = int(x / 50.0 + 18.0)
            elif x >= 10.0:
                l = int(x / 10.0 + 10.0)
            elif x > 1.0:
                l = int(x / 2.0 + 5.0)
            else:
                l = 5
            iii = int(x)
            kmax = n
            if iii > n:
                kmax = iii
            nm = kmax + l
            z = 1.0 / x if x != 0.0 else 0.0
            t3 = 0.0
            t2 = small
            sj = 0.0
            for _ in range(nm):
                k = nm - 1 - _
                t1 = (2.0 * k + 3.0) * z * t2 - t3
                if n == k:
                    sj = t1
                # rescale to avoid overflow
                if abs(t1) >= huge:
                    t1 /= huge
                    t2 /= huge
                    sj /= huge
                t3 = t2
                t2 = t1
            # t1 holds the last value computed
            # avoid division by zero
            if t1 == 0.0:
                return 0.0
            return w * sj / t1


    def _convol(self, t1: np.ndarray, tlast: np.ndarray, delta: float) -> np.ndarray:
        """
        Convolve t₁ with t_last to obtain t_next.

        This replicates the Fortran subroutine ``convol`` which
        computes the next term in the phonon expansion.  The
        convolution integral is evaluated on the discretised β grid
        using a trapezoidal rule.  The exponential factors in the
        original code are handled explicitly here.

        Parameters
        ----------
        t1 : ndarray
            Array of t₁(β) values.
        tlast : ndarray
            Array of tₙ₋₁(β) values.
        delta : float
            Grid spacing in β.

        Returns
        -------
        ndarray
            Array of tₙ(β) values.
        """
        n1 = len(t1)
        nl = len(tlast)
        # The length of the convolution result is (n1 + nl - 1)
        nn = n1 + nl - 1
        tnext = np.zeros(nn, dtype=float)
        tiny = 1.0e-30
        # perform discrete convolution with exponential weighting
        for k in range(nn):
            sum_val = 0.0
            for j in range(n1):
                i1 = k + j
                i2 = k - j
                f1 = 0.0
                be = j * delta
                if t1[j] > 0:
                    # first term: t_last(i1) * exp(-β)
                    if 0 <= i1 < nl:
                        f1 = tlast[i1] * math.exp(-be)
                    # second term: t_last(i2) without exponential if i2 >= 0
                    f2 = 0.0
                    if 0 <= i2 < nl:
                        f2 = tlast[i2]
                    elif i2 < 0 and (-i2) < nl:
                        be2 = -i2 * delta
                        f2 = tlast[-i2] * math.exp(-be2)
                    cc = t1[j] * (f1 + f2)
                    # trapezoidal end corrections for j=0 or j=n1-1
                    if j == 0 or j == n1 - 1:
                        cc *= 0.5
                    sum_val += cc
            tnext[k] = sum_val * delta
            # integrate for normalisation check (not returned here)
        # threshold small values
        tnext[tnext < tiny] = 0.0
        return tnext
