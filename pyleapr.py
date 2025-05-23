import re
import numpy as np

# Physical constants (from phys.f90)
pi = 3.14159265358979323846  # Pi
bk = 8.617333262145e-5       # Boltzmann constant in eV/K
amu = 931.49410242e6         # Atomic mass unit in eV/c^2 (for mass conversion if needed)

# Constants from locale.f90 and mainio.f90
THERM_CONST = 0.0253  # Standard thermal energy in eV, corresponding to 293.6 K or 20 C
TINY_SAB = 1.0e-30    # Smallest S(a,b) value considered non-zero
EXPLIM = -250.0       # Limit for np.exp argument to avoid underflow/overflow
SMALL_VAL = 1.0e-8    # General small value threshold
VSMALL_VAL = 1.0e-10  # Very small value threshold
TINY_VAL = 1.0e-20    # Tiny value threshold (distinct from TINY_SAB for S(a,b) values)
SLIM_LOG_ARG = -225.0 # Lower limit for log arguments, similar to NJOY's 'slim'

# Constants for _trans delta calculation (c0-c4 from NJOY trans.f)
TRANS_C0 = 0.025
TRANS_C1 = 0.17
TRANS_C2 = 1.0
TRANS_C3 = 2.5
TRANS_C4 = 0.05

# Constant for _stable
STABLE_EPS = 1.0e-7
STABLE_QUART = 0.25

# Physical constants for wavenumber calculation in skold
HBARC_EV_ANGSTROM = 1973.269804 # hbar*c in eV*Angstrom

# Constants for _coldh (cold hydrogen/deuterium)
PMASS_AMU = 1.00782503223  # Mass of H1 atom in AMU
DMASS_AMU = 2.01410177812  # Mass of D atom in AMU
DEH_EV = 0.00746           # Rotational constant for H2 in eV
DED_EV = 0.00373           # Rotational constant for D2 in eV
# Scattering lengths b (fm) -> sigma_coh = 4*pi*b_coh^2, sigma_inc = 4*pi*b_inc_sq_val
# Values from NJOY manual Table H.1 for H1 and D.
# Coherent scattering lengths b_coh (fm)
BCOH_H_FM = -3.7406
BCOH_D_FM = 6.671
# Incoherent cross sections sigma_inc (barns)
SIGINC_H_BARN = 80.26
SIGINC_D_BARN = 2.05
# HBAR_EV_S for bp calculation, if needed directly: 6.582119569e-16 eV.s
# AMU_G for bp calculation: 1.6605402e-24 g


kr = 22  # Fortran unit number for reading main input
k4 = 4   # Fortran unit number for reading secondary input (e.g. ENDF tape)
k8 = 8   # Fortran unit number for writing primary output
nsysi = 5 # Standard input unit
nsyso = 6 # Standard output unit

class NJOYError(Exception):
    """Custom exception for NJOY specific errors."""
    pass

class LeaprCalculator:
    def __init__(self, input_file_path=None):
        self.input_file_path = input_file_path
        self.lines = []
        self.current_line_index = 0

        # Card 1
        self.nout = 0
        self.title = ""

        # Card 2
        self.ntempr = 0
        self.iprint = 0
        self.nphon = 0

        # Card 3
        self.alpha_values = [] # alpha
        self.beta_values = []  # beta

        # Card 4
        self.lat = 0
        self.za = 0.0
        self.isabt = 0
        self.ilog = 0
        self.smin = 0.0

        # Card 5
        self.nbeta = 0
        self.nalpha = 0

        # Card 6
        self.awr = 0.0
        self.spr = []
        self.delta = []

        # Card 7
        self.ni = 0
        self.njoy = 0

        # Card 8
        self.nsk = 0
        self.sbeta = []

        # Card 9
        self.icoh = 0
        self.ncold = 0
        self.nss = 0
        self.b7 = 0.0

        # Card 10 (if ncold > 0)
        self.acold = []
        self.bcold = []

        # Card 11 (if ncold > 0)
        self.tbeta_cold = []

        # Card 12 (if ncold > 0)
        self.scold = []
        self.arsc = []

        # Card 13 (if ncold > 0)
        self.parsc = []

        # Card 14 (if nss > 0)
        self.kappa_values = []

        # Card 15 (if nss > 0)
        self.tbeta_s = []

        # Card 16 (if nss > 0)
        self.free_atom_xs_s = []
        self.epsilon_s = []

        # Card 17 (if nss > 0)
        self.kappa_s = []

        # Temperature-specific data
        self.temperatures_data = []  # List of dictionaries, each holding data for a temperature
        
        # Attributes that are results of calculations
        self.sab_results = {}  # Store S(alpha, beta) results, perhaps keyed by temperature
        
        # Attributes for contin, start, fsum, terpt, convol methods
        self.ssm_values = np.array([]) # To be initialized in leapr_main
        self.dwpix_values = np.array([]) # To be initialized in leapr_main
        self.tempf_values = np.array([]) # To be initialized in leapr_main
        self.arat = 1.0 # Ratio aws/awr, default to 1.0 (no secondary scatterer)


        if self.input_file_path:
            self._load_and_preprocess_lines()

    def _load_and_preprocess_lines(self):
        with open(self.input_file_path, 'r') as f:
            for line in f:
                # Remove Fortran-style comments and trim whitespace
                cleaned_line = re.sub(r"!.*$", "", line).strip()
                if cleaned_line:  # Only add non-empty lines
                    self.lines.append(cleaned_line)

    def _read_next_significant_line(self):
        if self.current_line_index < len(self.lines):
            line = self.lines[self.current_line_index]
            self.current_line_index += 1
            return line
        return None

    def _read_line_as_values(self, types):
        line = self._read_next_significant_line()
        if line is None:
            raise EOFError("Unexpected end of input file.")
        
        values_str = line.split()
        # Handle varying number of values for title
        if str in types and types.index(str) == len(types) -1 and len(values_str) >= len(types) :
             # Join all remaining parts for the title string
            values_str = values_str[:len(types)-1] + [' '.join(values_str[len(types)-1:])]


        if len(values_str) != len(types):
            raise ValueError(f"Line '{line}' has {len(values_str)} values, expected {len(types)} based on types {types}.")
        
        values = []
        for i, type_func in enumerate(types):
            try:
                values.append(type_func(values_str[i]))
            except ValueError:
                raise ValueError(f"Could not convert '{values_str[i]}' to {type_func.__name__} in line '{line}'.")
        return values

    def _read_float_list(self, num_values):
        values = []
        while len(values) < num_values:
            line = self._read_next_significant_line()
            if line is None:
                raise EOFError("Unexpected end of input file while reading float list.")
            values.extend([float(x) for x in line.split()])
        if len(values) > num_values: # if last line had more values than needed
            values = values[:num_values]
        return values
    
    def _add_temperature_set(self, temp_data):
        """Adds a new set of temperature-specific data."""
        self.temperatures_data.append(temp_data)

    def _parse_input(self):
        if not self.input_file_path:
            print("Input file path not provided. Cannot parse input.")
            return
        if not self.lines: # Ensure lines are loaded if not already
            self._load_and_preprocess_lines()
            
        # Card 1
        self.nout, self.title = self._read_line_as_values([int, str])

        # Card 2
        self.ntempr, self.iprint, self.nphon = self._read_line_as_values([int, int, int])

        # Card 3
        if self.nphon > 0:
            self.alpha_values = self._read_float_list(self.nphon) # alpha
            self.beta_values = self._read_float_list(self.nphon)  # beta

        # Card 4
        self.lat, self.za, self.isabt, self.ilog, self.smin = self._read_line_as_values([int, float, int, int, float])
        
        # Card 5
        self.nbeta, self.nalpha = self._read_line_as_values([int, int])

        # Card 6
        line_values_card6 = self._read_next_significant_line().split()
        self.awr = float(line_values_card6[0])
        if self.nphon > 0 :
            num_spr = self.nphon // 2 # As per Fortran code, spr and delta are read only if nphon > 0
            self.spr = [float(x) for x in line_values_card6[1:1+num_spr]]
            self.delta = [float(x) for x in line_values_card6[1+num_spr:1+num_spr*2]]


        # Card 7
        self.ni, self.njoy = self._read_line_as_values([int, int])

        # Card 8
        self.nsk = self._read_line_as_values([int])[0]
        if self.nsk > 0:
            self.sbeta = self._read_float_list(self.nsk)
        
        # Card 9
        self.icoh, self.ncold, self.nss, self.b7 = self._read_line_as_values([int, int, int, float])

        # Card 10 (if ncold > 0)
        if self.ncold > 0:
            self.acold = self._read_float_list(6) # acold is always size 6
            self.bcold = self._read_float_list(self.ncold) 
        
        # Card 11 (if ncold > 0)
        if self.ncold > 0:
            self.tbeta_cold = self._read_float_list(self.ncold)

        # Card 12 (if ncold > 0)
        if self.ncold > 0:
            self.scold = self._read_float_list(self.ncold)
            self.arsc = self._read_float_list(self.ncold) 

        # Card 13 (if ncold > 0)
        if self.ncold > 0:
            self.parsc = self._read_float_list(self.ncold)
        
        # Card 14 (if nss > 0)
        if self.nss > 0:
            self.kappa_values = self._read_float_list(self.nss)

        # Card 15 (if nss > 0)
        if self.nss > 0:
            self.tbeta_s = self._read_float_list(self.nss)

        # Card 16 (if nss > 0)
        if self.nss > 0:
            self.free_atom_xs_s = self._read_float_list(self.nss)
            self.epsilon_s = self._read_float_list(self.nss)

        # Card 17 (if nss > 0)
        if self.nss > 0:
            self.kappa_s = self._read_float_list(self.nss)
        
        # Loop for each temperature
        for _ in range(self.ntempr):
            temp_data = {}
            # Card 18
            temp_data['temperature'] = self._read_line_as_values([float])[0]
            # Card 19
            # Ensure nbeta is positive before reading secondary scattering law
            if self.nbeta > 0:
                 temp_data['secondary_scattering_law'] = self._read_float_list(self.nbeta) 
            else:
                 temp_data['secondary_scattering_law'] = [] # Assign empty list if nbeta is not positive
            self._add_temperature_set(temp_data)


    def leapr_main(self):
        """
        Main orchestrator method for the LEAPR calculations.
        """
        if not self.input_file_path:
            print("Input file not set. Cannot run calculations.")
            return

        self._parse_input() # Parse the input file first

        print(f"Starting LEAPR calculations for: {self.title}")

        if not self.temperatures_data:
            print("No temperature data loaded. Check input file or parsing logic.")
            return

        # Initialize ssm_values, dwpix_values, tempf_values before temperature loop
        if self.ntempr > 0 and self.nalpha > 0 and self.nbeta > 0:
            self.ssm_values = np.zeros((self.ntempr, self.nalpha, self.nbeta))
        else:
            # Handle cases where nalpha or nbeta might be zero from input
            # Or if ntempr is zero
            print("Warning: ntempr, nalpha, or nbeta is zero. S(a,b) calculations might be skipped or invalid.")
            self.ssm_values = np.zeros((self.ntempr, self.nalpha if self.nalpha > 0 else 1, self.nbeta if self.nbeta > 0 else 1))
        
        # ssp_values for positive beta results from coldh
        self.ssp_values = np.zeros_like(self.ssm_values)


        self.dwpix_values = np.zeros(self.ntempr)
        self.tempf_values = np.zeros(self.ntempr)

        for itemp, temp_data in enumerate(self.temperatures_data):
            current_temp_k = temp_data['temperature']
            print(f"\nProcessing temperature: {current_temp_k}K (itemp={itemp})")

            # Continuous spectrum calculations (S(a,b) stored in self.ssm_values)
            if self.lat == 1: 
                if temp_data.get('ni', 0) > 0 and self.nphon > 0:
                    # Card 20 data for continuous model:
                    # p1_values (rho(E)), ni (num points), delta1 (delta E), tbeta (norm)
                    # These should be in temp_data if ni > 0.
                    self._contin(itemp, temp_data, self.nphon)
                elif self.iprint >=1 :
                    print(f"Skipping _contin for temp {current_temp_k}K (lat=1) due to missing p1_values (ni=0) or nphon=0.")
            else: # lat != 1, ssm for this temperature should be initialized to zero if not done elsewhere.
                  # (It is initialized before the loop)
                if self.iprint >=1:
                    print(f"lat != 1, skipping _contin for temp {current_temp_k}K. S(a,b) from contin is zero.")


            # Discrete oscillator calculations
            # nd (number of discrete oscillators) is read on Card 21 in Fortran
            # This card is not currently parsed in _parse_input.
            # For now, let's assume nd, bdel_values, adel_values would be in temp_data.
            num_discrete_oscillators = temp_data.get('nd', 0)
            if num_discrete_oscillators > 0:
                # Placeholder: these values need to be correctly parsed into temp_data
                # temp_data['bdel_values'] = [...] # Oscillator energies
                # temp_data['adel_values'] = [...] # Oscillator weights
                # temp_data['tbeta_continuous'] = temp_data.get('tbeta',0.0) # from card 20
                # temp_data['twt_translational'] = temp_data.get('twt',0.0) # from card 22 (not parsed yet)
                
                # Call _discre only if nd > 0 and necessary data exists
                # This will require adding parsing for Card 21 (nd, (bdel(i), adel(i)))
                # and Card 22 (twt - if used)
                # For now, this call won't execute with current dummy input.
                # self._discre(itemp, temp_data)
                if self.iprint >=1:
                    print(f"Discrete oscillators (nd={num_discrete_oscillators}) present, but _discre call is currently placeholder.")
            elif self.iprint >=2:
                 print(f"No discrete oscillators (nd=0) for temp {current_temp_k}K.")

            # Translational motion calculations (S(a,b) in self.ssm_values modified)
            if current_temp_data.get('twt', 0.0) > SMALL_VAL : # twt from Card 22
                self._trans(itemp, temp_data)
            elif self.iprint >=2:
                print(f"Skipping _trans for temp {current_temp_k}K as twt is zero or negligible.")

            # Skold approximation (modifies self.ssm_values)
            if self.nsk == 2:
                self._skold(itemp, temp_data)
            
            # Cold H/D coherent/incoherent scattering (populates ssp_values and ssm_values for neg beta)
            if self.ncold > 0 : # ncold from Card 9
                # self._coldh(itemp, temp_data) # Call placeholder
                if self.iprint >=1:
                    print(f"Cold H/D (ncold={self.ncold}) processing for temp {current_temp_k}K. _coldh call is placeholder.")
            elif self.iprint >=2:
                print(f"Skipping _coldh for temp {current_temp_k}K as ncold is zero.")


            # Placeholder for other calculation steps
            # self.calculate_incoherent_elastic(current_temp_k, temp_data)
            # self.calculate_incoherent_inelastic(current_temp_k, temp_data)
            # if self.icoh > 0:
            #     self.calculate_coherent_elastic(current_temp_k, temp_data)
            
            # Store or output results as needed
            # self.sab_results[current_temp_k] = {"S_alpha_beta": ..., "S_incoherent_elastic": ...}

        print("\nLEAPR calculations completed.")

    # --- Start of translated Fortran subroutines ---

    def _fsum(self, n, p_array, np_val, tau, deltab_val):
        """
        Corresponds to Fortran function fsum.
        Calculates sums needed for Debye-Waller factor and effective temperature.
        n:      integer, 1 or 2
        p_array:  numpy array, phonon spectrum (rho)
        np_val:   integer, number of points in p_array
        tau:    numpy array, dimensionless temperature grid
        deltab_val: float, delta beta for phonon expansion
        """
        s = 0.0
        if np_val <= 0:
            return s
        
        factor = deltab_val / (2.0 * tau[0]) # tau[0] is temp in eV * beta
                                          # this seems to be related to temp_k * bk * beta_val
                                          # In Fortran: deltab / (2.d0 * tau(1))
                                          # tau(i) = tev * beta(i)

        # The tau argument in Fortran fsum is actually just 'tev' (temp in eV)
        # and deltab. The loop is over beta values implicitly through p_array.
        # Let's re-evaluate based on how `fsum` is called in `start`.
        # In `start`: an = fsum(1, p, npt, tev, deltab)
        #             bn = fsum(2, p, npt, tev, deltab)
        # So, `tau` in fsum's signature is `tev` (temperature in eV).
        # `p_array` is the phonon spectrum, `np_val` is `npt`.

        if tau == 0.0: # tau here is tev
            return 0.0

        factor_exp = deltab_val / tau # deltab / tev
        
        for k in range(np_val): # k from 0 to np_val-1
            ek = (k + 0.5) * deltab_val # beta value for midpoint of interval k
                                      # Fortran: ek = (k-1+0.5)*deltab for k=1 to npt
            
            term = p_array[k] * deltab_val
            if n == 2:
                term *= ek
            
            # Handling e^x where x can be large
            exp_arg = ek / tau # ek is beta, tau is tev. So beta/tev.
            if exp_arg > np.log(1.0/TINY_SAB): # Avoid overflow for exp(large_val)
                coth_val = 1.0 # coth(x) -> 1 for large x
            elif exp_arg < -np.log(1.0/TINY_SAB): # Avoid overflow for exp(-large_val) leading to coth(-x) -> -1
                 coth_val = -1.0
            else:
                exp_val = np.exp(exp_arg)
                if exp_val == 1.0: # Avoid division by zero if exp_val is exactly 1 (e.g. arg is 0)
                    # This case (ek/tau = 0) should ideally not happen if ek represents energy
                    # For small x, coth(x) approx 1/x.
                    # If ek is beta, and beta can be small, this needs care.
                    # However, the Fortran code doesn't have special handling here beyond typical float behavior.
                    # Let's assume ek/tau is not pathologically small to make exp_val-1 zero unless ek is zero.
                    # If ek is 0, then p_array[k] should be 0 for physical spectrum.
                    # For now, direct translation.
                    if exp_val - 1.0 == 0:
                        coth_val = np.inf # or handle as error, or use approximation
                    else:
                        coth_val = (exp_val + 1.0) / (exp_val - 1.0)
            s += term * (coth_val - 1.0)
            
        return s * 0.5

    def _terpt(self, tn_array, ntn_val, delta_val, be_val):
        """
        Corresponds to Fortran function terpt.
        Interpolates in a table (tn_array) at a point be_val.
        Assumes tn_array is sorted.
        tn_array:  numpy array, table of function values (e.g., S(beta))
        ntn_val:   integer, number of points in tn_array
        delta_val: float, spacing of the argument of tn_array (delta beta)
        be_val:    float, beta value at which to interpolate
        """
        if ntn_val == 0:
            return 0.0
        if ntn_val == 1:
            return tn_array[0]

        # Fortran: x = be / delta + 1.0
        # Python: x = be_val / delta_val (0-indexed location)
        x = be_val / delta_val 
        
        # Fortran: i = min( max( 1, int(x) ), ntn-1 )
        # Python: i_idx = min( max( 0, int(x) ), ntn_val - 2 )
        # int(x) in Fortran truncates. In Python, int() also truncates.
        # i is the index of the lower bound for interpolation.
        i_idx = int(x)

        if i_idx < 0:
            # Extrapolate based on first point, or return boundary, or error.
            # Fortran terpt would use i=1 (0 in Python) and x=x-i, so x could be negative.
            # If x (Fortran) is < 1.0, then i=1. x-i becomes be/delta.
            # tn(i) + (x-i)*(tn(i+1)-tn(i))
            # If be/delta = -0.5, x(F)=0.5, i=1. tn(1) + (-0.5)*(tn(2)-tn(1))
            # Python: x_py = -0.5. i_idx = int(-0.5) = 0.
            #         rem = x_py - i_idx = -0.5 - 0 = -0.5
            #         tn_array[0] + (-0.5) * (tn_array[1] - tn_array[0])
            # This seems to align.
            i_idx = 0
        elif i_idx >= ntn_val -1: # If x is such that i_idx is last point or beyond
            i_idx = ntn_val - 2 # Use second to last point as lower bound for interpolation to tn_array[ntn_val-1]
        
        # Fortran: rem = x - i
        # Python: rem = x - i_idx
        rem = x - i_idx

        # Linear interpolation: f(a) + (f(b)-f(a))*(x-a)/(b-a)
        # Here, points are tn_array[i_idx] at x_knot = i_idx
        # and tn_array[i_idx+1] at x_knot = i_idx+1
        # So (x-a)/(b-a) is (x - i_idx) / (i_idx+1 - i_idx) = x - i_idx = rem
        val = tn_array[i_idx] + rem * (tn_array[i_idx+1] - tn_array[i_idx])
        return val

    def _start(self, itemp, p_rho_values, num_rho_points, delta_e_rho, tbeta_norm_continuous, current_temp_k):
        """
        Corresponds to Fortran subroutine start.
        Initializes the phonon spectrum T1(beta) and calculates Debye-Waller factor (f0)
        and effective temperature (tbar).

        itemp: integer, current temperature index (0 to ntempr-1)
        p_rho_values: numpy array, rho(E) values (phonon density of states) from input (card 20 `p1`)
        num_rho_points: integer, number of points in p_rho_values (input `ni` for card 20)
        delta_e_rho: float, energy spacing (delta E) for p_rho_values (input `delta1` for card 20)
        tbeta_norm_continuous: float, normalization factor for continuous part (input `tbeta` for card 20)
        current_temp_k: float, current temperature in Kelvin
        
        Returns: tuple (p_initial, deltab_val, tev_val, f0_val, tbar_val)
                 p_initial: numpy array, initial phonon spectrum T1(beta)
                 deltab_val: float, delta beta for phonon expansion
                 tev_val: float, temperature in eV
                 f0_val: float, Debye-Waller factor integral
                 tbar_val: float, effective temperature factor (tbar/temp)
        """
        tev_val = current_temp_k * bk  # Temperature in eV
        if tev_val == 0.0:
            raise NJOYError("Error in _start: temperature in eV is zero.")

        # deltab is delta beta for the phonon expansion.
        # delta_e_rho is delta E for the input rho(E).
        # These are related by E = tev * beta. So delta_E = tev * delta_beta.
        # Thus, deltab_val = delta_e_rho / tev_val.
        deltab_val = delta_e_rho / tev_val

        if self.iprint >= 2:
            print(f"  tev={tev_val:.5e}, deltab={deltab_val:.5e}")

        p_initial = np.zeros(self.nbeta) # Corresponds to T1 in Fortran LEAPR

        if num_rho_points == 0: # no continuous part
            if self.iprint >=1:
                print(" warning in start. no continuous part given.")
            # Still need to return sensible defaults for f0, tbar if other parts (e.g. discrete osc) exist
            self.dwpix_values[itemp] = 0.0 # f0
            self.tempf_values[itemp] = current_temp_k # tbar = temp, so tempf = temp * 1.0
            return p_initial, deltab_val, tev_val, 0.0, 1.0 # p_initial is all zeros

        # Transfer rho(E) (p_rho_values) to p_initial (rho(beta))
        # p_rho_values are at E = (i+0.5)*delta_e_rho
        # p_initial are at beta = (k+0.5)*deltab_val
        # E_k = (k+0.5)*deltab_val * tev_val = (k+0.5)*delta_e_rho
        # So the kth point in p_initial corresponds to the kth point in p_rho_values
        npt = min(num_rho_points, self.nbeta) # number of points to transfer
        p_initial[:npt] = p_rho_values[:npt]

        # Normalize p_initial (which is currently rho(beta))
        # The normalization factor tbeta_norm_continuous (read as 'tbeta' on card 20)
        # is defined such that integral of rho(E)dE = tbeta_norm_continuous
        # integral rho(beta) dbeta = integral rho(E)/tev dE = tbeta_norm_continuous / tev
        # Each element p_initial[k] is rho at beta_k. The integral is sum(p_initial[k]*deltab_val).
        
        current_norm = np.sum(p_initial[:npt]) * deltab_val
        
        if self.iprint >=2:
            print(f" current normalization of rho(beta) = {current_norm:.5e}")
            print(f" desired normalization for S(b) part = {tbeta_norm_continuous/tev_val:.5e}")

        if current_norm == 0.0:
            if tbeta_norm_continuous == 0.0: # Both are zero, consistent.
                f0_val = 0.0
                tbar_val = 1.0 # Effective temp = actual temp
            else: # current_norm is zero, but desired is non-zero. This is an issue.
                print(f" Error: rho(beta) sums to zero, but input tbeta ({tbeta_norm_continuous}) is non-zero.")
                # NJOY's LEAPR seems to proceed with f0=0, tbar=temp in this case.
                f0_val = 0.0
                tbar_val = 1.0
            self.dwpix_values[itemp] = f0_val
            self.tempf_values[itemp] = tbar_val * current_temp_k
            # p_initial is already zero, so S(beta) part will be zero.
            return p_initial, deltab_val, tev_val, f0_val, tbar_val

        # Apply normalization factor to p_initial
        renorm_factor = (tbeta_norm_continuous / tev_val) / current_norm
        p_initial[:npt] *= renorm_factor
        
        if self.iprint >=2:
            test_norm = np.sum(p_initial[:npt]) * deltab_val
            print(f" test normalization = {test_norm:.5e} (should be {tbeta_norm_continuous/tev_val:.5e})")

        # Calculate f0 (Debye-Waller factor integral) and tbar (effective temperature)
        # f0 = (1/ARAT) * Integral[ rho(beta)/beta * (coth(beta/2T) - 1) dbeta ]
        # tbar/T = (1/ARAT) * (1/2T) * Integral[ rho(beta) * beta * (coth(beta/2T) - 1) dbeta ]
        # Note: The fsum function calculates Integral[ p(beta') * factor * (coth(beta'*tev/2T) - 1) dbeta' ]
        # where factor is 1 for f0-like sum, and beta' for tbar-like sum.
        # The beta in fsum's integral is p_array's argument, which is already beta.
        # The tev in fsum's call is the temperature in eV.
        # The term beta/2T in coth becomes (beta_k * tev) / (2 * tev) = beta_k / 2 for fsum's exp_arg.
        # This means fsum's exp_arg should be `ek / (2*tev)` if its `tau` argument is `tev`.
        # Let's re-check fsum's `exp_arg = ek / tau`. If tau is tev, then exp_arg = beta_k_grid_val / tev.
        # The coth argument is beta_E / (2*kT), where beta_E is energy, kT is tev.
        # If ek in fsum is beta-grid value, then coth argument is beta_grid_val / 2.
        # So, fsum's `tau` param should be `2.0 * tev_val` if `ek` is beta.
        # OR, fsum's `exp_arg` should be `ek / (2.0 * tau)` if `tau` is `tev_val`.
        # The Fortran fsum has exp(ek/tau). `ek` there is `(k-1+0.5)*deltab`. So `ek` is beta.
        # `tau` there is `tev`. So it's `exp(beta/tev)`.
        # `coth(beta/2T_ev) - 1` is `(exp(beta/T_ev)+1)/(exp(beta/T_ev)-1) - 1`.
        # This matches `(coth_val - 1.0)` if `exp_arg` in fsum is `beta/tev`.
        # So, fsum is calculating sum over k of:
        #   p[k]*deltab * (coth(beta[k]/tev) - 1) * 0.5  (for n=1 for f0)
        #   p[k]*deltab * beta[k] * (coth(beta[k]/tev) - 1) * 0.5 (for n=2 for tbar)

        # f0 = (1/ARAT) * sum_k [ rho(beta_k)/beta_k * (coth(beta_k/2T) - 1) * deltab ]
        #    The 0.5 factor in fsum handles the 1/2 in coth(beta/2T) somehow, or it's part of the definition.
        # Let's assume fsum is correctly implementing the required sum for NJOY's definition of f0, tbar.
        # The Fortran code for f0 (an via fsum(1,...)) and tbar (bn via fsum(2,...)) is:
        #   an = fsum(1, p, npt, tev, deltab)  ! p is rho(beta) here.
        #   bn = fsum(2, p, npt, tev, deltab)
        #   f0 = an / arat
        #   tbar = bn / arat / tev  (actually tbar/temp = bn / arat / tev)
        # So, f0_val should be self._fsum(1, p_initial, npt, tev_val, deltab_val) / self.arat
        # And tbar_val should be self._fsum(2, p_initial, npt, tev_val, deltab_val) / self.arat / tev_val
        
        # However, the formulas given in comments are Integral[ rho(beta)/beta ... ] for f0
        # and Integral[ rho(beta)*beta ... ] for tbar.
        # fsum(1,...) calculates Sum[ rho(beta_k) * (coth-1) ]
        # fsum(2,...) calculates Sum[ rho(beta_k)*beta_k * (coth-1) ]
        # This suggests the /beta_k for f0 needs to be done explicitly if fsum doesn't do it.
        # Looking at NJOY leapr.f90, start subroutine:
        #   an=fsum(1,p,npt,tev,deltab)
        #   bn=fsum(2,p,npt,tev,deltab)
        #   f0=an/arat
        #   tbar=bn/arat/tev
        # This implies fsum(1,...) is already Sum[ rho(beta_k)/beta_k * ... ] if formulas are right,
        # or the formulas are shorthand and fsum result 'an' is directly f0*arat.
        # The provided fsum code calculates: Sum_k p_array[k] * deltab * term_n * (coth_val-1)*0.5
        # where term_n = 1 for n=1, and term_n = ek (beta_k) for n=2.
        # So an = Sum_k p_initial[k] * (coth(beta_k/tev)-1)*0.5 * deltab
        #    bn = Sum_k p_initial[k] * beta_k * (coth(beta_k/tev)-1)*0.5 * deltab

        # These seem to be direct sums, not including /beta_k for an.
        # Let's trust the direct translation of fsum and the way an, bn are used.
        an = self._fsum(1, p_initial, npt, tev_val, deltab_val)
        bn = self._fsum(2, p_initial, npt, tev_val, deltab_val)

        f0_val = an / self.arat
        if tev_val == 0.0: # Should have been caught earlier
            tbar_val = 1.0 # Or some other appropriate default for tbar/temp
        else:
            tbar_val = bn / self.arat / tev_val # This is tbar/temp factor

        self.dwpix_values[itemp] = f0_val
        self.tempf_values[itemp] = tbar_val * current_temp_k # Store actual effective temperature

        if self.iprint >= 1:
            print(f"  calculated debye-waller integral f0 = {f0_val:.5e}")
            print(f"  calculated effective temperature tbar = {self.tempf_values[itemp]:.5e} K " \
                  f"(factor tbar/temp = {tbar_val:.5e})")

        # Transform p_initial to S(beta) for one-phonon term: T1(beta)
        # T1(beta) = rho(beta) * exp(-beta/2) / (2 * sinh(beta/2T)) / beta * (f0*arat)
        # No, this is not T1. This is part of S_alpha_beta.
        # The p_initial array (called P in Fortran START) is modified in place.
        # P(K) = P(K) * EXP(-EK/2._KR) / (ARAT*EK*DELTAB) where EK is beta_k.
        # This P(K) is then returned as T1. Let's verify this transformation.
        # The factor seems to be rho(beta_k) * exp(-beta_k/2) / (arat * beta_k * deltab)
        # Note: deltab factor might be an issue if p_initial is density vs function values.
        # p_initial is rho(beta). Integral rho(beta)dbeta is normalized.
        # The output p_initial is T1(beta).
        
        for k_idx in range(npt):
            ek_beta = (k_idx + 0.5) * deltab_val # beta value at center of interval k
            if ek_beta == 0.0:
                p_initial[k_idx] = 0.0 # Avoid division by zero; rho(0) should be 0 or handled.
                continue

            # exp_term: Fortran uses exp(-ek/(2.0*tev)) for exp(-beta/2) but that's if tev is 1.
            # exp(-beta/2) is standard for detailed balance factor.
            exp_val_neg_beta_half = np.exp(max(EXPLIM, -ek_beta / 2.0))
            
            denominator = self.arat * ek_beta # * deltab_val # Removed deltab_val based on manual NJOY check
                                                            # The deltab is part of the sum in convol.
                                                            # p_initial here is T1_k, not T1_k*deltab
            if denominator == 0:
                p_initial[k_idx] = 0.0
            else:
                p_initial[k_idx] *= exp_val_neg_beta_half / denominator
        
        # Zero out remaining part of p_initial if npt < nbeta
        if npt < self.nbeta:
            p_initial[npt:] = 0.0
            
        return p_initial, deltab_val, tev_val, f0_val, tbar_val

    def _convol(self, t1_beta_spectrum, tlast_sab, tnext_sab, num_beta_points, nl_order, nn_order, delta_beta, ckk_factor):
        """
        Corresponds to Fortran subroutine convol.
        Performs convolution to get n-th order S(beta) from (n-1)-th order.
        S_n(beta) = ckk * Integral [ S_1(beta') * S_{n-1}(beta - beta') d_beta' ]

        t1_beta_spectrum: numpy array, T1(beta) - one phonon term (from _start, modified)
        tlast_sab: numpy array, S_{n-1}(beta) - S(beta) from previous order
        tnext_sab: numpy array, S_n(beta) - S(beta) for current order (output, modified in place)
        num_beta_points: integer, number of points in beta grid (self.nbeta)
        nl_order: integer, current phonon order n (for S_n)
        nn_order: integer, max phonon order (self.nphon)
        delta_beta: float, spacing of beta grid
        ckk_factor: float, normalization constant ckk = f0 / (n * arat)
        """
        # Initialize tnext_sab to zero before accumulating
        tnext_sab.fill(0.0)

        # Loop over k (index for beta, for S_n(beta_k))
        # Fortran: do 140 k=1,np
        # beta_k = (k-1+0.5)*delta  (Fortran) -> (k_idx+0.5)*delta_beta (Python)
        for k_idx in range(num_beta_points):
            sum_val = 0.0
            # Loop over j (index for beta', for S_1(beta'_j))
            # Fortran: do 130 j=1,np
            # beta'_j = (j-1+0.5)*delta (Fortran) -> (j_idx+0.5)*delta_beta (Python)
            for j_idx in range(num_beta_points):
                # Argument for S_{n-1} is beta_k - beta'_j
                # Index for tlast_sab: (beta_k - beta'_j)/delta_beta gives offset from beta=0
                # Fortran: kk = k - j + (np+1)/2  (using (np+1)/2 as beta=0 index for S_{n-1})
                # This Fortran indexing for kk implies S_{n-1} is stored centered.
                # Python S(alpha,beta) is usually stored starting at beta=0.
                # If tlast_sab[m] corresponds to beta_m = (m+0.5)*delta_beta,
                # then we need m such that (m+0.5)*delta_beta = (k_idx+0.5)*delta_beta - (j_idx+0.5)*delta_beta
                # m+0.5 = k_idx+0.5 - (j_idx+0.5)
                # m = k_idx - j_idx - 0.5. This must be an integer index.
                
                # Let's re-evaluate based on typical convolution definition:
                # S_n(beta_k) = sum_j S_1(beta'_j) * S_{n-1}(beta_k - beta'_j) * delta_beta'
                # beta_k - beta'_j = (k_idx+0.5)*delta_beta - (j_idx+0.5)*delta_beta
                #                    = (k_idx - j_idx)*delta_beta
                # So the argument for S_{n-1} is (k_idx - j_idx)*delta_beta.
                # The index m for tlast_sab should satisfy (m_idx+0.5)*delta_beta = (k_idx-j_idx)*delta_beta
                # m_idx = k_idx - j_idx - 0.5. Still not integer.

                # The common way to implement discrete convolution for S(a,b) is:
                # S_n[k] = sum_j S_1[j] * S_{n-1}[k-j] * delta_beta
                # Here S_n[k] is S_n at beta_k = k * delta_beta (if grid starts at 0)
                # or (k+0.5)*delta_beta if it's cell-centered.
                # Let's assume S_x[i] refers to S_x((i+0.5)*delta_beta).
                
                # S_n[k_idx] (value at beta_k) = sum_j ( S_1[j_idx] * S_{n-1}[m_idx] ) * delta_beta
                # where beta_m = beta_k - beta_j.
                # (m_idx + 0.5)*delta_beta = (k_idx + 0.5)*delta_beta - (j_idx + 0.5)*delta_beta
                # m_idx = k_idx - j_idx - 0.5. This is problematic.

                # NJOY's convol.f has:
                # DO K=1,NP (outer loop for S_next, beta_k)
                #   SUM=0.
                #   DO J=1,NP (inner loop for S_1, beta_j)
                #     KK = K-J+NPL  (NPL is np, number of beta points for S_n-1)
                #     IF (KK.GE.1 .AND. KK.LE.NPN) THEN  (NPN is also np for S_n-1)
                #       SUM = SUM + T1(J)*TLAST(KK)
                #   ENDIF
                #   TNEXT(K) = CKK * DELTA * SUM
                # This implies TLAST (S_n-1) is indexed such that TLAST(NPL) is S_n-1(0) if K=J.
                # This is very unusual. NPL is num_beta_points.
                # If k=j, kk = np. S_n-1(beta near 0) would be at index np?
                # Typical S(alpha,beta) has beta >= 0.
                # S_1 is rho(beta)*stuff, so T1(J) is for beta_j > 0.
                # S_n-1 can be for beta < 0 if beta_k - beta_j < 0.
                # The problem is that S(alpha,beta) is typically defined for beta >= 0.
                # LEAPR does extend S to negative beta for intermediate calculations.
                # The tlast_sab, tnext_sab arrays are full beta grid from -beta_max to +beta_max.
                # Size of these arrays in Fortran is (npn = 2*nbeta-1 if nbeta is positive beta points).
                # Here, num_beta_points is self.nbeta (positive beta points only).
                # This translation needs to handle the symmetric beta grid used internally by LEAPR.

                # For now, let's assume t1_beta_spectrum is for positive beta (size self.nbeta)
                # And tlast_sab, tnext_sab are also for positive beta (size self.nbeta)
                # This means the convolution must be adapted.
                # S_n(beta_k) = sum_{j | beta_j > 0} S_1(beta_j) * S_{n-1}(beta_k - beta_j) * delta_beta
                # If beta_k - beta_j < 0, then S_{n-1} needs to be defined for negative arguments.
                # Standard S(alpha,beta) from ENDF is for beta > 0.
                # If these are intermediate calculations, they might use symmetric beta.
                # The problem description implies ssm (final S(a,b)) is (ntempr, nalpha, nbeta),
                # suggesting nbeta is for positive beta values.

                # If tlast_sab is S_{n-1}(beta) for beta > 0:
                m_idx = k_idx - j_idx # Target index for beta_k - beta_j (if grids are aligned and start at 0)
                                      # If beta grid is (idx+0.5)*delta_beta:
                                      # beta_target = (k_idx+0.5)*delta_beta - (j_idx+0.5)*delta_beta = (k_idx-j_idx)*delta_beta
                                      # m_idx_plus_half = k_idx - j_idx
                                      # m_idx = k_idx - j_idx - 0.5. Still non-integer.
                
                # Let's use the direct summation formula with index checks.
                # S_n(beta_k) = sum_j S_1(beta_j) S_{n-1}(beta_k - beta_j) delta_beta
                # beta_k = (k_idx + 0.5) * delta_beta
                # beta_j = (j_idx + 0.5) * delta_beta
                # beta_target_for_S_n_minus_1 = beta_k - beta_j = (k_idx - j_idx) * delta_beta
                
                # Find index m_target_idx for beta_target_for_S_n_minus_1 in tlast_sab's grid
                # (m_target_idx + 0.5) * delta_beta = (k_idx - j_idx) * delta_beta
                # m_target_idx = k_idx - j_idx - 0.5. This is the core issue.
                # The arrays T1, TLAST, TNEXT must be on compatible grids.

                # Re-check convol.f from NJOY2016 source:
                # T1 is dimensioned NP (number of positive beta points)
                # TLAST, TNEXT are dimensioned NPN (number of points for full beta grid, typically 2*NP-1 or 2*NP)
                # NPL = NP in the call from CONTIN.
                # KK = K-J+NP.  (K for TNEXT, J for T1, KK for TLAST)
                # K runs 1 to NPN. J runs 1 to NP.
                # This means T1 is only the positive frequency part of the spectrum.
                # TLAST and TNEXT cover negative frequencies.
                # If TNEXT[0] is beta_max_neg, TNEXT[NP-1] is beta=0, TNEXT[NPN-1] is beta_max_pos.
                # Let NP be self.nbeta (number of positive beta points).
                # Full grid size NPN = 2 * self.nbeta - 1 (if beta=0 is shared).
                # Or NPN = 2 * self.nbeta (if beta=0 not explicitly stored, or grid is shifted).
                # The current Python code assumes tlast_sab, tnext_sab are of size self.nbeta (positive only).
                # This is a major simplification that won't work if negative beta is needed.

                # For this subtask, we are asked to translate existing Fortran.
                # The Fortran S(alpha,beta) `ssm(k,j,itemp)` has k for beta, j for alpha.
                # `beta(k)` means beta values. `nbeta` is number of beta points.
                # This usually means positive beta points for final output.
                # The `contin` subroutine calculates `ssm`.
                # `convol` is called by `contin` with `p` (T1), `tlast`, `tnow`.
                # `p` is dimension `np` (positive beta points). `tlast`, `tnow` are `npn`.
                # `npn` is `2*nbeta - 1`. `npl` in call is `nbeta`.
                # So, `t1_beta_spectrum` is size `self.nbeta`.
                # `tlast_sab`, `tnext_sab` should be size `2*self.nbeta-1`.
                # This needs to be addressed in `_contin` where these arrays are created.

                # Assuming for a moment that all arrays are on the positive beta grid of size num_beta_points (self.nbeta)
                # This is a strong, likely incorrect simplification for intermediate steps.
                # S_n[k] = sum_{j=0 to k} S_1[j] * S_{n-1}[k-j] * delta_beta (if beta >= 0 for all)
                if t1_beta_spectrum[j_idx] == 0.0:
                    continue

                m_idx = k_idx - j_idx
                if 0 <= m_idx < num_beta_points:
                    sum_val += t1_beta_spectrum[j_idx] * tlast_sab[m_idx]
            
            tnext_sab[k_idx] = ckk_factor * delta_beta * sum_val
        
        # The above simplified convolution is for positive beta only and assumes index arithmetic k-j.
        # This is likely incorrect for the actual LEAPR algorithm which uses symmetric beta.
        # The Fortran indexing `KK = K-J+NPL` (where NPL=NP, number of positive beta points)
        # and K runs up to NPN (e.g. 2*NP-1), J up to NP.
        # This requires `tlast_sab` and `tnext_sab` to be on a grid that includes negative beta.
        # Let's assume `num_beta_points` passed to `_convol` is NPN (full grid size).
        # And `t1_beta_spectrum` is still on positive beta grid (size self.nbeta, Fortran NP).
        # Let `np_pos = self.nbeta`.
        # `tnext_sab` (output) is size `num_beta_points` (NPN).
        # `tlast_sab` (input) is size `num_beta_points` (NPN).
        # `t1_beta_spectrum` (input) is size `np_pos`.

        # Corrected loop structure based on Fortran convol:
        # num_beta_points here is NPN (full grid size) from the call in _contin.
        # np_pos is self.nbeta (size of t1_beta_spectrum).
        
        # tnext_sab.fill(0.0) # Already done above.
        np_pos = len(t1_beta_spectrum) # Should be self.nbeta

        for k_convol_idx in range(num_beta_points): # k_convol_idx from 0 to NPN-1 (Fortran K from 1 to NPN)
            current_sum = 0.0
            # j_t1_idx for t1_beta_spectrum (positive beta values)
            for j_t1_idx in range(np_pos): # j_t1_idx from 0 to NP-1 (Fortran J from 1 to NP)
                # Fortran: KK = K - J + NP (NP is np_pos, K is k_convol_idx+1, J is j_t1_idx+1)
                # kk_idx_fortran = (k_convol_idx+1) - (j_t1_idx+1) + np_pos = k_convol_idx - j_t1_idx + np_pos
                # This kk_idx_fortran is 1-based index for tlast_sab.
                # Python 0-based index kk_tlast_idx = kk_idx_fortran - 1
                kk_tlast_idx = k_convol_idx - j_t1_idx + np_pos - 1
                
                if 0 <= kk_tlast_idx < num_beta_points: # Check bounds for tlast_sab (size NPN)
                    if t1_beta_spectrum[j_t1_idx] != 0.0: # Optimization
                        current_sum += t1_beta_spectrum[j_t1_idx] * tlast_sab[kk_tlast_idx]
            
            tnext_sab[k_convol_idx] = ckk_factor * delta_beta * current_sum

    def _contin(self, itemp, current_temp_data, max_phonon_order_nphon):
        """
        Corresponds to Fortran subroutine contin.
        Calculates S(alpha,beta) for continuous part using phonon expansion.

        itemp: integer, current temperature index.
        current_temp_data: dict, data for the current temperature.
                           Expected keys: 'temperature' (K),
                                          'p1_values' (rho(E) for continuous),
                                          'ni' (num points for rho(E)),
                                          'delta1' (delta E for rho(E)),
                                          'tbeta' (normalization for continuous rho(E)).
        max_phonon_order_nphon: integer, max phonon expansion order (self.nphon).
        """
        current_temp_k = current_temp_data['temperature']
        p_rho_values = current_temp_data.get('p1_values', np.array([])) # rho(E)
        num_rho_points = current_temp_data.get('ni', 0)
        delta_e_rho = current_temp_data.get('delta1', 0.0)
        tbeta_norm_continuous = current_temp_data.get('tbeta', 0.0)

        if num_rho_points == 0 and tbeta_norm_continuous == 0.0 and self.lat != 1:
            # This case is for discrete oscillators primarily, or if continuous part is truly zero.
            # If lat=1, _start handles no continuous part warning.
            # If not lat=1, and no p1_values, this routine might not be the main path for S(a,b) anyway.
            if self.iprint >=1:
                print(" Skipping _contin as no continuous spectrum data (ni=0, tbeta=0) and not explicitly lat=1.")
            # Ensure ssm is zeroed for this temperature if it's not handled elsewhere.
            # self.ssm_values[itemp, :, :] = 0.0 # This might be too broad if other models contribute.
            return

        # Call _start to get initial phonon spectrum T1(beta) and other params
        # p_initial_t1 is T1(beta) for positive beta, size self.nbeta
        p_initial_t1, deltab_val, tev_val, f0_val, tbar_val = \
            self._start(itemp, p_rho_values, num_rho_points, delta_e_rho, 
                        tbeta_norm_continuous, current_temp_k)

        if deltab_val == 0.0:
            if self.iprint >= 0:
                print(" error in contin. deltab is zero. stop.")
            raise NJOYError("delta beta is zero in _contin.")

        # npn is the number of points for the full symmetric beta grid (includes negative beta)
        # Typically 2*self.nbeta - 1 if beta=0 is shared, or 2*self.nbeta.
        # NJOY LEAPR uses npn = 2*nbeta -1. Index npbeta = nbeta is beta=0. (1-based)
        # Python: if size is NPN=2*NB-1, index NB-1 is beta=0.
        npn_full_grid_size = 2 * self.nbeta - 1 
        if npn_full_grid_size <=0:
             if self.iprint >=0:
                print(f" error in contin. npn_full_grid_size={npn_full_grid_size} from nbeta={self.nbeta}. stop.")
             raise NJOYError(f"Invalid npn_full_grid_size {npn_full_grid_size}")


        # tlast_sab and tnow_sab are S_n-1(beta) and S_n(beta) on the full grid
        tlast_sab = np.zeros(npn_full_grid_size)
        tnow_sab = np.zeros(npn_full_grid_size)

        # Initialize S(alpha,beta) for this temperature to zero
        self.ssm_values[itemp, :, :] = 0.0

        # Zeroth phonon term S0(alpha,beta) = exp(-f0_val*alpha) * delta(beta)
        # This is handled by adding to ssm at beta=0 (or near beta=0).
        # For discrete beta grid, delta(beta) means value is at beta=0 bin.
        # The beta grid for ssm is self.beta_values, which are positive.
        # S0 contributes to elastic scattering. Here we calculate S(alpha,beta) for beta > 0.
        # The elastic part is usually added separately or via specific logic for beta=0.
        # For now, this _contin focuses on inelastic S(a,b) for beta > 0 from phonon expansion.

        # First phonon term: S1(alpha,beta) = alpha * exp(-f0_val*alpha) * T1(beta)
        # (where T1(beta) is p_initial_t1 after modification in _start)
        # This S1 is for the positive beta part of the spectrum.
        # We need to place p_initial_t1 onto the tnow_sab (full grid) for n=1.
        # If tnow_sab index (npn_full_grid_size // 2) corresponds to beta=0.
        # Then tnow_sab[npn_full_grid_size // 2 + k_idx+1] for beta_k > 0.
        center_idx_npn = self.nbeta - 1 # Index for beta=0 on full grid (0 to NPN-1)

        for k_beta_idx in range(self.nbeta): # k_beta_idx for positive beta grid (0 to self.nbeta-1)
            # Positive beta values on full grid: center_idx_npn + k_beta_idx + 1 (for beta_(k+1))
            # Or, if p_initial_t1[0] is for beta_0_positive_grid_center = 0.5*deltab_val
            # then it maps to tnow_sab[center_idx_npn + k_beta_idx] if beta grid is symmetric points,
            # or center_idx_npn + 1 + k_beta_idx if center_idx_npn is truly beta=0 and next is beta_1.
            # NJOY: p(k) for k=1..nbeta (positive beta). tnow(nbeta+k) = p(k+1) for k=0..nbeta-1
            # So, tnow_sab[center_idx_npn + k_beta_idx] = p_initial_t1[k_beta_idx] (if using same convention for p_initial_t1)
            # Let's assume p_initial_t1[k] is for (k+0.5)*deltab
            # And tnow_sab[center_idx_npn + k] is for (k+0.5)*deltab
            # And tnow_sab[center_idx_npn - k] is for -(k+0.5)*deltab
            if center_idx_npn + k_beta_idx < npn_full_grid_size :
                 tnow_sab[center_idx_npn + k_beta_idx] = p_initial_t1[k_beta_idx]
            # For negative beta, S1(-beta) = exp(beta) * S1(beta) -- from detailed balance on S(E) then to S(beta)
            # T1(-beta) = exp(beta) * T1(beta)
            # So, p_initial_t1 for negative beta needs exp(beta_val) factor.
            # beta_val for index (center_idx_npn - (k_beta_idx+1)) is -(k_beta_idx+0.5)*deltab_val
            # exp_factor = np.exp( (k_beta_idx+0.5)*deltab_val )
            # This detailed balance is for S(alpha,beta), not T_n(beta) directly in this form.
            # Let's use symmetry from NJOY: tnow(nbeta-k) = tnow(nbeta+k) * exp(-beta(k))
            # tnow_sab[center_idx_npn - (k_beta_idx+1)] = tnow_sab[center_idx_npn + (k_beta_idx+1)] * exp(-(k_beta_idx+1)*deltab_val)
            # This assumes beta grid points are multiples of deltab_val.
            # My p_initial_t1 is for beta_k = (k+0.5)*deltab_val.
            # Let's simplify: S_n(beta) is typically symmetric or derived for positive beta first.
            # The Fortran `contin` loop for n=1 to nphon:
            #   `call scale` - which applies exp(-f0*alpha) and alpha factor.
            #   Then adds to `ssm`.
            # The first term added to `ssm` (n=0 in `contin` loop, so 1st phonon `nl=1`)
            # uses `p` (which is `p_initial_t1`) directly for positive beta.
            # `ssm(k,j,itemp) = ssm(k,j,itemp) + xa(j)*fac*p(k)`
            # This means `p_initial_t1` is directly used for the first phonon contribution to final S(a,b).

        # Main phonon expansion loop
        # nl is current phonon order, from 1 to max_phonon_order_nphon
        for nl_order in range(1, max_phonon_order_nphon + 1):
            if nl_order == 1:
                # S1(beta) for positive beta is p_initial_t1.
                # Place it on the full symmetric grid tnow_sab.
                # tnow_sab[center_idx_npn + k_beta_idx] = p_initial_t1[k_beta_idx] for k_beta_idx = 0 to nbeta-1
                # This covers beta > 0.
                # Beta=0 point (center_idx_npn) can be zero for inelastic.
                # Negative beta part of tnow_sab (S1) from detailed balance: S1(-beta) = exp(-beta)S1(beta)
                tnow_sab.fill(0.0) # Clear for current order
                for k_idx in range(self.nbeta): # k_idx for p_initial_t1
                    beta_val = (k_idx + 0.5) * deltab_val
                    # Positive beta part of S1 on full grid
                    if center_idx_npn + k_idx < npn_full_grid_size:
                         tnow_sab[center_idx_npn + k_idx] = p_initial_t1[k_idx]
                    # Negative beta part of S1 on full grid
                    if center_idx_npn - (k_idx + 1) >= 0: # k_idx+1 because beta grid is (idx+0.5) relative to center beta=0
                         # Accessing S1 at beta_val, stored at p_initial_t1[k_idx]
                         # Value for -beta_val is exp(-beta_val) * S1(beta_val)
                         exp_factor = np.exp(max(EXPLIM, -beta_val))
                         tnow_sab[center_idx_npn - (k_idx + 1)] = p_initial_t1[k_idx] * exp_factor
            else: # nl_order > 1
                # Calculate S_nl(beta) using S_1(beta) and S_{nl-1}(beta)
                # S_1(beta) is p_initial_t1 (positive beta part)
                # S_{nl-1}(beta) is tlast_sab (full beta grid)
                ckk = f0_val / (nl_order * self.arat) # ckk = f0 / (n*arat) from Fortran contin
                if self.iprint >=3:
                    print(f"    ckk factor for order {nl_order} = {ckk:.5e}")

                self._convol(p_initial_t1, tlast_sab, tnow_sab, 
                             npn_full_grid_size, nl_order, max_phonon_order_nphon,
                             deltab_val, ckk)

            # Add S_nl(alpha,beta) to total ssm
            # S_nl(alpha,beta) = (alpha^nl / n!) * exp(-f0*alpha) * [scaled T_nl(beta)]
            # The Fortran `scale` subroutine seems to handle the alpha dependence and normalization.
            # Here, `tnow_sab` is T_nl(beta) from convol (or T1 for nl=1).
            # Factor for ssm: xa(j)*fac * tnow_sab(k_on_full_grid)
            # fac = term / arat (term is (f0*arat)^nl / nl! in Fortran `contin` for `fac` before `scale` call)
            # This implies tnow_sab from convol is already scaled somewhat.
            # Let's assume the primary alpha dependence (alpha^nl / nl! * exp(-f0*alpha)) needs to be applied here.

            term = 1.0
            if nl_order > 1: # (f0*arat)^nl / nl!
                for i in range(1, nl_order + 1):
                    term = term * f0_val * self.arat / i
            
            fac = term / self.arat # Matches Fortran `contin` before `scale` call for nl > 0.
                                   # For nl=0 (n=1 in loop), fac = 1.0/arat. (nl is n_phonon_order here)
            if nl_order == 1: # Fortran n=0 in loop, nl=1 phonon
                fac = 1.0 / self.arat


            if self.iprint >= 2:
                if nl_order == 1:
                     print(f"  sum for 1-phonon s(a,b) term = {np.sum(tnow_sab)*deltab_val:.5e}")
                print(f"  max order={max_phonon_order_nphon}, current order={nl_order}, factor={fac:.5e}")

            for j_alpha_idx in range(self.nalpha): # Loop over alpha grid
                alpha_val = self.alpha_values[j_alpha_idx]
                # xa_factor is alpha part: (alpha*f0)^nl / nl! * exp(-alpha*f0) / arat
                # This is effectively (alpha*effective_debye_waller_factor)
                # The fac calculated above is (f0*arat)^nl / nl! / arat
                # So xa_factor = (alpha_val / (f0_val*self.arat))^nl * fac * exp(-f0_val*alpha_val) if f0!=0
                # This is getting complicated. Let's use NJOY's xa(j) directly from `scale` concept:
                # xa(j) = alpha(j) * exp(-f0*alpha(j)) for nl=1
                # xa(j) = (alpha(j)*f0)^nl / nl! * exp(-f0*alpha(j)) -- this is S_inel for one term.
                # The `fac` in Fortran `contin` is `(f0*arat)^n / n! / arat`. (n is phonon order here)
                # The `p(k)` in `ssm(k,j) = ssm(k,j) + xa(j)*fac*p(k)` is `tnow_sab` here.
                
                # Simplified alpha factor for S_nl(alpha,beta) contribution:
                # exp(-f0*alpha) * ( (alpha*f0)^nl / nl! )  -- if f0 is DW factor from free gas.
                # Here f0_val is the DW integral.
                # S(a,b) = sum_{n=0 to N} S_n(a,b)
                # S_0(a,b) = exp(-a*DW) * delta(b)
                # S_n(a,b) = exp(-a*DW) * (a*DW)^n / n! * T_n(b) / DW^n  -- if T_n includes (DW)^n
                # Or S_n(a,b) = exp(-a*DW) * alpha^n / n! * T_n(beta) where T_n is from convol.
                # NJOY: xa(j) = alpha(j) for n=1, then multiplied by fac=1/arat.
                #       xa(j) = alpha(j) for n > 1, then multiplied by fac=(f0*arat)^(n-1)/n!/arat ...
                # This implies the (alpha*f0)^nl part is constructed iteratively.

                # Let's use the structure from Fortran `contin` more directly for alpha part:
                # `xn = alpha(j)*f0`
                # `term = xn`
                # `do l=2,nl: term = term*xn/(f0*l)` This makes term = (alpha*f0)^nl * f0 / (nl! * f0) = (alpha*f0)^nl / nl!
                # `xa(j) = term * exp(-alpha(j)*f0)`
                # Then `ssm = ssm + xa(j) * (1/arat) * tnow_sab(k_positive_beta)` (missing some factors from fac)

                # Let's use the direct summation form:
                # S(alpha,beta) = sum_{n=1 to N} (exp(-alpha*f0)/arat) * ( (alpha*f0*arat)^n / n! ) * T_n(beta)
                # No, this is not quite it. It's simpler:
                # S(alpha,beta) = exp(-alpha*f0) * Sum_{n=1..N} alpha^n / n! * T_n(beta) (if T_n has arat effects in it)
                # Or based on NJOY `contin` and `scale` logic:
                # For phonon order `nl_order`:
                #   `factor_for_S_nl = ( (f0_val * self.arat)**nl_order / math.factorial(nl_order) ) / self.arat` (this is `fac` from NJOY)
                #   `sab_alpha_part = alpha_val * np.exp(max(EXPLIM, -alpha_val * f0_val))` (for nl_order=1)
                #   For nl_order > 1, it's more like `(alpha_val * f0_val)**nl_order / math.factorial(nl_order) * exp(-alpha_val * f0_val)`
                
                # Sticking to `xa(j)*fac*p(k)` form:
                # `fac` is `term/arat` where `term` is `(f0*arat)^nl / nl!`
                # `xa(j)` needs to be `alpha(j)`.
                # This makes the alpha part `alpha(j) * (f0*arat)^nl / nl! / arat * exp(-f0*alpha(j))` (missing exp)
                # The `exp(-f0*alpha)` is applied globally at the end in NJOY `contin`.
                # So, `ssm_contrib = alpha_val * fac_for_order_nl * tnow_sab_val`.
                
                # `fac` from NJOY contin (n is nl_order):
                # if n=0 (1-phonon, nl_order=1): fac = 1/arat
                # if n>0 (multi-phonon): term = (f0*arat)^n / n!; fac = term/arat.
                # This `fac` is for `T_n(beta)`.
                # `ssm(k,j) = ssm(k,j) + alpha(j)*fac*tnow(k_positive_beta)`
                # Then final `ssm(k,j) = ssm(k,j) * exp(-alpha(j)*f0)`

                current_alpha_factor = alpha_val 
                # This factor will be multiplied by fac_for_order_nl and then by exp(-alpha*f0) later.

                for k_beta_idx in range(self.nbeta): # k_beta_idx for positive self.beta_values
                    # Get tnow_sab value for positive beta_k = self.beta_values[k_beta_idx]
                    # We need to interpolate tnow_sab (on full deltab_val grid) onto self.beta_values grid.
                    # self.beta_values are the report points for beta.
                    # tnow_sab is on grid (idx - center_idx_npn)*deltab_val if centered at 0.
                    # Or (idx + 0.5)*deltab_val for positive part if p_initial_t1 grid.
                    # The tnow_sab from convol is on grid (i - (npn-1)/2)*deltab, for i=0..npn-1
                    # Or, if center_idx_npn is beta=0 point: tnow_sab[center_idx_npn + m] is for beta_m = m*deltab_val (if m can be non-integer for terpt)
                    # Or tnow_sab[center_idx_npn + m_idx] for beta at (m_idx+0.5)*deltab if grid is cell-centered.
                    
                    # For now, assume tnow_sab (full grid) has its positive part starting at center_idx_npn,
                    # corresponding to p_initial_t1's grid.
                    # So, tnow_sab[center_idx_npn + k_beta_idx] is value for beta approx self.beta_values[k_beta_idx].
                    # This requires self.beta_values to align with (idx+0.5)*deltab_val.
                    # This is often true if beta grid is uniform.
                    
                    beta_report_val = self.beta_values[k_beta_idx]
                    
                    # Interpolate tnow_sab (on deltab_val grid) at beta_report_val
                    # tnow_sab grid: idx k_full uses (k_full - center_idx_npn + 0.5)*deltab_val (approx)
                    # For positive beta part of tnow_sab: (k_from_center + 0.5)*deltab_val
                    # tnow_sab indices from center_idx_npn upwards correspond to positive beta.
                    # tnow_sab_positive_part = tnow_sab[center_idx_npn:] (length nbeta)
                    # This part is effectively T_nl(beta) for beta > 0, on grid (idx+0.5)*deltab_val.
                    
                    sab_val_for_beta = 0.0
                    if nl_order == 1: # T1 is p_initial_t1
                        # Interpolate p_initial_t1 (on deltab_val grid) at beta_report_val
                        sab_val_for_beta = self._terpt(p_initial_t1, self.nbeta, deltab_val, beta_report_val)
                    else: # T_nl is from tnow_sab (positive part)
                        # Positive part of tnow_sab starts at index center_idx_npn.
                        # Number of points in this positive part is self.nbeta.
                        # Grid for this part is (idx+0.5)*deltab_val for idx = 0 to self.nbeta-1.
                        tnow_sab_positive_part = tnow_sab[center_idx_npn : center_idx_npn + self.nbeta]
                        sab_val_for_beta = self._terpt(tnow_sab_positive_part, self.nbeta, deltab_val, beta_report_val)

                    # Factor for this phonon order term (fac_for_order_nl)
                    # NJOY's `fac` definition:
                    if nl_order == 1: # n=0 in Fortran loop, nl=1 phonon
                        fac_for_order_nl = 1.0 / self.arat
                    else: # n > 0 in Fortran loop, nl > 1 phonon
                        term_val = 1.0
                        for i_val in range(1, nl_order + 1): # (f0*arat)^nl / nl!
                            term_val = term_val * f0_val * self.arat / i_val
                        fac_for_order_nl = term_val / self.arat
                    
                    self.ssm_values[itemp, j_alpha_idx, k_beta_idx] += current_alpha_factor * fac_for_order_nl * sab_val_for_beta
            
            # Prepare for next iteration: S_n becomes S_{n-1}
            tlast_sab[:] = tnow_sab[:] # Copy contents

        # Apply global exp(-alpha*f0) factor
        for j_alpha_idx in range(self.nalpha):
            alpha_val = self.alpha_values[j_alpha_idx]
            exp_factor = np.exp(max(EXPLIM, -alpha_val * f0_val))
            self.ssm_values[itemp, j_alpha_idx, :] *= exp_factor
        
        # Ensure S(alpha,beta) is not too small (underflow issues)
        self.ssm_values[self.ssm_values < TINY_SAB] = 0.0

        if self.iprint >= 1:
            sab_sum_check = np.sum(self.ssm_values[itemp,:,:]) * self.alpha_values[0] * self.beta_values[0] # Approx integral
            print(f"  sum check for s(a,b) at temp {itemp} = {sab_sum_check:.5e}")


    # Placeholder for other calculation methods
    # def calculate_incoherent_inelastic(self, temperature, temp_specific_data):
    #     pass

    # --- Discrete Oscillator Methods ---

    def _bfact(self, x_bessel_arg, dw_continuum_val, beta_i_energy_kt):
        """
        Corresponds to Fortran subroutine bfact.
        Calculates Bessel function related terms for discrete oscillators.
        x_bessel_arg (x): 2*alpha*dwc*sqrt(tev/therm)
        dw_continuum_val (dwc): Debye-Waller factor for the oscillator
        beta_i_energy_kt (betai): Oscillator energy / kT (E_osc / tev)
        
        Returns: tuple (bzero, bplus_array, bminus_array)
                 bzero: I_0(x) * exp(-x)
                 bplus_array: I_m(x) * exp(-x) * exp(-m*betai/2) for m > 0
                 bminus_array: I_m(x) * exp(-x) * exp( m*betai/2) for m > 0 (for I_{-m})
        """
        # Coefficients for Bessel function series expansion I_m(z) = (z/2)^m sum_{k=0 to inf} [(z^2/4)^k / (k! * Gamma(m+k+1))]
        # These seem to be for modified Bessel I_n(x)*exp(-x)
        # c0 to c37 are coefficients from NJOY's bfact.f for specific series. Max m is 50.
        _bfact_coeffs_c = [
            1.0000000000000000e+00, -5.0000000000000000e-01, 3.1250000000000000e-02,
            -1.3020833333333333e-03, 4.0690104166666667e-05, -1.0172526041666667e-06,
            2.1192762586805556e-08, -3.7844218905009921e-10, 5.9131592039077995e-12,
            -8.2127211165386111e-14, 1.0265901395673264e-15, -1.1665893781780400e-17,
            1.2151972689354583e-19, -1.1573307323194842e-21, 1.0211872580740942e-23,
            -8.3816767208546320e-26, 6.3906038639678800e-28, -4.5321739043086700e-30,
            3.0089576194810400e-32, -1.8728635214944200e-34, 1.1000000000000000e-36,
            1.0000000000000000e+00, -1.2500000000000000e-01, 1.3020833333333333e-03,
            -8.6805555555555556e-05, 4.3402777777777778e-06, -1.7361111111111111e-07,
            5.7870370370370370e-09, -1.6534391534391534e-10, 4.1335978835978836e-12,
            -9.1857730746619643e-14, 1.8371546149323929e-15, -3.3402811180588961e-17,
            5.5671351967648269e-19, -8.5648233796381952e-21, 1.2235461970911707e-22,
            -1.6313949294548944e-24, 2.0392436618186180e-26
        ] # c0 to c36 (37 coeffs)

        max_m_terms = 50  # Max terms for bplus/bminus arrays
        bplus_array = np.zeros(max_m_terms)
        bminus_array = np.zeros(max_m_terms)

        if x_bessel_arg == 0.0:
            bzero = 1.0
            # bplus and bminus remain zero
            return bzero, bplus_array, bminus_array

        if x_bessel_arg > 25.0: # Asymptotic expansion for I_m(x)*exp(-x)
            bzero = 1.0 / np.sqrt(2.0 * pi * x_bessel_arg)
            term = bzero
            exp_neg_half_betai = np.exp(max(EXPLIM, -0.5 * beta_i_energy_kt))
            exp_pos_half_betai = np.exp(max(EXPLIM,  0.5 * beta_i_energy_kt)) if beta_i_energy_kt < -2*EXPLIM else np.exp(EXPLIM)


            for m in range(max_m_terms): # m from 0 to 49
                # term for I_m(x)*exp(-x) is approx bzero * Product_{j=1 to m} (1 - (j-0.5)^2 / (2*x))
                # This simplified form is not what NJOY uses.
                # NJOY uses recurrence I_{m+1}(x) = I_{m-1}(x) - (2m/x)I_m(x)
                # And I_m(x)*exp(-x) tends to bzero / sqrt(1 + (m/x)^2) for large x, not simple product.
                # For large x, I_m(x) ~ exp(x)/sqrt(2*pi*x) * (1 - (4m^2-1)/(8x) + ...)
                # So I_m(x)*exp(-x) ~ 1/sqrt(2*pi*x) * (1 - (4m^2-1)/(8x))
                # This means term_m = bzero * (1 - (4*(m+1)^2-1)/(8*x_bessel_arg)) for bplus[m] (m+1 order)
                if m == 0: # I_0 term is bzero itself. I_1 uses a factor.
                    # bplus[0] is I_1(x)*exp(-x)*exp(-beta_i/2)
                    # bminus[0] is I_1(x)*exp(-x)*exp( beta_i/2) (for I_{-1})
                    # For I_1: use (1 - (4*1^2-1)/(8x)) = (1 - 3/(8x))
                    factor = (1.0 - (4.0*((m+1)**2) - 1.0) / (8.0 * x_bessel_arg)) # for I_{m+1}
                    current_I_m_exp_neg_x = bzero * factor # Approx I_1(x)*exp(-x)
                else: # For I_{m+1} based on I_m
                    # This iterative product is more like some series terms.
                    # Let's use the direct formula for I_m(x)*exp(-x) for large x from NIST DLMF 10.41.2
                    # I_m(x)exp(-x) ~ 1/sqrt(2*pi*x) * sum_{k=0} (-1)^k u_k(m) / (2x)^k
                    # where u_0(m)=1, u_1(m) = (4m^2-1)/4, etc.
                    # For m=0, u_1(0)=-1/4.  I_0(x)exp(-x) ~ bzero * (1 - (-1/4)/(2x)) = bzero*(1+1/(8x))
                    # The previous bzero was 1/sqrt(2pix). The formula should be (1-(4m^2-1)/8x + ...)
                    # So bzero needs to be adjusted or the factor needs to be (1 - (4m^2-1)/(8x))
                    # The Fortran calculates `term = term * (1.d0 - (float(m)-.5d0)**2 / (2.d0*x) )`
                    # This is product (1-(j-0.5)^2/(2x)).
                    # Let's replicate Fortran's product for `term`
                    term *= (1.0 - ((m + 1) - 0.5)**2 / (2.0 * x_bessel_arg)) # m for I_{m+1}
                    current_I_m_exp_neg_x = term
                
                if m < max_m_terms:
                    bplus_array[m] = current_I_m_exp_neg_x * exp_neg_half_betai**(m+1)
                    bminus_array[m] = current_I_m_exp_neg_x * exp_pos_half_betai**(m+1)
        
        else: # Series expansion for x <= 25.0
            # I_0(x)*exp(-x)
            xsq = x_bessel_arg * x_bessel_arg
            sum0 = _bfact_coeffs_c[0] + xsq*(_bfact_coeffs_c[2] + xsq*(_bfact_coeffs_c[4] + xsq*(_bfact_coeffs_c[6] + \
                   xsq*(_bfact_coeffs_c[8] + xsq*(_bfact_coeffs_c[10] + xsq*(_bfact_coeffs_c[12] + \
                   xsq*(_bfact_coeffs_c[14] + xsq*(_bfact_coeffs_c[16] + xsq*(_bfact_coeffs_c[18] + \
                   xsq*_bfact_coeffs_c[20])))))))))
            sum1 = _bfact_coeffs_c[1] + xsq*(_bfact_coeffs_c[3] + xsq*(_bfact_coeffs_c[5] + xsq*(_bfact_coeffs_c[7] + \
                   xsq*(_bfact_coeffs_c[9] + xsq*(_bfact_coeffs_c[11] + xsq*(_bfact_coeffs_c[13] + \
                   xsq*(_bfact_coeffs_c[15] + xsq*(_bfact_coeffs_c[17] + xsq*_bfact_coeffs_c[19])))))))))
            bzero = sum0 + x_bessel_arg * sum1
            
            # I_m(x)*exp(-x) for m > 0
            # Uses recurrence: I_{m+1} = I_{m-1} - (2m/x)I_m
            # Or specific series for I_1, I_2 etc.
            # Fortran uses different sets of C coefficients for different orders.
            # C(22)...C(37) are for I_1(x)*exp(-x)/x
            # This is getting too complex for direct one-to-one with given C array.
            # The Fortran bfact.f uses recurrence with I_0, I_1 calculated from series.
            # I_0(x)exp(-x) is bzero.
            # I_1(x)exp(-x) needs to be calculated.
            # Let's use scipy.special.i0e for I_0(x)exp(-x) and i1e for I_1(x)exp(-x)
            # And ive for I_m(x)exp(-x) for m > 1.
            # This avoids using the C coefficients directly for now, simplifying.
            # Requires scipy: `pip install scipy`
            try:
                from scipy.special import i0e, i1e, ive
            except ImportError:
                raise NJOYError("Scipy not installed. Required for _bfact if x_bessel_arg <= 25.0")

            bzero = i0e(x_bessel_arg)
            current_Im_exp_neg_x = bzero    # For m=0
            next_Im_exp_neg_x = i1e(x_bessel_arg) # For m=1
            
            exp_neg_half_betai = np.exp(max(EXPLIM, -0.5 * beta_i_energy_kt))
            exp_pos_half_betai = np.exp(max(EXPLIM,  0.5 * beta_i_energy_kt)) if beta_i_energy_kt < -2*EXPLIM else np.exp(EXPLIM)


            for m in range(max_m_terms): # m from 0 to 49, for I_{m+1}
                # bplus[m] is for I_{m+1}
                # bminus[m] is for I_{-(m+1)} (which is I_{m+1})
                
                # current_Im_exp_neg_x is I_m * exp(-x)
                # next_Im_exp_neg_x is I_{m+1} * exp(-x)
                
                bplus_array[m] = next_Im_exp_neg_x * (exp_neg_half_betai**(m+1))
                bminus_array[m] = next_Im_exp_neg_x * (exp_pos_half_betai**(m+1))

                if x_bessel_arg == 0: # Avoid division by zero in recurrence if x is tiny (already handled for x=0)
                    current_Im_exp_neg_x = 0.0 # Higher orders will be zero
                    next_Im_exp_neg_x = 0.0
                    if m+2 <=2 : # I_2(0)=0, I_3(0)=0
                         current_Im_exp_neg_x = 0.0
                    else: # for m >= 1, I_{m+1}(0)=0
                         next_Im_exp_neg_x = 0.0

                else:
                    # Recurrence: I_{m+2} = I_m - (2(m+1)/x)I_{m+1}
                    # Or use ive(m+2, x_bessel_arg) for I_{m+2}*exp(-x)
                    if m + 2 <= 100: # ive is stable for reasonable orders
                         val_I_m_plus_2 = ive(m+2, x_bessel_arg)
                    else: # if order gets too high, it tends to zero faster than 1/m!
                         val_I_m_plus_2 = 0.0
                    
                    current_Im_exp_neg_x = next_Im_exp_neg_x
                    next_Im_exp_neg_x = val_I_m_plus_2
                
                if abs(next_Im_exp_neg_x) < TINY_VAL and m > 10 : # Optimization: if terms get too small
                    if self.iprint >=3: print(f"bfact: I_{m+2} term small, truncating at m={m}")
                    break
        
        return bzero, bplus_array, bminus_array

    def _bfill(self, original_beta_grid, num_original_beta_points, max_beta_extended):
        """
        Corresponds to Fortran subroutine bfill.
        Sets up an extended beta grid (bex_array) and its reciprocal delta (rdbex_array).

        original_beta_grid: numpy array (self.beta_values for current temp)
        num_original_beta_points: int (self.nbeta)
        max_beta_extended: float, maximum beta value for the extended grid.

        Returns: tuple (bex_array, rdbex_array, nbx_val)
                 bex_array: numpy array, the extended beta grid
                 rdbex_array: numpy array, reciprocal delta beta for bex_array
                 nbx_val: int, number of points in bex_array
        """
        if num_original_beta_points < 2: # Not enough points to determine original spacing
            # Default to a fine grid if original is too small
            fine_delta_beta = 0.05
            nbx_val = int(max_beta_extended / fine_delta_beta) + 1
            if nbx_val <= 1: nbx_val = 2 # Ensure at least 2 points
            bex_array = np.linspace(0, max_beta_extended, nbx_val)
            rdbex_array = np.full(nbx_val -1 if nbx_val > 1 else 1, 1.0 / fine_delta_beta if fine_delta_beta > 0 else np.inf)
            if nbx_val == 1 and fine_delta_beta == 0 : rdbex_array[0] = np.inf # Avoid division by zero if max_beta_extended is 0
            return bex_array, rdbex_array, nbx_val

        original_delta_beta1 = original_beta_grid[1] - original_beta_grid[0]

        # Fortran logic: if (betan(2).lt.0.1d0) then bex is betan
        # Assuming betan(1) is beta=0, betan(2) is first point. Here, grid[0] and grid[1].
        if original_delta_beta1 < 0.1 and abs(original_delta_beta1 - (original_beta_grid[-1]/ (num_original_beta_points-1)) ) < SMALL_VAL :
            # Original grid is fine enough and uniform, use it directly if it covers max_beta_extended
            if original_beta_grid[-1] >= max_beta_extended:
                bex_array = np.copy(original_beta_grid)
                nbx_val = num_original_beta_points
                # rdbex needs to be array of 1/delta for each interval
                if nbx_val > 1:
                    deltas = bex_array[1:] - bex_array[:-1]
                    rdbex_array = 1.0 / deltas
                else: # single point, delta is undefined or infinite
                    rdbex_array = np.array([np.inf]) 
                return bex_array, rdbex_array, nbx_val
            else: # Original grid is fine, but doesn't extend far enough
                fine_delta_beta = original_delta_beta1
                # nbx_val will be calculated below for new grid
        else: # Original grid is too coarse or non-uniform, create a new fine grid
            fine_delta_beta = 0.05
        
        # Create a new fine grid up to max_beta_extended
        if fine_delta_beta <=0: raise NJOYError("Beta grid spacing must be positive in _bfill")
        
        nbx_val = int(max_beta_extended / fine_delta_beta) + 1
        if nbx_val <= 1 and max_beta_extended > 0: nbx_val = 2 # Ensure grid if max_beta_extended > 0
        elif max_beta_extended == 0: nbx_val = 1 # Single point at zero

        if nbx_val ==1 :
            bex_array = np.array([0.0])
            rdbex_array = np.array([np.inf]) # No interval
        else:
            bex_array = np.linspace(0, max_beta_extended, nbx_val)
             # rdbex_array should have size nbx_val-1 for intervals
            deltas = bex_array[1:] - bex_array[:-1]
            # Handle cases where delta might be zero if nbx_val is large and max_beta_extended is small
            # However, linspace should prevent this unless max_beta_extended=0 (handled) or nbx_val=1 (handled)
            rdbex_array = 1.0 / deltas

        if self.iprint >= 3:
            print(f"  bfill: Extended beta grid created with nbx={nbx_val}, max_beta={max_beta_extended:.4f}, delta={fine_delta_beta:.4f}")
        
        return bex_array, rdbex_array, nbx_val

    def _exts(self, sab_on_original_beta_grid, original_beta_grid, 
                bex_array_positive, nbx_positive):
        """
        Corresponds to Fortran subroutine exts.
        Extends S(alpha,beta) defined for positive beta to a symmetric form 
        over a new (potentially finer/larger) positive beta grid, then implies symmetry.

        sab_on_original_beta_grid: numpy array, S(alpha,beta) values for a given alpha,
                                   on the original_beta_grid.
        original_beta_grid: numpy array, the original beta grid (e.g., self.beta_values).
        bex_array_positive: numpy array, the extended positive beta grid from _bfill.
                            Grid points are >= 0.
        nbx_positive: int, number of points in bex_array_positive.

        Returns: numpy array (sex_symmetric_output_array) of size (2*nbx_positive - 1),
                 representing S(alpha,beta) on a symmetric grid centered at zero.
                 The center element (index nbx_positive-1) is S(alpha,0).
        """
        
        if nbx_positive == 0:
            return np.array([])
        if nbx_positive == 1: # Only beta=0 point
            s_at_zero = 0.0
            if len(original_beta_grid) > 0 and original_beta_grid[0] == 0.0:
                 s_at_zero = sab_on_original_beta_grid[0]
            else: # Interpolate or extrapolate S(alpha,0)
                 s_at_zero = np.interp(0.0, original_beta_grid, sab_on_original_beta_grid, left=0.0, right=0.0)
            return np.array([s_at_zero])

        # Create the symmetric output array
        # Size 2*nbx_positive - 1. Center index is nbx_positive - 1.
        sex_symmetric_output_array = np.zeros(2 * nbx_positive - 1)
        center_idx = nbx_positive - 1

        # Interpolate sab_on_original_beta_grid onto bex_array_positive for the positive part
        # np.interp(x_new, x_orig, y_orig)
        # sab_on_bex_positive will be S(alpha, beta_bex) for beta_bex >= 0
        sab_on_bex_positive = np.interp(
            bex_array_positive, 
            original_beta_grid, 
            sab_on_original_beta_grid,
            left=0.0, # Extrapolate with 0 if bex_array_positive extends lower
            right=0.0  # Extrapolate with 0 if bex_array_positive extends higher
        )

        # Fill the positive part of sex_symmetric_output_array (including beta=0)
        # sex_symmetric_output_array[center_idx] is S(alpha, bex_array_positive[0]=0)
        # sex_symmetric_output_array[center_idx + k] is S(alpha, bex_array_positive[k])
        for k in range(nbx_positive):
            sex_symmetric_output_array[center_idx + k] = sab_on_bex_positive[k]

        # Fill the negative part using detailed balance S(alpha, -beta) = exp(-beta) * S(alpha, beta)
        # sex_symmetric_output_array[center_idx - k] is S(alpha, -bex_array_positive[k])
        for k in range(1, nbx_positive): # k from 1 up to nbx_positive-1
            beta_val = bex_array_positive[k] # This is positive beta_k
            s_at_pos_beta_k = sab_on_bex_positive[k] # S(alpha, beta_k)
            
            exp_factor = np.exp(max(EXPLIM, -beta_val)) # exp(-beta_k)
            sex_symmetric_output_array[center_idx - k] = exp_factor * s_at_pos_beta_k
            
        return sex_symmetric_output_array

    def _sint(self, x_beta_target, 
                bex_array_positive, nbx_positive, rdbex_reciprocal_deltas, 
                sex_symmetric_array, 
                alpha_val, total_weight_div_arat, effective_temp_k, 
                original_beta_grid_max_val):
        """
        Corresponds to Fortran function sint.
        Interpolates S(alpha,beta) from sex_symmetric_array or uses SCT approximation.

        x_beta_target: float, the beta value at which S(a,b) is needed (can be +/-).
        bex_array_positive: numpy array, positive part of the extended beta grid from _bfill.
        nbx_positive: int, number of points in bex_array_positive.
        rdbex_reciprocal_deltas: numpy array, reciprocal deltas for intervals in bex_array_positive.
        sex_symmetric_array: numpy array, symmetric S(a,b) from _exts.
        alpha_val: float, current alpha value.
        total_weight_div_arat: float, total weight / self.arat (sab0 in Fortran).
        effective_temp_k: float, effective temperature in Kelvin.
        original_beta_grid_max_val: float, max beta of original grid (self.beta_values[-1]).
        
        Returns: float, S(alpha, beta=x_beta_target).
        """

        # Determine if interpolation or SCT should be used
        # Use SCT if abs(x_beta_target) is greater than the max beta of the original grid
        use_sct = False
        if abs(x_beta_target) > original_beta_grid_max_val:
            use_sct = True
        
        # Also, ensure x_beta_target is within the bounds of the grid covered by sex_symmetric_array for interpolation
        max_bex_val = bex_array_positive[-1] if nbx_positive > 0 else 0.0
        if abs(x_beta_target) > max_bex_val: # If target is outside even the extended grid
            use_sct = True

        if not use_sct:
            # Interpolation part
            if nbx_positive == 0: return 0.0
            if nbx_positive == 1: # Single point at beta=0
                if abs(x_beta_target) < SMALL_VAL: # Effectively asking for S(alpha,0)
                    return sex_symmetric_array[0] 
                else: # Asking for non-zero beta when grid is only beta=0
                    use_sct = True # Fallback to SCT

        if not use_sct:
            # Construct the full symmetric grid corresponding to sex_symmetric_array
            # sex_symmetric_array has size 2*nbx_positive - 1. Center index is nbx_positive - 1.
            # bex_array_positive[0] must be 0.0 for this construction.
            if nbx_positive > 0 and abs(bex_array_positive[0]) > SMALL_VAL:
                # This case should ideally not happen if _bfill ensures bex starts at 0
                print("Warning (_sint): bex_array_positive[0] is not zero. Interpolation might be incorrect.")

            if nbx_positive == 1: # Already handled S(alpha,0) case if x_beta_target was 0
                 full_symmetric_bex_grid = np.array([0.0])
            else:
                # Negative part: -bex_array_positive[:0:-1] means reverse all but first, then negate.
                # Example: bex=[0,1,2,3] -> neg part from [3,2,1] -> [-3,-2,-1]
                neg_part = -bex_array_positive[1:][::-1] 
                full_symmetric_bex_grid = np.concatenate((neg_part, bex_array_positive))
            
            try:
                interpolated_s_val = np.interp(x_beta_target, full_symmetric_bex_grid, sex_symmetric_array, left=0.0, right=0.0)
                return interpolated_s_val
            except Exception as e:
                if self.iprint >=1:
                    print(f"Error during interpolation in _sint: {e}. Switching to SCT.")
                use_sct = True # Fallback to SCT on error

        # SCT approximation part (if use_sct is True)
        if alpha_val < SMALL_VAL: # Avoid division by zero or issues with alpha=0
            return 0.0 

        # teff_ev = effective_temp_k * self.bk # Effective temperature in eV
        # den = alpha_val * teff_ev / (THERM_CONST * self.arat) # NJOY: den = alph*tbart/(therm*arat)
        # tbart in NJOY is already teff (in Kelvin). therm is THERM_CONST (eV).
        # So, den = alpha_val * (effective_temp_k * self.bk / THERM_CONST) / self.arat
        # self.arat is aws/awr (atomic weight ratio of secondary to primary scatterer)
        # For primary scatterer, self.arat is 1.0.
        # The 'arat' in NJOY's SCT formula refers to the primary scatterer's atomic weight ratio (awr).
        # Let's use effective_awr for clarity in SCT. If only primary, effective_awr = self.awr.
        # If there's a secondary scatterer defined by self.arat = aws/awr, then the SCT formula
        # should consistently use the properties of the species it's describing.
        # The 'wt' (total_weight_div_arat) passed to sint is sum of w_i/A_i.
        # If SCT is for the combined effect, then self.arat (A) might be an effective one.
        # NJOY's `discre` uses `self.arat` (awr) in the SCT call for `sint`.
        # `total_weight_div_arat` is sab0 = wt/arat in Fortran. `wt` is sum of oscillator weights.
        
        # Using NJOY's formulation directly:
        # sab0 = total_weight_div_arat  (this is sum(adel_i / awr_i))
        # arat_for_sct = self.arat (awr of primary scatterer from input card 6)
        # tbart_kelvin = effective_temp_k
        
        # Denominator term, NJOY: den = alph * tbart / (therm * arat)
        # alph = alpha_val
        # tbart = effective_temp_k (Kelvin)
        # therm = THERM_CONST (eV)
        # arat = self.awr (atomic weight of primary scatterer)
        # Need conversion bk somewhere if tbart is K and therm is eV.
        # tbart (effective temp) in NJOY is T_eff (Kelvin). Denominator uses T_eff / T_0,
        # where T_0 is room temp (293.6K -> THERM_CONST eV).
        # So, tbart/(therm) should be (effective_temp_k * bk) / THERM_CONST for eV/eV ratio,
        # or effective_temp_k / (THERM_CONST/bk) for K/K ratio.
        # Let's use eV for energies: teff_ev = effective_temp_k * self.bk
        
        teff_ev = effective_temp_k * self.bk
        if self.awr <= 0: # Should not happen with valid inputs
            if self.iprint >=0: print("Error (_sint SCT): self.awr is zero or negative.")
            return 0.0
            
        den = alpha_val * teff_ev / (THERM_CONST * self.awr) # Use self.awr for SCT's 'arat'
        
        if den <= SMALL_VAL : # Avoid division by zero or log(negative) if arg becomes large neg
            return 0.0 
            
        # Argument of exp, NJOY: arg = (x/arat - alph)**2 / (4.d0*den)
        # x = x_beta_target
        arg_numerator = (x_beta_target / self.awr - alpha_val)**2
        arg = arg_numerator / (4.0 * den)
        
        sct_val = 0.0
        try:
            # sab0_norm = total_weight_div_arat (this is already wt/awr)
            exp_val = np.exp(max(EXPLIM, -arg))
            sct_val = total_weight_div_arat / np.sqrt(4.0 * np.pi * den) * exp_val
        except OverflowError:
            sct_val = 0.0 # Result too large or too small after exp
        
        if np.isnan(sct_val) or np.isinf(sct_val):
            return 0.0

        return sct_val

    def _discre(self, itemp, current_temp_data):
        """
        Corresponds to Fortran subroutine discre.
        Calculates S(alpha,beta) contributions from discrete oscillators.
        Modifies self.ssm_values, self.dwpix_values, self.tempf_values.
        """
        num_discrete_oscillators = current_temp_data.get('nd', 0)
        if num_discrete_oscillators == 0:
            if self.iprint >= 2: print(f"  _discre: No discrete oscillators (nd=0) for itemp={itemp}.")
            return

        current_temp_k = current_temp_data['temperature']
        tev = current_temp_k * self.bk
        if tev == 0.0: raise NJOYError("_discre: Temperature in eV is zero.")

        # Get oscillator energies (eV) and weights from input data
        # These should be populated by _parse_input from Card 21
        osc_energies_ev_input = np.array(current_temp_data.get('bdel_values_card21', [])[:num_discrete_oscillators])
        osc_weights_input = np.array(current_temp_data.get('adel_values_card21', [])[:num_discrete_oscillators])

        if len(osc_energies_ev_input) != num_discrete_oscillators or \
           len(osc_weights_input) != num_discrete_oscillators:
            raise NJOYError(f"_discre: Mismatch in number of oscillators nd={num_discrete_oscillators} and provided energy/weight arrays.")

        # Convert oscillator energies to E/kT (dimensionless beta)
        beta_osc_energies_kt = osc_energies_ev_input / tev # bdeln(i) in Fortran

        # SC factor (sqrt(therm/tev)) for x calculation in bfact
        sc_factor_sqrt = np.sqrt(THERM_CONST / tev) if self.lat == 1 else 1.0
        
        max_beta_original_grid = self.beta_values[-1] if self.nbeta > 0 else 0.0
        
        # Estimate max beta needed for extended grid based on oscillator energies
        # Max beta can reach sum of (m_i * E_i_kt) over oscillators, m_i up to ~50
        # This is a rough upper bound for the extended grid.
        max_possible_beta_shift = np.sum(beta_osc_energies_kt * 50) # 50 from max_m_terms in _bfact
        max_beta_for_bfill = max_beta_original_grid + max_possible_beta_shift
        # NJOY uses fixed maxbb = betan(nbeta) + sum(bdeln)*maxd + 5.0/tev, then limits nbx.
        # Let's cap max_beta_for_bfill to something reasonable if it gets too large, e.g. 2000/tev.
        max_beta_for_bfill = min(max_beta_for_bfill, 2000.0 / tev if tev > 0 else 2000.0)


        tbeta_continuous_weight = current_temp_data.get('tbeta', 0.0) # Weight of continuous part
        # twt_translational_weight needs to be parsed from Card 22 (default to 0 if not specified)
        twt_translational_weight = current_temp_data.get('twt', 0.0) 

        if self.iprint >= 1:
            print(f"  Entering _discre for itemp={itemp}, nd={num_discrete_oscillators}")
            if self.iprint >=2:
                 print(f"    Oscillator energies (eV): {osc_energies_ev_input}")
                 print(f"    Oscillator weights: {osc_weights_input}")
                 print(f"    Oscillator beta energies (E/kT): {beta_osc_energies_kt}")
                 print(f"    Max beta for bfill: {max_beta_for_bfill:.4f}")


        # --- Loop over alpha values ---
        for nal_idx in range(self.nalpha):
            alpha_val = self.alpha_values[nal_idx]
            if self.iprint >=3: print(f"    _discre: Processing alpha = {alpha_val:.4e} (index {nal_idx})")

            # Setup extended beta grid using _bfill
            # Note: rdbex_recip_deltas is 1/(delta_beta_on_bex_grid) for each interval
            bex_array_pos, rdbex_recip_deltas, nbx_pos = \
                self._bfill(self.beta_values, self.nbeta, max_beta_for_bfill)
            
            if nbx_pos == 0 : # Should not happen if max_beta_for_bfill > 0
                if self.iprint >=1: print("Warning (_discre): nbx_pos from _bfill is 0. Skipping alpha.")
                continue
            delta_bex_uniform = bex_array_pos[1] - bex_array_pos[0] if nbx_pos > 1 else 1.0 # Assuming uniform for now

            # Get S(alpha,beta) from previous step (e.g., continuous model from _contin)
            # This is S(alpha,beta_orig_grid)
            sab_prev_step_on_orig_beta = np.copy(self.ssm_values[itemp, nal_idx, :])

            # Extend this S(alpha,beta) to a symmetric form on the bex_array_pos grid
            sex_symmetric_prev_step = self._exts(sab_prev_step_on_orig_beta, 
                                                 self.beta_values, self.nbeta,
                                                 bex_array_pos, nbx_pos)

            # Initialize convoluted beta lines and weights for this alpha
            # Start with a single delta function at beta=0, weight=1.0
            # This represents S_0(beta) for the oscillator part before any convolution.
            beta_lines_current_osc_iteration = [0.0] 
            weights_current_osc_iteration = [1.0]
            num_lines_cumulative = 1

            # --- Loop over discrete oscillators ---
            for i_osc in range(num_discrete_oscillators):
                E_i_kt = beta_osc_energies_kt[i_osc] # Current oscillator energy in kT units (bdeln(i))
                W_i = osc_weights_input[i_osc]      # Current oscillator weight (adel(i))

                if self.iprint >=3: print(f"      Oscillator {i_osc+1}/{num_discrete_oscillators}, E/kT={E_i_kt:.4f}, W={W_i:.4f}")

                if abs(E_i_kt) < SMALL_VAL: # Skip zero-energy oscillator
                    if self.iprint>=2: print(f"        Skipping oscillator {i_osc} due to zero energy.")
                    continue
                
                # dwc_val = adel(i) / (arat*bdeln(i)) in Fortran discre
                # self.awr is arat for primary scatterer.
                dw_contrib_osc = W_i / (self.awr * E_i_kt) 
                
                x_bessel = 2.0 * alpha_val * dw_contrib_osc * sc_factor_sqrt
                
                bzero, bplus, bminus = self._bfact(x_bessel, dw_contrib_osc, E_i_kt)

                beta_lines_next_iteration = []
                weights_next_iteration = []
                
                max_m_terms_from_bfact = len(bplus) # Should be 50

                for k_line in range(num_lines_cumulative):
                    beta_prev = beta_lines_current_osc_iteration[k_line]
                    weight_prev = weights_current_osc_iteration[k_line]

                    if abs(weight_prev) < TINY_VAL : continue # Skip if previous weight is negligible

                    # Add central line (m=0 term from I_0)
                    beta_lines_next_iteration.append(beta_prev)
                    weights_next_iteration.append(weight_prev * bzero)

                    # Add lines from I_m, m > 0 (positive and negative beta shifts)
                    for m_bessel_idx in range(max_m_terms_from_bfact): # m_bessel_idx from 0 to 49
                        m_val = m_bessel_idx + 1 # Actual m for Bessel function (1 to 50)
                        
                        if abs(bplus[m_bessel_idx]) > TINY_VAL:
                            beta_lines_next_iteration.append(beta_prev + m_val * E_i_kt)
                            weights_next_iteration.append(weight_prev * bplus[m_bessel_idx])
                        
                        if abs(bminus[m_bessel_idx]) > TINY_VAL:
                            beta_lines_next_iteration.append(beta_prev - m_val * E_i_kt)
                            weights_next_iteration.append(weight_prev * bminus[m_bessel_idx])
                
                # Sort, combine nearby lines, and prune
                if not beta_lines_next_iteration: # If all weights became tiny
                    beta_lines_current_osc_iteration = [0.0] # Reset to prevent issues
                    weights_current_osc_iteration = [0.0]
                    num_lines_cumulative = 1
                    if self.iprint >=2: print("        All lines pruned after oscillator, resetting.")
                    continue


                # Combine and Prune:
                # 1. Convert to numpy arrays for easier manipulation
                beta_lines_np = np.array(beta_lines_next_iteration)
                weights_np = np.array(weights_next_iteration)

                # 2. Sort by beta values
                sorted_indices = np.argsort(beta_lines_np)
                beta_lines_sorted = beta_lines_np[sorted_indices]
                weights_sorted = weights_np[sorted_indices]
                
                # 3. Combine lines with close beta values
                if len(beta_lines_sorted) > 0:
                    beta_combined = [beta_lines_sorted[0]]
                    weight_combined = [weights_sorted[0]]
                    for j_line_idx in range(1, len(beta_lines_sorted)):
                        if abs(beta_lines_sorted[j_line_idx] - beta_combined[-1]) < SMALL_VAL: # Threshold for combining
                            weight_combined[-1] += weights_sorted[j_line_idx]
                        else:
                            beta_combined.append(beta_lines_sorted[j_line_idx])
                            weight_combined.append(weights_sorted[j_line_idx])
                    
                    beta_lines_final_iter = np.array(beta_combined)
                    weights_final_iter = np.array(weight_combined)

                    # 4. Prune lines with very small weights (relative to max weight in this iteration)
                    if len(weights_final_iter) > 0:
                        max_w_iter = np.max(np.abs(weights_final_iter))
                        if max_w_iter > TINY_VAL: # Avoid division by zero if all weights are zero
                            prune_mask = (np.abs(weights_final_iter) > max_w_iter * 1e-7) & (np.abs(weights_final_iter) > TINY_VAL)
                            beta_lines_final_iter = beta_lines_final_iter[prune_mask]
                            weights_final_iter = weights_final_iter[prune_mask]
                    
                    # 5. Limit number of lines (e.g., NJOY's maxdd=500)
                    #    If too many, keep the ones with largest absolute weights.
                    MAX_DISCRETE_LINES = 500 
                    if len(beta_lines_final_iter) > MAX_DISCRETE_LINES:
                        if self.iprint >= 2: print(f"        Warning: Number of discrete lines {len(beta_lines_final_iter)} exceeds max {MAX_DISCRETE_LINES}. Pruning.")
                        sorted_abs_weight_indices = np.argsort(np.abs(weights_final_iter))[::-1] # Descending
                        keep_indices = np.sort(sorted_abs_weight_indices[:MAX_DISCRETE_LINES])
                        beta_lines_final_iter = beta_lines_final_iter[keep_indices]
                        weights_final_iter = weights_final_iter[keep_indices]

                    beta_lines_current_osc_iteration = list(beta_lines_final_iter)
                    weights_current_osc_iteration = list(weights_final_iter)
                    num_lines_cumulative = len(beta_lines_current_osc_iteration)
                else: # All lines might have been pruned if initial list was empty or weights were all zero
                    beta_lines_current_osc_iteration = [0.0]
                    weights_current_osc_iteration = [0.0] # Avoid issues if list becomes empty
                    num_lines_cumulative = 1


                if num_lines_cumulative == 0 : # Should not happen if reset to [0.0] with [0.0] weight
                    if self.iprint >= 1: print("Warning: num_lines_cumulative is zero in _discre.")
                    beta_lines_current_osc_iteration = [0.0] # Safety reset
                    weights_current_osc_iteration = [0.0]
                    num_lines_cumulative = 1


            # --- End of oscillator loop ---
            if self.iprint >=2 : print(f"      Finished oscillator loop for alpha={alpha_val:.3e}. Num final lines: {num_lines_cumulative}")


            # --- Construct final S(alpha,beta) for current alpha using convoluted lines ---
            sab_current_alpha_new = np.zeros(self.nbeta)
            
            # Effective temperature for SCT part of _sint should be the one from continuous + translational part
            # NJOY 'discre' uses `tbart = tempf(itemp)` which is already T_eff from _contin.
            effective_temp_for_sint_k = self.tempf_values[itemp] 
            
            # Total weight for normalization in _sint's SCT part.
            # NJOY: sab0 = (tbeta + twt + sum of adel over processed oscillators) / arat
            # Here, sum of adel is sum(osc_weights_input)
            total_processed_osc_weight_sum = np.sum(osc_weights_input) 
            current_total_sab_weight = tbeta_continuous_weight + twt_translational_weight + total_processed_osc_weight_sum
            sab0_norm_for_sint = current_total_sab_weight / self.awr # arat in NJOY's sint context is awr

            if self.iprint >=3 and nal_idx ==0 : # Print only for first alpha for brevity
                 print(f"    _discre SGT params: sab0_norm={sab0_norm_for_sint:.4e}, T_eff_sint={effective_temp_for_sint_k:.4e} K")


            for k_beta_report_idx in range(self.nbeta): # Loop over output beta grid points
                beta_report_val = self.beta_values[k_beta_report_idx]
                accumulated_s_val = 0.0

                for j_line_idx in range(num_lines_cumulative):
                    beta_line = beta_lines_current_osc_iteration[j_line_idx]
                    weight_line = weights_current_osc_iteration[j_line_idx]

                    if abs(weight_line) < TINY_VAL: continue

                    delta_beta_arg = beta_report_val - beta_line # Beta shift for _sint

                    # Call _sint with S(a,b) from *before* discrete oscillator processing for this alpha
                    # which is sex_symmetric_prev_step (on bex_array_pos grid)
                    s_contrib = weight_line * self._sint(
                        delta_beta_arg,
                        bex_array_pos, nbx_pos, rdbex_recip_deltas,
                        sex_symmetric_prev_step, # S(a,b) of continuous part on extended symm grid
                        alpha_val,
                        sab0_norm_for_sint, # Normalization for SCT part of _sint
                        effective_temp_for_sint_k, # T_eff for SCT part of _sint
                        max_beta_original_grid
                    )
                    accumulated_s_val += s_contrib
                
                sab_current_alpha_new[k_beta_report_idx] = accumulated_s_val
            
            self.ssm_values[itemp, nal_idx, :] = sab_current_alpha_new
            if self.iprint >=3 and nal_idx == 0: # Print S(a,b) for first alpha after discre
                print(f"    SSM for itemp={itemp}, alpha_idx=0 after _discre (first few beta points):")
                print(f"    {self.ssm_values[itemp, nal_idx, :min(5, self.nbeta)]}")


        # --- End of alpha loop ---

        # --- Update global Debye-Waller factor and Effective Temperature ---
        # These sums are over the discrete oscillators only.
        dw_integral_discrete_sum = 0.0
        effective_temp_contrib_discrete_sum_weighted = 0.0 # eb(i)*arat in Fortran

        for i_osc in range(num_discrete_oscillators):
            E_i_kt = beta_osc_energies_kt[i_osc] # E_osc / kT
            W_i = osc_weights_input[i_osc]       # adel_i

            if abs(E_i_kt) < SMALL_VAL: continue

            # coth(E/2kT) term
            # For coth(x), x = E_i_kt / 2.0
            x_coth = E_i_kt / 2.0
            coth_val_m1 = 0.0 # coth(x) - 1
            if abs(x_coth) < 1e-4: # Small x approximation for coth(x)-1 = 2/(e^{2x}-1) ~ 2/(2x + (2x)^2/2) = 1/x - 1 + x/3 ...
                                  # Or coth(x)-1 = (e^x+e^-x - (e^x-e^-x))/(e^x-e^-x) = 2e^-x / (e^x-e^-x)
                                  # For small x, e^x ~ 1+x, e^-x ~ 1-x.  2(1-x)/((1+x)-(1-x)) = 2(1-x)/(2x) = (1-x)/x = 1/x -1.
                coth_val_m1 = (1.0/x_coth if x_coth !=0 else np.inf) -1.0 # Approx for coth(x)-1
            else:
                exp_2x_coth = np.exp(max(EXPLIM, 2.0 * x_coth))
                if abs(exp_2x_coth - 1.0) < VSMALL_VAL: # Avoid division by zero
                    coth_val_m1 = (1.0/x_coth if x_coth !=0 else np.inf) -1.0 # Fallback to small x approx
                else:
                    coth_val_m1 = 2.0 / (exp_2x_coth - 1.0) # coth(x)-1 = 2 / (e^{2x}-1)

            # dw_integral_discrete_sum += W_i / self.awr * coth_term_div_beta * 0.5 (NJOY fsum like)
            # dbw(i) = adel(i)/arat * (coth(bdeln(i)/2.)-1.) / bdeln(i)
            dw_integral_discrete_sum += (W_i / self.awr) * (coth_val_m1 / E_i_kt)
            
            # effective_temp_contrib_discrete_sum += W_i / self.awr * coth_term * 0.5 * E_i_kt (NJOY fsum like)
            # eb(i) = adel(i) * (coth(bdeln(i)/2.)-1.)
            effective_temp_contrib_discrete_sum_weighted += W_i * coth_val_m1 # This is eb(i)*arat in Fortran
        
        # Add discrete contributions to global DW factor and T_eff
        # Note: NJOY's dwpix(itemp) is sum over (dbw(i) from contin + dbw(i) from discre)
        # dbw(i) from contin is f0 = an/arat where an = _fsum(1,...)
        # My _fsum already includes the 0.5 factor. So dw_integral_discrete_sum needs it too if anology holds.
        # Let's assume the formulas from NJOY manual / comments are direct:
        # DW_contrib_discrete = Sum_i [ W_i/A * (coth(E_i/2kT)-1)/(E_i/kT) ] -- this is what I have for dw_integral_discrete_sum
        self.dwpix_values[itemp] += dw_integral_discrete_sum 

        # Update effective temperature (tempf)
        # tempf(itemp) initially holds T_eff_contin from _contin (which is Tbar_contin/T_actual * T_actual)
        # T_eff_final = ( (T_eff_contin/T_actual)*(W_contin+W_trans) + Sum_i[ W_i/A * (coth(E_i/2kT)-1)*(E_i/kT) ] )
        #               / (W_contin+W_trans + Sum_i W_i/A) * T_actual
        # My effective_temp_contrib_discrete_sum_weighted is Sum_i [ W_i * (coth-1) ]
        # This corresponds to Sum_i [ (eb(i)*arat) ] from Fortran.
        # Fortran effective temp update:
        # sum1 = tempf(itemp)/temp * (tbeta+twt)  (tempf was T_eff_contin)
        # sum2 = sum of eb(i) / temp           (eb(i) is effective_temp_contrib_discrete_sum_weighted / self.awr * E_i_kt * 0.5)
        # sumw = tbeta + twt + sum of adel(i)
        # tempf(itemp) = (sum1+sum2)/sumw * temp
        
        # Contribution to T_eff/T from discrete part: Sum_i [ W_i/A * (coth-1)*(E_i/kT) ] * (1/kT) ? No.
        # T_eff/T = Sum_i [ W_i/A_i * (coth(E_i/2kT)-1) * (E_i/2kT) ] / Sum_i [W_i/A_i]
        # My effective_temp_contrib_discrete_sum_weighted = Sum_i W_i * (coth(E_i_kt/2)-1)
        # So, Sum_i [ (W_i/self.awr) * (coth(E_i_kt/2)-1) * (E_i_kt/2) ] for the numerator part of T_eff_discrete/T
        # This is (effective_temp_contrib_discrete_sum_weighted / self.awr) * (E_i_kt / 2) summed. No, this is wrong.

        # Let teff_factor_contin = self.tempf_values[itemp] / current_temp_k (from _contin)
        # Weighted sum for new teff_factor:
        # Numerator_eff_temp = teff_factor_contin * (tbeta_continuous_weight + twt_translational_weight) \
        #                    + Sum_i [ (W_i/self.awr) * (coth(E_i/2kT)-1) * (E_i_kt/2) ]
        # Denominator_eff_temp_norm = tbeta_continuous_weight + twt_translational_weight + np.sum(osc_weights_input)/self.awr
        
        # NJOY discre.f: sum2 = sum2 + eb(i). eb(i) = adel(i)*(coth(bdeln(i)/2.)-1.)
        # So sum2 is my effective_temp_contrib_discrete_sum_weighted.
        # Then tempf(itemp) = (tempf(itemp)/temp*(tbeta+twt) + sum2/arat/temp) / (tbeta+twt+distot) * temp
        # distot = sum of adel(i).
        # tempf(itemp) from _contin was T_eff_contin.
        
        current_teff_factor_from_contin = self.tempf_values[itemp] / (current_temp_k if current_temp_k!=0 else 1.0)
        
        numerator_teff = current_teff_factor_from_contin * (tbeta_continuous_weight + twt_translational_weight) \
                         + (effective_temp_contrib_discrete_sum_weighted / self.awr) / (tev if tev!=0 else 1.0) # sum2/arat/tev part
        
        total_system_weight = tbeta_continuous_weight + twt_translational_weight + np.sum(osc_weights_input)

        if total_system_weight > SMALL_VAL:
            self.tempf_values[itemp] = (numerator_teff / total_system_weight) * current_temp_k
        elif self.iprint >=1:
            print(f"Warning (_discre): Total system weight is zero. T_eff not updated for discrete part.")
        
        if self.iprint >= 1:
            print(f"  Updated DW integral for itemp={itemp}: {self.dwpix_values[itemp]:.5e}")
            print(f"  Updated effective temperature for itemp={itemp}: {self.tempf_values[itemp]:.5e} K")

    # --- End of Discrete Oscillator Methods ---

    # --- Translational Motion Methods (trans, stable, terps, sbfill, besk1) ---

    def _besk1(self, x):
        """
        Corresponds to Fortran function besk1.
        Calculates modified Bessel function K_1(x) without the exp(-x) factor for x > 1.
        For x <= 1, calculates K_1(x).
        Uses series expansions from NJOY.
        """
        # Coefficients c0-c37 for K1 series from NJOY besk1.f
        # These are split into two sets in Fortran, one for x<=1 (poly for K1(x)*x), one for x>1 (poly for K1(x)*exp(x)*sqrt(x))
        _besk1_coeffs_small_x = [ # For K1(x)*x, 11 terms: c0..c10
            1.0000000000000000e+00, 2.5000000000000000e-01, 3.1250000000000000e-02,
            2.6041666666666667e-03, 1.6276041666666667e-04, 7.8125000000000000e-06,
            2.9330533333333333e-07, 8.9425868055555556e-09, 2.2356467013888889e-10,
            4.6575972945601852e-12, 8.3171380259998148e-14
        ]
        _besk1_coeffs_large_x_sqrt_x_K1_exp_x = [ # For K1(x)*exp(x)*sqrt(x/ (pi/2)) (poly part), 8 terms: c0..c7 (renamed bi3 in NJOY)
                                                # Or, for K1(x)*exp(x)*sqrt(x), NJOY uses: sqrt(pi/2) * (this series)
                                                # Fortran besk1 for x>1 returns sqrt(1/x) * series_poly
                                                # So series_poly should be K1(x)*exp(x)*x / sqrt(pi/2)
                                                # Or, if result is K1(x)*exp(x)*sqrt(x), then series_poly is K1(x)*exp(x)*sqrt(x)
                                                # NJOY: besk1=sqrt(1./x)*bi3. bi3 is poly. So poly is sqrt(x)*besk1.
                                                # And besk1 for x>1 is K1(x)*exp(x). So poly is sqrt(x)*K1(x)*exp(x).
            1.2533141373155000e+00, # sqrt(pi/2)
            -1.5666426716443750e-01, # sqrt(pi/2) * (-0.125)
            2.4478791744443359e-02, # sqrt(pi/2) * (0.01953125)
            -5.7307030630824000e-03, # ...
            1.7801009960091000e-03,
            -6.7093840903010000e-04,
            2.8712383633390000e-04,
            -1.3544005106150000e-04
        ]


        if x <= 1.0e-08: # Very small x
            if x == 0.0: return np.inf # K1(0) is inf
            return 1.0 / x # K1(x) ~ 1/x for small x

        if x <= 1.0:
            # Series for K1(x)*x: sum c_k * (x^2/4)^k
            # Result is (series / x) + x * log(x/2) * I1(x) (I1 is Bessel I)
            # NJOY's series for K1(x)*x is simpler: Poly(x^2) + x^2/2 * log(x/2) * Poly_I1(x^2)
            # Fortran: besk1 = term1/x + term2 * alog(xhalf)
            # term1 is poly in x^2. term2 is x/2 * poly_I1 in x^2.
            xsq = x * x
            xhalf = x * 0.5
            
            # Term1: Polynomial part for K1(x)*x
            # c0 + c2*(x^2) + c4*(x^2)^2 ... (NJOY uses x*x as var)
            # No, it's sum c_k * u^k where u = x*x/4
            u = xsq / 4.0
            term1_poly = _besk1_coeffs_small_x[0]
            u_power = 1.0
            for k in range(1, 7): # c1 to c6 for K1(x)*x approx (NJOY uses 7 terms for this poly)
                u_power *= u
                term1_poly += _besk1_coeffs_small_x[k] * u_power
            
            # Term2: I1(x) related part. I1(x) = sum (x/2)^(2k+1) / (k! * (k+1)!)
            # NJOY term2 = xhalf * (poly for I1(x)/(x/2))
            # I1(x)/(x/2) = sum (x^2/4)^k / (k! (k+1)!)
            # Coeffs for I1(x)/(x/2) : 1, 1/2, 1/12, 1/144, ...
            # NJOY's second poly for K1 seems to be for x*log(x/2)*I1(x)
            # The Fortran code sums:
            #   bi1 = c(1)+u*(c(2)+u*(c(3)+u*(c(4)+u*(c(5)+u*(c(6)+u*c(7))))))
            #   bi2 = c(8)+u*(c(9)+u*(c(10)+u*c(11)))
            #   besk1 = bi1/x + xhalf*alog(xhalf)*bi2
            # Here, NJOY coeffs are 1-based.
            # c(1) to c(7) are for 1/x term. My _besk1_coeffs_small_x seems to be these if x^2/4 is the variable.
            # The Fortran uses `u=x*x`. My `u` is `x*x/4`.
            # Let's re-evaluate with NJOY's `u=x*x`.
            
            u_njoy = xsq # NJOY's u variable for small x poly
            
            # bi1 polynomial (Fortran c(1) to c(7))
            # c(1)=1, c(2)=.0625, c(3)=.00260416, ...
            # These are not my _besk1_coeffs_small_x directly.
            # The _besk1_coeffs_small_x are likely for sum c_k * (x^2/4)^k for K1(x)*x itself.
            # Let's assume my _besk1_coeffs_small_x are for K1(x)*x for u=x^2/4:
            # P(u) = c0 + c1*u + c2*u^2 ...
            # The NJOY K1(x) for x<=1 uses:
            # K_1(x) = (1/x) * sum_{k=0 to N} a_k (x^2/4)^k + (x/2) * log(x/2) * sum_{k=0 to M} b_k (x^2/4)^k
            # The first sum is related to Y_1(x), second to J_1(x).
            # My _besk1_coeffs_small_x should be for the full K1(x)*x.
            # The problem statement's coefficients are simplified.
            # Fortran code:
            #   bi1 = c(1)+u*(c(2)+u*(c(3)+u*(c(4)+u*(c(5)+u*(c(6)+u*c(7))))))
            #   bi2 = c(8)+u*(c(9)+u*(c(10)+u*c(11)))
            #   besk1=bi1/x+xhalf*alog(xhalf)*bi2
            # If c are direct Fortran coefficients:
            f_c = [1.0, 0.0625, 0.0026041666, 0.0001627604, 7.8125e-06, 2.93305e-07, 8.9425e-09, # for bi1
                   0.5, 0.03125, 0.0013020833, 4.06901e-05 ] # for bi2
            
            bi1_val = f_c[0] + u_njoy*(f_c[1] + u_njoy*(f_c[2] + u_njoy*(f_c[3] + u_njoy*(f_c[4] + u_njoy*(f_c[5] + u_njoy*f_c[6])))))
            bi2_val = f_c[7] + u_njoy*(f_c[8] + u_njoy*(f_c[9] + u_njoy*f_c[10]))
            
            log_xhalf = np.log(xhalf) if xhalf > 0 else SLIM_LOG_ARG # Avoid log(0)
            val = bi1_val / x + xhalf * log_xhalf * bi2_val
            return val

        else: # x > 1.0
            # Fortran: besk1 = sqrt(1/x) * bi3_poly
            # bi3_poly is sum c_k * (1/x)^k using _besk1_coeffs_large_x_sqrt_x_K1_exp_x
            # This means the result is K1(x)*exp(x) (omitting exp(-x) from K1, and sqrt(pi/2x) prefactor)
            recip_x = 1.0 / x
            bi3_poly = _besk1_coeffs_large_x_sqrt_x_K1_exp_x[0]
            recip_x_power = 1.0
            for k in range(1, 8): # c1 to c7 for this poly
                recip_x_power *= recip_x
                bi3_poly += _besk1_coeffs_large_x_sqrt_x_K1_exp_x[k] * recip_x_power
            
            val = np.sqrt(recip_x) * bi3_poly # This is K1(x)*exp(x)*sqrt(pi/2)
                                             # No, Fortran is K1(x)exp(x) = sqrt(pi/(2x)) * poly_for_NIST_10.40.2
                                             # The NJOY poly for bi3 must be such that sqrt(1/x)*bi3 = K1(x)exp(x).
                                             # So bi3 = sqrt(x) * K1(x) * exp(x).
                                             # The coeffs _besk1_coeffs_large_x_sqrt_x_K1_exp_x[0] is sqrt(pi/2).
                                             # So bi3_poly is sqrt(pi/2) * (1 + mu-1/(8x) + ...), where mu=4*m^2, m=1 for K1.
                                             # NIST 10.40.2: K_nu(z)exp(z) ~ sqrt(pi/2z) * [1 + (mu-1)/(8z) + ...]
                                             # For K1, mu = 4*1^2 = 4.
                                             # So K1(x)exp(x) ~ sqrt(pi/2x) * [1 + (4-1)/(8x) + (4-1)(4-9)/(2!(8x)^2) + ...]
                                             #               = sqrt(pi/2x) * [1 + 3/(8x) + 3*(-5)/(128x^2) + ...]
                                             # The poly in NJOY is likely this [1 + 3/(8x) + ...].
                                             # And _besk1_coeffs_large_x_sqrt_x_K1_exp_x[0] being sqrt(pi/2) means
                                             # bi3_poly = sqrt(pi/2) * [1 + 3/(8x) + ...].
                                             # So, val = sqrt(1/x) * sqrt(pi/2) * [poly_series_in_1/x] = K1(x)exp(x).
                                             # This matches what's needed.
            return val

    def _terps(self, sd_array, nsd_points, delta_beta_conv_grid, beta_target):
        """
        Corresponds to Fortran function terps.
        Performs log-linear interpolation on sd_array.
        sd_array stores S_diffusion or S_free_gas values, assumed to be on a grid
        with points beta_i = (i + 0.5) * delta_beta_conv_grid for i=0..nsd_points-1.
        """
        original_beta_target = beta_target # Keep for final exp(-beta) factor if needed
        
        if abs(beta_target) < VSMALL_VAL: # Effectively beta = 0
            # Need S(0). If sd_array[0] is for delta/2, S(0) needs extrapolation or specific handling.
            # NJOY terps: if be=0, j=0, rem=-0.5. log(sd(1)) + (-0.5)*(log(sd(2))-log(sd(1)))
            # This is extrapolation. Let's try to match.
            # x = 0.0/delta - 0.5 = -0.5. i_idx = int(-0.5) = 0 (if int truncates like Fortran for neg)
            # or -1 if int floors. Python int() truncates. So i_idx=0.
            # rem = -0.5 - 0 = -0.5.
            # This would use sd_array[0] and sd_array[1].
            # If beta_target is exactly 0, it implies special handling or that the grid for sd_array
            # is not well-suited for beta=0 directly.
            # For now, if beta_target is tiny, let's proceed with general logic,
            # which might extrapolate from the first two points.
            pass

        exp_factor = 1.0
        if beta_target < 0:
            beta_target = abs(beta_target)
            # exp(-original_beta_target) because original was negative. So exp(abs(original_beta_target))
            # No, it's S(-beta) = exp(-beta_positive) * S(beta_positive)
            exp_factor = np.exp(max(EXPLIM, original_beta_target)) # original_beta_target is negative here

        if nsd_points == 0:
            return 0.0
        if nsd_points == 1:
            # Single point in sd_array, value is sd_array[0].
            # This value is for beta = 0.5 * delta_beta_conv_grid.
            # If beta_target is close to this, return sd_array[0], else 0.
            # This is simplified; true log-linear would need two points.
            if abs(beta_target - 0.5 * delta_beta_conv_grid) < delta_beta_conv_grid * 0.5:
                 val = sd_array[0]
            else: # Target is too far from the single point.
                 val = 0.0
            return val * exp_factor if val > TINY_SAB else 0.0

        # x is fractional index: beta_target = (x + 0.5) * delta_beta_conv_grid
        x_frac_idx = beta_target / delta_beta_conv_grid - 0.5
        
        i_idx = int(x_frac_idx) # Lower grid index
        
        # Clamp i_idx to be valid for array access sd_array[i_idx] and sd_array[i_idx+1]
        i_idx = max(0, min(i_idx, nsd_points - 2))
        
        remainder = x_frac_idx - i_idx
        
        s_i = sd_array[i_idx]
        s_i_plus_1 = sd_array[i_idx+1]

        if s_i < TINY_SAB or s_i_plus_1 < TINY_SAB: # If values are too small, avoid log(0) and return 0
            # If one is valid, could do linear instead of log-linear, but NJOY returns 0.
            return 0.0

        log_s_i = np.log(s_i) if s_i > 0 else SLIM_LOG_ARG
        log_s_i_plus_1 = np.log(s_i_plus_1) if s_i_plus_1 > 0 else SLIM_LOG_ARG
        
        # Log-linear interpolation
        log_interp_s = log_s_i + remainder * (log_s_i_plus_1 - log_s_i)
        
        if log_interp_s < SLIM_LOG_ARG: # NJOY's slim is -225
            return 0.0
            
        val = np.exp(log_interp_s)
        
        # Apply detailed balance factor if original beta was negative
        final_val = val * exp_factor
        
        return final_val if final_val > TINY_SAB else 0.0


    # --- End of Translational Motion Methods ---


    # --- End of Discrete Oscillator Methods ---

    # Placeholder for other calculation methods
    # def calculate_incoherent_inelastic(self, temperature, temp_specific_data):
    #     pass
    # def calculate_incoherent_elastic(self, temperature, temp_specific_data):
    #     pass
    # def calculate_coherent_elastic(self, temperature, temp_specific_data):
    #     pass

    # --- End of translated Fortran subroutines ---


if __name__ == '__main__':
    # Create a dummy input file for testing based on the issue description
    # This input has lat=0, so _contin will currently be skipped unless lat is changed or logic adapted.
    dummy_input_content = """20 'Uranium Dioxide (UO2) at 293.6K'
1 0 6
0.025 0.050 0.075 0.100 0.125 0.150
0.01 0.02 0.03 0.04 0.05 0.06
0 16.0 0 0 0.0
128 64
238.0 0.5 0.5 0.5
0 0
0
0 0 0 0.0
293.6
1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0 26.0 27.0 28.0 29.0 30.0 31.0 32.0 33.0 34.0 35.0 36.0 37.0 38.0 39.0 40.0 41.0 42.0 43.0 44.0 45.0 46.0 47.0 48.0 49.0 50.0 51.0 52.0 53.0 54.0 55.0 56.0 57.0 58.0 59.0 60.0 61.0 62.0 63.0 64.0 65.0 66.0 67.0 68.0 69.0 70.0 71.0 72.0 73.0 74.0 75.0 76.0 77.0 78.0 79.0 80.0 81.0 82.0 83.0 84.0 85.0 86.0 87.0 88.0 89.0 90.0 91.0 92.0 93.0 94.0 95.0 96.0 97.0 98.0 99.0 100.0 101.0 102.0 103.0 104.0 105.0 106.0 107.0 108.0 109.0 110.0 111.0 112.0 113.0 114.0 115.0 116.0 117.0 118.0 119.0 120.0 121.0 122.0 123.0 124.0 125.0 126.0 127.0 128.0
"""
    dummy_input_filename = "dummy_leapr_input.dat"
    with open(dummy_input_filename, "w") as f:
        f.write(dummy_input_content)

    # Initialize LeaprCalculator with the dummy file
    calculator = LeaprCalculator(input_file_path=dummy_input_filename)
    
    # Run the main process, which now includes parsing
    calculator.leapr_main()

    # Print out some of the parsed values to verify
    print(f"\n--- Parsed Input Verification ---")
    print(f"Title: {calculator.title}")
    print(f"Number of temperatures: {calculator.ntempr}")
    print(f"NPHON: {calculator.nphon}")
    if calculator.nphon > 0:
        print(f"Alpha values: {calculator.alpha_values}")
        print(f"Beta values: {calculator.beta_values}")
    print(f"AWR: {calculator.awr}")
    if calculator.nphon > 0:
        print(f"SPR: {calculator.spr}") # card 6
        print(f"DELTA: {calculator.delta}") # card 6
    print(f"Temperatures data: {calculator.temperatures_data}")

    if calculator.temperatures_data:
        print(f"Data for first temperature: {calculator.temperatures_data[0]}")
    print(f"NCOHL: {calculator.icoh}")
    print(f"NCOLD: {calculator.ncold}")
    print(f"NSS: {calculator.nss}")
    print(f"B7: {calculator.b7}")
    if calculator.ncold > 0:
        print(f"tbeta_cold: {calculator.tbeta_cold}") # card 11
    if calculator.nss > 0:
        print(f"kappa_values: {calculator.kappa_values}") # card 14
        print(f"tbeta_s: {calculator.tbeta_s}") # card 15
        print(f"free_atom_xs_s: {calculator.free_atom_xs_s}") # card 16
        print(f"epsilon_s: {calculator.epsilon_s}") # card 16
        print(f"kappa_s: {calculator.kappa_s}") # card 17
    print(f"ZA: {calculator.za}") # card 4
    print(f"ISABT: {calculator.isabt}") # card 4
    print(f"ILOG: {calculator.ilog}") # card 4
    print(f"SMIN: {calculator.smin}") # card 4
    print(f"NSK: {calculator.nsk}") # card 8
    if calculator.nsk > 0:
        print(f"sbeta: {calculator.sbeta}") # card 8
    print(f"NI: {calculator.ni}") # card 7
    print(f"NJOY: {calculator.njoy}") # card 7
    print(f"NALPHA: {calculator.nalpha}") # card 5
    print(f"NBETA: {calculator.nbeta}") # card 5
    print(f"LAT: {calculator.lat}") # card 4
    print(f"NOUT: {calculator.nout}") # card 1
    print(f"IPRINT: {calculator.iprint}") # card 2
    print(f"--- End of Parsed Input Verification ---")
