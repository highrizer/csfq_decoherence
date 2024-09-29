from matplotlib import pyplot as plt
import numpy as np
from scipy import integrate, optimize
from scipy import constants

hbar = constants.hbar
Phi0 = constants.physical_constants["mag. flux quantum"][0]
kb = constants.k
e_charge = constants.e

def n_BE(omega, T):
    """
    Calculates the Bose-Einstein occupation number for a given frequency `omega` and temperature `T`.

    Parameters:
    omega (float): Angular frequency (rad/s).
    T (float): Temperature (K).

    Returns:
    float: Bose-Einstein occupation number.
    """
    return 1 / (np.exp(hbar * omega / kb / T) - 1)

def T_BE(omega, n):
    """
    Calculates the temperature `T` for a given frequency `omega` and occupation number `n` using the Bose-Einstein distribution.

    Parameters:
    omega (float): Angular frequency (rad/s).
    n (float): Bose-Einstein occupation number.

    Returns:
    float: Temperature (K).
    """
    return hbar * omega / kb / np.log(1 / n + 1)

def T_Noise_Q(omega, Ts, atts):
    """
    Calculates the noise temperature across multiple attenuation stages.

    Parameters:
    omega (float): Angular frequency (rad/s).
    Ts (array): Array of temperatures (K) at different stages.
    atts (array): Array of attenuation values (dB).

    Returns:
    array: Noise temperatures at each stage.

    Raises:
    ValueError: If the number of attenuation values is not one less than the number of temperature stages.
    """
    if not len(Ts) == len(atts) + 1:
        raise ValueError("number of attenuations needs to one less than number of temperature stages")
    n_noises = np.zeros_like(Ts)
    T_noises = np.zeros_like(Ts)
    T_noises[0] = Ts[0]
    n_noises[0] = n_BE(omega, Ts[0])
    for i in range(len(atts)):
        n_noises[i + 1] = (
            (1 - 10 ** (-atts[i] / 10)) * n_BE(omega, Ts[i + 1])
            + 10 ** (-atts[i] / 10) * n_noises[i]
        )
        T_noises[i + 1] = T_BE(omega, n_noises[i + 1])
    return T_noises

def S_V_Z(omega, Z, T):
    """
    Calculates the spectral density of voltage noise for a given impedance and temperature.

    Parameters:
    omega (float): Angular frequency (rad/s).
    Z (complex): Impedance (ohms).
    T (float): Temperature (K).

    Returns:
    float: Spectral density of voltage noise.
    """
    if np.abs(omega) < 1e-9:
        return 2 * kb * T * np.real(Z)
    else:
        return hbar * omega * np.real(Z) * (1 + 1 / np.tanh((hbar * omega) / (2 * kb * T)))

def S_flux_bias(omega, Zb, Lb, Mb, T):
    """
    Calculates the flux bias noise spectral density.

    Parameters:
    omega (float): Angular frequency (rad/s).
    Zb (float): Impedance (ohms).
    Lb (float): Inductance (H).
    Mb (float): Mutual inductance (H).
    T (float): Temperature (K).

    Returns:
    float: Flux bias noise spectral density.
    """
    if np.abs(omega) < 1e-9:
        return (Mb**2 / Zb) * (2 * kb * T) / Phi0**2
    else:
        Z = (Zb * 1j * omega * Lb) / (Zb + 1j * omega * Lb)
        return S_V_Z(omega, Z, T) / ((omega * Lb)**2) * (Mb**2) / Phi0**2

def S_simple_1f(omega, Af, alpha):
    """
    Calculates the spectral density of a simple 1/f^α noise model.

    Parameters:
    omega (float): Angular frequency (rad/s).
    Af (float): Amplitude of the noise.
    alpha (float): Exponent of frequency dependence.

    Returns:
    float: Spectral density.
    """
    return Af**2 / ((np.abs(omega) / (2 * np.pi))**alpha)

def S_thermal_1f(omega, A, alpha, T):
    """
    Calculates the noise spectral density for a 1/f^α noise with phenomenological thermal model,
    see Quintana 2017

    Parameters:
    omega (float): Angular frequency (rad/s).
    A (float): Amplitude of the noise.
    alpha (float): Exponent of frequency dependence.
    T (float): Temperature (K).

    Returns:
    float: Thermal noise spectral density.
    """
    return A * omega / (np.abs(omega)**alpha) * (1 + 1 / np.tanh(hbar * omega / (kb * T * 2)))

def S_ohmic(omega, B, gamma, T):
    """
    Calculates the ohmic noise spectral density.

    Parameters:
    omega (float): Angular frequency (rad/s).
    B (float): Amplitude of the noise.
    gamma (float): Exponent of frequency dependence.
    T (float): Temperature (K).

    Returns:
    float: Ohmic noise spectral density.
    """
    if gamma == 1 and np.abs(omega) < 1e-9:
        return B * kb * T * 2 / hbar
    elif np.abs(omega) < 1e-9:
        return B * np.abs(omega)**(gamma - 1) * kb * T * 2 / hbar
    else:
        return B * omega * np.abs(omega)**(gamma - 1) * (1 + 1 / np.tanh(hbar * omega / (kb * T * 2)))

def Gamma_1_qp_eq(coupling, omega, Ej, Delta, T):
    """
    Calculates the quasiparticle-induced decay rate in equilibrium.
    see SI section 3 of L. V. Abdurakhimov et al Appl. Phys. Lett. 115, 262601 (2019)

    Parameters:
    coupling (float): dimensionless coupling matrix element of sin(phi/2) operator
    omega (float): Angular frequency (rad/s).
    Ej (float): Josephson energy in Joule.
    Delta (float): Superconducting gap in Joule.
    T (float): Temperature (K).

    Returns:
    float: Quasiparticle-induced decay rate in equilibrium.
    """
    return coupling**2 * S_qp_eq(omega, Ej, Delta, T)

def Gamma_1_qp_neq(coupling, omega, chi, Ej, Delta):
    """
    Calculates the quasiparticle-induced decay rate out of equilibrium.
    see SI section 3 of L. V. Abdurakhimov et al Appl. Phys. Lett. 115, 262601 (2019)

    Parameters:
    coupling (float): dimensionless coupling matrix element of sin(phi/2) operator
    omega (float): Angular frequency (rad/s).
    chi (float): Non-equilibrium parameter.
    Ej (float): Josephson energy in Joule.
    Delta (float): Superconducting gap in Joule.

    Returns:
    float: Quasiparticle-induced decay rate out of equilibrium.
    """
    return coupling**2 * S_qp_neq(omega, chi, Ej, Delta)

def S_qp_eq(omega, Ej, Delta, T):
    """
    Calculates the quasiparticle spectral density in equilibrium.
    see SI section 3 of L. V. Abdurakhimov et al Appl. Phys. Lett. 115, 262601 (2019)

    Parameters:
    omega (float): Angular frequency (rad/s).
    Ej (float): Josephson energy in Joule.
    Delta (float): Superconducting gap in Joule.
    T (float): Temperature (K).

    Returns:
    float: Quasiparticle spectral density in equilibrium.
    """
    return 16 * Ej / np.pi / hbar * np.exp(-Delta / T / kb) * np.exp(hbar * omega / 2 / T / kb) * \
           np.i0(np.abs(hbar * omega / 2 / kb / T))

def S_qp_neq(omega, chi, Ej, Delta):
    """
    Calculates the quasiparticle spectral density out of equilibrium.
    see SI section 3 of L. V. Abdurakhimov et al Appl. Phys. Lett. 115, 262601 (2019)

    Parameters:
    omega (float): Angular frequency (rad/s).
    chi (float): Non-equilibrium parameter.
    Ej (float): Josephson energy in Joule.
    Delta (float): Superconducting gap in Joule.

    Returns:
    float: Quasiparticle spectral density out of equilibrium.
    """
    return chi * 8 * Ej / np.pi / hbar * np.sqrt(2 * Delta / hbar / omega)

def Gamma_1(transverse_coupling, noise_power):
    """
    Calculates the decay rate for transverse coupling to noise power.

    Parameters:
    transverse_coupling (float): Transverse coupling constant.
    noise_power (float): Noise power spectral density.

    Returns:
    float: Decay rate.
    """
    return (2 * np.pi)**2 * np.abs(transverse_coupling)**2 * noise_power

def Gamma_phi_white(longitudinal_coupling, noise_power):
    """
    Calculates the dephasing rate for longitudinal coupling to white noise power.

    Parameters:
    longitudinal_coupling (float): Longitudinal coupling constant (1/s).
    noise_power (float): Noise power spectral density.

    Returns:
    float: Dephasing rate.
    """
    return (2 * np.pi)**2 * np.abs(longitudinal_coupling)**2 * noise_power

def Gamma_purcell(g, kappa, omega_r, omega_q):
    """
    Calculates the Purcell decay rate.

    Parameters:
    g (float): Coupling constant (rad/s).
    kappa (float): Cavity decay rate (rad/s).
    omega_r (float): Resonator frequency (rad/s).
    omega_q (float): Qubit frequency (rad/s).

    Returns:
    float: Purcell decay rate.
    """
    return (kappa * (g / (omega_r - omega_q)) ** 2)

def A_simple_to_A_thermal(A, alpha, T):
    """
    Converts the amplitude of simple 1/f noise to thermal 1/f noise.

    Parameters:
    A (float): Amplitude of the noise.
    alpha (float): Exponent of frequency dependence.
    T (float): Temperature (K).

    Returns:
    float: Thermal noise amplitude.
    """
    return A / 2 * hbar / kb / T * (2 * np.pi) ** alpha

def eta_0(alpha, omega_low, tau):
    """
    Calculates the prefactor for Ramsey decay in the presence of 1/f^α noise.

    Parameters:
    alpha (float): Exponent of frequency dependence.
    omega_low (float): Low-frequency cutoff.
    tau (float): Time duration.

    Returns:
    float: Prefactor for Ramsey decay.
    """
    result = integrate.quad(
        lambda z: 1 / z ** alpha * (np.sin(z / 2) / z * 2) ** 2,
        omega_low * tau,
        np.inf
    )
    return result[0] * (2 * np.pi) ** (alpha - 1)

def eta_1(alpha):
    """
    Calculates the prefactor for echo decay in the presence of 1/f^α noise.

    Parameters:
    alpha (float): Exponent of frequency dependence.

    Returns:
    float: Prefactor for echo decay.
    """
    result = integrate.quad(
        lambda z: 1 / z ** alpha * (np.sin(z / 4) / z * 4) ** 2 * np.sin(z / 4) ** 2,
        0,
        np.inf
    )
    return result[0] * (2 * np.pi) ** (alpha - 1)

def Afz_to_non_prime(Afz_prime, Afx, c_fz_prime_fx):
    """
    Converts prime flux noise to non-prime flux noise.

    Parameters:
    Afz_prime (float): Prime flux noise amplitude.
    Afx (float): Flux noise amplitude.
    c_fz_prime_fx (float): Cross-correlation term.

    Returns:
    float: Non-prime flux noise amplitude.
    """
    return np.sqrt(Afz_prime ** 2 + c_fz_prime_fx * Afz_prime * Afx
                    + 1 / 4 * Afx ** 2)

def c_fz_fx_to_non_prime(Afz_prime, Afx, c_fz_prime_fx):
    """
    Converts the cross-correlation term from prime to non-prime flux noise.

    Parameters:
    Afz_prime (float): Prime flux noise amplitude.
    Afx (float): Flux noise amplitude.
    c_fz_prime_fx (float): Cross-correlation term.

    Returns:
    float: Non-prime cross-correlation term.
    """
    return (c_fz_prime_fx * Afz_prime * Afx + 1 / 2 * Afx ** 2) \
    / Afz_to_non_prime(Afz_prime, Afx, c_fz_prime_fx) / Afx

def dephasing_tunable_csfq_flux_noise(
        dE_dfz_0011,
        dE_dfx_0011,
        Afz_prime,
        Afx,
        c_fz_prime_fx,
        alpha=1,
        N=0,
        omega_low=None,
        tau=None
        ):
    """
    Calculates the dephasing rate for a tunable CSFQ (C-shunt flux qubit) under flux noise.

    Parameters:
    dE_dfz_0011 (float): Energy derivative with respect to flux `fz`.
    dE_dfx_0011 (float): Energy derivative with respect to flux `fx`.
    Afz_prime (float): Prime flux noise amplitude.
    Afx (float): Flux noise amplitude.
    c_fz_prime_fx (float): Cross-correlation term.
    alpha (float): Exponent of frequency dependence (default = 1).
    N (int): 0 for Ramsey, 1 for echo (default = 0).
    omega_low (float): Low-frequency cutoff for noise.
    tau (float): Time duration.

    Returns:
    float: Dephasing rate.
    """
    Afz = Afz_to_non_prime(Afz_prime, Afx, c_fz_prime_fx)
    c_fz_fx = c_fz_fx_to_non_prime(Afz_prime, Afx, c_fz_prime_fx)
    AE_0011 = np.sqrt(
        (Afz * dE_dfz_0011) ** 2
        + (Afx * dE_dfx_0011) ** 2
        + 2 * c_fz_fx * Afx * Afz * dE_dfz_0011 * dE_dfx_0011
    )
    if N == 0:
        eta = eta_0(alpha, omega_low, tau)
    elif N == 1:
        eta = eta_1(alpha)
    Gamma = eta ** (1 / (1 + alpha)) * (AE_0011) ** (2 / (1 + alpha))
    return 1 / Gamma


###dephasing calculation for arbitrary noise PSDs
def ramsey_func(omega, tau):
    """
    Ramsey filter function.

    Parameters:
    omega (float): Frequency (rad/s)
    tau (float): Time

    Returns:
    float: Value of the filter function for Ramsey decay.
    """
    return np.sinc(omega * tau / (2 * np.pi))**2

def echo_func(omega, tau):
    """
    Echo filter function.

    Parameters:
    omega (float): Frequency (rad/s)
    tau (float): Time

    Returns:
    float: Value of the filter function for Echo decay.
    """
    return np.sinc(omega * tau / (4 * np.pi))**2 * np.sin(omega * tau / 4)**2

def chi(tau, S_omega_01_func, filter_func, cutoff=2 * np.pi/(1000*100e-6)):
    """
    Computes the Chi integral for Ramsey or Echo decay.

    Parameters:
    tau (float): Time
    S_omega_01_func (callable): Noise PSD of omega_01,
    in (rad/s)/Hz
    filter_func (callable): Filter function
    cutoff (float): Lower integration limit
    Returns:
    float: Value of the Chi integral.
    """
    # Define the integrand
    def integrand(omega):
        return S_omega_01_func(omega) * filter_func(omega, tau)

    # Perform the integration from cutoff to infinity
    result, _ = integrate.quad(integrand, cutoff, np.inf)

    # Return the final result
    return (1 / (2 * np.pi)) * tau**2 * result

def t_phi_arb_psd(S_omega_01_func, t_phi_init=1e-7):
    """
    Calculates dephasing times T2Echo and T2Ramsey and their corresponding rates.
    
    Parameters:
    S_omega_01_func (callable): Noise PSD of omega_01,
    in (rad/s)/Hz
    in (rad/s)/Hz
    
    Returns:
    dict: Dictionary with dephasing times and rates for Ramsey and Echo decays.
    """
    # Find root for Ramsey
    ramsey_sol = optimize.root(lambda tau: chi(tau, S_omega_01_func, ramsey_func) - 1, t_phi_init)
    Tramsey = min(ramsey_sol.x) if ramsey_sol.success else np.nan

    # Find root for Echo
    echo_sol = optimize.root(lambda tau: chi(tau, S_omega_01_func, echo_func) - 1, t_phi_init)
    Techo = min(echo_sol.x) if echo_sol.success else np.nan

    return {
        "T2Echo": Techo,
        "T2Ramsey": Tramsey,
        "GammaEcho": 1 / Techo if Techo > 0 else np.nan,
        "GammaRamsey": 1 / Tramsey if Tramsey > 0 else np.nan
    }
    
if __name__ == "__main__":
    omega_low = 2 * np.pi / (1000 * 100e-6)
    tau = 100e-9
    # alphas = np.linspace(0.7, 1.2, 51)
    # eta_0s = [eta_0(alpha, omega_low, tau) for alpha in alphas]
    # eta_1s = [eta_1(alpha) for alpha in alphas]
    # #plt.plot(alphas, eta_0s, label="eta 0")
    # plt.plot(alphas, eta_1s, label="eta 1")
    # plt.legend()
    # plt.show()
    T2 = dephasing_tunable_csfq_flux_noise(
        40 * 2 * np.pi * 1e9,
        40 * 2 * np.pi * 1e9,
        Afz_prime=10e-6,
        Afx=6e-6,
        c_fz_prime_fx=0.5,
        omega_low=omega_low,
        tau=tau
    )
    print (T2)