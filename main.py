import os
import argparse
import numpy as np

from data_io import get_data_out_in_ratio_f_vsum, datacorrection_complex_leaking
from heat_diffusion import ss_heat
from deflection_model import theta_iso_free_thermal_expansion_model
from fitting import fit_inout, fit_ratio


def main():
    parser = argparse.ArgumentParser(description="D-TOPS Signal Analysis (Python)")
    parser.add_argument(
        "file_pos", type=str, nargs="?", default=None,
        help="Path to raw lock-in data text file (positional)"
    )
    parser.add_argument(
        "--file", type=str,
        # default="DTOPS_AlCaF2_10x_1p0-1p0_100k-100_313p2mV.txt",
        default = r"C:\Users\d-cahill\OneDrive - University of Illinois - Urbana\Documents\Data\TOPS data\aug2526\sic_dtops_highf.txt",
        help="Path to raw lock-in data text file"
    )
    parser.add_argument("--save", action="store_true", help="Save experimental and fitted data to test.dat")
    parser.add_argument("--no-plot", action="store_true", help="Disable plotting")
    args = parser.parse_args()

    data_path = args.file_pos if args.file_pos is not None else args.file
    if not os.path.exists(data_path) and not os.path.exists(f"{data_path}.txt"):
        # Fallback search in current directory
        base_name = os.path.basename(data_path)
        if os.path.exists(base_name):
            data_path = base_name
        elif os.path.exists(f"{base_name}.txt"):
            data_path = f"{base_name}.txt"

    print(f"Loading data from: {data_path}")

    # ======================= sample parameters ===============================
    # 
    lambda_down = [20.0, 0.2, 400]      # cross-plane thermal conductivity (W/m-K)
    eta_down = [1.0, 1.0, 1.0]           # anisotropy of thermal conductivity (kx/ky)
    C_down = [2.65e6, 0.1e6, 2.21e6]      # volumetric heat capacity (J/m^3-K); 2.73 for CaF2, 2.21 for SiC, 2.65 for NbV
    h_down = [60e-9, 1e-9, 1e-3]         # thickness (m)
    niu = 0.30                           # Poisson's ratio of the bulk material
    alpha_T = 2e-6                      # coefficient of thermal expansion (CTE) (K^-1)

    # Air
    lambda_up = 0.028
    eta_up = 1.0
    C_up = 1192.0
    h_up = 1e-3

    # ======================== experimental parameters ========================
    obj = 2.0                            # e.g., 10 for 10x lens
    lens_transmittance = 0.85            # 0.85 for use with window; 0.92 without window
    focal_length = 5.0 / obj * 40e-3     # focal length of objective lens

    r_rms = (5.0 / obj) * 12.8e-6        # Focused pump/probe beam 1/e^2 radius (m)
    xoffset = (5.0 / obj) * 13.5e-6      # Beam offset (m)
    C_probe = 0.90                      # Probe factor
    w_1_d = 0.92e-3                     # Probe beam 1/e^2 radius at detector

    incident_pump = 4e-3               # Average power of pump before lens (W)
    incident_probe = 1e-3              # Laser power of probe before lens (W)

    # n_metal = 2.9
    # k_metal = 8.2
    n_metal=2.63  # for NbV Nb0.43V0.57 at lamda 780nm or absorbance = 0.40
    k_metal=3.59
    sample_reflectance = (np.abs(n_metal - 1 + 1j * k_metal)**2) / (np.abs(n_metal + 1 + 1j * k_metal)**2)
    sample_absorbance=1-sample_reflectance
    #sample_absorbance=0.40
    A_pump = incident_pump * lens_transmittance * sample_absorbance * (4.0 / np.pi)
    A_dc_pump = incident_pump * lens_transmittance * sample_absorbance
    A_dc_probe = incident_probe * lens_transmittance * sample_absorbance

    # Steady-state heating estimates
    T_ss_heat_pump_est = 2.0 * np.pi * ss_heat(
        lambda_down, C_down, h_down, eta_down,
        lambda_up, C_up, h_up, eta_up,
        r_rms, A_dc_pump, xoffset
    )
    T_ss_heat_probe_est = 2.0 * np.pi * ss_heat(
        lambda_down, C_down, h_down, eta_down,
        lambda_up, C_up, h_up, eta_up,
        r_rms, A_dc_probe, 0.0
    )
    T_ss_heat_est = T_ss_heat_pump_est + T_ss_heat_probe_est
    print(f"Estimated initial steady-state heating: {T_ss_heat_est:.4f} K")

    # Options
    air_lens = {'enable': True, 'coef_air': 9e-7}
    reflection = {'enable': False, 'wavelength': 670e-9, 'dphidT': -6e-5}

    FDPBD_fitting1 = True
    FDPBD_fitting2 = False
    flag_save = args.save

    # ========================= signal processing =============================
    v_out_data, v_in_data, v_ratio_data, v_sum_data, f = get_data_out_in_ratio_f_vsum(data_path)

    amp_3 = -7.7176e-09
    amp_2 = 2.3877e-06
    amp_1 = -7.0848e-05
    amp_0 = 1.0

    delay_2 = 5.8742e-12
    delay_1 = -1.0728e-05
    delay_0 = 3.3137e-03

    complex_leaking = (amp_0 + amp_1 * np.sqrt(f) + amp_2 * (f) + amp_3 * (f**1.5)) * np.exp(
        1j * (delay_0 + delay_1 * f + delay_2 * (f**2))
    )

    v_corr_in, v_corr_out, v_corr_ratio = datacorrection_complex_leaking(v_out_data, v_in_data, complex_leaking)

    v_sum_avg = np.mean(v_sum_data) * 4.0

    # Frequency range reduction around out-of-phase peak fc
    abs_vout = np.abs(v_corr_out)
    peak_idx = np.argmax(abs_vout)
    fc = f[peak_idx]

    lowlim = np.sum(f < fc / 10.0)
    highlim = np.sum(f > fc * 10.0)

    length_raw_data = len(f)
    f = f[highlim : length_raw_data - lowlim]
    v_corr_in = v_corr_in[highlim : length_raw_data - lowlim]
    v_corr_out = v_corr_out[highlim : length_raw_data - lowlim]
    v_corr_ratio = v_corr_ratio[highlim : length_raw_data - lowlim]
    v_sum_data = v_sum_data[highlim : length_raw_data - lowlim]

    det_factor = np.sqrt(8.0 / np.pi) * (focal_length / w_1_d)
    theta_exp_in = v_corr_in * np.sqrt(2) / v_sum_avg / det_factor
    theta_exp_out = v_corr_out * np.sqrt(2) / v_sum_avg / det_factor
    theta_exp_ratio = -theta_exp_in / theta_exp_out

    # ============================= fitting ===================================
    fitted_k = lambda_down[2]
    uncertainty_k = np.nan
    fitted_alpha = alpha_T
    uncertainty_alpha = np.nan

    if FDPBD_fitting1:
        print("Running FDPBD Fitting 1 (in-phase + out-of-phase)...")
        x_guess = [lambda_down[2], alpha_T]
        bounds = ([0.0, -0.01], [1000.0, 0.01])

        x_sol, res, ci, perr = fit_inout(
            f, theta_exp_in, theta_exp_out,
            x_guess, bounds,
            niu, lambda_down, C_down, h_down, eta_down,
            lambda_up, C_up, h_up, eta_up,
            r_rms, C_probe, A_pump, xoffset,
            air_lens=air_lens, reflection=reflection
        )

        lambda_down[2] = x_sol[0]
        alpha_T = x_sol[1]

        fitted_k = x_sol[0]
        uncertainty_k = perr[0]
        fitted_alpha = x_sol[1]
        uncertainty_alpha = perr[1]

        print(f"Fitted bulk thermal conductivity (lambda_down[2]): {lambda_down[2]:.4f} +/- {uncertainty_k:.4f} W/m-K")
        print(f"Fitted bulk CTE (alpha_T): {alpha_T:.4e} +/- {uncertainty_alpha:.4e} 1/K")
        # if ci is not None:
            # print(f"95% Confidence Interval for lambda_down[2]: {ci[0]}")
            # print(f"95% Confidence Interval for alpha_T: {ci[1]}")

    if FDPBD_fitting2:
        print("Running FDPBD Fitting 2 (ratio only)...")
        x_guess = lambda_down[2]

        x_sol, res, perr_k = fit_ratio(
            f, v_corr_ratio,
            x_guess,
            niu, alpha_T, lambda_down, C_down, h_down, eta_down,
            lambda_up, C_up, h_up, eta_up,
            r_rms, C_probe, A_pump, xoffset,
            air_lens=air_lens, reflection=reflection
        )

        lambda_down[2] = x_sol
        fitted_k = x_sol
        uncertainty_k = perr_k
        print(f"Fitted bulk thermal conductivity (lambda_down[2]): {lambda_down[2]:.4f} +/- {uncertainty_k:.4f} W/m-K")

    # Appending record to results.dat
    results_file = "results.dat"
    file_exists = os.path.exists(results_file)
    filename = os.path.basename(data_path)

    with open(results_file, "a") as f_res:
        if not file_exists:
            f_res.write("# filename\tfitted_thermal_conductivity\tuncertainty_thermal_conductivity\tfitted_CTE\tuncertainty_CTE\n")
        f_res.write(f"{filename}\t{fitted_k:.6e}\t{uncertainty_k:.6e}\t{fitted_alpha:.6e}\t{uncertainty_alpha:.6e}\n")

    print(f"Appended fitting results for {filename} to {results_file}")

    # ============================= fitting results ===========================
    # print("Calculating final fitted model prediction...")
    theta_model = theta_iso_free_thermal_expansion_model(
        niu, alpha_T, f,
        lambda_down, C_down, h_down, eta_down,
        lambda_up, C_up, h_up, eta_up,
        r_rms, C_probe, A_pump, xoffset,
        air_lens=air_lens, reflection=reflection
    )

    theta_model_out = np.imag(theta_model)
    theta_model_in = np.real(theta_model)
    theta_model_ratio = -theta_model_in / theta_model_out

    T_ss_heat_pump = 2.0 * np.pi * ss_heat(
        lambda_down, C_down, h_down, eta_down,
        lambda_up, C_up, h_up, eta_up,
        r_rms, A_dc_pump, xoffset
    )
    T_ss_heat_probe = 2.0 * np.pi * ss_heat(
        lambda_down, C_down, h_down, eta_down,
        lambda_up, C_up, h_up, eta_up,
        r_rms, A_dc_probe, 0.0
    )
    T_ss_heat = T_ss_heat_pump + T_ss_heat_probe
    print(f"Final steady-state heating with fitted properties: {T_ss_heat:.4f} K")

    # Save results to test.dat
    if flag_save:
        output_data = np.column_stack([
            f, v_sum_data, theta_exp_in, theta_exp_out,
            theta_exp_ratio, theta_model_in, theta_model_out, theta_model_ratio
        ])
        np.savetxt("test.dat", output_data, fmt="%12.3e")
        print("Fitted results saved to test.dat")

    # Plot results
    if not args.no_plot:
        try:
            import matplotlib.pyplot as plt
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

            ax1.semilogx(f, 1e6 * theta_exp_in, 'ko', label='Exp In-phase')
            ax1.semilogx(f, 1e6 * theta_exp_out, 'kx', label='Exp Out-of-phase')
            ax1.semilogx(f, 1e6 * theta_model_in, 'r-', label='Model In-phase')
            ax1.semilogx(f, 1e6 * theta_model_out, 'r--', label='Model Out-of-phase')
            ax1.set_xlabel('f (Hz)')
            ax1.set_ylabel('In, out-of-phase (μrad)')
            ax1.legend()
            ax1.grid(True, which="both", ls="--", alpha=0.5)

            ax2.loglog(f, theta_exp_ratio, 'ko', label='Exp Ratio')
            ax2.loglog(f, theta_model_ratio, 'r-', label='Model Ratio')
            ax2.set_xlabel('f (Hz)')
            ax2.set_ylabel('Ratio (-in/out)')
            ax2.legend()
            ax2.grid(True, which="both", ls="--", alpha=0.5)

            plt.tight_layout()
            plt.show()
        except ImportError:
            print("matplotlib is not installed; skipping plot generation.")


if __name__ == "__main__":
    main()
