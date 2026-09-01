function theta = theta_iso_free_thermal_expansion_model( ...
    niu, alpha_T, f, ...
    lambda_down, C_down, h_down, eta_down, ...
    lambda_up, C_up, h_up, eta_up, ...
    r_rms, C_probe, A_pump, xoffset, ...
    air_lens, reflection)
% theta_iso_free_thermal_expansion_model
%
% Computes the probe beam deflection angle theta(f) using the isotropic
% free thermal expansion model. Two contributions are summed:
%
%   (1) Thermoelastic: surface deformation due to CTE (alpha_T)
%   (2) Air lens:      probe deflection due to dn_air/dT in the air above
%                      the surface (optional, controlled by air_lens.enable)
%
% T_s vs T_bs correction (thermoelastic path only):
%   BiFDTR_BO_TEMP returns T_s x S_probe (temperature at the air/metal
%   interface). The thermoelastic formula requires T_bs (temperature at the
%   top of the bulk substrate, z=0+). When interface conductance G_int is
%   finite, T_s != T_bs. The conversion uses:
%
%     conv_factor(k) = cosh(zeta2*L) + (lambda_metal*zeta2/G_int)*sinh(zeta2*L)
%     T_bs x S_probe = (T_s x S_probe) / conv_factor
%
%   The air lens path needs T_s (not T_bs), so BiFDTR_BO_TEMP is used
%   directly without any conv_factor correction.
%
% air_lens (optional; disabled if omitted):
%   .enable    - true/false
%   .coef_air  - -dn_air/dT (K^-1), approx +9e-7 at standard conditions

if nargin < 16
    air_lens.enable = false;
end

r_pump  = r_rms;
r_probe = r_rms;
kmax    = 1 / sqrt(r_pump^2 + r_probe^2) * 2;
q2      = (1i) * 2*pi*f * C_down(3) / lambda_down(3);

% Interface conductance and metal coating thickness for T_s -> T_bs conversion
G_int   = lambda_down(2) / h_down(2);
L_metal = h_down(1);

theta = zeros(length(f), 1);

for n = 1:length(f)

    % Scalar thermal wavenumber prefactor for the metal coating at frequency f(n)
    q2_metal_n = 1i * 2*pi*f(n) * C_down(1) / lambda_down(1);

    % ── (1) Thermoelastic contribution ───────────────────────────────────
    % Integrand: 2(1+nu)*alpha_T / (zeta_bulk + k) * T_bs * S_probe
    % T_bs x S_probe = BiFDTR_BO_TEMP / conv_factor(k)
    theta(n,1) = integral(@(kvect) ...
        thermoelastic_integrand(kvect, f(n), niu, alpha_T, q2(n), ...
            eta_down(3), lambda_down(1), q2_metal_n, G_int, L_metal, ...
            r_pump, r_probe, A_pump, xoffset, C_probe, ...
            lambda_up, C_up, h_up, eta_up, lambda_down, C_down, h_down, eta_down), ...
        0, kmax);

    % ── (2) Air lens contribution ─────────────────────────────────────────
    % Integrand: (coef_air / zeta_air) * T_s * S_probe
    % T_s x S_probe = BiFDTR_BO_TEMP directly (no conv_factor needed)
    if air_lens.enable
        q2_air_n = 1i * 2*pi*f(n) * C_up / lambda_up;
        theta(n,1) = theta(n,1) + integral(@(kvect) ...
            air_lens_integrand(kvect, f(n), q2_air_n, ...
                lambda_up, C_up, h_up, eta_up, ...
                lambda_down, C_down, h_down, eta_down, ...
                r_pump, r_probe, A_pump, xoffset, C_probe, air_lens.coef_air), ...
            0, kmax);
    end
        % ── (3) changes in the phase of the φr with temperature ─────────────────────────────────────────
    % Integrand: (coef_air / zeta_air) * T_s * S_probe
    % T_s x S_probe = BiFDTR_BO_TEMP directly (no conv_factor needed)
    if reflection.enable
            theta(n,1) = theta(n,1) + integral(@(kvect) ...
            reflection_integrand(kvect, f(n), reflection.dphidT, ...
        lambda_up, C_up, h_up, eta_up, ...
        lambda_down, C_down, h_down, eta_down, ...
        r_pump, r_probe, A_pump, xoffset, C_probe, reflection.wavelength), ...
            0, kmax);
    end

end
end

% =========================================================================
function val = thermoelastic_integrand(kvect, freq, niu, alpha_T, q2_bulk, ...
        eta_bulk, lambda_metal, q2_metal, G_int, L_metal, ...
        r_pump, r_probe, A_pump, xoffset, C_probe, ...
        lambda_up, C_up, h_up, eta_up, lambda_down, C_down, h_down, eta_down)

k            = 2*pi*kvect;
zeta_metal   = sqrt(q2_metal + k.^2);
zeta_metal_L = zeta_metal .* L_metal;
gamma_metal  = lambda_metal .* zeta_metal;

% conv and corr from transfer matrix + interface conductance
conv = cosh(zeta_metal_L) + (gamma_metal ./ G_int) .* sinh(zeta_metal_L);
corr = sinh(zeta_metal_L) ./ gamma_metal + cosh(zeta_metal_L) ./ G_int;

% BiFDTR_BO_TEMP returns T_s x S_probe
T_s_Sprobe = BiFDTR_BO_TEMP(kvect, freq, lambda_up, C_up, h_up, eta_up, ...
                             lambda_down, C_down, h_down, eta_down, ...
                             r_pump, r_probe, A_pump);

% pump heat flux q_s and probe sensitivity S_probe in k-space (kvect in cycles/m)
flx     = A_pump .* exp(-pi^2 .* r_pump^2  ./ 2 .* kvect.^2);
S_probe =          exp(-pi^2 .* r_probe^2 ./ 2 .* kvect.^2);

% Full expression (no approximation):
%   T_bs x S_probe = conv x (T_s x S_probe) - corr x q_s x S_probe
T_bs_Sprobe = conv .* T_s_Sprobe - corr .* flx .* S_probe;

val = -C_probe .* 8.*pi^2 .* kvect.^2 .* ...
      (-besselj(1, 2.*pi.*kvect.*xoffset)) .* ...
      (2.*(1+niu).*alpha_T ./ ...
       (sqrt(4.*pi^2.*eta_bulk.*kvect.^2 + q2_bulk) + 2.*pi.*kvect)) .* ...
      T_bs_Sprobe;
end

% =========================================================================
function val = air_lens_integrand(kvect, freq, q2_air, ...
        lambda_up, C_up, h_up, eta_up, ...
        lambda_down, C_down, h_down, eta_down, ...
        r_pump, r_probe, A_pump, xoffset, C_probe, coef_air)
% Integral kernel for air lens probe beam deflection.
%
% The air lens contribution arises from the refractive index change of air
% caused by surface heating: Z_air = (-dn_air/dT) / zeta_air * T_s.
%
% T_s x S_probe is obtained directly from BiFDTR_BO_TEMP (no conv_factor
% correction needed, since T_s is the air/metal interface temperature).

k        = 2*pi*kvect;
zeta_air = sqrt(q2_air + k.^2);

% BiFDTR_BO_TEMP returns T_s x S_probe at the surface (air/metal interface)
% This is exactly what the air lens term requires -- no T_bs correction needed
T_s_Sprobe   = BiFDTR_BO_TEMP(kvect, freq, lambda_up, C_up, h_up, eta_up, ...
                               lambda_down, C_down, h_down, eta_down, ...
                               r_pump, r_probe, A_pump);

% Z_air x S_probe = (coef_air / zeta_air) x T_s x S_probe
Z_air_Sprobe = coef_air ./ zeta_air .* T_s_Sprobe;

val = -C_probe .* 8.*pi^2 .* kvect.^2 .* ...
      (-besselj(1, 2.*pi.*kvect.*xoffset)) .* Z_air_Sprobe;
end

% =========================================================================
function val = reflection_integrand(kvect, freq, dphidT, ...
        lambda_up, C_up, h_up, eta_up, ...
        lambda_down, C_down, h_down, eta_down, ...
        r_pump, r_probe, A_pump, xoffset, C_probe, wavelength)

% fresnel reflection change with temperature
% T_s x S_probe is obtained directly from BiFDTR_BO_TEMP (no conv_factor
% correction needed, since T_s is the air/metal interface temperature).

k        = 2*pi*kvect;
% BiFDTR_BO_TEMP returns T_s x S_probe at the surface (air/metal interface)
% This is exactly what the air lens term requires -- no T_bs correction needed
T_s_Sprobe   = BiFDTR_BO_TEMP(kvect, freq, lambda_up, C_up, h_up, eta_up, ...
                               lambda_down, C_down, h_down, eta_down, ...
                               r_pump, r_probe, A_pump);

% Z_air x S_probe = (coef_air / zeta_air) x T_s x S_probe
Z_r_Sprobe = -wavelength .* dphidT .* T_s_Sprobe./(4*pi);

val = -C_probe .* 8.*pi^2 .* kvect.^2 .* ...
      (-besselj(1, 2.*pi.*kvect.*xoffset)) .* Z_r_Sprobe;
end