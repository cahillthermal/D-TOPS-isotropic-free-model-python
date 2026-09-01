% main code for signal analysis of D-TOPS updated by Jinchi on Nov 12, 2024
% sample: air/ metal coating/ bulk material to measure
% model: isotropic free thermal expansion model
% signal and fitted properties:
% in-phase and out-of-phase --> thermal conductivity and coefficient of thermal expansion of the bulk material
% or
% ratio --> thermal conductivity of the bulk material
%
% updated by Hyun in 2026:
% 1. Corrected the thermoelastic signal calculation to use T_bs instead of T_s
% 2. Updated to include the temperature-dependent contributions of the air mirage and Fresnel reflection effects

close all;
clearvars;

% define filenames for the data to be analyzed
FileNames_data = "DTOPS_AlCaF2_10x_1p0-1p0_100k-100_313p2mV";

% ======================= sample parameters ===============================

% down parameters:
% Al/CaF2

lambda_down=[150 0.2 10];            % cross-plane thermal conductivity (W/m-K) k update on 20250409 by Jenny
eta_down=[1 1 1];                     % anisotropy of thermal conductivity; eta=kx/ky; kx is in-plane; ky is cross-plane
C_down=[2.65 0.1 2.73]*1e6;           % volumetric heat capacity (J/m^3-K); 2.73 for CaF2
h_down=[60 1 1e6]*1e-9;               % thickness (m)
niu = 0.30;                           % Poisson's ratio of the bulk material
alpha_T = 20e-6;                     % coefficient of thermal expansion (CTE) of the bulk material (K^(-1))

% up parameter: air
lambda_up=0.028;
eta_up=1.0;
C_up=1192;
h_up=1e-3;


% ======================== experimental parameters ========================

% lens
obj=10;                   % e.g., 5 for 5x lens; 10 for 10x lens (5x new)
lens_transmittance = 0.85; % transmittance of lens for pump laser, 0.94 for 1x, 0.86 for 2x, 0.92 for both 5x for TOPS 2.0
focal_length = 5/obj*40e-3; % focal length of the objective lens

% laser spot size and beam offset
r_rms=(5/obj)*12.8e-6; % root-mean-square focused pump and probe beam 1/e^2 radius (m); 12.8 for TOPS 2.0
xoffset=(5/obj)*13.5e-6; % Beam offset (m); pump beam moved downwards by setting 13um move of the gimbal mount 13.4 for TOPS 2.0
C_probe = 0.90;                         % calibrated using CaF_2; depends on the ratio of xoffset and w_rms
w_1_d = 0.92e-3;                        % probe beam 1/e^2 radius at the detector (updated on 20250407)

% set laser power
incident_pump=1e-3;        % average power of digital power (square wave) pump before lens (W)
incident_probe=1e-3;      % laser power of cw probe  before lens(W)

% absorbance laser power
% % Index of refraction of metal coating at the wavelength of the pump laser
n_metal=2.9; % for Al at 780 nm (from Mathewson and Myers):
k_metal=8.2;
% n_metal=2.63; % for NbV Nb0.43V0.57 at lamda 780nm
% k_metal=3.59;
sample_reflectance=abs(n_metal-1+(1i)*k_metal)^2/abs(n_metal+1+(1i)*k_metal)^2; % reflectance of the sample surface
sample_absorbance=1-sample_reflectance;  % absorbance of sample surface
% sample_absorbance=0.4;  % 0.4 for NbV
A_pump=incident_pump*lens_transmittance*sample_absorbance*(4.0/pi); % Amplitude of the primary cosine component of the absorbed pump laser
A_dc_pump=incident_pump*lens_transmittance*sample_absorbance; % total DC component of the absorbed pump;
A_dc_probe=incident_probe*lens_transmittance*sample_absorbance; % total DC component of the absorbed probe;

% calculate the estimated steady state heating separately for pump and probe considering the offset
% conductivity of the bulk material to measure (as set above)
T_ss_heat_pump_est=2*pi*ss_heat(lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,A_dc_pump,xoffset);
T_ss_heat_probe_est=2*pi*ss_heat(lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,A_dc_probe,0.0);
T_ss_heat_est = T_ss_heat_pump_est+T_ss_heat_probe_est

% ====================== Air lens parameters ================================
%
%  Only Z_air is included (the Z_r reflection-phase contribution is NOT included)
%
%  air_lens.enable   : true  -> include the air optical-path-length change contribution (Z_air, Eq.31)
%                      false -> compute the thermoelastic contribution only
%  air_lens.coef_air : -dn_air/dT (K^-1), approx. +9e-7 under standard conditions
%                      (coef_air > 0 because dn_air/dT < 0)

air_lens.enable   = true;    % <- toggle on/off here
air_lens.coef_air = 9e-7;    % -dn_air/dT (K^-1)

% ====================== Fresnel reflection parameters ================================

reflection.enable   = false;    % <- toggle on/off here
reflection.wavelength = 670e-9;    % probe beam reflectance
reflection.dphidT = -6e-5;

% ======================== type of fitting ================================

% fitting1 fits lambda_down(3) and alpha_T using in-phase_and_out-of-phase
% fitting2 fits lambda_down(3) using ratio
% set one of them to 1 and the other to 0; set both to zero for no fitting
FDPBD_fitting1=1;
FDPBD_fitting2=0;

flag_save = 0; % 1 if need to save data (experimental and fitting curve)

% ========================= signal processing =============================

% load the arrays for the out-of-phase signal, in-phase signal,
% ratio (-in-phase/out-of-phase), frequency, and the detector SUM voltage
[V_out_data,V_in_data,V_ratio_data,V_SUM_data,f]=GetData_out_in_ratio_f_VSUM(FileNames_data);

% calculate the "leaking" data for correction of the frequency response due
% to the imperfection of pump modulation and detector response
Amplitude_corrected_3 = -7.7176e-09; % 3rd order
Amplitude_corrected_2 = 2.3877e-06; % 2nd order
Amplitude_corrected_1 = -7.0848e-05; % 1st order
Amplitude_corrected_0 = 1 ; % constant

delay_2 = 5.8742e-12; % 2nd order
delay_1 = -1.0728e-05; % 1st order
delay_0 = 3.3137e-03; % const
% 
% Amplitude_corrected_3 = -6.08e-09; % 3rd order
% Amplitude_corrected_2 = 1.18e-06; % 2nd order
% Amplitude_corrected_1 = 1.50e-04; % 1st order
% Amplitude_corrected_0 = 1 ; % constant
% 
% delay_2 = 7.33e-12; % 2nd order
% delay_1 = -1.10e-05; % 1st order
% delay_0 = 5.08e-03;  % const

% Amplitude_corrected_3 = -9.5095e-09; % 3rd order
% Amplitude_corrected_2 = 2.9421e-06; % 2nd order
% Amplitude_corrected_1 = -8.7298e-05; % 1st order
% Amplitude_corrected_0 = 1.2352e+00 ; % constant
% 
% delay_2 = 5.8742e-12; % 2nd order
% delay_1 = -1.0728e-05; % 1st order
% delay_0 = 3.3137e-03; % const
complex_leaking = (Amplitude_corrected_0 + Amplitude_corrected_1 * sqrt(f) + Amplitude_corrected_2 * sqrt(f).^2+ Amplitude_corrected_3 * sqrt(f).^3).* (exp(1i*(delay_0 + delay_1*f + delay_2*f.^2)));
% correct the measured data using the "leaking" data
[Vcorrected_in,Vcorrected_out,Vcorrected_ratio]=datacorrection_complex_leaking(V_out_data,V_in_data,complex_leaking);

% calculate the average detector SUM voltage
V_SUM_avg = mean(V_SUM_data)*4;

% finding out-of-phase peak fc and reduce the frequency range to fc/10 to fc*10
% assuming that the signal is collected starting from high frequency
length_raw_data = length(f);
Vout_and_f = sortrows([abs(Vcorrected_out) f]);
fc = Vout_and_f(length(f),2);
highlim = 0;
lowlim = 0;
for i_fr = 1:length(f)
    f_temp = f(i_fr);
    if f_temp<fc/10
        lowlim = lowlim + 1;
    else
        if f_temp>fc*10
            highlim = highlim + 1;
        end
    end
end
f = f(1+highlim:length_raw_data-lowlim);
Vcorrected_in=Vcorrected_in(1+highlim:length_raw_data-lowlim);
Vcorrected_out=Vcorrected_out(1+highlim:length_raw_data-lowlim);
Vcorrected_ratio=Vcorrected_ratio(1+highlim:length_raw_data-lowlim);

% convert signal (in unit of voltage) to deflection angle (in unit of rad)
det_factor=(8/pi)^0.5*(focal_length/w_1_d);
% det_factor is the factor such that:
% V_signal*sqrt(2)/V_SUM = det_factor*theta
theta_exp_in=Vcorrected_in*sqrt(2)/V_SUM_avg/det_factor;
theta_exp_out=Vcorrected_out*sqrt(2)/V_SUM_avg/det_factor;
theta_exp_ratio=-theta_exp_in./theta_exp_out;

% ============================= fitting ===================================

if FDPBD_fitting1==1

    % this fitting approach written by Jingyi Zhou in September 2022

    f2=[f,f];
    Xguess=[lambda_down(3),alpha_T];
    lb=[0,-100];   % set lower bounds to fitting parameters
    ub=[100,100];  % set upper bounds to fitting parameters

    max_theta_exp_in = max(abs(theta_exp_in));
    theta_exp_all=[theta_exp_out,theta_exp_in];
    theta_exp_all_norm=[theta_exp_out,theta_exp_in]./[max_theta_exp_in,3.0*max_theta_exp_in];
    % dgc note: changed weightings for in-phase and out-of-phase simultaneous
    % fitting to weight the out-of-phase signal a fixed factor of 3 more than the in-phase

    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    [Xsol,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@(X,f2) FIT_inout(X,f2,max_theta_exp_in,...
        niu,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,C_probe,A_pump,xoffset, air_lens,reflection),...
        Xguess,f2,theta_exp_all_norm,lb,ub,options);
    confidenceinterval95=nlparci(Xsol,residual,'jacobian',J)

    % reset the experimental parameters to the fitted values
    lambda_down(3)=Xsol(1);
    alpha_T = Xsol(2);
    % print the fitted values in command window
    lambda_measure = lambda_down(3)
    alpha_T_measure = alpha_T

end

if FDPBD_fitting2==1

    Xguess=lambda_down(3);
    Xsol=lsqcurvefit(@(X,f) FIT_ratio(X,f,...
        niu,alpha_T,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,C_probe,A_pump,xoffset,air_lens, reflection),...
        Xguess,f,Vcorrected_ratio);

    lambda_down(3)=Xsol;   % reset the experimental parameters to the fitted values
    lambda_measure = lambda_down(3)

end

% ============================= fitting results ===========================

% calculate the model prediction of theta using fitted properties
theta_model=theta_iso_free_thermal_expansion_model(niu,alpha_T,f,lambda_down,C_down,h_down,eta_down,...
    lambda_up,C_up,h_up,eta_up,r_rms,C_probe,A_pump,xoffset, air_lens,reflection);
theta_model;
theta_model_out=imag(theta_model);
theta_model_in=real(theta_model);
theta_model_ratio=-theta_model_in./theta_model_out;

% calculate the steady state heating using the fitted thermal
% conductivity of the bulk material to measure
T_ss_heat_pump=2*pi*ss_heat(lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,A_dc_pump,xoffset);
T_ss_heat_probe=2*pi*ss_heat(lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,A_dc_probe,0.0);
T_ss_heat = T_ss_heat_pump+T_ss_heat_probe

% plot the fitting of theta
figure(1)
subplot(1,2,1)
semilogx(f,1e6*theta_exp_in,'ko','linewidth',1.5); hold on
semilogx(f,1e6*theta_exp_out,'kx','linewidth',1.5); hold on
semilogx(f,1e6*theta_model_in,'r-','linewidth',1.5); hold on
semilogx(f,1e6*theta_model_out,'r--','linewidth',1.5); hold on
box on; axis tight;
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Arial');
xlabel('f (Hz)');
ylabel('in, out-of-phase (\murad)')
subplot(1,2,2)
loglog(f, theta_exp_ratio,'ko','linewidth',1.5); hold on
loglog(f, theta_model_ratio,'r-','linewidth',1.5); hold on
box on; axis tight;
set(gca,'linewidth',1.5,'fontsize',16,'fontname','Arial');
xlabel('f (Hz)');
ylabel('ratio')

% write data and model to a file test.dat
if flag_save == 1
    fileID = fopen('test.dat','w');
    fprintf(fileID,'%12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n', [f, V_SUM_data, theta_exp_in, theta_exp_out, ...
        theta_exp_ratio, theta_model_in, theta_model_out, theta_model_ratio]');
    fclose(fileID);
end