function [theta_model_ratio,fit_result]=FIT_ratio(X,f, ...
    niu,alpha_T,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,C_probe,A_pump,xoffset,air_lens, reflection)

lambda_down(3)=X;

theta_model=theta_iso_free_thermal_expansion_model(niu,alpha_T,f,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,C_probe,A_pump,xoffset,air_lens, reflection);
theta_model_out=imag(theta_model);
theta_model_in=real(theta_model);
theta_model_ratio=-theta_model_in./theta_model_out;

fit_result=X
