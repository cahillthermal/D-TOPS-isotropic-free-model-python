function [theta_model_all_norm,fit_result]=FIT_inout(X,f2,max_theta_exp_in, ...
    niu,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,C_probe,A_pump,xoffset, air_lens, reflection)

f=f2(:,1);
lambda_down(3)=X(1);
alpha_T=X(2);

theta_model=theta_iso_free_thermal_expansion_model(niu,alpha_T,f,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_rms,C_probe,A_pump,xoffset,air_lens,reflection);
theta_model_out=imag(theta_model);
theta_model_in=real(theta_model);
theta_model_all=[theta_model_out,theta_model_in];
theta_model_all_norm=theta_model_all./[max_theta_exp_in,3.0*max_theta_exp_in];

fit_result=X
