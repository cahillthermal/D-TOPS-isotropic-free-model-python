function [Vcorrected_in,Vcorrected_out,Vcorrected_ratio]=datacorrection_complex_leaking(Vout_data,Vin_data,complex_leaking)

Vcomplex_data=Vin_data+Vout_data.*(1i);
Vcorrected_complex=Vcomplex_data./complex_leaking;
Vcorrected_in=real(Vcorrected_complex);
Vcorrected_out=imag(Vcorrected_complex);
Vcorrected_ratio=-real(Vcorrected_complex)./imag(Vcorrected_complex);
