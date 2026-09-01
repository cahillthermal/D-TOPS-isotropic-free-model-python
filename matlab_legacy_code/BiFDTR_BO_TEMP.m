
function Integrand=BiFDTR_BO_TEMP(kvectin,freq,lambda_up,C_up,h_up,eta_up,lambda_down,C_down,h_down,eta_down,r_pump,r_probe,A_pump)
ii=sqrt(-1);
kvect=kvectin;
kvect2=kvect.^2;
%for liquid layer
alpha_up=lambda_up./C_up;
omega=2*pi*freq;
q2=(ii*omega./alpha_up);

un=sqrt(4*pi^2*eta_up*kvect2+q2);
gamman=lambda_up*un;

G_up=1./gamman; %The liquid G(k)

%for substrate

Nlayers=length(lambda_down); %# of layers
 
alpha_down=lambda_down./C_down;
q2=(ii*omega./alpha_down(Nlayers));
 
un=sqrt(4*pi^2*eta_down(Nlayers)*kvect2+q2);
gamman=lambda_down(Nlayers)*un;
Bplus=0;
Bminus=1;
if Nlayers~=1
    for n=Nlayers:-1:2
        q2=(ii*omega./alpha_down(n-1));
        unminus=sqrt(eta_down(n-1)*4*pi^2*kvect2+q2);
        gammanminus=lambda_down(n-1)*unminus;
        AA=gammanminus+gamman;
        BB=gammanminus-gamman;
        temp1=AA.*Bplus+BB.*Bminus;
        temp2=BB.*Bplus+AA.*Bminus;
        expterm=exp(unminus*h_down(n-1));
        Bplus=(0.5./(gammanminus.*expterm)).*temp1;
        Bminus=0.5./(gammanminus).*expterm.*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic=logical(h_down(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        Bplus(penetration_logic)=0;
        Bminus(penetration_logic)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        un=unminus;
        gamman=gammanminus;
    end
end
 
G_down=(Bplus+Bminus)./(Bminus-Bplus)./gamman; %The substrate G(k)

G = G_up.*G_down./(G_up+G_down);

S=exp(-pi^2*(r_probe^2)/2*kvect.^2); %coaxial S(k)
P=A_pump*exp(-pi^2*(r_pump^2)/2*kvect.^2);
Kernal= S.*P; %The rest of the integrand
Integrand=G.*Kernal;




