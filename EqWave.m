function EqWave(alpha,k_range)
%EqWave dispersion curves of anelastic equatorial waves
%   EqWave() plots dispersion curves of anelastic equatorial waves with 
%   (black) and without (red) the nontraditional Coriolis terms (NCTs) on a
%   wavenumber-frequency space. By default, the convective coupling
%   parameter is 0 (effective static stability is neutral), and the zonal
%   wavenumber range is [-20,20].
%
%   EqWave(alpha) further sets the convective coupling parameter, alpha.
%   Effective buoyancy frequency = sqrt(alpha*N2), where N2 = 1.3e-4 s^-2.
%   alpha = 1; % Dry, Extremely Stable
%   alpha = 1.6361e-1; % Strongly Stable
%   alpha = 1.6361e-2; % Mildly Stable
%   alpha = 1.6361e-4; % Weakly Stable
%   alpha = 1.6361e-6; % Slightly Stable
%   alpha = 0; % Neutral
%   alpha =-1.6361e-6; % Slightly Unstable
%
%   EqWave(alpha,k_range) further sets the zonal wavenumber range, 
%   [-k_range,k_range].

% Hing Ong, Mar 31, 2020

% Setup parameters
Omega = 7.292e-5; % planetary rotation rate, s^-1
N2 = 1.3e-4; % dry buoyancy frequency, s^-2
a = 6.38e6; % planetary radius, m
rH = 1.1e-4; % inverse density scale height, m^-1
m = 2.5e-4; % vertical wavelength, m^-1
n_omega = 65537; % number of frequency grid points
if (nargin < 1)
    alpha = 0;
end
if (nargin < 2)
    k_range = 20;
end

% Calculate parameters
beta = 2 * Omega / a;
Omega2 = Omega*Omega;
rH2 = rH*rH;
m2 = m*m;
m2_r4H2 = m2 + 0.25 * rH2;
aN2 = alpha * N2;
if (aN2+4*Omega2<=0)
    disp('Error: aN2+4*Omega2<=0')
    return
end

% Setup omega axis
omega = linspace(0,Omega,n_omega);
omega2 = omega.*omega;
aN2_om2 = aN2 - omega2;

% Calculate dispersion relation of zero-v waves w/ NCTs
k_Kelvin_a = (Omega*rH - sqrt(Omega2*rH2 + m2_r4H2.*aN2_om2)) .* omega ./ (-aN2_om2);
k_Kelvin_b = (Omega*rH + sqrt(Omega2*rH2 + m2_r4H2.*aN2_om2)) .* omega ./ (-aN2_om2);
k_Kelvin_a(k_Kelvin_a~=real(k_Kelvin_a)) = NaN;
k_Kelvin_b(k_Kelvin_b~=real(k_Kelvin_b)) = NaN;

% Calculate dispersion relation of nonzero-v waves w/ NCTs
A = - aN2_om2;
B = - (beta.*aN2_om2./omega + 2*Omega*rH*omega);
C1 = m2_r4H2*omega2 - Omega*beta*rH;
C2 = - beta * sqrt(m2*aN2_om2 + 0.25*rH2*(aN2_om2+4*Omega2));

n=0;
k_0_a = (-B + sqrt(B.^2 - 4*A.*(C1+C2*(2*n+1)))) *0.5./A;
k_0_b = (-B - sqrt(B.^2 - 4*A.*(C1+C2*(2*n+1)))) *0.5./A;
k_0_a(k_0_a~=real(k_0_a)) = NaN;
k_0_b(k_0_b~=real(k_0_b)) = NaN;
k_0_a(abs(k_0_a-k_Kelvin_a)*2^31<1) = NaN;
k_0_b(abs(k_0_b-k_Kelvin_a)*2^31<1) = NaN;
k_0_a(abs(k_0_a-k_Kelvin_b)*2^31<1) = NaN;
k_0_b(abs(k_0_b-k_Kelvin_b)*2^31<1) = NaN;

n=1;
k_1_a = (-B + sqrt(B.^2 - 4*A.*(C1+C2*(2*n+1)))) *0.5./A;
k_1_b = (-B - sqrt(B.^2 - 4*A.*(C1+C2*(2*n+1)))) *0.5./A;
k_1_a(k_1_a~=real(k_1_a)) = NaN;
k_1_b(k_1_b~=real(k_1_b)) = NaN;
k_1_a1 = k_1_a;
k_1_a2 = k_1_a;
k_1_a1(aN2_om2<=0) = NaN;
k_1_a2(aN2_om2>=0) = NaN;

k_Kelvin_a((aN2_om2+Omega*rH*omega./k_Kelvin_a).*k_Kelvin_a./((aN2_om2+4*Omega2).*omega)<=0) = NaN;
k_Kelvin_b((aN2_om2+Omega*rH*omega./k_Kelvin_b).*k_Kelvin_b./((aN2_om2+4*Omega2).*omega)<=0) = NaN;

% Calculate dispersion relation of zero-v waves w/o NCTs
nk_Kelvin_a = (- sqrt(m2_r4H2*omega2.*aN2_om2)) ./ (-aN2_om2);
nk_Kelvin_b = (+ sqrt(m2_r4H2*omega2.*aN2_om2)) ./ (-aN2_om2);
nk_Kelvin_a(nk_Kelvin_a~=real(nk_Kelvin_a)) = NaN;
nk_Kelvin_b(nk_Kelvin_b~=real(nk_Kelvin_b)) = NaN;

% Calculate dispersion relation of nonzero-v waves w/o NCTs
A = - aN2_om2;
B = - (beta.*aN2_om2./omega);
C1 = m2_r4H2*omega2;
C2 = - beta * sqrt(m2*aN2_om2 + 0.25*rH2*(aN2_om2));

n=0;
nk_0_a = (-B + sqrt(B.^2 - 4*A.*(C1+C2*(2*n+1)))) *0.5./A;
nk_0_b = (-B - sqrt(B.^2 - 4*A.*(C1+C2*(2*n+1)))) *0.5./A;
nk_0_a(nk_0_a~=real(nk_0_a)) = NaN;
nk_0_b(nk_0_b~=real(nk_0_b)) = NaN;
nk_0_a(abs(nk_0_a-nk_Kelvin_a)*2^31<1) = NaN;
nk_0_b(abs(nk_0_b-nk_Kelvin_a)*2^31<1) = NaN;
nk_0_a(abs(nk_0_a-nk_Kelvin_b)*2^31<1) = NaN;
nk_0_b(abs(nk_0_b-nk_Kelvin_b)*2^31<1) = NaN;

n=1;
nk_1_a = (-B + sqrt(B.^2 - 4*A.*(C1+C2*(2*n+1)))) *0.5./A;
nk_1_b = (-B - sqrt(B.^2 - 4*A.*(C1+C2*(2*n+1)))) *0.5./A;
nk_1_a(nk_1_a~=real(nk_1_a)) = NaN;
nk_1_b(nk_1_b~=real(nk_1_b)) = NaN;

nk_Kelvin_a(nk_Kelvin_a./omega<=0) = NaN;
nk_Kelvin_b(nk_Kelvin_b./omega<=0) = NaN;

% Plot
figure('Color','w')
hold on
plot(nk_Kelvin_a*a,omega/Omega,'r-')
plot(nk_Kelvin_b*a,omega/Omega,'r-')
plot(nk_0_a*a,omega/Omega,'r-.')
plot(nk_0_b*a,omega/Omega,'r-.')
plot(nk_1_a*a,omega/Omega,'r--')
plot(nk_1_b*a,omega/Omega,'r--')
plot(k_Kelvin_a*a,omega/Omega,'k-')
plot(k_Kelvin_b*a,omega/Omega,'k-')
plot(k_0_a*a,omega/Omega,'k-.')
plot(k_0_b*a,omega/Omega,'k-.')
plot(k_1_a1*a,omega/Omega,'k--')
plot(k_1_a2*a,omega/Omega,'k--')
plot(k_1_b*a,omega/Omega,'k--')
k_1_b(k_1_b*a>k_range) = NaN;
[~,ind] = max(k_1_b);
set(gca,'Box','on','XLim',[-k_range,k_range],'YLim',[0,omega(ind)/Omega])
xlabel('Zonal Wavenumber')
ylabel('Frequency (CPD)')
title(['{\alpha}N^2 / 4{\Omega}^2 = ',num2str(0.25*aN2/Omega2,'%0.3e')]);
hold off
end

