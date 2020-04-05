 Dispersion curves of anelastic equatorial waves
 
 
   EqWave() plots dispersion curves of anelastic equatorial waves with 
   (black) and without (red) the nontraditional Coriolis terms (NCTs) on a
   wavenumber-frequency space. By default, the convective coupling
   parameter is 0 (effective static stability is neutral), and the zonal
   wavenumber range is [-20,20].


   EqWave(alpha) further sets the convective coupling parameter, alpha.
   Effective buoyancy frequency = sqrt(alpha*N2), where N2 = 1.3e-4 s^-2.
   
   alpha = 1; % Dry, Extremely Stable
   
   alpha = 1.6361e-1; % Strongly Stable
   
   alpha = 1.6361e-2; % Mildly Stable
   
   alpha = 1.6361e-4; % Weakly Stable
   
   alpha = 1.6361e-6; % Slightly Stable
   
   alpha = 0; % Neutral
   
   alpha =-1.6361e-6; % Slightly Unstable


   EqWave(alpha,k_range) further sets the zonal wavenumber range, 
   [-k_range,k_range].


 Hing Ong, Mar 31, 2020
