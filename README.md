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

    The video, EqWave.mp4, displays dispersion curves of the waves with (black) and without (red) the NCTs. A unique effective buoyancy frequency is used to plot the curves on every frame of the animation. A unique musical pitch is played during every frame (except the last, explained later). The sound frequency played is proportional to the effective buoyancy frequency used. The ratio in the effective buoyancy frequency between any adjacent frames (except the last) is 2^(1/12), which is also the ratio in sound frequency between any adjacent keys in most of the modern musical instruments (sound of piano is used for this demonstration). The last frame used zero effective buoyancy frequency, so it comes with no sound.

 Hing Ong, Apr 5, 2020
