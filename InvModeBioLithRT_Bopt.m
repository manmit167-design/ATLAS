%--------------------------------------------------------------------------
%{
Mie theory Update by:
1) Dr. Manmit K Singh- Post Doctoral Fellow, 
Department of Earth and Environmental Geosciences,
University of Dayton, 300 College Park, Dayton, OH 45469
2) Dr. Umesh K Haritasya, 
Director, Sustainability Program and Graduate Certificate
Professor, Earth and Environmental Geosciences
University of Dayton, 300 College Park, Dayton, OH 45469
3) Dr. Jeffrey S. Kargel- Senior Scientist, Planetary Science Institute, 
   Tucson, AZ

Base Code:

1) Mr. Enrico Schiassi - PhD Student, Systems and Industrial Engineering, 
   The University of Arizona, Tucson, AZ
2) Dr. Roberto Furfaro - Professor, Systems and Industrial Engineering, 
   The University of Arizona, Tucson, AZ
3) Dr. Jeffrey S. Kargel- Senior Scientist, Planetary Science Institute, 
   Tucson, AZ

%}

%-------------------------------------------------------------------------- 
% InvModeBioLithRT_Bopt computes the objective function for inversion
% problem with Mie LUT interpolation for SPM absorption and backscattering.
%-------------------------------------------------------------------------- 
function [ss]= InvModeBioLithRT_Bopt(Fit,data)
%
%% Inputs

% Tune ph, CDOM, X
%global lambda S_CDOM S_X g_size type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs

% Tune CDOM, X, g_size 
 global C_ph lambda S_CDOM  type_Rrs_below type_case_water g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs

% Tune CDOM, X
%global C_ph lambda S_CDOM S_X g_size type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs

% Tune g_size, X
% global C_ph C_CDOM lambda S_CDOM S_X g_size type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs

% Tune X
%global C_ph C_CDOM lambda S_CDOM S_X g_size type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                    BioLith Model                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Water component concentrations
%% Water component concentratios

% Tune ph, CDOM, and X
%  C_ph= Fit(1);  % phytoplankton 
%  C_CDOM=Fit(2); % CDOM 
%  C_X=Fit(3);    % SPM

% Tune CDOM, X, g_size
 C_CDOM=Fit(1); % CDOM 
 C_X=Fit(2);    % SPM
 g_size=Fit(3); % grain size 

% Tune CDOM, and X
%  C_CDOM=Fit(1); % CDOM 
%  C_X=Fit(2);    % SPM

% Tune g_size, and X
 % g_size=Fit(1); % grain size 
 % C_X=Fit(2);    % SPM

% Tune X
%C_X=Fit;    % SPM

%% Pure water absorption (1/m)
load abs_W.A
wavelength_aw = abs_W(:,1);
absorpt_W = abs_W(:,2);
a_W = zeros(size(lambda));
for i = 1:length(lambda)
    a_W(i) = Absorp_W(lambda(i), wavelength_aw, absorpt_W);
end

%% Phytoplankton absorption (1/m)
load A0_A1_PhytoPlanc.dat
lam_p = A0_A1_PhytoPlanc(:,1);
a0_p = A0_A1_PhytoPlanc(:,2);
a1_p = A0_A1_PhytoPlanc(:,3);

a0 = zeros(size(lambda));
a1 = zeros(size(lambda));
for i=1:length(lambda)
    [a0(i), a1(i)]=a0_a1_phytoplanc1(lambda(i), lam_p, a0_p, a1_p);
end

aph_440 = 0.06 * (C_ph)^0.65;
abs_ph = zeros(size(lambda));
for i=1:length(lambda)
    abs_ph(i) = (a0(i) + a1(i)*log(aph_440)) * aph_440;
end

%% CDOM absorption coefficient [1/m]
Ga_CDOM = 1; Oa_CDOM = 0;
abs_CDOM_440 = (Ga_CDOM * C_CDOM) + Oa_CDOM;
abs_CDOM = abs_CDOM_440 * exp(-S_CDOM * (lambda - 440));

%% SPM absorption coefficient using Mie LUT interpolation
load('M:\ATLAS06March_2026_Final\ATLAS06March_2026_Final\Mie_LUT_g_wl.mat');
%%persistent gsizes wavelengths Qabs Qback
if isempty(gsizes) || isempty(wavelengths) || isempty(Qabs) || isempty(Qback)
    load('Mie_LUT_g_wl.mat', 'gsizes', 'wavelengths', 'Qabs', 'Qback');
end

density = 2600; % kg/m^3
vol_particle = (4/3)*pi*(g_size)^3;
mass_particle = density * vol_particle;
N_particles = (C_X * 1e-3) / mass_particle;

geom_cs = pi * g_size^2;

Qabs_vec = interp2(wavelengths, gsizes, Qabs, lambda, g_size * ones(size(lambda)), 'linear', 0);
abs_X = N_particles * Qabs_vec .* geom_cs;


%% Scattering and backscattering

switch type_case_water
    case 1
        b1 = 0.00144;
    case 2
        b1 = 0.00111;
end
lambda1 = 500;
bb_W = b1 * (lambda / lambda1) .^ (-4.32);

Qback_vec = interp2(wavelengths, gsizes, Qback, lambda, g_size*ones(size(lambda)), 'linear', 0);
bb_x = N_particles * Qback_vec .* geom_cs;

a_tot = a_W + abs_ph + abs_CDOM + abs_X;
bb = bb_W + bb_x;
ext = a_tot + bb;

ext_safe = ext; ext_safe(ext < 1e-10) = 1e-10;
omega_b = bb ./ ext_safe;

view_rad = view * (pi/180);
sun_rad = sun * (pi/180);

switch type_case_water
    case 1
        f_rs = 0.095;
    case 2
        f_rs = 0.0512 * (1 + 4.6659 * omega_b - 7.8387 * omega_b.^2 + 5.4571 * omega_b.^3) ...
            .* (1 + 0.1098 / cos(sun_w)) .* (1 + 0.4021 / cos(view_w));
end

Rrs_below_deep = f_rs .* omega_b;

switch type_Rrs_below
    case 0
        Rrs_below = Rrs_below_deep;
    case 1
        % Reflection factors of bottom surface [1/sr]
        B0=1/pi; B1=1/pi; B2=1/pi; B3=1/pi; B4=1/pi; B5=1/pi; BOTTOM=[B0,B1,B2,B3,B4,B5];
        % Bottom Albedo (constant)
            % wavelenght range [350;900] [nm]
            wavelength= (350:1:900)'; % [nm]
            load Bott0const.R
            Bott0= Bott0const(:,2);
            abott0=zeros(size(lambda)); %

            for i=1:size(lambda,1)
                lam=lambda(i);
                abott0(i)= BottomAlbedo0(lam,wavelength,Bott0);
            end
%             figure(1)
%             plot(lambda,abott0); grid on
%             xlim([400 700])
       
        % Bottom Albedo Sand
         % wavelenght range [350;1000] [nm]
            wavelength= (350:1:1000)'; % [nm]
            load Bott1SAND.R
            Bott1= Bott1SAND(:,2);
            abott1=zeros(size(lambda)); %

            for i=1:size(lambda,1)
                lam=lambda(i);
                abott1(i)= BottomAlbedo1(lam,wavelength,Bott1);
            end
%             figure(2)
%             plot(lambda,abott1); grid on
%             xlim([400 700])
           
        % Bottom Albedo of fine-grained sediment
         % wavelenght range [350;900] [nm]
            wavelength= (350:1:900)'; % [nm]
            load Bott2silt.R
            Bott2= Bott2silt(:,2);
            abott2=zeros(size(lambda));

            for i=1:size(lambda,1)
                lam=lambda(i);
                abott2(i)= BottomAlbedo2(lam,wavelength,Bott2);
            end
%             figure(3)
%             plot(lambda,abott2); grid on
%             xlim([400 700])

        % Bottom Albedo of green makrophyte "Chara contraria"
         % wavelenght range [350;900] [nm]
            wavelength= (350:1:900)'; % [nm]
            load Bott3chara.R
            Bott3= Bott3chara(:,2);
            abott3=zeros(size(lambda));

            for i=1:size(lambda,1)
                lam=lambda(i);
                abott3(i)= BottomAlbedo3(lam,wavelength,Bott3);
            end
%             figure(3)
%             plot(lambda,abott3); grid on
%             xlim([400 700])

        % Bottom Albedo of green makrophyte "Potamogeton perfoliatus"
         % wavelenght range [350;900] [nm]
            wavelength= (350:1:900)'; % [nm]
            load Bott4perfol.R
            Bott4= Bott4perfol(:,2);
            abott4=zeros(size(lambda));

            for i=1:size(lambda,1)
                lam=lambda(i);
                abott4(i)= BottomAlbedo4(lam,wavelength,Bott4);
            end
%             figure(4)
%             plot(lambda,abott4); grid on
%             xlim([400 700])

        % Bottom Albedo of green makrophyte "Potamogeton pectinatus"
         % wavelenght range [350;900] [nm]
            wavelength= (350:1:900)'; % [nm]
            load Bott5pectin.R
            Bott5= Bott5pectin(:,2);
            abott5=zeros(size(lambda));

            for i=1:size(lambda,1)
                lam=lambda(i);
                abott5(i)= BottomAlbedo5(lam,wavelength,Bott5);
            end
%             figure(5)
%             plot(lambda,abott5); grid on
%             xlim([400 700])

        abott=[abott0 abott1 abott2 abott3 abott4 abott5];
       
        Bottom= zeros(size(abott));
        Rrs_Bottom= zeros(size(abott)); % Bottom remote sensing reflectance [1/sr]
       
        for i=1:size(fA,2)
            Bottom(:,i)= fA(i)*abott(:,i);
            Rrs_Bottom(:,i)= BOTTOM(i)* Bottom(:,i); %fA(i)*abott(:,i);
        end
       
        Bottom=sum(Bottom,2);
        Rrs_Bottom=sum(Rrs_Bottom,2); % [1/sr]

        % Attenuation Coefficients
        switch type_case_water
            case 1
                k0=1.0395;
            case 2
                k0=1.0546;
        end
       
        Kd=k0*(ext/cos(sun_w));
        kuW=(ext/cos(view_w)).*((1+omega_b).^3.5421)*(1-(0.2786/cos(sun_w)));
        kuB=(ext/cos(view_w)).*((1+omega_b).^2.2658)*(1-(0.0577/cos(sun_w)));
        %
        Ars1=1.1576; Ars2=1.0389;
       
       
        Rrs_below_shallow= Rrs_below_deep.*(1-(Ars1*exp(-zB*(Kd+kuW)))) + Ars2*Rrs_Bottom.*exp(-zB*(Kd+kuB)) ;
        Rrs_below=Rrs_below_shallow;
end

%% Atmospheric corrections and surface reflectance calculations
load E0.txt;
load absO2.A;
load absO3.A;
load absWV.A;

E0_interp = zeros(size(lambda));
for i=1:length(lambda)
    E0_interp(i) = ExtraSun1(lambda(i), E0(:,1), E0(:,2));
end

abs_O2 = zeros(size(lambda));
abs_O3 = zeros(size(lambda));
abs_WV = zeros(size(lambda));
for i=1:length(lambda)
    abs_O2(i) = ExtraSun1(lambda(i), absO2(:,1), absO2(:,2));
    abs_O3(i) = ExtraSun1(lambda(i), absO3(:,1), absO3(:,2));
    abs_WV(i) = ExtraSun1(lambda(i), absWV(:,1), absWV(:,2));
end

M = 1 / (cos(sun_rad) + 0.50572 * (90 + 6.079975 - sun) ^ (-1.253));
M1 = (M * P) / 1013.25;
Moz = 1.0035 / sqrt(cos(sun_rad) ^ 2 + 0.007);

switch type_case_water
    case 1
        AM = 1;
    case 2
        AM = 1;
end

omega_a = ((-0.0032 * AM) + 0.972) * exp(RH * 3.06 * 10^-4);
Ha = 1; V = 15;
beta = 3.91 * (Ha / V);
tau_a = beta * (lambda / 550) .^ (-alpha);

Tr = exp(-M1 ./ ((115.6406 * lambda .^ 4) - (1.335 * lambda .^ 2)));
Taa = exp(-(1 - omega_a) * tau_a * M);
Tas = exp(-omega_a * tau_a * M);
Toz = exp(-abs_O3 * Hoz * Moz);
To = exp((-1.41 * abs_O2 * M1) ./ (1 + (118.3 * abs_O2 * M1)) .^ 0.45);
Twv = exp((-0.2385 * abs_WV * WV * M) ./ (1 + (20.07 * abs_WV * WV * M)) .^ 0.45);

B3 = 0.82 - 0.1417 * alpha;
B1 = B3 * (1.459 + B3 * (0.1595 + 0.4129 * B3));
B2 = B3 * (0.0783 + B3 * (-0.3824 - 0.5874 * B3));
Fa = 1 - 0.5 * exp((B1 + B2 * cos(sun_rad)) * cos(sun_rad));

Edd = E0_interp .* Tr .* Taa .* Tas .* Toz .* To .* Twv * cos(sun_rad);
Edsr = 0.5 * E0_interp .* (1 - Tr .^ 0.95) .* Taa .* Tas .* Toz .* To .* Twv * cos(sun_rad);
Edsa = E0_interp .* sqrt(Tr) .* Taa .* (1 - Tas) .* Toz .* To .* Twv .* Fa * cos(sun_rad);
Eds = Edsr + Edsa;

Ed = f_dd .* Edd + f_ds .* Eds;
Ls = g_dd .* Edd + g_dsr .* Edsr + g_dsa .* Edsa;

Rrs_above = rho_L .* (Ls ./ Ed);

sigma = 0.03; nW = 1.33; rho_U = 0.54; Q = 5;

den = 1 - rho_U * Q * Rrs_below;
den(den < 1e-10) = 1e-10;

Rrs = ((1 - sigma) .* (1 - rho_L) ./ (nW^2)) .* (Rrs_below ./ den) + Rrs_above;

ss = sum((Rrs_obs - Rrs) .^ 2);

end
