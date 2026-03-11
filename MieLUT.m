clc;clear;close all;
%addpath ('*\ATLAS06March_2026\functions_and_input-spectra\input-spectra')
% Define grids
%%%Choose a range matching your model’s application (typically, 
%%% grain size from ~0.1 µm to 30 µm, with fine resolution, 
%%% e.g., 100–250 bins, logarithmic spacing is recommended).
gsizes = logspace(log10(0.1e-6), log10(33.6e-6), 150); % [m], 0.1 µm to 33.6 µm
% Combine all required wavelengths for your sensors
% wl_hyperion1 = load('S:\Oct28_2025_Mie\Oct28_2025_Mie\Mie_LookUpTable\EMIT_Wavelength.mat'); % e.g., load('hyperion_bands.txt');
% wl_hyperion = wl_hyperion1.EMIT_Wavelength;
% 
% wl_emit1 = load('S:\Oct28_2025_Mie\Oct28_2025_Mie\Mie_LookUpTable\Hyperion_Wavelength.mat');
% wl_emit = wl_emit1.Hyperion_Wavelength;
% 
% wl_prisma1 = load('S:\Oct28_2025_Mie\Oct28_2025_Mie\Mie_LookUpTable\PRISMA_Wavelength.mat');
% wl_prisma = wl_prisma1.PRISMA_Wavelength;
% 
% all_wls = unique(round([wl_hyperion; wl_emit; wl_prisma])); % nm, integer grid
clear wl_prisma1 wl_prisma wl_emit1 wl_emit wl_hyperion1 wl_hyperion
% wavelengths = all_wls(400:900)'; % in nanometers
wavelengths = (300:1:2500)'; % in nanometers
num_gsizes = length(gsizes);
num_wls = length(wavelengths);

% Refractive index for silicates or SPM
m = 1.55 + 0.001i;

Qext = zeros(num_gsizes, num_wls);
Qsca = zeros(num_gsizes, num_wls);
Qabs = zeros(num_gsizes, num_wls);
Qback = zeros(num_gsizes, num_wls);

for ig = 1:num_gsizes
    g = gsizes(ig);
    for il = 1:num_wls
        lambda_nm = wavelengths(il);
        lambda_m = lambda_nm * 1e-9;
        x = 2 * pi * g / lambda_m;
        if x > 0
            abcd = Mie_abcd(m, x);
            an = abcd(1,:);
            bn = abcd(2,:);
            n = 1:length(an);
            cn = 2*n+1;
            Qext(ig, il) = (2 / x^2) * sum(cn .* real(an + bn));
            Qsca(ig, il) = (2 / x^2) * sum(cn .* (abs(an).^2 + abs(bn).^2));
            Qabs(ig, il) = Qext(ig, il) - Qsca(ig, il);
            s = sum(((-1).^n) .* cn .* (an - bn));
            Qback(ig, il) = (1 / x^2) * abs(s).^2;
        end
    end
end

% Save LUT to file
save('Mie_LUT_g_wl.mat', 'gsizes', 'wavelengths', 'Qext', 'Qsca', 'Qabs', 'Qback');
