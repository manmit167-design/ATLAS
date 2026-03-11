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

Base Code:

1) Mr. Enrico Schiassi - PhD Student, Systems and Industrial Engineering, 
   The University of Arizona, Tucson, AZ
2) Dr. Roberto Furfaro - Professor, Systems and Industrial Engineering, 
   The University of Arizona, Tucson, AZ
3) Dr. Jeffrey S. Kargel- Senior Scientist, Planetary Science Institute, 
   Tucson, AZ

%}
%--------------------------------------------------------------------------

clear ; clc; close all;
format long

%% Add necessary paths
%add path to the folder containing the functions and the spectra needed for
%the run
% addpath('/ATLAS Folder Location')
% addpath('/functions_and_input-spectra\functions')
% addpath('/functions_and_input-spectra\input-spectra')
addpath M:\ATLAS06March_2026_Final\ATLAS06March_2026_Final
addpath M:\ATLAS06March_2026_Final\ATLAS06March_2026_Final\functions_and_input-spectra\functions
addpath M:\ATLAS06March_2026_Final\ATLAS06March_2026_Final\functions_and_input-spectra\input-spectra

%% Input
%--------------------------------------------------------------------------

%% 1) Global Variables

% Warnings: 
% - 1: the names must be consistent in the main and the correspondent functions
% - 2: the tuned quantities are not part of the global variables
%

% Tune ph, CDOM, X
% global lambda S_CDOM S_X g_size type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs

% Tune CDOM, X, g_size 
 global C_ph lambda S_CDOM S_X type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs 

% Tune CDOM, X
% global C_ph lambda S_CDOM S_X g_size type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs

% Tune g_size, X
% global C_ph C_CDOM lambda S_CDOM S_X g_size type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs


% Tune X
%global C_ph C_CDOM lambda S_CDOM S_X g_size type_Rrs_below zB type_case_water fA g_dd g_dsr g_dsa f_dd f_ds alpha view view_w sun sun_w rho_L P RH Hoz WV  Rrs_obs


%% 2)  Wavelength

%{

Enter the range of wavelenght of interest.
The range allowed is the visible i.e. [400-2500] [nm]
The range MUST match the range used in the data set (if inverse mode is used)

%}

% Ask the user to choose mode (0 = inverse, 1 = forward)
FW_mode_only = input('Enter 0 for Inverse Modeling or 1 for Forward Modeling: ');

% Ask user to select the *.txt file
[filename, pathname] = uigetfile('*.txt', 'Select the input data file');
if isequal(filename,0)
    error('No file selected. Program terminated.');
else
    fullpath = fullfile(pathname, filename);
    disp(['Selected file: ', fullpath]);
end

% addpath S:\HMA_GLAM_BioLith-RT_5-master\Mie_LookUpTable\Mie_LookUpTable\functions_and_input-spectra\functions
% addpath S:\HMA_GLAM_BioLith-RT_5-master\Mie_LookUpTable\Mie_LookUpTable\functions_and_input-spectra\input-spectra
% Load the selected file
Spec_Data = load(fullpath);

switch FW_mode_only
    case 0 % -------- INVERSE MODE --------
        disp('User Chose Inverse Modeling...');
        data.RrsLam = Spec_Data;      % Assign spectral reflectance data
        lambda = data.RrsLam(:,1);    % Extract wavelength (first column)
        
    case 1 % -------- FORWARD MODE --------
        disp('User Chose Forward Modeling...');
        lambda = Spec_Data(:,1);      % Extract wavelength (first column)

    otherwise
        error('Invalid choice! Enter 0 or 1.');
end

% Display all available wavelengths with serial numbers
fprintf('\nAvailable wavelengths:\n');
for i = 1:length(lambda)
    fprintf('%3d: %8.2f nm\n', i, lambda(i));
end

% Ask user to provide wavelength range by serial number
prompt = '\nEnter start and end serial numbers (e.g., 10 60) or press Enter to use default [400–700 nm]: ';
userInput = input(prompt, 's');

if isempty(userInput)
    % Use default range of 400 nm to 700 nm
    defaultIdx = find(lambda >= 400 & lambda <= 700);
    if isempty(defaultIdx)
        error('No wavelengths found in the range 400–700 nm.');
    end
    selectedIdx = defaultIdx;
    fprintf('Using default wavelength range: 400–700 nm\n');
else
    % Parse user input
    rangeValues = str2num(userInput); %#ok<ST2NM> 
    if length(rangeValues) ~= 2 || any(rangeValues < 1) || any(rangeValues > length(lambda))
        error('Invalid input. Please enter two valid serial numbers within range.');
    end
    selectedIdx = rangeValues(1):rangeValues(2);
    fprintf('Using user-selected range: %.2f–%.2f nm\n', lambda(selectedIdx(1)), lambda(selectedIdx(end)));
end

% Final selected wavelength and spectral data
lambda_selected = lambda(selectedIdx);
if FW_mode_only == 0
    Rrs_selected = data.RrsLam(selectedIdx, :);
    data.RrsLam = Rrs_selected;
else
    Spec_Data = Spec_Data(selectedIdx, :);
end
lambda = lambda_selected;
clear lambda_selected


%% 3) Remote sensig reflectance observed/measured [1/sr]

switch FW_mode_only
    
    case 0 % inverse mode
        
        % Display available columns in data.RrsLam
        disp(['The data.RrsLam matrix has ', num2str(size(data.RrsLam, 2)), ' columns.']);
        disp('Column 1 is usually wavelength. Choose column >= 2 for Rrs.');

        % Ask user for column number (excluding wavelength column)
        x = input('Enter the column number for Rrs_obs (e.g., 2 for the first reflectance column): ');

        % Validate input
        if x < 2 || x > size(data.RrsLam, 2)
            error('Invalid column number! Must be between 2 and number of columns in data.RrsLam.');
        end

        % Extract observed Rrs values
        Rrs_obs = data.RrsLam(:, x);

        % Display confirmation
        disp(['Rrs_obs has been extracted from column ', num2str(x), '.']);

%         Rrs_obs= data.RrsLam(:,2); 
        clear x
        
    case 1 % forward mode
        
        Rrs_obs= zeros(size(lambda));
%         Rrs_obs= EMIT(:,2);
end

%% 4) Case Water Selection

type_case_water=2; % 1=case-1 water; 2=case-2 water
%% 5) Angles

% Enter the view (camera) and the sun angles [degrees] - zenith angles


% Ask user to input viewing and sun angles in degrees
view = input('Enter the viewing angle in degrees: ');  % e.g., 30
sun  = input('Enter the sun angle in degrees: ');      % e.g., 89.13

% Display confirmation
disp(['View angle: ', num2str(view), ' deg']);
disp(['Sun angle: ', num2str(sun), ' deg']);


%% 6) Water component concentrations

C_ph = eps    ;                 % [mg/m^3]  

% ---- Ask user for optical and physical parameters ----

% 1) CDOM concentration
C_CDOM = input('Enter CDOM concentration C_CDOM [mg/m^3]: ');

% 2) SPM concentration
C_X = input('Enter Suspended Particulate Matter C_X [g/m^3]: ');

% 3) Grain size in micrometers
g_size_microm = input('Enter grain size in micrometers [µm]: ');

% 4) Slope for CDOM absorption (automatically calculated)
S_CDOM = 0.0092+0.008*exp(-1*C_CDOM);

% 6) Convert grain size to meters
g_size = g_size_microm * 1e-6;

% Display all selected/derived values
disp('--------------------------------------------');
disp('User-selected and derived parameter values:');
disp(['C_CDOM = ', num2str(C_CDOM), ' [mg/m^3]']);
disp(['S_CDOM (derived) = ', num2str(S_CDOM)]);
disp(['C_X = ', num2str(C_X), ' [g/m^3]']);
disp(['S_X = ', num2str(S_X), ' (slope for a_X(lambda))']);
disp(['Grain size = ', num2str(g_size_microm), ' [µm] = ', num2str(g_size), ' [m]']);
disp('--------------------------------------------');

%% 7) Input for the Remote Sensing Reflectance [1/sr] 

disp('--- Input for Water Column Contribution and Atmospheric Parameters ---');

% Water column contribution
type_Rrs_below = input('Type of Rrs below water (0=deep water; 1=shallow water) [default 0]: ');
if isempty(type_Rrs_below)
    type_Rrs_below = 0;
end

% Bottom depth
zB = input('Bottom depth zB [m] [default 4.00]: ');
if isempty(zB)
    zB = 4.00;
end

% Areal fraction of bottom surface
disp('Enter areal fractions (sum should be 1):');
fA0 = input('fA0 (constant) [default 0]: '); if isempty(fA0), fA0 = 0; end
fA1 = input('fA1 (sand) [default 0]: '); if isempty(fA1), fA1 = 0; end
fA2 = input('fA2 (sediment) [default 1]: '); if isempty(fA2), fA2 = 1; end
fA3 = input('fA3 (Chara contraria) [default 0]: '); if isempty(fA3), fA3 = 0; end
fA4 = input('fA4 (Potamogeton perfoliatus) [default 0]: '); if isempty(fA4), fA4 = 0; end
fA5 = input('fA5 (Potamogeton pectinatus) [default 0]: '); if isempty(fA5), fA5 = 0; end
fA = [fA0, fA1, fA2, fA3, fA4, fA5];

% Atmospheric conditions
g_dd = input('g_dd [Irradiance intensity 1/sr] [default 0.05]: '); if isempty(g_dd), g_dd = 0.05; end
g_dsr = input('g_dsr [default 0]: '); if isempty(g_dsr), g_dsr = 0; end
g_dsa = input('g_dsa [default 0]: '); if isempty(g_dsa), g_dsa = 0; end

% Intensities of light sources
f_dd = input('f_dd [default 1]: '); if isempty(f_dd), f_dd = 1; end
f_ds = input('f_ds [default 1]: '); if isempty(f_ds), f_ds = 1; end

% Angstrom exponent
alpha = input('Alpha (Angstrom exponent) [default 1.317]: '); if isempty(alpha), alpha = 1.317; end

% Atmospheric pressure
P = input('Atmospheric pressure P [mbar] [default 1013.25]: '); if isempty(P), P = 1013.25; end

% Relative Humidity
RH = input('Relative Humidity RH [default 0.60]: '); if isempty(RH), RH = 0.60; end

% Scale height for ozone
Hoz = input('Scale height for ozone Hoz [cm] [default 0.300]: '); if isempty(Hoz), Hoz = 0.300; end

% Scale height of precipitable water
WV = input('Water vapor scale height WV [cm] [default 2.500]: '); if isempty(WV), WV = 2.500; end

% Display final chosen values
disp('--- Final Selected/Default Values ---');
fprintf('type_Rrs_below = %d\n', type_Rrs_below);
fprintf('zB = %.2f m\n', zB);
fprintf('fA = [%.2f %.2f %.2f %.2f %.2f %.2f]\n', fA);
fprintf('g_dd = %.3f, g_dsr = %.3f, g_dsa = %.3f\n', g_dd, g_dsr, g_dsa);
fprintf('f_dd = %.2f, f_ds = %.2f\n', f_dd, f_ds);
fprintf('alpha = %.3f\n', alpha);
fprintf('P = %.2f mbar\n', P);
fprintf('RH = %.2f\n', RH);
fprintf('Hoz = %.3f cm\n', Hoz);
fprintf('WV = %.3f cm\n', WV);

%% 8) Tuned quantities

%{

 All the fit quantities must be compotents of the following vector, not
 part of the global variables

%}

% Tune ph, CDOM, X
% Fit= [C_ph, C_CDOM, C_X];

% Tune CDOM, X, g_size 
 Fit= [C_CDOM, C_X, g_size];
  
% Tune CDOM, X
%Fit= [C_CDOM, C_X];

% Tune g_size, X
% Fit= [g_size, C_X];

% Tune X
%Fit= C_X;

%% 9) Geometry
% compute the view (camera) and sun angles in water [rad]; and the Fresnel
% coeff

[view_w, sun_w, rho_L] = Snell_law(view, sun);
%-------------------------------------------------------------------------------
switch FW_mode_only
    
    case 0 % inverse mode
        % Simulated remote sensing reflectance
        [Rrs0, Res0] = AOP_Rrs_Mie(Fit);
        Rrs0 = real(Rrs0);
%-----------------------------------------------------------------------------------  
%% Water component concentrations retrieval (constrained optimization framework)
% Linear Constraints: the water component cocentrations are non-negative
% (inequality const.)
        % Inequality and bound constraints for fmincon
        Aineq = -eye(length(Fit));
        bineq = zeros(length(Fit),1);
        Aeq = []; 
        beq = [];
        lb = []; % Lower bound
        ub = []; % Upper bound
        

        obj = @InvModeBioLithRT_Copt; % Objective function: Res = sum((Rrs_obs-Rrs).^2)
        % Solving linear constrained optimization problem
        [Fitret, Res] = fmincon(obj, Fit, Aineq, bineq, Aeq, beq, lb, ub);
        % Simulated remote sensing reflectance with retrieved water
        % component concentration via Constrained Optimization
        [Rrs_fit, Res_fit] = AOP_Rrs_Mie(Fitret);
   %% Water Component Concentrations retrieval (Bayesian Optimization framework)
        sample_size = length(lambda); % sample size
        p = 3;                        % number of parameters to be tuned-to be set by the user
        % as defined by the user in section 1 of Global Variable
        
       
        
        type_initial_values = 0;  % 0 = Raw guess for the tuned parameters will be the initial guess for the MCMC
                                  % 1 = parameters tuned via Constrained optimization will be the initial guess for the MCMC
        
        switch type_initial_values
            case 0
                mse = mean(Res0) / (sample_size - p); % estimate for the error variance
            case 1
                Fit = Fitret;
                mse = mean(Res_fit) / (sample_size - p); % estimate for the error variance
        end
        % Define the structure for the paraemters to be tuned:
        % {'name'- mandatory, initial value-mandatory, min value, max
        % value, prior_mu, prior_sig, targetflag, localflag
       % Tune ph, CDOM, X
%         params = {
%             {'C_{ph}',  Fit(1), 0}
%             {'C_{CDOM}',Fit(2), 0}
%             {'C_{X}',   Fit(3), 0}
%             };

        % Tune CDOM, X, g_size
         params = {
            {'C_{CDOM}',  Fit(1), 0,  }
            {'C_{X}',     Fit(2), 0, 300 }
            {'d [m]',     Fit(3), 0, 33.6e-06 }
            };
         
%         % Tune CDOM, X
%         params = {
%             %{'C_{ph}',  Fit(1), 0}
%             {'C_{CDOM}',Fit(1), 0, 1}
%             {'C_{X}',   Fit(2), 0, 200}
%             };

        % Tune g_size, X
        % params = {
        % 
        %     {'d [m]',   Fit(1), 0, 50e-05 }
        %     {'C_{X}',   Fit(2), 0, 4}
        %     };
        % 
        % Tune X
%         params = {
%         
%             {'C_{X}',   Fit, 0, 320}
%             };
        
% Set the functions and the options for the MCMC run.
% FUNCTIONS
%         model.ssfun            -2*log(likelihood) function
%         model.priorfunc        -2*log(prior) prior function
%         model.sigma2           initial error varinace
%         model.N                total number of obervations
%         model.S20              prior for sigma2
%         model.N0               prior accuracy for S20
%         model.nbatch           number of datasets
        model.ssfun = @InvModeBioLithRT_Bopt;
        model.sigma2 = mse;
        model.N = sample_size;
        %
        %    OPTIONS
        %    options.nsimu            number of simulations - mandatory
        %    options.qcov             proposal covariance
        %    options.method           'dram','am','dr', 'ram' or 'mh' - mandatory
        %    options.adaptint         interval for adaptation, if 'dram' or 'am' used
        %                             DEFAULT adaptint = 100
        %    options.drscale          scaling for proposal stages of dr
        %                             DEFAULT 3 stages, drscale = [5 4 3]
        %    options.updatesigma      update error variance. Sigma2 sampled with updatesigma=1
        %                             DEFAULT updatesigma=0
        %    options.verbosity        level of information printed
        %    options.waitbar          use graphical waitbar?
        %    options.burnintime       burn in before adaptation starts

        options.nsimu = 4000;
        options.method = 'dram';    % dram is default
        options.updatesigma = 1;
        
        %   mcmc run
        [res, chain, s2chain] = mcmcrun(model, data, params, options);
        
        %   Results
        chainstats(chain, res); % mean and std of the sampled posteriors
        [Rrs_B, Res_B] = AOP_Rrs_Mie(mean(chain)); % simulated Rrs with the Bayesian optimization tuned mean values
        
        % PLOTS
        % Chain Plots
        % mcmc sampling
        figure(1); clf;
        mcmcplot(chain, [], res, 'chainpanel');

        % Posteriors
        figure(3);
        mcmcplot(chain, [], res, 'denspanel', 2);
        
        % Remote sensing reflectance plot
        figure(5);
        plot(lambda, mean(Rrs0,2), lambda, mean(Rrs_fit,2), lambda, mean(Rrs_B,2), lambda, Rrs_obs, 'LineWidth', 3);
        grid on;
        title('Remote Sensing Reflectance', 'FontSize', 24);
        xlabel('Wavelength [nm]', 'FontSize', 20);
        ylabel('Rrs [1/sr]', 'FontSize', 20);
        legend('Rrs guess', 'Rrs fit C', 'Rrs fit B', 'Rrs measured');
        xlim([lambda(1), lambda(end)]);
        
    case 1 % forward mode
        
        % Simulated remote sensing reflectance
        [Rrs0, Res0] = AOP_Rrs_Mie(Fit);
        Rrs0 = real(Rrs0);
        
        % Uncomment next two lines to create synthetic reflectance data
        % noise = rand(size(lambda));
        % Rrs_obs_syn = Rrs0 + 0.1 * Rrs0 + 0.05 * noise .* Rrs0;
        
        figure(1);
        plot(lambda, Rrs0, '-b');
        grid on;
        title('Remote Sensing Reflectance', 'FontSize', 24);
        xlabel('Wavelength [nm]', 'FontSize', 20);
        ylabel('Rrs [1/sr]', 'FontSize', 20);
        xlim([lambda(1), lambda(end)]);
        hold on;
        plot(Spec_Data(:,1), Spec_Data(:,2), '-*r');
end
