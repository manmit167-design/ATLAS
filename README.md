# ATLAS
Advance Tool for Lake Assessment and Spectral Sensing
MATLAB	R2019b or later
Optimization Toolbox	Required for fmincon (constrained optimization in inverse mode)
MCMC Toolbox (Marko Laine)	Free external toolbox — required for mcmcrun, mcmcplot, chainstats.
Required Folder Structure  (do not rename subfolders)
ATLAS_root/
  ├── main.m          (main script)
  ├── AOP_Rrs_Mie.m
  ├── InvModeBioLithRT_Copt.m
  ├── InvModeBioLithRT_Bopt.m
  └── functions_and_input-spectra/
        ├── functions/        (all helper .m files)
        └── input-spectra/    (spectral library files)
Open the main script and replace the three addpath lines with your local paths:
addpath 'C:\YourPath\ATLAS_root'
addpath 'C:\YourPath\ATLAS_root\functions_and_input-spectra\functions'
addpath 'C:\YourPath\ATLAS_root\functions_and_input-spectra\input-spectra'
Prepare Your Input Data File  (.txt format)
Column	Content	Units / Notes	Mode
Column 1	Wavelength	nm  |  Valid range: 400–2500 nm	Both
Column 2+	Rrs spectrum	1/sr  |  Inverse mode only; script prompts for column number	Inverse
Parameters You Will Be Prompted For at Runtime
The script uses interactive input() prompts. Prepare these values before running:
Modeling mode	—	0 = Inverse  |  1 = Forward	Yes
Viewing angle	degrees	e.g., 30°	Yes
Solar zenith angle	degrees	e.g., 45°	Yes
C_CDOM	mg/m³	0 – 10 	Yes
C_X  (SPM)	g/m³	0 – 300 	Yes
Grain size (d)	µm	1 – 33.6 µm  	Yes
Water type	—	2 = Case-2 (default for glacial lakes)	Preset
Bottom depth (zB)	m	Default 4.00 m	If shallow
Atm. pressure (P)	mbar	Default 1013.25	Yes
Relative humidity	—	Default 0.60	Yes
Ozone height (Hoz)	cm	Default 0.300	Yes
Water vapour (WV)	cm	Default 2.500	Yes
Script paths are hardcoded — must be manually updated for every new machine.
Grain size is discretely binned in Mie LUT — not a continuous inversion.
MCMC runs 4,000 iterations by default — allow 5–15 min per spectrum.
Valid wavelength range: 400–2500 nm. Inputs outside this range will error.
p = 3 parameters hardcoded — change manually if tuning fewer/more variables.
