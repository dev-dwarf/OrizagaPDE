%% General Settings
HEAT = 0;
HEATLL = 1;
ALLENCAHN = 2;
if ~exist('generateData', 'var')
  dxMagnitude = 0;
  dtMagnitude = 1;
  equation = ALLENCAHN; % 0 = heat, 1 = heat w/ lateral loss, 2 = allen-cahn
end

N = 40*2^dxMagnitude; %% Number of grid points
explicit = true;

tSubd = 8; % Subdivisions for reference solutions - Allen-Cahn / Periodic
xSubd = 4;
dtMagnitudes = 6; % Magnitudes to generate data for
dxMagnitudes = 6;

save_outputs = false; % save outputs for data analysis
do_plots = ~save_outputs; % whether to generate figures

%% PDE parameters
loss_coeff = 0.15*10;
allen_coeff = (1/.15).^2;  % Coefficient needed for Allen-Cahn PDE.
alpha = 1.0; % diffusivity.

%% X Domain
a = 0;
b = pi;
dx = (b-a)/N; % N already includes 2^dxMagnitude term
dx2 = dx/2;

%% Time Domain
t0 = 0;
tf = 1.0;
dt =(dx2^2)/2;
dt = dt / (2^dtMagnitude);

%% Prelude
clf;
clc;
figures_so_far = 1;

switch (equation)
  case 0
    eqTitle = "Heat";
    eqLabel = "Heat";
  case 1
    eqTitle = "Heat w/ Lateral Loss";
    eqLabel = "HeatLateralLoss";
  case 2
    eqTitle = "Allen-Cahn";
    eqLabel = "AllenCahn";
end
