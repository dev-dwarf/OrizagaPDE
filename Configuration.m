%% General Settings
N = 25; %% Number of grid points

explicit = true;

HEAT = 0;
HEATLL = 1;
ALLENCAHN = 2;
equation = HEATLL; % 0 = heat, 1 = heat w/ lateral loss, 2 = allen-cahn
subdivide = 100; % Subdivisions for reference solutions

save_outputs = false; % save outputs for data analysis
do_plots = ~save_outputs; % whether to generate figures
dtMagnitudes = 6; % Magnitudes to generate data for
dxMagnitudes = 6;

%% PDE parameters
loss_coeff = 0.15*10;
allen_coeff = (1/.15).^2;  % Coefficient needed for Allen-Cahn PDE.
alpha = 1.0; % diffusivity.

%% X Domain
a = 0;
b = pi;
dx = (b-a)/N;
dx2 = dx/2;
dxMagnitude = 0;

%% Time Domain
t0 = 0;
tf = 1.0;
dt =(dx2^2)/2;
dtMagnitude = 0;

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
