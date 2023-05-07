clf;
clc;
clear all;
figures_so_far = 1;

addpath('./mole_MATLAB/');

%% Settings %%%%
N = 100; % number of grid points, x
M = N;

explicit = true;
equation = 2; % 0 = heat, 1 = heat w/ lateral loss, 2 = allen-cahn

%% Shared problem parameters

%% Heat Equation
loss_coeff = 0.15*10;
allen_coeff = (1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
alpha = 1.0; % thermal diffusivity.

%% Space Domain
a = -pi;
b = pi;
c = -pi;
d = pi;
dx = (b-a)/N;
dy = (d-c)/N;
dx2 = dx/2;
dy2 = dy/2;

%% Space discretization
xmesh = [a (a+dx2):dx:(b-dx2) b]'; 
ymesh = [c (c+dx2):dx:(d-dx2) d]';

[X, Y] = meshgrid(xmesh, ymesh);

%% Initial Condition
icfun = @(x,y) sin(pi*x).*sin(pi*y);
icfun = @(x,y) rand(size(x)).*rand(size(y)) - rand(size(x)).*rand(size(y));
U0 = icfun(X,Y);
e = min(min(U0));
f = max(max(U0));
e = -1;
f = 1;

%% Time Domain
t0 = 0;
tf = 10.0;
dt = (dx^2)/(4*alpha); % Von Neumann Stability Criterion
dt=(dx2^2)/2;

T = [t0:dt:(ceil(tf/dt)*dt)];

%% Mimetic Finite Differences Approach
%% Time discretization
t = t0;
timesteps = 1;

%% Initial Condition
u = flatten(U0, N, M);
unew = flatten(U0, N, M);

%% Mimetic Operator Matrix
order = 2;

L = lap2D(order, N, dx, M, dy);
R = robinBC2D(order, N, dx, M, dy, 0, 1);

if (explicit)
  MFD = speye(size(L)) + alpha*dt*(L-R);
else
  MFD = speye(size(L)) - alpha*dt*(L-R);
end

figure(1);
mesh(X, Y, unflatten(U0, N, M));


%% Integrate
while (t < tf+dt)
    %% update U
    if (explicit)
      unew = MFD * u;
    else
      unew = MFD \ u;
    end

    switch (equation)
      case 1
        unew = unew - loss_coeff*dt*u;
      case 2
        unew = unew - allen_coeff*dt*(u.^3-u);
    end        

    %% timestep
    timesteps = timesteps+1;
    t = timesteps*dt;
    u = unew;

    figure(2);
    colormap('cool');
   %% imagesc(unflatten(u, N, M));
    mesh(X, Y, unflatten(u, N, M));
    xlim([a b]);
    ylim([c d]);
    zlim([e f]);
    colorbar
    drawnow;
end

function v = flatten(m, N, M)
  v = reshape(m, (N+2)*(M+2), 1);
end

function m = unflatten(v, N, M)
  m = reshape(v, N+2, M+2);
end
