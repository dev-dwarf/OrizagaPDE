%% Allen-Cahn and Heat Equations with No-Flux Boundary conditions.
clf;
clc;
clear all;
figures_so_far = 1;

%addpath('./mole_MATLAB/');

%% Settings %%%%
N = 25; % number of grid points

explicit = true;
equation = 0; % 0 = heat, 1 = heat w/ lateral loss, 2 = allen-cahn

%% Shared problem parameters

%% Heat Equation
loss_coeff = 0.15*10;
allen_coeff = (1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
alpha = 1.0; % thermal diffusivity.

%% X Domain
a = 0;
b = pi;
dx = (b-a)/N;
dx2 = dx/2;

%% X discretization
X = [a (a+dx2):dx:(b-dx2) b]'; % solution domain is u(0) to u(N)

% Initial Condition
u0= 1.2*(rand(size(X))-1*rand(size(X)));
u0 = 0.15*(1.0+sin(3*2*3.14/(b-a) * X));
u0 = cos(X);

% Time Domain
t0 = 0;
tf = 1.0;
dt = (dx2^2)/(4*alpha); % Von Neumann Stability Criterion
dt=(dx2^2)/2;

T = [t0:dt:(ceil(tf/dt)*dt)];

%plot_frequency = (tf/dt)/N;
plot_frequency = 2; % for watching time evolution
plot_frequency = 1;

% Implicit "Standard" Finite Differences Approach
finiteDifference = zeros(length(T), length(X)); % to store solutions

% Time discretization
t = t0;
timesteps = 1;

% Initial Condition
u = u0;
unew = u0;
finiteDifference(1,:) = u0;

%% Finite Difference Operator Matrix
%% N
v = alpha*dt/dx^2;
a0 = v*ones(N+2, 1);
a1 = -2*v*ones(N+2, 1);
A = spdiags([a0 a1 a0], [-1 0 1], N+2, N+2);

%% Use Saulo's non-uniform spacing scheme near edges:
%% X = a + dx/2 
A(2, 1)   = (8/3)*v; A(2, 2)     = -4*v; A(2, 3)   = (4/3)*v;
%% X = b - dx/2 
A(N+1, N) = (4/3)*v; A(N+1, N+1) = -4*v; A(N+1, N+2) = (8/3)*v;

%% Enforce No-Flux B.Cs:
%% Really we are using the same centered difference laplacian over the boundary,
%% however we treat the A(1,0) point as equal to A(1,2), and likewise
%% for the other edge.
vedge = 4*v;
%% X = a
A(1, 1) = -2*vedge; A(1, 2) = 2*vedge;
%% X = b
A(N+2, N+1) = 2*vedge; A(N+2, N+2) = -2*vedge;

if (explicit)
  FD = speye(size(A)) + A;
else
  FD = speye(size(A)) - A;
end

%% Integrate
while (t < tf+dt)
    %% update U
    if (explicit)
      unew = FD * u;
    else
      unew = FD \ u;
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
    finiteDifference(timesteps,:) = u;
end

u_finiteDifference = u;

%% Implicit Mimetic Finite Differences Approach
mimetic = zeros(length(T), length(X));

%% Time discretization
t = t0;
timesteps = 1;

%% Initial Condition
u = u0;
unew = u0;
mimetic(1,:) = u0;

%% Mimetic Operator Matrix
order = 2;

L = lap(order, N, dx);
R = robinBC(order, N, dx, 1, 1);

if (explicit)
  MFD = speye(size(L)) + alpha*dt*(L-R);
else
  MFD = speye(size(L)) - alpha*dt*(L-R);
end

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
    mimetic(timesteps,:) = u;
end

u_mimetic = u;

%% Matlab Solver
matlab = pdepe(0, @(x, t, u, DuDx) pde(x, t, u, DuDx, equation, loss_coeff, allen_coeff), @(x) u0(find(X==x)), @bcfun, X, T);

%% Exact Solutions
has_exact = false;
if (equation == 0)
  has_exact = true;
  for t=1:size(T,2)
    exact(t, :) = (exp(-T(t))*cos(X))';
  end
elseif (equation == 1)
  has_exact = true;
  for t=1:size(T,2)
    exact(t, :) = exp(-T(t)*(1+loss_coeff))*cos(X)';
  end
end

%% Plots
figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,finiteDifference);
title("Finite Difference");
xlabel('x'); ylabel('t'); zlabel('u');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,mimetic);
title("Mimetic");
xlabel('x'); ylabel('t'); zlabel('u');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,matlab);
title("Matlab");
xlabel('x'); ylabel('t'); zlabel('u');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,abs(mimetic-finiteDifference));
title("Mimetic vs Finite Difference");
xlabel('x'); ylabel('t'); zlabel('diff');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,abs(mimetic-matlab));
title("Mimetic vs Matlab");
xlabel('x'); ylabel('t'); zlabel('diff');

if (has_exact)
  figure(figures_so_far); figures_so_far = figures_so_far + 1;
  mesh(X,T,abs(mimetic-exact));
  title("Mimetic vs Exact");
  xlabel('x'); ylabel('t'); zlabel('diff');
end

figure(figures_so_far); figures_so_far = figures_so_far + 1;
clf; hold on;
plot(X,u0, "k:", 'DisplayName', 'I.C');
plot(X, mimetic(end,:), "r-", 'DisplayName', 'Mimetic');
plot(X, finiteDifference(end,:), "b--", 'DisplayName', 'FiniteDifference');
plot(X, matlab(end,:), "k.", 'DisplayName', 'Matlab');

if (has_exact)
  plot(X, exact(end,:), "k-.", 'DisplayName', 'Exact');
end

switch (equation)
  case 0
  title("Final State (Heat)");
  case 1
  title("Final State (Heat w/ Lateral Loss)");
  case 2
  title("Final State (Allen-Cahn)");
end

xlabel('x');
ylabel('u');
xlim([a b]);
ylim([-2 2]);
legend();
hold off;

%% Function for matlab solver
function [c,f,s] = pde(x, t, u, DuDx, equation, loss_coeff, allen_coeff)
  c = 1;
  f = DuDx;
  switch (equation)
    case 0
      s = 0;
    case 1
      s = -loss_coeff*u;
    case 2
	  s = allen_coeff*(u-u.^3);
  end        
end

function [pl,ql,pr,qr] = bcfun(xl, ul, xr, ur, t)
  pl = 0;
  ql = 1;
  pr = 0;
  qr = 1;
end

