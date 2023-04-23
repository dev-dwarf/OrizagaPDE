%% Allen-Cahn and Heat Equations with Dirichlet (Zero) Boundary conditions.
clf;
clc;
clear all;
figures_so_far = 1;

%addpath('./mole_MATLAB/');

%% Settings %%%%
N = 50; % number of grid points

explicit = true;
equation = 1; % 0 = heat, 1 = heat w/ lateral loss, 2 = allen-cahn

%% Shared problem parameters

%% Heat Equation
loss_coeff = 0.15;
allen_coeff = (1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
alpha = 1.0; % thermal diffusivity.

%% X Domain
a = 0;
b = 2*pi;
dx = (b-a)/N;
dx2 = dx/2;

%% X discretization
X = [a (a+dx2):dx:(b-dx2) b]'; % solution domain is u(0) to u(N)

% Initial Condition
u0= 1.2*(rand(size(X))-1*rand(size(X)));
%% u0 = 0.15*(1*0+sin(3*2*3.14/(b-a) * X));

u0(1) = 0;
u0(end) = 0;

% Time Domain
t0 = 0;
tf = 5.0;
dt = (dx2^2)/(4*alpha); % Von Neumann Stability Criterion
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
v = alpha*dt/dx^2;
a0 = v*ones(N+2, 1);
a1 = -2*v*ones(N+2, 1);
A = spdiags([a0 a1 a0], [-1 0 1], N+2, N+2);

%% Use Saulo's non-uniform spacing scheme at boundaries
%% X = a + dx/2
A(2, 1)   = (8/3)*v; A(2, 2)     = -4*v; A(2, 3)   = (4/3)*v;
%% X = b - dx/2
A(N+1, N) = (4/3)*v; A(N+1, N+1) = -4*v; A(N+1, N+2) = (8/3)*v;

%% Enforce Dirichlet (Zero) B.Cs:
%% Remove the first and last rows, applying the rule U(a) == U(b) == 0.
A = A(2:end-1, 2:end-1);
u = u(2:end-1);
unew = u;

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
    finiteDifference(timesteps,:) = [0; u; 0];
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

%% Enforce Dirichlet (Zero) B.C:
%% Choice 1: Drop the first/last elements since they are 0.
L = L(2:end-1, 2:end-1);
u = u(2:end-1);
unew = u;

%% Choice 2: Use robinBC operator to enforce the B.C, where B.Cs
%%  are 1*u + 0*du' = g.
%% R = robinBC(order, N, dx, 1, 0);

%% Choice 3: Do nothing, which works because g in the above is 0.

if (explicit)
  MFD = speye(size(L)) + alpha*dt*(L);
else
  MFD = speye(size(L)) - alpha*dt*(L);
end
]
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
    mimetic(timesteps,:) = [0; u; 0];
end

u_mimetic = u;

%% Matlab Solver
matlab = pdepe(0, @(x, t, u, DuDx) pde(x, t, u, DuDx, equation, loss_coeff, allen_coeff), @(x) u0(find(X==x)), @bcfun, X, T);

%% Plots
figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,finiteDifference);
title("FiniteDifference");
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
xlabel('x'); ylabel('t'); zlabel('diff.');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,abs(mimetic-matlab));
title("Mimetic vs Matlab");
xlabel('x'); ylabel('t'); zlabel('diff.');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
clf; hold on;
plot(X,u0, "k:", 'DisplayName', 'I.C');
plot(X, mimetic(end,:), "r-", 'DisplayName', 'Mimetic');
plot(X, finiteDifference(end,:), "b--", 'DisplayName', 'FiniteDifference');
plot(X, matlab(end,:), "k.", 'DisplayName', 'Matlab');
switch (equation)
  case 0
  title("Final State (Heat)");
  case 1
  title("Final State (Heat w/ Lateral Loss)");
  case 2
  title("Final State (Allen-Cahn)");
end
xlabel('x')
ylabel('u')
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
      s = 0
    case 1
      s = -loss_coeff*u;
    case 2
	  s = allen_coeff*(u-u.^3);
  end        
end

function [pl,ql,pr,qr] = bcfun(xl, ul, xr, ur, t)
  pl = ul;
  ql = 0;
  pr = ur;
  qr = 0;
end
