%% Allen-Cahn and Heat Equations with No-Flux Boundary conditions.
clf;
clc;
clear all;
figures_so_far = 1;

%addpath('./mole_MATLAB/');

%% Settings %%%%
N = 50; % number of grid points

explicit = true;
allen_cahn = true;

%% Shared problem parameters

%% Heat Equation
coeff = (1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
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
u0 = 0.15*(1*0+sin(3*2*3.14/(b-a) * X));

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
%% N
v = alpha*dt/dx^2;
a0 = v*ones(N+2, 1);
a1 = -2*v*ones(N+2, 1);
A = spdiags([a0 a1 a0], [-1 0 1], N+2, N+2);

%% Use Saulo's non-uniform spacing scheme at boundaries
%% MODIFICATION FOR NO-FLUX B.Cs:
%% For periodic we had A(1, 0) = A(1, N+1)
%% Now we will have A(1, 0) = A(1, 2)
%% Hence we add the coefficient from A(1, N+1) to the entry for A(1,2)
%% X = a
A(1, 1) = -2*v; 
A(1, 2) = 2*v;

%% X = a + dx/2 (Same as for periodic case)
A(2, 1)   = (8/3)*v; A(2, 2)     = -4*v; A(2, 3)   = (4/3)*v;

%% X = b - dx/2 (N+2 no longer "wraps-around" to 1)
A(N+1, N) = (4/3)*v; A(N+1, N+1) = -4*v; A(N+1, N+2) = (8/3)*v;

%% X = b (new row, not needed for periodic because U(b) equaled U(a) before)
A(N+2, N+1) = 2*v; A(N+2, N+2) = -2*v;

if (explicit)
  FD = speye(size(A)) + A;
else
  FD = speye(size(A)) - A;
end

figure(figures_so_far); figures_so_far = figures_so_far + 1;
plot(X, u, 'o-');
str = sprintf('Initial Condition \t t = 0', t);
title(str)
xlabel('x')
ylabel('U')
xlim([a b]);
ylim([-2 2]);

%% Integrate
while (t < tf)
    %% update U
    if (explicit)
      unew = FD * u;
    else
      unew = FD \ u;
    end

    if (allen_cahn)
      unew = unew - coeff*dt*(u.^3-u);
    end

    %% timestep
    timesteps = timesteps+1;
    t = timesteps*dt;

    u = unew;
    finiteDifference(timesteps,:) = u;
end
finiteDifference(end,:) = u;

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,finiteDifference);
title("FiniteDifference");
xlabel('x'); ylabel('t'); zlabel('u');

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
R = robinBC(order, N, dx, 0, 1);

if (explicit)
  MFD = speye(size(L)) + alpha*dt*(L-R);
else
  MFD = speye(size(L)) - alpha*dt*(L-R);
end

%% Integrate
while (t < tf)
    %% update U
    if (explicit)
      unew = MFD * u;
    else
      unew = MFD \ u;
    end

    if (allen_cahn)
      unew = unew - coeff*dt*(u.^3-u);
    end

    %% timestep
    timesteps = timesteps+1;
    t = timesteps*dt;

    u = unew;
    mimetic(timesteps,:) = u;
end
mimetic(end,:) = u;

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,mimetic);
title("Mimetic");
xlabel('x'); ylabel('t'); zlabel('u');

u_mimetic = u;

%% Matlab Solver
matlab = pdepe(0, @(x, t, u, DuDx) pde(x, t, u, DuDx, allen_cahn, coeff), @(x) u0(find(X==x)), @bcfun, X, T);

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,matlab);
title("Matlab");
xlabel('x'); ylabel('t'); zlabel('u');


%% 3D plots
figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,abs(mimetic-finiteDifference));
title("Mimetic vs FiniteDifference");
xlabel('x'); ylabel('t'); zlabel('u');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,abs(mimetic-matlab));
title("Mimetic vs Matlab");
xlabel('x'); ylabel('t'); zlabel('u');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
hold on;
plot(X,u0, "c:", 'DisplayName', 'I.C');
plot(X, mimetic(end,:), "r-", 'DisplayName', 'Mimetic');
plot(X, finiteDifference(end,:), "b--", 'DisplayName', 'FiniteDifference');
plot(X, matlab(end,:), "k.", 'DisplayName', 'Matlab');
title("Final State")
xlabel('x')
ylabel('U')
xlim([a b]);
ylim([-2 2]);
legend();
hold off;

%% Function for matlab solver
function [c,f,s] = pde(x, t, u, DuDx, allen_cahn, coeff)
  c = 1;
  f = DuDx;
  if (allen_cahn)
	s = coeff*(u-u.^3);
  else
	s = 0;
  end
end

function [pl,ql,pr,qr] = bcfun(xl, ul, xr, ur, t)
  pl = 0;
  ql = 1;
  pr = 0;
  qr = 1;
end

  

