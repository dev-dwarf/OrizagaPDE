%% Allen-Cahn and Heat Equations with No-Flux Boundary conditions.
clf;
clc;
clear all;

addpath('./mole_MATLAB/');

%% Settings %%%%
N = 50; % number of grid points

explicit = true;
plots = true;
skip_standard = false;
allen_cahn = true;

%% Shared problem parameters

%% Heat Equation
coeff=(1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
alpha = 1.0; % thermal diffusivity.

%% X Domain
a = -2*pi;
b = 2*pi;
dx = (b-a)/N;
dx2 = dx/2;

%% X discretization
x = [a (a+dx2):dx:(b-dx2) b]'; % solution domain is u(0) to u(N)

% Initial Condition
u0= 1.2*(rand(size(x))-1*rand(size(x)));
%% u0 = 0.5*(1+sin(2*2*3.14/(b-a) * x));

% Time Domain
t0 = 0;
tf = 5.0;
dt = (dx2^2)/(4*alpha); % Von Neumann Stability Criterion

%plot_frequency = (tf/dt)/N;
plot_frequency = 2; % for watching time evolution
plot_frequency = 1;

% Implicit "Standard" Finite Differences Approach
standard = zeros(ceil(tf/dt), length(x)); % to store solutions

% Time discretization
t = t0;
timesteps = 1;

% Initial Condition
u = u0;
unew = u0;
standard(1,:) = u0;

%% Finite Difference Operator Matrix
%% N
v = alpha*dt/dx^2;
a0 = v*ones(N+2, 1);
a1 = -2*v*ones(N+2, 1);
A = spdiags([a0 a1 a0], [-1 0 1], N+2, N+2);

%% Use Saulo's non-uniform spacing scheme at boundaries
%% MODIFICATION FOR NO-FLUX B.Cs:
%% For periodic we had A(1, 0) = A(1, N+1)
%% Now we will have A(1, 0) = A(1, 1)
%% Hence we add the coefficient from A(1, N+1) to the entry for A(1,1)
%% X = a
A(1, 1) = -v; % = -2*v + v
A(1, 2) = v;

%% X = a + dx/2 (Same as for periodic case)
A(2, 1)   = (8/3)*v; A(2, 2)     = -4*v; A(2, 3)   = (4/3)*v;

%% X = b - dx/2 (N+2 no longer "wraps-around" to 1)
A(N+1, N) = (4/3)*v; A(N+1, N+1) = -4*v; A(N+1, N+2) = (8/3)*v;

%% X = b (new row, not needed for periodic because U(b) equaled U(a) before)
A(N+2, N+1) = v; A(N+2, N+2) = -v; % = -2*v + v

if (explicit)
  FD = speye(size(A)) + A;
else
  FD = speye(size(A)) - A;
end

if (plots)
  figure(1);
  plot(x, u, 'o-');
  str = sprintf('Initial Condition \t t = 0', t);
  title(str)
  xlabel('x')
  ylabel('U')
  xlim([a b]);
  ylim([-2 2]);
end

figures_so_far = 2;

%% Integrate
while (t <= tf)
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
	
    %% Make a plot every few timesteps and on the last timestep.
    if (plots && (mod(timesteps-1, plot_frequency) == 1 || t >= tf))
      % figure(figures_so_far); figures_so_far = figures_so_far + 1;
      figure(2);

      plot(x, u, 'o-');
      str = sprintf('Standard \t t = %.2f', t);
      title(str)
      xlabel('x')
      ylabel('U')
      xlim([a b]);
      ylim([-2 2]);
    end

    standard(timesteps,:) = u;
end

u_standard = u;

%% Implicit Mimetic Finite Differences Approach
mimetic = zeros(ceil(tf/dt), length(x));

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
L = L - R;

if (explicit)
  MFD = speye(size(L)) + alpha*dt*L;
else
  MFD = speye(size(L)) - alpha*dt*L;
end

%% Integrate
while (t <= tf)
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

    %% Make a plot every few timesteps and on the last timestep.
    if (plots && (mod(timesteps-1, plot_frequency) == 1 || t >= tf))
      %figure(figures_so_far); figures_so_far = figures_so_far + 1;
      figure(3);

      plot(x, u, 'o-');
      str = sprintf('Mimetic \t t = %.2f', t);
      title(str)
      xlabel('x')
      ylabel('U')
      xlim([a b]);
      ylim([-2 2]);
    end

    mimetic(timesteps,:) = u;
end

u_mimetic = u;

abs_diff = abs(u_standard - u_mimetic);

if (plots)
  figure(4);
  plot(x, abs_diff, 'o-');
  str = sprintf('Abs. Diff');
  title(str)
  xlabel('x')
  ylabel('dff')
  xlim([a b]);
  ylim([-2 2]);
end

%% 3D plots
if (plots)
  X = x;
  T = 0:dt:tf;

  %figure(figures_so_far); figures_so_far + 1;
  figure(5);
  mesh(X,T,standard);
  title("standard");
  xlabel('x'); ylabel('t'); zlabel('u');

  %figure(figures_so_far); figures_so_far + 1;
  figure(6);
  mesh(X,T,mimetic);
  title("mimetic");
  xlabel('x'); ylabel('t'); zlabel('u');

  %figure(figures_so_far); figures_so_far + 1;
  figure(7);
  mesh(X,T,abs(mimetic-standard));
  title("abs. difference");
  xlabel('x'); ylabel('t'); zlabel('diff');

end

