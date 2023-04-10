%% Allen-Cahn and Heat Equations with Periodic Boundary conditions.
clf;
clc;
clear all;

addpath('./mole_MATLAB/');

%% Settings %%%%
N = 50; % number of grid points

explicit = false;
plots = true;
skip_standard = false;
allen_cahn = true;

%% Shared problem parameters

%% Heat Equation
coeff=(1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
alpha = 1.0; % thermal diffusivity.

%% X Domain
a = 0;
b = 2*pi;
dx = (b-a)/N;
dx2 = dx/2;

%% X discretization
x = [a (a+dx2):dx:(b-dx2)]'; % solution domain is u(0) to u(N) {u(x(0)) = u(x(N+1))}
x_display = [x; b];

%% Initial Condition
if (allen_cahn)
  u0=1.2*(rand(size(x))-1*rand(size(x)));
  %u0 = 0.5*(1+sin(5*3.14/(b-a) * x));
else
  u0 = 0.5*(1+sin(2*3.14/(b-a) * x));
end

%% Time Domain
t0 = 0;
tf = 0.4;
dt = (dx2^2)/(4*alpha); % Von Neumann Stability Criterion

%plot_frequency = (tf/dt)/N;
plot_frequency = 2; % for watching time evolution
plot_frequency = 1;

%% Implicit "Standard" Finite Differences Approach
standard = zeros(ceil(tf/dt), length(x_display)); % to store solutions

%% Time discretization
t = t0;
timesteps = 1;

%% Initial Condition
u = u0;
unew = u0;
u_display = [u0; u0(1)];
standard(1,:) = u_display;

%% Finite Difference Operator Matrix
v = alpha*dt/dx^2;
a0 = v*ones(N+1, 1);
a1 = -2*v*ones(N+1, 1);
A = spdiags([a0 a1 a0], [-1 0 1], N+1, N+1);

%% Use Saulo's non-uniform spacing scheme at boundaries
%% X = a
A(1, N+1) = v;       A(1, 1)     = -2*v; A(1, 2)   = v;
%% X = a + dx/2
A(2, 1)   = (8/3)*v; A(2, 2)     = -4*v; A(2, 3)   = (4/3)*v;
%% X = b - dx/2
A(N+1, N) = (4/3)*v; A(N+1, N+1) = -4*v; A(N+1, 1) = (8/3)*v;

if (explicit)
  FD = speye(size(A)) + A;
else
  FD = speye(size(A)) - A;
end


if (plots)
  figure(1);
  plot(x_display, u_display, 'o-');
  str = sprintf('Initial Condition \t t = 0', t);
  title(str)
  xlabel('x')
  ylabel('U')
  xlim([a b]);
  ylim([-2 2]);
end

figures_so_far = 2;

if (skip_standard == 0)

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

    u = unew;

    %% timestep
    timesteps = timesteps+1;
    t = timesteps*dt;

    %% Make a plot every few timesteps and on the last timestep.
    u_display = [u; u(1)];
    if (plots && (mod(timesteps-1, plot_frequency) == 1 || t >= tf))
      % figure(figures_so_far); figures_so_far = figures_so_far + 1;
      figure(2);

      plot(x_display, u_display, 'o-');
      str = sprintf('Standard \t t = %.2f', t);
      title(str)
      xlabel('x')
      ylabel('U')
      xlim([a b]);
      ylim([-2 2]);
    end

    standard(timesteps,:) = u_display;
end
end

u_standard = u;

%% Implicit Mimetic Finite Differences Approach
mimetic = zeros(ceil(tf/dt), length(x_display));

%% X discretization
x = [a (a+dx/2):dx:(b-dx/2)]';
x_display = [x; b];

%% Time discretization
t = t0;
timesteps = 0;

%% Initial Condition
u = u0;
unew = u0;
u_display = [ u; u(1)];
mimetic(1,:) = u_display;

%% Mimetic Operator Matrix
order = 2;

%% Instead of using MOLE's lap(), construct laplacian from div and grad ourselves.
%% This lets us impose boundary conditions on D and G as needed.
D = div(order, N, dx); % 1D Mimetic divergence operator
G = grad(order, N, dx);

%% Periodic BC imposed on the divergence operator
D(1,2) = 1/(2*dx);
D(1,end-1) = -1/(2*dx);
D(end,2) = 1/(2*dx);
D(end,end-1) = -1/(2*dx);

%% Apply the rule that U(1) == U(end), and drop last row/column.
L = D*G;
L(:, 1) = L(:, 1) + L(:, end);
L = L(1:end-1,1:end-1);

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

    u_display = [u; u(1)];
    %% Make a plot every few timesteps and on the last timestep.
    if (plots && (mod(timesteps-1, plot_frequency) == 1 || t >= tf))
      %figure(figures_so_far); figures_so_far = figures_so_far + 1;
      figure(3);

      plot(x_display, u_display, 'o-');
      str = sprintf('Mimetic \t t = %.2f', t);
      title(str)
      xlabel('x')
      ylabel('U')
      xlim([a b]);
      ylim([-2 2]);
    end

    mimetic(timesteps,:) = u_display;
end

u_mimetic = u;

percent_diff = 100 * abs(u_standard - u_mimetic) ./ u_standard;

%% 3D plots
if (plots)
  X = x_display;
  T = 0:dt:tf;

  %figure(figures_so_far); figures_so_far + 1;
  figure(4);
  mesh(X,T,standard);
  title("standard");
  xlabel('x'); ylabel('t'); zlabel('u');

  %figure(figures_so_far); figures_so_far + 1;
  figure(5);
  mesh(X,T,mimetic);
  title("mimetic");
  xlabel('x'); ylabel('t'); zlabel('u');

  %figure(figures_so_far); figures_so_far + 1;
  figure(6);
  mesh(X,T,abs(mimetic-standard));
  title("abs. difference");
  xlabel('x'); ylabel('t'); zlabel('diff');

end

