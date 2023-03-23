% Heat Equation with periodic boundary conditions.
clf;
clc;
clear all;

%addpath('./mole_MATLAB/')

%Settings %%%%
N = 50; % number of grid points

explicit = 1;
plots = true;
plot_frequency = 20 * N;
skip_standard = false;
allen_cahn = true;

% Shared problem parameters

% Heat Equation
coeff=(1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
alpha = 1.0; % thermal diffusivity.

% X Domain
a = 0;
b = 2*pi;
dx = (b-a)/N;

% Time Domain
t0 = 0;
tf = 0.5;
dt = dx^2/(4*alpha); % Von Neumann Stability Criterion

% Implicit "Standard" Finite Differences Approach

% X discretization
dx = (b-a)/N;
x = [0:dx:(b-dx)]'; % solution domain is u(0) to u(N) {u(x(0)) = u(x(N+1))}
x_display = [x; b];

% Time discretization
t = t0;

timesteps = 0;

% Initial Condition
if (allen_cahn)
  u0=1.2*(rand(size(x))-1*rand(size(x)));
else
  u0 = sin(3.14 * x);
end

u_display = [u0; u0(1)];
u = u0;
unew = u0;

% Finite Difference Operator Matrix
v = alpha*dt/dx^2;
a0 = v*ones(N, 1);
a1 = -2*v*ones(N, 1);
A = spdiags([a0 a1 a0], [-1 0 1], N, N);
A(1, N) = v;
A(N, 1) = v;

if (explicit)
  FD = speye(size(A)) + A;
else
  FD = speye(size(A)) - A;
end

figures_so_far = 1;

if (skip_standard == 0)

% Integrate
while (t < tf)
    % update U
    if (explicit)
      unew = FD * u;
    else
      unew = FD \ u;
    end

    if (allen_cahn)
      unew = unew - coeff*dt*(u.^3-u);
    end

    % timestep
    timesteps = timesteps+1;
    t = timesteps*dt;

    % Make a plot every few timesteps and on the last timestep.
    if (plots && (mod(timesteps-1, plot_frequency) == 0 || t >= tf))
      u_display = [u; u(1)];
      figure(figures_so_far);
      figures_so_far = figures_so_far + 1;
      plot(x_display, u_display, 'o-');
      str = sprintf('Standard \t t = %.2f', t);
      title(str)
      xlabel('x')
      ylabel('U')
      xlim([a b]);
      ylim([-1 1]);
    end

    u = unew;
end
end


% Implicit Mimetic Finite Differences Approach

% X discretization
x = [a (a+dx/2):dx:(b-dx/2)]';
x_display = [x; b];

% Time discretization
t = t0;
timesteps = 0;

% Initial Condition
if (allen_cahn)
  u0=1.2*(rand(size(x))-1*rand(size(x)));
else
  u0 = sin(3.14 * x);
end
u_display = [u0; u0(1)];
u = u0;
unew = u0;

% Mimetic Operator Matrix
order = 2;

% Instead of using MOLE's lap(), construct laplacian from div and grad ourselves.
% This lets us impose boundary conditions on D and G as needed.
D = div(order, N, dx); % 1D Mimetic divergence operator
G = grad(order, N, dx);

% Periodic BC imposed on the divergence operator
D(1,2) = 1/(2*dx);
D(1,end-1) = -1/(2*dx);
D(end,2) = 1/(2*dx);
D(end,end-1) = -1/(2*dx);

% Periodic BC imposed on the gradient operator
temp = G(end, :);
G(end,:) = G(end,:) - G(1,:);
G(1,:) = G(1,:) - temp;

% Apply the rule that U(1) == U(end), and drop last row/column.
L = D*G;
L(:, 1) = L(:, 1) + L(:, end);
L = L(1:end-1,1:end-1);

if (explicit)
  MFD = speye(size(L)) + alpha*dt*L;
else
  MFD = speye(size(L)) - alpha*dt*L;
end

% Integrate
coeff=(1/.15).^2;  %Coefficient needed for Allen-Cahn PDE.
while (t < tf)
    % update U
    if (explicit)
      unew = MFD * u;
    else
      unew = MFD \ u;
    end

    if (allen_cahn)
      unew = unew - coeff*dt*(u.^3-u);
    end

    % timestep
    timesteps = timesteps+1;
    t = timesteps*dt;

    % Make a plot every few timesteps and on the last timestep.
    if (plots && (mod(timesteps-1, plot_frequency) == 0 || t >= tf))
      u_display = [u; u(1)];
      figure(figures_so_far);
      figures_so_far = figures_so_far + 1;
      plot(x_display, u_display, 'o-');
      str = sprintf('Mimetic \t t = %.2f', t);
      title(str)
      xlabel('x')
      ylabel('U')
      xlim([a b]);
      ylim([-1 1]);
    end

    u = unew;
end


