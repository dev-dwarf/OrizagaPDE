%% Allen-Cahn and Heat Equations with Periodic Boundary conditions.
clf;
clc;
clear all;
figures_so_far = 1;

addpath('./mole_MATLAB/');

%% Settings %%%%
N = 50; % number of grid points

explicit = false;
allen_cahn = false;

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
finiteDifference = zeros(ceil(tf/dt), length(x_display)); % to store solutions

%% Time discretization
t = t0;
timesteps = 1;

%% Initial Condition
u = u0;
unew = u0;
u_display = [u0; u0(1)];
finiteDifference(1,:) = u_display;

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

%% Enforce Periodic B.Cs:
%% Use centered difference laplacian.
vedge = 4*v;
%% X = a
A(1, N+1) = vedge;       A(1, 1)     = 2*vedge; A(1, 2)   = v*vedge;
%% Apply the rule that U(a) == U(b), and drop last row/column.
A(:, 1) = A(:, 1) + A(:, end);
A = A(1:end-1,1:end-1);


if (explicit)
  FD = speye(size(A)) + A;
else
  FD = speye(size(A)) - A;
end

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

    u_display = [u; u(1)];
    finiteDifference(timesteps,:) = u_display;
end

u_finiteDifference = u;

%% Implicit Mimetic Finite Differences Approach
mimetic = zeros(ceil(tf/dt), length(x_display));

%% X discretization
x = [a (a+dx/2):dx:(b-dx/2)]';
x_display = [x; b];

%% Time discretization
t = t0;
timesteps = 1;

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
%% Taken from MOLE example code: examples_MATLAB/hyperbolic1D.m
%% https://github.com/jcorbino/mole/blob/686bfe038c5717957b67496dc11c552f872609f2/examples_MATLAB/hyperbolic1D.m
D(1,2) = 1/(2*dx);
D(1,end-1) = -1/(2*dx);
D(end,2) = 1/(2*dx);
D(end,end-1) = -1/(2*dx);

%% Construct Mimetic Laplacian
L = D*G;

%% Apply the rule that U(a) == U(b), and drop last row/column.
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
    mimetic(timesteps,:) = u_display;
end

u_mimetic = u;

percent_diff = 100 * abs(u_finiteDifference - u_mimetic) ./ u_finiteDifference;

%% Plots
X = x_display;
T = 0:dt:tf;

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,finiteDifference);
title("Finite Difference");
xlabel('x'); ylabel('t'); zlabel('u');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,mimetic);
title("Mimetic");
xlabel('x'); ylabel('t'); zlabel('u');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
mesh(X,T,abs(mimetic-finiteDifference));
title("Mimetic vs Finite Difference");
xlabel('x'); ylabel('t'); zlabel('diff.');

figure(figures_so_far); figures_so_far = figures_so_far + 1;
clf; hold on;
plot(X, [u0; u0(1)], "c:", 'DisplayName', 'I.C');
plot(X, mimetic(end,:), "r-", 'DisplayName', 'Mimetic');
plot(X, finiteDifference(end,:), "b--", 'DisplayName', 'FiniteDifference');
if (allen_cahn)
  title("Final State (Allen-Cahn)");
else
  title("Final State (Heat)");
end
xlabel('x')
ylabel('u')
xlim([a b]);
ylim([-2 2]);
legend();
hold off;

