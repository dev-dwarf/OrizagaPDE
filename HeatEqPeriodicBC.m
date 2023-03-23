## Heat Equation with periodic boundary conditions.
clf;
clc;
clear all;

addpath('./mole_MATLAB/')

#### Settings #####
N = 50; # number of grid points

explicit = 1;
plots = true;
plot_frequency = 5 * N;
skip_standard = false;


##################################

## Shared problem parameters #####

## Heat Equation
alpha = 0.5; # thermal diffusivity.

# X Domain
a = 0;
b = 1;


# Time Domain
t0 = 0;
tf = 0.5;

##################################

## Implicit "Standard" Finite Differences Approach

# X discretization
dx = (b-a)/N;
x = [0:dx:(b-dx)]'; # solution domain is u(0) to u(N) {u(x(0)) = u(x(N+1))}
x_display = [x; b];

# Time discretization
t = t0;
dt = (dx^2) / (4 * alpha); # Use same timestep as mimetic code.

timesteps = 0;

# Initial Condition
u0 = sin(3.14 * x);
u_display = [u0; u0(1)];
u = u0;
unew = u0;

# Finite Difference Operator Matrix
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
endif

figures_so_far = 1;

if (skip_standard == 0)
# Integrate
while (t < tf)
    # update U
    if (explicit)
      unew = FD * u;
    else
      unew = FD \ u;
    endif

    # timestep
    timesteps = timesteps+1;
    t = timesteps*dt;

    # Make a plot every few timesteps and on the last timestep.
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
      ylim([0 1]);
    endif

    u = unew;
end
endif
##################################


## Implicit Mimetic Finite Differences Approach

# X discretization
N = (N ); # chosen so FD and MFD matrix have same dimensions
dx = (b-a)/N;
x = [a (a+dx/2):dx:(b-dx/2) b]';
x_display = [x; b];

# Time discretization
t = t0;
timesteps = 0;
dt = dx^2/(4*alpha); # Von Neumann Stability Criterion

# Initial Condition
u0 = sin(3.14 * x);
u0(1) = 0; # first/last index holds B.C, we want 0.
u0(end) = 0;
u_display = [u0; u0(1)];
u = u0;
unew = u0;

# Mimetic Operator Matrix
order = 2;

# Instead of using MOLE's lap(), construct laplacian from div and grad ourselves.
# This lets us impose boundary conditions on D and G as needed.
D = div(order, N, dx); % 1D Mimetic divergence operator
G = grad(order, N, dx);

% Periodic BC imposed on the divergence operator
D(1,2) = 1/(2*dx);
D(1,end-1) = -1/(2*dx);
D(end,2) = 1/(2*dx);
D(end,end-1) = -1/(2*dx);

# Periodic BC imposed on the gradient operator
temp = G(end, :);
G(end,:) = G(end,:) - G(1,:);
G(1,:) = G(1,:) - temp;

L = D*G;

if (explicit)
  MFD = speye(size(L)) + alpha*dt*L;
else
  MFD = speye(size(L)) - alpha*dt*L;
endif

# Integrate
while (t < tf)
    # update U
    if (explicit)
      unew = MFD * u;
    else
      unew = MFD \ u;
    endif

    # timestep
    timesteps = timesteps+1;
    t = timesteps*dt;

    # Make a plot every few timesteps and on the last timestep.
    if (plots && (mod(timesteps-1, plot_frequency) == 0 || t >= tf))
      figure(figures_so_far);
      figures_so_far = figures_so_far + 1;
      plot(x, u, 'o-');
      str = sprintf('Mimetic \t t = %.2f', t);
      title(str)
      xlabel('x')
      ylabel('U')
      xlim([a b]);
      ylim([0 1]);
    endif

    u = unew;
end

##################################
