clf;
clc;
clear all;
fig = 1;

%% Heat Equation
for n=0:8
  load(sprintf('./data/Heat_NoFlux_dt%d.mat', n));
  times(n+1) = log2(dt);
  
  %% Makes sense to average the norms by number of points
  mimeticError(n+1) = (norm(exact-mimetic, 2)/length(exact));
  finiteDifferenceError(n+1) = (norm(exact-finiteDifference, 2)/length(exact));

end
figure(fig); fig = fig + 1;
hold on;
plot(times, finiteDifferenceError, 'ko--', 'DisplayName', 'Finite Difference');
plot(times, mimeticError, 'ro--', 'DisplayName', 'Mimetic');
xlabel('log2(dt)'); ylabel('Normalized Mean Squared Error');
title('Error Comparison for Heat Eq. w/ No Flux B.Cs');
legend();
hold off;
save("./data/Error_Heat_NoFlux.mat", "dt", "mimeticError", "finiteDifferenceError");
