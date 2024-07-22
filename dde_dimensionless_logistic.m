%% Dimensionless Hutchinson Equation
% We study here the dynamics and behavior of the non-dimensionalized
% delayed differential equation (DDE) outlined in Beddington & May (1975),
% equation 13. Analysis of dynamic behavior comes from equation 14.

% The paper can be found here:
% https://www.sciencedirect.com/science/article/pii/0025556475900280

clc;
clearvars;
close all;

spacing = 0.01;

%% Run model

% Time
tf = 80;
t = linspace(1, tf, (tf - 1) / spacing);

% Lag
tau = [0 1.2 1.5 1.8];

% X for t <= 0; X = N / K
X_0 = 0.002;

% Solve & plot solutions
figure; 

for T = unique(tau)
    if T > 0
        sol = dde23(@dlhe, T, X_0, t);
    else
        sol = ode45(@he, t, X_0);
    end
    
    X = deval(sol, t);

    plot(t, smooth(X));

    i = find(tau == T, 1);
    legends{i} = sprintf('T = %g', tau(i));
    hold on
end

title('Dimensionless Hutchinson Equation Model of Population');
xlabel("Normalized Time (t')");
ylabel('Percent of Carrying Capacity (X)');
legend(legends, 'Location', 'southeast');

%% Plot dynamic behavior λ against delay τ
clearvars -except spacing;

syms lambda T
eqn = lambda == - exp(-lambda * T);
fcn = matlabFunction(rhs(eqn) - lhs(eqn));

T_max = 10;
lambdas = [];

for T = (0 + spacing):spacing:T_max

    soln = solve(fcn(T, lambda), lambda);
    lambdas = [lambdas soln];

end

T = linspace(0, T_max, T_max / spacing);

% plot stability behavior (stable or unstable)
figure; 
plot(T, real(lambdas));
xlabel('Normalized Delay (τ)');
ylabel('Re(λ)');
title('Dynamic Behahior w.r.t Normalized Delay (Stability)');
yline(0);

% plot dynamic behavior (oscillatory or monotonic)
figure;
plot(T, imag(lambdas));
xlabel('Normalized Delay (τ)');
ylabel('Im(λ)');
ylim([-0.5 2]);
title('Dynamic Behahior w.r.t Normalized Delay (Behavior)');
yline(0);

%% Define non-dimensionalized delayed logistic equations
% T > 0
function dX = dlhe(t, X, X_lag)

    dX = X_lag * (1 - X_lag);

end

% T = 0
function dX = he(t, X)

    dX = X * (1 - X);

end    