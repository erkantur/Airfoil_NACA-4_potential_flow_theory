% Define the NACA 4-digit airfoil
N = 200; % Number of panels
m = 0.02; % Max camber
p = 0.4;  % Location of max camber
t = 0.12; % Thickness-to-chord ratio
c = 1;    % Chord length

% Generate airfoil coordinates using NACA 4-digit formula
[x, y] = naca4(m, p, t, c, N);
x = x(:)';
y = y(:)';

% Angle of attack and freestream velocity
alpha = 5; % degrees
U_inf = 100; % m/s

% Panel method calculations
[theta, control_points] = panel_angles(x, y, N);
[A, RHS] = panel_system(theta, control_points, alpha, U_inf, N);
gamma = A \ RHS;

% Pressure coefficient
V_tangent = U_inf * sin(theta - deg2rad(alpha)) - gamma / (2*pi);
Cp = 1 - (V_tangent / U_inf).^2;

% Lift and empirical drag coefficients
Cl = 2 * sum(gamma) / (U_inf * c);
Cf = 0.074 / (Reynolds_number(U_inf, c, 1.225, 1.81e-5))^0.2; % Friction coefficient
Cd = Cf * 1.328 / sqrt(Reynolds_number(U_inf, c, 1.225, 1.81e-5)); % Empirical drag coefficient

% Plot the results
figure
subplot(2,1,1);
plot(x, y, 'b-');
axis equal
xlabel('X Coordinate')
ylabel('Y Coordinate')
title('NACA Airfoil Geometry')
grid on

subplot(2,1,2);
plot(control_points(:,1), Cp, 'r-');
xlabel('X Coordinate')
ylabel('Pressure Coefficient, Cp')
title('Pressure Distribution over Airfoil')
grid on

fprintf('Lift Coefficient, Cl: %.4f\n', Cl);
fprintf('Drag Coefficient, Cd: %.4f\n', Cd);

% Helper functions
function [theta, control_points] = panel_angles(x, y, N)
    control_points = zeros(N, 2);
    theta = zeros(N, 1);
    for i = 1:N-1
        control_points(i,:) = [0.5*(x(i)+x(i+1)), 0.5*(y(i)+y(i+1))];
        theta(i) = atan2((y(i+1)-y(i)),(x(i+1)-x(i)));
    end
    control_points(N,:) = [0.5*(x(N)+x(1)), 0.5*(y(N)+y(1))];
    theta(N) = atan2((y(1)-y(N)),(x(1)-x(N)));
end

function [A, RHS] = panel_system(theta, control_points, alpha, U_inf, N)
    A = zeros(N, N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                r = control_points(i,:) - control_points(j,:);
                A(i,j) = sin(theta(i)-theta(j))*(norm(r)^2 - dot(r,[cos(theta(j)) sin(theta(j))])) / ...
                          (norm(r)^2*pi) + cos(theta(i)-theta(j))*dot(r,[sin(theta(j)) -cos(theta(j))]) / ...
                          (norm(r)*pi);
            end
        end
    end
    RHS = U_inf * cos(theta - deg2rad(alpha));
end

function Re = Reynolds_number(U, c, rho, mu)
    Re = rho * U * c / mu;
end

function [x, y] = naca4(m, p, t, c, N)
    x = linspace(0, c, N);
    y_t = 5 * t * c * (0.2969 * sqrt(x/c) - 0.1260 * (x/c) - 0.3516 * (x/c).^2 + 0.2843 * (x/c).^3 - 0.1036 * (x/c).^4);
    if p == 0
        y_c = zeros(size(x));
        dyc_dx = zeros(size(x));
    else
        y_c = m * (2 * p * (x/c) - (x/c).^2) / p^2 .* (x <= p*c) + m * (1 - 2 * p + 2 * p * (x/c) - (x/c).^2) / (1 - p)^2 .* (x > p*c);
        dyc_dx = 2 * m * (1 / p - 1) / p^2 .* (x <= p*c) + 2 * m * (-1 / (1 - p) + 1) / (1 - p)^2 .* (x > p*c);
    end
    theta = atan(dyc_dx);
    x_u = x - y_t .* sin(theta);
    y_u = y_c + y_t .* cos(theta);
    x_l = x + y_t .* sin(theta);
    y_l = y_c - y_t .* cos(theta);
    x = [flip(x_u), x_l(2:end)];
    y = [flip(y_u), y_l(2:end)];
end
