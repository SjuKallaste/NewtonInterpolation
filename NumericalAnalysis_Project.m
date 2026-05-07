clc; clear; close all;

% <problem data>
fprintf('\n%s\n', repmat('-', 1, 52));
fprintf('%s\n', '  Problem Setup');
fprintf('%s\n', '  Press ENTER to accept the default value.');
fprintf('%s\n\n', repmat('-', 1, 52));

raw = input(sprintf('  x nodes  [default: [1 1.3 1.6 1.9 2.2]] : '), 's');
if isempty(raw), x_data = [1, 1.3, 1.6, 1.9, 2.2];
else, x_data = str2num(raw); end %#ok<ST2NM>

raw = input(sprintf('  f(x) vals [default: [0.1411 -0.6878 -0.9962 -0.5507 0.3115]] : '), 's');
if isempty(raw), y_data = [0.1411, -0.6878, -0.9962, -0.5507, 0.3115];
else, y_data = str2num(raw); end %#ok<ST2NM>

raw = input(sprintf('  f(t) expr [default: sin(3*t)] : '), 's');
if isempty(raw), raw = 'sin(3*t)'; end
f_expr = raw;
f_exact = str2func(['@(t) ' f_expr]);

raw = input(sprintf('  Eval point x [default: 1.5] : '), 's');
if isempty(raw), x_eval = 1.5;
else, x_eval = str2double(raw); end

raw = input(sprintf('  Integral lower bound [default: 1] : '), 's');
if isempty(raw), int_lo = 1;
else, int_lo = str2double(raw); end

raw = input(sprintf('  Integral upper bound [default: 1.6] : '), 's');
if isempty(raw), int_hi = 1.6;
else, int_hi = str2double(raw); end

n = length(x_data);
fprintf('\n');
% </problem data>

% <divided difference table construction>
D = zeros(n, n);
D(:, 1) = y_data(:);
for j = 2:n
    for i = 1:(n - j + 1)
        D(i, j) = (D(i+1, j-1) - D(i, j-1)) / (x_data(i+j-1) - x_data(i));
    end
end
% </divided difference table construction>

% <print divided difference table>
CW = 12;
TW = 8 + CW * n;
fprintf('\n%s\n', repmat('=', 1, TW));
fprintf('%s\n', sprintf('%*s', floor(TW/2) + 14, 'DIFF TABLE'));
fprintf('%s\n', repmat('=', 1, TW));
fprintf('%-8s', 'x');
for j = 0:n-1
    fprintf('%-*s', CW, sprintf('Diff #%d', j));
end
fprintf('\n%s\n', repmat('-', 1, TW));
for i = 1:n
    fprintf('%-8.4f', x_data(i));
    for j = 1:(n - i + 1)
        fprintf('%-*.4f', CW, D(i, j));
    end
    fprintf('\n');
end
fprintf('%s\n\n', repmat('=', 1, TW));
% </print divided difference table>

% <user menu>
choice = menu('Select The Method Used', ...
    'Linear (Forward)', ...
    'Quadratic (Forward)', ...
    'Cubic (Forward)', ...
    'Linear (Backward)', ...
    'Quadratic (Backward)', ...
    'Cubic (Backward)');
if choice == 0
    fprintf('%s\n', 'No method selected. Exiting...');
    return;
end
% </user menu>

% <method configuration>
degree_names = {'Linear', 'Quadratic', 'Cubic'};
if choice <= 3
    degree = choice;
    is_fwd = true;
else
    degree = choice - 3;
    is_fwd = false;
end
if is_fwd
    dir_str = 'Forward';
    used_x = x_data(1:degree+1);
    coeffs = D(1, 1:degree+1);
    omega_nodes = x_data(1:degree);
else
    dir_str = 'Backward';
    used_x = x_data(n-degree:n);
    coeffs = arrayfun(@(k) D(n - k + 1, k), 1:degree+1);
    omega_nodes = x_data(n:-1:n-degree+1);
end
method_str = sprintf('%s (%s)', degree_names{degree}, dir_str);
% </method configuration>

% <convert Newton form to standard polynomial>
p_std = zeros(1, degree + 1);
basis = 1;
for k = 1:degree + 1
    p_std = p_std + [zeros(1, degree + 1 - k), coeffs(k) * basis];
    if k <= degree
        basis = conv(basis, [1, -omega_nodes(k)]);
    end
end
% </convert Newton form to standard polynomial>

% <evaluate interpolation at x_eval>
P_x = polyval(p_std, x_eval);
f_x = f_exact(x_eval);
abs_err = abs(P_x - f_x);
rel_err = 100 * abs_err / abs(f_x);
% </evaluate interpolation at x_eval>

% <differentiate polynomial at x_eval>
p_d = polyder(p_std);
dP_x = polyval(p_d, x_eval);
h = 1e-7;
df_x = (f_exact(x_eval + h) - f_exact(x_eval - h)) / (2 * h);
% </differentiate polynomial at x_eval>

% <integrate polynomial over [int_lo, int_hi]>
p_i = polyint(p_std);
int_P = polyval(p_i, int_hi) - polyval(p_i, int_lo);
int_f = integral(f_exact, int_lo, int_hi);
% </integrate polynomial over [int_lo, int_hi]>

% <display results>
W = 64;
sep = repmat('=', 1, W);
thin = repmat('-', 1, W);
fprintf('%s\n', sep);
fprintf('%s\n', sprintf('  Method : %s Interpolation', method_str));
fprintf('%s\n', sprintf('  Interpolation Nodes: %s', mat2str(used_x)));
fprintf('%s\n', thin);
fprintf('%s\n', sprintf('  Newton Coefficients: %s', mat2str(round(coeffs, 5))));
fprintf('%s\n', sep);
fprintf('%s\n', sprintf('  P%d(%g) = %+.8f', degree, x_eval, P_x));
fprintf('%s\n', sprintf('  f(%g)  [exact] = %+.8f', x_eval, f_x));
fprintf('%s\n', sprintf('  Absolute Error = %.4e', abs_err));
fprintf('%s\n', sprintf('  Relative Error = %.4f%%', rel_err));
fprintf('%s\n', thin);
fprintf('%s\n', sprintf('  dP%d/dx at x = %g = %+.8f', degree, x_eval, dP_x));
fprintf('%s\n', sprintf('  f''(%g)  [exact] = %+.8f', x_eval, df_x));
fprintf('%s\n', thin);
fprintf('%s\n', sprintf('  Integral P%d [%g, %g] = %+.8f', degree, int_lo, int_hi, int_P));
fprintf('%s\n', sprintf('  Integral f  [%g, %g] = %+.8f', int_lo, int_hi, int_f));
fprintf('%s\n\n', sep);
% </display results>

% <plot>
x_plot = linspace(min(x_data) - 0.1, max(x_data) + 0.1, 600);
y_true = f_exact(x_plot);
y_poly = polyval(p_std, x_plot);

figure('Name', sprintf('Newton %s', method_str), 'NumberTitle', 'off', 'Color', 'w');
hold on; grid on; box on;
plot(x_plot, y_true, 'b-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('f(x) = %s', f_expr));
plot(x_plot, y_poly, 'r--', 'LineWidth', 2, ...
    'DisplayName', sprintf('P_{%d}(x) [%s]', degree, method_str));
plot(x_data, y_data, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Data Points');
plot(x_eval, P_x, 'rv', 'MarkerSize', 11, 'MarkerFaceColor', 'r', ...
    'DisplayName', sprintf('P_{%d}(%.1f) = %.4f', degree, x_eval, P_x));
plot(x_eval, f_x, 'b^', 'MarkerSize', 11, 'MarkerFaceColor', 'b', ...
    'DisplayName', sprintf('f(%.1f) = %.4f', x_eval, f_x));
xlabel('x', 'FontSize', 13);
ylabel('y', 'FontSize', 13);
title(sprintf('Newton %s Interpolation — EBS 207', method_str), 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);
hold off;
% </plot>