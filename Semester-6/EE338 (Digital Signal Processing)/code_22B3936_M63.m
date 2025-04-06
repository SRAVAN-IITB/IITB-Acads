%% Specifications and Initial Setup
clear;
close all;
clc;

% Sampling Parameters
Fs = 630e3; % Sampling frequency (Hz)

% Band Specifications (Hz)
f_pass1 = [65e3 95e3]; % First passband
f_pass2 = [210e3 240e3]; % Second passband
f_trans = 5e3; % Transition bandwidth
delta_1 = 0.15;
delta_2 = 0.15;

% Calculate stopband frequencies
f_stop1 = [f_pass1(1)-f_trans f_pass1(2)+f_trans];
f_stop2 = [f_pass2(1)-f_trans f_pass2(2)+f_trans];
fprintf('f_stop1 = %g, f_stop2 = %g\n', f_stop1/1000, f_stop2/1000);

omega_p1(1) = 2 * pi * f_pass1(1) / Fs; % Normalized passband edge (Group I)
omega_p1(2) = 2 * pi * f_pass1(2) / Fs; % Normalized passband edge (Group I)
omega_s1(1) = 2 * pi * f_stop1(1) / Fs; % Normalized stopband edge (Group I)
omega_s1(2) = 2 * pi * f_stop1(2) / Fs; % Normalized stopband edge (Group I)
fprintf('omega_p1 = %g to %g, omega_s1 = %g to %g \n', omega_p1(1), omega_p1(2), omega_s1(1), omega_s1(2));

omega_p2(1) = 2 * pi * f_pass2(1) / Fs; % Normalized passband edge (Group II)
omega_p2(2) = 2 * pi * f_pass2(2) / Fs; % Normalized passband edge (Group II)
omega_s2(1) = 2 * pi * f_stop2(1) / Fs; % Normalized stopband edge (Group II)
omega_s2(2) = 2 * pi * f_stop2(2) / Fs; % Normalized stopband edge (Group II)
fprintf('omega_p2 = %g to %g, omega_s2 = %g to %g \n', omega_p2(1), omega_p2(2), omega_s2(1), omega_s2(2));

% Calculating D1 and D2 to find out N
D1 = (1 / (1 - delta_1)^2) - 1;
D2 = (1 / delta_2^2) - 1;
fprintf('D1 = %g and D2 = %g\n', D1, D2);

% Computing Capital Omega
Omega_p1(1) = tan(omega_p1(1)/2);
Omega_p1(2) = tan(omega_p1(2)/2);
Omega_s1(1) = tan(omega_s1(1)/2);
Omega_s1(2) = tan(omega_s1(2)/2);
B_1 = Omega_p1(2)-Omega_p1(1);
Omega0_1 = sqrt(Omega_p1(2) * Omega_p1(1));

Omega_p2(1) = tan(omega_p2(1)/2);
Omega_p2(2) = tan(omega_p2(2)/2);
Omega_s2(1) = tan(omega_s2(1)/2);
Omega_s2(2) = tan(omega_s2(2)/2);
B_2 = Omega_p2(2)-Omega_p2(1);
Omega0_2 = sqrt(Omega_p2(2) * Omega_p2(1));

% Add a check to ensure the frequencies are valid
if any(omega_p1 > pi) || any(omega_s1 > pi) || any(omega_p2 > pi) || any(omega_s2 > pi)
    error('Normalized frequencies exceed Nyquist range (0 to pi).');
end

% Transforming to Low-pass
LP_Omega_p1(1) = -1;
LP_Omega_p1(2) = 1;
LP_Omega_s1(1) = ((Omega_s1(1))^2 - (Omega0_1)^2)/(B_1 * Omega_s1(1));
LP_Omega_s1(2) = ((Omega_s1(2))^2 - (Omega0_1)^2)/(B_1 * Omega_s1(2));

LP_Omega_p2(1) = -1;
LP_Omega_p2(2) = 1;
LP_Omega_s2(1) = ((Omega_s2(1))^2 - (Omega0_2)^2)/(B_2 * Omega_s2(1));
LP_Omega_s2(2) = ((Omega_s2(2))^2 - (Omega0_2)^2)/(B_2 * Omega_s2(2));

% For calculating N and Omega_c, we consider the min of either of the 
% LP_Omega_s(1) and LP_Omega_s(2) so as to enforce the more strigent
% condition
LP_Omega_s1_min = min(abs(LP_Omega_s1(1)), abs(LP_Omega_s1(2)));
LP_Omega_s2_min = min(abs(LP_Omega_s2(1)), abs(LP_Omega_s2(2)));
fprintf('LP_Omega_s1_min = %g and LP_Omega_s2_min = %g\n', LP_Omega_s1_min, LP_Omega_s2_min);
% Calculate order and LP_OmegaC for Group I
N1 = ceil(log10(D2 / D1) / (2 * log10(abs(LP_Omega_s1_min))));
LP_OMegaC_1_min = 1/((D1)^(0.5/N1));
LP_OMegaC_1_max = LP_Omega_s1_min/((D2)^(0.5/N1));
LP_OmegaC_1 = (LP_OMegaC_1_min + LP_OMegaC_1_max)/2;

% Calculate order and LP_OmegaC for Group II
N2 = ceil(log10(D2 / D1) / (2 * log10(abs(LP_Omega_s2_min))));
LP_OmegaC_2 = (1/(D1)^(0.5/N2) + LP_Omega_s2_min/((D2)^(0.5/N2)))/2;

fprintf('N1 = %d and N2 = %d\n', N1, N2);
fprintf('LP_OmegaC_1 = %g and LP_OmegaC_2 = %g\n', LP_OmegaC_1, LP_OmegaC_2);

% Use the higher order to ensure both bands meet specifications
N = max(N1, N2);
fprintf('Max required filter order: %g\n', N);

%% Print Design Specifications
fprintf('\n=== Filter Design Specifications ===\n');
fprintf('Sampling Frequency: %.2f kHz\n', Fs/1000);
fprintf('\nFirst Passband: %.2f - %.2f kHz\n', f_pass1(1)/1000, f_pass1(2)/1000);
fprintf('First Stopband: %.2f - %.2f kHz\n', f_stop1(1)/1000, f_stop1(2)/1000);
fprintf('\nSecond Passband: %.2f - %.2f kHz\n', f_pass2(1)/1000, f_pass2(2)/1000);
fprintf('Second Stopband: %.2f - %.2f kHz\n', f_stop2(1)/1000, f_stop2(2)/1000);
fprintf('\nTransition Width: %.2f kHz\n', f_trans/1000);
fprintf('Filter Order: %d\n', N);

%% Print Transfer Function Details
fprintf('\n=== Transfer Function Information ===\n');
fprintf('First Band:\n');
fprintf('Center Frequency (W01): %.4f\n', LP_OmegaC_1);
fprintf('Bandwidth (B1): %.4f\n', B_1);
fprintf('\nSecond Band:\n');
fprintf('Center Frequency (W02): %.4f\n', LP_OmegaC_2);
fprintf('Bandwidth (B2): %.4f\n', B_2);

%% Print Filter Equations
fprintf('\n=== Filter Equations ===\n');
fprintf('Analog Prototype: H(s) = 1/');
fprintf('(s^%d + 1)\n', N);
fprintf('\nFrequency Transformation:\n');
fprintf('s -> (s^2 + W0^2)/(Bs)\n');
fprintf('\nBilinear Transform:\n');
fprintf('s -> (z-1)/(z+1)\n');

%% Design First Bandpass Filter
% Calculate Butterworth poles for prototype LPF
fprintf('Calculated poles for the prototype LPF 1:\n');
poles1 = zeros(2*N1, 1);
for k = 1:2*N1
    poles1(k) = LP_OmegaC_1 * exp(1i*pi*(2*k + N1 - 1)/(2*N1));
    fprintf('Pole (%d) = %.4f + %.4fi\n', k, real(poles1(k)), imag(poles1(k)));

end

fprintf('\n');

% LHCP Poles
poles_LHCP_1 = [];
for k = 1:length(poles1)
    if real(poles1(k)) < 0
        poles_LHCP_1 = [poles_LHCP_1, poles1(k)];
        fprintf('Selected pole(%d) = %.4f + %.4fi\n', k, real(poles1(k)), imag(poles1(k)));
    end
end

% Plot the LHCP poles on the complex plane
figure;
plot(real(poles_LHCP_1), imag(poles_LHCP_1), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0, 0.5, 0]);
title('LHCP Poles of the Prototype Low-Pass Filter 1');
xlabel('Real');
ylabel('Imaginary');
grid on;
axis equal;

% Get transfer function of analog LPF
[num_lpf1, den_lpf1] = zp2tf([], poles_LHCP_1, LP_OmegaC_1^N1);

% Normalize coefficients
num_lpf1 = num_lpf1 / den_lpf1(1);
den_lpf1 = den_lpf1 / den_lpf1(1);

% Convert to bandpass using frequency transformation
syms s z;
analog_lpf1(s) = poly2sym(num_lpf1,s)/poly2sym(den_lpf1,s);
analog_bpf1(s) = analog_lpf1(((s*s) + (Omega0_1*Omega0_1))/((B_1)*s));

% Apply Bilinear Transform
digital_bpf1(z) = analog_bpf1((1-z)/(1+z));
fprintf('BPF_1 Digital Transfer Function:\n');
disp(digital_bpf1(z));

% Print coefficients of BPF_1
[num_bpf1_sym, den_bpf1_sym] = numden(digital_bpf1(z));
num_bpf1 = sym2poly(expand(num_bpf1_sym));
den_bpf1 = sym2poly(expand(den_bpf1_sym));

% Normalize coefficients
num_bpf1 = num_bpf1/den_bpf1(1);
den_bpf1 = den_bpf1/den_bpf1(1);

fprintf('\n=== First Filter Coefficients After Bilinear Transformation ===\n');
fprintf('Numerator coefficients (b1):\n');
for i = 1:length(num_bpf1)
    fprintf('b1(%d) = %.10g\n', i-1, num_bpf1(i));
end
fprintf('\nDenominator coefficients (a1):\n');
for i = 1:length(den_bpf1)
    fprintf('a1(%d) = %.10g\n', i-1, den_bpf1(i));
end

%% Design Second Bandpass Filter
% Calculate Butterworth poles for prototype LPF
fprintf('Calculated poles for the prototype LPF 2:\n');
poles2 = zeros(2*N2, 1);
for k = 1:2*N2
    poles2(k) = LP_OmegaC_2 * exp(1i*pi*(2*k + N2 - 1)/(2*N2));
    fprintf('Pole (%d) = %.4f + %.4fi\n', k, real(poles2(k)), imag(poles2(k)));

end

fprintf('\n');

% LHCP Poles
poles_LHCP_2 = [];
for k = 1:length(poles2)
    if real(poles2(k)) < 0
        poles_LHCP_2 = [poles_LHCP_2, poles2(k)];
        fprintf('Selected pole(%d) = %.4f + %.4fi\n', k, real(poles2(k)), imag(poles2(k)));
    end
end

% Plot the LHCP poles on the complex plane
figure;
plot(real(poles_LHCP_2), imag(poles_LHCP_2), 'x', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0, 0.5, 0]);
title('LHCP Poles of the Prototype Low-Pass Filter 2');
xlabel('Real');
ylabel('Imaginary');
grid on;
axis equal;

% Get transfer function of analog LPF
[num_lpf2, den_lpf2] = zp2tf([], poles_LHCP_2, LP_OmegaC_2^N2);

% Convert to bandpass using frequency transformation
syms s z;
analog_lpf2(s) = poly2sym(num_lpf2,s)/poly2sym(den_lpf2,s);
analog_bpf2(s) = analog_lpf2(((s*s) + (Omega0_2*Omega0_2))/((B_2)*s));

% Apply Bilinear Transform
digital_bpf2(z) = analog_bpf2((1-z)/(1+z));
fprintf('BPF_2 Digital Transfer Function:\n');
disp(digital_bpf2(z));

% Print coefficients of BPF_2
[num_bpf2_sym, den_bpf2_sym] = numden(digital_bpf2(z));
num_bpf2 = sym2poly(expand(num_bpf2_sym));
den_bpf2 = sym2poly(expand(den_bpf2_sym));

% Normalize coefficients
num_bpf2 = num_bpf2/den_bpf2(1);
den_bpf2 = den_bpf2/den_bpf2(1);

fprintf('\n=== Second Filter Coefficients After Bilinear Transformation ===\n');
fprintf('Numerator coefficients (b2):\n');
for i = 1:length(num_bpf2)
    fprintf('b2(%d) = %.10g\n', i-1, num_bpf2(i));
end
fprintf('\nDenominator coefficients (a2):\n');
for i = 1:length(den_bpf2)
    fprintf('a2(%d) = %.10g\n', i-1, den_bpf2(i));
end

%% Combine Discrete Filters
num_discrete_total = conv(num_bpf1, den_bpf2) + conv(num_bpf2, den_bpf1);
den_discrete_total = conv(den_bpf1, den_bpf2);

% Normalize combined filter coefficients
num_discrete_total = num_discrete_total/den_discrete_total(1);
den_discrete_total = den_discrete_total/den_discrete_total(1);

%% Print Coefficients After Bilinear Transformation
fprintf('\n=== Combined Filter Coefficients (NORMALIZED) After Bilinear Transformation ===\n');
fprintf('Numerator coefficients (b_total):\n');
for i = 1:length(num_discrete_total)
    fprintf('b_total(%d) = %.10g\n', i-1, num_discrete_total(i));
end
fprintf('\nDenominator coefficients (a_total):\n');
for i = 1:length(den_discrete_total)
    fprintf('a_total(%d) = %.10g\n', i-1, den_discrete_total(i));
end

% Print Transfer Function After Bilinear Transformation
fprintf('\n=== Transfer Function After Bilinear Transformation ===\n');
fprintf('H(z) = ');
for i = 1:length(num_discrete_total)
    fprintf('%.10g * z^%d ', num_discrete_total(i), i-1);
    if i < length(num_discrete_total)
        fprintf('+ ');
    end
end
fprintf('/\n');
for i = 1:length(den_discrete_total)
    fprintf('%.10g * z^%d ', den_discrete_total(i), i-1);
    if i < length(den_discrete_total)
        fprintf('+ ');
    end
end
fprintf('\n');
%% All Plots
fs = 630e3;  % Sampling frequency in Hz
f = linspace(0, fs/2, 1024);  % Frequency axis from 0 to Nyquist

[H1, w1] = freqz(num_bpf1, den_bpf1, 1024, fs); % Frequency response 
% Plot Magnitude Response
figure; hold on; grid on;
plot(w1, abs(H1), 'b', 'LineWidth', 2);  % Plot frequency response in blue
xlabel('Frequency (KHz)');
ylabel('Magnitude');
yline(0.15, '--k', 'Tolerance 0.15', 'LabelHorizontalAlignment', 'left'); 
yline(0.85, '--k', 'Tolerance 0.85', 'LabelHorizontalAlignment', 'left');

% Find intersection points for y=0.15
crossing_015_indices_1 = find(diff(sign(abs(H1) - 0.15))); % Find indices where the magnitude crosses 0.15
intersect_015_1 = w1(crossing_015_indices_1); % Frequencies at these indices

% Find intersection points for y=0.85
crossing_085_indices_1 = find(diff(sign(abs(H1) - 0.85))); % Find indices where the magnitude crosses 0.85
intersect_085_1 = w1(crossing_085_indices_1); % Frequencies at these indices

% Mark the intersection points
plot(intersect_015_1, 0.15*ones(size(intersect_015_1)), 'ro', 'MarkerFaceColor', 'r'); % Red markers for y=0.15 intersections
plot(intersect_085_1, 0.85*ones(size(intersect_085_1)), 'go', 'MarkerFaceColor', 'g'); % Green markers for y=0.85 intersections

% Label the intersection points for y = 0.15
for i = 1:length(intersect_015_1)
    text(intersect_015_1(i), 0.15, sprintf('%.1f kHz', intersect_015_1(i)/1e3), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Label the intersection points for y = 0.85
for i = 1:length(intersect_085_1)
    text(intersect_085_1(i), 0.85, sprintf('%.1f kHz', intersect_085_1(i)/1e3), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

% Labels and Legends
xlabel('Frequency (Hz)');
ylabel('Magnitude Response');
title('Magnitude Response of First Passband Butterworth Filter');
legend({'Magnitude Response', 'Tolerance of Passband', 'Tolerance of Passband'}, 'Location', 'Best');
hold off;
 

[H2, w2] = freqz(num_bpf2, den_bpf2, 1024, fs); % Frequency response
% Plot Magnitude Response
figure; hold on; grid on;
plot(w2, abs(H2), 'b', 'LineWidth', 2);  % Plot frequency response in blue
xlabel('Frequency (KHz)');
ylabel('Magnitude');
title('Magnitude Response of First Passband Butterworth Filter');
yline(0.15, '--k', 'Tolerance 0.15', 'LabelHorizontalAlignment', 'left'); 
yline(0.85, '--k', 'Tolerance 0.85', 'LabelHorizontalAlignment', 'left');

% Find intersection points for y=0.15
crossing_015_indices_2 = find(diff(sign(abs(H2) - 0.15))); % Find indices where the magnitude crosses 0.15
intersect_015_2 = w2(crossing_015_indices_2); % Frequencies at these indices

% Find intersection points for y=0.85
crossing_085_indices_2 = find(diff(sign(abs(H2) - 0.85))); % Find indices where the magnitude crosses 0.85
intersect_085_2 = w2(crossing_085_indices_2); % Frequencies at these indices

% Mark the intersection points
plot(intersect_015_2, 0.15*ones(size(intersect_015_2)), 'ro', 'MarkerFaceColor', 'r'); % Red markers for y=0.15 intersections
plot(intersect_085_2, 0.85*ones(size(intersect_085_2)), 'go', 'MarkerFaceColor', 'g'); % Green markers for y=0.85 intersections

% Label the intersection points for y = 0.15
for i = 1:length(intersect_015_2)
    text(intersect_015_2(i), 0.15, sprintf('%.1f kHz', intersect_015_2(i)/1e3), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Label the intersection points for y = 0.85
for i = 1:length(intersect_085_2)
    text(intersect_085_2(i), 0.85, sprintf('%.1f kHz', intersect_085_2(i)/1e3), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

% Labels and Legends
xlabel('Frequency (Hz)');
ylabel('Magnitude Response');
title('Magnitude Response of Second Passband Butterworth Filter');
legend({'Magnitude Response', 'Tolerance of Passband', 'Tolerance of Passband'}, 'Location', 'Best');
hold off;
 

[H, w] = freqz(num_discrete_total, den_discrete_total, 1024, fs); % Frequency response
% Plot Magnitude Response
figure; hold on; grid on;
plot(w, abs(H), 'b', 'LineWidth', 2);  % Plot frequency response in blue
xlabel('Frequency (KHz)');
ylabel('Magnitude');
title('Magnitude Response of Second Passband Butterworth Filter');
yline(0.15, '--k', 'Tolerance 0.15', 'LabelHorizontalAlignment', 'left'); 
yline(0.85, '--k', 'Tolerance 0.85', 'LabelHorizontalAlignment', 'left');

% Shading for Passband Regions
fill([65e3 95e3 95e3 65e3], [0.85 0.85 1 1], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([210e3 240e3 240e3 210e3], [0.85 0.85 1 1], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Shading for Stopband Regions
fill([0 60e3 60e3 0], [0 0 0.15 0.15], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([100e3 205e3 205e3 100e3], [0 0 0.15 0.15], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([245e3 fs/2 fs/2 245e3], [0 0 0.15 0.15], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Shading for Transition Bands
fill([60e3 65e3 65e3 60e3], [0 0 1 1], 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([95e3 100e3 100e3 95e3], [0 0 1 1], 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([205e3 210e3 210e3 205e3], [0 0 1 1], 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([240e3 245e3 245e3 240e3], [0 0 1 1], 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Find intersection points for y=0.15
crossing_015_indices = find(diff(sign(abs(H) - 0.15))); % Find indices where the magnitude crosses 0.15
intersect_015 = w(crossing_015_indices); % Frequencies at these indices

% Find intersection points for y=0.85
crossing_085_indices = find(diff(sign(abs(H) - 0.85))); % Find indices where the magnitude crosses 0.85
intersect_085 = w(crossing_085_indices); % Frequencies at these indices

% Mark the intersection points
plot(intersect_015, 0.15*ones(size(intersect_015)), 'ro', 'MarkerFaceColor', 'r'); % Red markers for y=0.15 intersections
plot(intersect_085, 0.85*ones(size(intersect_085)), 'go', 'MarkerFaceColor', 'g'); % Green markers for y=0.85 intersections

% Label the intersection points for y = 0.15
for i = 1:length(intersect_015)
    text(intersect_015(i), 0.15, sprintf('%.1f kHz', intersect_015(i)/1e3), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Label the intersection points for y = 0.85
for i = 1:length(intersect_085)
    text(intersect_085(i), 0.85, sprintf('%.1f kHz', intersect_085(i)/1e3), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
end

% Labels and Legends
xlabel('Frequency (Hz)');
ylabel('Magnitude Response');
title('Magnitude Response of Multi-Passband Butterworth Filter');
legend({'Magnitude Response', 'Tolerance of Passband', 'Tolerance of Passband', 'Passband', 'Transition Band', 'Stopband'}, 'Location', 'Best');
hold off;



%% Plot dB-Response in a Separate Window Using fvtool
fvtool(num_discrete_total, den_discrete_total, 'Fs', Fs);