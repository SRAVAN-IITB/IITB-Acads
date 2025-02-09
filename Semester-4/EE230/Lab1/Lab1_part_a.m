% Define the parameters of the RC filter
R1 = 1e3; % Resistance in ohms
C1 = 100e-9; % Capacitance in farads
fc1 = 1/(2*pi*R1*C1); % Cut-off frequency in hertz
wc1 = 2*pi*fc1; % Cut-off angular frequency in radians per second

% Define the parameters of the square wave input
A1 = 5/2; % Amplitude in volts (5 Vpp with 2.5 V DC offset)
f1 = 1/(2e-3); % Frequency in hertz (period of 2 milliseconds)
w1 = 2*pi*f1; % Angular frequency in radians per second
T1 = 1/f1; % Period in seconds
D1 = 0.5; % Duty cycle (fraction of period that is high)

% Define the time vector for simulation
t1 = 0:1e-5:3*T1; % Time in seconds from 0 to 3 cycles with 0.01 ms resolution

% Generate the square wave input
Vin1 = A1*square(w1*t1 + pi, D1*100) + 2.5; % Square wave function in MATLAB with 2.5 V DC offset

% Calculate the output of the RC filter
Vout1 = zeros(size(Vin1)); % Initialize the output vector
Vout1(1) = Vin1(1); % Initial condition
for i = 2:length(t1) % Loop over the time vector
    dt1 = t1(i) - t1(i-1); % Time step
    dV1 = dt1*(-Vout1(i-1)*wc1 + Vin1(i)*wc1); % Differential equation of RC filter
    Vout1(i) = Vout1(i-1) + dV1; % Update the output
end

% Find the index where the output reaches 63.2% of its max value
idx_632 = find(Vout1 >= 0.632*2*A1, 1, 'first');

% Plot the input waveform
figure;
subplot(2,1,1);
plot(t1, Vin1, 'b') % Plot the input in blue
xlabel('Time (s)') % Label the x-axis
ylabel('Voltage (V)') % Label the y-axis
title('Input Square Wave') % Add a title
grid on % Add a grid

% Plot both input and output waveforms
subplot(2,1,2);
plot(t1, Vin1, 'b', t1, Vout1, 'r') % Plot the input in blue and the output in red
hold on
plot(t1(idx_632), Vout1(idx_632), 'go') % Mark the point where it reaches 63.2%
hold off
xlabel('Time (s)') % Label the x-axis
ylabel('Voltage (V)') % Label the y-axis
legend('Input', 'Output', '63.2% Point') % Add a legend
title('RC Low Pass Filter Response to Square Wave Input') % Add a title
grid on % Add a grid
