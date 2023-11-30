% COMM.SYS.300 COMMUNICATION THEORY
% Matlab Exercise #1

%% 1 Sampling
% copy of instruction

Fs = 500;               % sampling frequency = 10*50Hz
Ts = 1/Fs;              % sampling interval
t = 0:Ts:0.1;           % sampling time instants [s]
x = sin(2*pi*50*t);     % signal vector
figure                  % Opens a new figure
plot(t,x)               % plot in time domain 

%% 2 Fourier transform
% copy of instruction

F_x = fft(x);                % DFT of X, saved to Fx
Nx = length(x);
Fo = 1/(Ts*Nx);              % frequency resolution.
disp(Nx)
disp(Fo);
freq1 = 0:Fo:(Nx-1)*Fo;      % One-sided frequency Axis
figure
plot(freq1,abs(F_x)/Nx)      % One-sided amplitude Spectrum

freq2 = -Nx/2*Fo:Fo:(Nx/2-1)*Fo;   % Two-sided frequency Axis
figure
plot(freq2,fftshift(abs(F_x)/Nx)); % Two-sided amplitude Spectrum


%% 3 SPECTRAL ILLUSTRATIONS

% 3.1 SIGNAL GENERATION

Fs = 16e9;            % sampling frequecy
Ts = 1/Fs;            % sampling interval
t = 0:Ts:(1000-1)*Ts; % sampling time instants [s]

% fc = 800 MHz sinusoidal‐signal 
% sampling frequency of Fs = 16 GHz (=16000MHz).

x=sin(2*pi*800e6*t);  % signal vector

figure('Name', '3-1. Plot in time domain')                
plot(t,x)             % plot in time domain 
title('Time domain Plot of x(t)')
xlabel('t [s]')
ylabel('Amplitude')
axis([0 20e-9 -1.2 1.2])

Nx=length(x);
Fo=1/(Ts*Nx);
freq2 = -Nx/2*Fo:Fo:(Nx/2-1)*Fo;
freq = freq2 ./1e6;

Fx = fft(x);

figure('Name', '3-1. Two-sided amplitude spectrum') 
plot(freq2, fftshift(abs(Fx)/Nx));
title('Two-sided amplitude spectrum of sinusoidal signal')
xlabel('frequency [Hz]')
ylabel('Amplitude [-]')
axis([-Fs/2 Fs/2 0 0.5]);

%% 3.2 MULTIPLICATION BETWEEN TWO SIGNALS

m = sin(2*pi*750e6*t);
s = x.*m;

figure('Name', '3-2. Multiplication')
plot(t,s)
title('Time domain Plot of multiplication of two signals s(t)')
xlabel('t [s]')
ylabel('Amplitude')
axis([0 60e-9 -1.2 1.2])

Fms = fft(s);
figure('Name','3-2. Multiplication') 
plot(freq2, fftshift(abs(Fms)/Nx));
title('Two-sided amplitude spectrum of multiplication of two signals')
xlabel('frequency [Hz]')
ylabel('Amplitude [-]')
axis([-Fs/2 Fs/2 0 0.25]);


%% 3.3 ADDING A NOISE SIGNAL

n = randn(size(s));
y = 10*s + n;

figure('Name', '3-3. Adding Noise') 
plot(t, y);
title('Time domain plot of sinusoidal signal with noise')
xlabel('t [s]')
ylabel('Amplitude [-]')

Fy = fft(y);
figure('Name', '3-3. Adding Noise') 
plot(freq, fftshift(abs(Fy)/Nx));
title('Two-sided amplitude spectrum of sinusoidal signal with noise')
xlabel('frequency [MHz]')
ylabel('Amplitude [-]')


%% 4 LINEAR FILTERING


% 4.1 LOW PASS BUTTERWORTH FILTER

order = 10;
f_cut = 200e6;
fr = f_cut/(Fs/2);           % Cut-off frequency normalized to 1.
[b,a] = butter(order, fr);   % Coefficients of the filter

% Plot the frequency response of the filter (help freqz)
freqz(b, a, Nx, Fs)
title('Frequency response of the Butterworth filter')

% Filter the signal y(t) with the generated Butterworth filter
y_filtered_butter = filter(b, a, y);
Fy_butter = fft(y_filtered_butter);

% Plot the filtered signal in time and frequency domain
figure('Name', '4-1. Butterworth filter')
plot(t, y)
hold on
plot(t, y_filtered_butter)
title('Butterworth, time domain')

figure('Name', '4-1. Butterworth filter')
plot(freq, fftshift(abs(Fy)/Nx))
hold on
plot(freq, fftshift(abs(Fy_butter)/Nx))
title('Butterworth, frequency domain')
xlabel('f [MHz]')
ylabel('amplitude [-]')

% low frequency part of the signal (50 Hz) is visible in Butterworth
% frequency and time domain plot (low frequency)

%% 4.2 BADN PASS FIR (Finite impulse response) FILTER

order = 60;    % filter order
f_filter = [0 0.8e9 1.3e9 1.8e9 2.3e9 Fs/2]/(Fs/2);
a_filter = [0 0 1 1 0 0];
b = firpm(order, f_filter, a_filter);

% Plot the impulse response of the filter by using stem‐function (help stem)
figure('Name','4-2. FIR filter')
stem(-order/2:order/2,b)
title('FIR, impulse response of the filter')

% Plot the amplitude response of the filter by using fft
F_b = fft(b, Nx);  % same length with the frequency vector for the fft

figure('Name', '4-2. FIR filter')
plot(freq, fftshift(abs(F_b)/Nx))
title('FIR, amplitude response of the filter')
xlabel('f [MHz]')
ylabel('amplitude [-]')

% Filter the signal y(t) (from task 3.3) with the generated FIR filter (help filter)
y_filtered_FIR = filter(b, 1, y);

% Plot the filtered signal in time and frequency domain
% Compare y(t) and the filtered y(t) with each other in time domain and
% frequency domain

Fy_fir = fft(y_filtered_FIR);

figure('Name','4-2. FIR filter')
plot(t, y)
hold on
plot(t, y_filtered_FIR)
title('FIR, time domain')

figure('Name', '4-2. FIR filter')
plot(freq, fftshift(abs(Fy)/Nx))
hold on
plot(freq, fftshift(abs(Fy_fir)/Nx))
title('FIR, frequency domain')
xlabel('f [MHz]')
ylabel('amplitude [-]')

% after band pass filter, the signal of summation of two frequency
% (800+750=1550 MHz) is visible in both time and frequency domain
