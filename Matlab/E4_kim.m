%% COMM.SYS.300 Exercise 4
% Kim Yukyeong
% 15.11.2023

%% 1 SYSTEM MODEL AND GENERATION OF QAM SYMBOLS
% 1.1 System model

clear; close all; clc;

% 1.2 System parameters
% useful functions: sqrt, bsxfun, repmat, mean, abs, randi
% YOUR CODE HERE:

% ** Try the above code by altering the SNR and the QAM constellation size (QPSK, 64-QAM, etc.)

% define alphabet size [number of symbols in QAM alphabet]
alphabet_size = 16;
% alphabet sizes are 4, 16, 64, 256, 1024, ... the possible values are given by the vector 2.^(2:2:N), for any N

% define the signal-to-noise power ratio [dB]
SNR = 20;                   

% System parameters
T = 1/10e6;                 % Symbol time interval [s]
fc = 75e6;                  % Carrier frequency
r = 20;                     % Oversampling factor (r samples per pulse)
N_symbols_per_pulse = 30;   % Duration of TX/RX-filters in numbers of symbols
alfa = 0.25;                % Roll-off factor (excess bandwidth)

% Sampling frequency and sampling time interval
Fs = r/T;                   % Sampling frequency
Ts = 1/Fs;                  % Sampling time interval

% Number of training symbols and data symbols
N_data_symbols     = 10000;  % Number of data symbols
N_training_symbols = 940;    % Number of training symbols

%% 1.3 GENERATION OF QAM SYMBOLS AND THE TRANSMITTED SYMBOL FRAME

% Generate the data symbols
%  - Create a QAM constellation (help bsxfun)
%  - Scale the constellation so that the expected average power of transmitted symbols equals to one
%  - Generate the specified number (N_data_symbols) of random data symbols
%  - Plot the transmitted symbols in complex plane (a constellation)

% qam_axis presents the symbol values in real/imaginary axis, for different alphabet/constellation sizes ():
qam_axis = -sqrt(alphabet_size)+1:2:sqrt(alphabet_size)-1;

% qam_axis = [-1 1];                % for QPSK
% qam_axis = [-3 -1 1 3];           % for 16-QAM
% qam_axis = [-7 -5 -3 -1 1 3 5 7]; % for 64-QAM

% generation of a complex constellation:
alphabet = bsxfun(@plus, qam_axis', 1j*qam_axis);
% **function bsxfun equivalent to 
% alphabet = repmat(qam_axis', 1, sqrt(alphabet_size)) + repmat(1j*qam_axis, sqrt(alphabet_size), 1);

alphabet = alphabet(:).';                             % alphabet symbols as a row vector


% Scaling the constellation, so that the mean power of a transmitted symbol is 1 
% (for QPSK ==> 1/sqrt(2), for 16-QAM ==> 1/sqrt(10)) 
alphabet_scaling_factor = 1/sqrt(mean(abs(alphabet).^2));
alphabet = alphabet*alphabet_scaling_factor;

% Random vector of symbol indices (i.e., numbers between 1...alphabet_size)
symbol_ind = randi(length(alphabet),1,N_data_symbols);
data_symbols = alphabet(symbol_ind);                   % Data symbols to be transmitted

figure ('Name', '1.3 Constellation')
plot(data_symbols,'bo')
xlabel('Re')
ylabel('Im')
title('Transmitted data symbols')


% Generate training (QAM) symbols and create the transmitted symbol frame
%  - Use the above-defined QAM constellation also for
%  - Generate the specified number (N_training_symbols) of random training symbols

% Generation of training symbols (similar to data symbols):
training_symbols = alphabet(randi(length(alphabet), 1, N_training_symbols));

% Concatenating the training and data symbols to get the overall transmitted symbols:
symbol_frame = [training_symbols data_symbols];


%% 2 TRANSMITTER STRUCTURE

% useful functions: rcosdesign, plot, legend, zeros, length, real, fft,
% fftshift
% YOUR CODE HERE:

% Implement the transit filter: Root-Raised-Cosine (RRC) and plot the pulse shape
p = rcosdesign(alfa,N_symbols_per_pulse, r,'sqrt');
figure ('Name', '2 TX')
plot(-N_symbols_per_pulse*r/2*Ts:Ts:N_symbols_per_pulse*r/2*Ts,p,'b')
hold on
plot(-N_symbols_per_pulse*r/2*Ts:T:N_symbols_per_pulse*r/2*Ts, p(1:r:end),'ro')
xlabel('time [s]')
ylabel('Amplitude')
title('Transmit/receive RRC filter (pulse shape)')
legend('Pulse shape', 'Ideal symbol-sampling locations')


% Filter the transmitted symbol frame
%  - upsample the symbol sequence rate to match with sampling rate of the filter/pulse:

symbols_upsampled = zeros(size(1:r*length(symbol_frame)));   
symbols_upsampled(1:r:r*length(symbol_frame)) = symbol_frame;


% the up-sampled sequence looks like {a1 0 0... a2 0 0... a3 0 0...}
x_LP = filter(p, 1, symbols_upsampled);           % Transmitter filtering
x_LP = x_LP(1+(length(p)-1)/2:end);               % Filter delay correction

% x_LP is the complex-valued lowpass equivalent signal of the transmitted real-valued bandpass signal x_BP

% Implement the upconversion (modulation = frequency translation to the carrier frequency fc)
%  - Define the time vector for the oscillator signal based on the reference clock time in the TX
%  - Generate the complex-exponential carrier signal using the above-defined time vector
%  - Multiply (="mix") the carrier signal with the low-pass equivalent QAM signal (x_LP)


% Time vector for the TX oscillator signal:
t_TX_oscillator = 0:Ts:Ts*(length(x_LP)-1);

% TX oscillator signal:
TX_oscillator_signal = exp(1j*2*pi*fc*t_TX_oscillator);

% Carrier modulation / upconversion (complex valued):
x_BP_complex = x_LP.*TX_oscillator_signal;

% finalize the TX process (lowpass-to-bandpass transformation) take the real part of the signal (and scale with sqrt(2))
% Taking the real value to finalize the lowpass-to-bandpass transformation:
x_BP = sqrt(2)*real(x_BP_complex);


% Plot the following figures
% - The complex-valued low-pass equivalent signal x_LP in time and frequency domain
% - The complex-valued bandpass signal x_BP_complex in frequency domain
% - The real-valued bandpass signal x_BP in time and frequency domain

figure ('Name', '2 TX')                 
plot(t_TX_oscillator, abs(x_LP))
xlabel('Time [s]')
ylabel('Amplitude (of a complex signal)')
title('Complex-valued Lowpass signal in time domain')

figure ('Name', '2 TX')                 
plot(t_TX_oscillator, x_BP)                         % notice no abs needed
xlabel('Time [s]')
ylabel('Amplitude')
title('Real-valued Bandpass signal in time domain')

NFFT = 2^14;                                        % FFT size
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts);             % frequency vector

figure ('Name', '2 TX')
plot(f/1e6, fftshift(abs(fft(x_LP,NFFT))))
xlabel('Frequency [MHz]')
ylabel('Amplitude ')
title('Amplitude spectrum of the lowpass signal')

figure ('Name', '2 TX')
plot(f/1e6, fftshift(abs(fft(x_BP_complex,NFFT)))) 
xlabel('Frequency [MHz]')
ylabel('Amplitude')
title('Amplitude spectrum of the Modulated (complex-valued) BP signal')

figure ('Name', '2 TX')
plot(f/1e6, fftshift(abs(fft(x_BP,NFFT))))          
xlabel('Frequency [MHz]')
ylabel('Amplitude')
title('Amplitude spectrum of the real-valued BP signal')



%% 3 CHANNEL MODEL
% useful functions: randnm varm 
% YOUR CODE HERE:

% simple AWGN channel model.
%  - create white random noise,
%  - scale it with the proper scaling factor to obtain the desired SNR, 
%  - add it on top of the transmitted signal (x_BP).

% Generate the noise
n = randn(size(x_BP));              % White Gaussian random noise
P_x_BP = var(x_BP);                 % Signal power
P_n = var(n);                       % Noise power

% Defining noise scaling factor based on the desired SNR:
noise_scaling_factor = sqrt(P_x_BP/P_n/10^(SNR/10)*(r/(1+alfa)));
% disp("Noise scaling factor = " + noise_scaling_factor);


% Q. Why r/(1+alfa) part?
% --> Because desired(given) SNR value is based on the signal-to-noise where there is no oversampling,
% but our signal is oversampled by factor of fs/W, so the data sample power is reduced about oversampling rate,
% SNR = Es(no oversampling)/N0,   W = (1+alpa)/T
% Es_oversampled/N0 = (Es/N0)/(fs/W)
% Es_oversampled/N0 = SNR /(fs/W) = SNR * W/fs = SNR * (1+alpa)/T * (T/r)


% Noisy signal
y_BP = x_BP + noise_scaling_factor*n;


% Plot the amplitude spectrum of the noisy bandpass signal
figure ('Name', '3 Channel')
hold on
plot(f/1e6, fftshift(abs(fft(x_BP, NFFT))), 'LineWidth', 1) 
plot(f/1e6, fftshift(abs(fft(y_BP,NFFT))), 'Color', [1, 0, 0, 0.15])          
xlabel('Frequency [MHz]')
ylabel('Amplitude')
legend('Noiseless', 'Noisy signal')
title('Amplitude spectrum of the bandpass signal')
hold off

% ** add propagation delay in part 5 before adding noise ==> implemented in Section 5


%% 4 RECEIVER STRUCTURE

% 4.1 Downconversion and channel filtering

% useful functions: length, exp, filter
% YOUR CODE HERE:

% Implement the downconversion (i.e. demodulation back to the baseband) 
% by using the same principle as in the upconversion in TX, remember to use the RX reference clock time:

% Time vector for the RX oscillator signal:
t_RX_oscillator = 0:Ts:Ts*(length(y_BP)-1);

% RX oscillator signal (notice the minus-sign compared to TX oscillator):
RX_oscillator_signal = exp(-1j*2*pi*fc*t_RX_oscillator);

% Carrier demodulation / downconversion (signal becomes complex again)
y_BP_downconverted = y_BP.*RX_oscillator_signal;


% Plot the amplitude spectrum of the downconverted signal
figure ('Name', '4.1 RX')
plot(f/1e6, fftshift(abs(fft(y_BP_downconverted, NFFT))))          
xlabel('Frequency [MHz]')
ylabel('Amplitude')
title('Amplitude spectrum of the demodulated signal')

% Filter the received signal with the receive filter (RRC)
X_LP_received = sqrt(2)*filter(p, 1, y_BP_downconverted);   % Receiver filtering
X_LP_received = X_LP_received(1+(length(p)-1)/2:end);       % Filter delay correction


% Plot the amplitude spectrum of the downconverted and filtered signal
figure ('Name', '4.1 RX')
plot(f/1e6, fftshift(abs(fft(X_LP_received, NFFT))))          
xlabel('Frequency [MHz]')
ylabel('Amplitude')
title('Amplitude spectrum of demodulated & filtered signal')


% t_RX_oscillator = ... 
% **add clock timing offset in part 5 instead of perfect synchronization

%% 4.2 Signal sampling
% useful functions: 
% YOUR CODE HERE:

% Sampling the received signal in order to get symbol samples
RX_symbol_frame = X_LP_received(1:r:end);

% Take user data symbols and training symbols in separate vectors:
RX_training_symbols = RX_symbol_frame(1:N_training_symbols);
RX_data_symbols = RX_symbol_frame(N_training_symbols+1:end);


% Plot the received symbol samples in complex plane (constellation)
figure ('Name', '4.2 Signal sampling')
plot(RX_data_symbols,'bo')
hold on
plot(alphabet,'rs')
hold off
xlabel('Re')
ylabel('Im')
title('Received data symbols of perfect synchronization')


%% 4.3 Obtaining symbol decisions
% useful functions: abs, min, mean
% YOUR CODE HERE:

% Calculate the Euclidian distance between each symbol sample and each alphabet symbol:
alphabet_error_matrix = abs(bsxfun(@minus, alphabet.', RX_data_symbols));

% rows represent the alphabet symbol indices and columns represent the received symbol indices 
% Euclidian distance between the 5th received symbol and the 3rd symbol in the alphabet is given as
%   alphabet_error_matrix(3,5)

% Find the indices of the alphabet symbols, which have the minimum distance to the symbol samples:
[~, estimated_symbol_ind] = min(alphabet_error_matrix);


% Finding out which symbols were estimated incorrecly:
symbol_errors = ...
 estimated_symbol_ind ~= symbol_ind(1:length(estimated_symbol_ind));


% due to filtering transitions, we lose some of the last symbols in the symbol frame. 
% In practice, continue taking a few samples after the frame to try to get all the symbols. 
% However, filter transitions in the beginning and end of the frame are always creating non-idealities to the transmission 
% (the same is also happening in frequency domain: compare data in the middle and in the edge of the used band).

% Symbol Error Rate (SER) (0 means 0% of errors, 1 means 100% of errors)
SER = mean(symbol_errors);

disp('Case of perfect synchronization')
fprintf('SNR = %d, alphabet size = %d --> Symbol error rate: %.7f\n\n', SNR, alphabet_size, SER)

% YOUR ANSWERS HERE:
% Q. What is the worst value for SER? Why?
% --> worst value of BER (bit error rate) is 0.5,
% 100% bit error can be easily corrected to 0%, 0.6% BER is same as 0.4%  
% --> relation between BER and SER is ==> BER = M/(2M-2) * SER  
% --> SER can be considered as worst when BER is 0.5, 

% --> SER = BER * (2M-2)/M = 0.5 * (2M-2)/M = (M-1)/M
% --> for BPSK, worst SER == 1/2   = 0.5 
% --> 4-QAM,  worst SER == 3/4     = 0.75  
% --> 16-QAM, worst SER == 15/16   = 0.9375
% --> 64-QAM, worst SER == 63/64   = 0.9844


% Try the above code by altering the SNR and the QAM constellation size (QPSK, 64-QAM, etc.)
% higher SNR -> lower SER, smaller alphabet size -> lower SER

% SNR = 10 db
% --> QPSK   : 0.012136 ( 1 % error)
% --> 16-QAM : 0.348947 (35 %)
% --> 64-QAM : 0.763189 (76 %)

% SNR = 20 db
% --> QPSK   : 0.000000 (no error)
% --> 16-QAM : 0.000401 (0.05 %)
% --> 64-QAM : 0.143330 (  14 %)

% SNR = 100 db
% --> QPSK   : 0.000000 (no error)
% --> 16-QAM : 0.000000 
% --> 64-QAM : 0.000000 

%% 5. TIMING SYNCHRONIZATION AND PHASE CORRECTION
% useful functions: randi
% YOUR CODE HERE:

% Add an unknown propagation delay to the AWGN channel model. 
% In contrast to, the RX filter delay which is known, 
% the channel delay is random and must be estimated in the RX. 
% This is usually referred as timing synchronization. 

% NOTE: This should be added to the channel part before the noise is added to X_BP (see Sect. 3).

% Lastly, to consider the fact that the oscillator clocks are not synchronized, 
% we define separate reference clock times for the TX and RX. 
% define the TX and RX clocks having random offsets uniformly distributed between 0…1s.

% Time vector for the TX and RX oscillator signal:
TX_clock_start_time = rand;             % Clock start time in the TX oscillator
RX_clock_start_time = rand;             % Clock start time in the RX oscillator

% TX oscillator time vectors are redefined to use these offsets:
t_TX_oscillator = TX_clock_start_time + (0:Ts:Ts*(length(x_LP)-1));

% TX oscillator signal:
TX_oscillator_signal = exp(1j*2*pi*fc*t_TX_oscillator);

% Carrier modulation / upconversion (still complex valued):
x_BP_complex = x_LP.*TX_oscillator_signal;

% Taking the real value to finalize the lowpass-to-bandpass transformation:
x_BP = sqrt(2)*real(x_BP_complex);

% unknown propagation delay
p_delay = randi(940);                  % Delay in samples
x_BP = [zeros(1, p_delay), x_BP];      % Add delay zeros

% addding AWGN noise
n = randn(size(x_BP));              % White Gaussian random noise
P_x_BP = var(x_BP);                 % Signal power
P_n = var(n);                       % Noise power

% Defining noise scaling factor based on the desired SNR:
noise_scaling_factor = sqrt(P_x_BP/P_n/10^(SNR/10)*(r/(1+alfa)));

% Delayed noisy signal
y_BP = x_BP + noise_scaling_factor*n;

% Plot the amplitude spectrum of the bandpass signal
NFFT = 2^14;                                        
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts);             
figure ('Name', '5 Propagation Delay')
hold on
plot(f/1e6, fftshift(abs(fft(x_BP, NFFT)))) 
plot(f/1e6, fftshift(abs(fft(y_BP,NFFT))), 'Color', [1, 0, 0, 0.1])          
xlabel('Frequency [MHz]')
ylabel('Amplitude')
title('Amplitude spectrum of the propagation delayed BP signal')
legend('noiseless', 'noisy')
hold off

% RX oscillator time vectors are redefined to use these offsets:
t_RX_oscillator = RX_clock_start_time + (0:Ts:Ts*(length(y_BP)-1));

% RX oscillator signal (notice the minus-sign compared to TX oscillator):
RX_oscillator_signal = exp(-1j*2*pi*fc*t_RX_oscillator);

% Carrier demodulation / downconversion (signal becomes complex again)
y_BP_downconverted = y_BP.*RX_oscillator_signal;

% Filter the received signal with the receive filter (RRC similar to TX)
X_LP_received = sqrt(2)*filter(p, 1, y_BP_downconverted);  % Receiver filtering
X_LP_received = X_LP_received(1+(length(p)-1)/2:end);      % Filter delay correction

RX_symbol_frame = X_LP_received(1:r:end);
RX_training_symbols = RX_symbol_frame(1:N_training_symbols);
RX_data_symbols = RX_symbol_frame(N_training_symbols+1:end);

alphabet_error_matrix = abs(bsxfun(@minus, alphabet.', RX_data_symbols));

[~, estimated_symbol_ind] = min(alphabet_error_matrix);

symbol_errors = ...
 estimated_symbol_ind ~= symbol_ind(1:length(estimated_symbol_ind(1:length(N_data_symbols))));

SER = mean(symbol_errors);

disp('Case of random channel delay')
fprintf('SNR = %d, alphabet size = %d --> Symbol error rate: %.7f\n\n', SNR, alphabet_size, SER)

%% 5.1 Timing synchronization
% useful functions: xcorr, zeros, max
% YOUR CODE HERE:

% To perform the correlation between the known training symbols and the received QAM signal, 
% should upsample the training symbols in order match the sampling rates 

training_signal = zeros(size(1:r*length(training_symbols))); 
% Zero vector initilized for Up-sampled symbol sequence

training_signal(1:r:r*length(training_symbols)) = training_symbols;
% the up-sampled sequence looks like {a1 0 0 a2 0 0 a3 0 0 a4 0 ...}

% Calculate the cross-correlation as function time-delay between received signal and upsampled training symbols 
[corr_fun, delay] = xcorr(X_LP_received, training_signal);

% Plot the amplitude of cross-correlation function
figure ('Name', '5.1 Timing synchronization')
plot(delay,abs(corr_fun))
xlabel('Delay [samples]')
ylabel('Correlation')
title('Cross-correlation between transmitted and received training symbols')

% Find the correlation peak to obtain the timing delay (help max):
% Find the sample index with the maximum correlation value:
[~, max_ind] = max(abs(corr_fun));

% Estimated delay: The 1st signal sample should be taken at "timing_ind+1"
timing_ind = delay(max_ind);


% Now, the sampling can be performed as
RX_symbol_frame = X_LP_received(timing_ind+1:r:end);

% where after the training symbols and data symbols are separated, as above.
RX_training_symbols = RX_symbol_frame(1:N_training_symbols);
RX_data_symbols = RX_symbol_frame(N_training_symbols+1:end);


%% 5.2 Phase error estimation and phase compensation
% useful functions: 
% YOUR CODE HERE:

% Plot the received symbol samples in complex plane (constellation) 
% before phase error compensation to observe the effects of the error.

figure ('Name', '5.2 Phase Compensation')
plot(RX_data_symbols,'bo')
hold on
plot(alphabet,'rs')
hold off
xlabel('Re')
ylabel('Im')
title('Received data symbols with phase error')

% Without going to details here, we use the following expression to calculate the so called channel estimate; 
% Here RX_training_symbols and training_symbols should be row-vectors

channel_estimate = ...                   
  RX_training_symbols*training_symbols'/norm(training_symbols)^2;       

% the phase error can be compensated by simply dividing the signal with the channel estimate
RX_data_symbols = RX_data_symbols/channel_estimate;


% plot the received symbol samples in complex plane (constellation) to observe the compensation result.
figure ('Name', '5.2 Phase Compensation')
plot(RX_data_symbols,'bo')
hold on
plot(alphabet,'rs')
hold off
xlabel('Re')
ylabel('Im')
title('Received data symbols after phase compensation')


%% 5.3 Symbol decisions and calculating the symbol error rate (SER)
% useful functions: 
% YOUR CODE HERE:

alphabet_error_matrix = abs(bsxfun(@minus, alphabet.', RX_data_symbols));

% Find the indices of the alphabet symbols, which have the minimum distance to the symbol samples:
[~, estimated_symbol_ind] = min(alphabet_error_matrix);

% Finding out which symbols were estimated incorrecly:
symbol_errors = ...
 estimated_symbol_ind ~= symbol_ind(1:length(estimated_symbol_ind));

% Symbol Error Rate (SER) (0 == no errors, 1 == 100% of errors)
SER = mean(symbol_errors);

disp('Case of the estimated delay')
fprintf('SNR = %d, alphabet size = %d --> Symbol error rate: %.7f\n', SNR, alphabet_size, SER)


% Q. Compare the realization of random channel delay to the estimated delay. 
% What can you say about the accuracy?
% --> Before time synchronization, the accuracy (based on the SER value) is very low 
% even in case of the very high SNR (~100db) and small size of alphabet (4),
% SER value is close to the worst one
% --> After time synchronization, accuracy of estimated delay is almost
% same as where there is no random delay

% Q. How the phase error effects the received constellation? How well is this error compensated?
% --> Received signal has certain degree of phase error which can cause symbol error, and shown as
% rotated by the phase error(certain angle for all symbols) in constellation 
% --> phase error is compensated very well by channel estimation using training symbol
