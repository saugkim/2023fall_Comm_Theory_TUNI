% COMM.SYS.300 COMMUNICATION THEORY
% Matlab Exercise #2
% Kim Yukyeong

% Random variables and random processes

%% 1. RANDOM VARIABLES

% 1.1 ROLLING A FAIR 6‐FACED DICE (DISCRETE VALIABLE)

% Generate random samples for rolling a fair 6‐faced dice for N times 

clc; close all; clear;

N_faces = 6;                            % Number of faces in the dice
N_trials = 10;    %100000               % Number of trials
trials = randi(N_faces, 1, N_trials);   % Getting random integers 1...6

% Plot the histogram of the outcome values (help histogram, help hist, help bar, help histcounts)

x_histogram_centers = 1:6;              % x-axis for the bin center coordinates
histogram_edges = 0.5:1:6.5;            % x-axis for the bin edge coordinates
[hist_count, ~]= histcounts(trials, histogram_edges);

figure('Name', '1.1 Rolling histogram')
bar(x_histogram_centers, hist_count)
title('Dice rolling histogram')
xlabel('Rolled number')
ylabel('occurrence number')
grid on                                 

% Histograms can be plotted directly as ** 
% h = histogram(trials,'BinLimits', [0.5 6.5], 'BinMethod', 'integers') **

trials_ = [10, 100, 1000, 100000];
figure('Name', '1.1 Rolling dice histogram')
for idx = 1:4
  subplot(4, 1, idx)
  histogram(randi(N_faces, 1, trials_(idx)),'BinLimits', [0.5 6.5], 'BinMethod', 'integers');
  grid on
  one_face_occurrence = trials_(idx)/N_faces;      
  hold on                                    
  plot([0.5 6.5], [one_face_occurrence one_face_occurrence], 'r') 
  hold off
  title(['Number of trial = ' num2str(trials_(idx))])
end
xlabel('Rolled number')
sgtitle('Dice rolling histogram (Number of occurrence)');
legend_ = legend('Experimental','Analytic');
legend_.Position(1) = 0.80;
legend_.Position(2) = 0.85;  


% Normalize the histogram values with N_trials (number of trials), 
% so that we get the experimental probability density function (pdf) for the dice. 
% After the normalization the sum of the histogram bin area ("integral") is equal to one, 
% as required with pdfs (when talking about discrete pdf, often refer to pmf-probability mass function).

pdf_experimental = hist_count/N_trials;

% Plot the normalized histogram, and define the true (analytic) pdf for the fair dice and
% plot it on top of the experimental pdf

figure ('Name', '1.1 Normalized histogram')
bar(x_histogram_centers, pdf_experimental)
xlabel('Rolled number')
ylabel('pdf')
grid on
title('Dice rolling normalized histogram')
one_face_probability = 1/N_faces;           % probability of one face of the dice
hold on                                     % avoids removing the previous plot
plot([0.5 6.5], [one_face_probability one_face_probability], 'r')  % plotting true pdf
hold off
legend('Experimental pdf','True(analytic) pdf')       % Using legend we can name different data curves in the plot (in order of appearance)

% Are the true pdf and the experimental pdf perfectly equal?
%   -> NO, not at all. (in case of number of trial = 10)

%  Repeat the process for a few of times and compare the outcomes:
%  Is the experimental pdf varying between simulations? 
%   -> YES, it varies for each attempt. 

%  Increase the number of trials to N=100000 and repeat the process for a few times
%  Is the experimental pdf varying between simulations now? 
%   -> NO, experimental pdf is very stable(similar) for each trial and values are quite close(equal) to true pdf 


%% 1.2

% Generate 1000 samples for a normally distributed random variable X~(μ,σ2), 
% where the mean and variance are given as μ=3, and σ^2=4

clear; clc; close all;

N1 = 100;
N2 = 100000;

N_samples = N2;

mu = 3;
sigma = sqrt(4);             

X = mu + sigma * randn(1, N_samples);

% Calculate the following statistics of the observed samples
%  -> Mean, Standard deviation, Variance 

fprintf('Mean value of generated samples = %6.4f\n', mean(X));
fprintf('Standard deviation of generated samples = %6.4f\n', std(X));
fprintf('Variance of generated samples = %6.4f\n\n', var(X));


% Define the histogram bin 
bin_width = 0.5;                        % bin width (in the x-axis)
bin_centers = -7:bin_width:13;          % x-axis for the bin center coordinates
bin_edges = (bin_centers(1)-bin_width/2):bin_width:(bin_centers(end)+bin_width/2);


% Plot the experimental pdf (i.e. the normalized histogram so that the sum of histogram bin area is equal to one)

[hist_count, ~] = histcounts(X, bin_edges);
pdf_experimental = hist_count/sum(hist_count*bin_width);

figure('Name', '1.2 Gaussian distributed random samples')
bar(bin_centers, pdf_experimental, 1)
hold on
title('Generated samples normalized histogram')
xlabel('X')
ylabel('pdf')
grid on

% Plot the true (analytic) pdf of X~(μ,σ2) on top of the experimental pdf

pdf_true = 1/(sqrt(2*pi)*sigma)*exp(-(mu-bin_edges).^2/(2*sigma^2));

plot(bin_edges, pdf_true, 'r','LineWidth', 2.5)   % defines a specific line width
hold off
legend('Experimental pdf','True pdf')             % Using legend for different data curves in the plot


% Repeat the experiment with different number of samples: N=100 and N=100000
% See the difference in the fitting between the experimental pdf and the analytic pdf as the number of samples changes
%  -> when number of samples is small (=100), difference between experimental and
%  analytic pdf is huge, but increasing number of samples, experimental pdf
%  getting similar to the analytic pdf


% Based on the histogram with N=100000 samples, define the probability P(X>5.25)
% Integrate (i.e. sum) histogram bins, whose x‐axis center values are larger than 5.25

b = 5.25;
indices_with_bin_center_larger_than_b = bin_centers > b;
considered_bin_values = pdf_experimental(indices_with_bin_center_larger_than_b);

probability_X_larger_than_b = sum(considered_bin_values*bin_width);

% Check your answer with the Q‐function P(X>b)= Q((b-μ)/σ)
analytic_probability = qfunc((b-mu)/sigma);

fprintf('Sum of histogram (experimental) of X>5.25 = %.5f\n', probability_X_larger_than_b)
fprintf('Probability (analytic) of X>5.25 using Q-function = %.5f\n', analytic_probability)


%% 2. RANDOM PROCESSES

% 2.1 WHITE NOISE VS. COLORED NOISE

% Generate a zero mean white Gaussian noise signal with variance σ2=3

clc; clear; close all;

N = 10000;                              % Number of generated samples
noise_var = 3;                          % Desired noise variance
noise = sqrt(noise_var)*randn(1, N);    % Noise signal generation

% Plot the noise signal
figure ('Name', '2.1 White noise signal')
plot(noise)
xlabel('sample index')
ylabel('noise amplitude')
title('White noise')
xlim([0 100])                           % define the x-axis limits

% Plot the histogram of the noise signal to see that it is Gaussian distributed
figure('Name', '2.1 White noise histogram')
histogram(noise, 40)
xlabel('noise amplitude')
ylabel('histogram count')
title('White noise histogram')


% Implement a low pass FIR filter (help firpm, recap from the previous exercise)
% Use filter order of 60 and a transition band from 0.1*Fs/2 to 0.2*Fs/2  

N_filter = 60;
h = firpm(N_filter, [0 0.1 0.2 1], [1 1 0 0]);

% Plot the impulse response and amplitude (frequency) response of the filter
% define the amplitude response as a function of normalized frequency (–Fs/2…Fs/2 → ‐1…1)

N_freq = length(noise);
freq_vec_filter = -1:2/N_freq:(1-2/N_freq);
figure('Name', '2.1 Impulse response of filter')
stem(-N_filter/2:N_filter/2, h);
title('Impulse response of the filter')

%frequency vector values normalized between -1 and 1
figure('Name', '2.1 Low pass FIR filter')
plot(freq_vec_filter, 10*log10(fftshift(abs(fft(h, N_freq)))))
xlabel('Normalized frequency (F_s/2=1)')
ylabel('Amplitude')
title('Amplitude response of the filter')

% Filter the noise signal using the above filter
filtered_noise = filter(h, 1, noise);

% Plot the white noise and filtered noise signals in the same figure and compare
% Remember the filter delay in order to align the signals properly in time
filtered_noise = filtered_noise(N_filter/2+1:end);          % remove the delay

N_samples_plot = length(filtered_noise);

figure('Name', '2.1 Filtering noise signal')
plot(noise(1:N_samples_plot))
hold on
plot(filtered_noise(1:N_samples_plot),'r')
legend('White noise','Colored noise')
xlabel('sample index')
ylabel('noise amplitude')
title('White noise and filtered (colored) noise')
xlim([0 100])    
hold off


% Plot the histogram of the filtered noise signal
figure('Name', '2.1 Filtered noise histogram')
histogram(filtered_noise, 40)
xlabel('noise amplitude')
ylabel('histogram count')
title('Filtered noise histogram')

fprintf('Mean of filtered noise = %.5f\n',  mean(filtered_noise));
fprintf('Standard deviation of filtered noise = %.5f\n', std(filtered_noise));

% How it is distributed? 
% What is the distribution of the output signal of an LTI system, if the input signal is Gaussian distributed?
% -> Output signal is also gaussian distribution, with standard deviation ~= 0.65 and mean of 0


% Compute and plot the autocorrelation function of the original white noise signal 
[corr_fun, lags] = xcorr(noise);

% we normalize the max-value to 1 and use stem-function in order to emphasize the impulse-like nature of the outcome
figure('Name', '2.1 Autocorrelation White noise')
stem(lags, corr_fun/max(corr_fun))
xlabel('\tau')          
ylabel('R(\tau)')
title('Autocorrelation of white noise signal')
xlim([-30 30])


% What type of autocorrelation function should the white noise have?
% -> peak is shown only when there is no(=zero) time difference (tau = 0), 
% -> there is no correlation at all delays, it means that this process is completely uncorrelated.


% Compute and plot the autocorrelation function of the filtered (i.e., colored) noise signal

[corr_fun, lags] = xcorr(filtered_noise);
figure('Name', '2.1 Autocorrelation Filtered noise')
stem(lags, corr_fun/max(corr_fun))
title('Autocorrelation of filtered noise signal')
xlabel('\tau')          
ylabel('R(\tau)')
xlim([-60 60])


% Compare this with the previously plotted impulse response of the filter
%  -> filtered noise signal has similar(same) shape of the impulse response
%  of the filter, correlation between samples is found.


% Plot the power spectra (amplitude spectrum squared) of the two noise signals in the same figure 
% use the decibel scale 20*log10(abs(.))

noise_abs_spec = 20*log10(abs(fft(noise(1:length(filtered_noise)))));
filtered_noise_abs_spec = 20*log10(abs(fft(filtered_noise)));

% Define the frequency vector values (normalized between -1 and 1):
freq_vec = -1:2/length(noise_abs_spec):1-2/length(noise_abs_spec);
figure('Name', '2.1 Noise spectra')
plot(freq_vec,fftshift(noise_abs_spec))
hold on
plot(freq_vec,fftshift(filtered_noise_abs_spec),'r')
hold off
xlabel('Normalized frequency (F_s/2=1)')
ylabel('power [dB]')
title('Noise spectra')
legend('White noise','Filtered (coloured) noise')


% Notice that both signals are still Gaussian distributed, but the power spectra are different.  
% Try to remember what is the connection between the correlation function and the power spectral density function?

% Make sure you understand the difference between the following concepts: Gaussian noise and white noise
% I.e., here both noise signals are Gaussian, but only the other one is white


%% 2.2 RANDOM WALK MODEL

% Let's consider a random sequence X[n], the so called random walk model

% where W[i] are binary i.i.d.(independent and identically distributed) random variables with
% probabilities P[W[i]=s]=p and P[W[i]=(–s)]=1–p

% Generate a random process of 2000 samples and 5000 realizations(=repeatitions, ensemble size)
% Use first p=0.5 and s=1

clear; clc; close all;

ensembles = [50, 500, 5000, 10000];
probs = [0.5, 0.51];
stepsizes = [1, 2];

% Number of samples for each realization and number of realizations
N_samples  = 2000;                
N_ensemble = ensembles(3);         

% Step probability and step size:
p = probs(1);                           % P(Wi=s) = p, P(Wi=-s) = 1-p
s = stepsizes(1);                       % step length
n = 1:N_samples;                        % vector of samples indices

% Generating matrix of randomly generated steps:
W = rand(N_ensemble, N_samples);        % (i.e. uniformly distributed random values between 0 and 1)
indices_with_positive_s = W < p;        % find out steps going "up"
W(indices_with_positive_s) = s;         % Define steps for going "up"
W(~indices_with_positive_s) = -s;       % Define steps for going "down"

% The overall "random walk" is achieved by taking the cumulative sum over the steps:
X = cumsum(W, 2);  % each row describes one random walk realization, so the sum is taken over the 2nd dimension


% Plot five example realizations of the random walk process
figure('Name', '2.2 Random Walk process')
for ind = 1:5
  subplot(5, 1, ind)
  plot(n, X(ind,:))
  ylabel('X(n)')
  grid on
  title(['Realization #' num2str(ind)])   %num2str converts a numerical value into a character value
end
xlabel('n')

% Here is a handy way to get a full screen figure (otherwise the figure might be to unclear):
set(gcf,'units','normalized','outerposition', [0 0 1 1])

% Calculate the ensemble mean and the ensemble variance of the process

% Plot the calculated ensemble mean E[X[n]] and variance E[(X[n]-E[X[n]])2] 
% with the theoretical values given as 
%   Theoretical mean:       E[X[n]] = ns(2p-1)
%   Theoretical variance:   E[(X[n] - E[X[n]])2] = np(2s)2(1-p)

% Notice that the ensemble mean and variance are now functions of the sequence index n ("=time"),
% so there will be a specific mean and variance for each sequence index

mean_theory = n*s*(2*p-1);          % Theoretical mean
var_theory = n*(2*s)^2*p*(1-p);     % Theoretical variance
mean_observed = mean(X);            % Empirical mean (i.e., what we observe)
var_observed = var(X);              % Empirical variance (i.e., what we observe)


figure('Name', ' 2.2 Mean Random Walk')
plot(n, mean_observed,'b','LineWidth',3)
hold on
plot(n, mean_theory,'r:','LineWidth',2)
hold off
legend('observed mean','theoretical mean')
ylim([-2 2])                        % set the axis limits in y-direction only
xlabel('n')
ylabel('Mean')
title('Mean over the sample index')


figure('Name', '2.2 Variance Random Walk')
plot(n, var_observed,'b','LineWidth',3)
hold on
plot(n, var_theory,'r:','LineWidth',2)
hold off
legend('observed variance','theoretical variance')
xlabel('n')
ylabel('Variance')
title('Variance over the sample index')


% Re‐run the process with different values of p, s, and try different number of realizations  
% p = 0.5, 0.51
% s = 1, 2
% N_ensembles = 50, 500, 5000, 10000

% What if the number of realizations is very low?
%  -> observed values are deviated a lot from theoretical values

% Is the process stationary (in wide sense)? Why?/Why not?
%  -> It is NOT wide sense stationary process, because variance depends on time, 
%     variance of random walk process increases linearly with n.
