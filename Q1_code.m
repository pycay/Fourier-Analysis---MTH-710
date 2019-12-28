clear; close all; clc;
% Load data from file
load sample009.mat
% Compute the number of data
n = length(sample);
% Set the sample frequency
Fs = 8000;
% Compute the sample period
T = 1/Fs;
% Define an array with sample times
time = 1/Fs:1/Fs:n*T;
% Plotting data points
figure(1)
plot(time,sample)
title('Sampled data')
xlabel('time[s]')
ylabel('X(t)')
% Compute the fourier transform of the sampled signal
Y = fft(sample);
% Compute the two sided spectrum
P2 = abs(Y/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1); %convert two sided spectrum to single sided spectrum
% Define the frequency domain
f = Fs*(0:(n/2))/n;
figure(2)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('amp')
% Order the spectrum by descending amplitude
% B stores the frequency components in order of Amplitude (first component
% is the most important and so on).
[B,I] = sort(P1,'descend'); %Arrangement I indicates which position each component is in, once ordered
%[B,I] = sort(___)  returns a collection of index vectors 
%for any of the previous syntaxes. I is the same size as A and 
%describes the arrangement of the elements of A into B 
%along the sorted dimension. For example, if A is a vector, then B = A(I).
B = sort(P1, 'descend');

% N, the number of term in the sum, is the number of entries in the array
% P1. Then, N = length(P1) = 16000.
N_number_of_terms_in_sum=length(P1)
% Each point/bin in the FFT output array is spaced by the frequency
% resolution: delta_f = Fs/n = 8000/32000 = 0.25 Hz
 delta_f = Fs/n
% Then fn = n*delta_f
fn=n*delta_f;
fn
Highest_Amp=B(1)
freq_of_highest_amp=I(1)*delta_f

% Amplitude An of each term is in the n-th position of P1.