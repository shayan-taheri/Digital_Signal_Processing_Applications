%% Shayan Taheri --- Project 3 of EEE 5513: Digital FIR Filter Design

%% Project Description

% Design of linear phase digital FIR (finite impulse response) filters
% using the Fourier series method.

% Center Wc of the ideal filter, "where applicable", in the center
% of the transition period to account for the smearing of the window.

% Use the minimum number of samples, unless specified otherwise.

%% Preparing Environment
clear all;
close all;
clc;
 
%% Question 1

% Design a Low Pass Filter (LPF):
% Number of h(n) samples is odd.
% fs = 48 (kHz), f1 = 9.2 (kHz), fp = 5.2 (kHz).
% Stopband Attenuation >= 50 (dB).
% Passband Ripple <= 0.3 (dB).
% Use the Hamming window if feasible.
% Print the coeffcient values of h(n).
% Plot the frequency response (magnitude and phase).
% Hamming Window: w(n) = 0.54 + [0.46 * cos(2*pi*n/L)] , |n| <= L/2

fs = 48e3; % Sampling Rate (Sampling Frequency)

% Wp = 2 * pi * F_passband/F_sampling --> Passband Angular Frequency
% Ws = 2 * pi * F_stopband/F_sampling --> Stopband Angular Frequency
% h[n] = sin(Wc * (n ? a))/(pi * (n ? a)) --> Linear Phase Finite Impulse Response

T = 1/fs; % Sampling Period
f_LPL = 9.2e3; % Lower Passband Limit (for Transition End Point)
% f_UPL = fn (Nyquist Limit) - f_LPL; % Upper Passband Limit
fp = 5.2e3; % Passband Edge (Starting Point for Transition)
passband_ripple = 0.3; % Passband Ripple
wc = 2*pi*(f_LPL + fp)/2; % Angular Cutoff Frequency (Middle Point in Transition)
transition_inter = 2*pi*(f_LPL-fp)/fs; % Transition Band Interval (Delta W)
L = ceil(6.6*pi/transition_inter-1) + 1; % Window Length
n = 0:1:L; % Number of Window Analysis Points
w = hamming(L+1); % Hamming Window Construction

% Creating a Sinc Function
h = sin(T * wc * (n-L/2)) ./ (pi * (n-L/2)); % Linear Phase Finite Impulse Response
h(L/2+1) = sinc(0)*(wc*T)/pi;

window_func = w'; % Window Function

H_trans = window_func .* h % Transfer Function

[H_out,W_out] = freqz(H_trans,1); % Frequency Response of Digital Filter

figure();

subplot(2,1,1);
plot((W_out/pi)*(fs/2),abs(H_out),'k');
xlabel('Frequency (Hz)','FontSize',14);
ylabel('|H_D(w)|','FontSize',16);
title('Filter Frequency Response','FontSize',14);

subplot(2,1,2);
plot((W_out/pi)*(fs/2),angle(H_out)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',14);
ylabel('\theta','FontSize',16);

%% Question 2

% Design a Band Pass Filter (BPF):
% Passband: 300 (Hz) --> 500 (Hz)
% Transition Band: 100 (Hz)
% Passband Ripple: 0.1 (dB)
% Stopband after: 60 (dB)
% fs = 2 (kHz)
% Use the Blackman window if feasible
% Print the coeffcient values of h(n), n = 0 --> L
% Plot the frequency response (magnitude and phase)
% Blackman Window: w(n) = 0.42 + [0.5 * cos((2*pi).*(n/L))] + [0.08 * cos((4*pi).*(n/L))]
% in Blackman Window: |n| <= L/2

fs = 2e3; % Sampling Rate (Sampling Frequency)
T = 1/fs; % Sampling Period
f_LPL = 600; % Lower Passband Limit (for Transition End Point)
fp = 500; % Passband Edge (Starting Point for Transition)
passband_ripple = 0.1; % Passband Ripple
f_cut_1 = (200+300)/2; % First Regular Cutoff Frequency (Middle Point in Transition)
w_cut_1 = 2*pi*f_cut_1; % First Angular Cutoff Frequency
f_cut_2 = (500+600)/2; % Second Regular Cutoff Frequency (Middle Point in Transition)
w_cut_2 = 2*pi*f_cut_2; % Second Angular Cutoff Frequency
transition_inter = 2*pi*(f_LPL-fp)/fs; % Transition Band Interval (Delta W)
L = round(11*pi/transition_inter-1) + 1; % Window Length
n = 0:1:L; % Number of Window Analysis Points
w = blackman(L+1); % Blackman Window Construction

% Creating a Sinc Function
h = (sin(T * w_cut_2 * (n-L/2)) ./ (pi * (n-L/2))) - (sin(T * w_cut_1 * (n-L/2)) ./ (pi * (n-L/2)));
h(L/2+1) = (sinc(0)*(w_cut_2*T)/pi) - (sinc(0)*(w_cut_1*T)/pi);

window_func = w'; % Window Function
H_trans = window_func .* h % Transfer Function
[H_out,W_out] = freqz(H_trans,1); % Frequency Response of Digital Filter

figure();

subplot(2,1,1);
plot((W_out/pi)*(fs/2),abs(H_out),'k');
xlabel('Frequency (Hz)','FontSize',14);
ylabel('|H_D(w)|','FontSize',16);
title('Filter Frequency Response','FontSize',14);

subplot(2,1,2);
plot((W_out/pi)*(fs/2),angle(H_out)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',14);
ylabel('\theta','FontSize',16);

%% Question 3

% Design a Hilbert transformer h(n), in which n = 0 to 40.
% Use a Kaiser window with alpha = 3 and Ws = 10 rad/sec
% Plot the magnitude and phase of Hd(W) for "W" from 0 to Ws/2.

L = 80; % Window Length
n = 0:1:L/2; % Number of Window Analysis Points
alpha = 3; % Coefficient in Kaiser Window
ws = 10; % Angular Sampling Frequency
fs = ws/(2*pi); % Regular Sampling Frequency
w = kaiser(L+1,alpha); % Kaiser Window Construction

% Impulse Response of Bandpass Hilbert Transform Filter
h = (2-2*cos(n.*pi)) ./ (2*pi.*n);
h(1) = 0;
h_part2 = h;

for k = 1:(L/2+1)
    h(k) = -h_part2(L/2+2-k);
end

h(L/2+1:L+1) = h_part2;

window_func = w'; % Window Function

H_trans = window_func .* h % Transfer Function

[H_out,W_out] = freqz(H_trans); % Frequency Response of Digital Filter

figure();

subplot(2,1,1);
plot((W_out/pi)*(ws/2),abs(H_out));
xlabel('\omega (rad/sec)','FontSize',16);
ylabel('|H_D(w)|','FontSize',16);
title('Filter Frequency Response','FontSize',14);

subplot(2,1,2);
plot((W_out/pi)*(ws/2),angle(H_out)*180/pi);
xlabel('\omega (rad/sec)','FontSize',16);
ylabel('\theta','FontSize',16);

%% Question 4

% Design a digital differentiator for "W" from 0 to Ws/2.
% Length (L) = 20 and fs = 1 (Hz)
% (a) Using rectangular window
% (b) Using Hamming window
% Obtain the weights for these cases.
% Plot the magnitude and phase from "W = 0" to "W = Ws/2" for (a), (b), and Ideal cases.

clear all;


L = 20; % Window Length
n = 0:1:L/2; % Number of Window Analysis Points
fs = 1; % Sampling Rate (Regular Sampling Frequency)
ws = 2*pi*fs; % Angular Sampling Frequency
T = 1/fs; % Sampling Period

% Impulse Response of Digital Differentiator Filter
h = cos(n*pi) ./ (n*T);
h(1) = 0;

h_part2 = h;

for k = 1:(L/2+1)
    h(k) = -h_part2(L/2+2-k);
end

h(L/2+1:L+1) = h_part2;

w_ideal = -ws/2:0.01:ws/2; % Ideal Window Construction
H_ideal = j*w_ideal % Transfer Function for Ideal Window

w_ham = hamming(L+1); % Hamming Window Construction
window_ham_func = w_ham'; % Hamming Window Function
H_ham = window_ham_func .* h % Transfer Function for Hamming Window
[Hd_ham,W_ham] = freqz(H_ham); % Frequency Response for Hamming Window

w_rec = rectwin(L+1); % Rectangular Window Construction
window_rec_func = w_rec'; % Rectangular Window Function
H_rec = window_rec_func .* h % Transfer Function for Rectangular Window
[Hd_rec,W_rec] = freqz(H_rec); % Frequency Response for Rectangular Window

figure();

subplot(2,1,1);
plot(w_ideal,abs(H_ideal));  
xlabel('\omega (rad/sec)','FontSize',16);
ylabel('|H_D(w)|','FontSize',16);
title('Filter Frequency Response for Ideal Window','FontSize',14);
 
subplot(2,1,2);
plot(w_ideal,angle(H_ideal)*180/pi);
xlabel('\omega (rad/sec)','FontSize',16);
ylabel('\theta','FontSize',16);

figure();

subplot(2,1,1);
plot((W_ham/pi)*(ws/2),abs(Hd_ham));  
xlabel('\omega (rad/sec)','FontSize',16);
ylabel('|H_D(w)|','FontSize',16);
title('Filter Frequency Response for Hamming Window','FontSize',14);

subplot(2,1,2);
plot((W_ham/pi)*(ws/2),angle(Hd_ham)*180/pi);
xlabel('\omega (rad/sec)','FontSize',16);
ylabel('\theta','FontSize',16);

figure();

subplot(2,1,1);
plot((W_rec/pi)*(ws/2),abs(Hd_rec));  
xlabel('\omega (rad/sec)','FontSize',16);
ylabel('|H_D(w)|','FontSize',16);
title('Filter Frequency Response for Rectangular Window','FontSize',14);

subplot(2,1,2);
plot((W_rec/pi)*(ws/2),angle(Hd_rec)*180/pi);
xlabel('\omega (rad/sec)','FontSize',16);
ylabel('\theta','FontSize',16);
