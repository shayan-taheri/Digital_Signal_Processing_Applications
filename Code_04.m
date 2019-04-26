%% Shayan Taheri --- Project 4 of EEE 5513: Fourier Analysis

%% Project Description

%{

1) Compute 16-point DFT and IDFT (magnitude and phase) of x1(n) and x2(n),
where:

x1(n) = cos(2*pi*n/16) + 2*cos(4*pi*n/16)
x2(n) = sin(2*pi*n/16) + 2*sin(4*pi*n/16)

a. Using the defined equations for DFT and IDFT.
b. Your own FFT (Decimation in Time)

2) For n = 0 to n = 15, compute y(n) using:

a. Difference Equation
b. Convolution Summation
c. FFT (by overlap and add method in your own 16-point FFT)

* x(n) = 4*sin(2*pi*n/8)
* h(n) = alpha ^ n, which is a decreasing exponential discrete signal

3) Given x(n) = 3*sin(6*pi*n/16), compute and compare the energy in the
signal from its 8-point discrete time-domain and DFT representations.

4) x(t) is a periodic triangle signal with amplitude range of 0 to 1 and
period of 1 us. Obtain the Fourier series coefficients (complex
exponential) from D.C. to the tenth harmonic.

%}

%% Preparing Environment
clear all;
close all;
clc;
 
%% Question 1

N = 16; % Total number of sample points for the DFT/FFT analysis

n = 0:N-1; % Range of the number of sample points for DFT/FFT expression

x1 = cos((2*pi.*n)./16) + 2*cos((4*pi.*n)./16); % Making x1(n) signal
x2 = sin((2*pi.*n)./16) + 2*sin((4*pi.*n)./16); % Making x2(n) signal

X1 = zeros(N); % Vector for DFT of x1(n) using its formula
X2 = zeros(N); % Vector for DFT of x2(n) using its formula

x1_inv = zeros(N); % Vector for IDFT of X1(k) using its formula
x2_inv = zeros(N); % Vector for IDFT of X1(k) using its formula

% DFT Formula: X(k) = Sigma{n=0,N-1} x(n).exp(-2*pi*j*k*n/N)
% N = Number of sample points
% n = Current sample for analysis
% x(n) = Signal value at discrete time of "n"
% k = Frequency of result signal
% X(k) = Represents amplitude and phase of the frequency-domain signal

% Computing the forward discrete Fourier transform (DFT) of signals
% "k+1" and "n+1" --> Since the first array element starts from "1" in MATLAB
for k = 0:N-1 % The range of sum operation in DFT (Theory)
    X1(k+1) = sum(x1(n+1).*exp(-j*(2*pi./N)*k.*n));
    X2(k+1) = sum(x2(n+1).*exp(-j*(2*pi./N)*k.*n));
end

% Computing the inverse discrete Fourier transform (IDFT) of signals
% It is similar to DFT with the difference of having negative exponential
% weight and division of the expression by the total number of samples
for k = 0:N-1
    x1_inv(k+1) = (1/N)*sum(X1(n+1).*exp(j*(2*pi./N)*k.*n));
    x2_inv(k+1) = (1/N)*sum(X2(n+1).*exp(j*(2*pi./N)*k.*n));
end

% Part A: Using the defined equations for DFT and IDFT

% Figure 1
figure();
subplot(2,1,1)
stem(n,abs(x1), 'k');
title('Analysis of x_1 Signal','FontSize',12);
ylabel('Magnitude of x_1','FontSize',14);
xlabel('n','FontSize',16);
subplot(2,1,2)
stem(n, angle(x1), 'b');
ylabel('Phase of x_1','FontSize',14);
xlabel('n','FontSize',16);

% Figure 2
figure();
subplot(2,1,1)
stem(n,abs(x1_inv), 'k');
title('Inverse Fourier Analysis of X_1 Signal','FontSize',12);
ylabel('Magnitude of IDFT of X_1','FontSize',14);
xlabel('n','FontSize',16);
subplot(2,1,2)
stem(n, angle(x1_inv), 'b');
ylabel('Phase of IDFT of X_1','FontSize',14);
xlabel('n','FontSize',16);

% Figure 3
figure();
subplot(2,1,1)
stem(n,abs(x2), 'k');
title('Analysis of x_2 Signal','FontSize',12);
ylabel('Magnitude of x_2','FontSize',14);
xlabel('n','FontSize',16);
subplot(2,1,2)
stem(n,angle(x2),'b');
ylabel('Phase of x_2','FontSize',14);
xlabel('n','FontSize',16);

% Figure 4
figure();
subplot(2,1,1)
stem(n, abs(x2_inv), 'k');
title('Inverse Fourier Analysis of X_2 Signal','FontSize',12);
ylabel('Magnitude of IDFT of X_2','FontSize',14);
xlabel('n','FontSize',16);
subplot(2,1,2)
stem(n, angle(x2_inv), 'b');
ylabel('Phase of IDFT of X_2','FontSize',14);
xlabel('n','FontSize',16);

% Figure 5
figure();
subplot(2,1,1)
stem(n, abs(X1), 'k');
title('Fourier Analysis of X_1 Signal','FontSize',12);
ylabel('Magnitude of X_1','FontSize',14);
xlabel('k','FontSize',16);
subplot(2,1,2)
stem(n, angle(X1), 'b');
ylabel('Phase of X_1','FontSize',14);
xlabel('k','FontSize',16);

% Figure 6
figure();
subplot(2,1,1)
stem(n, abs(X2), 'k');
title('Fourier Analysis of X_2 Signal','FontSize',12);
ylabel('Magnitude of X_2','FontSize',14);
xlabel('k','FontSize',16);
subplot(2,1,2)
stem(n, angle(X2), 'b');
ylabel('Phase of X_2','FontSize',14);
xlabel('k','FontSize',16);

% Part B: Your own FFT (Decimation in Time)

X1_FFT = FFT_16(x1); % Getting FFT of x1(n) signal using own function
X2_FFT = FFT_16(x2); % Getting FFT of x2(n) signal using own function

% Figure 7
figure();
subplot(2,1,1)
stem(n, abs(X1_FFT), 'k');
title('Fast Fourier Analysis of X_1 Signal','FontSize',12);
ylabel('Magnitude of X_1','FontSize',14);
xlabel('k','FontSize',16);
subplot(2,1,2)
stem(n, angle(X1_FFT), 'b');
ylabel('Phase of X_1','FontSize',14);
xlabel('k','FontSize',16);

% Figure 8
figure();
subplot(2,1,1)
stem(n, abs(X2_FFT), 'k');
title('Fast Fourier Analysis of X_2 Signal','FontSize',12);
ylabel('Magnitude of X_2','FontSize',14);
xlabel('k','FontSize',16);
subplot(2,1,2)
stem(n, angle(X2_FFT), 'b');
ylabel('Phase of X_2','FontSize',14);
xlabel('k','FontSize',16);

%% Question 2

% M: The memory module for the input signal, x[n].
% It is used to perform "Difference Equations" or "Convolution Sum" filters.
M = zeros(3,16);

i = 0:3; % The power of scaling factor in the decreasing exponential function
h = 0.5 .^ i; % Discrete time-domain version of the filter, containing 4 elements.

% Difference Equation Filter: y[n] = x[n] + y[n-1]

% Convolution Sum Filter: y[n] = x[n]*h[n] = Sigma{k = -inf,+inf} x[k].h[n-k]
% ... y[n] = Sigma{k = -inf,+inf} x[n-k].h[k]

for n=0:N-1
    
    % Making 4 copies of the input signal corresponding to 4 elements of
    % the "H" filter
    x(n+1) = 4*sin((2*pi*n)/8);
    M(1,n+2) = x(n+1);
    M(2,n+3) = x(n+1);
    M(3,n+4) = x(n+1);
    
    % Computing the output signal using the difference equation filter
    % y_last = [0.5 * x(n)] + [0.5^2 * x(n-1)] + [0.5^3 * x(n-2)]
    y_last = 0.5*M(1,n+1)+0.5^2*M(2,n+1)+ 0.5^3*M(3,n+1);
    y(n+1) = x(n+1) + y_last;
    
    % Computing the output signal using the convolution summation filter
    % y_convolution(n+1) = h * [x(n+1);x(n);x(n-1);x(n-2)] --> Convolution
    y_convolution(n+1) = h*[x(n+1);M(1,n+1);M(2,n+1); M(3,n+1)];

end

% Appending zeros to the filter vector to have the same size as
% the input signal
h = [h,zeros(1,12)];

X_fft = FFT_16(x); % Getting fast Fourier transform of the input signal
H_fft = FFT_16(h); % Getting fast Fourier transform of the filter
y_ifft = IFFT_16(1/N * (X_fft .* H_fft)); % Getting inverse FFT (IFFT) of the output

figure();
subplot(3,1,1);
stem(1:N,y,'k');
title('Computing y(n) using Difference Equation Filter','FontSize',12);
ylabel('Amplitude','FontSize',14);
xlabel('n','FontSize',16);
subplot(3,1,2);
stem(1:N,y_convolution,'b');
title('Computing y(n) using Convolution Summation Filter','FontSize',12);
ylabel('Amplitude','FontSize',14);
xlabel('n','FontSize',16);
subplot(3,1,3);
stem(1:N,y_ifft,'r');
title('Computing y(n) using Fast Fourier Transform','FontSize',12);
ylabel('Amplitude','FontSize',14);
xlabel('n','FontSize',16);

%% Question 3

% Considering 8 sample points for Fourier analysis
n = 0:7; % Range of the number of sample points for DFT/FFT representations

x_time = 3 * sin(6*pi*n/16); % Discrete time-domain signal
Energy_Time = sum(x_time.^2); % Signal energy calculation in time-domain

X_freq = fft(x_time,length(n)); % Getting FFT of the time-domain signal
Energy_Freq = (1/length(n))*sum(abs(X_freq).^2); % Signal energy calculation in frequency-domain

figure();
subplot(2,1,1);
stem(n,abs(x_time),'k');
title('Energy Analysis of Time-Domain Signal','FontSize',12);
ylabel('Magnitude','FontSize',14);
xlabel('n','FontSize',16);
subplot(2,1,2);
stem(n,angle(x_time),'b');
ylabel('Phase','FontSize',14);
xlabel('n','FontSize',16);

figure();
subplot(2,1,1);
stem(n,X_freq,'k');
title('Energy Analysis of Frequency-Domain Signal','FontSize',12);
ylabel('Magnitude','FontSize',14);
xlabel('k','FontSize',16);
subplot(2,1,2);
stem(n,angle(X_freq),'b');
ylabel('Phase','FontSize',14);
xlabel('k','FontSize',16);

%% Question 4

% Triangle Function (with zero DC): x[n] = 1 - (1 - |n|/N); |n| < N

K = 0.1:0.05:10; % Range of harmonics for the Fourier series of x(t)

Exact_FSeries = (cos(K*pi)-1) ./ ((pi*K).^2);

n_b = 0:63; % Range "b" of the number of sample points for DFT/FFT expression
Series_64_Time = 1 - abs(1 - n_b./32); % Making the 64-point time-domain series
Series_64_Freq = fft(Series_64_Time)/64; % Getting FFT of the 64-point time series

n_c = 0:255; % Range "c" of the number of sample points for DFT/FFT expression
Series_128_Time = 1 - abs(1 - n_c./128); % Making the 128-point time-domain series
Series_128_Freq = fft(Series_128_Time)/256; % Getting FFT of the 128-point time series

figure();
subplot(2,1,1);
plot(1:length(K),abs(Exact_FSeries),'k');
title('Analysis of Signal Fourier Series in Exact Mode','FontSize',12);
ylabel('Magnitude','FontSize',14);
ylim([-0.2,0.8]);
xlabel('n','FontSize',16);
subplot(2,1,2);
plot(1:length(K),angle(Exact_FSeries),'b');
ylabel('Phase','FontSize',14);
xlabel('n','FontSize',16);

figure();
subplot(2,1,1)
plot(n_b,abs(Series_64_Freq),'k');
title('Analysis of 64-Point DFT of Signal Fourier Series','FontSize',12);
ylabel('Magnitude','FontSize',14);
ylim([-0.2,0.8]);
xlabel('k','FontSize',16);
subplot(2,1,2)
plot(n_b, angle(Series_64_Freq),'b');
ylabel('Phase','FontSize',14);
xlabel('k','FontSize',16);

figure();
subplot(2,1,1)
plot(n_c,abs(Series_128_Freq),'k');
title('Analysis of 128-Point DFT of Signal Fourier Series','FontSize',12);
ylabel('Magnitude','FontSize',14);
ylim([-0.2,0.8]);
xlabel('k','FontSize',16);
subplot(2,1,2)
plot(n_c,angle(Series_128_Freq),'b');
ylabel('Phase','FontSize',14);
xlabel('k','FontSize',16);
