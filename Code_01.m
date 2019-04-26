%% Shayan Taheri --- Project 1 of EEE 5513: A Signal Processing System
%% Preparing Environment
clear all;
close all;
clc;
 
f_3dB = 2000; % Regular Cut-Off Frequency
w_3dB = 2*pi*f_3dB; % Angular Cut-Off Frequency
fs = 8000;
Ts = 1/fs; % Sampling Time and Regular Sampling Frequency
ws = 2*pi/Ts; % Angular Sampling Frequency
 
f_var = 1;  % Regular Frequency Variation Size
w_var = 2*pi*f_var; % Angular Frequency Variation Size
f_sp = 0:f_var:fs+2000; % Specified Range of Regular Frequencies
w = 2*pi*f_sp; % Specified Range of Angular Frequencies
 
n = 1; % The initial order of the Butterworth Low-Pass Filter
SNR = 0; % Signal to Noise Ratio
 
while (SNR < 60)
    
    [z,p,k] = buttap(n);   % Building the Butterworth Low-Pass Filter
    [b,a] = zp2tf(z,p,k);  % Forming Transfer Function Polynomials
    [b,a] = lp2lp(b,a,w_3dB); % Transforming the Transfer Function to the 3dB Cut-Off Frequency 
    f_sp = -fs:f_var:fs; % The specified regular sampling frequency range for the transfer function
    w = 2*pi*f_sp;

    % Returning the complex frequency response of an analog filter specified at an angular frequencies
    % NODE 2 Signal with abs(X1(w)) = 1, Angle = 0
    X2 = freqs(b,a,w);

    % Shifted Versions of X2(w)
    X2_Shift = zeros(1,length(X2));
    
    for NPoints = 1:8000;
        X2_Shift(1,NPoints) = X2(1,NPoints+8001);
        X2_Shift(1,NPoints+8000) = X2(1,NPoints);
    end
    
    X3 = X2+X2_Shift; % Signal Spectrum at Node 3
 
    % Calculate the SNR
    S2 = 0;
    N2 = 0;

    for NPoints=8001:f_3dB+8000
        
        S2=S2+abs(X2(NPoints))^2*w_var; % Signal Power at Node 2
        N2=N2+abs(X2_Shift(NPoints))^2*w_var; % Noise Power at Node 2
        
    end
    
    SNR = 10*log10(S2/N2);

    n = n+1; % Updating the filter order

end
 
%% Part A
disp (['LPF Order = ',num2str(n-1),'; Signal to Noise Ratio = ',num2str(SNR)]);
 
%% Part B
figure();
subplot(2,1,1);
plot(w/(2*pi),abs(X2),'k');
title('X2(\omega) at Butterworth LPF Output -- f_{3dB} = 2 (KHz), Order = 6','FontSize',12);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H(f)|','FontSize',14);
subplot(2,1,2);
plot(w/(2*pi),angle(X2)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

%% Part C
figure();
subplot(2,1,1);
plot(w/(2*pi),abs(X3),'k');
title('X3(\omega) at Ideal Sampler Output','FontSize',14);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|X3|','FontSize',14);
subplot(2,1,2);
plot(w/(2*pi),angle(X3)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);
 
%% Part D 
DSP_H = sawtooth(w/8000,0.5); % Constructing DSP Transfer Function, H(z) -- Step 1
DSP_H = (1-DSP_H) ./ 2; % Constructing DSP Transfer Function, H(z) -- Step 2

X4 = DSP_H .* X3; % Signal Spectrum at Node 4

figure();
subplot(2,1,1);
plot(w/(2*pi),abs(X4),'k');
title('X4(\omega) at DSP Output','FontSize',14);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|X4|','FontSize',14);
subplot(2,1,2);
plot(w/(2*pi),angle(X4)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

%% Part E
H_sh=sin(w.*Ts./2)./(w.*Ts./2); % Sample/Hold Transfer Function = Sinc(pi*f_sp*Ts)

X5 = X4 .* H_sh; % Signal Spectrum at Node 5
 
figure();
subplot(2,1,1);
plot(w/(2*pi),abs(X5),'k');
title('X5(\omega) at Sample/Hold Output','FontSize',14);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|X5|','FontSize',14);
subplot(2,1,2);
plot(w/(2*pi),angle(X5)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

%% Parts F and G
T_signal = 1/f_3dB; % Time Period (Minimum) for Signal Power

SNR_Max = 60; % Maximum Value of Signal to Noise Ratio

syms tau; % Defining the Tau Variable

tau = solve(10*log10(3/(8*pi^2)*(T_signal/tau)^2)-SNR_Max); % Solving Equation
tau = double(tau); % Making the result double precision

DAC_Bits = 12; % Number of Bits for Digital to Analog Converter

disp (['Tau (Max) = ',num2str(abs(tau(1))),'; Signal to Noise Quantization Ratio (SNQR) = ',num2str(6*DAC_Bits)]);
 
%% Part H

% f_3dB = 2000; --> Regular Cut-Off Frequency
% fs = 8000; --> Regular Sampling Frequency
% f_var = 1;  --> Regular Frequency Variation Size
% f_sp = -fs:f_var:fs; % The specified regular sampling frequency range for the transfer function
% w = 2*pi*f_sp; --> Specified Range of Angular Frequencies

Hr = zeros(1,length(w)); % Interpolation Filter Transfer Function

% Tau Hold = 125 us --> Frequency Hold = 8KHz

% Making the intermediate points equal to "One", and the rest equal to "Zero".
for NPoints=1:2*f_3dB
    
    Hr(NPoints+6001) = 1;

end
 
figure();
plot(w/(2*pi),Hr,'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_{R}|','FontSize',14);
title('Ideal Hr(\omega) -- Interpolation Filter','FontSize',14);

%% Part I
X6 = Hr .* X5; % Signal Spectrum at Node 6
 
figure();
subplot(2,1,1);
plot(w/(2*pi),abs(X6),'k');
title('X6(\omega) at Interpolation Filter Output','FontSize',14);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|X6|','FontSize',14);
subplot(2,1,2);
plot(w/(2*pi),angle(X6)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);
