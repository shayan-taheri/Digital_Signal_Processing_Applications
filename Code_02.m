%% Shayan Taheri --- Project 2 of EEE 5513: Digital Filter Design
%% Preparing Environment
clear all;
close all;
clc;
 
%% Question 1

fs = 10000; % Sampling Frequency
f_var = 1;  % Regular Frequency Variation Size
w_var = 2*pi*f_var; % Angular Frequency Variation Size
fsp = 0:f_var:fs/2; % Specified Range of Regular Frequencies
w = 2*pi*fsp; % Specified Range of Angular Frequencies
 
f_cut = 2000; % Regular (3dB) Cutoff Frequency

[z,p,k] = cheb1ap(4,0.5); % Constructing Chebyshev Type I Analog Lowpass Filter
[b,a] = zp2tf(z,p,k); % Forming transfer function polynomials from the zeros, poles, and gains
[b,a] = lp2lp(b,a,2*pi*f_cut); % Transforming an analog low-pass filter into a low-pass filter with cutoff angular frequency
[bz,az] = impinvar(b,a,fs); % Creating a digital filter with the specified numerator and denominator coefficients

H = freqs(b,a,w); % Returning the complex frequency response of an analog filter specified at the specified angular frequencies

[H_D,w_D] = freqz(bz,az); % Returning the frequency response and its corresponding angular frequency vector

[z1,p1,k1] = tf2zp(bz,az); % Finding the zeros, poles, and gains of a transfer function

% Showing the magnitude and phase of the poles

disp (['Pole 1 Magnitude = ',num2str(abs(p1(1))),'; Pole 1 Phase = ',num2str(angle(p1(1))*180/pi)]);
disp (['Pole 2 Magnitude = ',num2str(abs(p1(2))),'; Pole 2 Phase = ',num2str(angle(p1(2))*180/pi)]);
disp (['Pole 3 Magnitude = ',num2str(abs(p1(3))),'; Pole 3 Phase = ',num2str(angle(p1(3))*180/pi)]);
disp (['Pole 4 Magnitude = ',num2str(abs(p1(4))),'; Pole 4 Phase = ',num2str(angle(p1(4))*180/pi)]);

figure();
subplot(2,1,1);
plot((w_D/pi)*(fs/2),abs(H_D),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D(e^{jwT})|','FontSize',14);
subplot(2,1,2);
plot((w_D/pi)*(fs/2),angle(H_D)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w/(2*pi),abs(H/H(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H(jw)/H(0)|','FontSize',14);
subplot(2,1,2);
plot(w/(2*pi),angle(H)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot((w_D/pi)*(fs/2), abs(H_D/H_D(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D(e^{jwT})/H_D(1)|','FontSize',14);
subplot(2,1,2);
plot((w_D/pi)*(fs/2),angle(H_D/H_D(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

%% Question 2

fs_2 = 20000; % Sampling Frequency for Question 2
f_var2 = 1; % Regular Frequency Variation Size
w_var2 = 2*pi*f_var2; % Angular Frequency Variation Size
fsp2 = 0:f_var2:fs_2/2; % Specified Range of Regular Frequencies
w2 = 2*pi*fsp2; % Specified Range of Angular Frequencies
 
f_cut2 = 2000; % Regular (3dB) Cutoff Frequency

[z2,p2,k2] = cheb1ap(4,0.5); % Constructing Chebyshev Type I Analog Lowpass Filter
[A,B,C,D] = zp2ss(z2,p2,k2); % Converting a zero-pole-gain representation of a given system to an equivalent state-space representation
[A1,B1,C1,D1] = lp2lp(A,B,C,D,2*pi*f_cut2); % Transforming an analog low-pass filter into a low-pass filter with cutoff angular frequency
[b2,a2] = ss2tf(A1,B1,C1,D1); % Converting a state-space representation of a system into an equivalent transfer function
[bz2,az2] = impinvar(b2,a2,fs_2); % Creating a digital filter with the specified numerator and denominator coefficients
 
H2 = freqs(b2,a2,w2); % Returning the complex frequency response of an analog filter specified at the specified angular frequencies

[H_D2,w_D2] = freqz(bz2,az2);

[z2,p2,k2] = tf2zp(bz2,az2);

% Showing the magnitude and phase of the poles

disp (['Pole 1 Magnitude = ',num2str(abs(p2(1))),'; Pole 1 Phase = ',num2str(angle(p2(1))*180/pi)]);
disp (['Pole 2 Magnitude = ',num2str(abs(p2(2))),'; Pole 2 Phase = ',num2str(angle(p2(2))*180/pi)]);
disp (['Pole 3 Magnitude = ',num2str(abs(p2(3))),'; Pole 3 Phase = ',num2str(angle(p2(3))*180/pi)]);
disp (['Pole 4 Magnitude = ',num2str(abs(p2(4))),'; Pole 4 Phase = ',num2str(angle(p2(4))*180/pi)]);

% Comparison of the Zpmax

disp (['Zpmax for Pole 1 = ',num2str(max(abs(p1))),'; Zpmax for Pole 2 = ',num2str(max(abs(p2)))]);

figure();
subplot(2,1,1);
plot(w_D2/pi*fs_2/2, abs(H_D2),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D2(e^{jwT})|','FontSize',14);
subplot(2,1,2);
plot(w_D2/pi*fs_2/2,angle(H_D2)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w2/(2*pi), abs(H2/H2(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_{2}(jw)/H_{2}(0)|','FontSize',14);
subplot(2,1,2);
plot(w2/(2*pi),angle(H2/H2(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);
 
figure();
subplot(2,1,1);
plot(w_D2/pi*fs_2/2, abs(H_D2/H_D2(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D2(e^{jwT})/H_D2(1)|','FontSize',14);
subplot(2,1,2);
plot(w_D2/pi*fs_2/2,angle(H_D2/H_D2(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

% Comparing the Appriximation Accuracy

figure();
subplot(2,1,1);
plot(w/pi/2, abs(H/H(1)),'b');
hold on;
plot(w_D/pi*fs/2, abs(H_D/H_D(1)),'k:');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Magnitude','FontSize',12);
title('|H(jw)/H(0)| and  |H_D(e^{jwT})/H_D(1)|','FontSize',14);
legend('|H(jw)/H(0)|','|H_D(e^{jwT})/H_D(1)|');
hold off;
subplot(2,1,2);
plot(w/pi/2,angle(H/H(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);
hold on;
plot(w_D/pi*fs/2,angle(H_D/H_D(1))*180/pi,'k:');
legend('|H(jw)/H(0)|','|H_D(e^{jwT})/H_D(1)|');
hold off;
 
figure();
subplot(2,1,1);
plot(w2/(2*pi), abs(H2/H2(1)),'b');
hold on;
plot(w_D2/pi*fs_2/2, abs(H_D2/H_D2(1)),'k:');
hold off;
title('|H_{2}(jw)/H_{2}(0)| and |H_D2(e^{jwT})/H_D(1)|','FontSize',14);
ylabel('Magnitude','FontSize',12);
xlabel('Frequency (Hz)','FontSize',12);
legend('|H_{2}(jw)/H_{2}(0)|','|H_D2(e^{jwT})/H_D(1)|');
subplot(2,1,2);
plot(w2/(2*pi),angle(H2/H2(1))*180/pi,'b');
hold on;
plot(w_D2/pi*fs_2/2,angle(H_D2/H_D2(1))*180/pi,'k:');
ylabel('\theta','FontSize',14);
xlabel('Frequency (Hz)','FontSize',12);
legend('|H_{2}(jw)/H_{2}(0)|','|H_D2(e^{jwT})/H_D(1)|');
hold off;

%% Question 3

fs3 = 15000;  % Sampling Frequency for Question 3
f_var3 = 1; % Regular Frequency Variation Size
w_var3 = 2*pi*f_var3; % Angular Frequency Variation Size
fsp3 = 0:f_var3:fs3/2; % Specified Range of Regular Frequencies
w3 = 2*pi*fsp3; % Specified Range of Angular Frequencies
 
f_cut3 = 3000; % Regular (3dB) Cutoff Frequency

[z3,p3,k3] = cheb1ap(4,0.5);
[A,B,C,D] = zp2ss(z3,p3,k3);
[A1,B1,C1,D1] = lp2lp(A,B,C,D,2*pi*f_cut3);
[b3,a3] = ss2tf(A1,B1,C1,D1);
[bz3,az3] = impinvar(b3,a3,fs3);
 
H3 = freqs(b3,a3,w3);

[H_D3,w_D3] = freqz(bz3,az3);

[z3,p3,k3] = tf2zp(bz3,az3);

disp (['Pole 1 Magnitude = ',num2str(abs(p3(1))),'; Pole 1 Phase = ',num2str(angle(p3(1))*180/pi)]);
disp (['Pole 2 Magnitude = ',num2str(abs(p3(2))),'; Pole 2 Phase = ',num2str(angle(p3(2))*180/pi)]);
disp (['Pole 3 Magnitude = ',num2str(abs(p3(3))),'; Pole 3 Phase = ',num2str(angle(p3(3))*180/pi)]);
disp (['Pole 4 Magnitude = ',num2str(abs(p3(4))),'; Pole 4 Phase = ',num2str(angle(p3(4))*180/pi)]);

figure();
subplot(2,1,1);
plot(w_D3/pi*fs3/2, abs(H_D3),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D3(e^{jwT})|','FontSize',14);
subplot(2,1,2);
plot(w_D3/pi*fs3/2,angle(H_D3)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w3/(2*pi),abs(H3/H3(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H3(jw)/H3(0)|','FontSize',14);
subplot(2,1,2);
plot(w3/(2*pi),angle(H3/H3(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w_D3/pi*fs3/2, abs(H_D3/H_D3(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D3(e^{jwT})/H_D3(1)|','FontSize',14);
subplot(2,1,2);
plot(w_D3/pi*fs3/2,angle(H_D3/H_D3(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

%% Question 4

fs4 = 40000;  % Sampling Frequency for Question 4
f_var4 = 1; % Regular Frequency Variation Size
w_var4 = 2*pi*f_var4; % Angular Frequency Variation Size
fsp4 = 0:f_var4:fs4/2; % Specified Range of Regular Frequencies
w4 = 2*pi*fsp4; % Specified Range of Angular Frequencies
 
f_cut4 = 4000; % Regular (3dB) Cutoff Frequency

[z4,p4,k4] = cheb1ap(4,0.5);
[A,B,C,D] = zp2ss(z4,p4,k4);
[A1,B1,C1,D1] = lp2lp(A,B,C,D,2*pi*f_cut4);
[b4,a4] = ss2tf(A1,B1,C1,D1);
[bz4,az4] = impinvar(b4,a4,fs4);
 
H4 = freqs(b4,a4,w4);

[H_D4,w_D4] = freqz(bz4,az4);

[z4,p4,k4] = tf2zp(bz4,az4);

disp (['Pole 1 Magnitude = ',num2str(abs(p4(1))),'; Pole 1 Phase = ',num2str(angle(p4(1))*180/pi)]);
disp (['Pole 2 Magnitude = ',num2str(abs(p4(2))),'; Pole 2 Phase = ',num2str(angle(p4(2))*180/pi)]);
disp (['Pole 3 Magnitude = ',num2str(abs(p4(3))),'; Pole 3 Phase = ',num2str(angle(p4(3))*180/pi)]);
disp (['Pole 4 Magnitude = ',num2str(abs(p4(4))),'; Pole 4 Phase = ',num2str(angle(p4(4))*180/pi)]);

figure();
subplot(2,1,1);
plot(w_D4/pi*fs4/2, abs(H_D4),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D4(e_{jwT})|','FontSize',14);
subplot(2,1,2);
plot(w_D4/pi*fs4/2,angle(H_D4)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w4/(2*pi),abs(H4/H4(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H4(jw)/H4(0)|','FontSize',14);
subplot(2,1,2);
plot(w4/(2*pi),angle(H4/H4(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);
 
 
figure();
subplot(2,1,1);
plot(w_D4/pi*fs4/2,abs(H_D4/H_D4(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D4(e^{jwT})/H_D4(1)|','FontSize',14);
subplot(2,1,2);
plot(w_D4/pi*fs4/2,angle(H_D4/H_D4(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

%% Question 5
 
fs5 = 10000;   % Sampling Frequency for Question 5
f_var5 = 1; % Regular Frequency Variation Size
w_var5 = 2*pi*f_var5; % Angular Frequency Variation Size
fsp5 = 0:f_var5:fs5/2; % Specified Range of Regular Frequencies
w5 = 2*pi*fsp5; % Specified Range of Angular Frequencies
 
f_cut5 = 2000; % Regular (3dB) Cutoff Frequency

[z,p,k] = cheb1ap(4,0.5);
[A,B,C,D] = zp2ss(z,p,k);
[A1,B1,C1,D1] = lp2lp(A,B,C,D,2*pi*f_cut5);
[b,a] = ss2tf(A1,B1,C1,D1);
[bz5,az5] = bilinear(b,a,fs5); % Converting an s-domain transfer function to a discrete equivalent
 
H5 = freqs(b,a,w5);

[H_D5,w_D5] = freqz(bz5,az5);

f_cut5_new = atan(f_cut*2*pi/fs/2)*2*fs/(2*pi);

[z5,p5,k5] = tf2zp(bz5,az5);

disp (['New Regular Cutoff Frequency of H_D5 = ',num2str(f_cut5_new)]);

disp (['Pole 1 Magnitude = ',num2str(abs(p5(1))),'; Pole 1 Phase = ',num2str(angle(p5(1))*180/pi)]);
disp (['Pole 2 Magnitude = ',num2str(abs(p5(2))),'; Pole 2 Phase = ',num2str(angle(p5(2))*180/pi)]);
disp (['Pole 3 Magnitude = ',num2str(abs(p5(3))),'; Pole 3 Phase = ',num2str(angle(p5(3))*180/pi)]);
disp (['Pole 4 Magnitude = ',num2str(abs(p5(4))),'; Pole 4 Phase = ',num2str(angle(p5(4))*180/pi)]);

figure();
subplot(2,1,1);
plot(w_D5/pi*fs5/2, abs(H_D5),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D5(e^{jwT})|','FontSize',14);
subplot(2,1,2);
plot(w_D5/pi*fs5/2,angle(H_D5)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w5/(2*pi),abs(H5/H5(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H5(jw)/H5(0)|','FontSize',14);
subplot(2,1,2);
plot(w5/(2*pi),angle(H5/H5(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w_D5/pi*fs5/2,abs(H_D5/H_D5(1)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D5(e^{jwT})/H_D5(1)|','FontSize',14);
subplot(2,1,2);
plot(w_D5/pi*fs5/2,angle(H_D5/H_D5(1))*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

%% Question 6

w_3dB = 10000; % Angular (3dB) Cutoff Frequency
w_s = 4 * w_3dB; % Angular Sampling Frequency
f_s = w_s/(2*pi); % Regular Sampling Frequency
f_bw_var = 1; % Butterworth Regular Frequency Variation Size
f_bw = 0:f_bw_var:f_s/2; % Butterworth Regular Frequency
w_bw = 2*pi*f_bw; % Butterworth Angular Frequency
b = [0 0 0 0 0 0 1]; % The coefficients of the filter numerator
a = [1 3.8637033 7.464106 9.1416202 7.464106 3.8637033 1]; % The coefficients of the filter denominator
[bn,an] = lp2hp(b,a,w_3dB); % Transforming an analog low-pass filter prototype given by polynomial coefficients into a low-pass filter with cutoff angular frequency Wo
[bt,at] = bilinear(bn,an,f_s); % Converting an s-domain transfer function to a discrete equivalent
H_bw = freqs(bn,an,w_bw); % Returning the complex frequency response of an analog filter specified at the specified angular frequencies
[Z, P, k] = tf2zp(bt, at); % Finding the zeros, poles, and gains of a transfer function
[H_Dbw, w_Dbw] = freqz(bt, at); % Returning the frequency response and its corresponding angular frequency vector

w_3db_new = atan(w_3dB/(2*f_s))* (2*f_s);

disp (['New Angular Cutoff Frequency of Butterworth-based Digital Filter = ',num2str(w_3db_new)]);

disp (['Pole 1 Magnitude = ',num2str(abs(P(1))),'; Pole 1 Phase = ',num2str(angle(P(1))*180/pi)]);
disp (['Pole 2 Magnitude = ',num2str(abs(P(2))),'; Pole 2 Phase = ',num2str(angle(P(2))*180/pi)]);
disp (['Pole 3 Magnitude = ',num2str(abs(P(3))),'; Pole 3 Phase = ',num2str(angle(P(3))*180/pi)]);
disp (['Pole 4 Magnitude = ',num2str(abs(P(4))),'; Pole 4 Phase = ',num2str(angle(P(4))*180/pi)]);

figure();
subplot(2,1,1);
plot(w_Dbw * f_s, 20*log10(abs(H_Dbw)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel({'H_D6(e^{jwT})'},'FontSize',14);
subplot(2,1,2);
plot(w_Dbw * f_s,angle(H_Dbw)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w_bw/(2*pi),abs(H_bw/H_bw(end)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H6(jw)/H6(max)|','FontSize',14);
subplot(2,1,2);
plot(w_bw/(2*pi),angle(H_bw/H_bw(end)*180/pi),'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot((w_Dbw)/pi*(f_s/2),abs(H_Dbw/H_Dbw(end)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_D6(e^{jwT})/H_D6(max)|','FontSize',14);
subplot(2,1,2);
plot((w_Dbw)/pi*(f_s/2),angle(H_Dbw/H_Dbw(end)*180/pi),'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();

subplot(2,1,1);
plot((w_Dbw)/pi*(f_s/2),abs(H_Dbw/H_Dbw(end)),'b');
hold on;
plot(w_bw/(2*pi),abs(H_bw/H_bw(end)),':k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Magnitude','FontSize',14);
legend('|H_D6(e^{jwT})/H_D6(max)|','|H6(jw)/H6(max)|');

subplot(2,1,2);
plot((w_Dbw)/pi*(f_s/2),angle(H_Dbw/H_Dbw(end)*180/pi),'b');
hold on;
plot(w_bw/(2*pi),angle(H_bw/H_bw(end)*180/pi),':k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);
legend('|H_D6(e^{jwT})/H_D6(max)|','|H6(jw)/H6(max)|');

%% Question 7

w_3dB = 10000; % Angular (3dB) Cutoff Frequency
w_s = 4 * w_3dB; % Angular Sampling Frequency
f_s = w_s/(2*pi); % Regular Sampling Frequency
f_bw_var = 1; % Butterworth Regular Frequency Variation Size
f_bw = 0:f_bw_var:f_s/2; % Butterworth Regular Frequency
w_bw = 2*pi*f_bw; % Butterworth Angular Frequency
w_warp = 2 * f_s * tan(w_3dB/(2*f_s)); % Angular Frequency for Warping
f_warp = w_warp/(2*pi); % Regular Frequency for Warping
b = [0 0 0 0 0 0 1]; % The coefficients of the filter numerator
a = [1 3.8637033 7.464106 9.1416202 7.464106 3.8637033 1]; % The coefficients of the filter denominator
[bn,an] = lp2hp(b,a,w_3dB); % Transforming an analog low-pass filter prototype given by polynomial coefficients into a low-pass filter with cutoff angular frequency Wo

% Converting an s-domain transfer function to a discrete equivalent
% with considering the optional matching frequency for warping
[bt, at] = bilinear(bn,an,f_s,f_warp);

[Z, P, k1] = tf2zp(bt, at); % Finding the zeros, poles, and gains of a transfer function
H_warp = freqs(bn,an,w_bw); % Returning the complex frequency response of an analog filter specified at the specified angular frequencies
[H_Dwarp, w_Dwarp] = freqz(bt,at); % Returning the frequency response and its corresponding angular frequency vector

disp (['Pole 1 Magnitude = ',num2str(abs(P(1))),'; Pole 1 Phase = ',num2str(angle(P(1))*180/pi)]);
disp (['Pole 2 Magnitude = ',num2str(abs(P(2))),'; Pole 2 Phase = ',num2str(angle(P(2))*180/pi)]);
disp (['Pole 3 Magnitude = ',num2str(abs(P(3))),'; Pole 3 Phase = ',num2str(angle(P(3))*180/pi)]);
disp (['Pole 4 Magnitude = ',num2str(abs(P(4))),'; Pole 4 Phase = ',num2str(angle(P(4))*180/pi)]);

figure();
subplot(2,1,1);
plot(w_Dwarp * f_s, 20*log10(abs(H_Dwarp)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel({'H_{D-warp}(e^{jwT})'},'FontSize',14);
subplot(2,1,2);
plot(w_Dwarp * f_s,angle(H_Dwarp)*180/pi,'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot(w_bw/(2*pi),abs(H_warp/H_warp(end)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_{warp}(jw)/H_{warp}(max)|','FontSize',14);
subplot(2,1,2);
plot(w_bw/(2*pi),angle(H_warp/H_warp(end)*180/pi),'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();
subplot(2,1,1);
plot((w_Dwarp)/pi*(f_s/2),abs(H_Dwarp/H_Dwarp(end)),'k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('|H_{D-warp}(e^{jwT})/H_{D-warp}(max)|','FontSize',14);
subplot(2,1,2);
plot((w_Dwarp)/pi*(f_s/2),angle(H_Dwarp/H_Dwarp(end)*180/pi),'b');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);

figure();

subplot(2,1,1);
plot((w_Dwarp)/pi*(f_s/2),abs(H_Dwarp/H_Dwarp(end)),'b');
hold on;
plot(w_bw/(2*pi),abs(H_warp/H_warp(end)),':k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Magnitude','FontSize',14);
legend('|H_{D-warp}(e^{jwT})/H_{D-warp}(max)|','|H_{warp}(jw)/H_{warp}(max)|');

subplot(2,1,2);
plot((w_Dwarp)/pi*(f_s/2),angle(H_Dwarp/H_Dwarp(end)*180/pi),'b');
hold on;
plot(w_bw/(2*pi),angle(H_warp/H_warp(end)*180/pi),':k');
xlabel('Frequency (Hz)','FontSize',12);
ylabel('\theta','FontSize',14);
legend('|H_{D-warp}(e^{jwT})/H_{D-warp}(max)|','|H_{warp}(jw)/H_{warp}(max)|');
