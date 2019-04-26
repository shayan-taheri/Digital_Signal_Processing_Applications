%% Shayan Taheri --- Project 5 of EEE 5513: Adaptive Digital Signal Processing

%% Preliminary Material

%{

*************** Adaptive Filtering System ***************

# General Description:

W(n+1) = W(n) + u(n).e(n).X(n)
e(n) = d(n) - W_transpose(n).X(n)

W(n) = [W_0(n);W_1(n);...;W_{L-1}(n)] --> Coefficient Vector

# The optimal version of the coefficient vector contains impulse responses.

X(n) = [X(n);X(n-1);...;X(n-L+1)] --> Input Signal Vector
d(n) --> Desired Output (Response) Signal
u(n) --> Step Size
eta(n)/e(n) --> Noise/Error Signal
L --> Vector Length

# Options for the Input Signal and the Desired Output Signal:

A. Collection from the physical environment.
B. Getting from the mathematical/statistical model of the physical
environment.

# How to specify "u(n) and L"? --> Using "Trial and Error" Method

The parameters are found after a number of tries (iterations).

# Analytical Model for Characterization:

The average behavior (statistical expectation value) of signal is used.

Example: Average Behavior of Coefficient Vector:
E{W_i(n)} = Integral{-inf,+inf} W.P_{W_i}(W,n).dW

# There is "joint probability distribution" between the input signal and
the output signal.

In a "joint probability distribution", two (or more) random processes have
their values within their specified range.

*************** Wiener Filter ***************

# Description:

- It is a stationary filter initially with considering the signals to be
stationary linear random processes.

- The filter becomes adaptive by recomputing and updating the coefficients
periodically.

- The filter coefficients are computed based on minimization of the average
squared difference between the actual filter output and the desired filter
output.

# Wiener Filter Representation:

y_estimate(m) = Sigma{k=0,k=P-1} W_k.X(m-k) = W_transpose.X

W_k = argmin E{e^2(m)} --> Minimization of the Mean Square Error

e(m) = y(m) - y_estimate(m) --> The Wiener filter error signal

y_estimate --> The least mean square (LMS) error-based estimation of the
desired filter output.
W_K --> Coefficient Vector
X --> Filter Input Signal
y --> Desired Filter Output Signal
m --> Discrete-Time Index
P --> Filter Length
P-1 --> Filter Order

# Z-Transform of Wiener Filter:

Y_estimate(z) = Sigma{k=0,k=P-1} W_k.X(z).z^{-k}

*************** Adaptive Filtering by Least-Mean-Square (LMS) Algorithms ***************

# LMS Algorithms:

A class of adaptive filters that can function alike to a desired filter
when its coefficients cause getting the least mean square of the error
signal (which is the difference between the desired output signal
and the actual output signal).

Its adaptation mechanism works based on the error at the current time.

# Adaptive Mechanism of LMS Algorithms:

- Function: Getting the optimal filter weights by updating the filter weights
in a way to reach as close as possible to the desired filter output.

- Process:

A. Initializing Weights by Zero (i.e. Weights = 0).

B. Gradient of E{e^2(n)} = ?

C. Updating Weights.

# Updating Weights:

if "Gradient Value" > 0 --> Decrease Weights

else if "Gradient Value" < 0 --> Increase Weights

Model: W_{n+1} = W_n - u.Gradient(E{e^2(n)})

# Convergence Factor (u)

Getting the optimal weights through updating them iteratively will not be
achieved completely, but a convergence is possible in mean.

The accuracy and the speed of convergence process is determined by the
convergence factor (u).

# LMS Algorithms <--> Wiener Filter

The LMS algorithms can be transformed into the Wiener filter by relaxing
the error criterion for reduction of the current sample error (instead of
the total error).

------------------

# LMS Adaptive Filter Definition (Using Average Behavior Analysis):

A. V(n+1) = V(n) - u(n).X(n).X_transpose(n).V(n) + u(n).eta(n).X(n)
B. E{V(n+1)} = E{V(n)} - u(n).E{X(n).X_transpose(n).V(n)} + u(n).E{eta(n).X(n)}

[V(n) --> Coefficient Error Vector] ~ [Err(n) --> Mean-Squared Error]
X(n) --> Input Signal
Err(n) --> E{e^2(n)}
u(n) --> Step Size --> Does not depend on X(n),d(n), and W(n).
V(n) = W(n) - W_{optimal}
eta(n) --> Noise Signal
u (Convergence Factor) --> Controls the accuracy and the speed of
convergence and adaptation.

- Characterization of System Parameters --> Obtaining Adequate Performance

*************** Adaptive Filtering by Homogeneous Systems ***************

# Homogeneity Condition

A change in the system input signal causes a corresponding change in the
system output signal.

X[n] <--> Y[n] ====> K.X[n] <--> K.Y[n]

X[n] --> Input Signal
Y[n] --> Output Signal
K --> Constant Value

# Homogeneous Function

A homogeneous function has multiplicative scaling behavior.

F(a.X_1,a.X_2,...,a.X_N) = (a^k).F(X_1,X_2,...,X_N)

K --> Degree of Homogeneity

Both "a" and "K" are constant numbers.

# Homogeneous Polynomial

A homogeneous polynomial is a polynomial whose all of its non-zero terms have the same total (sum) degree.

# Homogeneous Distribution

For Y = t.X,

Integral{R^n} f(t.X).p(X).dX = (t^K).Integral{R^n} f(X).p(X).dX <-->

<--> (t^(-n)).Integral{R^n} f(Y).p(Y/t).dY = (t^k).Integral{R^n} f(Y).p(Y).dY

p() --> Test Function (Usage for evaluation of the characteristics of
optimization algorithms, such as convergence rate).

# Homogeneous System

In a homogeneous system of linear algebraic equations, the right side of
all equations are filled by zeros.

A.X = 0

A --> A matrix of M * N size.
X = [X_1;X_2;...;X_N] --> A m vector of N * 1 size.

An always (or trivial) solution for the homogeneous system is "X = 0".

In another definition, a system of linear equations is "homogenous" if all
of the constant terms in their right side are equal to zero.

*************** White Gaussian Noise ***************

# Additive White Gaussian Noise (AWGN):

- A basic noise model for mathematical modeling of the effects of existing
random processes in nature.

- The model has constant power spectral density.

- It can be described as a sequence of uncorrelated random values.

# Properties of Additive White Gaussian Noise:

- Additive = It is added to the inherent system noise.

- White = It has uniform power across the system frequency band.

- Gaussian = It has a normal distribution in the time domain with an
average value of zero.

%}

%% Project Description

%{

Question 1. Obtain the Wiener weights for the following system:

x(n) --> Input White Gaussian Noise Signal

eta(n) = 0 --> Input Signal

F(z) = c_0 + c_1.z^(-1) + c_2.z^(-2) --> Unknown Channel --> Output: d(n)

For F(z), we have c_0 = 1, c_1 = 2, and c_2 = 3.

A(z) = W_0 + W_1.z^(-1) + W_2.z^(-2) --> Modeling System --> Output: y(n)

e(n) = d(n) + eta(n) - y(n) --> Output Signal

W* = [R]^(-1).[P]

Case A) Use 100 samples of the input signals.
Case B) USe 1000 samples of the input signals.

Question 2. Repeat the first question when the modeling system becomes
adaptive using the LMS algorithm.

Question 3. Repeat the first question when the modeling system becomes
adaptive using the homogeneous algorithm.

Question 4. Repeat the first question when:

* eta(n) = 0.2 * sin((2*pi*n)/16)
* The modeling system becomes adaptive using the LMS algorithm.

%}

%% Preparing Environment
clear all;
close all;
clc;
 
%% Question 1

% W* = [R]^(-1).[P] --> The Complex Conjugate of the Weight Vector

c = [1,2,3]; % The Coefficients of the Unknown Channel Model

n1 = 100; % For Processing "100" Number of Samples
x1 = randn(1,n1); % Generating the noise signal.
% Convolution of the noise signal and the unknown channel transfer function
d1 = conv(x1,c);
d1 = d1(1:n1); % Choosing 100 number of samples of the convoltuion result
% Calculation of "R" and "P" Matrices for Getting Weight Vector
for i = 1:3 % According to the number of coefficients.
   
   x1_shifted = [zeros(1,i-1),x1(1:n1-i+1)];
   r1(i) = (x1 * x1_shifted')/n1;
   P1(i) = (d1 * x1_shifted')/n1;
   
end

for k = 1:3

    for j = 1:3
        
       R1(k,j) = r1(abs(k-j)+1);
    
    end 

end

W1 = R1\P1';
display('Q1: The Wiener Filter Weights for 100 Number of Samples:');
display(W1);

n2 = 1000; % For Processing "1000" Number of Samples
x2 = randn(1,n2);
d2 = conv(x2,c);
d2 = d2(1:n2);
% Calculation of "R" and "P" Matrices for Getting Weight Vector
for i = 1:3 % According to the number of coefficients.
   
   x2_shifted = [zeros(1,i-1),x2(1:n2-i+1)];
   r2(i) = (x2 * x2_shifted')/n2;
   P2(i) = (d2 * x2_shifted')/n2;
   
end

for k = 1:3

    for j = 1:3
        
       R2(k,j) = r2(abs(k-j)+1);
    
    end 

end

W2 = R2\P2';
display('Q1: The Wiener Filter Weights for 1000 Number of Samples:');
display(W2);

%% Question 2

c = [1,2,3]; % The Coefficients of the Unknown Channel Model

n1 = 100; % For Processing "100" Number of Samples
x1 = randn(1,n1); % Generating the noise signal.
% Convolution of the noise signal and the unknown channel transfer function
d1 = conv(x1,c);
d1 = d1(1:n1); % Choosing 100 number of samples of the convoltuion result
mu=0.02; % The step size for updating weights.

% Initializing Filter Weights, which is done in accordance with 
% the order of the channel/filter transfer function.
for i = 1:3
   
    w0(i) = 0;
    w1(i) = 0;
    w2(i) = 0;
    
end

% Updating Filter Weights
for j = 3:n1
    
    y1(j) = w0(j)*x1(j) + w1(j)*x1(j-1) + w2(j)*x1(j-2); % Calculate the filter output.
    e1(j) = d1(j) - y1(j); % Calculate the error signal.
    
    % Calculate the weights (i.e. coefficients).
    w0(j+1) = w0(j) + 2*mu*x1(j)*e1(j);
    w1(j+1) = w1(j) + 2*mu*x1(j-1)*e1(j);
    w2(j+1) = w2(j) + 2*mu*x1(j-2)*e1(j);

end

display('Q2: The Adaptive Wiener Filter Weights for 100 Number of Samples:');
display([w0(n1+1),w1(n1+1),w2(n1+1)]);

figure();
plot(1:n1,w0(1:n1),'^:');
hold on;
plot(1:n1,w1(1:n1),'*:');
plot(1:n1,w2(1:n1),'o:');
hold off;
legend('w0(n)','w1(n)','w2(n)');
xlabel('Iteration Index','FontSize',16);
ylabel('Coefficient (Weight) Value','FontSize',14);
title('Q2: Analysis of the LMS-Based Adaptive Wiener Filter','FontSize',12);

% ------------------

n2 = 1000; % For Processing "1000" Number of Samples
x2 = randn(1,n2);
d2 = conv(x2,c);
d2 = d2(1:n2);
mu=0.02; % The step size for updating weights.

% Initializing Filter Weights, which is done in accordance with 
% the order of the channel/filter transfer function.
for i = 1:3
   
    w0(i) = 0;
    w1(i) = 0;
    w2(i) = 0;
    
end

% Updating Filter Weights
for j = 3:n2
    
    y2(j) = w0(j)*x2(j) + w1(j)*x2(j-1) + w2(j)*x2(j-2); % Calculate the filter output.
    e2(j) = d2(j) - y2(j); % Calculate the error signal.
    
    % Calculate the weights (i.e. coefficients).
    w0(j+1) = w0(j) + 2*mu*x2(j)*e2(j);
    w1(j+1) = w1(j) + 2*mu*x2(j-1)*e2(j);
    w2(j+1) = w2(j) + 2*mu*x2(j-2)*e2(j);

end

display('Q2: The Adaptive Wiener Filter Weights for 1000 Number of Samples:');
display([w0(n2+1),w1(n2+1),w2(n2+1)]);

figure();
plot(1:n2,w0(1:n2),'^:');
hold on;
plot(1:n2,w1(1:n2),'*:');
plot(1:n2,w2(1:n2),'o:');
hold off;
legend('w0(n)','w1(n)','w2(n)');
xlabel('Iteration Index','FontSize',16);
ylabel('Coefficient (Weight) Value','FontSize',14);
title('Q2: Analysis of the LMS-Based Adaptive Wiener Filter','FontSize',12);

%% Question 3

c = [1,2,3]; % The Coefficients of the Unknown Channel Model

n1 = 100; % For Processing "100" Number of Samples
x1 = randn(1,n1); % Generating the noise signal.
% Convolution of the noise signal and the unknown channel transfer function
d1 = conv(x1,c);
d1 = d1(1:n1); % Choosing 100 number of samples of the convoltuion result
mu=0.02; % The step size for updating weights.

% Initializing Filter Weights, which is done in accordance with 
% the order of the channel/filter transfer function.
for i = 1:3
   
    w0(i) = 0;
    w1(i) = 0;
    w2(i) = 0;
    
end

% Updating Filter Weights
for j = 3:n1
    
    y1(j) = w0(j)*x1(j) + w1(j)*x1(j-1) + w2(j)*x1(j-2); % Calculate the filter output.
    e1(j) = d1(j) - y1(j); % Calculate the error signal.
    mu = 1/(2*(x1(j)^2 + x1(j-1)^2 + x1(j-2)^2)); % Calculate the step size parameter (u).
    
    % Calculate the weights (i.e. coefficients).
    w0(j+1) = w0(j) + 2*mu*x1(j)*e1(j);
    w1(j+1) = w1(j) + 2*mu*x1(j-1)*e1(j);
    w2(j+1) = w2(j) + 2*mu*x1(j-2)*e1(j);

end

display('Q3: The Adaptive Wiener Filter Weights for 100 Number of Samples:');
display([w0(n1+1),w1(n1+1),w2(n1+1)]);

figure();
plot(1:n1,w0(1:n1),'^:');
hold on;
plot(1:n1,w1(1:n1),'*:');
plot(1:n1,w2(1:n1),'o:');
hold off;
legend('w0(n)','w1(n)','w2(n)');
xlabel('Iteration Index','FontSize',16);
ylabel('Coefficient (Weight) Value','FontSize',14);
title('Q3: Analysis of the Homogeneous-Based Adaptive Wiener Filter','FontSize',12);

% ------------------

n2 = 1000; % For Processing "1000" Number of Samples
x2 = randn(1,n2);
d2 = conv(x2,c);
d2 = d2(1:n2);
mu=0.02; % The step size for updating weights.

% Initializing Filter Weights, which is done in accordance with 
% the order of the channel/filter transfer function.
for i = 1:3
   
    w0(i) = 0;
    w1(i) = 0;
    w2(i) = 0;
    
end

% Updating Filter Weights
for j = 3:n2
    
    y2(j) = w0(j)*x2(j) + w1(j)*x2(j-1) + w2(j)*x2(j-2); % Calculate the filter output.
    e2(j) = d2(j) - y2(j); % Calculate the error signal.
    mu = 1/(2*(x2(j)^2 + x2(j-1)^2 + x2(j-2)^2)); % Calculate the step size parameter (u).
    
    % Calculate the weights (i.e. coefficients).
    w0(j+1) = w0(j) + 2*mu*x2(j)*e2(j);
    w1(j+1) = w1(j) + 2*mu*x2(j-1)*e2(j);
    w2(j+1) = w2(j) + 2*mu*x2(j-2)*e2(j);

end

display('Q3: The Adaptive Wiener Filter Weights for 1000 Number of Samples:');
display([w0(n2+1),w1(n2+1),w2(n2+1)]);

figure();
plot(1:n2,w0(1:n2),'^:');
hold on;
plot(1:n2,w1(1:n2),'*:');
plot(1:n2,w2(1:n2),'o:');
hold off;
legend('w0(n)','w1(n)','w2(n)');
xlabel('Iteration Index','FontSize',16);
ylabel('Coefficient (Weight) Value','FontSize',14);
title('Q3: Analysis of the Homogeneous-Based Adaptive Wiener Filter','FontSize',12);

%% Question 4

c = [1,2,3]; % The Coefficients of the Unknown Channel Model

n1 = 100; % For Processing "100" Number of Samples
x1 = randn(1,n1); % Generating the noise signal.
% Convolution of the noise signal and the unknown channel transfer function
d1 = conv(x1,c);
d1 = d1(1:n1); % Choosing 100 number of samples of the convoltuion result
mu=0.02; % The step size for updating weights.

% Initializing Filter Weights, which is done in accordance with 
% the order of the channel/filter transfer function.
for i = 1:3
   
    w0(i) = 0;
    w1(i) = 0;
    w2(i) = 0;
    
end

% Updating Filter Weights
for j = 3:n1
    
    eta(j) = 0.2*sin((2*pi*j)/16); % The non-zero input signal
    y1(j) = w0(j)*x1(j) + w1(j)*x1(j-1) + w2(j)*x1(j-2); % Calculate the filter output.
    e1(j) = d1(j) + eta(j) - y1(j); % Calculate the error signal.
    
    % Calculate the weights (i.e. coefficients).
    w0(j+1) = w0(j) + 2*mu*x1(j)*e1(j);
    w1(j+1) = w1(j) + 2*mu*x1(j-1)*e1(j);
    w2(j+1) = w2(j) + 2*mu*x1(j-2)*e1(j);

end

display('Q4: The Adaptive Wiener Filter Weights for 100 Number of Samples:');
display([w0(n1+1),w1(n1+1),w2(n1+1)]);

figure();
plot(1:n1,w0(1:n1),'^:');
hold on;
plot(1:n1,w1(1:n1),'*:');
plot(1:n1,w2(1:n1),'o:');
hold off;
legend('w0(n)','w1(n)','w2(n)');
xlabel('Iteration Index','FontSize',16);
ylabel('Coefficient (Weight) Value','FontSize',14);
title('Q4: Analysis of the LMS-Based Adaptive Wiener Filter for Non-Zero Input Signal','FontSize',12);

figure();
plot(1:n1,eta(1:n1),'^:');
hold on;
plot(1:n1,e1(1:n1),'*:');
hold off;
legend('\eta','e');
xlabel('Iteration Index','FontSize',16);
ylabel('Amplitude','FontSize',14);
title('Q4: Analysis of the LMS-Based Adaptive Wiener Filter for Non-Zero Input Signal','FontSize',12);

% ------------------

n2 = 1000; % For Processing "1000" Number of Samples
x2 = randn(1,n2);
d2 = conv(x2,c);
d2 = d2(1:n2);
mu=0.02; % The step size for updating weights.

% Initializing Filter Weights, which is done in accordance with 
% the order of the channel/filter transfer function.
for i = 1:3
   
    w0(i) = 0;
    w1(i) = 0;
    w2(i) = 0;
    
end

% Updating Filter Weights
for j = 3:n2
    
    eta(j) = 0.2*sin((2*pi*j)/16); % The non-zero input signal
    y2(j) = w0(j)*x2(j) + w1(j)*x2(j-1) + w2(j)*x2(j-2); % Calculate the filter output.
    e2(j) = d2(j) + eta(j) - y2(j); % Calculate the error signal.
    
    % Calculate the weights (i.e. coefficients).
    w0(j+1) = w0(j) + 2*mu*x2(j)*e2(j);
    w1(j+1) = w1(j) + 2*mu*x2(j-1)*e2(j);
    w2(j+1) = w2(j) + 2*mu*x2(j-2)*e2(j);

end

display('Q4: The Adaptive Wiener Filter Weights for 1000 Number of Samples:');
display([w0(n2+1),w1(n2+1),w2(n2+1)]);

figure();
plot(1:n2,w0(1:n2),'^:');
hold on;
plot(1:n2,w1(1:n2),'*:');
plot(1:n2,w2(1:n2),'o:');
hold off;
legend('w0(n)','w1(n)','w2(n)');
xlabel('Iteration Index','FontSize',16);
ylabel('Coefficient (Weight) Value','FontSize',14);
title('Q4: Analysis of the LMS-Based Adaptive Wiener Filter for Non-Zero Input Signal','FontSize',12);

figure();
plot(1:n2,eta(1:n2),'^:');
hold on;
plot(1:n2,e2(1:n2),'*:');
hold off;
legend('\eta','e');
xlabel('Iteration Index','FontSize',16);
ylabel('Amplitude','FontSize',14);
title('Q4: Analysis of the LMS-Based Adaptive Wiener Filter for Non-Zero Input Signal','FontSize',12);
