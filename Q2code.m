clear all, close all;
clc
%% part 2 (c) plot the real coefficients (call this the original function)
a = pi/3;
b = 2*pi/3;
n = [-32:-1 1:32];
w = -pi:0.01:pi;
c0 = (-2*a+2*b)/(2*pi);
cn = (-sin(a*n)+sin(b*n))./(n*pi);
s = c0 + real(sum(diag(cn)*exp(i*transpose(n)*w))); %diag(cn) so that 1st row multiplied by c1, 2nd row multiplied by c2, etc... 
%and transpose(n) since matrix multiplication operations need to be obeyed
figure(1);
plot(w, s)
xlim([-pi, pi])
title('s(t) vs t')
%% part (d) plot log of original
log_original = 20*log10(abs(c0 + sum(diag(cn)*exp(i*transpose(n)*w))));
figure(2);
plot(w, log_original)
xlim([-pi, pi])
title('Log(s(t)) vs t')
%% part (e) barlet with original
NN = 32;
w10 = 1;
w1n = 1 - abs(n/NN);
sw1 = c0*w10 + real(sum(diag(cn)*diag(w1n)*exp(i*transpose(n)*w))); 
figure(3);
plot(w, s), hold on
plot(w, sw1)
xlim([-pi, pi])
title('Bartlet window')
%% part(e) blackman with original
w20 = 0.42 + 0.5 + 0.08;
w2n = 0.42 + 0.5*cos(pi*n/NN) + 0.08*cos(2*pi*n/NN);
sw2 = c0*w20 + real(sum(diag(cn)*diag(w2n)*exp(i*transpose(n)*w)));
figure(4);
plot(w, s), hold on
plot(w, sw2)
xlim([-pi, pi])
title('Blackman window')
%% part (e) Hanning with original
w30 = 0.5 + 0.5;
w3n = 0.5 + 0.5*cos(pi*n/NN);
sw3 = c0*w30 + real(sum(diag(cn)*diag(w3n)*exp(i*transpose(n)*w)));
figure(5);
plot(w, s), hold on
plot(w, sw3)
xlim([-pi, pi])
title('Hanning window')
%% part (e) plot log of 4 functions
log_Barlet = 20*log10(c0*w10 + sum(diag(cn)*diag(w1n)*exp(i*transpose(n)*w)));
log_Blackman = 20*log10(c0*w20 + sum(diag(cn)*diag(w2n)*exp(i*transpose(n)*w)));
log_Hanning = 20*log10(c0*w30 + sum(diag(cn)*diag(w3n)*exp(i*transpose(n)*w)));
figure(6);
plot(w, log_Barlet), hold on
plot(w, log_Blackman)
plot(w, log_Hanning)
plot(w, log_original)
xlim([-pi, pi])
title('Log(s(t)) vs t')
legend({'Bartlet', 'Blackman', 'Hanning', 'Original'})
%% part (g)
fs = 8820;
a = 2*pi*600/fs;
b = 2*pi*1662/fs;
NN = 17; %since matlab index starts at 1, we want coefficients h[0] to h[16]
h = zeros(NN+1,1);
h(1) = w20*(-2*a+2*b)/(2*pi);
for j=1:NN
  w2 = 0.42 + 0.5*cos(pi*j/NN) + 0.08*cos(2*pi*j/NN);
  h(j+1) = w2*(-sin(a*j)+sin(b*j))/(j*pi);
end

disp('h[j] =')
disp(h)
