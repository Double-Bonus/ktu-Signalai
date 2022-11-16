%%
clc, clear, close all;
set(0,'defaultAxesFontName','TimesNewRoman')
set(0,'defaultAxesFontSize',10)
   
load("signalai/EKG_13")

f_SL_Hz = 1; % projektuojamo zemu dazniu filtro slopinimo juostos ribinis daznis
f_pr_Hz = 0.4; % pralaidumo juostos ribinis daznis
f_d_Hz = 500;

time_s = 11;             
time_n =(0:(time_s*f_d_Hz)-1)*1/f_d_Hz;

% Solve equation to find M
syms x
eqn = (f_SL_Hz^2)*x^3 - (f_pr_Hz^2 -((f_SL_Hz+f_pr_Hz)^2))*x^2 + (2*f_d_Hz*(f_SL_Hz+f_pr_Hz))*x -f_d_Hz^2 == 0;
M_opt = vpasolve(eqn, x);

M = 50; % Decimacijos koeficientas parenkamas toks, kad diskretizavimo daznio ir koeficiento M dalybos
        % rezultatas buti sveikas skaicius

F = (f_SL_Hz - f_pr_Hz)/f_SL_Hz;
  
D1_s = (2*M*(1 - sqrt( ((M*F) / (2-F)))));
D1_d = 2 - F*(M*1);

D1_opt = D1_s / D1_d;

D2_opt = M/D1_opt;

D1 = 10;
D2 = 5;
if D1*D2 ~= M
    disp("error");
end  

%% pries paduodant naudoiti ta pati filtra kai p1 uzduotyje
Fs_low = 500;  % Sampling Frequency

N_low = 90;  % Order
Fpass = 45;  % Passband Frequency
Fstop = 47;  % Stopband Frequency
Wpass = 1;   % Passband Weight
Wstop = 10;   % Stopband Weight

% Calculate the coefficients using the FIRLS function.
[b_low, a_low] = firls(N_low, [0 Fpass Fstop Fs_low/2]/(Fs_low/2), [1 1 0 0], [Wpass Wstop]);
ekg = filter(b_low, a_low, ekg);
freqz(b_low, 1, 15000, f_d_Hz);


% Decimation
if 0 % using matlab API
    ekg_1 = decimate(ekg, D1, 'fir');
    ekg_2 = decimate(ekg_1, D2, 'fir');
    current_FD = f_d_Hz / (D1 * D2)
else
    b_safety_H1 = fir1(45,(24.9/(f_d_Hz/2))); % Naudojamas panasus filtras i decimate (filtro eile = 30; fp() < (500 / 2*10) 
    [b_H1_fqH, b_H1_fqW] =  freqz(b_safety_H1, 1, 15000, f_d_Hz);
    ekg_1 = filter(b_safety_H1, 1, ekg);
    ekg_1 = downsample(ekg_1, D1);
    current_FD = f_d_Hz / D1;
    
    b_safety_H2 = fir1(60,(3.5/(current_FD/2))); % Naudojamas panasus filtras i decimate (filtro eile = 30; fp() < (50 / 2*5) 
    [b_H2_fqH, b_H2_fqW] =  freqz(b_safety_H2, 1, 15000, current_FD);
    ekg_2 = filter(b_safety_H2, 1, ekg_1);
    ekg_2 = downsample(ekg_2, D2);
    current_FD = current_FD / D2;
end

% Low-pass filtras
b_lowPass = fir1(40, (0.6/(current_FD/2)));
[B_H_fqH, B_H_fzW] =  freqz(b_lowPass, 1, 15000, current_FD);
ekg_3 = filter(b_lowPass, 1, ekg_2);


if 0 % usisng MATLAB API
    ekg_4 = interp(ekg_3, D2);
    ekg_5 = interp(ekg_4, D1); % dreifas
else % make filter ourself
    ekg_4 = upsample(ekg_3, D2);
    ekg_4 = filter(b_safety_H2, 1, ekg_4);
    ekg_4 = ekg_4 * D2;
    
    ekg_5 = upsample(ekg_4, D1);
    ekg_5 = filter(b_safety_H1, 1, ekg_5);
    ekg_5 = ekg_5 * D1;
    
end
ekg_6 = ekg - ekg_5;

% check this API:
% grpdelay()


figure() % Plot filters magnitute responces
subplot(311); 
plot(b_H1_fqW, 20*log10(abs(b_H1_fqH)))
xlabel('f, Hz'); ylabel('A, dB'); ylim([-100 3]);
grid on; title('Amplitudes daznine charakteristika, H1 filtras');
subplot(312); 
plot(b_H2_fqW, 20*log10(abs(b_H2_fqH)))
xlabel('f, Hz'); ylabel('A, dB'); ylim([-100 3]);
grid on; title('Amplitudes daznine charakteristika, H2 filtras');
subplot(313); 
plot(B_H_fzW, 20*log10(abs(B_H_fqH)))
xlabel('f, Hz'); ylabel('A, dB'); ylim([-100 3]);
grid on; title('Amplitudes daznine charakteristika, H filtras');
saveas(gca,"outFigs/extra-H1-3-mag.jpg");


figure() % Plot EKG signals
subplot(311); 
plot(time_n, ekg); 
title('Laiko sritis: pradinis EKG');
xlabel('t, s'); ylabel('A, mV'); grid on;


subplot(312); 
plot(time_n, ekg_5); 
title('Laiko sritis: Dreifas');
xlabel('t, s'); ylabel('A, mV'); grid on;


subplot(313); 
plot(time_n, ekg_6); 
title('Laiko sritis: po multirat- RIR');
xlabel('t, s'); ylabel('A, mV'); grid on;


% delay of filter
velinimas = length(b_safety_H1) + length(b_safety_H2)*D1 + length(b_lowPass)*D1*D2/2


ekg_withZeros = [zeros(1, velinimas), ekg];
dreif_withZeros = [ekg_5, zeros(1, velinimas)];


ekg_noDreif = ekg_withZeros - dreif_withZeros;
figure()
subplot(211); 
plot(time_n, ekg_withZeros(velinimas+1:end)); 
xlabel('t, s'); ylabel('A, mV');
grid on; title('Laiko sritis: padinis EKG');

subplot(212); 
plot(time_n, ekg_noDreif(velinimas+1:end)); 
xlabel('t, s'); ylabel('A, mV'); 
grid on; title('Laiko sritis: EKG po filtro');
saveas(gca,"outFigs/extra-EkgAmpl.jpg");

% figure()
% plot(time_n, dreif_withZeros(velinimas+1:end)); 
% title('Laiko sritis: dreifas');
% xlabel('t, s'); ylabel('A, mV'); grid on;



figure()
subplot(211)
[x1, y1] = f_getFreqOfSignal(ekg_withZeros(velinimas+1:end), f_d_Hz);
plot(x1, y1); ylim([-100 0]);
title('Dazniu sritis: padinis EKG');  grid on;
xlim([0 f_d_Hz/2])

subplot(212); 
[x2, y2] = f_getFreqOfSignal(ekg_noDreif(velinimas+1:end), f_d_Hz);
plot(x2, y2); 
title('Dazniu sritis: EKG po filtro');  grid on;
xlim([0 f_d_Hz/2]); ylim([-100 0]);
saveas(gca,"outFigs/extra-Spekkt.jpg");



function [x_fq, y_fq] = f_getFreqOfSignal(sig, fd)
    nfft = length(sig);
    sig_fq = abs(fft(sig))/nfft;    
    sig_fq = 20*log10(sig_fq/max(sig_fq));   
    k = 0:1:nfft-1;
    
    x_fq = k*fd/nfft;
    y_fq = sig_fq;
end
