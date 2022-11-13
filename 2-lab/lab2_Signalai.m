clc
clear all
close all
   
load("signalai/EKG_13")
fd = 500; %add _Hz            % Diskretizavimo freq
ts = 11;             
t =(0:(ts*fd)-1)*1/fd;

%Noramlius font size ir db iki -80
figure_cnt = 1;

%% 3.3.1 RIR filtras

% Filtro koeficientai
b = [0.0181, 0.0319, -0.0084, -0.0803, -0.0449, 0.1709, 0.4008, 0.4008, ...
    0.1709, -0.0449, -0.0803, -0.0084, 0.0319, 0.0181]; 
a = zeros(1,length(b)); 
a(1) = 1;

ekg_filtered = filter(b, a, ekg);
figure(figure_cnt); figure_cnt=figure_cnt+1;
plot(t, ekg_filtered)

%% 3.3.2
figure(figure_cnt); figure_cnt=figure_cnt+1;
impz(b, a);

%% 3.3.3
n = 15000;
figure(figure_cnt); figure_cnt=figure_cnt+1;
freqz(b, a, n, fd)
datacursormode('on')


%% 3.3.4
figure(figure_cnt); figure_cnt=figure_cnt+1;
zplane(b,a)


%% 3.3.5
%% plot for time
figure(figure_cnt); figure_cnt=figure_cnt+1;
subplot(211); 
plot(t, ekg); 
title('Laiko sritis prie? RIR');
% xlabel('t, s'); ylabel('A, mV'); xlim(T_DISP); grid on;
subplot(212); 
plot(t, ekg_filtered);
title('Laiko sritis po RIR');
% xlabel('t, s'); ylabel('A, mV'); xlim(T_DISP); grid on;


%% plot freq: plot signal in freqency domain
% (Same as in 1 Lab) - use function
figure(figure_cnt); figure_cnt=figure_cnt+1;

nfft = length(ekg);
ekg_freq = abs(fft(ekg))/nfft;    
ekg_freq = 20*log10(ekg_freq/max(ekg_freq));   
k = 0:1:nfft-1;                 
f = k*fd/nfft;                  


subplot(211)
plot(f, ekg_freq, 'k');
xlabel('f, Hz')
ylabel('Sa, dB')
title('Signalo spektras pries apdorojima')
xlim([0 fd/2])

ekgFilter_freq = abs(fft(ekg_filtered))/nfft;    
ekgFilter_freq = 20*log10(ekgFilter_freq/max(ekgFilter_freq));  
subplot(212)
plot(f, ekgFilter_freq, 'k');
xlabel('f, Hz')
ylabel('Sa, dB')
title('Signalo spektras po RIR apdorojimo')
xlim([0 fd/2])


%% 3.4 NIR Filtras

fDelta_Hz = 1.5;

% constans:
f0_Hz = 50;
K0 = 1;
L = 3;
S = mag2db(150); % double check

%3.4.2
K = K0 * 10^(-S/20);
Kr = K0 * 10^(-L/20);

%3.4.3
N = fd/f0_Hz;

%3.4.4
beta = sqrt((Kr^2 - K0^2)/(K^2 - Kr^2)) * tan((N*pi*fDelta_Hz) / (2*fd));

k1 = (K0 + (K*beta)) / (1 + beta);
k2 = (K0 - (K*beta)) / (1 + beta);
k3 = (1 - beta) / (1 + beta);


bNIR = [k1, zeros(1,N-1), -k2];
aNIR = [1, zeros(1,N-1), -k3];

%3.4.6
ekg_afterNIR = filter(bNIR, aNIR, ekg_filtered);

% analyzsis
figure(figure_cnt); figure_cnt=figure_cnt+1;

%3.4.7
impz(bNIR, aNIR);
grid on;
xlim([0 100]);
title('Impusline charakteristika');

%3.4.8
figure(figure_cnt); figure_cnt=figure_cnt+1;
freqz(bNIR, aNIR, n, fd)
title('NIR analize');
ylim([-50 5])

% 3.4.9
figure(figure_cnt); figure_cnt=figure_cnt+1;
zplane(bNIR, aNIR)
title('NIR zplane');


% 3.4.10
    % plot for time
figure(figure_cnt); figure_cnt=figure_cnt+1;
subplot(211); 
plot(t, ekg); 
title('Laiko sritis neapdoroto ekg');
% xlabel('t, s'); ylabel('A, mV'); xlim(T_DISP); grid on;
subplot(212); 
plot(t, ekg_afterNIR);
title('Laiko sritis po NIR');
% xlabel('t, s'); ylabel('A, mV'); xlim(T_DISP); grid on;



% plot freq: plot signal in freqency domain
% (Same as in 1 Lab) - use function
figure(figure_cnt); figure_cnt=figure_cnt+1;

nfft = length(ekg);
ekg_freq = abs(fft(ekg))/nfft;    
ekg_freq = 20*log10(ekg_freq/max(ekg_freq));   
k = 0:1:nfft-1;                 
f = k*fd/nfft;                  


subplot(211)
plot(f, ekg_freq, 'k');
xlabel('f, Hz')
ylabel('Sa, dB')
title('Signalo spektras neapdoroto ekg')
xlim([0 fd/2])

ekg_afterNIR_fq = abs(fft(ekg_afterNIR))/nfft;    
ekg_afterNIR_fq = 20*log10(ekg_afterNIR_fq/max(ekg_afterNIR_fq));  
subplot(212)
plot(f, ekg_afterNIR_fq, 'k');
xlabel('f, Hz')
ylabel('Sa, dB')
title('Signalo spektras po NIR apdorojimo')
xlim([0 fd/2])

figure(figure_cnt); figure_cnt=figure_cnt+1;
f_plotAllSignalsTime(t, ekg, ekg_filtered, ekg_afterNIR)


figure(figure_cnt); figure_cnt=figure_cnt+1;
f_plotAllSignalsFreq(ekg, ekg_filtered, ekg_afterNIR, fd)

function f_plotAllSignalsTime(t, initial, RIR, NIR)

    subplot(311)
    plot(t, initial); 
    title('Laiko sritis neapdoroto ekg');
    subplot(312); 
    plot(t, RIR);
    title('Laiko sritis po RIR');

    subplot(313); 
    plot(t, NIR);
    title('Laiko sritis po RIR ir NIR');
    
end

function f_plotAllSignalsFreq(initial, RIR, NIR, fd)

    subplot(311)
    [x1, y1] = getFreqOfSignal(initial, fd);
    plot(x1, y1); 
    title('Dazniu sritis neapdoroto ekg');
    xlim([0 fd/2])
    
    subplot(312); 
    [x2, y2] = getFreqOfSignal(RIR, fd);
    plot(x2, y2); 
    title('Dazniu sritis po RIR');
    xlim([0 fd/2])

    subplot(313); 
    [x3, y3] = getFreqOfSignal(NIR, fd);
    plot(x3, y3); 
    title('Dazniu sritis po RIR ir NIR');
    xlim([0 fd/2])
    
end

function [x_fq, y_fq] = getFreqOfSignal(sig, fd)

    nfft = length(sig);
    sig_fq = abs(fft(sig))/nfft;    
    sig_fq = 20*log10(sig_fq/max(sig_fq));   
    k = 0:1:nfft-1;
    
    x_fq = k*fd/nfft;
    y_fq = sig_fq;

end
