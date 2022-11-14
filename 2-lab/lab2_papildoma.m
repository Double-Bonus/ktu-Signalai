clc, clear, close all;
   
load("signalai/EKG_13")

f_SL_Hz = 1; % projektuojamo iemeji dazni? filtro slopinimo juostos ribinis dazni
f_pr_Hz = 0.4; % pralaidumo juostos ribinis daznis
f_d_Hz = 500;

time_s = 11;             
time_n =(0:(time_s*f_d_Hz)-1)*1/f_d_Hz;


m = [f_SL_Hz^2 - f_pr_Hz^2; -((f_SL_Hz+f_pr_Hz)^2); 2*f_d_Hz*(f_SL_Hz+f_pr_Hz); -f_d_Hz^2]; % isprende lygti randame optimalia verte
M_opt = roots(m)

M = 50; % Decimacijos koeficientas parenkamas toks, kad diskretizavimo da?nio ir koeficiento M dalybos
  % rezultatas b?t? sveikas skai?ius


F = (f_SL_Hz - f_pr_Hz)/f_SL_Hz;
  
D1_s = (2*M*(1 - sqrt( ((M*F) / (2-F)))));
D1_d = 2 - F*(M*1);

D1_opt = D1_s / D1_d

D2_opt = M/D1_opt

D1 = 10
D2 = 5
if D1*D2 ~= M
    disp("error");
end  

%%
%pries paduodant naudoiti ta pati filtra kai p1 uzduotyje??
% low-pass filter before analysis
% firls

% Decimation
if 0 % using matlab API
    ekg_1 = decimate(ekg, D1, 'fir');
    ekg_2 = decimate(ekg_1, D2, 'fir');
    current_FD = f_d_Hz / (D1 * D2)
else
    b_safety_H1 = fir1(30,(20/(f_d_Hz/2))); % Naudojamas panasus filtras i decimate (filtro eile = 30; fp(20) < (500 / 2*10) 
    [b_H1_fqH, b_H1_fqW] =  freqz(b_safety_H1, 1, 15000, f_d_Hz);
    ekg_1 = filter(b_safety_H1, 1, ekg);
    ekg_1 = downsample(ekg_1, D1);
    figure()
    freqz(b_safety_H1, 1, 20000, f_d_Hz);  % should be atleast -50dB
    current_FD = f_d_Hz / D1;
    
    b_safety_H2 = fir1(30,(4/(current_FD/2))); % Naudojamas panasus filtras i decimate (filtro eile = 30; fp(4) < (50 / 2*5) 
    [b_H2_fqH, b_H2_fqW] =  freqz(b_safety_H2, 1, 15000, current_FD);
    ekg_2 = filter(b_safety_H2, 1, ekg_1);
    ekg_2 = downsample(ekg_2, D2);
    figure()
    freqz(b_safety_H2, 1, 20000, current_FD);  % should be atleast -50dB
    current_FD = current_FD / D2;
end

% Low-pass filtras
b_lowPass = fir1(40, (0.67/(current_FD/2)));
[B_H_freqz,B_Hw_freqz] =  freqz(b_lowPass, 1, 15000, current_FD);
figure()
freqz(b_safety_H2, 1, 20000, current_FD); % should be atleast -50dB
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


figure()
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
subplot(311); 
plot(time_n, ekg_withZeros(velinimas+1:end)); 
title('Laiko sritis: padinis EKG');
xlabel('t, s'); ylabel('A, mV'); grid on;


subplot(312); 
plot(time_n, ekg_noDreif(velinimas+1:end)); 
title('Laiko sritis: nuemus dreifa');
xlabel('t, s'); ylabel('A, mV'); grid on;



subplot(313); 
plot(time_n, dreif_withZeros(velinimas+1:end)); 
title('Laiko sritis: dreifas');
xlabel('t, s'); ylabel('A, mV'); grid on;