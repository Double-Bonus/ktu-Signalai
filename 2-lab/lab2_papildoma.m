clc
clear all
close all
   
load("signalai/EKG_13")


f_SL_Hz = 1; % projektuojamo ?em?j? da?ni? filtro slopinimo juostos ribinis da?ni
f_pr_Hz = 0.4; % pralaidumo juostos ribinis da?nis
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

%pries paduodant naudoiti ta pati filtra kai p1 uzduotyje??

%%

% low-pass filter before analysis
% firls


% Decimation
ekg_1 = decimate(ekg, D1, 'fir');
ekg_2 = decimate(ekg_1, D2, 'fir');

% Low-pass filtras
current_FD = f_d_Hz / (D1 * D2)
b_lowPass = fir1(51,(0.67/(current_FD/2)),'low');
[B_H_freqz,B_Hw_freqz] =  freqz(b_lowPass, 1, 15000, f_d_Hz);
freqz(b_lowPass, 1, 15000, current_FD);
ekg_3 = filter(b_lowPass, 1, ekg_2);


% LEVEL 1 UPSAMPLE
ekg_4 = interp(ekg_3, D2);
ekg_5 = interp(ekg_4, D1); % dreifo linija

ekg_6 = ekg - ekg_5;

% pacioje pradziojevelinnimas????
% grpdelay()


figure(2)
subplot(311); 
plot(time_n, ekg); 
title('Laiko sritis prie? RIR');
xlabel('t, s'); ylabel('A, mV'); grid on;




subplot(312); 
plot(time_n, ekg_5); 
title('Laiko sritis po RIR');
xlabel('t, s'); ylabel('A, mV'); grid on;



subplot(313); 
plot(time_n, ekg_6); 
title('Laiko sritis po RIR');
xlabel('t, s'); ylabel('A, mV'); grid on;



% try with delay
velinimas = length(b_lowPass)*D1*D2/2 %TODO get legit value here


ekg_withZeros = [zeros(1, velinimas), ekg];
dreif_withZeros = [ekg_5, zeros(1, velinimas)];


ekg_noDreif = ekg_withZeros - dreif_withZeros;
figure(3)
subplot(311); 
plot(time_n, ekg_withZeros(velinimas+1:end)); 
title('Laiko sritis prie? RIR');
xlabel('t, s'); ylabel('A, mV'); grid on;


subplot(312); 
plot(time_n, ekg_noDreif(velinimas+1:end)); 
title('Laiko sritis po RIR');
xlabel('t, s'); ylabel('A, mV'); grid on;



subplot(313); 
plot(time_n, dreif_withZeros(velinimas+1:end)); 
title('Laiko sritis po RIR');
xlabel('t, s'); ylabel('A, mV'); grid on;





