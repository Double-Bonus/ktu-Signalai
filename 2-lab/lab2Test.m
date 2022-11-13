clc, clear, close all;

fd = 500;

t_sec = 15;
t = 0: 1/fd: t_sec-1/fd;

% Test with square
t_sq = length(t)/t_sec/10; % 100ms
y_sq = zeros(1, fd);
for i = 1:t_sq
    y_sq(i) = 3; 
end
y_sq = repmat(y_sq, 1, t_sec);
plot(t, y_sq)


[y_multi, vel] = f_getMultirate(y_sq, fd)

figure()
plot(t, y_multi(vel+1:end))


yNir = f_getNirFilter(y_sq, fd);
figure()
plot(t, yNir)


%% Test wih triangle

fd = 500;
t_sec = 15;
t = 0: 1/fd: t_sec-1/fd;

t_triangle_sec = 0.202;
t_temp = 0: 1/fd: t_triangle_sec-1/fd;

for ii = 1:length(t_temp)
    if t_temp(ii) <= 0.1
        y_triangle(ii) = 15*t_temp(ii);   
    else
        y_triangle(ii) = -15*t_temp(ii) +3;
    end
end
figure()
plot(t_temp, y_triangle)

y_triangle = [y_triangle, zeros(1, fd - length(t_temp))];
y_triangle = repmat(y_triangle, 1, t_sec);
figure()
plot(t, y_triangle)

% get filters
[y_multi, vel] = f_getMultirate(y_triangle, fd)

figure()
plot(t, y_multi(vel+1:end))
title('Multirate filtras')

yNir = f_getNirFilter(y_triangle, fd);
figure()
plot(t, yNir)
title('NIR suku tipo filtras filtras')



function [out_y, velinimas] = f_getMultirate(sig_y, f_d_Hz)
    D1 = 10
    D2 = 5

    % Decimation
    b_safety_H1 = fir1(30,(20/(f_d_Hz/2))); % Naudojamas panasus filtras i decimate (filtro eile = 30; fp(20) < (500 / 2*10) 
    ekg_1 = filter(b_safety_H1, 1, sig_y);
    ekg_1 = downsample(ekg_1, D1);
    current_FD = f_d_Hz / D1;

    b_safety_H2 = fir1(30,(4/(current_FD/2))); % Naudojamas panasus filtras i decimate (filtro eile = 30; fp(4) < (50 / 2*5) 
    ekg_2 = filter(b_safety_H2, 1, ekg_1);
    ekg_2 = downsample(ekg_2, D2);
    current_FD = current_FD / D2;

    % Low-pass filtras
    b_lowPass = fir1(40, (0.67/(current_FD/2)));
    ekg_3 = filter(b_lowPass, 1, ekg_2);

    ekg_4 = interp(ekg_3, D2);
    ekg_5 = interp(ekg_4, D1); % dreifo linija
    ekg_6 = sig_y - ekg_5;

    velinimas = length(b_safety_H1) + length(b_safety_H2)*D1 + length(b_lowPass)*D1*D2/2

    ekg_withZeros = [zeros(1, velinimas), sig_y];
    dreif_withZeros = [ekg_5, zeros(1, velinimas)];

    out_y = ekg_withZeros - dreif_withZeros;
end

function out_y = f_getNirFilter(sig_y, fd)
    fDelta_Hz = 1.5;
    f0_Hz = 50;
    K0 = 1;
    L = 3;
    S = mag2db(150);

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
    out_y = filter(bNIR, aNIR, sig_y);




end  