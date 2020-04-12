clc
clear all
close all
nfig = 1;

% передаем сообщение из 32 нулей
k = 32;
Message = zeros(1, k);

% добавляем контрольную сумму
G = [1 0 0 1 1 1 1 0 1 0 1 1 0 0 1 0 1];
size_g = length(G);
length_C = k + size_g - 1;
CRC = [Message, zeros(1, size_g - 1)];
[b, r] = deconv(CRC, G);
r_t = mod(r, 2);
CRC = CRC + r_t; % сообщение с контрольной суммой

tblen1 = 4;
trellis1 = poly2trellis(tblen1, [13 17]); % получили решетку для скорости 1/2
Cod1 = convenc(CRC, trellis1); % закодировали
Model1 = Cod1*(-2) + 1; % замодулировали

num_mes = 200000; % число сообщений, которые мы хотим передать
SNRdB = -12:2:0;
SNR = 10.^(SNRdB/10);
sigma = sqrt((1 ./ (2 * SNR)));
T1 = zeros(1, length(SNRdB)); % пропускная способность HARQ1
Pdec1 = zeros(1, length(SNRdB)); % вероятность ошибки декодирования HARQ1
for n = 1:length(SNRdB)
    n_exp = 0; % число передач
    corr_mes = 0; % число правильно переданных сообщений
    bufer = zeros(1, 96); 
    while corr_mes < num_mes
        rnd_x = randn();
        rnd_y = randn();
        H = (rnd_x*rnd_x + rnd_y*rnd_y)^0.5; % замирание
        noise = sigma(n)*randn(1, length(Model1));
        r = Model1*H + noise; % выход канала
        
        % декодирование
        decod = vitdec(r, trellis1, tblen1, 'trunc', 'unquant'); % мягкое декодирование
  
        [p, s] = deconv(decod, G);
        s = mod(s, 2);
        q = sum(s);
        b = sum(CRC ~= decod);
        if q ~= 0 % синдром не равен нулю, пробуем декодировать данные из буфера
            bufer = bufer + r; % аккумулируем сообщения в буфере
            new_decod = vitdec(bufer, trellis1, tblen1, 'trunc', 'unquant'); % декодируем
            [p, s] = deconv(new_decod, G);
            s = mod(s, 2);
            q = sum(s);
            b = sum(CRC ~= new_decod);
            if q == 0 % если синдром равен нулю - сообщение успешно декодировалось
                corr_mes = corr_mes + 1;
                bufer = zeros(1, 96); % очищаем буфер
            end
        else
            corr_mes = corr_mes + 1;
        end
        Pdec1(n) = Pdec1(n) + ((q == 0)&(b > 0)); % если синдром равен нулю, а число неправильно декодированных бит больше нуля, то появляется ошибка декодирования
        n_exp = n_exp + 1;
    end
    T1(n) = k * num_mes/(length(Cod1)*n_exp); % считаем пропускную способность
    Pdec1(n) = Pdec1(n)/n_exp;
end

tblen2 = [5 4];
trellis2_7 = poly2trellis(tblen2, [23 35 0 32 07 16 14; 0 5 13 07 13 05 10]); % решетка для кода со скоростью 2/7
trellis2_5 = poly2trellis(tblen2, [23 35 0 32 07; 0 5 13 07 13]); % решетка для кода со скоростью 2/5
trellis2_3 = poly2trellis(tblen2, [23 35 0; 0 5 13]); % решетка для кода со скоростью 2/3
Cod2_7 = convenc(CRC, trellis2_7); % сначала кодируем сообщение с самой низкой скоростью

% здесь по сути пакет1, пакет2, пакет3 нужны для замодулированных
% сообщений, они будут формироваться в процессе передачи заново
P1 = zeros(1, 3*(length(Cod2_7)/7)); % создаем первый пакет, который будем пытаться передать с максимальной скоростью 2/3
j = 1;
count = 0;
for i = 1: length(P1)
    if count < 3
        P1(1, i) = Cod2_7(1, j);
        j = j + 1;
        count = count + 1;
    else
        j = j + 4;
        count = 0;
    end
end

P2 = zeros(1, 2*(length(Cod2_7)/7)); % добавочный пакет данных, для передачи со скоростью 2/5
j = 1;
count = 0;
for i = 1: length(P2)
    if count < 2
        P2(1, i) = Cod2_7(1, j + 3);
        j = j + 1;
        count = count + 1;
    else
        j = j + 5;
        count = 0;
    end
end

P3 = zeros(1, 2*(length(Cod2_7)/7)); % пакет 3, оставшиеся биты для передачи со скоростью 2/7
ModelP1 = P1*(-2) + 1;
ModelP2 = P2*(-2) + 1;
ModelP3 = P3*(-2) + 1;

T2 = zeros(1, length(SNRdB)); % пропускная способность HARQ2
Pdec2 = zeros(1, length(SNRdB)); % вероятность ошибки декодирования HARQ2
av_len = 0; % средняя длина закодированного сообщения, которое декодировалось с синдромом равным 0
for n = 1:length(SNRdB)
    n_exp = 0;
    corr_mes = 0;
    len = 0; % длина закодированного сообщения, которое декодировалось
    while corr_mes < num_mes
        q = -1;
        rnd_x = randn();
        rnd_y = randn();
        H = (rnd_x*rnd_x + rnd_y*rnd_y)^0.5; % замирание
        noise2_3 = sigma(n)*randn(1, length(ModelP1));
        r2_3 = ModelP1*H + noise2_3;
        % декодирование
        decod2_3 = vitdec(r2_3, trellis2_3, 5, 'trunc', 'unquant');
        [p, s3] = deconv(decod2_3, G);
        s3 = mod(s3, 2);
        q3 = sum(s3);
        b = sum(CRC ~= decod2_3);
        if q3 ~= 0 % если не смогло успешно передаться со скоростью 2/3, то передаем со скоростью 2/5
            noiseP2 = sigma(n)*randn(1, length(ModelP2));
            rP2 = ModelP2*H + noiseP2;
            r2_5 = zeros(1, length(r2_3)+length(rP2)); % формируем пакет, который содержит в себе предыдущий пакет, который не передался со скоростью 2/3, и новый, который дошел
            count = 0;
            j = 1;
            t = 1;
            for i = 1:length(r2_5)
               if count < 3
                   r2_5(1, i) = r2_3(1, j);
                   j = j + 1;
                   count  = count + 1;
               elseif count < 5
                   r2_5(1, i) = rP2(1, t);
                   t = t + 1;
                   if count + 1 < 5
                       count = count + 1;
                   else
                       count = 0;
                   end
               end
            end
            decod2_5 = vitdec(r2_5, trellis2_5, 5, 'trunc', 'unquant'); % пытаемся декодировать
            [p, s5] = deconv(decod2_5, G);
            s5 = mod(s5, 2);
            q5 = sum(s5);
            b = sum(CRC ~= decod2_5);
            if q5 ~= 0 % если не передалось со скоростью 2/5, пытаемся передать со скоростью 2/7
                noiseP3 = sigma(n)*randn(1, length(ModelP3));
                rP3 = ModelP3*H + noiseP3;
                r2_7 = zeros(1, length(Cod2_7)); % содержит в себе пакет для скорости 2/5 и новые дошедшие данные
                count = 0;
                j = 1;
                t = 1;
                for i = 1:length(r2_7)
                   if count < 5
                       r2_7(1, i) = r2_5(1, j);
                       j = j + 1;
                       count  = count + 1;
                   elseif count < 7
                       r2_7(1, i) = rP3(1, t);
                       t = t + 1;
                       if count + 1 < 7
                           count = count + 1;
                       else
                           count = 0;
                       end
                   end
                end
                decod2_7 = vitdec(r2_7, trellis2_7, 5, 'trunc', 'unquant'); % пытаемся декодировать
                [p, s7] = deconv(decod2_7, G);
                s7 = mod(s7, 2);
                q7 = sum(s7);
                b = sum(CRC ~= decod2_7);
                if q7 == 0
                    corr_mes = corr_mes + 1;
                    q = q7;
                    len = length(r2_7);
                end
            else 
                corr_mes = corr_mes + 1;
                q = q5;
                len = length(r2_5);
            end
        else
            corr_mes = corr_mes + 1;
            q = q3;
            len = length(r2_3);
        end
        Pdec2(n) = Pdec2(n) + ((q == 0)&(b > 0)); % считаем вероятность ошибки декодирования
        n_exp = n_exp + 1;
        av_len = av_len + len; 
    end
    av_len = av_len / num_mes;
    T2(n) = k * num_mes/(av_len*n_exp); % считаем пропускную способность
    Pdec2(n) = Pdec2(n)/n_exp;
end

figure(1);
plot(SNRdB, T1, 'r', SNRdB, T2, 'b', 'LineWidth', 2)
grid on
legend('Пропускная способность T_1', 'Пропускная способность T_2');
title('Пропускная способность от SNRdB');
xlabel('SNRdB');
ylabel('Пропускная способность, T');

figure(2);
semilogy(SNRdB, Pdec1, 'r', SNRdB, Pdec2, 'b', 'LineWidth', 2)
grid on
legend('Вероятность ошибки декодирования Type1', 'Вероятность ошибки декодирования Type2');
title('Вероятности ошибки декодирования от SNRdB');
xlabel('SNRdB');
ylabel('Вероятность ошибки декодирования');
