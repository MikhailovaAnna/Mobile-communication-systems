clc
clear all
close all

k = 8;
mes_book = zeros(2^k, k);
for i = 0 : 2^k - 1
    mes_book(i + 1,:) = de2bi(i,k);
end

check_bits = [9, 13, 15, 16];
length_cod = k + length(check_bits) + 4; % 4 проверочных символа и 4 нуля
code_book = zeros(2^k, length_cod);
for i = 0: 2^k - 1
%    добавляем проверочные символы
      q = 1;
      w = 1;
      for j = 5: length_cod
          if j == check_bits(1, q)
              code_book(i + 1, j) = 2;
              q = q + 1;
          else
              code_book(i + 1, j) = mes_book(i + 1, w);
              w = w + 1;
          end
      end
%     считаем ксор тех позиций в двоичной системе, где единички
    contr_sum = zeros(1, length(check_bits));
    for j = 4:length_cod - 1
        if code_book(i + 1, j + 1) == 1
            pos = length_cod - j;
            z = (num2str(str2double(dec2bin(pos))) - '0');
            z = [zeros(1, length(check_bits) - length(z)), z];
            contr_sum = contr_sum + z;
        end
    end
    contr_sum = mod(contr_sum, 2);
    for j = 0: length(check_bits) - 1
        code_book(i + 1, check_bits(1, j + 1)) = contr_sum(1, j + 1);
    end
end

% модулятор
mod_book = zeros(2^k, length_cod - 4);
for i = 0:2^k - 1
    mod_book(i + 1,:) = code_book(i + 1, 5:end)*(-2) + 1;
end

num_mes = 10000;
SNRdB = -15:3:0;
SNR = 10.^(SNRdB/10);
sigma = sqrt((1 ./ (2 * SNR)));
Pb_hard = zeros(1, length(SNRdB));
Pb_theor = erfc(sqrt(2 * 10 .^ (SNRdB / 10)) ./ sqrt(2)) * 0.5;
for n = 1:length(SNRdB)
    corr_mes = 0;
    while corr_mes < num_mes
        i = randi([0, 2^k - 1]);
        mes = mes_book(i + 1,:);
        % кодириуем
        cod = code_book(i + 1,:);
        % модулируем
        mod_mes = mod_book(i + 1,:);
        % добавляем шум
        noise = sigma(n)*randn(1, length(mod_mes));
        r_noise = mod_mes + noise;
        % демодуляция
        c_demod = r_noise < 0;
        c_demod = [[0 0 0 0], c_demod];
        % декодирование
        decod = zeros(1, length_cod);
        %добавляем проверочные символы, пишем в нужные места 2
        q = 1;
        w = 1;
        for j = 0: length_cod - 1
            if j + 1 ~= check_bits(1, q)
                decod(1, j + 1) = c_demod(1, w);
                w = w + 1;
            else
                q = q + 1;
            end
        end
        %пересчитываем контрольную сумму
        contr_sum = zeros(1, length(check_bits));
        for j = 0:length(decod) - 1
            if decod(j + 1) == 1
                z = (num2str(str2double(dec2bin(length(decod) - j))) - '0');
                z = [zeros(1, length(check_bits) - length(z)), z];
                contr_sum = contr_sum + z;
            end
        end
        contr_sum = mod(contr_sum, 2);
        % записываем контрольную сумму
        for j = 0:length(check_bits) - 1
            decod(1, check_bits(j + 1)) = contr_sum(j + 1);
        end
        % сравниваем контрольную сумму до и после
        pos = 0;
        flag = 0;
        for j = 0:length(check_bits) - 1
            if c_demod(1, check_bits(j + 1)) ~= decod(1, check_bits(j + 1))
                pos = pos + length(cod) - check_bits(j + 1) + 1;
            end
        end
        % ищем позицию, где была ошибка и исправляем ее
        if pos ~= 0 %pos > -1 && flag == 1
            c_demod(1, pos) = mod(decod(1, pos) + 1, 2);
        end
        % убираем проверочные биты
        new_decod_hard = zeros(1, length(decod) - length(check_bits));
        w = 1;
        h = 1;
        for j = 0: length(decod) - 1
            if j + 1 ~= check_bits(1, w)
                new_decod_hard(1, h) = c_demod(1, j + 1);
                h = h + 1;
            else
                w = w + 1;
            end 
        end
        new_decod_hard = new_decod_hard(1, 5:end);
        b = sum(mes ~= new_decod_hard);
        Pb_hard(n) = Pb_hard(n)+ b/length(mes);
        corr_mes = corr_mes + 1;
    end
    Pb_hard(n) = Pb_hard(n)/num_mes;
end

Pb_soft = zeros(1, length(SNRdB));
for n = 1:length(SNRdB)
    corr_mes = 0;
    while corr_mes < num_mes
        i = randi([0, 2^k - 1]);
        mes = mes_book(i + 1,:);
        % кодириуем
        cod = code_book(i + 1,:);
        % модулируем
        mod_mes = mod_book(i + 1,:);
        % добавляем шум
        noise = sigma(n)*randn(1, length(mod_mes));
        r_noise = mod_mes + noise;
        [~, I] = min(pdist2(mod_book(1:end,:), r_noise));
        soft_res = mes_book(I,:);
        b = sum(mes ~= soft_res);
        Pb_soft(n) = Pb_soft(n)+ b/length(mes);
        corr_mes = corr_mes + 1;
    end
    T(n) = k * num_mes/(length_cod*corr_mes);
    Pb_soft(n) = Pb_soft(n)/num_mes;
end

figure(1);
semilogy(SNRdB, Pb_hard, 'r', SNRdB, Pb_soft, 'g', SNRdB, Pb_theor, 'ro', 'LineWidth', 2)
grid on
legend('Вероятность ошибки на бит жесткого декодера', 'Вероятность ошибки на бит мягкого декодера', 'Теоретическая вероятность ошибки на бит');
title('Вероятности ошибки на бит от SNRdB');
xlabel('SNRdB');
ylabel('Вероятность ошибки на бит');