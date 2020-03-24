clc
clear all
close all
nfig = 1;

k = 8;
% все слова
M = zeros(2^k, k);
for i = 0 : 2^k - 1
    M(i + 1,:) = de2bi(i,k);
end

G = [1 0 0 1 1 1 1 0 1 0 1 1 0 0 1 0 1];
size_g = length(G);
length_C = k + size_g - 1;
C = zeros(2^k, length_C);
for i = 0: 2^k - 1
    C(i + 1,:) = [M(i + 1,:), zeros(1, size_g - 1)];
    [q, r] = deconv(C(i + 1,:), G);
    r_t = mod(r, 2);
    C(i + 1,:) = C(i + 1,:) + r_t;
end

weight_word = sum(C(2:end, :), 2);
d_min = min(weight_word);

check_bits = [16, 24, 28, 30, 31];
length_cod = length_C + 7; % 5 проверочных символов и 2 нуля
Cod = zeros(2^k, length_cod);
for i = 0: 2^k - 1
%    добавляем проверочные символы, пишем в нужные места 2
      q = 1;
      w = 1;
      for j = 3: length_cod
          if j == check_bits(1, q)
              Cod(i + 1, j) = 2;
              q = q + 1;
          else
              Cod(i + 1, j) = C(i + 1, w);
              w = w + 1;
          end
      end
%     считаем ксор тех позиций в двоичной системе, где единички
    c = zeros(1, 5);
    for j = 0:length_cod - 1
        if Cod(i + 1, j + 1) == 1
            pos = length_cod - j;
            z = (num2str(str2double(dec2bin(pos))) - '0');
            z = [zeros(1, 5 - length(z)), z];
            c = c + z;
        end
    end
    contr_sum = mod(c, 2);
    for j = 0: length(check_bits) - 1
        Cod(i + 1, check_bits(1, j + 1)) = contr_sum(1, j + 1);
    end
end

% модулятор
R = zeros(2^k, length_cod - 2);
for i = 0:2^k - 1
    R(i + 1,:) = Cod(i + 1, 3:length_cod)*(-2) + 1;
end

num_mes = 1000;
SNRdB = -8:2:0;
SNR = 10.^(SNRdB/10);
sigma = sqrt((1 ./ (2 * SNR)));
Pb = zeros(1, length(SNRdB));
Pb_theor = erfc(sqrt(2 * 10 .^ (SNRdB / 10)) ./ sqrt(2)) * 0.5;
Pdec = zeros(1, length(SNRdB));
Ped_th = zeros(1, length(SNRdB));
Ped_pr = zeros(1, length(SNRdB));
T = zeros(1, length(SNRdB));
for n = 1:length(SNRdB)
    n_exp = 0;
    corr_mes = 0;
    while corr_mes < num_mes
        q = 2; % количество единичек в синдроме
        i = randi([0, 2^k - 1]);
        per = 0;
        m = M(i + 1,:);
        % добавляем контрольную сумму
        c = C(i + 1,:);
        % кодириуем
        cod = Cod(i + 1,:);
        % модулируем
        r = R(i + 1,:);
        while q ~= 0
            % добавляем шум
            noise = sigma(n)*randn(1, length(r));
            r_noise = r + noise;
            % демодуляция
            c_demod = r_noise < 0;
            % добавление нулей
            c_demod = [[0 0], c_demod];
            % декодирование
            decod = c_demod;
            %пересчитываем контрольную сумму
            e = zeros(1, 5);
            for j = 0:length(decod) - 1
                if decod(j + 1) == 1
                    z = (num2str(str2double(dec2bin(length(decod) - j))) - '0');
                    z = [zeros(1, 5 - length(z)), z];
                    e = e + z;
                end
            end
            contr_sum = mod(e, 2);
      
            pos = bi2de(contr_sum(end:-1:1));
            if pos == 0
                new_decod = zeros(1, length(decod) - 5);
                w = 1;
                h = 1;
                for j = 0: length(decod) - 1
                    if j + 1 ~= check_bits(1, w)
                        new_decod(1, h) = c_demod(1, j + 1);
                        h = h + 1;
                    else
                        w = w + 1;
                    end 
                end
                % убираем добавленные нули
                new_decod = new_decod(1, 3:length(new_decod));
                [p, s] = deconv(new_decod, G);
                s = mod(s, 2);
                q = sum(s);
               
                per = per + 1;
                b = sum(c ~= new_decod);
                Pb(n) = Pb(n)+ b/length(new_decod);
                Pdec(n) = Pdec(n) + ((q == 0)&(b > 0));
                q = 0;
            else
            % исправляем ошибку 
                pos = length(cod) - pos + 1;
                c_demod(1, pos) = mod(decod(1, pos) + 1, 2);
                % убираем проверочные биты
                new_decod = zeros(1, length(decod) - 5);
                w = 1;
                h = 1;
                for j = 0: length(decod) - 1
                    if j + 1 ~= check_bits(1, w)
                        new_decod(1, h) = c_demod(1, j + 1);
                        h = h + 1;
                    else
                        w = w + 1;
                    end 
                end
                % убираем добавленные нули
                new_decod = new_decod(1, 3:length(new_decod));
                % считаем синдром
                [p, s] = deconv(new_decod, G);
                s = mod(s, 2);
                q = sum(s);
                per = per + 1;
                b = sum(c ~= new_decod);
                Pb(n) = Pb(n)+ b/length(new_decod);
                Pdec(n) = Pdec(n) + ((q == 0)&(b > 0));
            end
        end
        [SNRdB(n), i, per]
        n_exp = n_exp + per;
        corr_mes = corr_mes + 1;
    end
    Pb(n) = Pb(n)/n_exp;
    Pdec(n) = Pdec(n)/n_exp;
    Ped_th(n) = (2^k - 1) * Pb_theor(n)^d_min;
    Ped_pr(n) = (2^k - 1) * Pb(n)^d_min;
    T(n) = k * num_mes/(length_cod*n_exp);
end

figure(1);
semilogy(SNRdB, Pdec, 'r', SNRdB, Ped_th, 'g', SNRdB, Ped_pr, 'bo', 'LineWidth', 2)
% semilogy(SNRdB, Ped_pr, 'r', SNRdB, Ped_th, 'g', SNRdB, Ped_pr, 'bo', 'LineWidth', 2)
grid on
legend('Вероятность ошибки декодирования(при жестком декодировании)', 'Вероятность ошибки декодирования(при мягком декодировании)', 'Теоретическая вероятность ошибки декодирования');
title('Вероятности ошибки декодирования от SNRdB');
xlabel('SNRdB');
ylabel('Вероятность ошибки декодирования');

figure(2);
plot(SNRdB, T, 'r', 'LineWidth', 2)
grid on
title('Пропускная способность от SNRdB');
xlabel('SNRdB');
ylabel('Пропускная способность, T');

figure(3);
plot(SNRdB, Pb_theor, 'r', SNRdB, Pb, 'g', 'LineWidth', 2)
grid on
legend('Вероятность ошибки на бит(теория)', 'Вероятность ошибки на бит(практика)');
title('Вероятности ошибки на бит от SNRdB');
xlabel('SNRdB');
ylabel('Вероятность ошибки на бит');