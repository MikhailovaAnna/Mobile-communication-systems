clc
clear all
close all
nfig = 1;

k = 4;
% все слова
M = zeros(2^k, k);
for i = 0 : 2^k - 1
    str = dec2bin(i,k);
    for j = 1: k
        M(i + 1, j) = str(j) - 48;
    end
end

% добавление контрольной суммы
G = [1 0 0 1 1 1 1 0 1 0 1 1 0 0 1 0 1]; % порождающий многочлен 16-й степени, младщие разряды справа
size_g = length(G);
t = zeros(1, size_g - 1);
C = zeros(2^k, k + size_g - 1);
for i = 0: 2^k - 1
    C(i + 1,:) = [M(i + 1,:), t]; % умножаем сообщение на х в степени порождающего многочлена
    [q, r] = deconv(C(i + 1,:), G); % находим частное и остаток от деления
    r_t = mod(r, 2);
    C(i + 1,:) = C(i + 1,:) + r_t;
end

% создание книги весов
weight_word = sum(C(2:end,:), 2);
d_min = min(weight_word); % минимальное расстояние кода
weight_book = zeros(1, 2^k);
for i = 1:2^k
    weight_book(i) = sum(weight_word == (i));
end

% создание книги модуляций
R = zeros(2^k, k + size_g - 1);
R = C.*(-2) + 1;

exp = 200000;
SNRdB = -20:5:10;
Peb = zeros(1, length(SNRdB)); % практическая вероятность ошибки на бит
Pdec = zeros(1, length(SNRdB)); % практическая ошибка декодирования
Ped = zeros(1, length(SNRdB)); % более точная верхняя граница декодирования
Ped_t = zeros(1, length(SNRdB)); % точная верхняя граница декодирования 
Peb_theor = erfc(sqrt(2 * 10 .^ (SNRdB / 10)) ./ sqrt(2)) * 0.5; % теоретическая вероятность ошибки на бит
SNR = 10.^(SNRdB/10);
for n = 1:length(SNRdB)
    n_exp = 0;
    A_i = zeros(exp, k + size_g - 1);
    for n_exp = 0:exp
        i = randi([0, 2^k - 1]);
        m = M(i + 1,:);
        c = C(i + 1,:);
        r = R(i + 1,:);
        A_i(n_exp + 1, :) = c;
        % начинаем шумить
        sigma = sqrt((1 / (2 * SNR(n))));
        noise = sigma*randn(1, k + size_g - 1);
        r_noise = r + noise;
        % демоделирование
        c_demod = r_noise < 0;
        % вычисление синдрома и вероятностей ошибки
        [p, s] = deconv(c_demod, G);
        s = mod(s, 2);
        q = sum(s);
        b = sum(c ~= c_demod); % mod(c + c_demod, 2);
        Peb(n) = Peb(n)+ b/length(c_demod);
        Pdec(n) = Pdec(n) + ((q == 0)&(b > 0));
    end
    Peb(n) = Peb(n)/exp;
    Pdec(n) = Pdec(n)/exp;
    Ped(n) =(2^length(m) - 1) * Peb_theor(n)^d_min;
    for i = d_min:length(weight_book)
        Ped_t(n) = Ped_t(n) + weight_book(i)*(Peb_theor(n)^i)*(1 - Peb_theor(n)^(length(weight_book) - i));
    end
    Ped_t(n) = Ped_t(n)/exp;
end

figure(1);
semilogy(SNRdB, Peb, 'ro', SNRdB, Peb_theor, 'b--', 'LineWidth', 2)
grid on
legend('Вероятность ошибки на бит(практика)', 'Вероятность ошибки на бит(теория)');
title('Вероятность ошибки на бит от SNRdB');
xlabel('SNRdB');
ylabel('Pb');

figure(2);
plot(SNRdB, Pdec, 'ro', SNRdB, Ped, 'b', SNRdB, Ped_t, 'g', 'LineWidth', 2)
grid on
legend('Практическая вероятность ошибки декодирования', 'Точная верхняя граница вероятности ошибки декодирования',  'Точная вероятность ошибки декодирования');
title('Вероятность ошибки декодирования от SNRdB');
xlabel('SNRdB');
ylabel('Pdec');