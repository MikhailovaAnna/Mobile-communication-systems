clc
clear all
close all
nfig = 1;

% �������� ��������� �� 32 �����
k = 32;
Message = zeros(1, k);

% ��������� ����������� �����
G = [1 0 0 1 1 1 1 0 1 0 1 1 0 0 1 0 1];
size_g = length(G);
length_C = k + size_g - 1;
CRC = [Message, zeros(1, size_g - 1)];
[b, r] = deconv(CRC, G);
r_t = mod(r, 2);
CRC = CRC + r_t; % ��������� � ����������� ������

tblen1 = 4;
trellis1 = poly2trellis(tblen1, [13 17]); % �������� ������� ��� �������� 1/2
Cod1 = convenc(CRC, trellis1); % ������������
Model1 = Cod1*(-2) + 1; % ��������������

num_mes = 200000; % ����� ���������, ������� �� ����� ��������
SNRdB = -12:2:0;
SNR = 10.^(SNRdB/10);
sigma = sqrt((1 ./ (2 * SNR)));
T1 = zeros(1, length(SNRdB)); % ���������� ����������� HARQ1
Pdec1 = zeros(1, length(SNRdB)); % ����������� ������ ������������� HARQ1
for n = 1:length(SNRdB)
    n_exp = 0; % ����� �������
    corr_mes = 0; % ����� ��������� ���������� ���������
    bufer = zeros(1, 96); 
    while corr_mes < num_mes
        rnd_x = randn();
        rnd_y = randn();
        H = (rnd_x*rnd_x + rnd_y*rnd_y)^0.5; % ���������
        noise = sigma(n)*randn(1, length(Model1));
        r = Model1*H + noise; % ����� ������
        
        % �������������
        decod = vitdec(r, trellis1, tblen1, 'trunc', 'unquant'); % ������ �������������
  
        [p, s] = deconv(decod, G);
        s = mod(s, 2);
        q = sum(s);
        b = sum(CRC ~= decod);
        if q ~= 0 % ������� �� ����� ����, ������� ������������ ������ �� ������
            bufer = bufer + r; % ������������ ��������� � ������
            new_decod = vitdec(bufer, trellis1, tblen1, 'trunc', 'unquant'); % ����������
            [p, s] = deconv(new_decod, G);
            s = mod(s, 2);
            q = sum(s);
            b = sum(CRC ~= new_decod);
            if q == 0 % ���� ������� ����� ���� - ��������� ������� ��������������
                corr_mes = corr_mes + 1;
                bufer = zeros(1, 96); % ������� �����
            end
        else
            corr_mes = corr_mes + 1;
        end
        Pdec1(n) = Pdec1(n) + ((q == 0)&(b > 0)); % ���� ������� ����� ����, � ����� ����������� �������������� ��� ������ ����, �� ���������� ������ �������������
        n_exp = n_exp + 1;
    end
    T1(n) = k * num_mes/(length(Cod1)*n_exp); % ������� ���������� �����������
    Pdec1(n) = Pdec1(n)/n_exp;
end

tblen2 = [5 4];
trellis2_7 = poly2trellis(tblen2, [23 35 0 32 07 16 14; 0 5 13 07 13 05 10]); % ������� ��� ���� �� ��������� 2/7
trellis2_5 = poly2trellis(tblen2, [23 35 0 32 07; 0 5 13 07 13]); % ������� ��� ���� �� ��������� 2/5
trellis2_3 = poly2trellis(tblen2, [23 35 0; 0 5 13]); % ������� ��� ���� �� ��������� 2/3
Cod2_7 = convenc(CRC, trellis2_7); % ������� �������� ��������� � ����� ������ ���������

% ����� �� ���� �����1, �����2, �����3 ����� ��� ����������������
% ���������, ��� ����� ������������� � �������� �������� ������
P1 = zeros(1, 3*(length(Cod2_7)/7)); % ������� ������ �����, ������� ����� �������� �������� � ������������ ��������� 2/3
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

P2 = zeros(1, 2*(length(Cod2_7)/7)); % ���������� ����� ������, ��� �������� �� ��������� 2/5
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

P3 = zeros(1, 2*(length(Cod2_7)/7)); % ����� 3, ���������� ���� ��� �������� �� ��������� 2/7
ModelP1 = P1*(-2) + 1;
ModelP2 = P2*(-2) + 1;
ModelP3 = P3*(-2) + 1;

T2 = zeros(1, length(SNRdB)); % ���������� ����������� HARQ2
Pdec2 = zeros(1, length(SNRdB)); % ����������� ������ ������������� HARQ2
av_len = 0; % ������� ����� ��������������� ���������, ������� �������������� � ��������� ������ 0
for n = 1:length(SNRdB)
    n_exp = 0;
    corr_mes = 0;
    len = 0; % ����� ��������������� ���������, ������� ��������������
    while corr_mes < num_mes
        q = -1;
        rnd_x = randn();
        rnd_y = randn();
        H = (rnd_x*rnd_x + rnd_y*rnd_y)^0.5; % ���������
        noise2_3 = sigma(n)*randn(1, length(ModelP1));
        r2_3 = ModelP1*H + noise2_3;
        % �������������
        decod2_3 = vitdec(r2_3, trellis2_3, 5, 'trunc', 'unquant');
        [p, s3] = deconv(decod2_3, G);
        s3 = mod(s3, 2);
        q3 = sum(s3);
        b = sum(CRC ~= decod2_3);
        if q3 ~= 0 % ���� �� ������ ������� ���������� �� ��������� 2/3, �� �������� �� ��������� 2/5
            noiseP2 = sigma(n)*randn(1, length(ModelP2));
            rP2 = ModelP2*H + noiseP2;
            r2_5 = zeros(1, length(r2_3)+length(rP2)); % ��������� �����, ������� �������� � ���� ���������� �����, ������� �� ��������� �� ��������� 2/3, � �����, ������� �����
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
            decod2_5 = vitdec(r2_5, trellis2_5, 5, 'trunc', 'unquant'); % �������� ������������
            [p, s5] = deconv(decod2_5, G);
            s5 = mod(s5, 2);
            q5 = sum(s5);
            b = sum(CRC ~= decod2_5);
            if q5 ~= 0 % ���� �� ���������� �� ��������� 2/5, �������� �������� �� ��������� 2/7
                noiseP3 = sigma(n)*randn(1, length(ModelP3));
                rP3 = ModelP3*H + noiseP3;
                r2_7 = zeros(1, length(Cod2_7)); % �������� � ���� ����� ��� �������� 2/5 � ����� �������� ������
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
                decod2_7 = vitdec(r2_7, trellis2_7, 5, 'trunc', 'unquant'); % �������� ������������
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
        Pdec2(n) = Pdec2(n) + ((q == 0)&(b > 0)); % ������� ����������� ������ �������������
        n_exp = n_exp + 1;
        av_len = av_len + len; 
    end
    av_len = av_len / num_mes;
    T2(n) = k * num_mes/(av_len*n_exp); % ������� ���������� �����������
    Pdec2(n) = Pdec2(n)/n_exp;
end

figure(1);
plot(SNRdB, T1, 'r', SNRdB, T2, 'b', 'LineWidth', 2)
grid on
legend('���������� ����������� T_1', '���������� ����������� T_2');
title('���������� ����������� �� SNRdB');
xlabel('SNRdB');
ylabel('���������� �����������, T');

figure(2);
semilogy(SNRdB, Pdec1, 'r', SNRdB, Pdec2, 'b', 'LineWidth', 2)
grid on
legend('����������� ������ ������������� Type1', '����������� ������ ������������� Type2');
title('����������� ������ ������������� �� SNRdB');
xlabel('SNRdB');
ylabel('����������� ������ �������������');
