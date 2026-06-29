%==========================================================================
% Simülasyon (kodlu, rate-1/2 BCC+Viterbi) ile literatürdeki kodsuz AWGN
% teorik BER eğrilerinin karşılaştırması.
%
% Eksen: alt-taşıyıcı başına Es/N0 (dB). Bizim simülasyonun "paket örnek
% SNR'ı" Es/N0'a çevriliyor:  Es/N0 = packetSNR * (Nfft/numUsed) / Ppacket
%==========================================================================
clear; clc; close all;
Nfft = 64; numUsed = 52;
Ppkt = 0.99;                                  % ~bütün-paket ort. gücü (DATA=1)
offset_dB = 10*log10(Nfft/numUsed) - 10*log10(Ppkt);   % ~ +1.0 dB
fprintf('packetSNR -> Es/N0 ofseti = %.2f dB\n', offset_dB);

snr = 0:2:24;                                 % simülasyonun paket örnek SNR'ı
% --- Simüle CODED BER (15 paket x 60 OFDM sembolü üzerinden ortalanmış) ---
berQPSK = [0.47455 0.32489 0.13377 0.0077708 0 0 0 0 0 0 0 0 0];
ber16   = [0.50124 0.49824 0.47756 0.41201 0.20969 0.041108 0.0018422 0.000336 0 0 0 0 0];
ber64   = [0.49935 0.49893 0.49361 0.48594 0.46148 0.39375 0.26721 0.096734 0.012787 0.0012354 2.3164e-5 0 0];
EsN0sim = snr + offset_dB;

% --- Teorik kodsuz AWGN BER (Gray-kodlu) vs Es/N0 ---
Q  = @(x) 0.5*erfc(x/sqrt(2));
es = 0:0.2:32; esl = 10.^(es/10);
T_QPSK = Q(sqrt(esl));                         % QPSK
T_16   = 0.75 * Q(sqrt(0.2*esl));              % 16-QAM
T_64   = (7/12) * Q(sqrt(esl/21));             % 64-QAM

figure('Position',[100 100 760 560]); hold on; grid on;
c1=[0 0.45 0.74]; c2=[0.85 0.33 0.1]; c3=[0.47 0.67 0.19];
semilogy(es,T_QPSK,'--','Color',c1,'LineWidth',1.6);
semilogy(es,T_16,  '--','Color',c2,'LineWidth',1.6);
semilogy(es,T_64,  '--','Color',c3,'LineWidth',1.6);
psim(EsN0sim,berQPSK,c1); psim(EsN0sim,ber16,c2); psim(EsN0sim,ber64,c3);
set(gca,'YScale','log'); ylim([1e-5 1]); xlim([0 30]);
xlabel('E_s/N_0  (alt-taşıyıcı başına, dB)'); ylabel('BER');
title('Simülasyon (kodlu, rate 1/2)  vs  Teori (kodsuz AWGN)');
legend({'QPSK teori (kodsuz)','16-QAM teori (kodsuz)','64-QAM teori (kodsuz)', ...
        'QPSK sim (kodlu)','16-QAM sim (kodlu)','64-QAM sim (kodlu)'}, ...
        'Location','southwest','FontSize',9);
print(gcf,'compare_BER.png','-dpng','-r120');
fprintf('compare_BER.png kaydedildi.\n');

function psim(x,b,c)
    m = b>0;
    semilogy(x(m),b(m),'o-','Color',c,'LineWidth',1.6,'MarkerFaceColor',c,'MarkerSize',6);
end
