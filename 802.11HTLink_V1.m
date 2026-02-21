clear; clc;
%ege alper uzun 23050211026
%heavily referenced from https://www.mathworks.com/help/wlan/ug/802-11n-packet-error-rate-simulation-for-2x2-tgn-channel.html
%ht config

cfgHT = wlanHTConfig;

%disp(cfgHT) %uncomment to see ht configs
%sim parameters, mcs stands for modulation scheme

snr = 10:1:22; 
%i removed bpsk and qpsk from this simulation since they wouldnt
%get any meaningful errors in these SNR ranges with our current simulation

mcsRange = [3,4,5,6,7];%16QAM 1/2CR 3/4 CR, 64QAM 2/3 3/4 5/6 CR respectively

%i removed MCS 0 1 2 they wouldn't generate any bit errors(close to 0)
%higher mcs could be used with VHT but since openofdm supports up to MCS 7
maxNumBits = 1e6; 
maxNumPackets = 1e3;
bitErrorRate = zeros(numel(snr), numel(mcsRange));

%snr steps could be increased, max num bits and maxnumpackets can be
%reduced for better performance.

% process
for m = 1:numel(mcsRange)
    
    cfgHT.MCS = mcsRange(m);
    ofdmInfo = wlanHTOFDMInfo('HT-Data', cfgHT);
    ind = wlanFieldIndices(cfgHT);
    
    fprintf('SNR(dB)\t\tBER\n');
    fprintf('-----------------------------\n');
    
    for i = 1:numel(snr)
        
        packetSNR = snr(i) - 10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
        
        totalBitErrors = 0;
        totalBitsTransmitted = 0;
        numPackets = 0;
        
        while totalBitsTransmitted < maxNumBits && numPackets < maxNumPackets
            
            % waveform gen
            txPSDU = randi([0 1], cfgHT.PSDULength*8, 1);
            tx = wlanWaveformGenerator(txPSDU, cfgHT);
            % awgn
            rx = awgn(tx, packetSNR, 'measured');
            % siso, not 2x2 for now
            chanEst = ones(ofdmInfo.NumTones, 1, 1);      
            htdata = rx(ind.HTData(1):ind.HTData(2), :);
            % noise variance
            nVarHT = 10^(-snr(i)/10);
            rxPSDU = wlanHTDataRecover(htdata, chanEst, nVarHT, cfgHT);
            % bit error
            numErrors = biterr(txPSDU, rxPSDU);
            totalBitErrors = totalBitErrors + numErrors;
            totalBitsTransmitted = totalBitsTransmitted + numel(txPSDU);
            numPackets = numPackets + 1;
        end
        
        % BER calc and print
        bitErrorRate(i, m) = totalBitErrors / totalBitsTransmitted;
        
        fprintf('%2d\t\t%e\n', snr(i), bitErrorRate(i,m));
    end
end
% display
disp('BER Table');
T = array2table([snr.', bitErrorRate]);
varNames = {'SNR_dB'};
for k = 1:numel(mcsRange)
    varNames{end+1} = sprintf('BER_MCS%d', mcsRange(k));
end
T.Properties.VariableNames = varNames;
disp(T);

% plot
figure;
semilogy(snr, bitErrorRate, '-o', 'LineWidth', 1.5);
grid on;
axis([min(snr) max(snr) 1e-6 1]);
xlabel('SNR (dB)');
ylabel('BER');
title('802.11n BER Simulation (AWGN)');

% legend
legendEntries = cell(1, numel(mcsRange));
for k = 1:numel(mcsRange)
    legendEntries{k} = sprintf('MCS %d', mcsRange(k));
end
legend(legendEntries, 'Location', 'southwest');
