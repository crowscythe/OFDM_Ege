%==========================================================================
% AI IS USED TO MAKE COMMENTS ETC.
% OFDM transceiver using only base MATLAB functions (no toolboxes)
%
% TX chain:
%   random bits -> scrambler -> FEC (Conv Encoder) -> block interleaver 
%   -> QAM mapper -> pilot insertion -> IFFT -> cyclic prefix -> AWGN
%
% RX chain:
%   remove CP -> FFT -> subcarrier extraction -> channel estimate
%   -> equalization -> QAM demapper -> deinterleaver -> FEC (Viterbi) 
%   -> descrambler
%
% Supported modulations: 'QPSK', '16QAM', '64QAM'
% FEC: Rate 1/2 Convolutional Code (K=3), Polynomials [7, 5] in octal.
% Noise: AWGN defined by Eb/N0 sweep.
%==========================================================================
clear; clc; close all;
rng(1);  % repeatable data and noise

%% ========================= User parameters ==============================
modType = '64QAM';              % 'QPSK', '16QAM', or '64QAM'
Nfft = 64;                      % FFT/IFFT size (must be even)
cpLen = 16;                     % cyclic prefix length

% INCREASED symbol count to resolve lower BER (10^-5 to 10^-6 range)
numOFDMSymbols = 5000;          

% INCREASED resolution of the sweep
EbNoDbVec = 0:1:22;             

channelType = 'AWGN';           % currently implemented: 'AWGN'

% Carrier allocation (Zero-centered)
usedCarriers  = [-26:-1 1:26]; 
pilotCarriers = [-21 -7 7 21]; 

%% =======================================================================
modType = upper(char(modType));

% Input validation
if Nfft < 2 || mod(Nfft,2) ~= 0, error('Nfft must be a positive even integer.'); end
if cpLen < 0 || cpLen >= Nfft, error('cpLen must be 0 <= cpLen < Nfft.'); end
usedCarriers = sort(usedCarriers(:).');
pilotCarriers = sort(pilotCarriers(:).');
validCarriers = -Nfft/2:(Nfft/2-1); validCarriers(validCarriers == 0) = [];
dataCarriers = usedCarriers(~ismember(usedCarriers, pilotCarriers));

%% Modulation setup
switch modType
    case 'QPSK'
        bitsPerSym = 2;
        normFactor = sqrt(2);
    case '16QAM'
        bitsPerSym = 4;
        normFactor = sqrt(10);
    case '64QAM'
        bitsPerSym = 6;
        normFactor = sqrt(42);
    otherwise
        error('Unsupported modulation.');
end

%% Derived parameters & FEC Dimensions
fftCenter = Nfft/2 + 1;                 
dataBins = dataCarriers + fftCenter;    
pilotBins = pilotCarriers + fftCenter;  
numDataSubc = length(dataCarriers);
numPilotSubc = length(pilotCarriers);

bitsPerOFDMSymbol = numDataSubc * bitsPerSym;
numDataSymbols = numDataSubc * numOFDMSymbols;
numTxCodedBits = bitsPerOFDMSymbol * numOFDMSymbols;
ofdmSymbolLen = Nfft + cpLen;

% Interleaver setup
candidateRows = [16 12 8 6 4 3 2];
interRows = 1;
for k = 1:length(candidateRows)
    if mod(bitsPerOFDMSymbol, candidateRows(k)) == 0
        interRows = candidateRows(k); break;
    end
end
interCols = bitsPerOFDMSymbol / interRows;

% FEC Setup: Rate 1/2
numUncodedBits = (numTxCodedBits / 2) - 2; 

fprintf('Modulation                : %s\n', modType);
fprintf('Nfft / CP length          : %d / %d\n', Nfft, cpLen);
fprintf('Used / Pilot / Data bins  : %d / %d / %d\n', length(usedCarriers), numPilotSubc, numDataSubc);
fprintf('Payload coded bits/frame  : %d\n', numTxCodedBits);
fprintf('Payload data bits/frame   : %d (after FEC)\n\n', numUncodedBits);

%% Generate transmit bits
txUncodedBits = double(rand(numUncodedBits,1) > 0.5);

%% Scrambler
pnSeq = zeros(numUncodedBits,1);
lfsr = [1 0 1 1 0 1 1]; 
for n = 1:numUncodedBits
    pnSeq(n) = lfsr(end);
    feedbackBit = mod(lfsr(7) + lfsr(4), 2);
    lfsr(2:end) = lfsr(1:end-1);
    lfsr(1) = feedbackBit;
end
scrambledUncodedBits = mod(txUncodedBits + pnSeq, 2);

%% FEC: Rate-1/2 Convolutional Encoder
% Generator polynomials: [7, 5] in octal -> 111 and 101 binary
paddedBits = [scrambledUncodedBits; 0; 0]; % Append 2 flush bits
txCodedBits = zeros(numTxCodedBits, 1);
reg = [0 0];
for k = 1:length(paddedBits)
    b = paddedBits(k);
    out1 = mod(b + reg(1) + reg(2), 2);
    out2 = mod(b + reg(2), 2);
    txCodedBits(2*k - 1) = out1;
    txCodedBits(2*k)     = out2;
    reg = [b reg(1)];
end

%% Block interleaver
scrambledBitBlocks = reshape(txCodedBits, bitsPerOFDMSymbol, numOFDMSymbols);
interleavedBitBlocks = zeros(size(scrambledBitBlocks));
for symIdx = 1:numOFDMSymbols
    temp = reshape(scrambledBitBlocks(:,symIdx), interRows, interCols);
    interleavedBitBlocks(:,symIdx) = reshape(temp.', [], 1);
end
interleavedBits = interleavedBitBlocks(:);

%% QAM mapping
txBitGroups = reshape(interleavedBits, bitsPerSym, []).';
switch modType
    case 'QPSK'
        I = 1 - 2*txBitGroups(:,1); Q = 1 - 2*txBitGroups(:,2);
        txDataSymbols = (I + 1j*Q) / normFactor;
    case '16QAM'
        bI1 = txBitGroups(:,1); bI2 = txBitGroups(:,2);
        bQ1 = txBitGroups(:,3); bQ2 = txBitGroups(:,4);
        I = -(1 - 2*bI1) .* (3 - 2*bI2); Q = -(1 - 2*bQ1) .* (3 - 2*bQ2);
        txDataSymbols = (I + 1j*Q) / normFactor;
    case '64QAM'
        gray8_binaryOrder = [-7 -5 -1 -3 7 5 1 3];
        idxI = txBitGroups(:,1)*4 + txBitGroups(:,2)*2 + txBitGroups(:,3) + 1;
        idxQ = txBitGroups(:,4)*4 + txBitGroups(:,5)*2 + txBitGroups(:,6) + 1;
        I = gray8_binaryOrder(idxI); Q = gray8_binaryOrder(idxQ);
        txDataSymbols = (I + 1j*Q) / normFactor;
end

%% Pilot insertion, IFFT, and CP
pilotPattern = ones(numPilotSubc,1);
pilotPattern(2:2:end) = -1;      
pilotMatrix = repmat(pilotPattern, 1, numOFDMSymbols);
txDataMatrix = reshape(txDataSymbols, numDataSubc, numOFDMSymbols);

txGrid = zeros(Nfft, numOFDMSymbols);
txGrid(dataBins, :) = txDataMatrix;
txGrid(pilotBins, :) = pilotMatrix;

txTime = ifft(ifftshift(txGrid, 1), Nfft, 1);
if cpLen > 0
    txWithCP = [txTime(end-cpLen+1:end,:); txTime];
else
    txWithCP = txTime;
end
txSerial = txWithCP(:);
txSignalPower = mean(abs(txSerial).^2);

%% BER simulation (Eb/N0 Sweep)
ber = zeros(size(EbNoDbVec));
totalSamples = numOFDMSymbols * ofdmSymbolLen;

for snrIdx = 0:length(EbNoDbVec)
    if snrIdx == 0
        currentEbNoDb = Inf; % No-noise sanity check first
    else
        currentEbNoDb = EbNoDbVec(snrIdx);
    end
    
    %% AWGN Channel Setup using Eb/N0
    if isinf(currentEbNoDb)
        noise = zeros(size(txSerial));
    else
        % SNR conversion: accounts for CP, empty bins, and coding rate
        snrLin = (10^(currentEbNoDb/10)) * (numUncodedBits / totalSamples);
        noiseVar = txSignalPower / snrLin;
        noise = sqrt(noiseVar/2) * (randn(size(txSerial)) + 1j*randn(size(txSerial)));
    end
    rxSerial = txSerial + noise;

    %% Receiver: remove CP, FFT
    rxWithCP = reshape(rxSerial, ofdmSymbolLen, numOFDMSymbols);
    rxNoCP = rxWithCP(cpLen+1:end, :);
    rxGrid = fftshift(fft(rxNoCP, Nfft, 1), 1);

    %% Channel estimation
    rxPilots = rxGrid(pilotBins, :);
    rxData = rxGrid(dataBins, :);
    HhatPilots = rxPilots ./ pilotMatrix;
    Hhat = mean(HhatPilots, 1);                          
    rxDataEq = rxData ./ repmat(Hhat, numDataSubc, 1);  

    %% QAM demapper
    rxDataSymbols = rxDataEq(:);
    switch modType
        case 'QPSK'
            rxBitGroups = zeros(numDataSymbols, 2);
            rxBitGroups(:,1) = real(rxDataSymbols) < 0;
            rxBitGroups(:,2) = imag(rxDataSymbols) < 0;
        case '16QAM'
            rxScaled = rxDataSymbols * normFactor; levels = [-3 -1 1 3];
            bitsAxis = [0 0; 0 1; 1 1; 1 0];
            [~, idxIhat] = min(abs(repmat(real(rxScaled), 1, 4) - repmat(levels, numDataSymbols, 1)), [], 2);
            [~, idxQhat] = min(abs(repmat(imag(rxScaled), 1, 4) - repmat(levels, numDataSymbols, 1)), [], 2);
            rxBitGroups = [bitsAxis(idxIhat,:), bitsAxis(idxQhat,:)];
        case '64QAM'
            rxScaled = rxDataSymbols * normFactor; levels = [-7 -5 -3 -1 1 3 5 7];
            bitsAxis = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 1 1; 1 0 1; 1 0 0];
            [~, idxIhat] = min(abs(repmat(real(rxScaled), 1, 8) - repmat(levels, numDataSymbols, 1)), [], 2);
            [~, idxQhat] = min(abs(repmat(imag(rxScaled), 1, 8) - repmat(levels, numDataSymbols, 1)), [], 2);
            rxBitGroups = [bitsAxis(idxIhat,:), bitsAxis(idxQhat,:)];
    end
    rxInterleavedBits = reshape(rxBitGroups.', [], 1);

    %% Deinterleaver
    rxInterleavedBitBlocks = reshape(rxInterleavedBits, bitsPerOFDMSymbol, numOFDMSymbols);
    rxDeinterleavedBitBlocks = zeros(size(rxInterleavedBitBlocks));
    for symIdx = 1:numOFDMSymbols
        temp = reshape(rxInterleavedBitBlocks(:,symIdx), interCols, interRows).';
        rxDeinterleavedBitBlocks(:,symIdx) = temp(:);
    end
    rxCodedBits = rxDeinterleavedBitBlocks(:);

    %% FEC
    N_trellis = numTxCodedBits / 2;
    nextState = [1 2; 3 4; 1 2; 3 4];  % State transitions
    expectedOut = zeros(4, 2, 2);      % Expected output bits
    expectedOut(1,1,:) = [0 0]; expectedOut(1,2,:) = [1 1];
    expectedOut(2,1,:) = [1 0]; expectedOut(2,2,:) = [0 1];
    expectedOut(3,1,:) = [1 1]; expectedOut(3,2,:) = [0 0];
    expectedOut(4,1,:) = [0 1]; expectedOut(4,2,:) = [1 0];

    PM = [0; Inf; Inf; Inf];
    survivor = zeros(4, N_trellis);
    inputHistory = zeros(4, N_trellis); 

    for k = 1:N_trellis
        r1 = rxCodedBits(2*k - 1);
        r2 = rxCodedBits(2*k);
        newPM = Inf(4,1);
        for state = 1:4
            if isinf(PM(state)), continue; end
            for inBit = 0:1
                nState = nextState(state, inBit+1);
                e1 = expectedOut(state, inBit+1, 1);
                e2 = expectedOut(state, inBit+1, 2);
                bm = (r1 ~= e1) + (r2 ~= e2); % Hamming branch metric
                if PM(state) + bm < newPM(nState)
                    newPM(nState) = PM(state) + bm;
                    survivor(nState, k) = state;
                    inputHistory(nState, k) = inBit;
                end
            end
        end
        PM = newPM;
    end

    % Traceback
    rxScrambledUncodedBits = zeros(numUncodedBits, 1);
    currState = 1; % Force finish at state 00 due to flush zeros
    for k = N_trellis:-1:1
        prevState = survivor(currState, k);
        if k <= numUncodedBits
            rxScrambledUncodedBits(k) = inputHistory(currState, k);
        end
        currState = prevState;
    end

    %% Descrambler
    rxUncodedBits = mod(rxScrambledUncodedBits + pnSeq, 2);

    %% BER Calculation
    bitErrors = sum(rxUncodedBits ~= txUncodedBits);
    if snrIdx == 0
        if bitErrors ~= 0
            error('No-noise sanity check FAILED: %d bit errors.', bitErrors);
        else
            fprintf('No-noise sanity check PASSED: 0 bit errors.\n\n');
        end
    else
        ber(snrIdx) = bitErrors / numUncodedBits;
        fprintf('Eb/N0 = %5.1f dB, BER = %12.6e, bit errors = %d / %d\n', ...
                currentEbNoDb, ber(snrIdx), bitErrors, numUncodedBits);
    end
end

%% Plot BER vs Eb/N0
berForPlot = ber;
berForPlot(berForPlot == 0) = 0.5 / numUncodedBits; % Protect log scale
figure;
semilogy(EbNoDbVec, berForPlot, 'o-','LineWidth',1.5,'MarkerSize',6);
grid on;
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title(sprintf('OFDM, BER vs Eb/N0 (%s)', modType));
