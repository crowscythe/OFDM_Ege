```matlab
%==========================================================================
% Simplified OFDM transceiver using only base MATLAB functions (no toolboxes)
%
%==========================================================================

clear; clc; close all;
rng(1);  % repeatable data and noise

%% ========================= User parameters ==============================
modType = '16QAM';              % 'QPSK', '16QAM', or '64QAM'
Nfft = 64;                      % FFT/IFFT size (default 802.11-style: 64)
cpLen = 16;                     % cyclic prefix length
numOFDMSymbols = 500;           % HT-DATA OFDM symbols
snrDbVec = 0:2:24;              % SNR sweep in dB (packet sample SNR)
channelType = 'AWGN';           % currently implemented: 'AWGN'

% Zero-centered subcarrier numbering:
% valid bins for even Nfft are $-Nfft/2 ... -1, 1 ... Nfft/2-1$
% Default below is a common 64-point OFDM layout with 52 used carriers.
usedCarriers  = [-26:-1 1:26];  % all occupied bins except DC
pilotCarriers = [-21 -7 7 21];  % fixed pilot locations, subset of usedCarriers

% Simplified HT packet field configuration
numHTSTFRepeats = 10;           % repeated short symbols
numHTLTFSymbols = 1;            % simplified single-stream style
numHTSIGSymbols = 2;            % HT-SIG as OFDM symbol(s), default 2
%% =======================================================================

modType = upper(char(modType));
channelType = upper(char(channelType));

%% Basic parameter checks
if ~isvector(snrDbVec)
    error('snrDbVec must be a vector.');
end

if Nfft < 2 || Nfft ~= floor(Nfft) || mod(Nfft,2) ~= 0
    error('Nfft must be a positive even integer.');
end

% Simplified HT-STF generation below uses every 4th carrier so that the
% short training waveform repeats with period $Nfft/4$.
if mod(Nfft,4) ~= 0
    error('Simplified HT-STF generation requires Nfft divisible by 4.');
end

if cpLen ~= floor(cpLen) || cpLen < 0 || cpLen >= Nfft
    error('cpLen must be an integer with 0 <= cpLen < Nfft.');
end

if numOFDMSymbols < 1 || numOFDMSymbols ~= floor(numOFDMSymbols)
    error('numOFDMSymbols must be a positive integer.');
end

if numHTSTFRepeats < 1 || numHTSTFRepeats ~= floor(numHTSTFRepeats)
    error('numHTSTFRepeats must be a positive integer.');
end

if numHTLTFSymbols < 1 || numHTLTFSymbols ~= floor(numHTLTFSymbols)
    error('numHTLTFSymbols must be a positive integer.');
end

if numHTSIGSymbols < 1 || numHTSIGSymbols ~= floor(numHTSIGSymbols)
    error('numHTSIGSymbols must be a positive integer.');
end

usedCarriers = usedCarriers(:).';
pilotCarriers = pilotCarriers(:).';

if any(usedCarriers ~= round(usedCarriers))
    error('usedCarriers must contain integer bin indices.');
end

if any(pilotCarriers ~= round(pilotCarriers))
    error('pilotCarriers must contain integer bin indices.');
end

if length(unique(usedCarriers)) ~= length(usedCarriers)
    error('usedCarriers contains duplicates.');
end

if length(unique(pilotCarriers)) ~= length(pilotCarriers)
    error('pilotCarriers contains duplicates.');
end

usedCarriers = sort(usedCarriers);
pilotCarriers = sort(pilotCarriers);

validCarriers = -Nfft/2:(Nfft/2-1);
validCarriers(validCarriers == 0) = [];

if ~all(ismember(usedCarriers, validCarriers))
    error('Some entries in usedCarriers are outside the valid FFT-bin range.');
end

if any(usedCarriers == 0)
    error('usedCarriers must not include the DC carrier (0).');
end

if ~all(ismember(pilotCarriers, usedCarriers))
    error('pilotCarriers must be a subset of usedCarriers.');
end

dataCarriers = usedCarriers(~ismember(usedCarriers, pilotCarriers));

if isempty(dataCarriers)
    error('No data subcarriers remain after removing pilotCarriers from usedCarriers.');
end

if isempty(pilotCarriers)
    error('At least one pilot carrier is required for channel estimation.');
end

%% Modulation setup and exact constellation-power sanity check
% Normalization factors are chosen so the average constellation power is 1:
% QPSK   : divide by sqrt(2)
% 16-QAM : divide by sqrt(10)
% 64-QAM : divide by sqrt(42)
switch modType
    case 'QPSK'
        bitsPerSym = 2;
        normFactor = sqrt(2);
        sanityConst = [1+1j; 1-1j; -1+1j; -1-1j] / normFactor;

    case '16QAM'
        bitsPerSym = 4;
        normFactor = sqrt(10);
        axisLevels = [-3 -1 1 3];
        [II, QQ] = meshgrid(axisLevels, axisLevels);
        sanityConst = (II(:) + 1j*QQ(:)) / normFactor;

    case '64QAM'
        bitsPerSym = 6;
        normFactor = sqrt(42);
        axisLevels = [-7 -5 -3 -1 1 3 5 7];
        [II, QQ] = meshgrid(axisLevels, axisLevels);
        sanityConst = (II(:) + 1j*QQ(:)) / normFactor;

    otherwise
        error('Unsupported modulation. Use ''QPSK'', ''16QAM'', or ''64QAM''.');
end

if abs(mean(abs(sanityConst).^2) - 1) > 1e-12
    error('Constellation normalization sanity check failed.');
end

%% Derived parameters
fftCenter = Nfft/2 + 1;                 % DC index after fftshift
usedBins = usedCarriers + fftCenter;    % MATLAB indices of all used bins
dataBins = dataCarriers + fftCenter;    % MATLAB indices of data bins
pilotBins = pilotCarriers + fftCenter;  % MATLAB indices of pilot bins

numUsedSubc = length(usedCarriers);
numDataSubc = length(dataCarriers);
numPilotSubc = length(pilotCarriers);

bitsPerOFDMSymbol = numDataSubc * bitsPerSym;
numDataSymbols = numDataSubc * numOFDMSymbols;
numTxBits = bitsPerOFDMSymbol * numOFDMSymbols;
ofdmSymbolLen = Nfft + cpLen;

% Choose a block interleaver size that exactly fits one OFDM payload block.
candidateRows = [16 12 8 6 4 3 2];
interRows = 1;
for k = 1:length(candidateRows)
    if mod(bitsPerOFDMSymbol, candidateRows(k)) == 0
        interRows = candidateRows(k);
        break;
    end
end
interCols = bitsPerOFDMSymbol / interRows;

% Fixed pilot symbols inserted at known positions in each OFDM symbol.
pilotPattern = ones(numPilotSubc,1);
pilotPattern(2:2:end) = -1;  % fixed BPSK pilot values
pilotMatrix = repmat(pilotPattern, 1, numOFDMSymbols);

fprintf('Modulation                         : %s\n', modType);
fprintf('Nfft / CP length                   : %d / %d\n', Nfft, cpLen);
fprintf('Used / Pilot / Data bins           : %d / %d / %d\n', numUsedSubc, numPilotSubc, numDataSubc);
fprintf('Bits per QAM symbol                : %d\n', bitsPerSym);
fprintf('Bits per OFDM payload symbol       : %d\n', bitsPerOFDMSymbol);
fprintf('Interleaver size                   : %d x %d (per HT-DATA OFDM symbol)\n', interRows, interCols);
fprintf('HT-STF repeats / HT-LTF / HT-SIG   : %d / %d / %d\n', numHTSTFRepeats, numHTLTFSymbols, numHTSIGSymbols);
fprintf('Total HT-DATA payload bits/packet  : %d\n\n', numTxBits);

%% Generate transmit payload bits
txBits = double(rand(numTxBits,1) > 0.5);

%% Scrambler: pseudo-random PN sequence from a 7-bit LFSR
% Polynomial used: $x^7 + x^4 + 1$
%
% Scrambling is modulo-2 addition:
%   scrambledBits = mod(txBits + pnSeq, 2)
%
% Descrambling uses the SAME pnSeq and the SAME modulo-2 addition because
% XOR is its own inverse.
pnSeq = zeros(numTxBits,1);
lfsr = [1 0 1 1 0 1 1];  % any non-zero seed is valid

for n = 1:numTxBits
    pnSeq(n) = lfsr(end);
    feedbackBit = mod(lfsr(7) + lfsr(4), 2);
    lfsr(2:end) = lfsr(1:end-1);
    lfsr(1) = feedbackBit;
end

scrambledBits = mod(txBits + pnSeq, 2);

%% Block interleaver
% Interleaving is performed per HT-DATA OFDM payload block.
%
% TX operation:
%   1) reshape one bit block to [interRows x interCols]
%   2) transpose
%   3) vectorize
%
% This spreads adjacent bits apart. The receiver applies the inverse
% reshape/transpose to deinterleave after QAM demapping.
scrambledBitBlocks = reshape(scrambledBits, bitsPerOFDMSymbol, numOFDMSymbols);
interleavedBitBlocks = zeros(size(scrambledBitBlocks));

for symIdx = 1:numOFDMSymbols
    temp = reshape(scrambledBitBlocks(:,symIdx), interRows, interCols);
    interleavedBitBlocks(:,symIdx) = reshape(temp.', [], 1);
end

interleavedBits = interleavedBitBlocks(:);

%% QAM mapping for HT-DATA (Gray coded, unit average constellation power)
txBitGroups = reshape(interleavedBits, bitsPerSym, []).';

switch modType
    case 'QPSK'
        % Gray-coded QPSK: each axis carries 1 bit.
        I = 1 - 2*txBitGroups(:,1);
        Q = 1 - 2*txBitGroups(:,2);
        txDataSymbols = (I + 1j*Q) / sqrt(2);

    case '16QAM'
        % 2-bit Gray mapping on each axis:
        % 00 -> -3, 01 -> -1, 11 -> +1, 10 -> +3
        bI1 = txBitGroups(:,1);
        bI2 = txBitGroups(:,2);
        bQ1 = txBitGroups(:,3);
        bQ2 = txBitGroups(:,4);

        I = -(1 - 2*bI1) .* (3 - 2*bI2);
        Q = -(1 - 2*bQ1) .* (3 - 2*bQ2);
        txDataSymbols = (I + 1j*Q) / sqrt(10);

    case '64QAM'
        % 3-bit Gray mapping on each axis.
        % Monotonic levels [-7 -5 -3 -1 1 3 5 7] use Gray labels:
        % [000 001 011 010 110 111 101 100]
        %
        % The table below is indexed by the natural binary value of [b1 b2 b3].
        gray8_binaryOrder = [-7 -5 -1 -3 7 5 1 3];

        idxI = txBitGroups(:,1)*4 + txBitGroups(:,2)*2 + txBitGroups(:,3) + 1;
        idxQ = txBitGroups(:,4)*4 + txBitGroups(:,5)*2 + txBitGroups(:,6) + 1;

        I = gray8_binaryOrder(idxI);
        Q = gray8_binaryOrder(idxQ);
        txDataSymbols = (I + 1j*Q) / sqrt(42);
end

fprintf('Measured average mapped-symbol power = %.4f\n', mean(abs(txDataSymbols).^2));

%% ========================= Generate HT-STF ==============================
% HT-STF is a known repetitive short training field used for AGC and coarse
% timing. A sparse frequency-domain pattern on every 4th carrier creates a
% periodic time-domain waveform with period $Nfft/4$.

stfCarriers = usedCarriers(mod(usedCarriers,4) == 0);
if isempty(stfCarriers)
    error('HT-STF generation failed: no used carriers are multiples of 4.');
end
stfBins = stfCarriers + fftCenter;

htStfGrid = zeros(Nfft,1);

% Constant-modulus QPSK-like phases help keep the pattern well behaved.
stfPhasePattern = [1+1j; 1-1j; -1+1j; -1-1j] / sqrt(2);
htStfPattern = repmat(stfPhasePattern, ceil(length(stfBins)/length(stfPhasePattern)), 1);
htStfGrid(stfBins) = htStfPattern(1:length(stfBins));

htStfBase = ifft(ifftshift(htStfGrid, 1), Nfft, 1);

% Because only every 4th subcarrier is active, the waveform repeats every
% $Nfft/4$ samples. Extract one short period and repeat it.
htStfShortLen = Nfft/4;
htStfShort = htStfBase(1:htStfShortLen);

% Normalize the repeated short symbol to unit average power.
htStfShortPower = mean(abs(htStfShort).^2);
if htStfShortPower <= 0
    error('HT-STF normalization failed: zero or negative power.');
end
htStfShort = htStfShort / sqrt(htStfShortPower);

if abs(mean(abs(htStfShort).^2) - 1) > 1e-10
    error('HT-STF normalization sanity check failed.');
end

% Repeat the short symbol multiple times.
ht_stf = repmat(htStfShort, numHTSTFRepeats, 1);

%% ========================= Generate HT-LTF ==============================
% HT-LTF carries known symbols for channel estimation.
% Here it is implemented as one or more OFDM symbols with known BPSK values
% on the FULL set of used subcarriers and no payload data.

htLtfPattern = ones(numUsedSubc,1);
htLtfPattern(2:2:end) = -1;  % deterministic known BPSK pattern

htLtfGrid = zeros(Nfft, numHTLTFSymbols);
htLtfGrid(usedBins, :) = repmat(htLtfPattern, 1, numHTLTFSymbols);

htLtfTime = ifft(ifftshift(htLtfGrid, 1), Nfft, 1);

% Normalize EACH OFDM symbol to unit average power BEFORE CP.
htLtfPow = mean(abs(htLtfTime).^2, 1);
if any(htLtfPow <= 0)
    error('HT-LTF normalization failed: zero or negative power.');
end
htLtfTime = htLtfTime ./ repmat(sqrt(htLtfPow), Nfft, 1);

if max(abs(mean(abs(htLtfTime).^2, 1) - 1)) > 1e-10
    error('HT-LTF normalization sanity check failed.');
end

if cpLen > 0
    htLtfWithCP = [htLtfTime(end-cpLen+1:end, :); htLtfTime];
else
    htLtfWithCP = htLtfTime;
end

ht_ltf = htLtfWithCP(:);

%% ========================= Generate HT-SIG ==============================
% HT-SIG carries simplified control information such as modulation and
% packet length. This is NOT bit-exact to the standard, but the structure
% is present and BPSK-modulated on OFDM data subcarriers with pilots placed
% exactly like the payload OFDM symbols.

switch modType
    case 'QPSK'
        htSigModBits = [0; 0];
    case '16QAM'
        htSigModBits = [0; 1];
    case '64QAM'
        htSigModBits = [1; 0];
end

payloadLengthBytes = ceil(numTxBits/8);
payloadLengthBytesClipped = min(payloadLengthBytes, 2^16 - 1);

htSigLengthStr = dec2bin(payloadLengthBytesClipped, 16);
htSigLengthBits = double(htSigLengthStr(:) - '0');

htSigBpscStr = dec2bin(min(bitsPerSym, 2^4 - 1), 4);
htSigBpscBits = double(htSigBpscStr(:) - '0');

htSigNsymStr = dec2bin(min(numOFDMSymbols, 2^12 - 1), 12);
htSigNsymBits = double(htSigNsymStr(:) - '0');

% Simplified control information payload.
htSigCoreBits = [htSigModBits; htSigLengthBits; htSigBpscBits; htSigNsymBits];
htSigParityBit = mod(sum(htSigCoreBits), 2);
htSigPayloadBits = [htSigCoreBits; htSigParityBit];

% Fill the available HT-SIG data subcarriers across all HT-SIG symbols.
numHTSigDataBits = numDataSubc * numHTSIGSymbols;
htSigBits = zeros(numHTSigDataBits, 1);
copyLen = min(length(htSigPayloadBits), numHTSigDataBits);
htSigBits(1:copyLen) = htSigPayloadBits(1:copyLen);

if copyLen < numHTSigDataBits
    fillerPattern = repmat([0; 1; 1; 0], ceil((numHTSigDataBits - copyLen)/4), 1);
    htSigBits(copyLen+1:end) = fillerPattern(1:(numHTSigDataBits - copyLen));
end

% BPSK mapping: 0 -> +1, 1 -> -1
htSigDataSymbols = 1 - 2*htSigBits;

htSigGrid = zeros(Nfft, numHTSIGSymbols);
htSigGrid(dataBins, :) = reshape(htSigDataSymbols, numDataSubc, numHTSIGSymbols);
htSigGrid(pilotBins, :) = repmat(pilotPattern, 1, numHTSIGSymbols);

htSigTime = ifft(ifftshift(htSigGrid, 1), Nfft, 1);

% Normalize EACH OFDM symbol to unit average power BEFORE CP.
htSigPow = mean(abs(htSigTime).^2, 1);
if any(htSigPow <= 0)
    error('HT-SIG normalization failed: zero or negative power.');
end
htSigTime = htSigTime ./ repmat(sqrt(htSigPow), Nfft, 1);

if max(abs(mean(abs(htSigTime).^2, 1) - 1)) > 1e-10
    error('HT-SIG normalization sanity check failed.');
end

if cpLen > 0
    htSigWithCP = [htSigTime(end-cpLen+1:end, :); htSigTime];
else
    htSigWithCP = htSigTime;
end

ht_sig = htSigWithCP(:);

%% ========================= Generate HT-DATA =============================
% HT-DATA reuses the EXISTING OFDM payload chain:
% bits -> scrambler -> interleaver -> QAM -> pilot insertion -> IFFT -> CP

txDataMatrix = reshape(txDataSymbols, numDataSubc, numOFDMSymbols);

% Frequency-domain OFDM grid kept in fftshift order so carrier numbering is
% intuitive around DC. Unused bins remain zero.
txGrid = zeros(Nfft, numOFDMSymbols);
txGrid(dataBins, :) = txDataMatrix;
txGrid(pilotBins, :) = pilotMatrix;

htDataTime = ifft(ifftshift(txGrid, 1), Nfft, 1);

% IMPORTANT POWER NORMALIZATION:
% Normalize EACH HT-DATA OFDM symbol to unit average power BEFORE CP.
% This keeps packet fields comparable in power and makes the packet-level
% sample SNR definition consistent across modulation choices.
htDataPow = mean(abs(htDataTime).^2, 1);
if any(htDataPow <= 0)
    error('HT-DATA normalization failed: zero or negative power.');
end
htDataTime = htDataTime ./ repmat(sqrt(htDataPow), Nfft, 1);

if max(abs(mean(abs(htDataTime).^2, 1) - 1)) > 1e-10
    error('HT-DATA normalization sanity check failed.');
end

if cpLen > 0
    htDataWithCP = [htDataTime(end-cpLen+1:end, :); htDataTime];
else
    htDataWithCP = htDataTime;
end

% Serialize the HT-DATA field for packet assembly.
ht_data_serial = htDataWithCP(:);

%% =========================== Packet assembly ============================
% Concatenate all time-domain fields in the correct HT packet order.
txPacket = [ht_stf ; ht_ltf ; ht_sig ; ht_data_serial];
txPacket = txPacket(:);

% Save field lengths and packet indices.
lenHTSTF  = length(ht_stf);
lenHTLTF  = length(ht_ltf);
lenHTSIG  = length(ht_sig);
lenHTDATA = length(ht_data_serial);

idxHtDataStart = lenHTSTF + lenHTLTF + lenHTSIG + 1;
idxHtDataEnd = idxHtDataStart + lenHTDATA - 1;

if idxHtDataEnd ~= length(txPacket)
    error('Packet assembly indexing mismatch.');
end

% Field and packet power checks.
pwrHTSTF = mean(abs(ht_stf).^2);
pwrHTLTF = mean(abs(ht_ltf).^2);
pwrHTSIG = mean(abs(ht_sig).^2);
pwrHTDATA = mean(abs(ht_data_serial).^2);
txPacketPower = mean(abs(txPacket).^2);

fprintf('HT field sample lengths and average powers:\n');
fprintf('  HT-STF   : %6d samples, avg power = %.4f\n', lenHTSTF,  pwrHTSTF);
fprintf('  HT-LTF   : %6d samples, avg power = %.4f\n', lenHTLTF,  pwrHTLTF);
fprintf('  HT-SIG   : %6d samples, avg power = %.4f\n', lenHTSIG,  pwrHTSIG);
fprintf('  HT-DATA  : %6d samples, avg power = %.4f\n', lenHTDATA, pwrHTDATA);
fprintf('  Packet   : %6d samples, avg power = %.4f\n', length(txPacket), txPacketPower);
fprintf('  HT-DATA starts at packet sample index %d\n\n', idxHtDataStart);

%% ============================ BER simulation ============================
% BER is computed ONLY on the HT-DATA payload bits for now.

ber = zeros(size(snrDbVec));

% First run a no-noise sanity check. Then sweep SNR.
for snrIdx = 0:length(snrDbVec)

    if snrIdx == 0
        currentSnrDb = Inf;
    else
        currentSnrDb = snrDbVec(snrIdx);
    end

    %% Channel: AWGN
    % SNR is defined per complex received time-domain sample:
    %   $SNR = P_{packet} / P_{noise}$
    % where $P_{packet}$ is measured over the WHOLE transmitted packet.
    switch channelType
        case 'AWGN'
            if isinf(currentSnrDb)
                noise = zeros(size(txPacket));
            else
                noiseVar = txPacketPower / (10^(currentSnrDb/10));
                noise = sqrt(noiseVar/2) * (randn(size(txPacket)) + 1j*randn(size(txPacket)));
            end
            rxPacket = txPacket + noise;

        otherwise
            error('Unsupported channel type. Only ''AWGN'' is implemented.');
    end

    %% Receiver: ideal extraction of HT-DATA from the received packet
    rxHtDataSerial = rxPacket(idxHtDataStart:idxHtDataEnd);

    %% Remove CP and take FFT
    rxWithCP = reshape(rxHtDataSerial, ofdmSymbolLen, numOFDMSymbols);
    rxNoCP = rxWithCP(cpLen+1:end, :);
    rxGrid = fftshift(fft(rxNoCP, Nfft, 1), 1);

    %% Subcarrier extraction and pilot-based channel estimation
    % Known pilots estimate the per-OFDM-symbol complex gain:
    %   Hhat_pilot = received_pilot / known_pilot
    %
    % Because only an AWGN channel is implemented here, the true channel is
    % flat and equal to 1. The per-symbol mean pilot ratio therefore also
    % captures the TX-side symbol normalization factor applied before CP.
    %
    % If a frequency-selective channel is added later, this section can be
    % extended by estimating the channel from HT-LTF and/or interpolating
    % pilot-based estimates across frequency.
    rxPilots = rxGrid(pilotBins, :);
    rxData = rxGrid(dataBins, :);

    HhatPilots = rxPilots ./ pilotMatrix;
    Hhat = mean(HhatPilots, 1);                          % 1 x numOFDMSymbols
    rxDataEq = rxData ./ repmat(Hhat, numDataSubc, 1);  % equalize all data bins

    %% QAM demapper (hard decision)
    % After demapping we are back in the bit domain, so the receiver can
    % deinterleave and then descramble, which reverses the TX bit processing.
    rxDataSymbols = rxDataEq(:);

    switch modType
        case 'QPSK'
            rxBitGroups = zeros(numDataSymbols, 2);
            rxBitGroups(:,1) = real(rxDataSymbols) < 0;
            rxBitGroups(:,2) = imag(rxDataSymbols) < 0;

        case '16QAM'
            % Scale back to raw PAM levels before nearest-level decisions.
            rxScaled = rxDataSymbols * sqrt(10);
            levels = [-3 -1 1 3];
            bitsAxis = [0 0; 0 1; 1 1; 1 0];

            distI = abs(repmat(real(rxScaled), 1, length(levels)) - repmat(levels, numDataSymbols, 1));
            [~, idxIhat] = min(distI, [], 2);

            distQ = abs(repmat(imag(rxScaled), 1, length(levels)) - repmat(levels, numDataSymbols, 1));
            [~, idxQhat] = min(distQ, [], 2);

            rxBitGroups = [bitsAxis(idxIhat,:), bitsAxis(idxQhat,:)];

        case '64QAM'
            % Scale back to raw PAM levels before nearest-level decisions.
            rxScaled = rxDataSymbols * sqrt(42);
            levels = [-7 -5 -3 -1 1 3 5 7];
            bitsAxis = [0 0 0; ...
                        0 0 1; ...
                        0 1 1; ...
                        0 1 0; ...
                        1 1 0; ...
                        1 1 1; ...
                        1 0 1; ...
                        1 0 0];

            distI = abs(repmat(real(rxScaled), 1, length(levels)) - repmat(levels, numDataSymbols, 1));
            [~, idxIhat] = min(distI, [], 2);

            distQ = abs(repmat(imag(rxScaled), 1, length(levels)) - repmat(levels, numDataSymbols, 1));
            [~, idxQhat] = min(distQ, [], 2);

            rxBitGroups = [bitsAxis(idxIhat,:), bitsAxis(idxQhat,:)];
    end

    rxInterleavedBits = reshape(rxBitGroups.', [], 1);

    %% Deinterleaver
    % Inverse of the TX block interleaver:
    % RX block -> reshape to [interCols x interRows] -> transpose -> vectorize
    rxInterleavedBitBlocks = reshape(rxInterleavedBits, bitsPerOFDMSymbol, numOFDMSymbols);
    rxDeinterleavedBitBlocks = zeros(size(rxInterleavedBitBlocks));

    for symIdx = 1:numOFDMSymbols
        temp = reshape(rxInterleavedBitBlocks(:,symIdx), interCols, interRows).';
        rxDeinterleavedBitBlocks(:,symIdx) = temp(:);
    end

    rxScrambledBits = rxDeinterleavedBitBlocks(:);

    %% Descrambler
    % Same PN sequence and same modulo-2 addition as the scrambler.
    rxBits = mod(rxScrambledBits + pnSeq, 2);

    %% BER / sanity check
    % BER compares original HT-DATA payload bits to the recovered bits.
    bitErrors = sum(rxBits ~= txBits);

    if snrIdx == 0
        if bitErrors ~= 0
            error('No-noise sanity check FAILED: %d bit errors found.', bitErrors);
        else
            fprintf('No-noise sanity check PASSED: 0 bit errors.\n\n');
        end
    else
        ber(snrIdx) = bitErrors / numTxBits;
        fprintf('SNR = %5.1f dB, BER = %12.6e, bit errors = %d / %d\n', ...
                currentSnrDb, ber(snrIdx), bitErrors, numTxBits);
    end
end

%% Plot BER vs SNR
% Exact zero BER cannot be shown on a log axis, so for plotting only,
% replace zeros by half an error event. The printed BER values remain exact.
berForPlot = ber;
berForPlot(berForPlot == 0) = 0.5 / numTxBits;

figure;
semilogy(snrDbVec, berForPlot, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Packet sample SNR (dB)');
ylabel('Bit Error Rate (BER)');
title(sprintf(['Simplified 802.11n HT packet OFDM BER vs sample SNR ' ...
               '(%s, Nfft=%d, CP=%d, Data=%d, Pilots=%d)'], ...
      modType, Nfft, cpLen, numDataSubc, numPilotSubc));
```
