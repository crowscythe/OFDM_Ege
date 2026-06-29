%==========================================================================
% Simplified IEEE 802.11n HT-packet OFDM link in base MATLAB (no toolboxes)
%
% TX : info bits -> scrambler -> BCC (K=7) -> puncturing -> interleaver
%      -> QAM map -> pilot insertion -> IFFT -> CP
% Packet : [ HT-STF | HT-LTF | HT-SIG | HT-DATA ]  (common IFFT scaling)
% RX : STF autocorrelation packet detect -> LTF matched-filter fine timing
%      -> LTF channel estimation -> per-subcarrier equalization
%      -> HT-SIG decode (CRC-checked MCS/length)
%      -> QAM demap -> deinterleave -> depuncture -> Viterbi -> descramble
%
% Everything is hand-rolled (encoder, Viterbi, sync, CRC) so the script runs
% on a bare MATLAB install. Parts that are intentionally simplified relative
% to the standard are flagged with "SIMPLIFIED" comments and listed in the
% README.
%==========================================================================

clear; clc; close all;
rng(1);  % repeatable data and noise

%% ========================= User parameters ==============================
modType   = '16QAM';            % 'QPSK', '16QAM', or '64QAM'
codeRate  = '1/2';              % '1/2', '2/3', or '3/4' (BCC + puncturing)
Nfft      = 64;                 % FFT/IFFT size
cpLen     = 16;                 % cyclic prefix length
numOFDMSymbols = 100;           % HT-DATA OFDM symbols per packet
numPackets = 20;                % packets averaged per SNR point (Monte-Carlo)
snrDbVec  = 0:2:24;             % SNR sweep in dB (packet sample SNR)
channelType = 'AWGN';           % currently implemented: 'AWGN'

% Zero-centered subcarrier numbering (64-point, 52 used carriers).
usedCarriers  = [-26:-1 1:26];
pilotCarriers = [-21 -7 7 21];

% Preamble / packet field configuration.
numHTSTFRepeats = 10;           % repeated short symbols (no CP)
numHTLTFSymbols = 2;            % HT-LTF OFDM symbols used for channel est.
numHTSIGSymbols = 2;            % HT-SIG OFDM symbols (BPSK), standard = 2

% Receiver sync configuration.
syncPadSamples  = 120;          % silence (noise-only) samples around packet
stfDetectThresh = 0.6;          % autocorrelation detection threshold
fineSearchWin   = 160;          % +/- search range for LTF fine timing (~STF len)

%% =======================================================================
modType     = upper(char(modType));
channelType = upper(char(channelType));

%% ----------------------------- validation ------------------------------
if mod(Nfft,4) ~= 0,            error('Nfft must be divisible by 4 (HT-STF).'); end
if cpLen < 0 || cpLen >= Nfft,  error('Need 0 <= cpLen < Nfft.'); end
if numHTSIGSymbols ~= 2,        error('This version assumes numHTSIGSymbols = 2.'); end

usedCarriers  = sort(usedCarriers(:).');
pilotCarriers = sort(pilotCarriers(:).');
validCarriers = -Nfft/2:(Nfft/2-1); validCarriers(validCarriers==0) = [];
if ~all(ismember(usedCarriers,validCarriers)), error('usedCarriers out of range.'); end
if ~all(ismember(pilotCarriers,usedCarriers)), error('pilotCarriers must be subset of usedCarriers.'); end
dataCarriers = usedCarriers(~ismember(usedCarriers,pilotCarriers));

switch codeRate
    case '1/2', rNum=1; rDen=2; punMask=logical([1 1]);
    case '2/3', rNum=2; rDen=3; punMask=logical([1 1 1 0]);
    case '3/4', rNum=3; rDen=4; punMask=logical([1 1 1 0 0 1]);
    otherwise, error('codeRate must be ''1/2'', ''2/3'' or ''3/4''.');
end
punPeriod = numel(punMask);      % mother-code bits per puncture period
punKeep   = sum(punMask);        % transmitted bits per puncture period

%% --------------------------- modulation setup --------------------------
switch modType
    case 'QPSK',  bitsPerSym=2; qamNorm=sqrt(2);
    case '16QAM', bitsPerSym=4; qamNorm=sqrt(10);
    case '64QAM', bitsPerSym=6; qamNorm=sqrt(42);
    otherwise, error('Unsupported modulation.');
end

%% --------------------------- derived params ----------------------------
fftCenter   = Nfft/2 + 1;
usedBins    = usedCarriers  + fftCenter;
dataBins    = dataCarriers  + fftCenter;
pilotBins   = pilotCarriers + fftCenter;
numUsedSubc = numel(usedCarriers);
numDataSubc = numel(dataCarriers);
numPilotSubc= numel(pilotCarriers);

ofdmSymbolLen = Nfft + cpLen;
Ncbps         = numDataSubc * bitsPerSym;       % coded bits per HT-DATA symbol
totalCoded    = Ncbps * numOFDMSymbols;         % coded bits over whole payload

% Size the BCC payload so that, after rate-1/2 encoding + puncturing, the
% transmitted coded stream exactly fills the available data subcarriers
% (any shortfall is zero-padded and stripped at the RX).
numTailBits = 6;
nPunPeriods = floor(totalCoded / punKeep);
codedLenExact = nPunPeriods * punKeep;          % real coded bits (<= totalCoded)
motherLen   = nPunPeriods * punPeriod;          % mother-code bits before puncturing
Lin         = motherLen / 2;                    % BCC input length (info + tail)
numInfoBits = Lin - numTailBits;
if numInfoBits < 8
    error('numOFDMSymbols too small for the chosen rate; increase it.');
end
codedPadLen = totalCoded - codedLenExact;       % zero padding to fill the grid

% Block interleaver geometry (SIMPLIFIED: not the standard two-permutation
% 802.11 interleaver; a plain row/column block interleaver per OFDM symbol).
candidateRows = [16 12 8 6 4 3 2];
interRows = 1;
for k = 1:numel(candidateRows)
    if mod(Ncbps,candidateRows(k))==0, interRows = candidateRows(k); break; end
end
interCols = Ncbps / interRows;

% Fixed BPSK pilots (SIMPLIFIED: static polarity, not the 127-length sign LFSR).
pilotPattern = ones(numPilotSubc,1); pilotPattern(2:2:end) = -1;

% Standard rate-1/2 BCC trellis (K=7, g0=133, g1=171 octal).
trellis = buildTrellis();

mcsIndex = mcsFromMode(modType, codeRate);

fprintf('Modulation / code rate    : %s , %s  (MCS %d)\n', modType, codeRate, mcsIndex);
fprintf('Nfft / CP                 : %d / %d\n', Nfft, cpLen);
fprintf('Used/Data/Pilot subc      : %d / %d / %d\n', numUsedSubc,numDataSubc,numPilotSubc);
fprintf('Coded bits/OFDM symbol    : %d\n', Ncbps);
fprintf('Payload info bits/packet  : %d (+%d tail), coded->%d, pad %d\n', ...
        numInfoBits, numTailBits, codedLenExact, codedPadLen);
fprintf('Interleaver (per symbol)  : %d x %d  [SIMPLIFIED block]\n', interRows, interCols);
fprintf('Preamble STF/LTF/SIG syms : %d reps / %d / %d\n\n', numHTSTFRepeats,numHTLTFSymbols,numHTSIGSymbols);

%% ========================= HT-DATA: TX setup ===========================
% The HT-DATA payload is (re)generated per packet inside the Monte-Carlo loop
% (see genHTData below). Here we only fix the things that stay constant across
% packets: the scrambler PN sequence and the common IFFT scale.
scramblerSeed = [1 0 1 1 0 1 1];
pnData = lfsrPN(numInfoBits, scramblerSeed);

%% ===================== common frequency->time scaling ==================
% A SINGLE scale is shared by HT-LTF, HT-SIG and HT-DATA so the channel
% estimate obtained from HT-LTF applies directly to the other fields. The
% value cancels on equalization, so we use its analytic expectation: the IFFT
% of numUsedSubc unit-power subcarriers has mean power numUsedSubc/Nfft^2.
commonScale = Nfft / sqrt(numUsedSubc);   % -> unit average time power

%% ============================ HT-STF ===================================
% SIMPLIFIED: sparse pattern on every 4th subcarrier -> period-(Nfft/4)
% time waveform, repeated. Used for AGC, packet detection, coarse timing.
stfCarriers = usedCarriers(mod(usedCarriers,4)==0);
stfBins = stfCarriers + fftCenter;
stfPhase = [1+1j;1-1j;-1+1j;-1-1j]/sqrt(2);
stfVals  = repmat(stfPhase, ceil(numel(stfBins)/4), 1);
htStfGrid = zeros(Nfft,1); htStfGrid(stfBins) = stfVals(1:numel(stfBins));
htStfBase = ifft(ifftshift(htStfGrid,1), Nfft, 1);
htStfShort = htStfBase(1:Nfft/4);
htStfShort = htStfShort / sqrt(mean(abs(htStfShort).^2));   % unit power
ht_stf = repmat(htStfShort, numHTSTFRepeats, 1);

%% ============================ HT-LTF ===================================
% Known BPSK pattern on all used subcarriers -> used for channel estimation.
htLtfFreq = ones(numUsedSubc,1); htLtfFreq(2:2:end) = -1;
htLtfGrid = zeros(Nfft,numHTLTFSymbols);
htLtfGrid(usedBins,:) = repmat(htLtfFreq,1,numHTLTFSymbols);
htLtfTime = ifft(ifftshift(htLtfGrid,1), Nfft, 1) * commonScale;
ht_ltf = addCP(htLtfTime, cpLen);
ht_ltf = ht_ltf(:);

%% ============================ HT-SIG ===================================
% Closer to the standard: field bits + CRC-8 + 6 tail bits, rate-1/2 BCC,
% block-interleaved, BPSK on the data subcarriers of 2 OFDM symbols.
% SIMPLIFIED: field layout and CRC polynomial are not bit-exact to 802.11.
sigCodedCap = numDataSubc * numHTSIGSymbols;     % BPSK -> 1 coded bit/carrier
sigLin      = sigCodedCap / 2;                   % rate-1/2 input length
sigCrcBits  = 8;
sigFieldLen = sigLin - numTailBits - sigCrcBits; % usable field bits

% Pack control fields (MCS, HT-Length in bytes), rest reserved/zero.
htLengthBytes = ceil(numInfoBits/8);
sigFields = zeros(sigFieldLen,1);
sigFields(1:7)  = de2biMSB(mcsIndex, 7);
sigFields(8:23) = de2biMSB(min(htLengthBytes,2^16-1), 16);
sigCrc   = crc8(sigFields);
sigInfo  = [sigFields; sigCrc];                  % CRC-protected payload
sigCoderIn  = [sigInfo; zeros(numTailBits,1)];
sigCoded    = bccEncode(sigCoderIn, trellis);    % rate 1/2 (no puncturing)
sigInter    = blockInterleave(sigCoded, 2, sigCodedCap/2);
sigSymsBPSK  = 1 - 2*sigInter;                   % BPSK: 0->+1, 1->-1
htSigGrid = zeros(Nfft,numHTSIGSymbols);
htSigGrid(dataBins,:)  = reshape(sigSymsBPSK, numDataSubc, numHTSIGSymbols);
htSigGrid(pilotBins,:) = repmat(pilotPattern,1,numHTSIGSymbols);
htSigTime = ifft(ifftshift(htSigGrid,1), Nfft, 1) * commonScale;
ht_sig = addCP(htSigTime, cpLen); ht_sig = ht_sig(:);

%% ===================== fixed preamble bookkeeping ======================
% HT-STF/LTF/SIG are identical for every packet (HT-DATA varies). Field sample
% lengths are constant, so the RX field offsets can be precomputed.
lenHTSTF=numel(ht_stf); lenHTLTF=numel(ht_ltf);
lenHTSIG=numel(ht_sig); lenHTDATA = ofdmSymbolLen*numOFDMSymbols;
trueStart = syncPadSamples + 1;        % silence is padded around each packet

% Parameter bundle passed to the per-packet TX/RX helpers.
P = struct('numInfoBits',numInfoBits,'numTailBits',numTailBits,'pnData',pnData, ...
    'trellis',trellis,'punMask',punMask,'punPeriod',punPeriod,'motherLen',motherLen, ...
    'codedPadLen',codedPadLen,'codedLenExact',codedLenExact,'Ncbps',Ncbps, ...
    'numOFDMSymbols',numOFDMSymbols,'interRows',interRows,'interCols',interCols, ...
    'modType',modType,'Nfft',Nfft,'cpLen',cpLen,'ofdmSymbolLen',ofdmSymbolLen, ...
    'commonScale',commonScale,'numDataSubc',numDataSubc,'dataBins',dataBins, ...
    'pilotBins',pilotBins,'pilotPattern',pilotPattern,'dataCarriers',dataCarriers, ...
    'pilotCarriers',pilotCarriers,'usedCarriers',usedCarriers,'lenHTDATA',lenHTDATA);

fprintf('Field samples  STF/LTF/SIG/DATA : %d / %d / %d / %d\n\n', ...
        lenHTSTF,lenHTLTF,lenHTSIG,lenHTDATA);

%% ===================== no-noise sanity check (1 packet) ================
% Build one noiseless packet and verify the full chain: exact sync, HT-SIG
% CRC pass, and zero payload bit errors.
[ht_data, infoBits] = genHTData(P);
rxStream = [zeros(syncPadSamples,1); ht_stf; ht_ltf; ht_sig; ht_data; zeros(syncPadSamples,1)];
[rxBits, startHat, Hhat] = receivePacket(rxStream, ht_ltf, htLtfFreq, ...
    usedBins, lenHTSTF, lenHTLTF, lenHTSIG, numHTLTFSymbols, stfDetectThresh, ...
    fineSearchWin, trueStart, P);
[sigOk, sigMcs, sigLenBytes] = decodeHTSIG(rxStream, startHat+lenHTSTF+lenHTLTF, ...
    cpLen, Nfft, dataBins, pilotBins, pilotPattern, Hhat, dataCarriers, usedCarriers, ...
    numDataSubc, numHTSIGSymbols, trellis, numTailBits, sigFieldLen, sigCrcBits);
if ~sigOk, error('HT-SIG CRC failed at infinite SNR.'); end
if sum(rxBits~=infoBits)~=0
    error('No-noise sanity FAILED: %d errors (sync off by %d).', ...
          sum(rxBits~=infoBits), startHat-trueStart);
end
fprintf('HT-SIG decoded : MCS=%d, length=%d bytes, CRC OK (true MCS=%d, %d bytes)\n', ...
        sigMcs, sigLenBytes, mcsIndex, htLengthBytes);
fprintf('No-noise sanity: PASSED (0 errors, start=%d/%d)\n\n', startHat, trueStart);

%% ============== BER simulation (Monte-Carlo averaged over packets) ======
ber = zeros(size(snrDbVec));
for snrIdx = 1:numel(snrDbVec)
    snrDb = snrDbVec(snrIdx);
    totErr = 0; totBits = 0;
    for pk = 1:numPackets
        % Fresh payload + fresh noise for each packet realization.
        [ht_data, infoBits] = genHTData(P);
        txPacket = [ht_stf; ht_ltf; ht_sig; ht_data];
        txPacketPower = mean(abs(txPacket).^2);
        txStream = [zeros(syncPadSamples,1); txPacket; zeros(syncPadSamples,1)];

        noiseVar = txPacketPower / (10^(snrDb/10));
        rxStream = txStream + ...
            sqrt(noiseVar/2)*(randn(size(txStream))+1j*randn(size(txStream)));

        rxBits = receivePacket(rxStream, ht_ltf, htLtfFreq, usedBins, ...
            lenHTSTF, lenHTLTF, lenHTSIG, numHTLTFSymbols, stfDetectThresh, ...
            fineSearchWin, trueStart, P);

        totErr  = totErr + sum(rxBits ~= infoBits);
        totBits = totBits + numInfoBits;
    end
    ber(snrIdx) = totErr / totBits;
    fprintf('SNR=%5.1f dB | BER=%11.4e | errors %d/%d (%d packets)\n', ...
            snrDb, ber(snrIdx), totErr, totBits, numPackets);
end

%% ============================== plot ===================================
berPlot = ber; berPlot(berPlot==0) = 0.5/(numInfoBits*numPackets);
figure;
semilogy(snrDbVec, berPlot, 'o-','LineWidth',1.5,'MarkerSize',6); grid on;
xlabel('Packet sample SNR (dB)'); ylabel('BER');
title(sprintf('802.11n HT OFDM, %s rate %s, Nfft=%d CP=%d (%d pkts/pt)', ...
      modType, codeRate, Nfft, cpLen, numPackets));

%% =======================================================================
%% ============================ local functions =========================
%% =======================================================================

function [ht_data, infoBits] = genHTData(P)
% Generate one HT-DATA field from fresh random payload bits:
% info -> scramble -> BCC -> puncture -> interleave -> QAM -> IFFT -> CP.
    infoBits = double(rand(P.numInfoBits,1) > 0.5);
    scrambled = mod(infoBits + P.pnData, 2);
    coderIn = [scrambled; zeros(P.numTailBits,1)];
    mother = bccEncode(coderIn, P.trellis);
    punc = mother(repmat(P.punMask(:), P.motherLen/P.punPeriod, 1));
    grid = [punc; zeros(P.codedPadLen,1)];
    blocks = reshape(grid, P.Ncbps, P.numOFDMSymbols);
    inter = zeros(size(blocks));
    for s = 1:P.numOFDMSymbols
        inter(:,s) = blockInterleave(blocks(:,s), P.interRows, P.interCols);
    end
    sym = qamMap(inter(:), P.modType);
    txGrid = zeros(P.Nfft, P.numOFDMSymbols);
    txGrid(P.dataBins,:)  = reshape(sym, P.numDataSubc, P.numOFDMSymbols);
    txGrid(P.pilotBins,:) = repmat(P.pilotPattern, 1, P.numOFDMSymbols);
    tt = ifft(ifftshift(txGrid,1), P.Nfft, 1) * P.commonScale;
    ht_data = reshape(addCP(tt, P.cpLen), [], 1);
end

function [rxBits, startHat, Hhat] = receivePacket(rxStream, ht_ltf, htLtfFreq, ...
        usedBins, lenHTSTF, lenHTLTF, lenHTSIG, numHTLTFSymbols, thresh, win, trueStart, P)
% Full HT-DATA receiver: STF detect -> LTF fine timing -> LTF channel est ->
% equalize + pilot phase track -> QAM demap -> deinterleave -> depuncture ->
% Viterbi -> descramble. Returns the recovered payload info bits.
    Nfft=P.Nfft; cpLen=P.cpLen;
    coarse = stfDetect(rxStream, Nfft/4, Nfft/4, thresh);
    if isempty(coarse), coarse = trueStart; end
    startHat = ltfFineTiming(rxStream, ht_ltf, coarse, lenHTSTF, win);

    % Safety clamp: keep the whole packet inside the stream even if a noisy
    % realization mis-locks the detector (such a packet just decodes to errors
    % instead of crashing the simulation).
    maxStart = numel(rxStream) - (lenHTSTF+lenHTLTF+lenHTSIG+P.lenHTDATA) + 1;
    startHat = min(max(startHat,1), maxStart);

    ltfStart  = startHat + lenHTSTF;
    dataStart = ltfStart + lenHTLTF + lenHTSIG;
    Hhat = ltfChannelEst(rxStream, ltfStart, cpLen, Nfft, usedBins, htLtfFreq, numHTLTFSymbols);

    rxData = rxStream(dataStart : dataStart + P.lenHTDATA - 1);
    rxNoCP = reshape(rxData, P.ofdmSymbolLen, P.numOFDMSymbols);
    rxNoCP = rxNoCP(cpLen+1:end,:);
    rxFreq = fftshift(fft(rxNoCP, Nfft, 1), 1);

    Hd = Hhat(ismember(P.usedCarriers,P.dataCarriers));  Hd=Hd(:);
    Hp = Hhat(ismember(P.usedCarriers,P.pilotCarriers)); Hp=Hp(:);
    eqD = rxFreq(P.dataBins,:)  ./ repmat(Hd, 1, P.numOFDMSymbols);
    eqP = rxFreq(P.pilotBins,:) ./ repmat(Hp, 1, P.numOFDMSymbols);
    phErr = angle(sum(eqP .* repmat(conj(P.pilotPattern),1,P.numOFDMSymbols),1));
    eqD = eqD .* repmat(exp(-1j*phErr), P.numDataSubc, 1);

    coded = qamDemap(eqD(:), P.modType);
    blocks = reshape(coded, P.Ncbps, P.numOFDMSymbols);
    deint = zeros(size(blocks));
    for s = 1:P.numOFDMSymbols
        deint(:,s) = blockDeinterleave(blocks(:,s), P.interRows, P.interCols);
    end
    grid = deint(:);
    mother = depuncture(grid(1:P.codedLenExact), P.punMask, P.motherLen);
    coderIn = viterbiDecode(mother, P.trellis);
    rxBits = mod(coderIn(1:P.numInfoBits) + P.pnData, 2);
end

function pn = lfsrPN(N, seed)
% Additive scrambler PN sequence, polynomial x^7 + x^4 + 1.
    pn = zeros(N,1); lfsr = seed(:).';
    for n = 1:N
        pn(n) = lfsr(end);
        fb = mod(lfsr(7)+lfsr(4),2);
        lfsr(2:end) = lfsr(1:end-1); lfsr(1) = fb;
    end
end

function t = buildTrellis()
% Rate-1/2, K=7 convolutional code, generators g0=133, g1=171 (octal).
% State = 6 memory bits [s1..s6], s1 = most recent. Output pair [A B].
    nS = 64;
    t.nextState = zeros(nS,2);   % next state for input 0/1
    t.outA = zeros(nS,2); t.outB = zeros(nS,2);
    for st = 0:nS-1
        s1=bitget(st,6); s2=bitget(st,5); s3=bitget(st,4);
        s4=bitget(st,3); s5=bitget(st,2); s6=bitget(st,1);
        for u = 0:1
            % g0 = [1 0 1 1 0 1 1] over [u s1 s2 s3 s4 s5 s6]
            A = mod(u + s2 + s3 + s5 + s6, 2);
            % g1 = [1 1 1 1 0 0 1]
            B = mod(u + s1 + s2 + s3 + s6, 2);
            ns = u*32 + s1*16 + s2*8 + s3*4 + s4*2 + s5;  % memory shift
            t.nextState(st+1,u+1) = ns;
            t.outA(st+1,u+1) = A; t.outB(st+1,u+1) = B;
        end
    end
    % Predecessor table for the Viterbi butterfly.
    t.predState = zeros(nS,2); t.predIn = zeros(nS,2);
    t.predA = zeros(nS,2); t.predB = zeros(nS,2);
    cnt = ones(nS,1);
    for st = 0:nS-1
        for u = 0:1
            ns = t.nextState(st+1,u+1);
            c = cnt(ns+1);
            t.predState(ns+1,c) = st; t.predIn(ns+1,c) = u;
            t.predA(ns+1,c) = t.outA(st+1,u+1);
            t.predB(ns+1,c) = t.outB(st+1,u+1);
            cnt(ns+1) = c+1;
        end
    end
end

function coded = bccEncode(inBits, t)
% Encode column bit vector; output [A0 B0 A1 B1 ...]. Start state 0.
    N = numel(inBits); coded = zeros(2*N,1); st = 0;
    for n = 1:N
        u = inBits(n);
        coded(2*n-1) = t.outA(st+1,u+1);
        coded(2*n)   = t.outB(st+1,u+1);
        st = t.nextState(st+1,u+1);
    end
end

function info = viterbiDecode(rxCoded, t)
% Hard-decision Viterbi for the terminated rate-1/2 code. rxCoded holds the
% mother-code bit pairs; punctured positions are marked -1 (erasure).
    nStages = numel(rxCoded)/2; nS = 64;
    INF = 1e9; metric = INF*ones(nS,1); metric(1) = 0;  % start state 0
    survPrev = zeros(nS,nStages,'uint8');
    survIn   = zeros(nS,nStages,'uint8');
    pA = t.predA; pB = t.predB; pS = t.predState; pI = t.predIn;
    for k = 1:nStages
        r1 = rxCoded(2*k-1); r2 = rxCoded(2*k);
        % branch metric (Hamming) for each next state's two predecessors
        bm1 = (r1~=-1).*(pA(:,1)~=r1) + (r2~=-1).*(pB(:,1)~=r2);
        bm2 = (r1~=-1).*(pA(:,2)~=r1) + (r2~=-1).*(pB(:,2)~=r2);
        c1 = metric(pS(:,1)+1) + bm1;
        c2 = metric(pS(:,2)+1) + bm2;
        useFirst = c1 <= c2;
        metric = min(c1,c2);
        survPrev(:,k) = uint8(pS(:,2)); survPrev(useFirst,k) = uint8(pS(useFirst,1));
        survIn(:,k)   = uint8(pI(:,2)); survIn(useFirst,k)   = uint8(pI(useFirst,1));
    end
    % traceback from terminated state 0
    info = zeros(nStages,1); st = 0;
    for k = nStages:-1:1
        info(k) = survIn(st+1,k);
        st = double(survPrev(st+1,k));
    end
end

function out = depuncture(rxBits, mask, motherLen)
% Insert erasures (-1) at punctured positions to rebuild the mother stream.
    out = -ones(motherLen,1);
    keepIdx = find(repmat(mask(:), motherLen/numel(mask), 1));
    out(keepIdx) = rxBits;
end

function y = blockInterleave(x, rows, cols)
    y = reshape(reshape(x, rows, cols).', [], 1);
end
function y = blockDeinterleave(x, rows, cols)
    y = reshape(reshape(x, cols, rows).', [], 1);
end

function sym = qamMap(bits, modType)
    g = reshape(bits, [], 1);
    switch modType
        case 'QPSK'
            b = reshape(g,2,[]).';
            I = 1-2*b(:,1); Q = 1-2*b(:,2);
            sym = (I+1j*Q)/sqrt(2);
        case '16QAM'
            b = reshape(g,4,[]).';
            I = -(1-2*b(:,1)).*(3-2*b(:,2));
            Q = -(1-2*b(:,3)).*(3-2*b(:,4));
            sym = (I+1j*Q)/sqrt(10);
        case '64QAM'
            b = reshape(g,6,[]).';
            gray = [-7 -5 -1 -3 7 5 1 3];
            iI = b(:,1)*4+b(:,2)*2+b(:,3)+1;
            iQ = b(:,4)*4+b(:,5)*2+b(:,6)+1;
            sym = (gray(iI).'+1j*gray(iQ).')/sqrt(42);
    end
end

function bits = qamDemap(sym, modType)
    sym = sym(:);
    switch modType
        case 'QPSK'
            b = [real(sym)<0, imag(sym)<0];
            bits = reshape(b.',[],1);
        case '16QAM'
            r = sym*sqrt(10); lvl=[-3 -1 1 3]; map=[0 0;0 1;1 1;1 0];
            bits = pamBits(r,lvl,map);
        case '64QAM'
            r = sym*sqrt(42); lvl=[-7 -5 -3 -1 1 3 5 7];
            map=[0 0 0;0 0 1;0 1 1;0 1 0;1 1 0;1 1 1;1 0 1;1 0 0];
            bits = pamBits(r,lvl,map);
    end
end
function bits = pamBits(r, lvl, map)
    n=numel(r); L=numel(lvl);
    [~,iI]=min(abs(repmat(real(r),1,L)-repmat(lvl,n,1)),[],2);
    [~,iQ]=min(abs(repmat(imag(r),1,L)-repmat(lvl,n,1)),[],2);
    grp=[map(iI,:) map(iQ,:)];
    bits=reshape(grp.',[],1);
end

function y = addCP(timeSyms, cpLen)
    if cpLen>0, y=[timeSyms(end-cpLen+1:end,:); timeSyms]; else, y=timeSyms; end
end

function start = stfDetect(stream, D, L, thresh)
% Delay-autocorrelation packet detector. Returns the first sample of the first
% sustained run (>= D samples) above the threshold, which marks the rising edge
% of the STF plateau. Requiring a sustained run rejects isolated noise spikes.
    N = numel(stream); start = [];
    M = zeros(N,1);
    for n = 1:N-2*D-L
        seg  = stream(n:n+L-1);
        segD = stream(n+D:n+D+L-1);
        P = sum(seg .* conj(segD));
        R = sum(abs(segD).^2);
        if R>0, M(n) = (abs(P)^2)/(R^2); end
    end
    above = M > thresh; runStart = []; runLen = 0;
    for n = 1:N
        if above(n)
            if runLen==0, runStart = n; end
            runLen = runLen + 1;
            if runLen >= D, start = runStart; return; end
        else
            runLen = 0;
        end
    end
end

function startHat = ltfFineTiming(stream, ltfTemplate, coarseStart, lenSTF, win)
% Matched-filter fine timing: slide the known HT-LTF over a window around the
% expected LTF position and pick the correlation peak; map back to packet start.
    expected = coarseStart + lenSTF;
    lo = max(1, expected-win); hi = expected+win;
    L = numel(ltfTemplate); best=-1; bestS=expected;
    tpl = ltfTemplate(:);
    for s = lo:hi
        if s+L-1 > numel(stream), break; end
        c = abs(sum(conj(tpl).*stream(s:s+L-1)));
        if c>best, best=c; bestS=s; end
    end
    startHat = bestS - lenSTF;
end

function Hhat = ltfChannelEst(stream, ltfStart, cpLen, Nfft, usedBins, knownFreq, nSym)
% Per-subcarrier LS channel estimate averaged over the HT-LTF symbols.
    acc = zeros(numel(usedBins),1);
    for m = 0:nSym-1
        base = ltfStart + m*(Nfft+cpLen) + cpLen;
        sym = stream(base : base+Nfft-1);
        F = fftshift(fft(sym,Nfft),1);
        acc = acc + F(usedBins)./knownFreq;
    end
    Hhat = acc / nSym;
end

function [ok, mcs, lenBytes] = decodeHTSIG(stream, sigStart, cpLen, Nfft, ...
        dataBins, pilotBins, pilotPattern, Hhat, dataCarriers, usedCarriers, ...
        numDataSubc, nSym, t, numTail, fieldLen, crcBits)
% Equalize, BPSK-demap, deinterleave, Viterbi-decode HT-SIG and check CRC.
    Hd = Hhat(ismember(usedCarriers,dataCarriers)); Hd=Hd(:);
    Hp = Hhat(~ismember(usedCarriers,dataCarriers)); Hp=Hp(:);
    eqData = zeros(numDataSubc,nSym);
    for m = 0:nSym-1
        base = sigStart + m*(Nfft+cpLen) + cpLen;
        sym = stream(base:base+Nfft-1);
        F = fftshift(fft(sym,Nfft),1);
        d = F(dataBins)./Hd;
        p = F(pilotBins)./Hp;
        ph = angle(sum(p.*conj(pilotPattern)));
        eqData(:,m+1) = d*exp(-1j*ph);
    end
    rxBits = double(real(eqData(:))<0);              % BPSK demap
    coded  = blockDeinterleave(rxBits, 2, numel(rxBits)/2);
    % Rate-1/2, no puncturing: the deinterleaved stream is already [A B A B...].
    decoded = viterbiDecode(coded, t);
    info = decoded(1:end-numTail);
    fields = info(1:fieldLen);
    crcRx  = info(fieldLen+1:fieldLen+crcBits);
    ok = isequal(crc8(fields), crcRx(:));
    mcs = biMSB2de(fields(1:7));
    lenBytes = biMSB2de(fields(8:23));
end

function c = crc8(bits)
% CRC-8 (poly x^8+x^2+x+1, init 0). SIMPLIFIED vs. the 802.11 HT-SIG CRC.
    reg = zeros(8,1); poly = [0 0 0 0 0 1 1 1]';  % x^2+x+1 feedback taps
    for i = 1:numel(bits)
        fb = mod(reg(1)+bits(i),2);
        reg(1:7) = reg(2:8);
        reg(8) = fb;
        if fb, reg = mod(reg + poly,2); end
    end
    c = reg;
end

function b = de2biMSB(val, n)
    b = zeros(n,1);
    for k = 1:n, b(k) = bitget(val, n-k+1); end
end
function v = biMSB2de(b)
    n=numel(b); v=0; for k=1:n, v=v+b(k)*2^(n-k); end
end

function m = mcsFromMode(modType, codeRate)
    key = [modType '-' codeRate];
    switch key
        case 'QPSK-1/2', m=1; case 'QPSK-3/4', m=2;
        case '16QAM-1/2',m=3; case '16QAM-3/4',m=4;
        case '64QAM-2/3',m=5; case '64QAM-3/4',m=6;
        otherwise, m=0;   % BPSK-1/2 / unspecified combos
    end
end
