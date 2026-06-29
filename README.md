# OFDM_Ege — Simplified IEEE 802.11n HT Packet (Base MATLAB)

A from-scratch **IEEE 802.11n HT** OFDM transmitter/receiver written in **pure MATLAB** (no toolboxes). It builds a full HT packet, sends it through an AWGN channel, then recovers it with **real synchronization, channel estimation and Viterbi decoding** — and measures BER vs SNR. The scrambler, convolutional encoder, Viterbi decoder, CRC, packet detector and timing recovery are all hand-rolled.

**Signal chain:**
```
bits → scrambler → conv. code (BCC) + puncturing → interleaver → QAM → IFFT + CP
     → [ HT-STF | HT-LTF | HT-SIG | HT-DATA ] → AWGN
     → STF detect → LTF timing + channel est → equalize → QAM demap
     → deinterleave → depuncture → Viterbi → descramble → BER
```

Run **`802.11_V4.m`** in MATLAB. Compare against theory with **`compare_theory.m`**.

---

## What's inside

| Block | Implementation |
|---|---|
| **Modulation** | QPSK / 16-QAM / 64-QAM, Gray-coded |
| **FEC** | Convolutional K=7 (g0=133, g1=171), rates 1/2, 2/3, 3/4 + hard-decision Viterbi |
| **Preamble** | HT-STF (detect/timing), HT-LTF (channel est.), HT-SIG (CRC-checked MCS/length) |
| **Sync** | STF autocorrelation detection + LTF matched-filter fine timing |
| **Channel est.** | LS estimate from HT-LTF + pilot phase tracking |
| **Measurement** | BER vs SNR, Monte-Carlo averaged over many packets |

Packet order: **`HT-STF → HT-LTF → HT-SIG → HT-DATA`** (common IFFT scaling, so the LTF channel estimate applies to every field).

---

## Results

Coded BER vs packet SNR (rate 1/2, AWGN, averaged over 15 packets/point):

| QPSK | 16-QAM | 64-QAM |
|:---:|:---:|:---:|
| ![QPSK](QPSK_BER.png) | ![16-QAM](16QAM_BER.png) | ![64-QAM](64QAM_BER.png) |

Error-free around **~7 dB (QPSK)**, **~13 dB (16-QAM)**, **~19 dB (64-QAM)** — higher-order modulation = more throughput, but needs more SNR.

### vs. theory

Simulated **coded** BER overlaid on textbook **uncoded** AWGN curves (`compare_theory.m`):

![Sim vs theory](compare_BER.png)

- Modulation ordering matches theory exactly → the chain is correct.
- Coded curves cross the uncoded ones: a bit worse at very low SNR, then a steep waterfall with **~3–4 dB coding gain** — the classic hard-decision convolutional-code behaviour.

---

## Configure & run

Edit the parameters at the top of `802.11_V4.m`, then run it (no toolboxes needed):

```matlab
modType   = '16QAM';   % 'QPSK', '16QAM', '64QAM'
codeRate  = '1/2';     % '1/2', '2/3', '3/4'
Nfft = 64;  cpLen = 16;
numOFDMSymbols = 100;  % OFDM symbols per packet
numPackets     = 20;   % packets averaged per SNR point
snrDbVec  = 0:2:24;
```

A no-noise sanity check (0 errors + HT-SIG CRC + exact timing) runs before the sweep.

---

## Simplified / not yet implemented

- Block interleaver (not the standard two-permutation 802.11 interleaver)
- Hard-decision Viterbi only (no soft-decision / LLR)
- HT-SIG structure is close to but not bit-exact with the standard
- AWGN only (no CFO / multipath), single spatial stream
