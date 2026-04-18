# 📡 Base-MATLAB OFDM Transceiver with FEC

A complete, first-principles implementation of an Orthogonal Frequency-Division Multiplexing (OFDM) transceiver. This script models a full physical layer (PHY) communication link in pure MATLAB. It does not require the Communications Toolbox or any other external add-ons.

-----

## 🛠️ System Architecture

The simulation models a complete transmit (TX) and receive (RX) chain over an Additive White Gaussian Noise (AWGN) channel.

-----

## 📡 HT Packet Structure

The transmitter generates a simplified IEEE 802.11n High Throughput (HT) packet. The time-domain waveform is formed by concatenating:

**HT-STF → HT-LTF → HT-SIG → HT-DATA**

-----

## 🧱 HT Block Descriptions

### HT-STF (Short Training Field)

**Purpose:**

  * AGC convergence
  * Coarse timing synchronization

**Generation:**

  * A sparse frequency-domain pattern is placed on every 4th used subcarrier.
  * IFFT produces a periodic time-domain waveform.
  * A segment of length $N_{fft}/4$ is extracted and repeated.
  * No cyclic prefix is used.

### HT-LTF (Long Training Field)

**Purpose:**

  * Channel estimation

**Generation:**

  * Known BPSK symbols are mapped onto all used subcarriers.
  * IFFT is applied.
  * Each OFDM symbol is normalized.
  * Cyclic prefix (CP) is added.

### HT-SIG (Signal Field)

**Purpose:**

  * Carries control information (modulation, payload length, symbol count).

**Generation:**

  * Control bits are constructed and BPSK-mapped onto data subcarriers.
  * Pilot symbols are inserted at fixed pilot locations.
  * IFFT is applied.
  * Each symbol is normalized.
  * CP is added.
    *(Note: The current receiver does not decode HT-SIG; it acts as a placeholder).*

### HT-DATA (Payload)

**Purpose:**

  * Payload transmission

**Generation:**

  * `bits → scrambler → interleaver → QAM mapper → pilot insertion → IFFT → normalization → CP`
  * Supports QPSK, 16-QAM, and 64-QAM.
  * Gray-coded constellations.
  * Fixed BPSK pilot symbols.

-----

## ⚡ Power Normalization

All OFDM symbols are normalized after IFFT and before CP insertion to ensure consistent average power across all packet fields. This is applied independently to HT-LTF, HT-SIG, and each HT-DATA symbol (HT-STF is normalized directly in the time domain):

$$P = \frac{1}{N} \sum |x[n]|^2 = 1$$

-----

## 📶 SNR Definition

SNR is defined per received complex time-domain sample:

$$SNR = \frac{P_{packet}}{P_{noise}}$$

Where $P_{packet}$ is the average power over the entire transmitted packet (HT-STF + HT-LTF + HT-SIG + HT-DATA + all CPs).

AWGN is added as:

$$n \sim \mathcal{CN}(0, \sigma^2)$$

With the noise variance defined as:

$$\sigma^2 = \frac{P_{packet}}{10^{SNR/10}}$$

### ⚠️ Important Note

  * CP and pilot symbols consume transmit energy but do not carry new payload bits.
  * Therefore, the effective $E_b/N_0$ is lower than the plotted sample SNR.
  * BER is computed only on HT-DATA payload bits.

-----

## 🔁 TX/RX Processing Chains

### Transmit Chain (TX)

1.  Data generation
2.  Scrambling (LFSR)
3.  Block interleaving
4.  QAM mapping (Gray-coded, unit power)
5.  Pilot insertion
6.  OFDM modulation (IFFT)
7.  Power normalization
8.  Cyclic prefix (CP)
9.  HT packet assembly

### Receive Chain (RX)

1.  Ideal HT-DATA extraction (no sync yet)
2.  CP removal & FFT
3.  Pilot-based channel estimation (per-symbol)
4.  Equalization (Zero-Forcing)
5.  QAM demapping (hard decision)
6.  Deinterleaving
7.  Descrambling
8.  BER computation

-----

## 🔬 Mathematical & Technical Specifications

  * **Scrambler (PN Sequence):** Polynomial $x^7 + x^4 + 1$
  * **OFDM Parameters:**
      * FFT Size: 64
      * Active Subcarriers: 52 (48 data, 4 pilots)
      * Subcarrier indices: `[-26:-1, 1:26]` (DC excluded)
  * **Channel Model:**
      * AWGN only
      * No multipath, CFO, or timing offset (ideal conditions)

-----

## 📊 Simulation Results

The Bit Error Rate (BER) is evaluated versus SNR. Higher-order modulations increase spectral efficiency but require higher SNR for the same BER.

| QPSK Performance | 16-QAM Performance | 64-QAM Performance |
| :---: | :---: | :---: |
| \<img src="QPSK\_BER.png" width="300" alt="QPSK BER"\> | \<img src="16QAM\_BER.png" width="300" alt="16-QAM BER"\> | \<img src="64QAM\_BER.png" width="300" alt="64-QAM BER"\> |

-----

## 🚀 Usage Guide

### Running the Simulation

Configure parameters at the top of the script:

```matlab
modType = '16QAM';      % Options: 'QPSK', '16QAM', '64QAM'
Nfft = 64;              % FFT size
cpLen = 16;             % Cyclic prefix length
numOFDMSymbols = 500;   % Number of HT-DATA symbols
snrDbVec = 0:2:24;      % SNR sweep array
```
