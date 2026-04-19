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

## 🟢 HT-STF (High Throughput Short Training Field)

**Purpose:**
- Automatic Gain Control (AGC) convergence  
- Packet detection  
- Coarse timing synchronization  
- Carrier Frequency Offset (CFO) estimation  

**Structure:**
- Constructed from a predefined frequency-domain sequence  
- Designed to produce short periodic repetitions (~16 samples in time domain)  
- No cyclic prefix (CP)  

---

## 🔵 HT-LTF (High Throughput Long Training Field)

**Purpose:**
- Channel estimation (including MIMO channel estimation)  

**Structure:**
- Known BPSK-modulated sequence mapped onto subcarriers  
- One or more OFDM symbols depending on the number of spatial streams  
- Cyclic prefix (CP) included  

---

## 🟣 HT-SIG (High Throughput Signal Field)

**Purpose:**
- Carries PHY control information:
  - Modulation and Coding Scheme (MCS)  
  - Payload length  
  - Other transmission parameters  

**Structure:**
- Composed of **two OFDM symbols**  
- BPSK modulation with **rate-1/2 convolutional coding**  
- Standard OFDM processing chain:
  - Subcarrier mapping  
  - Pilot insertion  
  - IFFT  
  - Cyclic prefix (CP)  

---

## 🔴 HT-DATA (Payload)

**Purpose:**
- Transmission of user data  

**Processing Chain:**
bits → scrambler → FEC encoder → interleaver → QAM mapper → pilot insertion → IFFT → CP

**Details:**
- Supports modulation schemes:
  - QPSK  
  - 16-QAM  
  - 64-QAM  
- Uses Gray-coded constellations  
- Pilot symbols:
  - BPSK-modulated  
  - Follow a predefined polarity sequence (time-varying)  

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
