# UAV-IRS-Communications

## Overview
This project consists of two MATLAB scripts that simulate and analyze RIS-aided wireless communication systems under Nakagami-m fading and Inverse Gamma shadowing. The simulations evaluate the impact of reflecting elements, phase shifts, and swarm-based modifications on system performance. The two key functionalities implemented in this project are:

1. **Performance comparison of RIS-aided systems** - Analyzes the spectral efficiency and outage probability for different system configurations.
2. **Swarm-based modification** - Introduces a swarm-based approach where multiple UAV-mounted RIS elements are utilized to enhance signal reliability.

## Functionalities
### 1. Performance comparison of RIS-aided Systems
This script simulates a standard RIS-aided communication system and evaluates:
- End-to-end channel fading and shadowing models
- Phase shift optimization for improved signal quality
- Outage probability analysis under varying spectral efficiency thresholds
- Analytical approximations of fading distributions

#### Key features:
- Uses Nakagami-m fading and Inverse Gamma shadowing to model channel impairments.
- Implements moment matching techniques to approximate the distribution of the received signal.
- Evaluates the impact of different phase-shift configurations on outage probability.
- Performs Laplace transform analysis for verifying analytical approximations.

### 2. Swarm-based modification
This script extends the original RIS model by incorporating a swarm-based approach where multiple UAV-mounted RIS elements work collaboratively to enhance transmission reliability.

#### Key features:
- Introduces a swarm of UAVs, each carrying an IRS with a predefined number of reflecting elements.
- Implements independent fading and shadowing models for each IRS within the swarm.
- Simulates various transmit SNR scenarios to evaluate system performance.
- Compares the swarm-based RIS system against the conventional single RIS setup.

## Simulation parameters
- **Number of Simulation Trials**: `1e5`
- **Number of Reflecting Elements (N)**: Configurable (`5` by default)
- **Transmit SNR**: Varies from `0 dB` to `20 dB`
- **Spectral Efficiency Threshold (R_th)**: Ranges from `1` to `4` b/s/Hz (comparison script) and fixed at `1` b/s/Hz (swarm script)
- **Path-Loss Exponent (PLE)**: `2.7`

## Methodology
1. **Nakagami-m Fading**: Models small-scale fading in wireless channels using shape (`m`) and spread (`Ω`) parameters.
2. **Inverse Gamma Shadowing**: Models large-scale fading effects using shape (`α`) and spread (`β`) parameters.
3. **Moment Matching**: Used to approximate the statistical properties of the received signal.
4. **Laplace Transform Analysis**: Validates the analytical approximation of the received signal distribution.
5. **Optimal Phase-Shift Configuration**: Ensures constructive interference at the receiver.
6. **Generalized Gaussian Quadrature**: Approximates the probability density function (PDF) of the aggregated received signal.

## Dependencies
- MATLAB (Tested on R2023a, but should be compatible with earlier versions)
- Statistics and Machine Learning Toolbox

## References

1. T. N. Do, G. Kaddoum, T. L. Nguyen, D. B. da Costa, and Z. J. Haas, **"Aerial Reconfigurable Intelligent Surface-Aided Wireless Communication Systems,"** *2021 IEEE 32nd Annual International Symposium on Personal, Indoor and Mobile Radio Communications (PIMRC)*, Helsinki, Finland, 2021, pp. 525-530. [DOI: 10.1109/PIMRC50174.2021.9569450](https://doi.org/10.1109/PIMRC50174.2021.9569450)  

2. B. Shang, E. Bentley, and L. Liu, **"UAV Swarm-Enabled Aerial Reconfigurable Intelligent Surface: Modeling, Analysis, and Optimization,"** *IEEE Transactions on Communications*, vol. PP, pp. 1-1, 2022. [DOI: 10.1109/TCOMM.2022.3173369](https://doi.org/10.1109/TCOMM.2022.3173369)
