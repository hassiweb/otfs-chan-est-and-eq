# Channel Estimation and Equalization for CP-OFDM-based OTFS in Fractional Doppler Channels

This is a code package of the following scientific paper:

- N. Hashimoto, N. Osawa, K. Yamazaki and S. Ibi, "[*Channel Estimation and Equalization for CP-OFDM-based OTFS in Fractional Doppler Channels*](https://ieeexplore.ieee.org/abstract/document/9473532)," 2021 IEEE International Conference on Communications Workshops (ICC Workshops), 2021, pp. 1-7, doi: 10.1109/ICCWorkshops50388.2021.9473532.
- Noriyuki HASHIMOTO, et al., "[*Channel Estimation and Equalization for CP-OFDM-based OTFS in Fractional Doppler Channels*](https://arxiv.org/abs/2010.15396)," arXiv:2010.15396v3 [cs.IT], Jan. 2021.

## Abstract
Orthogonal time frequency and space (OTFS) modulation is a promising technology that satisfies high Doppler requirements for future mobile systems. OTFS modulation encodes information symbols and pilot symbols into the two-dimensional (2D) delay-Doppler (DD) domain. The received symbols suffer from inter-Doppler interference (IDI) in the fading channels with fractional Doppler shifts that are sampled at noninteger indices in the DD domain. IDI has been treated as an unavoidable effect because the fractional Doppler shifts cannot be obtained directly from the received pilot symbols. In this paper, we provide a solution to channel estimation for fractional Doppler channels. The proposed estimation provides new insight into the OTFS input-output relation in the DD domain as a 2D circular convolution with a small approximation. According to the input-output relation, we also provide a low-complexity channel equalization method using the estimated channel information. We demonstrate the error performance of the proposed channel estimation and equalization in several channels by simulations. The simulation results show that in high-mobility environments, the total system utilizing the proposed methods outperforms orthogonal frequency division multiplexing (OFDM) with ideal channel estimation and a conventional channel estimation method using a pseudo sequence. 

## Content of Code Package
Main funcitons of this code package are `OTFS.m` and `OFDM.m`.  Fig. 3 in the paper is generated using these codes.

These codes are frameworks for the transceiver of OTFS and OFDM, respectively.  Most of algorithms are definded in class files located in `classes` directory.

The proposed channel estimation algorithm is defined in [`classes/OtfsPilotResponseBasedPathParameterEstimator.m`](classes/OtfsPilotResponseBasedPathParameterEstimator.m), and the proposed channel equalization algorithm is defined in [`classes/OtfsDeconvolutionalEqualizer.m`](classes/OtfsDeconvolutionalEqualizer.m).

## License and Referencing
This code package is licensed under the MIT license.  If you in any way use this code for research that results in publications, please cite our paper refered above.
