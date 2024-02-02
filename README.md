[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=emilbjornson/RIS-massive-MIMO)

Is Channel Estimation Necessary to Select Phase-Shifts for RIS-Assisted Massive MIMO?
==================

This is a code package is related to the following scientific article:

Özlem Tuğfe Demir and Emil Björnson, “[Is Channel Estimation Necessary to Select Phase-Shifts for RIS-Assisted Massive MIMO?](https://ieeexplore.ieee.org/document/9786573
),” IEEE Transactions on Wireless Communications, vol. 21, no. 11, pp. 9537-9552, November 2022, doi: 10.1109/TWC.2022.3177700.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. We encourage you to also perform reproducible research!


## Abstract of Article

Reconfigurable intelligent surfaces (RISs) consist of many passive elements of metamaterials whose impedance can be controllable to change the characteristics of wireless signals impinging on them. Channel estimation is a critical task when it comes to the control of a large RIS when having a channel with a large number of multipath components. In this paper, we derive Bayesian channel estimators for two RIS-assisted massive multiple-input multiple-output (MIMO) configurations: i) the short-term RIS configuration based on the instantaneous channel estimates; ii) the long-term RIS configuration based on the channel statistics. The proposed methods exploit spatial correlation characteristics at both the base station and the planar RISs, and other statistical characteristics of multi-specular fading in a mobile environment. Moreover, a novel heuristic for phase-shift selection at the RISs is developed. A computationally efficient fixed-point algorithm, which solves the max-min fairness power control optimally, is proposed. Simulation results demonstrate that the proposed uplink RIS-aided framework improves the spectral efficiency of the cell-edge mobile user equipments substantially in comparison to a conventional single-cell massive MIMO system. The impact of several channel effects are studied to gain insight about when the channel estimation, i.e., the short-term configuration, is preferable in comparison to the long-term RIS configuration to boost the spectral efficiency.

## Content of Code Package

The article contains 8 simulation figures, numbered 5-12. Figures 5, 6, and 7 are generated by the Matlab script Figs5_6_7.m. Figures 8, 9, 10, 11, and 12 are generated by the Matlab scripts Fig8.m, Fig9.m, Fig10.m, Fig11.m, and Fig12.m, respectively. The package also contains other Matlab functions used by some of the scripts.



See each file for further documentation.

## Acknowledgements

The work of Ö. T. Demir and E. Björnson were supported by the Swedish Foundation for Strategic Research under Grant FFL18-0277. 

## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
