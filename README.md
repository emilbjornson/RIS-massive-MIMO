Is Channel Estimation Necessary to Select Phase-Shifts for RIS-Assisted Massive MIMO?
==================

This is a code package is related to the following scientific article:

Özlem Tuğfe Demir and Emil Björnson, “[Is Channel Estimation Necessary to Select Phase-Shifts for RIS-Assisted Massive MIMO?](https://ieeexplore.ieee.org/document/9786573
),” IEEE Transactions on Wireless Communications, vol. 21, no. 11, pp. 9537-9552, November 2022, doi: 10.1109/TWC.2022.3177700.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. We encourage you to also perform reproducible research!


## Abstract of Article

Reconfigurable intelligent surfaces (RISs) consist of many passive elements of metamaterials whose impedance can be controllable to change the characteristics of wireless signals impinging on them. Channel estimation is a critical task when it comes to the control of a large RIS when having a channel with a large number of multipath components. In this paper, we derive Bayesian channel estimators for two RIS-assisted massive multiple-input multiple-output (MIMO) configurations: i) the short-term RIS configuration based on the instantaneous channel estimates; ii) the long-term RIS configuration based on the channel statistics. The proposed methods exploit spatial correlation characteristics at both the base station and the planar RISs, and other statistical characteristics of multi-specular fading in a mobile environment. Moreover, a novel heuristic for phase-shift selection at the RISs is developed. A computationally efficient fixed-point algorithm, which solves the max-min fairness power control optimally, is proposed. Simulation results demonstrate that the proposed uplink RIS-aided framework improves the spectral efficiency of the cell-edge mobile user equipments substantially in comparison to a conventional single-cell massive MIMO system. The impact of several channel effects are studied to gain insight about when the channel estimation, i.e., the short-term configuration, is preferable in comparison to the long-term RIS configuration to boost the spectral efficiency.

## Content of Code Package

The article contains 12 simulation figures, numbered 2-13. Figures 2 and 3 are generated by the Matlab script Fig2_3.m. Figures 4, 5, 12, and 13 are generated by the Matlab script Fig4_5_12_13.m. Figures 6 and 7 are generated by the Matlab script Fig6_7.m. Figures 8, 9, 10, and 11 are generated respectively by the Matlab scripts Fig8.m, Fig9.m, Fig10.m, and Fig11.m. The package also contains the Matlab scripts functionExampleSetup.m, functionExpectations_lmmse.m, functionExpectations_lmmse_random_pilot_sequence.m, functionExpectations_ls.m, optimize_power_coh_lin.m, optimize_power_coh_nonlin.m, optimize_power_noncoh_lin.m, and optimize_power_noncoh_nonlin.m that are MATLAB functions used by some of the scripts.

For the optimization problems, the convex programming solver CVX http://cvxr.com/cvx/ is used. Please adjust the default solver of CVX to SDPT3 not to encounter any issues. The version we have tested was SDPT3 4.0. 

See each file for further documentation.

## Acknowledgements

The work of Ö. T. Demir and E. Björnson was partially supported by ELLIIT and the Wallenberg AI, Autonomous Systems and Software Program (WASP) funded by the Knut and Alice Wallenberg Foundation. 

## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
