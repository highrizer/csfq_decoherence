# csfq_decoherence
This repository contains data on the decoherence rate of a tunable capacitively-shunted flux qubit (CSFQ), and code used to analyze them. The corresponding manuscript can be found on [arXiv:2307.13961](https://doi.org/10.48550/arXiv.2307.13961).

The respository is organized as follows
- Data: contains experimental data
- Simulation: contains simulated decoherence rates
- DeviceParameters: contains circuit model parameters and qubit sensitivities to various noise channels, as simulated using a generic circuit Hamiltonian builder (private code, methodoligies explained in detail in [arXiv:2010.14929](http://arxiv.org/abs/2010.14929)
- TunableResonator.nb: notebook for finding the symbolic Hamiltonian of a rf-SQUID terminated CPW resoantor, and numerical computation of its interaction strength with the CSDQ
- TunableResonator.wl: functions for numerically finding the interaction rates of the tunable resonator with the CSFQ
- decoherence.py functions for evaluating relaxation and dephasing from various noise channels
- awg_dephasing.ipynb: notebook for evaluating the impact of AWG noise on qubit dephasing
- junction_noise_dephasing.ipynb: notebook for evaluating the impact of junction critical current noise on qubit dephasing
- fit_dephasing_and_relaxation.ipynb: notebook for fitting experimental data to extract noise amplitudes

	