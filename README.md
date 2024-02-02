# INVDFT: Java Program to optimize parameter for Exchange Correlation Functional in  Density Functional Theory
INVDFT, a Java based software, employs the SCF technique to solve the KSDFT equation. The optimization of XC function parameters in INVDFT is accomplished using the Newton-Raphson and Monte Carlo methods
## Features
- Parameter optimization of Exchange Correlation Functional
- SCF with electron density mixing
- DFT: Neural Network, LDA and GGA functionals
- Hartree-Fock
## Requirements
- Linux OS
- Java (https://www.java.com/en/download/) version >= 11
-  MySQL for training ANN Exchange
## Example of parameter optimization

The LSDA model in equation  can be simplified to equation. We will demonstrate how to use INVDFT to find the parameters $\beta, \gamma$ in equation and then compare it with the LSDA model, where the LSDA model has $\beta=-0.930525$ and $\gamma=1.3333$. 
$\epsilon_{xc}=\beta \int(\rho_{\-}^\gamma+\rho_+^\gamma) dr^3$

------------
	{
	"name": ["H","OH,O,H}","H2,H,H}","Na2,Na,Na}","H2S,S,H,H}"],
	"con": 0.0015936254980079682,
	"re": [
	[-0.5],[-0.162231075697211],[-0.164621513944223],[-0.0264541832669323],[-0.27601593625498]
                          ],
	"simpul": [1,1],
	"name_ann": "ann_new_dft"
	}
------------
## Example of DFT Calculation
------------
	H2O;
	Water;
	{"xyz": 
	{"1": [0.0000, 0.0000, 0.1173], 
	"2": [0.0000, 0.7572, -0.4692],
	"3": [0.0000, -0.7572, -0.4692]}, 
	"atom": ["O", "H","H"],
	"Spin_dn": 5, "Spin_up": 5};
------------
## Example of Hartree Fock Calculation
