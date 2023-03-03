# dispersion

## Overview

The real and imaginary part of certain functions, called analytic functions, are related through the Hilbert transform. The Kramers-Kronig relation is another name for the Hilbert transform. The anomalous structure factor corrections, f' and f", are related through the Hilbert transform/Kramers-Kronig relation. This code implements a penalty function that takes f' and f" as input, and then computes the corresponding f" and f', respectively. We can calculate f" from f' with the Hilbert transform, and we can calculate f' from f" by taking the Hilbert transform and negating the result. The Hilbert transform is implemented with the Fast Fourier Transform (FFT), directly copying the method in [scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html), except re-implementing in PyTorch. The mean squared error comparing the computed f' and f" with the input f' and f", respectively, is calculated, and can be used as a penalty term in an optimization problem. The code supports automatic differentiation through implementation in PyTorch, so gradients of the penalty term can easily be obtained. The main use case for this code is to optimize the anomalous structure factor corrections for serial crystallography data to determine the electronic state of atoms.

One problem in implementing the Hilbert transform is that the values of the function at all points needs to be known. We get around this by realizing that the Hilbert transform is linear, i.e. the Hilbert transform of a sum of two functions is the sum of the Hilbert transform of each individual function. Generally, we have measurements over a large measurement range of energies for atoms in their ground state. We realize that the change in the anomalous structure factor for an atom from the ground state to a different electronic state is very small, and the change is only around the band-edge. That allows us to optimize *only* the change of the anomalous structure factor from the ground state, a function that is zero except around the band-edge. This problem is much more computationally tractable as we only have to optimize the values of f' and f" in a small window around the band-edge. As this function goes to zero at both ends, this method also prevents numerical edge effects due to the implementation of the Hilbert transform with the FFT.

## Requirements

Install PyTorch with the following command:

`conda install pytorch==1.13.1 torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia`

## References

[Scipy Documentation of the Hilbert Transform](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.hilbert.html)

[Wikipedia](https://en.wikipedia.org/wiki/Hilbert_transform)
