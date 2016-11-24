This is a Matlab code for convolutional sparse coding, implementing the method proposed in 

*Michal Sorel, Filip Sroubek, "Fast convolutional sparse coding using matrix inversion lemma", Digital Signal Processing,
 vol. 55, pp. 44-51, 2016* (<http://www.sciencedirect.com/science/article/pii/S1051200416300276>)

Convolutional sparse coding is an alternative to standard sparse coding better suited for modelling shift-invariant signals.
The most time-consuming part of both kernel learning and feature extraction is inversion of certain linear operator 
related to convolution. In this work we show how these inversions can be computed non-iteratively in 
the Fourier domain using the matrix inversion lemma even for multiple training signals. 
This greatly speeds up computation and  makes convolutional sparse coding computationally feasible even for large problems.

The matrix inversion lemma to speed up the convolutional sparse coding was already independently used in 
recent papers of B. Wohlberg and Heide et al. (see bellow) 
**The main problem of these methods is that they work efficiently only with one input image.
In practice we would like the kernels to represent well not one particular signal/image but a 
set of examples. The main advantage of the proposed method is that it is fast even 
in such cases.**

For a more detailed description see <http://zoi.utia.cas.cz/convsparsecoding> or the paper itself.

**Files**

*example_trainkernels.m* - example of training kernels

*example_featureextraction.m* - example of computing feature maps, once kernels are trained

*convsparseF.m* - main function implementing several variants of the proposed algorithm
and our implementation of other methods (for details, see help included in the file). 

+ 3D variant of the proposed algorithm (method = -1)
+ Consensus variant of the proposed algorithm (method = -3)
+ Tiling variant of the proposed algorithm (method = -4)
+ *M.D. Zeiler, D. Krishnan, G.W. Taylor, R. Fergus, "Deconvolutional networks", 2010* (method =  3)
+ *H. Bristow, A. Eriksson, S. Lucey, "Fast convolutional sparse coding", 2013* (method = 1 and method = 2) 
+ *B. Wohlberg, "Efficient Algorithms for Convolutional Sparse Representations", 2016* (method = 0)

*out_mesto_4.mat* - data used as input images for learning kernels (images of Prague
Old Town and Prague Castle). They are stored as a 4D array.

*exp_feature_learning_E2time.m* comparison of convergence for the proposed algorithm and two 
	older algorithms

*exp_feature_learning_convergence_f.m*  kernel learning experiment for multi-scale kernels. 

**Terms of Use**

This code can be freely used for research purposes by academic and other non-profit organizations. 
If you use it in your research, please cite our paper *Michal Sorel, Filip Sroubek, "Fast convolutional
sparse coding using matrix inversion lemma", Digital Signal Processing, vol. 55, pp. 44-51, 2016*.