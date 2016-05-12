This is a Matlab code for convolutional sparse coding, implementing the method proposed in 

*Michal Sorel, Filip Sroubek, "Fast convolutional sparse coding using matrix inversion lemma, Digital Signal Processing", 2016*
<http://www.sciencedirect.com/science/article/pii/S1051200416300276>

Convolutional sparse coding is an alternative to standard sparse coding better suited for modelling shift-invariant signals.
The most time-consuming part of both kernel learning and feature extraction is inversion of certain linear operator 
related to convolution. In this work we show how these inversions can be computed non-iteratively in 
the Fourier domain using the matrix inversion lemma even for multiple training signals. 
This greatly speeds up computation and  makes convolutional sparse coding computationally feasible even for large problems.

The matrix inversion lemma to speed up the convolutional sparse coding was already used in 
recent papers 
*B. Wohlberg, "Efficient Convolutional Sparse Coding", 2014*, *F. Heide, W. Heidrich, G. Wetzstein, "Fast and flexible convolutional sparse coding", 2015*
and *B. Wohlberg, "Efficient Algorithms for Convolutional Sparse Representations", 2016*.

**The main problem of these methods is that they work efficiently only with one input image.
In practice we would like the kernels to represent well not one particular signal/image but a large
database of examples. Our main contribution is that the proposed method is fast even 
for large training sets.**

For a more detailed description see <http://zoi.utia.cas.cz/convsparsecoding>.

**Files**

The code consists of the main function *convsparseF* implementing several variants of
the proposed algorithm and our implementation of other methods 
*M.D. Zeiler, D. Krishnan, G.W. Taylor, R. Fergus, "Deconvolutional networks", 2010* (method 3),
*H. Bristow, A. Eriksson, S. Lucey, "Fast convolutional sparse coding", 2013* (methods 1 and 2) and
*B. Wohlberg, "Efficient Algorithms for Convolutional Sparse Representations", 2016* (method 0).

*convsparseF.m* - main function (see help included in the file)

*out_mesto_4.mat* - data used as input images for learning kernels (images of Prague
Old Town and Prague Castle). They are stored as a 4D array.

*exp_feature_learning_E2time.m* comparison of convergence for the proposed algorithm and two 
	older algorithms

*exp_feature_learning_convergence_f.m*  kernel learning experiment for multi-scale kernels. 
To get really beautiful kernels, algorithm requires many iterations (1000-2000). It also consumes much
memory so be sure you have a machine with at least 16GB of memory and preferably more.

**Terms of Use**

This code can be freely used for research purposes by academic and other non-profit organizations. 
If you find it useful, please cite our paper *Michal Sorel, Filip Sroubek, "Fast convolutional
sparse coding using matrix inversion lemma", Digital Signal Processing, 2016*.