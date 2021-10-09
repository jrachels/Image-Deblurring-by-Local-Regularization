# Image-Deblurring-by-Local-Regularization

The code implements the local regularization method for deblurring in 2 dimensions. The paper it's based on can be found here:

https://users.math.msu.edu/users/lamm/Preprints/Cui_Lamm_Scofield_images/paper.pdf .

The contents of this repository were created by a team of students at the 2016 Mathematics REU at Willamette. I haven't touched any of this in years. 

To see it in action, just run alpha_test.m. With the default parameters, it will take about an hour to run. The output should be the three images shown below. 

The first image is the "real image" that the user/researcher wants to construct. 

![original](https://user-images.githubusercontent.com/92210470/136655995-93917dfe-00fc-4796-9081-d404843bc163.jpg)

The second image is a Gaussian blur of the real image with some added Poison noise. This is the input the local regularization algorithm receives.

![noisyblurred](https://user-images.githubusercontent.com/92210470/136655988-28f07dce-4822-4d22-a2fd-980c355f3108.jpg)


The third image is the attempted reconstruction of the original image. This is the output of the local regularization algorithm. Our measure of success is how closely this image resembles the first image. In particular, we are looking at the "sharpness" of the reconstructed image.

![reconstruction](https://user-images.githubusercontent.com/92210470/136656004-6d462a58-0a3f-455d-8313-730ebb83436f.jpg)

## References
<a id="1">[1]</a> 
C. Cui, P.K. Lamm and T.L. Scofield, *Local regularization for n-dimensional integral equations  with applications to image processing*, Inverse Problems **23** (2007), 1611-1633.
