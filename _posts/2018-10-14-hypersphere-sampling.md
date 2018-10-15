---
layout: post
title: Sampling in the hypersphere
author: Adrian Baez-Ortega
date: 2018/10/14
---


This post is the first of a short series where I will be presenting some of the most remarkable and surprising concepts and techniques from the field of statistical physics, all of which are beautifully featured in Werner Krauth's marvelous book [*Statistical Mechanics: Algorithms and Computations*](https://global.oup.com/academic/product/statistical-mechanics-algorithms-and-computations-9780198515364) (Oxford University Press). (As a side note, I positively recommend Krauth's related *free* [online course](https://www.coursera.org/learn/statistical-mechanics), which is every bit as delightful as his book, and more accessible to those lacking a physics background.)

The problem in focus today is deceptively simple: how to efficiently sample random points inside a sphere. Such sphere, however, may be defined not in the familiar three dimensions, but in any number of dimensions. We shall see that an intuitive but naive sampling strategy, even though at first provides a solution, is notwithstanding bound to become highly inefficient as the number of dimensions, *D*, increases; and, in our way towards a superior sampling algorithm, we shall encounter some amazingly clever and powerful mathematical concepts, including sample transformations and isotropic probability distributions.

Let us begin by defining the objects in question. We use the term *hypersphere* to refer to the general version of the sphere in an arbitrary number of dimensions. Thus, the two-dimensional hypersphere and the three-dimensional hypersphere correspond, respectively, to what we normally call a circle and a sphere. Hyperspheres with many more dimensions can also be defined, even if our capacity for picturing them is impaired by our everyday experience of three-dimensional space (try for a moment to imagine the aspect of a sphere in four dimensions!). In order to save us from having to read many long words, however, we shall refer to hyperspheres simply as spheres. In addition, we only consider, for simplicity, the case of the *unit sphere*, which has a radius *r* = 1.

Analogously, the *hypercube* is the general version of the cube in any dimensions, and we will simply call it the *D*-dimensional cube. We should note that, for *D* = 1, both the one-dimensional sphere and one-dimensional cube correspond to a straight line segment.

The reason why cubes are of interest here, is that they are instrumental in the simplest and most intuitive way of sampling random points inside a sphere. Let us consider, for example, the sphere in *two dimensions* — that is, the circle. The most obvious technique for sampling points inside a circle is to first sample random points (pairs of *x* and *y* coordinates) inside a *square* containing the circle, and then *accepting* each sampled point only if it is indeed inside the circle (that is, if its distance to the centre is smaller than the circle radius).

We can easily write a short algorithm for this naive sampling method, and plot the resulting points. We will plot any points lying outside the circle (rejected samples) in red, and any points inside the circle (accepted samples) in blue.

``` r
# Naive sampling of N random points (each a vector of coordinates 
# [x_1, ..., x_D]) inside a D-dimensional sphere of radius r
naive.hypersphere = function(N, D, r=1) {
    set.seed(1)
    samples = list(accepted = matrix(NA, nrow=N, ncol=D, 
                                     dimnames=list(NULL, paste0("x_", 1:D))),
                   rejected = matrix(NA, nrow=1e6, ncol=D,
                                     dimnames=list(NULL, paste0("x_", 1:D))))
    i = j = 0
    while (i < N) {
        
        # For each dimension (1, 2, ..., D),
        # sample a random coordinate in [-1, 1]
        point = runif(n=D, min=-1, max=1)
        
        # Check if the point is inside the sphere,
        # i.e., if (x_1^2 + x_2^2 + ... + x_D^2) < r
        accept = sum(point ^ 2) < r
        
        # Add to the table of accepted/rejected samples
        if (accept) {
            i = i + 1
            samples$accepted[i, ] = point
        }
        else {
            j = j + 1
            samples$rejected[j, ] = point
        }
        
    }
    samples$rejected = samples$rejected[1:j, ]
    samples
}
```

``` r
# Sample 1000 random points inside the circle of radius 1
points.circle = naive.hypersphere(N=5000, D=2)

# Plot accepted points in blue and rejected points in red
plot(points.circle$rejected, pch=16, col="firebrick", 
     xlab="x", ylab="y", main="5000 random points")
points(points.circle$accepted, pch=16, col="cornflowerblue")
```

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/unnamed-chunk-2-1.png)

This method can be applied to the sphere in any number of dimensions, as illustrated below for the cases *D* = 1, 2, 3.

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/sampling.png)

As you might have already noticed, the problem with this sampling algorithm is that the region in which the points are sampled (the cube) and the region where they are accepted (the sphere) rapidly grow more distinct from each other as we go to higher dimensions. While for *D* = 1 the sampling is perfect (as both sphere and cube correspond to the same line segment), we already notice that the difference between the 3-dimensional cube and the sphere within is substantially larger than the difference between the 2-dimensional square and its circle. In general, the ratio of the volume of the sphere to that of the cube gives the algorithm's *acceptance rate* (the fraction of accepted samples) for *D* dimensions, as shown below.

<table style="width:100%;">
<colgroup>
<col width="11%" />
<col width="15%" />
<col width="13%" />
<col width="20%" />
<col width="23%" />
<col width="15%" />
</colgroup>
<thead>
<tr class="header">
<strong>
<th align="center">Dimensions</th>
<th>Sphere analogue</th>
<th>Cube analogue</th>
<th>Volume of unit sphere</th>
<th>Volume of cube of side 2</th>
<th>Acceptance rate</th>
</strong>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center"><strong><em>D</em> = 1</strong></td>
<td>Line segment</td>
<td>Line segment</td>
<td><em>V</em><sub>S</sub> = length = 2</td>
<td><em>V</em><sub>C</sub> = length = 2</td>
<td><em>V</em><sub>S</sub> / <em>V</em><sub>C </sub> = 1</td>
</tr>
<tr class="even">
<td align="center"><strong><em>D</em> = 2</strong></td>
<td>Circle</td>
<td>Square</td>
<td><em>V</em><sub>S</sub> = area = 𝜋</td>
<td><em>V</em><sub>C</sub> = area = 4</td>
<td><em>V</em><sub>S</sub> / <em>V</em><sub>C </sub> = 𝜋 / 4 ≈ 0.79</td>
</tr>
<tr class="odd">
<td align="center"><strong><em>D</em> = 3</strong></td>
<td>Sphere</td>
<td>Cube</td>
<td><em>V</em><sub>S</sub> = volume = 4/3 𝜋</td>
<td><em>V</em><sub>C</sub> = volume = 8</td>
<td><em>V</em><sub>S</sub> / <em>V</em><sub>C </sub> = 𝜋 / 6 ≈ 0.52</td>
</tr>
<tr class="even">
<td align="center">...</td>
<td>...</td>
<td>...</td>
<td>...</td>
<td>...</td>
<td>...</td>
</tr>
<tr class="odd">
<td align="center"><strong>Large <em>D</em></strong></td>
<td>Hypersphere</td>
<td>Hypercube</td>
<td><em>V</em><sub>S</sub> = 𝜋<sup> (<i>D</i> / 2)</sup> <em>r<sup> D</sup></em> / Γ((<em>D</em> / 2) + 1)</td>
<td><em>V</em><sub>C</sub> = 2<em><sup> D</sup></em></td>
<td><em>V</em><sub>S</sub> / <em>V</em><sub>C </sub> ≪ 1</td>
</tr>
</tbody>
</table>

The acceptance rate already drops dramatically for *D* = 3, foreshadowing the method's fate in higher dimensions. For large values of *D*, this method is definitely too inefficient to be of any use, as nearly all of the proposed samples are rejected.

As an interesting side note, notice that the number 𝜋 is involved in the acceptance rate; this implies that, as we sample and accept/reject points with this algorithm, we are inadvertently computing the value of 𝜋 to an increasing degree of accuracy. For example, we can approximate 𝜋 from the empirical acceptance rate of our previous sampling exercise on the circle.

``` r
# For D = 2:
# Acceptance rate = (accepted samples) / (total samples) = 𝜋 / 4
accept.rate = nrow(points.circle$accepted) / 
    (nrow(points.circle$accepted) + nrow(points.circle$rejected))
accept.rate
```

    ## [1] 0.7818608

``` r
# Approximation of 𝜋
accept.rate * 4
```

    ## [1] 3.127443

It seems that we should sample quite a bit longer if we are to obtain an acceptable estimate of 𝜋!

We have seen that the acceptance-rate problem derives from the fact that we are sampling from a distribution which does not exactly match the region we are interested in. Let us now turn to a slightly more sophisticated approach, where the coordinates are sampled not from a uniform distribution between –1 and 1, but from a standard Gaussian (or normal) distribution. We will see that, by using a simple transformation, a vector of *D* Gaussian random numbers can be converted into a random point inside the *D*-dimensional unit hypersphere, or even on its surface. This method achieves perfect sampling, meaning that no samples have to be rejected, and thus the acceptance rate is 1. Hence, this algorithm will not suffer from the efficiency concerns of the previous one.

Let us begin by sampling points as vectors of *D* independent Gaussian coordinates. In two dimensions, we sample two coordinates, *x* and *y*, which gives the following distribution of samples.

``` r
# Plot 1000 random points in two dimensions,
# each composed of two independent Gaussian random coordinates
plot(x = rnorm(n=1000),
     y = rnorm(n=1000),
     pch=16, col="cornflowerblue", xlab="x", ylab="y", main="1000 random points")
```

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/unnamed-chunk-4-1.png)

So far, this looks more like a cloud than a circle. The crucial detail here, however, is that the points in this cloud are uniformly distributed in *all directions*. This is because the combination of two Gaussian distributions gives an *isotropic* distribution in two dimensions. [Isotropic](https://en.wikipedia.org/wiki/Isotropy) is just a fancy Greek name for something that is uniform in all orientations; in other words, something whose aspect does not change when it is rotated. An orange is normally isotropic, while a banana never is.

The Gaussian distribution has the *unique* property that the combination of *D* independent distributions generates an isotropic distribution in *D* dimensions. *This is the key insight* that enables perfect sampling in a sphere of any dimensionality.

(The proof that the Gaussian distribution is isotropic is itself deeply elegant, but will not be shown here. Briefly, it involves some clever variable transformations on an integral over the Gaussian distribution, resulting in two independent integrals on the angular and radial variables. It then becomes evident that the integral over the angular variable spans all angles uniformly, which is the definition of isotropy.)

Let us now climb from two to three dimensions, where the combination of three independent Gaussian variables (for *x*, *y* and *z*) again produces an isotropic, globular cluster of points.

``` r
# Plot 1000 random points in three dimensions (using package plot3D),
# each composed of three independent Gaussian random coordinates
plot3D::points3D(x = rnorm(n=1000),
                 y = rnorm(n=1000),
                 z = rnorm(n=1000),
                 col="cornflowerblue", pch=16, alpha=0.7,
                 bty="g", phi=20, main="1000 random points")
```

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/unnamed-chunk-5-1.png)

Now, all we have to do is reshaping this cloud into a perfect sphere. This requires transforming the distribution of distances, or *radii*, of the points, so that they are uniformly distributed across all radii. Because the volume increases exponentially with distance from the centre, however, longer radii demand more points than shorter ones (imagine that each point were a small sphere in itself, and consider how many would be needed to fill the space just beneath the sphere's surface, in contrast with the amount needed to fill the space surrounding the sphere centre).

To be precise, achieving uniform sampling in terms of radii requires a probability distribution proportional to the (*D* – 1)-th power of the radius *r* (with 0 &lt; *r* &lt; 1):

<center>
𝜋(<i>r</i>) ∝ <i>r</i><sup><i> D</i> – 1</sup>
</center>

Fortunately, sampling from this distribution is easy: we only need to sample from a uniform distribution between 0 and 1, and take the *D*-th root of each sample; the transformed samples will then follow the exponential distribution shown above.

To assign the radii from the distribution above to the points that we sampled from the isotropic Gaussian distribution, we need to divide the value of each coordinate by the point's original distance (given by the square root of the sum of its squared coordinates), and then multiply the coordinates again by the new radius, which is sampled from the distribution above. After this correction, the points will uniformly cover not just all angles, but also all distances between 0 and 1.

Below is a function that ties all these ideas together into an algorithm for perfect sampling of points in a *D*-dimensional sphere.

``` r
# Perfect sampling of N random points (vectors of coordinates 
# [x_1, ..., x_D]) in/on a D-dimensional sphere of radius r
perfect.hypersphere = function(N, D, r=1, surface=FALSE) {
    
    # Sample D vectors of N Gaussian coordinates
    set.seed(1)
    samples = sapply(1:D, function(d) {
        rnorm(n=N)
    })
    
    # Normalise all distances (radii) to 1
    radii = apply(samples, 1, function(s) {
        sqrt(sum(s ^ 2))
    })
    samples = samples / radii
    
    # Sample N radii with exponential distribution
    # (unless points are to be on the surface)
    if (!surface) {
        new.radii = runif(n=N) ^ (1 / D)
        samples = samples * new.radii
    }
    
    samples
}
```

Note that this function includes an additional argument, `surface`, which allows us to skip the sampling of new radii after the points have been normalised to have radius equal to 1. This results in points which are uniformly distributed *on the surface* of the sphere, instead of inside it.

Below are two examples of 1000 random points, *inside* and *on* the unit sphere, respectively.

``` r
# Sample 1000 random points inside the unit sphere (D=3)
points.sphere = perfect.hypersphere(N=1000, D=3)
plot3D::points3D(x = points.sphere[, 1],
                 y = points.sphere[, 2],
                 z = points.sphere[, 3],
                 col="cornflowerblue", pch=16, alpha=0.7, bty="g", phi=20, 
                 main="1000 random points inside the unit sphere (D = 3)")
```

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
# Sample 1000 random points on the surface of the unit sphere (D=3)
points.sphere = perfect.hypersphere(N=1000, D=3, surface=TRUE)
plot3D::points3D(x = points.sphere[, 1],
                 y = points.sphere[, 2],
                 z = points.sphere[, 3],
                 col="cornflowerblue", pch=16, alpha=0.7, bty="g", phi=20, 
                 main="1000 random points on the unit sphere surface (D = 3)")
```

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/unnamed-chunk-7-2.png)

This perfect sampling algorithm, if formally quite simple, provides a powerful example of the importance of understanding the topology of the space in which sampling is to be done, and of the beauty and elegance which so often characterises the mathematical tools, such as variable and sample transformations, that are ubiquitous in statistical physics.

Finally, it is very important to realise that, although I have chosen to present this sampling approach mainly because of the mathematical ingenuity that it encapsulates, hypersphere sampling is not just simply beautiful, but also of prime importance in statistical mechanics and other branches of physics. To mention one example, the velocities of *N* classical particles in a gas with constant [kinetic energy](https://en.wikipedia.org/wiki/Kinetic_energy) are distributed as random points on the surface of a 3*N*-dimensional hypersphere. This is because the fixed kinetic energy is proportional to the sum of the squared velocities, and therefore any valid set of particle velocities (each with *x*, *y* and *z* components) must satisfy that the sum of their squares is proportional to the kinetic energy — just as any point on the surface of a sphere has coordinates whose sum-of-squares is, by definition, equal to the square of the sphere radius.

Interestingly, the distribution of points in a hypersphere is also connected with [J. J. Thomson](https://en.wikipedia.org/wiki/J._J._Thomson)'s landmark 1904 [atomic model](https://en.wikipedia.org/wiki/Plum_pudding_model), which described the atom as a positively charged sphere with negatively charged electrons embedded in it. More generally, the problem of carefully defining the sampling region is of foremost importance in [Markov chain Monte Carlo](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampling, being the motivation of some of the most advanced methods in the field, such as Hamiltonian Monte Carlo and Adiabatic Monte Carlo.

­

*NB. This post was largely inspired by the contents of Werner Krauth's book* [Statistical Mechanics: Algorithms and Computations](https://global.oup.com/academic/product/statistical-mechanics-algorithms-and-computations-9780198515364) *(Oxford University Press), and its companion [online course](https://www.coursera.org/learn/statistical-mechanics).*
