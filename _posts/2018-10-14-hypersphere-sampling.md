---
layout: post
title: Sampling in the hypersphere
author: Adrian Baez-Ortega
date: 2018/10/14
---


This is the first of a short series of posts where I will present some of the most beautiful and surprising concepts from the field of statistical physics, all of which are featured in Werner Krauth's marvelous book [*Statistical Mechanics: Algorithms and Computations*](https://global.oup.com/academic/product/statistical-mechanics-algorithms-and-computations-9780198515364) (Oxford University Press). (Incidentally, I recommend taking a look at the related *free* [online course](https://www.coursera.org/learn/statistical-mechanics) available on Coursera, which is every bit as delightful as Krauth's book, and more accessible to those without a physics background.)

The problem in focus today is deceptively simple: how to efficiently sample random points inside a sphere in *D* dimensions. We will see that the more intuitive, naive sampling strategy is bound to become unbearably inefficient as the number of dimensions, *D*, increases; and, in our way to a much better sampling algorithm, we shall be witness to some amazingly clever and powerful mathematical concepts, including sample transformations and isotropic compound probability distributions.

Let us begin by defining the objects in question. We use the term *hypersphere* to refer to the general version of the sphere in an arbitrary number of dimensions. Thus, the two-dimensional hypersphere and the three-dimensional hypersphere correspond, respectively, to what we normally call *circle* and *sphere*. Hyperspheres can also be defined for much higher dimensions, even if our capacity to picture them is impaired by our everyday experience of space in three dimensions (try for a moment to imagine the shape of a sphere in four dimensions!). In order to save us from reading long words, however, we will simply refer to hyperspheres as D*-dimensional spheres*. Also, we will only consider, for simplicity, the case of the sphere whose radius is equal to one, also known as the *unit sphere*.

Analogously, the *hypercube* is the general version of a cube in *D* dimensions, and we will call it simply the *D*-dimensional cube. We should note that, for *D* = 1, both the one-dimensional sphere and the one-dimensional cube (and everything one-dimensional, really), correspond to a straight line (or segment).

The reason why cubes are of interest here is that they are required for the simplest and most intuitive way of sampling random points inside a sphere. Let us consider, for example, the sphere in two dimensions — the circle. The most obvious technique for sampling points inside a circle is to first sample random points (pairs of *x* and *y* coordinates) inside the square that contains the circle, and then checking if the points are inside the circle (that is, if their distance to the centre is smaller than 1), in which case we accept them as samples.

We can easily write a small algorithm for such naive sampling and plot the result. (Points that were sampled and rejected are plotted red, and accepted samples are shown in blue.)

``` r
# Naive sampling of N random points (vectors of coordinates 
# [x_1,...,x_D]) inside a D-dimensional sphere of radius r
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
        
        # Check if the sample point is within the sphere:
        # (x_1^2 + x_2^2 + ... + x_D^2) < r ?
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

As you have probably noticed already, the problem with this naive sampling method is that the volume where we sample the points (the cube) and the volume where we accept the points (the sphere) rapidly grow more different from each other as we go to higher dimensions. While for *D* = 1 the sampling is perfect (as both sphere and cube correspond to the same line segment), we already notice that the difference between the 3-dimensional cube and the sphere within is much larger than the difference between the circle and the square for *D* = 2. In general, the ratio of the volume of the sphere to that of the cube gives the *acceptance rate* (the fraction of accepted points) of the algorithm for *D* dimensions, as shown below.

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

The acceptance rate already drops dramatically for *D* = 3, foreshadowing the inefficiency of the method for higher dimensions. For large values of *D*, this method is indeed to slow to be of any use, as nearly all of the proposed samples are rejected.

One interesting detail is that the acceptance rate includes the number 𝜋; this means that, as we sample random points with this algorithm, we are unadvertently computing the value of 𝜋 to an increasing degree of accuracy. For example, we can approximate 𝜋 from the empirical acceptance rate of our previous sampling exercise on the circle.

``` r
# For D=2: acceptance rate = (accepted samples) / (total samples) = pi / 4
accept.rate = nrow(points.circle$accepted) / 
    (nrow(points.circle$accepted) + nrow(points.circle$rejected))
accept.rate
```

    ## [1] 0.7818608

``` r
# Approximation of pi
accept.rate * 4
```

    ## [1] 3.127443

It seems that we should have sampled longer to get a decent estimate of 𝜋!

We have seen that the problem with the acceptance rate derives from the fact that we are sampling from an unsuitable distribution. Let us now turn to a more clever approach, where coordinates are sampled not from a uniform distribution between –1 and 1, but from a standard Gaussian (or normal) distribution. We will see that, by using a simple transformation, a vector of *D* Gaussian numbers can be converted into a random point sampled from the *D*-dimensional unit hypersphere, or on its surface. This method achieves perfect sampling, meaning that there are no samples are rejected and the acceptance rate is 1.

Let us begin by sampling points as a vectors of *D* independent Gaussian coordinates. In two dimensions, we sample two coordinates, *x* and *y*, which gives the following distribution of samples.

``` r
# Plot 1000 random points in two dimensions,
# each composed of two independent Gaussian random coordinates
plot(x = rnorm(n=1000),
     y = rnorm(n=1000),
     pch=16, col="cornflowerblue", xlab="x", ylab="y", main="1000 random points")
```

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/unnamed-chunk-4-1.png)

So far, this looks more like a cloud than a circle. The crucial detail to notice, however, is that the points in this cloud are uniformly distributed in all directions. This is because the combination of two Gaussian distributions gives an *isotropic* distribution in two dimensions. [Isotropic](https://en.wikipedia.org/wiki/Isotropy) is just a fancy Greek term for something that is uniform in all orientations; in other words, something whose aspect does not change when it is rotated. Clearly, the shape of the cloud of points above is not going to be different if we rotate the plot to the left or to the right.

The Gaussian distribution is *unique*, in that the combination of *D* independent distributions produce an isotropic distribution in *D* dimensions. This is the key insight that enables perfect sampling in a sphere of any number of dimensions.

(The proof that the Gaussian distribution is isotropic in *D* dimensions is itself very beautiful, but will not be shown here. Briefly, it involves clever variable transformations on the integral that goes over the Gaussian distribution, in order to obtain two independent integrals on the angular and radial variables. It is then evident that the integral over the angular variable spans all angles uniformly, which is the definition of isotropy.)

Let us now move to three-dimensional space, where the combination of three independent Gaussian distributions (for *x*, *y* and *z*) again produces an isotropic cloud of points.

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

Now all we have left to do is to give this cloud a 'haircut' to fit within a three-dimensional sphere. This requires transforming the distribution of distances (radii) of the points, so that all radii are covered uniformly. Because volume increases exponentially with the radius, greater distances from the centre of the sphere require larger numbers of points than shorter distances. (Imagine that each point were a small sphere in itself, and think of how many would be needed to fill the space just beneath the surface of the sphere, in contrast with the number needed to fill the space surrounding the sphere centre.) To be precise, to achieve uniform sampling in terms of the radii requires a distribution of radii proportional to the (*D* – 1)-th power of each radius *r* (with 0 &lt; *r* &lt; 1):

<center>
𝜋(<i>r</i>) ∝ <i>r</i><sup><i>D</i> – 1</sup>
</center>
Fortunately, sampling from this distribution is easy: we only need to sample from a uniform distribution and take the *D*-th root of each sample; the transformed samples will follow the exponential distribution above.

To assign the radii from the distribution above to the points that we have sampled from the isotropic Gaussian distribution, we need to divide the value of each coordinate by the point's current distance (the square root of the sum of its squared coordinates), and then multiply the coordinates by the new radius, sampled from the distribution above. Once this adjustment is done, the points will not only cover all angles uniformly, but also all distances between 0 and 1.

Below is a function that ties all these ideas together into a method for perfect sampling of points in a *D*-dimensional sphere.

``` r
# Perfect sampling of N random points (vectors of coordinates 
# [x_1,...,x_D]) in/on a D-dimensional sphere of radius r
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

Note that this function includes a new argument, `surface`. This argument allows us to skip the sampling of new radii after the points have been normalised to have radius equal to 1, and therefore to sample points uniformly distributed *on the surface* of the sphere, instead of inside it.

Below are two examples of 1000 random points, *inside* and *on* the sphere, respectively.

``` r
# Sample 1000 random points inside the unit sphere (D=3)
points.sphere = perfect.hypersphere(N=1000, D=3)
plot3D::points3D(x = points.sphere[,1],
                 y = points.sphere[,2],
                 z = points.sphere[,3],
                 col="cornflowerblue", pch=16, alpha=0.7, bty="g", phi=20, 
                 main="1000 random points inside the unit sphere (D = 3)")
```

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
# Sample 1000 random points on the surface of the unit sphere (D=3)
points.sphere = perfect.hypersphere(N=1000, D=3, surface=TRUE)
plot3D::points3D(x = points.sphere[,1],
                 y = points.sphere[,2],
                 z = points.sphere[,3],
                 col="cornflowerblue", pch=16, alpha=0.7, bty="g", phi=20, 
                 main="1000 random points on the unit sphere surface (D = 3)")
```

![]({{ site.baseurl }}/images/hypersphere-sampling_files/figure-markdown_github/unnamed-chunk-7-2.png)

This perfect sampling algorithm, if formally very simple, provides an example of the importance of understanding the topology of the space in which sampling is to be done, and of the elegant and powerful mathematical tricks, such as sample transformation, that are ubiquitious across all of statistical physics.

Finally, it is very important to realise that, although I have presented on this sampling algorithm mainly because of the mathematical ingenuity and beauty that it encapsulates, hypersphere sampling is a technique of enormous importance in statistical mechanics and particle physics. To mention one example, the distribution of velocities of *N* classical particles in a gas with a fixed [kinetic energy](https://en.wikipedia.org/wiki/Kinetic_energy) can be sampled as random points on the surface of a hypersphere in 3*N* dimensions. This is because the kinetic energy, which is fixed, is proportional to the sum of the squared velocities; and so the particle velocities that are sampled (each of which has *x*, *y* and *z* components) must satisfy that the sum of their squares is proportional to the kinetic energy, just as all the points on the surface of a sphere have coordinates whose sum of squares is, by definition, equal to the square of the sphere radius. Interestingly, the distribution of points in a hypersphere is also conceptually connected with [J. J. Thomson](https://en.wikipedia.org/wiki/J._J._Thomson)'s landmark 1904 [atomic model](https://en.wikipedia.org/wiki/Plum_pudding_model), which described the atom as a positively charged sphere with negative-charged electrons embedded in it. More generally, the problem of carefully defining the space to sample from is of foremost importance in [Markov chain Monte Carlo](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) sampling algorithms, being the motivation of some of the most advanced techniques of the field, such as Hamiltonian Monte Carlo and Adiabatic Monte Carlo.

**\[NB. This post is largely based on the contents of Werner Krauth's book [*Statistical Mechanics: Algorithms and Computations*](https://global.oup.com/academic/product/statistical-mechanics-algorithms-and-computations-9780198515364) (Oxford University Press), and its companion [online course](https://www.coursera.org/learn/statistical-mechanics).\]**
