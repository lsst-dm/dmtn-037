



:tocdepth: 1

.. Please do not modify tocdepth; will be fixed when a new Sphinx theme is shipped.

.. note::

   **This technote is not yet published.**

   A description of the algorithm to calculate a model of the true sky and forward model atmospheric effects to create matched templates.

.. sectnum::

Introduction
============

.. figure:: /_static/StarFast_1.3_airmass_10deg_diff.png
   :name: fig-intro_dipoles
   :target: http://target.link/url

   A difference of two simulated g-band images, observed at elevation angles 10 degrees apart, with no PSF matching to show the dipoles characteristic of mis-subtraction.

Differential Chromatic Refraction (DCR) stretches sources along the zenith direction, so that the sources appear to shift position and the PSF elongates.
If all sources had the same spectrum we could correct for the elongation of the PSF with standard PSF-matching tools.
In practice, treating DCR in this way can reduce the severity of most of the dipoles in a difference image, since many sources have roughly the same spectrum, but it will make others worse when they don't match the set selected for calibration.
In :numref:`fig-intro_dipoles` most of the dipoles point in the same direction since the majority of sources in the simulation have similar (though not identical) spectra, but there are clearly several sources with different spectra where the diploes point in the opposite direction.
These differences can't be accounted for simply with a spatially-varying PSF matching kernel, because the distribution of sources with different spectra is not a simple spatially-varying function; a quasar may be surrounded by a very red galaxy, or a cool, red, star may appear in front of a cluster of young hot stars.

Fundamentally, the problem with correcting DCR through any type of astrometric calibration or PSF matching is that a given pixel of a CCD collects photons from a different solid angle of the sky, and that solid angle depends on wavelength and the observing conditions.
Coadding or differencing two images taken under different conditions mixes flux from (slightly) different directions in every pixel, and if the conditions are different enough the quality of the images will degrade (see :numref:`fig-max_dcr` below for an estimate).
Instead, we must use our knowledge of DCR to build deep estimates of the un-refracted sky at sub-filter wavelength resolution, and forward model template images to use for image differencing.

Below, I will first briefly discuss Differential Chromatic Refraction (DCR) in the context of optical telescopes, and LSST in particular.
Next I describe my `DCR iterative forward modeling`_ approach for solving DCR with software, including the `Mathematical framework`_, `Implementation`_, `Examples with simulated images`_, and `Examples with DECam images`_.
I wrap up with a discussion of `The DCR Sky Model`_, a new data product made possible by this approach that allows me to recover `Simulated source spectra`_.



The origins of Differential Chromatic Refraction
------------------------------------------------

As with any medium, light entering Earth's atmosphere will refract.
If the incident light arrives from any direction other than straight overhead, then refraction will bend its path towards an observer such that it appears to come from a location closer to zenith (see :numref:`fig-atmos_refraction`).
While the bulk refraction is easily corrected during astrometric calibration, the index of refraction increases for bluer light and leads to an angular refraction that is wavelength-dependent.
Since optical telescopes have a finite bandwidth, this means that a point-like source will get smeared along the direction towards zenith, possibly across multiple pixels.
In :numref:`fig-max_dcr`, I calculate the absolute worst case DCR in each of the LSST bands by calculating the difference in refraction for two laser beams of light from the same position but each tuned to one of the extreme wavelengths of the filter.
The DCR in this scenario is far more severe than a realistic source spectrum, but it illustrates how the effect will scale with wavelength and zenith angle.
Please see `Appendix: Refraction calculation`_ for details on the calculation of refraction used throughout this note.

.. figure:: /_static/refraction.gif
   :name: fig-atmos_refraction
   :target: http://www.astro.caltech.edu/~mcs/CBI/pointing/

   Light incident from a source outside the Earth's atmosphere refracts and it's path bends towards zenith.
   To an observer on Earth, the source appears to be at a higher elevation angle :math:`\Delta E` than it truly is.
   For a quick and dirty approximation, :math:`\Delta E` (in arcseconds) :math:`\sim 90 - E` (in degrees) = zenith angle.
   (Image credit: Martin Shepherd, Caltech)


.. figure:: /_static/DCR_ZA-wavelength.png
   :name: fig-max_dcr
   :target: https://dmtn-017.lsst.io

   Calculation of the maximum DCR in each of the LSST filters, as a function of zenith angle. 


DCR iterative forward modeling
==============================

Because refraction is wavelength-dependent :eq:`eqn-refraction`, the smearing of astronomical sources will depend on their intrinsic spectrum.
In practice, this appears as an elongation of the measured PSF in an image, and a jitter in the source location leading to mis-subtractions in image differencing when looking for transient or variable sources.
This effect can be properly corrected when designing a telescope by incorporating an atmospheric dispersion corrector (ADC) in the optical path in front of the detector, but it cannot be removed in software once the photons have been collected.
Instead, we can mitigate the effect by building DCR-matched template images for subtraction, and use our knowledge of DCR to build a deep model of the sky from a collection of images taken under a range of observing conditions. 

In the sections below, I begin by describing the notation and mathematical framework I will use to solve the problem.
I will then lay out the details of my particular implementation of the algorithm, along with a brief discussion of approaches that I have not implemented but might be of interest in the future.
Finally, I will give examples using data simulated using `StarFast <http://dmtn-012.lsst.io/en/latest/>`_, as well as images from the `DECam HiTS <https://arxiv.org/abs/1609.03567>`_ survey.

Mathematical framework
----------------------

The point-spread function (PSF) of a point-like source is a combination of many effects arising from the telescope and its environment.
In this note, I will assume that all instrumental effects have been properly measured and accounted for, though this is of course a simplification.
Even without instrumental effects, though, we still have the atmosphere to contend with.
Turbluence in the atmosphere leads to a blurring of images (seeing), and the severity can vary significantly between nights, or even over the same night.
In my initial implementation I will ignore the effect of variable seeing, and use only observations with comparable PSFs, but in order to make full use of all the data available it will have to be addressed in the future (see `Extension to variable seeing`_).
Thus, for this initial investigation I will assume that the only effect that changes the shape of the PSF over a set of observations is DCR.

:numref:`fig-subband_diagram` below illustrates the approach of this algorithm. Since DCR arises from the change in the index of refraction of the atmosphere across a filter bandpass, if we can build a model sky in smaller sub-bands the effect is greatly reduced.
Further, if DCR is the only effect on the PSF, then we can construct these sub-band models with a small enough bandwidth such that the model to be the same for all observations regardless of airmass and parallactic angle, except for a bulk shift of the entire image.
We only measure the combined image from all sub-bands, however, so those shifted sub-band images must be added together, which results in an apparent elongation of the PSF.

.. figure:: /_static/DCR_subband_diagram.png
   :name: fig-subband_diagram

   A star observed at different parallactic and zenith angles may appear slightly elongated in the zenith direction, but this is due to DCR.
   If the full filter were split into two sub-bands (here **A** and **B**), the star would appear round(er) but shifted slightly from the position measured across the full band.
   If we assume that the 'true' image within each sub-band does not change, then the changing elongation in the **N** observed images along the zenith direction can be attributed to these shifts, which we can calculate with :eq:`eqn-refraction`.
   
Notation
^^^^^^^^

I will use matrix notation throughout this note, with vectors written as lower case letters and 2D matrices written as upper case letters.
In this context, images are written as vectors, with all of the pixel values unwrapped.
To emphasize this point, I have added arrows over the vectors, though these don't convey any additional meaning.
Image data is written as :math:`\overrightarrow{s_i}`, where the subscript :math:`i` loops over the input images, while model data is written as :math:`\overrightarrow{y_\alpha}` with the subscript :math:`\alpha` looping over sub-bands.
While it is often convenient for :math:`\overrightarrow{y_\alpha}` to have the same resolution and overall pixelization as :math:`\overrightarrow{s_i}`, in general even a non-gridded pixelization such as `HEALPix <https://healpix.jpl.nasa.gov>`_ could be used.

The matrix :math:`B_{i\alpha}` encodes the transformation due to DCR of model plane :math:`\alpha` to image :math:`i`, and the reverse transformation is written as :math:`B_{\alpha i}^\star`.
Since the sub-bands have a narrow bandwidth (but see `Finite bandwidth considerations`_ below), the effect of DCR is a uniform shift of all pixels, so:

.. math::
   :label: eqn-BB_identity

   B_{\alpha i}^\star B_{i\alpha} = \mathbb{1}

Finally, the measured PSF of each image :math:`i` is given by :math:`Q^{(i)}`, which is a matrix that does not change the size of the image.
Or, to put it in more familiar terms, it represents the convolution of any given image with the measured PSF of image :math:`i`.
Ignoring effects such as von Karmen turbulence that may lead to a wavelength-dependent PSF size, one fiducial PSF is used for all models without any index: :math:`P`.

Iterative solution derivation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The image :math:`\overrightarrow{s_i}` is the sum of all of the sub-band models (see :numref:`fig-subband_diagram`), each shifted by the appropriate amount of DCR relative to the effective wavelength of the full filter from :eq:`eqn-DCR`:

.. math::
   :label: eqn-basic_sum

    \sum_\alpha B_{i\alpha}  \overrightarrow{y_\alpha} =  \overrightarrow{s_i}

Applying the reverse shift for one sub-band :math:`\gamma`, we can re-write :eq:`eqn-basic_sum` as:

.. math::
   :label: eqn-iterative_sum

    \overrightarrow{y_\gamma} = B_{\gamma i}^\star  \overrightarrow{s_i} - B_{\gamma i}^\star  \sum_{\alpha  \neq \gamma} B_{i\alpha}  \overrightarrow{y_\alpha}  

While this may not at first appear to help, we can now solve this problem iteratively.
In each iteration, we can solve for a new set of sub-band models :math:`\overrightarrow{y_\gamma}` using the solutions from the last iteration as fixed input.
I will discuss the solutions :math:`\overrightarrow{y_\gamma}` in a later section, `The DCR Sky Model`_.

Once we have a set of model :math:`\overrightarrow{y_\gamma}`, we can use that to predict the template for a future observation :math:`k`:

.. math::
   :label: eqn-basic_template

    \parallel \overrightarrow{s_k}\!\!\parallel  = \sum_\alpha B_{k\alpha}  \overrightarrow{y_\alpha}

Finite bandwidth considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using three sub-bands for the DCR model the approximation that there is negligible DCR within a sub-band will break down for high airmass observations in the LSST g- or u-band.
For example, the differential refraction between 420nm and 460nm under typical observing conditions and at airmass 1.3 is 0.27 arcseconds, or about one LSST pixel.
With that amount of variation across a sub-band we clearly cannot expect the simple shift of :eq:`eqn-basic_sum` to work for both low and high airmass observations.
One option would be to increase the number of sub-bands of the model, but this introduces additional degrees of freedom that may not be well constrained if we have only a few high airmass observations.
Another option would be to exclude high airmass observations, but that would be very unfortunate because the large lever arm of DCR in those observations has the potential to better constrain the model (note that we could still build DCR-matched templates for high airmass observations, even if they are not used to calculate the model).
Instead, we can modify :math:`B_{i\alpha}` to include the effective smearing caused by finite bandwidth:

.. math::
   :label: eqn-modified_B

    {B}'_{i\alpha} = \int_{\alpha_0}^{\alpha_1} f(\lambda)B_{i\lambda}\mathrm{d}\lambda

And :eq:`eqn-basic_sum` becomes:

.. math::
   :label: eqn-finite_sum

    \sum_\alpha {B}'_{i\alpha}  \overrightarrow{y_\alpha} =  \overrightarrow{s_i}

This transformation is the integral of the simple shift :math:`B_{i\alpha}` across the sub-band, optionally weighted by the filter profile :math:`f(\lambda)`.
The identity :eq:`eqn-BB_identity` no longer holds, so we must either accept an approximation or attempt a deconvolution to obtain a modified :eq:`eqn-iterative_sum`.
We have studiously avoided performing any outright deconvolutions, so in my implementation I instead neglect finite bandwidth effects in the reverse transformation and set :math:`{B'}_{i\alpha}^\star = {B}_{i\alpha}^\star`.
Now :eq:`eqn-iterative_sum` becomes:

.. math::
   :label: eqn-iterative_finite_sum

    \overrightarrow{y_\gamma} = B_{\gamma i}^\star  \overrightarrow{s_i} - B_{\gamma i}^\star  \sum_{\alpha  \neq \gamma} {B}'_{i\alpha}  \overrightarrow{y_\alpha} 

And new DCR-matched template images are calculated from the resulting model with:

.. math::
   :label: eqn-finite_template

    \parallel \overrightarrow{s_k}\!\!\parallel  = \sum_\alpha {B}'_{k\alpha}  \overrightarrow{y_\alpha}

While the above approximations may be unnecessary for low airmass observations, they are also not hurt by making it.
It is possible to use equations :eq:`eqn-basic_sum` - :eq:`eqn-basic_template` for observations where the DCR within sub-bands is small and equations :eq:`eqn-finite_sum` - :eq:`eqn-finite_template` otherwise, but since the approximation improves as the amount of DCR across a sub-band decreases it seems safe to use for all observations.
We may still end up wishing to use the simple shift for low airmass observations if there is a performance difference, however.

Extension to variable seeing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While not implemented yet, there is a fairly clear path forward to extend the iterative solution from :eq:`eqn-iterative_sum` to the case where additional effects beyond DCR introduce changes to the PSF
Variable seeing is the primary concern in this case, but in principle instrumental and other effects could be accounted for in this manner as well. 

Now, we need to convolve the model with the measured PSF of the image :math:`Q^{(i)}`, and convolve the image with the fiducial PSF used for the model :math:`P`.
This modifies :eq:`eqn-basic_sum` above:

.. math::
   :label: eqn-psf_sum

   \sum_\alpha B_{i\alpha}  Q^{(i)} \overrightarrow{y_\alpha}  = P  \overrightarrow{s_i} 

Now we can once again apply the reverse shift for one sub-band, and re-write :eq:`eqn-psf_sum` as:

.. math::
   :label: eqn-psf_iterative_sum

    Q^{(i)}\overrightarrow{y_\gamma} = B_{\gamma i}^\star  P \overrightarrow{s_i} - B_{\gamma i}^\star  \sum_{\alpha  \neq \gamma} B_{i\alpha}  Q^{(i)} \overrightarrow{y_\alpha}  

Unfortunately, we now have improved estimates for :math:`Q^{(i)}\overrightarrow{y_\gamma}` when what we really want is :math:`y_\gamma`.
This problem is identical to the standard problem of image co-addition, however, so at this point we would hook in an existing algorithm for combining images with variable PSFs.

Implementation
--------------

There are four main factors to consider when turning :eq:`eqn-iterative_sum` into an effective algorithm: 
what initial solution to use as the starting point for iterations,
what conditioning to apply to the new solution found in each iteration,
how to detect and down-weight contaminated data,
and how to determine when to exit the loop.
These factors will each be described in a subsection below.


Finding the initial solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming we have no prior spectral information, the best initial guess is that all sub-bands have the same flux in all pixels.
If all model planes are equal at the start, a good guess for the flux distribution within a sub-band is the standard co-add of the input images, divided by the number of model planes being used (since those will be summed).
A proper inverse-variance weighting of the input images as part of the coaddition will help make the best estimate, and if there are many input images we could restrict the coaddition to use only those observed near zenith (with negligible DCR).
An advantage of selecting the simple coadd as the starting point, is that the solution should immediately converge if the input data exhibits no actual DCR effects, such as redder bands (*i*-band or redder for LSST) or zenith observations.
However, since this image is only the starting point of an iterative process, the final solution should not be sensitive to small errors at this stage.

Conditioning the iterative solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A common failure mode of iterative forward modeling algorithms is oscillating solutions.
In these cases, :eq:`eqn-iterative_sum` may produce intermediate solutions for :math:`\overrightarrow{y_\gamma}` with very large amplitude in one iteration, leading to very small amplitude (or negative) solutions in the next iteration, for example.
Conditioning of the solution can mitigate this sort of failure, and also help reach convergence faster.
Some useful types of conditioning include:

* Instead of taking the current solution from :eq:`eqn-iterative_sum` directly, use a weighted average of the current and last solutions.
  This eliminates most instances of oscillating solutions, since it restricts the relative change of the solution between iterations.
  In the current implementation, the weights are chosen to be the convergence metrics of the two solutions, which allows the overall solution to converge rapidly when possible but cautiously if the solutions oscillate.
  While the increased rate towards convergence is helpful for well behaved data, the greatest benefit appears when using larger numbers of subfilters.
  Adding more subfilters beyond the standard three increases the number of degrees of freedom of the problem and increases the susceptibility to unstable and diverging solutions.
  Using the dynamic weights calculated from the convergence metric allows the solution to make small improvements 

* Threshold the solutions.
  Instead of solutions diverging through solutions oscillating between iterations, the solution might 'oscillate' between model planes.
  While it is possible for all of the flux in an image to come from one single model plane, with zero from all others, there are limits.
  Solutions with more flux near a source in a single plane than the initial coadd are likely to be unphysical, and also likely to be paired with deeply negative pixels in the other planes.
  Care must be taken to avoid overly strict thresholds that impair convergence (such as applying the preceding test to even noise-like pixels), but reasonable restrictions can eliminate extreme outliers.

* Frequency regularization.
  In addition to comparing the current solution to the last or initial solutions, we could also apply restrictions on variations between model planes.
  For example, we could calculate the slope (or higher derivatives) of the spectrum for every pixel in the model across the sub-bands, and apply a threshold.
  Any values deviating more than a set amount from the line (or higher order curve) fit by that slope could be fixed to the fit instead, and minimum and maximum slopes could be set.
  While I have written an option within the DCR modeling code to enforce this sort of regularization, in practice I have found the additional benefit to be negligible when combined with the preceding forms of conditioning, and leave it turned off by default.

Weighting the input data
^^^^^^^^^^^^^^^^^^^^^^^^

Weighting of the input data takes two forms:
weights that are static properties of the image (such as the variance plane),
and dynamic weights that may change between iterations.

* Static weights.
  In most cases the static weights will be just the inverse of the variance planes of the images, and best practice is to maintain separate arrays of inverse-variance weighted image values and the corresponding inverse variance values.
  All transformations and convolutions are applied to both equally, and the properly weighted solution is the transformed weighted-image sum divided by the transformed weights sum.

* Dynamic weights.
  The simplest form of dynamic weights is a flag, which indicates whether a particular image is to be used in calculating the new solution with :eq:`eqn-iterative_sum` or not.
  If an estimated template is made for each image using the new solution and equation :eq:`eqn-basic_template`, then those templates should become a better fit to the images with every iteration.
  While it is possible to have a catastrophic failure where the model performs worse for *all* images, it might also improve for most and degrade for a few.
  For example, if there are astrometric errors for one image, the pixel-based model may be misaligned to that image and the subtraction residuals may increase between iterations.
  In that case, that image would hurt the calculation of the overall model more than the additional data was helping it, and that image should be excluded from the next iteration.
  However, in case the apparent divergence was a fluke, convergence should still be tested for that image on all subsequent iterations in case the fit improves with a better model.
  It might be possible to re-calibrate images that are flagged in this way, with the hope that an improved astrometric solution would also improve the fit to the model.

Determining the end condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Iterative forward modeling does not have an end condition that can be predetermined, and without setting a limit it would run indefinitely.
Possible end conditions include:

* Fixed time / number of iterations.
  The simplest option is to set an upper bound on the number of iterations, to ensure that the loop does exit within a finite time.
  However, the limit should be set high enough that it does not get hit in typical useage.

* Test for convergence.
  There are two types of tests that check for convergence; one that is fast and one that is accurate.
  The fast check simply compares the rate of change of the solution, and if the difference between the new and the last solution is less than a specified threshold fraction of the average of the two the solution has converged.
  The accurate check, on the other hand, creates a template for each image using :eq:`eqn-basic_template`, and calculates a convergence metric from the difference of each image from its template.
  Once the convergence metric changes by less than a specified level between iterations we can safely exit the loop since additional iterations will provide insignificant improvement.
  This is slower since templates must be created for each image, for every iteration, but that has the advantage of allowing convergence to be checked for each individual image at no extra cost (see "Dynamic weights" above).
  If the extra computational cost is not considered prohibitive, then the second test of convergence is far superior, since it is more accurate and enables additional tests and weighting.
  The level of convergence reached will clearly impact the quality of the results; see `Simulated source spectra`_ below for more analysis.

* Test for divergence.
  If a template is made for testing convergence, we should naturally test also for divergence.
  If the convergence metric actually degrades with a new iteration, that is a clear sign that some feature of the image is being modeled incorrectly and more iterations will only exacerbate the problem.
  One special case is if only some of the images degrade, because it is possible that they contain astrometric or other errors, and we could choose to continue with those observations flagged (see "Dynamic weights" above).
  Otherwise, it is safest to exit immediately and discard the solution from the current iteration, using the last solution from before it started to diverge instead.

  * One possible modification is to calculate a spatially-varying convergence metric, and mask regions that degrade in future iterations.
    This allows an improved solution to be found even if one very bright feature (such as an improperly-masked saturated star or cosmic ray) is modeled incorrectly.


Note that if a convergence test is used, it should only be allowed to exit the loop after a minimum number of iterations have passed.
It will depend on how the convergence metric is calculated and the choice of initial solution, but the first iteration can show a slight degradation of convergence.

Examples with simulated images
------------------------------

To test the above algorithm, I first ran it on images simulated using `StarFast <http://dmtn-012.lsst.io/en/latest/>`_.
As shown in :numref:`fig-sim_image` below, these images contain a moderately crowded field of stars (no galaxies) with Kolmogorov double-gaussian PSFs, realistic SEDs, photon shot noise, and no other effects other than DCR.
In this example, the model is built using three frequency planes and eight input simulations of the field, with airmass ranging between 1.0 and 2.0 (*not* including the simulated observation shown in :numref:`fig-sim_image`).
I use :eq:`eqn-basic_template` to build a DCR-matched template for the simulated observation (:numref:`fig-sim_template`), and subtract this template to make the difference image (:numref:`fig-sim_template_diff`).
For comparison, in :numref:`fig-sim_image_diff` I subtract a second simulated image generated with the same field 10 degrees closer to zenith.
Note that the noise in :numref:`fig-sim_template_diff` is :math:`\sim \sqrt{2}` lower than in :numref:`fig-sim_image_diff`, because the template is a form of coadded image with significantly lower noise than an individual exposure.

.. figure:: /_static/simulations/simulated_108_image.png
   :name: fig-sim_image

   Simulated g-band image with airmass 1.3

.. figure:: /_static/simulations/simulated_image_108_template.png
   :name: fig-sim_template

   The DCR-matched template for the simulated image in :numref:`fig-sim_image`.

.. figure:: /_static/simulations/simulated_image_108_template_difference.png
   :name: fig-sim_template_diff

   Image difference of the simulated image :numref:`fig-sim_image` with its template from :numref:`fig-sim_template`.
   The image difference is taken as a direct pixel subtraction, rather than an Alard & Lupton or ZOGY style subtraction.
   While those more sophisticated styles of subtraction are used in production, here the emphasis is on comparing the raw differences, before PSF matching or other corrections.

.. figure:: /_static/simulations/simulated_image_108_112_difference.png
   :name: fig-sim_image_diff

   Image difference of :numref:`fig-sim_image` with another simulation of the same field 10 degrees closer to zenith (airmass 1.22).
   As in :numref:`fig-sim_template_diff` above, this image difference is a direct pixel subtraction and not an Alard & Lupton or ZOGY style subtraction.

Dipole mitigation
^^^^^^^^^^^^^^^^^^

A quantitative assessment of the quality of the template from :numref:`fig-sim_template` is the number of false detection in the image difference, :numref:`fig-sim_template_diff`.
Since there are no moving or variable sources in these simulated images, any source detected in the image difference is a false detection.
From comparing  :numref:`fig-sim_template_diff` and :numref:`fig-sim_image_diff` it is clear by eye that the DCR-matched template produces fewer dipoles, and we can see in :numref:`fig-dcr_dipoles` that this advantage holds regardless of airmass or elevation angle.
The DCR-matched template appears to perform about as well as an exposure taken within 5 degrees of the science image, but recall that the noise - and therefore the detection limit - is higher when using a single image as a template.


.. figure:: /_static/simulations/DCR_dipole_plot.png
   :name: fig-dcr_dipoles

   Plot of the number of false detections for simulated images as a function of airmass, for different templates.

Examples with DECam images
--------------------------

For a more rigorous test, I have also built DCR-matched templates for `DECam HiTS <https://arxiv.org/abs/1609.03567>`_ observations, which were calibrated and provided by Francisco Förster. 
For these images I used the implementation outlined above using the simplified equation :eq:`eqn-iterative_sum`, despite the images having variable seeing.
Because I have not yet implemented :eq:`eqn-psf_iterative_sum` using measured PSFs for each image, I have simply excluded observations with PSF FWHWs greater than 4 pixels (2.5 - 4 pixel widths are common). 
However, it should be noted that the images with PSFs at the larger end of that range are not as well matched by their templates, so restricting the input images going into the model to those with good seeing may not be sufficient in the future.

.. figure:: /_static/Decam/0410998_image.png
   :name: fig-decam_image

   Decam observation 410998, 'g'-band with airmass 1.33.

.. figure:: /_static/Decam/0410998_template.png
   :name: fig-decam_template

   The DCR-matched template for 410998, constructed from 12 observations with airmass ranging between 1.13 and 1.77.
   Note that the noise level is significantly decreased.

.. figure:: /_static/Decam/0410998-template_difference.png
   :name: fig-decam_template_diff

   Difference image of :numref:`fig-decam_image` - :numref:`fig-decam_template`.

.. figure:: /_static/Decam/0410998-0411232_difference.png
   :name: fig-decam_image_diff

   Image difference of :numref:`fig-decam_image` with a second DECam observation taken an hour later and approximately 10 degrees closer to zenith (at airmass 1.23).

From the above images we can make several observations.
First is that the DCR-matched template in :numref:`fig-decam_template_diff` performs at least as well as a well-matched reference image taken within an hour and 10 degrees of our science image (:numref:`fig-decam_image_diff`).
The DCR-matched template, however, has significantly reduced noise and artifacts (such as cosmic rays), since we are able to coadd 12 observations because we are not constrained by needing to match observing conditions.
Thus, a DCR-matched template calculated at zenith should provide the cleanest and deepest estimate of the static sky possible from the LSST survey, since we can coadd all images without sacrificing image quality or resolution.

While constructing the DCR sky model for a full DECam CCD was time consuming on my laptop (taking about 10 minutes on one core), forward modeling a DCR-matched template from the model takes only 1-2 seconds on the same machine.
Adding variable the PSFs from :eq:`eqn-psf_iterative_sum` will increase the time required to calculate the DCR sky model, but should not increase the time to forward model templates.

The DCR Sky Model
=================

The DCR sky model :math:`\overrightarrow{y_\alpha}` from :eq:`eqn-iterative_sum`  consists of a deep coadd in each subfilter :math:`\alpha`.
As with other coadds, the images are defined on instrument-agnostic tracts and patches of the sky, and must be warped to the WCS of the science image after constructing templates with :eq:`eqn-basic_template`.
While the sky model was designed for quickly building matched templates for image differencing, it is an interesting data product in its own right.
For example, we can run source detection and measurement on each sub-filter image :math:`\overrightarrow{y_\alpha}` and measure the spectra of sources within a single band.
An example visualization of the sub-bands of the DCR sky model is in :numref:`fig-filled_footprints`, below, while a more detailed look at spectra is in `Simulated source spectra`_.
This view can help identify sources with steep or unusual spectra, such as quasars with high emission in a narrow band, and could be used to help deblending and star-galaxy separation.
However, because of the inherent assumption that the true sky is static it cannot estimate the spectrum of transient or variable sources.

.. figure:: /_static/sim_filled_footprints_color2.png
   :name: fig-filled_footprints

   Source measurements in three sub-bands of the DCR sky model are converted to RGB values and used to fill the footprints of detected sources.
   The combined full-band model is displayed behind the footprint overlay.

Simulated source spectra
------------------------

As mentioned above, we can look in more detail at the spectra of individual sources.
The simulated images contain roughly 2500 stars ranging from type M to type B, with a distribution that follows local abundances.
Each star is simulated at high frequency resolution using Kurucz SED profiles and propagated through a model of the LSST g-band filter bandpass, so it is straightforward to compare the measured spectrum of a source to its input spectrum once it is matched.
To measure the subfilter source spectra, I run a modified version of multiband photometry from the LSST software stack.
This performs source detection on each subfilter coadd image, merges the detections from all subfilters, and performs forced photometry on each.
A few example comparison spectra for a range of stellar types are in :numref:`fig-starspectrum_bright1` - :numref:`fig-starspectrum_faint` below.

.. figure:: /_static/spectra/new/star_spectra_sim436_bright.png
   :name: fig-starspectrum_bright1

   Example input spectrum for a type F star with surface temperature ~7130K (solid blue).
   The flux measured in each sub-band is marked with a with red '+', and the average values of the simulated spectrum across each subfilter is marked with a blue 'x' for comparison.

.. figure:: /_static/spectra/new/star_spectra_sim423_bright.png
   :name: fig-starspectrum_bright2

   Example input spectrum for a type F star with surface temperature ~6370K (solid blue).
   The symbols are as :numref:`fig-starspectrum_bright1`.

.. figure:: /_static/spectra/new/star_spectra_sim324_mid.png
   :name: fig-starspectrum_mid

   Example input spectrum for a type G  star with surface temperature ~5420K (solid blue).
   The symbols are as :numref:`fig-starspectrum_bright1`.

.. figure:: /_static/spectra/new/star_spectra_sim262_mid.png
   :name: fig-starspectrum_mid2

   Example input spectrum for a type K  star with surface temperature ~4610K (solid blue).
   The symbols are as :numref:`fig-starspectrum_bright1`.

.. figure:: /_static/spectra/new/star_spectra_sim006_faint.png
   :name: fig-starspectrum_faint

   Example input spectrum for a type M star with surface temperature ~3620K (solid blue).
   The symbols are as :numref:`fig-starspectrum_bright1`.

While the above spectra are representative of the typical stars in the simulation, it is helpful to look at the entire set.
For this comparison, in :numref:`fig-colorcolor01` below I plot the simulated color between the blue and the red subfilters against the measured color, ignoring the center subfilter.
For this example, I let the forward modeling proceed until it reached 1% convergence, which took 8 iterations.
While there is a clear correlation between the measured and simulated spectra the slope is clearly off, with the measurements flatter than the simulations.
If a 0.1% convergence threshold is used, however, then after 24 iterations the agreement improves (:numref:`fig-colorcolor001`).
The example spectra in :numref:`fig-starspectrum_bright1` - :numref:`fig-starspectrum_faint` above used the 0.1% threshold.

.. figure:: /_static/spectra/color-color_sim01.png
   :name: fig-colorcolor01

   Simulated - measured color of detected sources within g-band with a 1% convergence threshold.

.. figure:: /_static/spectra/color-color_sim001.png
   :name: fig-colorcolor001

   Simulated - measured color of detected sources within g-band with a 0.1% convergence threshold.

.. Comment out this figure for now
  .. figure:: /_static/spectra/color-color_sim003.png
     :name: fig-colorcolor003

     Simulated - measured color of detected sources within g-band with a 0.3% convergence threshold.

Appendix: Refraction calculation
==================================

While the true density and index of refraction of air varies significantly with altitude, I will follow :cite:`Stone1996` in approximating it as a simple exponential profile in density that depends only on measured surface conditions.
While this is an approximation, it is reportedly accurate to better than 10 milliarcseconds for observations within 65 degrees of zenith, which should be sufficient for normal LSST operations.

The refraction of monochromatic light is given by

.. math::
   :label: eqn-refraction

   R(\lambda) &= r_0 n_0(\lambda) \sin z_0 \int_1^{n_0(\lambda)} \frac{dn}{n \left(r^2n^2 -r_0^2n_0(\lambda)^2\sin^2z_0\right)^{1/2}} \nonumber\\
    &\simeq \kappa (n_0(\lambda) - 1) (1 - \beta) \tan z_0 - \kappa (1 - n_0(\lambda)) \left(\beta - \frac{n_0(\lambda) - 1}{2}\right) \tan^3z_0

where :math:`n_0(\lambda)`, :math:`\kappa`, and :math:`\beta` are given by equations :eq:`eqn-n_lambda`, :eq:`eqn-kappa`, and :eq:`eqn-beta` below. 
The differential refraction relative to a reference wavelength is simply:

.. math::
   :label: eqn-DCR

   \Delta R(\lambda) = R(\lambda) - R(\lambda_{ref})

The index of refraction as a function of wavelength :math:`\lambda` (in Angstroms) can be calculated from the relative humidity (:math:`RH`, in percent), surface air temperature (:math:`T`, in Kelvin), and pressure (:math:`P_s` in millibar):

.. math::
   :label: eqn-n_lambda

   n_0( \lambda ) &=\:& 1 + \Delta n_s + \Delta n_w \\
   \\
   \Delta n_s &=\:& \bigg(2371.34 + \frac{683939.7}{130 -\sigma(\lambda)} + \frac{4547.3}{38.9 - \sigma(\lambda)^2}\bigg) D_s \times 10^{-8} \\
   \\
   \Delta n_w &=\:& \big(6487.31 + 58.058 \sigma(\lambda)^2 - 0.71150\sigma(\lambda)^4 + 0.08851\sigma(\lambda)^6\big) D_w \times 10^{-8} \\
   \\
   \sigma(\lambda) &=\:& 10^4/\lambda \;\;\;( \mu m^{-1})
   

Where the density factors for water vapor :math:`D_w` and dry air :math:`D_s` are given by :eq:`eqn-D_w` and :eq:`eqn-D_s` (from  :cite:`Owens67`), and the water vapor pressure :math:`P_w` is calculated from the relative humidity :math:`RH` with :eq:`eqn-P_w`.

.. math::
   :label: eqn-D_w

   D_w = \bigg[1+P_w (1+3.7\times10^{-4}P_w)\bigg(-2.37321\times 10^{-3} + \frac{2.23366}{T} - \frac{710.792}{T^2} + \frac{7.75141\times 10^4}{T^3}\bigg)\bigg] \frac{P_w}{T} 

.. math::
   :label: eqn-D_s

   D_s = \bigg[1 + (P_s - P_w) \bigg( 57.90 \times 10^{-8} -  \frac{9.3250\times 10^{-4}}{T} + \frac{0.25844}{T^2}\bigg)\bigg] \frac{P_s - P_w}{T}

.. math::
   :label: eqn-P_w

   P_w = RH\times 10^{-4}\times e^{(77.3450 + 0.0057 T - 7235.0/T)}/T^{8.2}


The ratio of local gravity at the observing site to :math:`g= 9.81 m/s^2` is given by

.. math::
   :label: eqn-kappa  

    \kappa = g_0/g = 1 + 5.302\times 10^{-3} \sin^2\phi - 5.83\times 10^{-6} \sin^2(2\phi) - 3.15\times 10^{-7} h \label{eqn:kappa}

By assuming an exponential density profile for the atmosphere, the ratio :math:`\beta` of the scale height of the atmosphere to radius of the observing site from the Earth's core can be approximated by:

.. math::
   :label: eqn-beta

   \beta &= \frac{1}{R_\oplus}\int_{0}^\infty \frac{\rho}{\rho_0} dh \nonumber \\
    &\simeq \frac{P_s}{\rho_0g_0 R_\oplus} = \frac{k_BT}{m g_0 R_\oplus} \nonumber \\
    &=  4.5908\times 10^{-6} T 

where :math:`m` is the average mass of molecules in the atmosphere, :math:`R_\oplus` is the radius of the Earth, :math:`k_B` is the Boltzmann constant, and :math:`g_0` is the acceleration due to gravity at the Earth's surface.

References
==========

.. bibliography:: DCR_references.bib
