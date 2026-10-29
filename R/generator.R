# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Joel Sjoberg; Copyright (c) 2026 Bogdan Oancea

# Seeded calls restore the caller's RNG state; stream environments retain their
# own state between batches. R's generator intentionally differs from NumPy.
.radar_with_rng <- function(rng, expr) {
  if (is.null(rng)) return(force(expr))
  stream <- is.environment(rng) && inherits(rng,"radar_rng")
  if(!stream && (!is.numeric(rng) || length(rng)!=1L || !is.finite(rng)))
    stop("rng must be NULL, a finite numeric seed or a RADAR random stream",call.=FALSE)
  had <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had) old <- get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  on.exit({
    if(stream && exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rng$state <- get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
    if(had) assign(".Random.seed",old,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
  })
  if(stream) assign(".Random.seed",rng$state,envir=.GlobalEnv) else set.seed(as.integer(rng))
  force(expr)
}
.radar_rng <- function(seed=NULL) {
  e <- new.env(parent=emptyenv()); class(e)<-"radar_rng"
  if(is.null(seed)) seed <- sample.int(.Machine$integer.max,1L)
  e$state <- .radar_with_rng(seed,get(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
  e
}
.gen_integer <- function(lo,hi,n=1L) {
  if(hi<=lo) stop("the upper sampling bound must exceed the lower bound",call.=FALSE)
  as.integer(floor(stats::runif(n,lo,hi)))
}
.gen_length <- function(length) {
  if(!is.numeric(length)||length(length)!=1L||!is.finite(length)||length<2L||length!=floor(length))
    stop("length must be an integer of at least 2",call.=FALSE)
  as.integer(length)
}
.gen_normalize <- function(x,eps=1e-11) (x-min(x))/(max(x)-min(x)+eps)

# Rational approximation to the complex probability function, evaluated using
# a Fourier expansion. The Voigt profile is Re(w(z))/(sigma*sqrt(2*pi)).
.gen_wofz <- local({
  n <- 32L; m <- 2L*n; k <- seq.int(-m+1L,m-1L); ell <- sqrt(n/sqrt(2))
  t <- ell*tan(k*pi/(2*m)); f <- c(0,exp(-t*t)*(ell*ell+t*t))
  f <- f[c(seq.int(2*n+1L,4*n),seq_len(2*n))]
  a <- rev((Re(stats::fft(f))/(2*m))[seq.int(2L,n+1L)])
  function(z) {
    zz <- (ell+1i*z)/(ell-1i*z); p <- rep(0i,length(z))
    for(ai in a) p <- p*zz+ai
    2*p/(ell-1i*z)^2 + 1/sqrt(pi)/(ell-1i*z)
  }
})
.gen_voigt <- function(x,sigma,gamma) {
  .gen_normalize(Re(.gen_wofz((x+1i*gamma)/(sigma*sqrt(2))))/(sigma*sqrt(2*pi)))
}

#' Draw elementary spectral backgrounds
#' @param length Number of samples.
#' @param sin_max_periodicity Maximum sine periodicity.
#' @param pol_max_degree Maximum polynomial degree.
#' @param eps Stabilizer for the Lorentzian baseline.
#' @param rng Optional integer seed or internal random stream.
#' @return A list of linear, sinusoidal, polynomial, Gaussian and Lorentzian curves.
#' @export
make_baselines <- function(length,sin_max_periodicity=8,pol_max_degree=10,eps=1e-8,rng=NULL) {
  length <- .gen_length(length)
  if(pol_max_degree<2L) stop("pol_max_degree must be at least 2",call.=FALSE)
  .radar_with_rng(rng,{
    scalers <- stats::runif(4); x <- seq(0,1,length.out=length)
    ends <- stats::runif(2); lin <- seq(ends[1],ends[2],length.out=length)
    offset <- stats::runif(1,0,2*pi); stretch <- stats::runif(1,.1,pi*sin_max_periodicity)
    sin <- .gen_normalize(sin(x*stretch+offset))*scalers[1]
    sin <- sin+stats::runif(1,0,max(0,1-max(sin)))
    degree <- .gen_integer(2,pol_max_degree+1L)
    pol <- numeric(length)
    for(i in 0:degree) pol <- pol+stats::runif(1,-1,1)*x^i
    pol <- .gen_normalize(pol)*scalers[2]
    pol <- pol+stats::runif(1,0,max(0,1-max(pol)))
    a <- stats::runif(1,-1,1)
    b <- if(stats::runif(1)>.5) stats::runif(1,.8,1) else stats::runif(1,0,.2)
    cc <- stats::runif(1,.3,1)
    gauss <- .gen_normalize(a*exp(-(x-b)^2/(2*cc^2)))*scalers[3]
    gauss <- gauss+stats::runif(1)*max(0,1-max(gauss))
    a <- stats::runif(1,.1,.5); if(stats::runif(1)>.5) a <- -a
    b <- if(stats::runif(1)>.5) stats::runif(1,.8,1) else stats::runif(1,0,.2)
    lor <- .gen_normalize((1/pi)*a/(a*a+(x-b)^2+eps))*scalers[4]
    lor <- lor+stats::runif(1)*max(0,1-max(lor))
    list(linear=lin,sinusoidal=sin,polynomial=pol,gaussian=gauss,lorentzian=lor)
  })
}

#' Draw a random spectral baseline
#' @param length Number of samples.
#' @param rng Optional integer seed or internal random stream.
#' @return Numeric baseline vector.
#' @export
generate_baseline <- function(length,rng=NULL) {
  .radar_with_rng(rng,{
    bls <- make_baselines(length)
    mask <- .gen_integer(0,2,length(bls))
    Reduce(`+`,Map(`*`,bls,mask))
  })
}

#' Draw a random Raman band
#' @param length Width in samples, at least one.
#' @param extremity Half-width of the band coordinate interval.
#' @param rng Optional integer seed or internal random stream.
#' @return Gaussian, Lorentzian or Voigt band of random height.
#' @export
get_random_peak <- function(length=16L,extremity=.5,rng=NULL) {
  if(length<1 || length!=floor(length)) stop("peak length must be a positive integer",call.=FALSE)
  .radar_with_rng(rng,{
    line <- seq(-extremity,extremity,length.out=length)
    rbf_epsilon <- .gen_integer(5,10); lorentz_epsilon <- stats::runif(1,.1,.3)
    sigma <- stats::runif(1,.01,.3); gamma <- stats::runif(1,0,.1); u <- stats::runif(1)
    if(u<.33) p <- exp(-(line*rbf_epsilon)^2) else if(u<.66) {
      p <- (1/pi)*lorentz_epsilon/(lorentz_epsilon^2+line^2); p <- p-min(p)
      if(max(p)>0) p <- p/max(p) else p[] <- 1
    } else p <- .gen_voigt(line,sigma,gamma)
    p*stats::runif(1)
  })
}

#' Draw a normalized peak-only spectrum
#' @param length Number of samples.
#' @param max_num_peaks Exclusive upper bound on the number of bands.
#' @param min_peak_width,max_peak_width Band width bounds, upper exclusive.
#' @param eps Stabilizing denominator.
#' @param rng Optional integer seed or internal random stream.
#' @return A numeric spectrum whose maximum is at most one.
#' @export
generate_peaks <- function(length=1000L,max_num_peaks=10L,min_peak_width=16L,max_peak_width=100L,eps=1e-7,rng=NULL) {
  length <- .gen_length(length)
  if(min_peak_width<1L || max_peak_width<=min_peak_width) stop("invalid peak width bounds",call.=FALSE)
  .radar_with_rng(rng,{
    spectrum <- numeric(length)
    n <- if(max_num_peaks<=1L) 1L else .gen_integer(1,max_num_peaks)
    for(i in seq_len(n)) {
      width <- .gen_integer(min_peak_width,max_peak_width); mid <- width%/%2L
      loc <- .gen_integer(mid,max(length-mid,mid+1L))
      peak <- get_random_peak(width)
      start <- loc-mid+1L; end <- min(start+width-1L,length)
      if(start<=end) spectrum[start:end] <- spectrum[start:end]+peak[seq_len(end-start+1L)]
    }
    spectrum[c(1L,length)] <- 0
    spectrum/(max(spectrum)+eps)
  })
}

#' Draw cosmic-ray spikes
#' @param length Number of samples.
#' @param max_num Exclusive upper bound on the number of spikes.
#' @param rng Optional integer seed or internal random stream.
#' @return Numeric vector of single-pixel spikes.
#' @export
generate_cosmic_rays <- function(length=1000L,max_num=30L,rng=NULL) {
  length <- .gen_length(length)
  .radar_with_rng(rng,{
    spectrum <- numeric(length); n <- if(max_num>0L) .gen_integer(0,max_num) else 0L
    if(n) {
      locs <- .gen_integer(0,length,n)+1L; heights <- stats::runif(n)
      # Assign in order so repeated positions behave like NumPy assignment.
      for(i in seq_len(n)) spectrum[locs[i]] <- heights[i]
    }
    spectrum
  })
}

#' Draw white spectral noise
#' @param length Number of samples.
#' @param scaling_factor Maximum Gaussian noise standard deviation.
#' @param rng Optional integer seed or internal random stream.
#' @return Gaussian noise with zero endpoints.
#' @export
generate_noise <- function(length=1000L,scaling_factor=.1,rng=NULL) {
  length <- .gen_length(length)
  .radar_with_rng(rng,{
    spectrum <- stats::rnorm(length,0,stats::runif(1)*scaling_factor)
    spectrum[c(1L,length)] <- 0; spectrum
  })
}

#' Draw four unscaled spectral components
#' @param length Number of samples.
#' @param min_width,max_width Band width bounds, upper exclusive.
#' @param num_peaks Exclusive upper bound on the number of bands.
#' @param rng Optional integer seed or internal random stream.
#' @return Named list containing all four component vectors.
#' @export
generate_spectrum <- function(length,min_width=16L,max_width=256L,num_peaks=10L,rng=NULL) {
  length <- .gen_length(length)
  .radar_with_rng(rng,{
    baseline <- generate_baseline(length); rays <- generate_cosmic_rays(length); noise <- generate_noise(length)
    if(length<500L) { min_width <- max(as.integer(length/32),1L); max_width <- max(as.integer(length/2),min_width+1L) }
    peaks <- generate_peaks(length,num_peaks,min_width,max_width)
    list(baseline=baseline,cosmic_rays=rays,noise=noise,peaks=peaks)
  })
}

#' Draw one synthetic spectrum and its targets
#' @param length Number of samples.
#' @param min_width,max_width Band width bounds, upper exclusive.
#' @param max_num_peaks Exclusive upper bound on band count.
#' @param increase_noise Amplify noise and min-max normalize, as in the test generator.
#' @param rng Optional integer seed or internal random stream.
#' @param legacy_scaling Reproduce the inconsistent original training input scaling.
#' @return List with spectrum and targets. Except with legacy_scaling=TRUE,
#'   the spectrum equals the sum of its four targets.
#' @export
generate_sample <- function(length=1024L,min_width=5L,max_width=256L,max_num_peaks=10L,
                            increase_noise=FALSE,rng=NULL,legacy_scaling=FALSE) {
  length <- .gen_length(length)
  .radar_with_rng(rng,{
    repeat {
      t <- generate_spectrum(length,min_width,max_width,max_num_peaks)
      scalers <- .gen_integer(0,10,3L)
      if(sum(scalers)==0) next
      scalers <- scalers/sum(scalers)
      bl <- t$baseline*scalers[1]; ray <- t$cosmic_rays*scalers[2]; peak <- t$peaks*scalers[3]
      noise <- t$noise; s <- bl+ray+peak
      if(max(s)>0) break
    }
    if(!increase_noise) {
      maxim <- max(s); bl <- bl/maxim; ray <- ray/maxim; peak <- peak/maxim
      if(stats::runif(1)<.05) noise[] <- 0
      noise <- noise*scalers[3]
      s <- (if(legacy_scaling) s else bl+ray+peak)+noise
    } else {
      noise <- noise*.gen_integer(3,10)*scalers[3]; s <- s+noise
      minimum <- min(s); s <- s-minimum; bl <- bl-minimum; maxim <- max(s)
      s <- s/maxim; bl <- bl/maxim; ray <- ray/maxim; peak <- peak/maxim; noise <- noise/maxim
    }
    list(spectrum=s,targets=list(baseline=bl,cosmic_rays=ray,noise=noise,peaks=peak))
  })
}

#' Construct a synthetic spectral dataset
#' @param X Input matrix with one spectrum per row.
#' @param baseline,cosmic_rays,noise,peaks Ground-truth component matrices.
#' @param params Generator parameters.
#' @return A synthetic_dataset object with X, spectra, the components, and params.
#' @export
synthetic_dataset <- function(X,baseline,cosmic_rays,noise,peaks,params=list()) {
  X <- .radar_matrix(X)
  vals <- list(baseline=baseline,cosmic_rays=cosmic_rays,noise=noise,peaks=peaks)
  for(c in COMPONENTS) {
    vals[[c]] <- .radar_matrix(vals[[c]],c)
    if(!identical(dim(vals[[c]]),dim(X))) stop("all targets must match X",call.=FALSE)
  }
  structure(c(list(X=X,spectra=X),vals,list(params=params)),class="synthetic_dataset")
}
#' @export
`$.synthetic_dataset` <- function(x,name) {
  v <- unclass(x)
  switch(name,length=ncol(v$X),n_spectra=nrow(v$X),targets=v[COMPONENTS],v[[name]])
}
#' @export
length.synthetic_dataset <- function(x) nrow(x$X)
#' @export
print.synthetic_dataset <- function(x,...) {
  cat(sprintf("RADAR synthetic dataset: %d spectra x %d samples\n",nrow(x$X),ncol(x$X))); invisible(x)
}
#' @export
`[.synthetic_dataset` <- function(x,i,...) {
  if(missing(i)) return(x)
  vals <- lapply(c("X",COMPONENTS),function(c) x[[c]][i,,drop=FALSE]); names(vals)<-c("X",COMPONENTS)
  do.call(synthetic_dataset,c(vals,list(params=x$params)))
}

#' Generate synthetic Raman spectra
#' @param n Number of spectra.
#' @param length Samples per spectrum.
#' @param min_width,max_width Band width bounds, upper exclusive.
#' @param max_num_peaks Exclusive upper bound on band count.
#' @param increase_noise Use amplified-noise test settings.
#' @param seed Optional integer random seed, independent of the caller's RNG state.
#' @param legacy_scaling Reproduce the original inconsistent training scaling.
#' @return A synthetic_dataset containing spectra, components, and generation parameters.
#' @details Seeded generation is reproducible within R.
#' All arrays use R double precision. Voigt bands use a native complex rational approximation.
#' @export
generate_dataset <- function(n,length=1024L,min_width=10L,max_width=300L,max_num_peaks=10L,
                             increase_noise=TRUE,seed=NULL,legacy_scaling=FALSE) {
  length <- .gen_length(length)
  if(!is.numeric(n)||length(n)!=1L||!is.finite(n)||n<1L||n!=floor(n)) stop("n must be a positive integer",call.=FALSE)
  .radar_with_rng(seed,{
    X <- matrix(0,n,length); vals <- stats::setNames(lapply(COMPONENTS,function(c) matrix(0,n,length)),COMPONENTS)
    for(i in seq_len(n)) {
      sample <- generate_sample(length,min_width,max_width,max_num_peaks,increase_noise,legacy_scaling=legacy_scaling)
      X[i,] <- sample$spectrum
      for(c in COMPONENTS) vals[[c]][i,] <- sample$targets[[c]]
    }
    params <- list(n=n,length=length,min_width=min_width,max_width=max_width,max_num_peaks=max_num_peaks,
                   increase_noise=increase_noise,legacy_scaling=legacy_scaling,
                   seed=if(is.environment(seed)) "stream" else seed)
    do.call(synthetic_dataset,c(list(X=X),vals,list(params=params)))
  })
}
#' @rdname generate_dataset
#' @export
generate <- generate_dataset

#' Generate independent training and test datasets
#' @param n_train,n_test Numbers of training and test spectra.
#' @param length Samples per spectrum.
#' @param seed Optional integer seed.
#' @param train_kwargs,test_kwargs Lists of generation parameters overriding the defaults.
#' @return Named list with train and test datasets.
#' @export
train_test_datasets <- function(n_train,n_test,length=1024L,seed=0L,train_kwargs=NULL,test_kwargs=NULL) {
  seeds <- .radar_with_rng(seed,sample.int(.Machine$integer.max,2L))
  tr <- utils::modifyList(list(min_width=5L,max_width=256L,max_num_peaks=5L,increase_noise=FALSE),if(is.null(train_kwargs)) list() else train_kwargs)
  te <- utils::modifyList(list(min_width=10L,max_width=300L,max_num_peaks=10L,increase_noise=TRUE),if(is.null(test_kwargs)) list() else test_kwargs)
  tr$n <- n_train; tr$length <- length; tr$seed <- seeds[1]
  te$n <- n_test; te$length <- length; te$seed <- seeds[2]
  list(train=do.call(generate_dataset,tr),test=do.call(generate_dataset,te))
}

#' Convert synthetic data to native model training arrays
#' @param x A synthetic_dataset object.
#' @return List with X of dimensions n x length x 1 and targets in model-output order.
#' @export
as_training_data <- function(x) {
  if(!inherits(x,"synthetic_dataset")) stop("x must be a synthetic_dataset",call.=FALSE)
  list(X=array(x$X,dim=c(nrow(x$X),ncol(x$X),1L)),targets=x$targets)
}

#' Create a stream of synthetic training batches
#' @param length Samples per spectrum.
#' @param batch_size Spectra per batch.
#' @param min_width,max_width Band width bounds, upper exclusive.
#' @param max_num_peaks Exclusive upper bound on band count.
#' @param increase_noise Use amplified-noise test settings.
#' @param seed Optional integer seed.
#' @param legacy_scaling Reproduce original inconsistent training input scaling.
#' @return A function. Each invocation returns a list with a three-dimensional
#'   input array X, target matrices, and the corresponding dataset.
#' @export
batch_generator <- function(length=750L,batch_size=100L,min_width=5L,max_width=256L,
                            max_num_peaks=5L,increase_noise=FALSE,seed=NULL,legacy_scaling=FALSE) {
  stream <- .radar_rng(seed)
  force(length);force(batch_size);force(min_width);force(max_width);force(max_num_peaks);force(increase_noise);force(legacy_scaling)
  function() {
    data <- generate_dataset(batch_size,length,min_width,max_width,max_num_peaks,increase_noise,stream,legacy_scaling)
    out <- as_training_data(data); out$dataset <- data; out
  }
}
