test_that("core checks pass", {
  
  # Flat spectra scale to zero and the baseline absorbs the original offset.
  x <- matrix(c(1,3,5,2,2,2),nrow=2,byrow=TRUE)
  z <- minmax_normalize(x)
  expect_true(all(identical(as.numeric(z$normalized[1,]),c(0,.5,1)),all(z$normalized[2,]==0)))
  back <- denormalize_components(list(baseline=z$normalized/2,peaks=z$normalized/2),z$scaling)
  expect_true(all(max(abs(back$baseline+back$peaks-x))<1e-12))
  expect_error(minmax_normalize(1:5))
  
  # Data retain all four targets, repeat deterministically, and isolate seeded calls.
  set.seed(42); old_rng <- .Random.seed
  first <- generate(10,length=128,seed=7)
  expect_true(all(identical(old_rng,.Random.seed)))
  second <- generate_dataset(10,length=128,seed=7)
  third <- generate(10,length=128,seed=8)
  expect_true(all(identical(first$X,second$X),!identical(first$X,third$X),
            identical(dim(first$X),c(10L,128L)),length(first)==10L,first$length==128L,
            max(abs(first$X-Reduce(`+`,first$targets)))<1e-12))
  expect_true(all(identical(dim(first[1:3]$X),c(3L,128L)),all(vapply(first$targets,is.matrix,logical(1)))))
  for(increased in c(FALSE,TRUE)) for(i in seq_len(10)) {
    sample <- generate_sample(64,increase_noise=increased,rng=i)
    expect_true(all(all(is.finite(sample$spectrum)),max(abs(sample$spectrum-Reduce(`+`,sample$targets)))<1e-12))
    if(increased) expect_true(all(abs(min(sample$spectrum))<1e-12,abs(max(sample$spectrum)-1)<1e-12))
  }
  legacy <- generate(10,length=128,increase_noise=FALSE,legacy_scaling=TRUE,seed=3)
  expect_true(all(max(abs(legacy$X-Reduce(`+`,legacy$targets)))>.001))
  for(b in make_baselines(100,rng=1)) expect_true(all(min(b)>=-1e-12,max(b)<=1+1e-12))
  expect_error(generate(2,length=1)); expect_error(generate(0)); expect_error(generate_peaks(min_peak_width=3,max_peak_width=3))
  
  # The native Voigt calculation matches independently computed SciPy values.
  v <- getFromNamespace(".gen_voigt","radar.raman")(c(-.5,-.1,0,.1,.5),.03,.01)
  expect_true(all(max(abs(v-c(0,0.046838011255973346,0.9999999999990361,0.046838011255973346,0)))<1e-7))
  
  # Successive batches advance their own stream and leave the caller's RNG intact.
  stream <- batch_generator(length=32,batch_size=2,seed=3)
  stream_again <- batch_generator(length=32,batch_size=2,seed=3)
  old_rng <- .Random.seed; b1 <- stream(); b2 <- stream()
  expect_true(all(identical(old_rng,.Random.seed),identical(dim(b1$X),c(2L,32L,1L)),
            identical(b1$X,stream_again()$X),!identical(b1$X,b2$X),length(b1$targets)==4L))
  expected_batch <- generate(2,length=32,min_width=5,max_width=256,max_num_peaks=5,increase_noise=FALSE,seed=3)
  expect_true(all(identical(b1$X,array(expected_batch$X,c(2L,32L,1L)))))
  tt <- train_test_datasets(3,2,length=64,seed=5)
  expect_true(all(length(tt$train)==3L,length(tt$test)==2L,!tt$train$params$increase_noise,tt$test$params$increase_noise))
  
  # Derived spectra, regular map ordering, descending axes and one-point bands.
  s <- matrix(seq_len(24),nrow=4L,ncol=6L)
  comps <- radar_components(s,s/4,s/8,s/8,s/2,wavenumbers=6:1,meta=list(offset=1:4,scale=2:5),map_shape=c(2L,2L))
  expect_true(all(max(abs(comps$reconstruction-s))<1e-12,max(abs(comps$residual))<1e-12,
            identical(comps$corrected,s/2),identical(comps$denoised,s*7/8),length(comps)==4L))
  m <- as_map(comps)
  expect_true(all(identical(dim(m),c(2L,2L,6L)),identical(as.numeric(m[2,1,]),as.numeric(s[2,]/2)),
            identical(as_map(comps,"denoised"),array(s*7/8,c(2L,2L,6L)))))
  area <- band_image(comps,2,5)
  expected <- rowSums((s[,2:4]+s[,3:5])/4)
  expect_true(all(max(abs(as.numeric(area)-expected))<1e-12,all(band_image(comps,3,3)==0)))
  expect_error(as_map(comps,"typo")); expect_error(band_image(comps,2,5,"typo"))
  part <- comps[c(2,4)]
  expect_true(all(identical(part$meta$offset,c(2L,4L)),identical(part$meta$scale,c(3L,5L)),is.null(part$map_shape)))
  
  # Positioned maps use physical coordinates when file rows arrive x-fast.
  px <- matrix(seq_len(36),nrow=6,ncol=6)
  pos <- as.matrix(expand.grid(x=0:2,y=0:1))
  positioned <- radar_components(px,px/4,px/8,px/8,px/2,wavenumbers=1:6,
    positions=pos,position_names=c("x","y"),map_shape=c(2L,3L))
  expect_true(all(identical(as.numeric(as_map(positioned)[2,1,]),as.numeric(px[4,]/2)),
            identical(as.numeric(as_map(positioned)[1,2,]),as.numeric(px[2,]/2))))
  position_area <- band_image(positioned,1,6)
  raw_area <- rowSums((px[,-1,drop=FALSE]+px[,-6,drop=FALSE])/4)
  expect_true(all(position_area[2,1]==raw_area[4],position_area[1,2]==raw_area[2],
            isTRUE(all.equal(positioned$spectra,px))))
  
  # Implicit C-order maps also preserve arbitrary-dimensional geometry.
  cube_rows <- matrix(as.numeric(seq_len(72)),nrow=12,ncol=6)
  cmap <- radar_components(cube_rows,cube_rows/4,cube_rows/8,cube_rows/8,cube_rows/2,
    wavenumbers=1:6,map_shape=c(2L,3L,2L),meta=list(map_order="C"))
  expect_true(all(identical(as.numeric(as_map(cmap)[2,1,1,]),as.numeric(cube_rows[7,]/2)),
            identical(as.numeric(as_map(cmap)[1,2,2,]),as.numeric(cube_rows[4,]/2))))
  
  # Metrics follow the source's population summaries and pooled R squared.
  t <- matrix(c(0,1,2),nrow=1L); p <- matrix(c(.5,1,1),nrow=1L)
  expect_true(all(abs(rmse(t,p)-sqrt(1.25/3))<1e-12,max_error(t,p)==1,r2(t,t)==1,
            is.nan(r2(matrix(1,1,3),matrix(1,1,3))),snr(rep(1,10),rep(0,10))>1000))
  pred <- lapply(first$targets,function(z) z+.1)
  err <- evaluate_components(first,pred); stat <- error_statistics(first,pred)
  for(c in COMPONENTS) expect_true(all(abs(err[[c]]$rmse-.1)<1e-12,abs(err[[c]]$max_error-.1)<1e-12,
                               abs(stat[[c]]$rmse_median-.1)<1e-12,stat[[c]]$rmse_variance<1e-24))
  expect_error(rmse(1:3,1:4))
  
  # Plotting works without an optional graphics dependency.
  plot_file <- tempfile(fileext=".pdf")
  grDevices::pdf(plot_file)
  plot(comps,index=1)
  truth <- synthetic_dataset(s,s/4,s/8,s/8,s/2)
  plot_comparison(truth,comps,index=2)
  grDevices::dev.off()
  expect_true(all(file.info(plot_file)$size>0))
  unlink(plot_file)
})

