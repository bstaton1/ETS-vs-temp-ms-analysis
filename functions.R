
##### FUNCTION TO CREATE A FIGURE FILE OF EITHER PDF OR PNG TYPE ####

make_figure = function(file, file_type, height, width) {

  file_type = tolower(file_type)

  if (file_type == "png") {
    f = function(file, height, width) png(file, height = height * 600, width = width * 600, res = 600)
  }

  if (file_type == "pdf") {
    f = function(file, height, width) pdf(file, height = height, width = width)
  }

  f(file = paste(file, file_type, sep = "."), height = height, width = width)
}

##### FUNCTION TO FIT A GAM TO AN INDIVIDUAL FISH OF A GIVEN SPECIES #####

fit_gam = function(spp, fish) {
  list(
    gam = gam(ETS ~ s(Temp), data = subset(dat, Species == spp & Fish == fish)),
    Species = spp,
    Fish = fish
  )
}

##### FUNCTION TO CREATE A PREDICTION DATA SET FOR ONE FISH

create_pred_data = function(fit) {
  expand.grid(Species = fit$Species,
              Fish = fit$Fish,
              Temp = seq(min(fit$gam$model$Temp), max(fit$gam$model$Temp), length = 100)
  )
}

##### FUNCTION TO GET PREDICTED ACTIVITY PROFILE, INCLUDING BOOTSTRAP ITERS #####

# bootstrap procedure:
#  1.) Fit original gam
#  2.) Extract coefficient estimates and their vcov matrix
#  3.) Draw n_boot replicates of coefficient estimates
#  4.) Produce the predicted curve for each set of coefficients

boot_gam = function(fit, n_boot) {

  # extract fitted coefs and their variance covariance matrix
  beta = coef(fit$gam)
  beta_vcov = vcov(fit$gam)

  # produce n_boot random replicate coefficient estimates
  beta_rand = mvtnorm::rmvnorm(n_boot, beta, beta_vcov)

  # create prediction data
  pred_data = create_pred_data(fit)

  # FIXME: figure out what this is called
  # extract the Linear Predictor Matrix (see ?predict.gam)
  Xp = predict(object = fit$gam, newdata = pred_data, type = "lpmatrix")

  # obtain predicted curve for original estimates
  pred_data$boot_0 = Xp %*% beta

  # obtain predicted curve for each bootstrap iteration
  boot_preds = apply(beta_rand, 1, function(betas) Xp %*% betas)

  # combine with rest of output
  colnames(boot_preds) = paste0("boot_", 1:n_boot)
  pred_data = cbind(pred_data, as.data.frame(boot_preds))

  # return the output
  return(pred_data)
}

##### FUNCTION TO CALCULATE THE CRITICAL TEMPERATURE VALUES OF INTEREST #####
# from a single activity profile (i.e., one bootstrap iteration)

get_Tvals = function(pred_data, iter) {

  # build the name of the bootstrap iter to calculate Tvals for
  boot_name = paste0("boot_", iter)

  # subset out only the needed info for this boot iter
  pred = pred_data[,c("Temp", boot_name)]; colnames(pred)[2] = "ETS"

  # calculate temperature at max ETS for this boot iter
  out = data.frame(Type = "max", Temp = pred$Temp[which.max(pred$ETS)], ETS = max(pred$ETS))

  # find the indices that result in ETS within 90% of the maximum
  indices_within90perc = which(pred$ETS > 0.9 * out$ETS)

  # get the temperatures associated with the bounds of the 90% of max ETS range
  out = rbind(out, cbind(Type = c("lwr", "upr"), pred[range(indices_within90perc),c("Temp", "ETS")]))

  # calculate the difference between upper and lower bound
  out = rbind(out, data.frame(Type = "breadth", Temp = diff(out$Temp[2:3]), ETS = NA))

  # format output
  rownames(out) = NULL
  out = cbind(Species = unique(pred_data$Species), Fish = unique(pred_data$Fish), iter = iter, out)

  # return output
  return(out)
}

##### FUNCTION TO SUMMARIZE THE CRITICAL TEMPERATURE VALUES OF INTEREST #####
# from bootstrapped activity profiles

summarize_Tvals = function(Tvals) {

  # extract the point estimates from the best fit curve
  pt_est = subset(Tvals, iter == 0)[,c("Species", "Fish", "Type", "ETS", "Temp")]; colnames(pt_est)[5] = "temp_mean"

  # calculate bootstrap summaries of uncertainty: 95% CI and SD
  q2.5 = aggregate(Temp ~ Species + Fish + Type, data = subset(Tvals, iter != 0), FUN = function(x) quantile(x, 0.025)); colnames(q2.5)[4] = "temp_q2.5"
  q10 = aggregate(Temp ~ Species + Fish + Type, data = subset(Tvals, iter != 0), FUN = function(x) quantile(x, 0.1)); colnames(q10)[4] = "temp_q10"
  q25 = aggregate(Temp ~ Species + Fish + Type, data = subset(Tvals, iter != 0), FUN = function(x) quantile(x, 0.25)); colnames(q25)[4] = "temp_q25"
  q50 = aggregate(Temp ~ Species + Fish + Type, data = subset(Tvals, iter != 0), FUN = function(x) quantile(x, 0.5)); colnames(q50)[4] = "temp_q50"
  q75 = aggregate(Temp ~ Species + Fish + Type, data = subset(Tvals, iter != 0), FUN = function(x) quantile(x, 0.75)); colnames(q75)[4] = "temp_q75"
  q90 = aggregate(Temp ~ Species + Fish + Type, data = subset(Tvals, iter != 0), FUN = function(x) quantile(x, 0.9)); colnames(q90)[4] = "temp_q90"
  q97.5 = aggregate(Temp ~ Species + Fish + Type, data = subset(Tvals, iter != 0), FUN = function(x) quantile(x, 0.975)); colnames(q97.5)[4] = "temp_q97.5"
  stdev = aggregate(Temp ~ Species + Fish + Type, data = subset(Tvals, iter != 0), FUN = function(x) sd(x)); colnames(stdev)[4] = "temp_stdev"

  # combine each summary with the others
  out_list = list(pt_est, q2.5, q10, q25, q50, q75, q90, q97.5, stdev)
  out = Reduce(function(x, y) merge(x, y, all = TRUE), out_list)

  # return the output
  return(out)
}

##### FUNCTION TO PLOT FITTED GAM OVER DATA WITH BOOTSTRAP CURVES AS WELL #####

plot_gam = function(pred_data, Tval_ests) {

  # set graphics parameters
  par(mar = c(3,3,2,1), mgp = c(2,0.35,0), tcl = -0.15, cex.axis = 0.8, lend = "square", ljoin = "mitre")

  # draw the bootstrapped profiles
  matplot(x = pred_data$Temp, y = pred_data[,paste0("boot_", 1:n_boot)], type = "l",
          col = scales::alpha("grey", 0.5), lty = 1, las = 1,
          xlab = "Temperature (C)", ylab = "ETS Activity",
          main = paste(unique(pred_data$Species), unique(pred_data$Fish), sep = "-"))

  # draw the critical Tvalue estimates with uncertainty
  add_Tval_ests(Tval_ests, "max")
  add_Tval_ests(Tval_ests, "lwr")
  add_Tval_ests(Tval_ests, "upr")
  add_Tval_ests(Tval_ests, "breadth")

  # draw the best fit profile
  lines(x = pred_data$Temp, y = pred_data$boot_0, col = "royalblue", lwd = 4)

  # draw the data points
  points(ETS ~ Temp,
         data = subset(dat, Species == levels(pred_data$Species) & Fish == unique(pred_data$Fish)),
         col = scales::alpha("royalblue", 0.5), pch = 16, cex = 1.5
  )
}

##### FUNCTION TO DRAW THE CRITICAL TEMPERATURE VALUES ON TOP OF PLOT #####

add_Tval_ests = function(Tval_ests, type) {

  # get plot coordinates
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])

  if (type == "breadth") {
    at_x = usr[1]
    at_y = usr[4] - ydiff * 0.05
    ests = round(Tval_ests[Tval_ests$Type == type,c("temp_mean", "temp_q2.5", "temp_q97.5")], 1)
    label = paste0("Breadth: ", ests[1], " (", ests[2], " - ", ests[3], ")")
    text(x = at_x, y = at_y, label, cex = 0.8, pos = 4)

  } else {
    # set the vertical offset -- so CIs don't overlap
    off = ifelse(type == "max", 0.05, ifelse(type == "lwr", 0.025, 0.075))

    # add the point estimate
    points(x = Tval_ests[Tval_ests$Type == type,"temp_mean"], y = usr[3] + ydiff * off,
           col = "salmon", pch = 16, cex = 1.5)

    # add a vertical segment from point estimate up to point on mean curve
    segments(Tval_ests[Tval_ests$Type == type,"temp_mean"], usr[3] + ydiff * off,
             Tval_ests[Tval_ests$Type == type,"temp_mean"], Tval_ests[Tval_ests$Type == type,"ETS"], col = "salmon", lty = 2)

    # add the 95% CIs for the estimate
    segments(Tval_ests[Tval_ests$Type == type,"temp_q2.5"], usr[3] + ydiff * off,
             Tval_ests[Tval_ests$Type == type,"temp_q97.5"], usr[3] + ydiff * off, col = "salmon")
  }
}
