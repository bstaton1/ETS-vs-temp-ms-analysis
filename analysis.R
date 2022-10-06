
# clear the workspace
rm(list = ls(all = TRUE))

# for fitting gams
library(mgcv)

# number of bootstrap iters per fish
n_boot = 1000

# acquire the data
dat = read.csv("ETS-vs-temp-data.csv")

# what do they look like
head(dat)

# better species names:
Species_pretty = c(
  "black" = "Blacktail Shiner",
  "bluegill" = "Bluegill",
  "creek" = "Creek Chub",
  "rough" = "Rough Shiner",
  "sculpin" = "Banded Sculpin",
  "poosa" = "Tallapoosa Shiner",
  "stone" = "Largescale Stoneroller"
)
dat$Species = unname(Species_pretty[dat$Species])
dat$Species = as.factor(dat$Species)
dat = dat[order(dat$Species),]

# Drop these individuals: activity thermal profile looks highly questionable based on other individuals of the same species
dat = dat[-which(dat$Species == "Tallapoosa Shiner" & dat$Fish == 1),]
dat = dat[-which(dat$Species == "Rough Shiner" & dat$Fish == 10),]
dat = dat[-which(dat$Species == "Blacktail Shiner" & dat$Fish == 5),]

# quick count of records by species and fish
table(dat$Species, dat$Fish)

# CHOOSE WHETHER TO PERFORM ANOVA WITH INVERSE VARIANCE WEIGHTS
# on the individual T-values using bootstrap standard devs
do_weight = TRUE

# name of output directory
# contents not tracked by git if called "output"
out_dir = "output"

# if the output directory doesn't exist, create it
if (!dir.exists(out_dir)) dir.create(out_dir)

# load the functions
source("functions.R")

# decide on the file type for figure files
file_type = "pdf"

##### PART 1: FIT ALL GAMS AND PERFORM BOOTSTRAPPING #####

# container objects
all_Tval_ests = NULL  # stores summaries critical temperature values for each fish
all_pred_data = NULL  # stores bootstrap output of ETS curves

make_figure(file.path(out_dir, "boot-gams"), "pdf", h = 4, w = 6)
for (spp in levels(dat$Species)) {

  # print a progress indicator for species
  cat("\rSpecies: ", spp, "\n", sep = "")

  # subset only the data for this species
  dat_spp = subset(dat, Species == spp)

  for (fish in unique(dat_spp$Fish)) {

    # print a progress indicator for individual
    cat("\r", "Individual:", fish)

    # fit the gam for this fish
    fit = fit_gam(spp, fish)

    # produce bootstrap activity profiles for this fish
    pred_data = boot_gam(fit, n_boot)

    # obtain the critical Tvalues for each bootstrapped profile
    Tvals = do.call(rbind, lapply(0:n_boot, function(i) get_Tvals(pred_data, i)))

    # summarize the bootstrapped Tvalues
    Tval_ests = summarize_Tvals(Tvals)

    # make the plot for this fish
    plot_gam(pred_data, Tval_ests)

    # combine the Tval_ests with others
    all_Tval_ests = rbind(all_Tval_ests, Tval_ests)

    # combine the pred_data with others
    all_pred_data = rbind(all_pred_data, pred_data)
  }
}; cat("\n")
dev.off()

##### PART 2: PLOT ALL GAMS OVER TOP EACHOTHER #####

# function to make plot for one species
plot_gams = function(spp, side, ylim, letter) {

    # plot the raw data
    if (side == "left") {
      mar = c(2,2,0.25,0.4)
    } else {
      mar = c(2,0.4,0.25,2)
    }
    par(yaxs = "i", mar = mar)

    plot(ETS ~ Temp, data = subset(dat, Species == spp), cex = 1.4, xpd = TRUE, yaxt = "n", pch = 16, col = scales::alpha("grey25", 0.4), ylim = ylim)
    usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
    text(x = usr[1], y = usr[4] - ydiff * 0.075, labels = paste0("(", letter, ") ", spp), pos = 4, cex = 1.2)

    # lines for each fish's fit
    junk = lapply(unique(subset(all_pred_data, Species == spp)$Fish), function(x) {
      lines(boot_0 ~ Temp, data = subset(all_pred_data, Fish == x & Species == spp), col = "black")
    })

    # put the yaxis on the right side
    axis(side = ifelse(side == "left", 2, 4), las = 2)

}

make_figure(file.path(out_dir, "FIG1-ETS-curves"), file_type, h = 7, w = 5)

# which side should the yaxis go on for each panel
side = rep(c("left", "right"), each = 4)

# the order of the species
ord_spp = c("Largescale Stoneroller", "Banded Sculpin",
            "Blacktail Shiner", "Tallapoosa Shiner",
            "Bluegill", "Rough Shiner", "Creek Chub")

# graphical parameters
par(mfcol = c(4,2), mgp = c(2,0.35,0), tcl = -0.15, oma = c(1,1.5,0,0))

# loop over species, creating the plot for each one
for (s in 1:7) {
  plot_gams(spp = ord_spp[s], side = side[s],
            ylim = c(0,3.25),
            letter = letters[s])
}

# outer axes labels
mtext(side = 1, outer = TRUE, latex2exp::TeX("$Temperature\\,(\\degree C)$"))
mtext(side = 2, outer = TRUE, latex2exp::TeX("$Complex\\,III\\,Activity\\,(mL\\,O_2\\cdot\\,gWW^{-1}\\cdot\\,hr^{-1})$"), line = -0.25)

dev.off()

##### PART 3: FIT ANOVA MODELS COMPARING T-VALUES AMONG SPECIES #####

# function to fit the anova model, extract estimates, and perform multiple comparisons
fit_anova = function(type, do_weight = TRUE) {

  # subset the estimates to fit model to
  fit_data = subset(all_Tval_ests, Type == type)

  # fit model depending on whether data points are weighted by bootstrap SD
  if (do_weight) {
    fit = lm(temp_mean ~ Species, data = fit_data, weights = 1/temp_stdev^2)
  } else {
    fit = lm(temp_mean ~ Species, data = fit_data)
  }

  # produce predictions
  pred_data = data.frame(Species = unique(all_Tval_ests$Species))
  preds = predict(fit, newdata = pred_data, se.fit = TRUE)

  # get the fitted values and 95% CIs
  pred_data$estimate = preds$fit
  pred_data$lwr95ci = preds$fit + qt(0.025, preds$df) * preds$se.fit
  pred_data$upr95ci = preds$fit + qt(0.975, preds$df) * preds$se.fit
  pred_data$type = type
  pred_data = pred_data[order(pred_data$estimate),]

  # perform multiple comparisons and CLD (compact letter display)
  tukey = as.data.frame(emmeans:::pairs.emmGrid(emmeans::emmeans(fit, ~Species, type = "response")))
  p_vals = tukey$p.value; names(p_vals) = stringr::str_replace(tukey$contrast, " - ", "-")
  cld = multcompView::multcompLetters(p_vals, reversed = TRUE)$Letters
  cld = data.frame(Species = names(cld), cld = unname(cld))

  # combine clds with rest of predictions
  pred_data = merge(pred_data, cld, by = "Species")
  pred_data = pred_data[order(pred_data$estimate),]

  # return output
  list(pred_data = pred_data, fit = fit, type = type, tukey = tukey)

}

# function to plot the information for an anova
anova_plot = function(anova_out, xlim, label, x_axis) {

  # y-axis location of each species
  at_y = 1:7
  names(at_y) = as.character(anova_out$pred_data$Species)

  # blank plot with appropriate dimensions
  plot(at_y ~ anova_out$pred_data$estimate,
       xlim = xlim, ylim = range(at_y + 0.7, at_y - 0.2),
       type = "n", yaxt = "n", xlab = "", ylab = "", xaxt = "n")

  # draw y-axis
  axis(side = 2, at = at_y, labels = anova_out$pred_data$Species, las = 2, tcl = 0)

  # draw x-axis if requested
  if (x_axis) {
    axis(side = 1)
  }

  # add error bars and point estimates
  segments(anova_out$pred_data$lwr95ci, at_y, anova_out$pred_data$upr95ci, at_y)
  points(at_y ~ anova_out$pred_data$estimate, pch = 16, cex = 1.5)

  # draw cld letters
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  text(x = usr[2], y = at_y, labels = anova_out$pred_data$cld, pos = 2)

  # add a panel label
  text(x = usr[1] - xdiff * 0.015, y = usr[4] - ydiff * 0.065, labels = label, font = 2, cex = 1.1, pos = 4)

  # draw individual fish estimates (treated as data)
  fit_data = subset(all_Tval_ests, Type == anova_out$type)
  fit_data$at_y_i = at_y[as.character(fit_data$Species)]
  fit_data$at_y_ir = fit_data$at_y_i + runif(nrow(fit_data), -0.3, 0.3)
  with(fit_data, segments(temp_q25, at_y_ir, temp_q75, at_y_ir, col = scales::alpha("grey25", 0.25)))
  with(fit_data, points(x = temp_mean, y = at_y_ir, pch = 16, cex = 1.2, col = scales::alpha("grey25", 0.25)))
}

# create the plot
set.seed(5) # set random seed for replicable jitter of points

make_figure(file.path(out_dir, paste0("FIG2-summary-plot-wt-", do_weight)), file_type, h = 7, w = 3.5)
par(mfrow = c(4,1), mar = c(1,8.75,0.5,0.25), mgp = c(2,0.2,0), tcl = -0.15, oma = c(1.5,0,0,0), lend = "square", ljoin = "mitre")
anova_plot(fit_anova("upr", do_weight), c(18, 47.5), latex2exp::TeX("(a) T_{up}"), TRUE)
anova_plot(fit_anova("max", do_weight), c(18, 47.5), latex2exp::TeX("(b) T_{max}"), TRUE)
anova_plot(fit_anova("lwr", do_weight), c(18, 47.5), latex2exp::TeX("(c) T_{low}"), TRUE)
anova_plot(fit_anova("breadth", do_weight), c(2, 24), latex2exp::TeX("(d) T_{breadth}"), TRUE)
mtext(side = 1, latex2exp::TeX("$Temperature\\,(\\degree C)$"), line = 1.4, cex = 0.8)
dev.off()

# print the p-values and f-values for overall ANOVAs (reported in-text)
summary(fit_anova("upr", do_weight)$fit)
summary(fit_anova("max", do_weight)$fit)
summary(fit_anova("lwr", do_weight)$fit)
summary(fit_anova("breadth", do_weight)$fit)

##### PART 4: PLOT COMPARING T VALUES TO LITERATURE CT_MAX #####

# how far away from the single estimate should the endpoints
# of the bar be when there is only one CT_max value?
single_off = 0.05

# fit the anova models, extract estimates
lwr = fit_anova("lwr", do_weight)$pred_data[,c("Species", "estimate")]; colnames(lwr)[2] = "lwr"
max = fit_anova("max", do_weight)$pred_data[,c("Species", "estimate")]; colnames(max)[2] = "max"
upr = fit_anova("upr", do_weight)$pred_data[,c("Species", "estimate")]; colnames(upr)[2] = "upr"
pub = data.frame(
  Species = c("Rough Shiner",
              "Blacktail Shiner",
              "Bluegill",
              "Banded Sculpin",
              "Tallapoosa Shiner",
              "Largescale Stoneroller",
              "Creek Chub"
  ),
  pub_lwr = c(31.8,
              36.4 - single_off,
              31.5,
              29.4,
              34.0,
              35.8,
              35.7 - single_off),
  pub_upr = c(38.6,
              36.4 + single_off,
              37.5,
              30.9,
              36.2,
              38.0,
              35.7 + single_off)
)


out = merge(merge(lwr, merge(max, upr)), pub)
out = out[order(out$max),]

# y-axis location of each species
at_y = 1:7
names(at_y) = as.character(out$Species)

make_figure(file.path(out_dir, paste0("FIG3-all-T-ests-wt-", do_weight)), file_type, h = 3.5, w = 6)


par(mar = c(1,8.75,0.5,0.5), mgp = c(2,0.2,0), tcl = -0.15, oma = c(1.5,0,0,0), lend = "square", ljoin = "mitre")

# blank plot with appropriate dimensions
plot(at_y ~ out$max,
     xlim = c(20,45), ylim = range(at_y + 0.2, at_y - 0.2),
     type = "n", yaxt = "n", xlab = "", ylab = "", xaxt = "n")

# draw y-axis
axis(side = 2, at = at_y, labels = names(at_y), las = 2, tcl = 0)

# draw x-axis if requested
axis(side = 1)

segments(out$pub_lwr, at_y, out$pub_upr, at_y, lwd = 10, col = "grey")

points(at_y ~ out$lwr, pch = 18, col = "blue", cex = 1.5)
points(at_y ~ out$max, pch = 16, col = "black", cex = 1.5)
points(at_y ~ out$upr, pch = 17, col = "red", cex = 1.5)

legend("topleft", x.intersp = 0.5, pt.cex = 1.5, legend = c(latex2exp::TeX("T_{up}"), latex2exp::TeX("T_{opt}"), latex2exp::TeX("T_{low}"), latex2exp::TeX("CT_{max}")),
       pch = c(17, 16, 18, 15), col = c("red", "black", "blue", "grey"), bty = "n")

mtext(side = 1, latex2exp::TeX("$Temperature\\,(\\degree C)$"), line = 1.3, cex = 1)

dev.off()

##### TABLES #####

# fit all anova models
anovas = lapply(c("lwr", "max", "upr", "breadth"), function(t) fit_anova(type = t, do_weight = do_weight))

# export numerical output from anova analyses
write.csv(do.call(rbind, lapply(anovas, function(t) out = cbind(type = t$type, t$tukey))), file.path(out_dir, "tukey.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(anovas, function(t) out = cbind(type = t$type, t$pred_data))), file.path(out_dir, "anova_ests.csv"), row.names = FALSE)

# export individual-level critical temperature values
write.csv(all_Tval_ests, file.path(out_dir, "individuals.csv"), row.names = FALSE)
