.libPaths("/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library")

library(ggplot2)
library(parallel)
library(gamlss)
library(tidyverse)
library(rcompanion)


fit_braincharts_model <- function(roi, dat, sw_dir, out_dir, plotting, ncores = 1){
  #' fit_braincharts_model(roi, data, sw_dir, out_dir)
  #'
  #' Fit a braincharts model to out-of sample data
  #' 
  #' @param roi The Desikan roi name to which the model should be applied
  #' @param NewData.all Df with out-of sample data and covariates in format derived from https://brainchart.shinyapps.io/brainchart/
  #' @param sw_dir Location of the Braincharts software and pretrained models
  #' @param out_dir Location where output should be saved
  #' @param plotting Whether histograms of centile scores should be plotted
  #' @param ncores How many cores will be used to compute model refit. Do not overload the computer!!
  
  
  print(paste("MODELLING: ", roi, "...................."))

  #------------------------- SLICE DATA --------------------------------
  # slice dat to include only current roi and meta information
  NewData <- dat[,1:10]
  NewData[[roi]] <- dat[[roi]]
  
  
  #------------------------- ADAPT MODEL TO NEW STUDY --------------------------------
  # load FIT and BOOT object per roi from pretrained braincharts model, distinguish model directory accordingly
  # if(grepl('CT_', roi, fixed = TRUE)) { 
  #   mdl_dir <- file.path(sw_dir, "model_new", "Model_CT_Refit")
  # } else if(grepl('SA_', roi, fixed=TRUE)) {
  #   mdl_dir <- file.path(sw_dir, "BrainCharts_refit_28022024")
  # } else {
  #   mdl_dir <- file.path(sw_dir, "model_new")
  # }
  mdl_dir <- file.path(sw_dir, "BrainCharts_refit_28022024")
  print(mdl_dir)
  tryCatch({
    FitObj <- readRDS( file=file.path(mdl_dir, paste("FIT_", roi, ".rds" , sep="")))
    BootObj <- readRDS( file=file.path(mdl_dir, paste("BOOT_", roi, ".rds" , sep="")))
  }, error = function(e){
    print(paste("FAILED TO LOAD BRAINCHARTS FILES", e$message))
    next
  })

  # calculate random effects for new study
  RESULT <- Calc.Novel( NewData, Fit=FitObj, BootFit=BootObj, Apply=TRUE, NumberCores = ncores )
  
  # save output
  outname <- file.path(out_dir, paste("result_gamlss_", roi, ".rds", sep=""))
  saveRDS(RESULT, file = outname)


  #------------------------- SUMMARIZE OUTPUT INTO A CSV FILE --------------------------------
  # summarize output into a df
  N <- nrow(NewData)
  NewDataDev <- data.frame(participant=character(N),
                            centile_score=double(N), #q.wre
                            bootu=double(N), #upper CI
                            bootl=double(N)) #lower CI

  # store relevant output in df
  NewDataDev$participant <- RESULT[["data"]][["participant"]]
  NewDataDev$age <- RESULT[["data"]][['age']]
  NewDataDev$age_days <- RESULT[["data"]][['age_days']]
  NewDataDev$sex <- RESULT[["data"]][['sex']]
  NewDataDev$study <- RESULT[["data"]][['study']]
  NewDataDev$fs_version <- RESULT[["data"]][['fs_version']]
  NewDataDev$country <- RESULT[["data"]][['country']]
  NewDataDev$run <- RESULT[["data"]][['run']]
  NewDataDev$session <- RESULT[["data"]][['session']]
  NewDataDev$dx <- RESULT[["data"]][['dx']]
  levels(NewDataDev$dx) <- c(levels(NewDataDev$dx), "preterm")
  NewDataDev$dx[is.na(NewDataDev$dx)] <- "preterm"

  centile_score_name <- paste(roi,"Transformed.q.wre", sep="") #centile score
  NewDataDev$centile_score <- RESULT[["data"]][[centile_score_name]]
  bootu_name <- paste(roi,"Transformed.q.wre.bootu", sep="")
  NewDataDev$bootu <- RESULT[["data"]][[bootu_name]]
  bootl_name <- paste(roi,"Transformed.q.wre.bootl",sep="")
  NewDataDev$bootl <- RESULT[["data"]][[bootl_name]]

  outnameCSV <- file.path(out_dir, paste("result_deviation_", roi, ".csv", sep=""))
  write.csv(NewDataDev, file=outnameCSV, row.names=FALSE)



  #------------------------- PLOT HISTOGRAMS OF RESULTING CENTILE SCORES PER GROUP --------------------------------
  if (plotting == TRUE) {
    print("PLOTTING HISTOGRAMS........")

    # calculate median per group
    med <- plyr::ddply(NewDataDev, "dx", summarise, grp.median=median(centile_score))

    plt <- ggplot(NewDataDev, aes(x=centile_score, color=dx)) +
      geom_histogram(fill="white", position="dodge",binwidth=0.05, alpha=0.5) +
      geom_vline(data=med, aes(xintercept=grp.median),color=c('red','blue'),linetype="dashed")+
      theme(legend.position="right")

    plot_dir <- file.path(out_dir, 'dev_score_distribution')
    dir.create(plot_dir, showWarnings = FALSE, recursive=TRUE)
    file_name=file.path(plot_dir, paste("dev_score_", roi, ".png",sep=""))
    png(file_name)
    print(plt)
    dev.off()
  }
    
  return(RESULT)
} 


plot_pop_curve_plotting <- function(dat, term_data_best, preterm_data_best, r, roi) {
  factor <- 10000
  ylabel_name <- substring(roi, 4)  # remove the prefix
  label_size <- 14
  
  ggplot(dat, aes(x = AgeNormal)) +
    geom_line(aes(y = factor * PRED.u975.pop), color = "azure4") +
    geom_line(aes(y = factor * PRED.l025.pop), color = "azure4") +
    geom_line(aes(y = factor * PRED.m500.pop), color = "black") +
      
    geom_point(data = term_data_best, aes(x = Age, y = factor * term_data_best[[r]]), pch = 19, cex = 0.5, color = "black") +
    geom_point(data = preterm_data_best, aes(x = Age, y = factor * preterm_data_best[[r]]), pch = 19, cex = 0.5, color = "red") +
      
    labs(x = "Age [years]", y = paste(ylabel_name, "[mm]", sep=" "))+
    theme_minimal() +
    theme(axis.title.y = element_text(size = label_size), axis.title.x = element_text(size = label_size),
          axis.text.x = element_text(size = label_size-2), axis.text.y = element_text(size = label_size-2),
          legend.position = "none")
}


plot_pop_curves <- function(roi, sw_dir, out_dir, best_analysis_dir) {
  print(paste("PLOTTING", roi, "........."))
  
  mdl_dir <- file.path(sw_dir, "BrainCharts_refit_28022024")
  
  FIT <- readRDS(file = file.path(mdl_dir, paste0("FIT_", roi, ".rds")))
  RESULT.BEST <- readRDS(file.path(best_analysis_dir, paste0("result_gamlss_", roi, ".rds")))
  
  # plot percentile curves with simulated data based on model parameters
  POP.CURVE.LIST <- list(AgeTransformed = seq(log(90), log(365 * 45), length.out = 2^12), sex = c("Female", "Male"))
  POP.CURVE.RAW <- do.call(expand.grid, POP.CURVE.LIST)
  CURVE <- Apply.Param(NEWData = POP.CURVE.RAW, FITParam = FIT$param)
  CURVE <- transform(CURVE, AgeNormal = (exp(AgeTransformed) - 280) / 365.245)
  CURVE.DF <- as.data.frame(CURVE)
  
  # look for ROI
  r <- paste0(roi, "Transformed.normalised")
  
  # Helper function for plotting actual dataset
  plot_and_save <- function(data_sex, roi, term_data_best, preterm_data_best, sex_label) {
    plot <- plot_pop_curve_plotting(
      dat = data_sex,
      term_data_best = term_data_best,
      preterm_data_best = preterm_data_best,
      r = r,
      roi = roi
    )
    file_name <- paste0("pop_curve_", roi, "_", sex_label, ".tiff")
    ggsave(filename = file_name, plot = plot, device = "tiff", path = out_dir, width = 7, height = 7, dpi = 500, bg = "white")
  }
  
  # Subset and plot data for females and males
  plot_and_save(
    data_sex = subset(CURVE.DF, sex == "Female"),
    roi = roi,
    term_data_best = subset(RESULT.BEST[["data"]], dx == "CN" & sex == "Female"),
    data_preterm_fem_best <- subset(RESULT.BEST[["data"]], (is.na(dx) | dx != "CN") & sex == "Female"),
    sex_label = "fem"
  )
  
  plot_and_save(
    data_sex = subset(CURVE.DF, sex == "Male"),
    roi = roi,
    term_data_best = subset(RESULT.BEST[["data"]], dx == "CN" & sex == "Male"),
    preterm_data_best = subset(RESULT.BEST[["data"]], (is.na(dx) | dx != "CN") & sex == "Male"),
    sex_label = "male"
  )
}