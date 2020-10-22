##### Functions used for correlation plot #####

# Histogram function to show univariate distribution
hist_func <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram() +
    theme(panel.background = element_blank())
}

# Function to retrieve correlation
# and color according to spearman's rho
cor_func <- function(data, mapping, ...){
  
  # Retrieve data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr1 <- cor(x, y, method = 'spearman', use = 'complete.obs')
  corr2 <- cor(x, y, method = 'pearson', use = 'complete.obs')
  spear <- paste('\u03C1 =', as.character(round(corr1, 2)))
  pearsr <-  paste('r =' , as.character(round(corr2, 2)))
  
  colFn <- colorRampPalette(c("coral1", "white", "palegreen3"),
                            interpolate ='spline')
  fill <- colFn(100)[findInterval(corr1, seq(-1, 1, length=100))]
  
  ggally_text(
    label = paste(spear, pearsr, sep = '\n'), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...) + 
    theme(panel.background = element_rect(fill=fill)) 
}

# Function to plot bivariate distributions with fitted OLS regression
smooth_func <- function(data, mapping, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  ggplot(data = data, aes (x = x, y = y))+
    geom_point(shape = 16, colour = 'black', size = 0.5) +
    geom_smooth(method = 'lm', se = FALSE, colour = 'red') +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = 'black'))
}

# Function to combine the above functions into a correlation plot
corr_plot <- function(data, title, ...){
  
  ggpairs(data, columns = 1: ncol(data), title = title, switch = 'both',
          upper = list(continuous = wrap(cor_func)),
          lower = list(continuous = smooth_func),
          diag = list(continuous = hist_func),
          axisLabels = "none") + 
          theme(strip.background = element_rect(fill = "white"),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank())
}

##### Functions for percentile box-plot for bootstrapped coefficients #####
percentil_func <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}

# Check if the range between the 2.5th and 97.5th percentile ever crosses zero
range_func <- function(x){
  r <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
  if(min(r) < 0 && max(r) < 0){return(TRUE)}
  if(min(r) > 0 && max(r) > 0){return(TRUE)}
  else{return(FALSE)}}

# Function to plot coefficient and accuracy scores for bootstrapped logistic
# regression models. 
boxplot_func <- function(input_data, accuracy = TRUE, alt_acc = NULL, ...){
  bri <- mean(input_data[, names(input_data) == 'brier'])
  auc <- mean(input_data[, names(input_data) == 'auc'])
  
  # If alternative accuracy scores are provided, the accuracy scores from 
  # running with only island variables are subtracted from the accuracy of the
  # full model.
  if(!(is.null(alt_acc))){
    alt_bri <- mean(alt_acc[, names(alt_acc) == 'brier'])
    bri_score <- paste('Brier ', format(round(bri, 2), nsmall = 2), ' - ',
                       format(round(alt_bri, 2), nsmall = 2), ' = ',
                       format(round(bri, 2) - round(alt_bri, 2), nsmall = 2))
    alt_auc <- mean(alt_acc[, names(alt_acc) == 'auc'])
    auc_score <- paste('AUC ', format(round(auc, 2), nsmall = 2), ' - ',
                       format(round(alt_auc, 2), nsmall =2), ' = ', 
                       format(round(auc, 2) - round(alt_auc, 2), nsmall =2))
    score_lab <- paste(bri_score, auc_score, sep = '\n')
  }
  else{
    bril <- paste("Brier ", format(round(bri,2), nsmall = 2))
    aucl <- paste("AUC ", format(round(auc,2), nsmall = 2))
    score_lab <- paste(bril, aucl, sep = "\n")
  }
  
  # Exlude intercept, brier, auc and vif scores for plotting
  d <- melt(input_data[!(names(input_data) %in% c('(Intercept)', 'brier',
                      'auc', grep('vif_', names(input_data), value = TRUE)))])
  
  # Identify coefs above/below 0
  d$col <- 0
  a <- apply(input_data, 2, range_func)
  d$col[d$variable %in% names(a)[which(a == TRUE)]] <- 1
  con_f <- ifelse(unique(d$variable) %in% names(a)[which(a == TRUE)], 
                  'bold.italic', 'plain')
  
  # Call plot
  plt <- ggplot(d, aes(x = variable , y = value, fill = as.factor(col))) + 
    stat_summary(fun.data = percentil_func, geom = 'errorbar') +
    stat_summary(fun.data = percentil_func, geom = 'boxplot') +
    theme_bw() + geom_hline(yintercept = 0, cex = 0.5,
                            linetype = 'dashed') +
    xlab('') +
    ylab('') +
    scale_fill_manual(values = c('coral1','palegreen3')) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = con_f),
          legend.position = 'none', plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Accuracy scores are optional, and with two positional choices
  if(accuracy == TRUE){
  plt + annotate('text', -Inf, Inf, hjust = -0.1, vjust = 1.2, 
                 label = score_lab, size = 2.5)
  }
  else if(accuracy == 'bottom'){
    plt + annotate('text', -Inf, -Inf, hjust = -0.1, vjust = -0.2,
                   label = score_lab, size = 2.5)
  }
  else {
    plt
  }
}

##### Function to plot random forest results #####
rfplot_func <- function (input_data, accuracy = TRUE, impute = TRUE, 
                         alt_acc = NULL, exclude = NULL, x_digits = NULL){

  # The results with imputation are returned as list of lists of lists and 
  # require a bit of unpacking. 
  if(impute == TRUE){
    brier <- round(mean(melt(lapply(1:ncol(input_data[[1]]),
                    function(i) input_data[[i]][1,]))[,1]), 2)
    auc <- round(mean(melt(lapply(1:ncol(input_data[[1]]),
                    function(i) input_data[[i]][2,]))[,1]), 2)
    vars <- melt(lapply(1:ncol(input_data[[1]]),
                    function(i) input_data[[i]][3,]))
    varimp <- as.data.frame(lapply(1:length(unique(vars[,1])),
                    function(x) mean(vars[vars == 
                                as.character(unique(vars[x,1])), 3])))
    names(varimp) <- as.character(unique(vars[,1]))
  }
  # Without imputation the results are returned differently
  else{
    brier <- round(mean(melt(input_data[1,])[,1]), 2)
    auc <-    round(mean(melt(input_data[2,])[,1]), 2)
    varimp <- as.data.frame(lapply(1:length(unique(melt(input_data[3,])[,1])), 
                    function(x) mean(melt(input_data[3,])[melt(input_data[3,]) 
                    == as.character(unique(melt(input_data[3,])[x,1])),][,3]))) 
    names(varimp) <- as.character(unique(melt(input_data[3,])[,1]))
  }
  
  # If the accuracy scores are to have a separate accuracy score subtracted.
  # In this case the two island variables are passed in as the 'alt' accuracy
  # scores.
  if(!(is.null(alt_acc))){
    alt_brier <- mean(melt(alt_acc[1,])[,1])
    alt_auc <- mean(melt(alt_acc[2,])[,1])
    
    bri_score <- paste('Brier ', format(round(brier, 2), nsmall = 2), ' - ',
                       format(round(alt_brier, 2), nsmall = 2), ' = ',
                       format(round(brier, 2) - round(alt_brier,2), nsmall = 2))
    auc_score <- paste('AUC ', format(round(auc, 2), nsmall = 2), ' - ',
                       format(round(alt_auc, 2), nsmall =2), ' = ', 
                       format(round(auc, 2) - round(alt_auc, 2), nsmall =2))
    score_lab <- paste(bri_score, auc_score, sep = '\n')

    varimp <- varimp[,!(names(varimp) %in% exclude)]
  }
  else{
    # Establish label for accuracy scores
    b_lab <- paste("Brier ", format(round(brier, 2), nsmall = 2))
    a_lab <- paste("AUC ", format(round(auc, 2), nsmall = 2))
    score_lab <- paste(b_lab, a_lab, sep = "\n") 
  }
  
  # Call to plot
  plt <- ggplot(melt(varimp), aes(x = reorder(variable, value), y = value)) +
    geom_bar(stat = 'identity', fill = 'slategray1', color = 'black') +
    ylab('Variable importance') +
    xlab('') +
    theme_bw() +
    geom_hline(yintercept = 0) +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.title = element_text(size = 8),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0, unit = "pt")) +
    scale_y_continuous(labels = function(x) round(x, 2)) +
    coord_flip()
  
  # If number of digits on x-axis is specified (coord_flip() above switches x/y)
  if(!(is.null(x_digits))){
    plt <- plt + 
      scale_y_continuous(labels = number_format(accuracy = x_digits))
  }
  
  # Accuracy scores for the plot are optional 
  if(accuracy == TRUE){
    plt + annotate('text', -Inf, Inf, hjust = 1.1, vjust = -0.2,
                   label = score_lab, size = 2.5)
  }
  else{
    plt
  }
}

# Utility function to retrieve a legend for having one legend with
# multiple ggplots (used for the island histograms)
g_legend <- function(plt){
  tmp <- ggplot_gtable(ggplot_build(plt))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
