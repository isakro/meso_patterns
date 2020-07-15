library(ggplot2)
library(GGally)
library(gridExtra)
library(cowplot)
library(png)
library(grid)
library(reshape2)
library(scales)
library(car)
library(boot)
library(mice)
library(ModelMetrics)
library(caret)
library(randomForest)

# Load functions
source('functions/plot_functions.R')
source('functions/functions.R')

# For reproducibility
set.seed(123)

# Time script
start <- Sys.time()

# Preprocessing ---------------------------------

# Read in the raw data.
sitdat <- read.csv('../gis_output/site_data.csv')
sampdat <- read.csv('../gis_output/sample_data.csv')

# Set column 'class' to hold whether site/nonsite, and class of 
# sample (from hull or buffer constraint).
sitdat$class <- 'site'
sampdat$class[sampdat$sample_type == 'buffer'] <- 'buff'
sampdat$class[sampdat$sample_type == 'hull'] <- 'hull'

# Transform average aspect to deviation from south (note that the
# -n flag used with r.slope.aspect in GRASS returns north as
# 360 instead of the default 90)
sitdat$dev_south <- abs(sitdat$aspect_average - 180)  
sampdat$dev_south  <- abs(sampdat$aspect_average - 180) 

# Retrieve the columns of interest. Including the response but leaving that out
# of consideration for now. Elevation is only kept for imputation (see below). 
# As it was used to limit sampling, it is not used as a predictor.
icols <- c('class', 'phase', 'location', 'island_size','infiltration_mode', 
           'dev_south', 'viewsize_5k', 'viewsize_10k', 'avg_fetch', 
           'emergence_50', 'emergence_500', 'emergence_1k', 'elev_minimum') 
sites <- sitdat[, icols]
nonsites <- sampdat[, icols]

# Combine the data
dt <- rbind.data.frame(sites, nonsites)

# Rename most variables
names(dt)[names(dt) == 'elev_minimum'] <- 'elev' 
names(dt)[names(dt) == 'location'] <- 'loc'
names(dt)[names(dt) == 'island_size'] <- 'isl_si'
names(dt)[names(dt) == 'infiltration_mode'] <- 'infil'
names(dt)[names(dt) == 'viewsize_5k'] <- 'view_5'
names(dt)[names(dt) == 'viewsize_10k'] <- 'view_10'
names(dt)[names(dt) == 'avg_fetch'] <- 'fetch'
names(dt)[names(dt) == 'emergence_50'] <- 'emerg_50'
names(dt)[names(dt) == 'emergence_500'] <- 'emerg_500'
names(dt)[names(dt) == 'emergence_1k'] <- 'emerg_1k'

# Set infiltration class 5 to NA as this is 'not classified' by the 
# Geological Survey. Although ordinal, it is left as a numeric variable for
# correlation plots below. 
dt$infil[dt$infil == 5] <- NA

# Set island size to 0 in case of NA
dt$isl_si[is.na(dt$isl_si)] <- 0

# Make location a binary numerical variable where 0 = island, 1 = mainland
dt$loc <- as.numeric(dt$loc) - 1

# Inspect independent variables for correlation (cor |0.8| across all phases 
# is the defined threshold). Correlation among viewshed, fetch and shoreline
# emergence variables is problematic for the EM. In addition, some extremely
# skewed distributions. Colors follow spearmans rho.
corr_plot(dt[dt$phase == 'em', !(colnames(dt) %in% 
          c('class', 'phase', 'elev'))],
          'Independent variables - Early Mesolithic') 

# Same for Middle Mesolithic.
corr_plot(dt[dt$phase == 'mm', !(colnames(dt) %in% 
           c('class', 'phase', 'elev'))],
          'Independent variables - Middle Mesolithic')

# Correlations between short distance and longer distance emergence gone in 
# Late Mesolithic and fetch is just below the 0.8 threshold with viewshed 5km.
corr_plot(dt[dt$phase == 'lm', !(colnames(dt) %in% 
          c('class', 'phase', 'elev'))],
          'Independent variables - Late Mesolithic')

# Viewshed 10km has the strongest correlation with both fetch and viewshed 5km,
# and also intuitively seems like it is the variable resembling the other
# two the most. Viewshed was therefore retained at 5 km. Similar reasoning for
# shoreline emergence.  Using the 50 m and 1 km buffers seems most likely to 
# capture any long/short distance difference. Following Dorman et al. 2012, the
# variables are renamed to clarify that they are functioning as stand-ins for 
# several. 
dt$view <- dt$view_5
dt$emerg_shdist <- dt$emerg_50
dt$emerg_lgdist <- dt$emerg_1k
dt <- dt[, !(names(dt) %in% c('view_5', 'view_10', 'emerg_50',
                              'emerg_500', 'emerg_1k'))]

# Make infiltration level ordinal for imputation with random forest (see below,
# as well as in the functions constructed around random forest in functions.r).
# The variable is imputed as an ordinal variable to only replace missing with 
# observed values. The variable is then included in the models as numerical
# instead of categorical to not loose information on the order.
dt$infil <- factor(dt$infil, order = TRUE, levels = c(1, 2, 3, 4))

# Log transform the skewed variables. Having gone through alternatives, log
# transform seemed most reasonable. 0s are causing trouble for island size,
# so this is treated seperately below
cols <- c('view', 'fetch', 'emerg_shdist', 'emerg_lgdist')
dt_log <- log(dt[cols])
colnames(dt_log) <- paste('log', colnames(dt_log), sep = '_')
dt_log <- data.frame(dt[, !(names(dt) %in% cols)], dt_log) 

# As island size contains a lot of zero values, it was transformed by 
# log(isl_si), and then retaining the binary variable for whether
# or not the location is on an island (following Hosmer et al.
# 2013:106-107)
dt_log$isl_si <- log(dt_log$isl_si)
dt_log$isl_si[dt_log$isl_si == '-Inf'] <- 0
names(dt_log)[names(dt_log) == 'isl_si'] <- 'log_isl_si'

# A few randomly generated points are situated at elevations without
# any surrounding emergence values. These are given NA instead of 0.
# Another option would have been to retain these as zero, but would 
# caused issues for log transform as with the island size variable above.
dt_log$log_emerg_shdist[dt_log$log_emerg_shdist == '-Inf'] <- NA
dt_log$log_emerg_lgdist[dt_log$log_emerg_lgdist == '-Inf'] <- NA

#  Transform continous variables to take on values between 0 and 1
dt_log[,c('log_isl_si','dev_south', 'log_view', 'log_fetch',
          'log_emerg_shdist','log_emerg_lgdist')] <- 
            as.data.frame(apply(dt_log[,c('log_isl_si','dev_south', 'log_view',
            'log_fetch', 'log_emerg_shdist', 'log_emerg_lgdist')], 2,
            scale_func))

# Make infiltration numeric (ppm used for imputation with logistic regression
# below only 'donate' already observed values to missing values, instead of 
# using averages as rfImpute() does).
dt_log$infil <- as.numeric(as.character(dt_log$infil))

# Logistc regression - Sites/non-sites ---------------------------------

# Number of bootstrap samples
boot_n <- 9999

# Select Early Mesolithic data, samples from the hull constraints. 
# Make response factor.
em_hull_dat <- dt_log[dt_log$phase == 'em' & dt_log$class != 'buff',]
em_hull_dat$class <- as.factor(em_hull_dat$class)

# Call bootstrap of logistic regression. 
# logistic_reg() is defined in functions.r
em_hull_boot <- boot(data = em_hull_dat[, names(em_hull_dat) != 'phase'], 
                    dependent = 'class', extraPar = TRUE,
                    statistic = logistic_reg, R = boot_n)

# Retrieve data 
em_hull_mod <- as.data.frame(em_hull_boot$t)
names(em_hull_mod) <- names(em_hull_boot$t0)

# Call boxplot_func() defined in plot_functions.r which
# plots the results.
em_hull_log <- boxplot_func(em_hull_mod)

# Repeat for Early Mesolithic - Buffer sample
em_buff_dat <- dt_log[dt_log$phase == 'em' & dt_log$class != 'hull',]
em_buff_dat$class <- as.factor(em_buff_dat$class)

em_buff_boot <- boot(data = em_buff_dat[names(em_buff_dat) != 'phase'], 
                     dependent = 'class', extraPar = TRUE, 
                     statistic = logistic_reg, R = boot_n)

em_buff_mod <- as.data.frame(em_buff_boot$t)
names(em_buff_mod) <- names(em_buff_boot$t0)
em_buff_log <- boxplot_func(em_buff_mod)

# Middle Mesolithic - Hull sample
mm_hull_dat <- dt_log[dt_log$phase == 'mm' & dt_log$class != 'buff',]
mm_hull_dat$class <- as.factor(mm_hull_dat$class)

mm_hull_boot <- boot(data = mm_hull_dat[, names(mm_hull_dat) != 'phase'], 
                     extraPar = TRUE, dependent = 'class', 
                     statistic = logistic_reg, R = boot_n)

mm_hull_mod <- as.data.frame(mm_hull_boot$t)
names(mm_hull_mod) <- names(mm_hull_boot$t0)
mm_hull_log <- boxplot_func(mm_hull_mod)

# Middle Mesolithic - Buffer sample
mm_buff_dat <- dt_log[dt_log$phase == 'mm' & dt_log$class != 'hull',]
mm_buff_dat$class <- as.factor(mm_buff_dat$class)

mm_buff_boot <- boot(data = mm_buff_dat[, names(mm_buff_dat) != 'phase'], 
                     extraPar = TRUE, dependent = 'class',
                     statistic = logistic_reg, R = boot_n)

mm_buff_mod <- as.data.frame(mm_buff_boot$t)
names(mm_buff_mod) <- names(mm_buff_boot$t0)
mm_buff_log <- boxplot_func(mm_buff_mod)

# Late Mesolithic - Hull sample
lm_hull_dat <- dt_log[dt_log$phase == 'lm' & dt_log$class != 'hull',]
lm_hull_dat$class <- as.factor(lm_hull_dat$class)

lm_hull_boot <- boot(data = lm_hull_dat[, names(lm_hull_dat) != 'phase'], 
                    extraPar = TRUE, statistic = logistic_reg, 
                    dependent = 'class', R = boot_n)

lm_hull_mod <- as.data.frame(lm_hull_boot$t)
names(lm_hull_mod) <- names(lm_hull_boot$t0)
lm_hull_log <- boxplot_func(lm_hull_mod)

# Late Mesolithic - Buffer sample.
lm_buff_dat <- dt_log[dt_log$phase == 'lm' & dt_log$class != 'hull',]
lm_buff_dat$class <- as.factor(lm_buff_dat$class)

lm_buff_boot <- boot(data = lm_buff_dat[,!(names(lm_buff_dat) %in% c('phase',
                     'log_emerg_500', 'log_emerg_1k'))], extraPar = TRUE, 
                     dependent = 'class', statistic = logistic_reg,
                     R = boot_n)

lm_buff_mod <- as.data.frame(lm_buff_boot$t)
names(lm_buff_mod) <- names(lm_buff_boot$t0)
lm_buff_log <- boxplot_func(lm_buff_mod)

# Random forest - Sites/non-sites ---------------------------------

# Number of trees is set to a large arbitrary number
tree_n <- 5000

# Random forest. Early Mesolithic - Hull sample
rfd_em_hull <- dt[dt$phase == 'em' & dt$class != 'buff', 
                  names(dt) != 'phase']
rfd_em_hull$class <- as.factor(rfd_em_hull$class)

# Nested cross-validation of random forest for five imputed datasets defined
# in functions.r
em_hull_rf <- imp_nest_rf(rfd_em_hull, nimp = 5, ntrees = tree_n,
                          depend = 'class', nindep = ncol(rfd_em_hull) - 2)

# Pass results to plot function defined in plot_functions.r
em_hull_rfp <- rfplot_func(em_hull_rf)

# Random forest. Early Mesolithic - Buffer sample
rfd_em_buff <- dt[dt$phase == 'em' & dt$class != 'hull', 
                  names(dt) != 'phase']
rfd_em_buff$class <- as.factor(rfd_em_buff$class)

em_buff_rf <- imp_nest_rf(rfd_em_buff, nimp = 5, ntrees = tree_n, 
                          depend = 'class', nindep = ncol(rfd_em_buff) - 2)

em_buff_rfp <- rfplot_func(em_buff_rf)

# Random forest. Middle Mesolithic - Hull sample
rfd_mm_hull <- dt[dt$phase == 'mm' & dt$class != 'buff', 
                  names(dt) != 'phase']
rfd_mm_hull$class <- as.factor(rfd_mm_hull$class)

mm_hull_rf <- imp_nest_rf(rfd_mm_hull, nimp = 5, ntrees = tree_n, 
                          depend = 'class', nindep = ncol(rfd_mm_hull) - 2)

mm_hull_rfp <- rfplot_func(mm_hull_rf)

# Random forest. Middle Mesolithic - Buffer sample
rfd_mm_buff <- dt[dt$phase == 'mm' & dt$class != 'hull', 
                  names(dt) != 'phase']
rfd_mm_buff$class <- as.factor(rfd_mm_buff$class)

mm_buff_rf <- imp_nest_rf(rfd_mm_buff, nimp = 5, ntrees = tree_n, 
                          depend = 'class', nindep = ncol(rfd_mm_buff) - 2)

mm_buff_rfp <- rfplot_func(mm_buff_rf)

# Random forest. Late Mesolithic - Hull sample
rfd_lm_hull <- dt[dt$phase == 'lm' & dt$class != 'buff', 
                  names(dt) != 'phase']
rfd_lm_hull$class <- as.factor(rfd_lm_hull$class)

lm_hull_rf <- imp_nest_rf(rfd_lm_hull, nimp = 5, ntrees = tree_n, 
                          depend = 'class', nindep = ncol(rfd_lm_hull) - 2)

lm_hull_rfp <- rfplot_func(lm_hull_rf)

# Random forest. Late Mesolithic - Buffer sample
rfd_lm_buff <- dt[dt$phase == 'lm' & dt$class != 'hull', 
                  names(dt) != 'phase']
rfd_lm_buff$class <- as.factor(rfd_lm_buff$class)

lm_buff_rf <-imp_nest_rf(rfd_lm_buff, nimp = 5, ntrees = tree_n, 
                         depend = 'class', nindep = ncol(rfd_lm_buff) - 2)

lm_buff_rfp <- rfplot_func(lm_buff_rf)

# Plot and save site/non-site results --------------------------------- 
em_hull <- readPNG('../gis_output/maps/em_hull.png')
plot_row <- plot_grid(rasterGrob(em_hull), em_hull_log, em_hull_rfp, 
                     labels = 'AUTO', hjust = c(-3.7, -0.5, -0.5), ncol = 3,
                     rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/em_hull_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

em_buff <- readPNG('../gis_output/maps/em_buff.png')
plot_row <- plot_grid(rasterGrob(em_buff), em_buff_log, em_buff_rfp, 
                      labels = 'AUTO', hjust = c(-3.7, -0.5, -0.5), ncol = 3,
                      rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/em_buff_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

mm_hull <- readPNG('../gis_output/maps/mm_hull.png')
plot_row <- plot_grid(rasterGrob(mm_hull), mm_hull_log, mm_hull_rfp, 
                      labels = 'AUTO', hjust = c(-3.7, -0.5, -0.5), ncol = 3,
                      rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/mm_hull_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

mm_buff <- readPNG('../gis_output/maps/mm_buff.png')
plot_row <- plot_grid(rasterGrob(mm_buff), mm_buff_log, mm_buff_rfp, 
                      labels = 'AUTO', hjust = c(-3.7, -0.5, -0.5), ncol = 3,
                      rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/mm_buff_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

lm_hull <- readPNG('../gis_output/maps/lm_hull.png')
plot_row <- plot_grid(rasterGrob(lm_hull), lm_hull_log, lm_hull_rfp, 
                      labels = 'AUTO', hjust = c(-3.7, -0.5, -0.5), ncol = 3,
                      rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/lm_hull_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

lm_buff <- readPNG('../gis_output/maps/lm_buff.png')
plot_row <- plot_grid(rasterGrob(lm_buff), lm_buff_log, lm_buff_rfp, 
                      labels = 'AUTO', hjust = c(-3.7, -0.5, -0.5), ncol = 3,
                      rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/lm_buff_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

# Logistic Regression - Comparison across phases -------------------------------

# First checking the frequencies of samples on islands in the different phases.
# 0 = island, 1 = mainland
isls_h <- rbind.data.frame(dt[dt$phase == 'em' & dt$class == 'hull',],
                           dt[dt$phase == 'mm' & dt$class == 'hull',],
                           dt[dt$phase == 'lm' & dt$class == 'hull',])
isls_h$phase <- factor(isls_h$phase, levels(isls_h$phase)[c(1,3,2)])

# Plot
hull_isl <- ggplot(isls_h, aes(x = phase, fill = factor(loc, 
                          labels = c('Island', 'Mainland')))) +
  geom_bar(position = 'dodge', alpha = 0.5, colour = 'black') + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  xlab('') +
  ggtitle('Hull samples') +
  scale_fill_manual('legend', 
                    values = c('Island' = '#E69F00', 'Mainland' = '#56B4E9')) + 
  theme(legend.position = 'bottom', legend.spacing.x = unit(0.2, 'cm'))

# Retrieve data - Buffer sample
isls_b <- rbind.data.frame(dt[dt$phase == 'em' & dt$class == 'buff',],
                           dt[dt$phase == 'mm' & dt$class == 'buff',],
                           dt[dt$phase == 'lm' & dt$class == 'buff',])
isls_b$phase <- factor(isls_b$phase, levels(isls_b$phase)[c(1,3,2)])

# Plot
buff_isl <- ggplot(isls_b, aes(x = phase, fill = factor(loc,
                          labels = c('Island', 'Mainland')))) +
  geom_bar(position = 'dodge', alpha = 0.5, colour = 'black') + 
  theme_classic() +
  theme(legend.title = element_blank()) +
  xlab('') +
  ggtitle('Buffer samples') +
  scale_fill_manual('legend', values = c('Island' = '#E69F00',
                                         'Mainland' = '#56B4E9')) 

# Retireve legend to have one for both plots
leg <- g_legend(hull_isl)

island_hist <- grid.arrange(arrangeGrob(hull_isl + 
                  theme(legend.position = 'none'),
                  buff_isl + theme(legend.position = 'none'), nrow = 1),
                            leg, nrow = 2, heights = c(10,1))
ggsave('../figures/island_hist.png', island_hist, width = 10, height = 7,
       units = 'cm', dpi = 600)

# Retrieve sites per phase with corresponding hull samples
em_s <- dt_log[dt_log$phase == 'em' & dt_log$class == 'site',]
em_h <- dt_log[dt_log$phase == 'em' & dt_log$class == 'hull',]
mm_s <- dt_log[dt_log$phase == 'mm' & dt_log$class == 'site',]
mm_h <- dt_log[dt_log$phase == 'mm' & dt_log$class == 'hull',]
lm_s <- dt_log[dt_log$phase == 'lm' & dt_log$class == 'site',]
lm_h <- dt_log[dt_log$phase == 'lm' & dt_log$class == 'hull',]

# Columns for which to find percentile rank compared to hull sample.
# (Note that this is also done for the ordinal infiltration variable). 
pcrank_cols <- c('log_isl_si','infil', 'dev_south', 'log_fetch', 'log_view',
                 'log_emerg_shdist', 'log_emerg_lgdist')

# Find percentile rank for each of the variables, and
# reinsert loc, elevation and phase.
em_p <- as.data.frame(prcr_func(em_s, em_h, pcrank_cols))
em_p$loc <- em_s$loc
em_p$elev <- em_s$elev
em_p$phase <- 'em'

mm_p <- as.data.frame(prcr_func(mm_s, mm_h, pcrank_cols))
mm_p$loc <- mm_s$loc
mm_p$elev <- mm_s$elev
mm_p$phase <- 'mm'

lm_p <- as.data.frame(prcr_func(lm_s, lm_h, pcrank_cols))
lm_p$loc <- lm_s$loc
lm_p$elev <- lm_s$elev
lm_p$phase <- 'lm'

# Combine the datasets
p_dat <- rbind(em_p, mm_p, lm_p)

# Reorder
p_dat <- p_dat[,c("phase", "loc", "log_isl_si",
                  "infil", "dev_south", "elev", "log_fetch",
                  "log_view", "log_emerg_shdist", "log_emerg_lgdist")] 

# Compare Early Mesolithic with Middle Mesolithic
em_mm <- p_dat[p_dat$phase != 'lm',]
em_mm$phase <- factor(em_mm$phase, levels = c('mm', 'em'))

# First run using only island variables
em_mm_islboot <- boot(data = em_mm, extraPar = TRUE, statistic = logistic_reg,
                  dependent = 'phase', independent = c('loc', 'log_isl_si'),
                  R = boot_n)

em_mm_isl <- as.data.frame(em_mm_islboot$t)
names(em_mm_isl) <- names(em_mm_islboot$t0)
em_mm_islp <- boxplot_func(em_mm_isl)

# Second run to get difference in performance. 
em_mm_boot <- boot(data = em_mm, extraPar = TRUE, statistic = logistic_reg,
                   dependent = 'phase', R = boot_n)

em_mm_mod <- as.data.frame(em_mm_boot$t)
names(em_mm_mod) <- names(em_mm_boot$t0)
em_mm_log <- boxplot_func(em_mm_mod[, !(names(em_mm_mod) %in%
              c('loc', 'log_isl_si'))], alt_acc = em_mm_isl, 
              accuracy = 'bottom')


# Compare Early Mesolithic with Late Mesolithic
em_lm <- p_dat[p_dat$phase != 'mm',]
em_lm$phase <- factor(em_lm$phase, levels = c('lm', 'em'))

em_lm_islboot <- boot(data = em_lm, extraPar = TRUE, statistic = logistic_reg,
                      dependent = 'phase',
                      independent = c('loc', 'log_isl_si'), R = boot_n)
em_lm_isl <- as.data.frame(em_lm_islboot$t)
names(em_lm_isl) <- names(em_lm_islboot$t0)
em_lm_islp <- boxplot_func(em_lm_isl)

em_lm_boot <- boot(data = em_lm, extraPar = TRUE, statistic = logistic_reg,
                   dependent = 'phase', R = boot_n)
em_lm_mod <- as.data.frame(em_lm_boot$t)
names(em_lm_mod) <- names(em_lm_boot$t0)
em_lm_log <- boxplot_func(em_lm_mod[,!(names(em_lm_mod) %in%
              c('loc', 'log_isl_si'))], alt_acc = em_lm_isl,
              accuracy = 'bottom')

# Compare Middle Mesolithic with Late Mesolithic
mm_lm <- p_dat[p_dat$phase != 'em',]
mm_lm$phase <- as.factor(mm_lm$phase)

mm_lm_islboot <- boot(data = mm_lm, extraPar = TRUE, statistic = logistic_reg,
                      dependent = 'phase',
                      independent = c('loc', 'log_isl_si'), R = boot_n)

mm_lm_isl <- as.data.frame(mm_lm_islboot$t)
names(mm_lm_isl) <- names(mm_lm_islboot$t0)
mm_lm_islp <- boxplot_func(mm_lm_isl)

mm_lm_boot <- boot(data = mm_lm, extraPar = TRUE, statistic = logistic_reg,
                   dependent = 'phase', R = boot_n)

mm_lm_mod <- as.data.frame(mm_lm_boot$t)
names(mm_lm_mod) <- names(mm_lm_boot$t0)
mm_lm_log <- boxplot_func(mm_lm_mod[, !(names(mm_lm_mod) %in%
                              c('loc', 'log_isl_si'))], alt_acc = mm_lm_isl)

# Random forest - Comparison accross phases ---------------------------------

# Retrieve sites per phase with corresponding hull samples
emrf_s <- dt[dt$phase == 'em' & dt$class == 'site',]
emrf_h <- dt[dt$phase == 'em' & dt$class == 'hull',]
mmrf_s <- dt[dt$phase == 'mm' & dt$class == 'site',]
mmrf_h <- dt[dt$phase == 'mm' & dt$class == 'hull',]
lmrf_s <- dt[dt$phase == 'lm' & dt$class == 'site',]
lmrf_h <- dt[dt$phase == 'lm' & dt$class == 'hull',]

# Columns for which to find percentile rank (everything except loc)
pcrank_cols <- c('infil', 'isl_si', 'dev_south', 'fetch',
                 'view', 'emerg_shdist', 'emerg_lgdist')

# Find percentile rank for each of the variables, and
# reinsert loc, elevation and phase.
emrf_p <- as.data.frame(prcr_func(emrf_s, emrf_h, pcrank_cols))
emrf_p$loc <- emrf_s$loc
emrf_p$elev <- emrf_s$elev
emrf_p$phase <- 'em'

mmrf_p <- as.data.frame(prcr_func(mmrf_s, mmrf_h, pcrank_cols))
mmrf_p$loc <- mmrf_s$loc
mmrf_p$elev <- mmrf_s$elev
mmrf_p$phase <- 'mm'

lmrf_p <- as.data.frame(prcr_func(lmrf_s, lmrf_h, pcrank_cols))
lmrf_p$loc <- lmrf_s$loc
lmrf_p$elev <- lmrf_s$elev
lmrf_p$phase <- 'lm'

# Combine the datasets
prf_dat <- rbind(emrf_p, mmrf_p, lmrf_p)

# Random forest. Early Mesolithic - Middle Mesolithic
emrf_mm <- prf_dat[prf_dat$phase != 'lm',]
emrf_mm$phase <- factor(emrf_mm$phase, levels = c('mm', 'em'))

em_mm_rfisl <- imp_nest_rf(emrf_mm, indep = c('loc', 'isl_si'),
                           ntrees = tree_n, depend = 'phase', nindep = 2)
rfplot_func(em_mm_rfisl, impute = FALSE)

emrf_mm_rf <- imp_nest_rf(emrf_mm, nimp = 5, ntree = tree_n, depend = 'phase',
                          nindep = ncol(emrf_mm)-3)
em_mm_rfp <- rfplot_func(emrf_mm_rf, alt_acc = em_mm_rfisl,
                         exclude = c('loc', 'isl_si'), x_digits = 0.001)

# Random forest. Early Mesolithic - Late Mesolithic
emrf_lm <- prf_dat[prf_dat$phase != 'mm',]
emrf_lm$phase <- factor(emrf_lm$phase, levels = c('lm', 'em'))

em_lm_rfisl <- imp_nest_rf(emrf_lm, indep = c('loc', 'isl_si'),
                           ntrees = tree_n, depend = 'phase', nindep = 2)
rfplot_func(em_lm_rfisl, impute = FALSE)

emrf_lm_rf <- imp_nest_rf(emrf_lm, nimp = 5, ntree = tree_n, depend = 'phase',
                          nindep = ncol(emrf_lm)-3)
em_lm_rfp <- rfplot_func(emrf_lm_rf, alt_acc = em_lm_rfisl,
                         exclude = c('loc', 'isl_si'), x_digits = 0.001)

# Random forest. Middle Mesolithic - Late Mesolithic
mmrf_lm <- prf_dat[prf_dat$phase != 'em',]
mmrf_lm$phase <- factor(mmrf_lm$phase, levels = c('lm', 'mm'))

mm_lm_rfisl <- imp_nest_rf(mmrf_lm, indep = c('loc', 'isl_si'),
                           ntrees = tree_n, depend = 'phase', nindep = 2)
rfplot_func(mm_lm_rfisl, impute = FALSE)

mmrf_lm_rf <- imp_nest_rf(mmrf_lm, nimp = 5, ntree = tree_n, depend = 'phase',
                          nindep = ncol(mmrf_lm)-3)
mm_lm_rfp <- rfplot_func(mmrf_lm_rf, alt_acc = mm_lm_rfisl, 
                         exclude = c('loc', 'isl_si'), x_digits = 0.001)

# Plot results of comparison between phases --------------------------------- 

em_lm <- readPNG('../gis_output/maps/em_lm.png')
plot_grid(rasterGrob(em_lm), em_lm_log, em_lm_rfp,
          labels = 'AUTO', hjust = c(-3.7, -0.5, -0.5), ncol = 3,
          rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/em_lm_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

em_mm <- readPNG('../gis_output/maps/em_mm.png')
plot_grid(rasterGrob(em_mm), em_mm_log, em_mm_rfp,
          labels = 'AUTO', hjust = c(-3.7, -0.5, -0.5), ncol = 3,
          rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/em_mm_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

mm_lm <- readPNG('../gis_output/maps/mm_lm.png')
plot_grid(rasterGrob(mm_lm), mm_lm_log, mm_lm_rfp,
                  labels = 'AUTO',  hjust = c(-3.7, -0.5, -0.5), ncol = 3,
          rel_widths = c(4/10, 3/10, 3/10))
ggsave('../figures/mm_lm_fig.png', width = 30, height = 10,
       units = 'cm', dpi = 600)

# Find elapsed time
end <- Sys.time() - start
print(end)
