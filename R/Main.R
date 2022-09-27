library(devtools)
install_github("josephlewis/leastcostpath@dev")
library(leastcostpath)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)

################
################
#### SET-UP ####
################
################

set.seed(1)
nsims <- 1000
neigh <- 4

source("./R/Functions.R")

#############################################
#############################################
#### PREPROCESSS DIGITAL ELEVATION MODEL ####
#############################################
#############################################

dem <- terra::rast("./Data/FS_240.tif")
dem <- normalise_raster(dem)
dem <- dem * 10

####################################
####################################
#### CREATE COST FUNCTION TABLE ####
####################################
####################################

cfs <- c("tobler", "tobler offpath", "modified tobler", "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female", 'rees', "davey", 'garmy', 'kondo-saino', 'naismith', 'campbell', "campbell 2019", "herzog", "llobera-sluckin", 'minetti', "wheeled transport")

index_table <- expand.grid(1:length(cfs), 1:length(cfs))

##############################
##############################
#### CREATE COST SURFACES ####
##############################
##############################

cs_list <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem, cost_function = x, neighbours = neigh, crit_slope = 12, percentile = 0.5)
  return(cs)
  }
)

###########################################################################
###########################################################################
#### CALCULATE AND COMPARE LCPs CREATED USING DIFFERENT COST FUNCTIONS ####
###########################################################################
###########################################################################

lcps <- list()

for (i in 1:nsims) {
  
  print(i)
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(dem)), 2, type = "random")
  
  compare_lcp <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp <- do.call(rbind, compare_lcp)
  compare_lcp$sim_no <- i
  
  lcps[[i]] <- compare_lcp  

}

lcps2 <- do.call(rbind, lcps)
lcps3 <- sf::st_drop_geometry(lcps2)

#FS = fractal surface, 240 = 2.40 fractal dimension, 4 = neighbours, 1000 = nsims 
write.csv(lcps3, file = "./Output/Data/FS_240_4_1000.csv", row.names = FALSE)

############################
#### METHODS PLOT ONE ####
############################

cf_df <- list()

for (i in 1:length(cfs)) {
  
  print(i)
  print(cfs[i])
  
  time_cfs <- c("tobler", "tobler offpath", "rees", "davey", "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female", "modified tobler", "garmy", "kondo-saino", "naismith", "campbell", "campbell 2019")
  
  energy_cfs <- c("herzog", "llobera-sluckin", "minetti")
  
  if(cfs[i] %in% time_cfs) { 
    type <- "Time-based"
  } else if (cfs[i] %in% energy_cfs) { 
    type <- "Energy-based"  
  } else { 
    type <- "Wheel-based"  
  }
  
  cf <- cost(cost_function = cfs[i], crit_slope = 12, percentile = 0.5)
  
  slope <- seq(-0.4, 0.4, 0.01)
  vals <- 1/cf(slope)
  vals <- (vals - min(vals)) / (max(vals) - min(vals))
  
  cf_df[[i]] <- data.frame(slope = slope, cf_vals = vals, cf = cfs[i], type = type)
  
}

cf_df <- do.call(rbind, cf_df)
cf_df$type <- factor(cf_df$type, levels = c("Time-based", "Energy-based", "Wheel-based"))
cf_df$cf <- factor(cf_df$cf, levels = cfs)

methods_plot1 <- ggplot(cf_df) + 
  geom_vline(xintercept = 0, lty = 2, colour = "grey80", size = 0.2) +
  geom_line(aes(x = slope, y = cf_vals, group = cf, colour = type)) +
  facet_wrap(~cf, nrow = 6) + 
  scale_colour_manual(values = c("#003f5c", "#bc5090", "#ffa600")) + 
  labs(x = "Mathematical Slope", y = "Relative cost", colour = "Theory type") + 
  theme_classic() + 
  theme(legend.position = "bottom")

ggplot2::ggsave(filename = "methods_plot_01.png", plot = methods_plot1, path = "./Output/Figures/", dpi = 300, width = 6, height = 6)

###################
###################
#### ANALYSIS #####
###################
###################

############################
##### RESULTS PLOT ONE #####
############################

lcps4 <- lcps3[lcps3$normalised_pdi == 0,]

lcp_df1 <- as.data.frame(table(lcps4$cf1_type, lcps4$cf2_type))

lcp_df1 <- lcp_df1 %>%
  group_by(Var1) %>%
  mutate(perc = Freq /sum(Freq) * 100, Var2 = Var2) %>%
  arrange(Var1)

lcp_df1$Var1 <- factor(lcp_df1$Var1, levels = rev(c("Time-based", "Energy-based", "Wheel-based")))
lcp_df1$Var2 <- factor(lcp_df1$Var2, levels = c("Time-based", "Energy-based", "Wheel-based"))

results_plot1 <- ggplot(lcp_df1, aes(x = Var1, y = perc, fill = Var2)) +
  geom_bar(position="fill", stat="identity", show.legend = TRUE) +
  labs(x = "Route Type", y = "Percentage (%)", fill = "Theory Type") +
  geom_text(aes(label = paste0(round(perc, 2))),
            position = position_fill(vjust = 0.5), size = 4) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(expand = expansion(mult = c(0.3, 0))) + 
  scale_fill_manual(values = c("#003f5c", "#bc5090", "#ffa600")) + 
  coord_flip() + 
  theme_classic() + 
  theme(axis.line.y =element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom") + 
  theme(legend.justification = "right", legend.margin = margin(t = 2, r = 4, b = 2, l = 2, unit = "mm"))

ggplot2::ggsave(plot = results_plot1, filename = "./Output/Figures/results_plot_01.png", width = 10, height = 4, dpi = 300)

############################
##### RESULTS PLOT TWO #####
############################

lcp_df2 <- as.data.frame(table(lcps4$cf1_type, lcps4$cf2))

lcp_df2 <- lcp_df2 %>%
  group_by(Var1) %>%
  mutate(perc = Freq /sum(Freq) * 100, Var2 = Var2) %>%
  filter(Var1 == "Wheel-based") %>%
  arrange(Freq)

colnames(lcp_df2) <- c("cf1_type", "cf2", "Freq", "perc")
lcp_df2 <- left_join(x = lcp_df2, y = unique(lcps4[c("cf2", "cf2_type")]), by = "cf2")
lcp_df2$cf2_type <- factor(lcp_df2$cf2_type, levels = c("Time-based", "Energy-based", "Wheel-based"))
lcp_df2$Var2 <- factor(lcp_df2$cf2, levels = as.character(unique(lcp_df2$cf2)))

results_plot2 <- ggplot(lcp_df2, aes(x = Var2, y = perc, fill = cf2_type)) +
  geom_bar(stat="identity", show.legend = TRUE) +
  labs(x = "Cost Function", y = "Percentage (%)", fill = "Theory Type") +
  geom_text(aes(label = paste0(round(perc, 2))), hjust = "left", size = 4, nudge_y = 2) + 
  scale_fill_manual(values = c("#003f5c", "#bc5090", "#ffa600")) + 
  lims(y = c(0, 100)) + 
  coord_flip() + 
  theme_classic() + 
  theme(legend.position = "bottom") + 
  theme(legend.title = element_text( size=8), legend.text=element_text(size=8))

ggplot2::ggsave(plot = results_plot2, filename = "./Output/Figures/results_plot_02.png", width = 6, height = 4, dpi = 300)