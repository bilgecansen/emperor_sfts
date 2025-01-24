
library(foreach)
library(tidyverse)
library(patchwork)
library(MCMCvis)
library(rstan)
library(doSNOW)
library(ggrepel)
library(rstanarm)
library(sf)
library(ggspatial)


# Install with devtools::install_github("hrbrmstr/ggalt")
library(ggalt)

empe_sites <- read.csv("data/empe_sites_new.csv")
sites_update <- read.csv("data/colony_attributes_update.csv")

data_pop <- readRDS("data/data_pop_empe.rds")
dat_fled <- read.csv("data/count_fecundity_updated.csv")
dat_sara <- read.csv("data/Breeding.csv")

theme_set(theme_bw())


# Load environmental data -------------------------------------------------

data_esm <- readRDS("data/data_env_empe.rds") %>%
  filter(site_id %in% unique(data_pop$sat$site_id)) %>%
  filter(year >= 2009 & year <= 2018) %>%
  select(site_id, year, contains(c("aice"))) %>%
  arrange(site_id)

data_faice <- readRDS("data/data_faice.rds")[[6]] %>%
  filter(site_id %in% unique(data_pop$sat$site_id)) %>%
  filter(year >= 2009 & year <= 2018) %>%
  select(site_id, year, contains(c("rearing"))) %>%
  arrange(site_id)

data_fdice <- readRDS("data/data_fdice.rds") %>%
  filter(site_id %in% unique(data_pop$sat$site_id)) %>%
  filter(year >= 2009 & year <= 2018) %>%
  select(site_id, year, contains(c("rearing"))) %>%
  arrange(site_id)

data_ftice <- readRDS("data/data_ftice.rds")[[6]] %>%
  filter(site_id %in% unique(data_pop$sat$site_id)) %>%
  filter(year >= 2009 & year <= 2018) %>%
  arrange(site_id)

data_env_avg1 <- data_esm %>%
  group_by(site_id) %>%
  summarise(across(-year, mean)) %>%
  ungroup()

data_env_avg2 <- left_join(data_faice, data_fdice, 
                           by = c("site_id", "year")) %>%
  left_join(data_ftice, by = c("site_id", "year")) %>%
  group_by(site_id) %>%
  summarise(across(-year, mean)) %>%
  ungroup()

data_env_avg <-left_join(data_env_avg1, data_env_avg2, by = "site_id")

env_mat <- data_env_avg[,-1] %>%
  as.matrix()

env_mat_std <- apply(env_mat, 2, function(x) (x - mean(x))/sd(x))


# Fit horshoe regression --------------------------------------------------

# Output from LaRue et al. (2024)
N_data <- readRDS("data/N_data.rds")

## check site order is the same for abundance and env datasets
all.equal(N_data$site_id, data_env_avg$site_id)

dat_finn <- list(y = N_data$mean,
                 y_sd = N_data$sd,
                 X = env_mat_std,
                 N = nrow(env_mat_std),
                 M = ncol(env_mat_std),
                 scale = 1)

res_finn <- stan(file = 'lm_finn.stan', 
     data = dat_finn,
     iter = 3000,
     control = list(adapt_delta = 0.9999, 
                    max_treedepth = 20))

# Plot slopes
beta <-  abs(MCMCsummary(res_finn, params = "beta")[,1])
idx_beta <- order(beta, decreasing = T)
var_names <- c("SIC (Laying)", "SIC (Incubation)", "SIC (Rearing)", 
               "SIC (Nonbreed)","FIA (Rearing)", "NOW (Rearing)", "TOE", "TOB")
  
theme_set(theme_bw())
g_imp1 <- ggplot() +
  geom_col(aes(x = factor(var_names[idx_beta], 
                          levels = var_names[idx_beta][8:1]), 
               y = beta[idx_beta]), 
           fill = "#FF800E", alpha = 0.8, color = "black") +
  #geom_point(mapping = aes(x = factor(var_names[idx_beta], 
  #                                    levels = var_names[idx_beta][8:1]), 
  #                         y = beta[idx_beta]), col = "darkorange", size = 2) +
  #geom_segment(aes(y = 0, 
  #                 x = var_names[idx_beta], 
  #                 yend = beta[idx_beta], 
  #                 xend = var_names[idx_beta]), 
  #             color = "darkorange", linewidth = 1, alpha = 0.8) +
  labs(y = "Slope Estimate (absolute)", title = "Colony Abundance") +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 9),
        axis.title.x = element_blank(),
        panel.border = element_blank(), 
        #panel.grid.major = element_blank()) +
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10)) + 
        #axis.line.y = element_line(colour = "black"),
        #axis.ticks.x = element_blank()) + 
  coord_flip()


# Fit linear regression ---------------------------------------------------

x1 <- env_mat_std[,1]
x2 <- env_mat_std[,6]

dat_lm1 <- list(y = N_data$mean,
                y_sd = N_data$sd,
                X1 = x1,
                X_pred1 = seq(min(x1), max(x1), length.out = 100),
                X2 = x2,
                X_pred2 = seq(min(x2), max(x2), length.out = 100),
                N = 50,
                L = 100,
                S_pred = (dat_sara$climate - 
                          mean(dat_sara$climate))/sd(dat_sara$climate),
                S = length(dat_sara$climate))

dat_lm2 <- list(y = N_data$mean,
                y_sd = N_data$sd,
                X = x2,
                N = 50,
                L = 100,
                X_pred = seq(min(x2), max(x2), length.out = 100),
                S_pred = (dat_sara$climate - 
                            mean(dat_sara$climate))/sd(dat_sara$climate),
                S = length(dat_sara$climate))

res_lm_multi <- stan(file = 'lm_multi.stan', 
                     data = dat_lm1,
                     iter = 3000)

res_lm_faice <- stan(file = 'lm.stan', 
                     data = dat_lm2,
                     iter = 10000,
                     thin = 2)

# Model fit
MCMCsummary(res_lm_multi, params = "Rsq")

# Univariate plots
y1 <- MCMCsummary(res_lm_multi, params = "mu_pred1")
y1_mean <- y1[,1]
y1_min <- y1[,3]
y1_max <- y1[,5]

y2 <- MCMCsummary(res_lm_faice, params = "mu_pred")
y2_mean <- y2[,1]
y2_min <- y2[,3]
y2_max <- y2[,5]

x_pred1 <- seq(min(x1), max(x1), length.out = 100)*sd(env_mat[,1]) + 
  mean(env_mat[,1])

g1 <- ggplot() +
  geom_ribbon(aes(x = x_pred1, ymin = y1_min, ymax = y1_max), alpha = 0.8, 
              fill = "grey") +
  geom_segment(aes(x = env_mat[,1], xend = env_mat[,1], y = N_data$min, 
                   yend = N_data$max)) +
  geom_line(aes(x = x_pred1, y = y1_mean), col = "#FF800E", linewidth = 1.5, 
            linetype = 2) +
  geom_point(aes(x = env_mat[,1], y = N_data$mean), alpha= 0.8, 
             color = "#FF800E", size = 3) +
  geom_point(aes(x = env_mat[,1], y = N_data$mean), shape = 1, size = 3, 
             stroke = 1.1) +
  labs(y = "Colony Abundance (log)", x = "SIC (Laying)") +
  theme(#axis.title.y = element_blank(),
        panel.border = element_blank(), 
        #panel.grid.major = element_blank()) +
        panel.grid.minor = element_blank())

x_pred2 <- seq(min(x2), max(x2), length.out = 100)*sd(env_mat[,6]) +
  mean(env_mat[,6])

g2 <- ggplot() +
  geom_ribbon(aes(x = x_pred2, ymin = y2_min, ymax = y2_max), alpha = 0.8, 
              fill = "grey") +
  geom_segment(aes(x = env_mat[,6], xend = env_mat[,6], y = N_data$min, 
                   yend = N_data$max)) +
  geom_line(aes(x = x_pred2, y = y2_mean), col = "#FF800E", size = 1.5, 
            linetype = 2) +
  geom_point(aes(x = env_mat[,6], y = N_data$mean), alpha= 0.8, 
             color = "#FF800E", size = 2) +
  geom_point(aes(x = env_mat[,6], y = N_data$mean), shape = 1, size = 2, 
             stroke = 1.1) +
  labs(y = "Colony Abundance (log)", x = "NOW (Rearing)") +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        plot.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.minor = element_blank())

g1
ggsave("fig_S1.jpeg", width = 6, height = 5, units = "in") 


# Maps --------------------------------------------------------------------

empe_coord <- data_pop$sat %>%
  distinct(site_id, img_long, img_lat)
empe_coord <- empe_coord[order(empe_coord$site_id),]
empe_coord <- empe_coord[-c(9, 34, 35, 36, 46),]

#world <- map_data("world")

z <- N_data$mean
z2 <- c()
for (i in 1:length(z)) {
  if (z[i] < 5) z2[i] <- "4"
  else
    if (z[i] < 6) z2[i] <- "5"
    else
      if (z[i] < 7) z2[i] <- "6"
      else
        if (z[i] < 8) z2[i] <- "7"
        else
          if (z[i] < 9) z2[i] <- "8"
          else
            if (z[i] < 10) z2[i] <- "9"
            else
              if (z[i] < 11) z2[i] <- "10"
}
z2 <- factor(z2, levels = 4:10)

antarctic <- st_as_sf(maps::map('world', plot = FALSE, fill = TRUE))
idx <- which(antarctic$ID == "Antarctica")
antarctic <- antarctic[idx,]
antarctic <- st_transform(antarctic, 3031)

col_dat <- st_multipoint(as.matrix(empe_coord[,2:3])) %>%
  st_sfc()
st_crs(col_dat) <- 4326
col_dat <- st_transform(col_dat, 3031)

z <- foreach(h = 1:nrow(empe_coord)) %do% {
  st_point(as.matrix(empe_coord[h,2:3]))
}
z_sfc <- st_sfc(z)
st_crs(z_sfc) <- 4326
z_sfc <- st_transform(z_sfc, 3031)

d <- st_sf(data.frame(color = pred_fled, 
                      size = z2, 
                      geom = z_sfc))

m1 <- ggplot() +
  geom_sf(data = antarctic, alpha = 0.5, fill="lightblue") +
  geom_sf(data = d, alpha = 0.9, size = 0.1) +
  theme(#panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) 

m2 <- ggplot() +
  geom_sf(data = antarctic, alpha = 0.5, fill="lightblue") +
  geom_sf(data = d[34,], alpha = 0.9, size = 1) +
  theme(#panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA)) 

ggplot() +
  geom_sf(data = antarctic, alpha = 0.5, fill="lightblue") +
  geom_sf(data = d[34,], alpha = 0.9, size = 4) +
  theme(#panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA))

ggsave("map1.pdf", width = 2, height = 2, units = "in") 


# Comparison with Labrouse et al. predictions -----------------------------

# Environmental data
data_esm_pgeo <- readRDS("data/data_env_empe.rds") %>%
  filter(site_id %in% "PGEO") %>%
  filter(year >= 1979 & year <= 2017) %>%
  select(site_id, year, contains(c("aice"))) %>%
  arrange(site_id)

data_esm_pgeo$NOW <- dat_sara$climate

env_mat_pgeo <- data_esm_pgeo[,-c(1:2)] %>%
  as.matrix()

env_mat_pgeo_std <- apply(env_mat_pgeo, 2, function(x) (x - mean(x))/sd(x))

# Horseshoe regression
dat_finn_pgeo <- list(y = dat_sara$yvar,
                      X = env_mat_pgeo_std,
                      N = nrow(env_mat_pgeo_std),
                      M = ncol(env_mat_pgeo_std),
                      scale = 1)

res_finn_pgeo <- stan(file = 'lm_finn_alt.stan', 
                      data = dat_finn_pgeo,
                      iter = 3000,
                      control = list(adapt_delta = 0.9999, 
                                     max_treedepth = 20))

beta_pgeo <-  abs(MCMCsummary(res_finn_pgeo, params = "beta")[,1])
idx_beta_pgeo <- order(beta_pgeo, decreasing = T)
var_names_pgeo <- c("SIC (Laying)", "SIC (Incubation)", "SIC (Rearing)", 
                    "SIC (Nonbreed)",
               "NOW (Rearing)")
g_imp2 <- ggplot() +
  geom_col(aes(x = factor(var_names_pgeo[idx_beta_pgeo], 
                          levels = var_names_pgeo[idx_beta_pgeo][5:1]), 
               y = beta_pgeo[idx_beta_pgeo]), 
           fill = "#006BA4", alpha = 0.8, color = "black") +
  #geom_point(aes(x = factor(var_names_pgeo[idx_beta_pgeo], 
  #                          levels = var_names_pgeo[idx_beta_pgeo][5:1]), 
  #               y = beta_pgeo[idx_beta_pgeo]), col = "darkorange", size = 2) +
  #geom_segment(aes(y = 0, 
  #                 x = var_names_pgeo[idx_beta_pgeo], 
  #                 yend = beta_pgeo[idx_beta_pgeo], 
  #                 xend = var_names_pgeo[idx_beta_pgeo]), 
  #             color = "darkorange", size = 1, alpha = 0.8) +
  labs(y = "Slope Estimate (absolute)", title = "Breeding Success (PGEO)") +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 9),
        axis.title.x = element_blank(),
        panel.border = element_blank(), 
        #panel.grid.major = element_blank()) +
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10)) + 
  #axis.line.y = element_line(colour = "black"),
  #axis.ticks.x = element_blank()) + 
  coord_flip()

# Univariate Plots
res_sara <- stan_lm(yvar ~ climate, data = dat_sara, prior = NULL)
y_sara <- predict(res_sara, newdata = data.frame(climate = dat_sara$climate))
cor(y_sara, dat_sara$yvar)^2

x_pred3 <- seq(min(dat_sara$climate), max(dat_sara$climate), length.out = 100)
y3 <- predict(res_sara, newdata = data.frame(climate = x_pred3), se.fit = T)
y3_mean <- y3$fit
y3_min <- y3$fit - 1.96*y3$se.fit
y3_max <- y3$fit + 1.96*y3$se.fit

g3 <- ggplot() +
  geom_ribbon(aes(x = x_pred3*1000, ymin = y3_min, ymax = y3_max), 
              alpha = 0.8, fill = "grey") +
  geom_line(aes(x = x_pred3*1000, y = y3_mean), col = "#006BA4", 
            size = 1.5, linetype = 2) +
  geom_point(aes(x = dat_sara$climate*1000, y = dat_sara$yvar), alpha= 0.8, 
             color = "#006BA4", size = 2) +
  geom_point(aes(x = dat_sara$climate*1000, y = dat_sara$yvar), shape = 1, 
             size = 2, stroke = 1.1) +
  labs(y = "Breeding Success", x = "NOW (Rearing)") +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        plot.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.minor = element_blank())

# Time for Space substitution
pred_sara <- predict(res_sara, 
                     newdata = data.frame(climate = data_env_avg$fdice_rearing/1000))

cor(pred_sara, N_data$mean)

#g4 <- ggplot() +
  #geom_point(aes(x = pred_sara, y = log(N_mean)), alpha= 0.8, 
             #color = "darkorange", size = 3) +
  #geom_point(aes(x = pred_sara, y = log(N_mean)), shape = 1, size = 3, 
             #stroke = 1.1) +
  #geom_smooth(aes(x = pred_sara, y = log(N_mean)), method = "lm", se = F, 
              #col = "black", size = 1.5, linetype = 2) +
  #labs(y = "Colony Abundance (log)", x = "Breeding Success Prediction") +
  #theme(panel.border = element_blank(), 
        #panel.grid.minor = element_blank())

#m3 <- ggplot() + 
  #geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               #color = "black", fill = "white") + 
  #geom_path(data = world, aes(x = long, y = lat, group = group)) + 
  #scale_y_continuous(breaks = (-2:2) * 30) + 
  #scale_x_continuous(breaks = (-4:4) * 45) + 
  #theme(axis.title.x = element_blank(), 
        #axis.text.x = element_blank(), 
        #axis.ticks = element_blank(), 
        #axis.title.y = element_blank(), 
        #axis.text.y = element_blank(),
        #panel.border = element_blank()) +
  #coord_map("ortho", orientation = c(-90, 0, 0), ylim = c(-90, -60)) +
  #geom_point(data = empe_coord, 
             #aes(x = img_long, y = img_lat, col = pred_sara, size = z2),
             #alpha = 0.9) +
  #geom_point(data = empe_coord, 
             #aes(x = img_long, y = img_lat, size = z2), 
             #shape = 1, 
             #stroke = 1.1,
             #alpha = 0.9) +
  #scale_size_manual(values = c("4" = 1, 
                               #"5" = 2, 
                               #"6" = 3,
                               #"7" = 5,
                               #"8" = 7,
                               #"9" = 10,
                               #"10" = 13)) +
  #geom_label_repel(data = empe_coord, 
                   #aes(x = img_long, y = img_lat, 
                       #label = site_id), 
                   #size = 2,
                   #nudge_x = 7,
                   #nudge_y = 7) +
  #scale_color_gradientn(colors = terrain.colors(10, rev = T)) +
  #labs(size = "Abundance (log)", col = "B.S. Predictions")

#(g4 | m3) + 
  #plot_annotation(tag_levels = "a")
#ggsave("fig_space.pdf", width = 12, height = 6, units = "in") 


# Comparison with fledgling abundance predictions ------------------------

# Horseshoe regression
dat_finn_fled <- list(y = log(dat_fled[28:66,3]),
                      X = env_mat_pgeo_std,
                      N = nrow(env_mat_pgeo_std),
                      M = ncol(env_mat_pgeo_std),
                      scale = 1)

res_finn_fled <- stan(file = 'lm_finn_alt.stan', 
                      data = dat_finn_fled,
                      iter = 3000,
                      control = list(adapt_delta = 0.9999, 
                                    max_treedepth = 20))

beta_fled <-  abs(MCMCsummary(res_finn_fled, params = "beta")[,1])
idx_beta_fled <- order(beta_fled, decreasing = T)
var_names_fled <- c("SIC (Laying)", "SIC (Incubation)", "SIC (Rearing)", 
                    "SIC (Nonbreed)", "NOW (Rearing)")
g_imp3 <- ggplot() +
  #geom_point(aes(x = factor(var_names_fled[idx_beta_fled], 
                            #levels = var_names_fled[idx_beta_fled][5:1]), 
                 #y = beta_fled[idx_beta_fled]), 
             #col = "darkorange", size = 2) +
  #geom_segment(aes(y = 0, 
                   #x = var_names_fled[idx_beta_fled], 
                   #yend = beta_fled[idx_beta_fled], 
                   #xend = var_names_fled[idx_beta_fled]), 
               #color = "darkorange", size = 1, alpha = 0.8) +
  geom_col(aes(x = factor(var_names_fled[idx_beta_fled], 
                          levels = var_names_fled[idx_beta_fled][5:1]), 
               y = beta_fled[idx_beta_fled]), 
           fill = "#006BA4", alpha = 0.8, color = "black") +
  labs(y = "Slope Estimate (absolute)", title = "Fledgling Abundance (PGEO)") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text = element_text(size = 9),
        panel.border = element_blank(), 
        #panel.grid.major = element_blank()) +
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10)) + 
  #axis.line.y = element_line(colour = "black"),
  #axis.ticks.x = element_blank()) + 
  coord_flip()

(g_imp1 /  #inset_element(m1, left = 0.7, bottom = 0, right = 1, top = 0.8,
                          #align_to = "full")) / 
g_imp2 /  #inset_element(m2, left = 0.7, bottom = 0.05, right = 1, top = 0.8,
                          #align_to = "full")) / 
g_imp3) +  #inset_element(m2, left = 0.7, bottom = 0.05, right = 1, top = 0.8,
                          #align_to = "full"))) +
  plot_annotation(tag_levels = list(c("a", "b", "c", ""))) +
  plot_layout(axes = "collect_x",
              heights = c(4,3,3))

ggsave("fig2.pdf", width = 90, height = 180, units = "mm", dpi = 600)
ggsave("fig2.jpeg", width = 90, height = 180, units = "mm", dpi = 600)

# Univariate plots
dat_lm_fled <- data.frame(z = log(dat_fled[28:66,3]),
                          climate = dat_sara$climate)
res_fled <- stan_lm(z ~ climate, data = dat_lm_fled, prior = NULL)
y_fled <- predict(res_fled, newdata = data.frame(climate = dat_sara$climate))
cor(y_fled, dat_lm_fled$z)

x_pred4 <- seq(min(dat_lm_fled$climate), max(dat_lm_fled$climate), 
               length.out = 100)
y4 <- predict(res_fled, newdata = data.frame(climate = x_pred4), se.fit = T)
y4_mean <- y4$fit
y4_min <- y4$fit - 1.96*y4$se.fit
y4_max <- y4$fit + 1.96*y4$se.fit

g6 <- ggplot() +
  geom_ribbon(aes(x = x_pred4*1000, ymin = y4_min, ymax = y4_max), 
              alpha = 0.8, fill = "grey") +
  geom_line(aes(x = x_pred4*1000, y = y4_mean), col = "#006BA4", 
            size = 1.5, linetype = 2) +
  geom_point(aes(x = dat_lm_fled$climate*1000, y = dat_lm_fled$z), 
             alpha= 0.8, color = "#006BA4", size = 2) +
  geom_point(aes(x = dat_lm_fled$climate*1000, y = dat_lm_fled$z), 
             shape = 1, size = 2, stroke = 1.1) +
  labs(y = "Fledling Abundance (log)", x = "NOW (Rearing)") +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        plot.title = element_text(size = 10),
        panel.border = element_blank(), 
        panel.grid.minor = element_blank())

(g2 / g3 / g6) +
  plot_annotation(tag_levels = "a") +
  plot_layout(axes = "collect_x",
              heights = c(4,3,3))
ggsave("fig3.pdf", width = 90, height = 180, units = "mm", dpi = 600)
ggsave("fig3.jpeg", width = 90, height = 180, units = "mm", dpi = 600)

# Time for Space substitution
pred_fled <- predict(res_fled, 
                     newdata = data.frame(climate = 
                                            data_env_avg$fdice_rearing/1000))

cor(pred_fled, N_data$mean)

g7 <- ggplot() +
  geom_point(aes(x = pred_fled, y = N_data$mean, col = pred_fled, size = z2), 
             alpha = 0.8) +
  geom_point(aes(x = pred_fled, y =  N_data$mean, size = z2), shape = 1, 
             stroke = 1.1) +
  geom_smooth(aes(x = pred_fled, y =  N_data$mean), method = "lm", se = F, 
              col = "black", size = 1.5, linetype = 2) +
  scale_size_manual(values = c("4" = 2, 
                               "5" = 2, 
                               "6" = 3,
                               "7" = 4,
                               "8" = 6,
                               "9" = 8,
                               "10" = 10)) +
  scale_colour_gradient2_tableau(
    palette = "Orange-Blue Diverging") +
  labs(y = "Colony Abundance (log)", 
       x = "Fledgling Abundance Prediction (log)") +
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

antarctic <- st_as_sf(maps::map('world', plot = FALSE, fill = TRUE))
idx <- which(antarctic$ID == "Antarctica")
antarctic <- antarctic[idx,]
antarctic <- st_transform(antarctic, 3031)

col_dat <- st_multipoint(as.matrix(empe_coord[,2:3])) %>%
  st_sfc()
st_crs(col_dat) <- 4326
col_dat <- st_transform(col_dat, 3031)

z <- foreach(h = 1:nrow(empe_coord)) %do% {
  st_point(as.matrix(empe_coord[h,2:3]))
}
z_sfc <- st_sfc(z)
st_crs(z_sfc) <- 4326
z_sfc <- st_transform(z_sfc, 3031)

d <- st_sf(data.frame(color = pred_fled, 
                      size = z2, 
                      geom = z_sfc))

m4 <- ggplot() +
  geom_sf(data = antarctic, alpha = 0.5) +
  annotation_scale() +
  geom_sf(data = d, aes(col = color, size = size), alpha = 0.9) +
  geom_sf(data = d, aes(size = size), alpha = 0.9, shape = 1, 
          stroke = 1.1) +
  scale_size_manual(values = c("4" = 1, 
                               "5" = 2, 
                               "6" = 3,
                               "7" = 5,
                               "8" = 7,
                               "9" = 10,
                               "10" = 13)) +
  scale_colour_gradient2_tableau(
    palette = "Orange-Blue Diverging") +
  labs(size = "Colony Abundance (log)", 
       col = "Fledling Abundance\nPrediction (log)") +
  coord_sf(label_axes = list(bottom = "E", left = "E", 
                             top = "E", right = "E")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

g7 + m4 +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 20))
ggsave("fig5.pdf", width = 360, height = 240, units = "mm", dpi = 600)

ggplot() +
  geom_sf(data = antarctic, alpha = 0.5, fill="lightblue") +
  geom_sf(data = d, aes(size = size), alpha = 0.6, 
          stroke = 1.1) +
  scale_size_manual(values = c("4" = 1, 
                               "5" = 2, 
                               "6" = 3,
                               "7" = 5,
                               "8" = 7,
                               "9" = 10,
                               "10" = 13)) +
  theme(#panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none")
ggsave("map2.pdf", width = 6, height = 6, units = "in") 


# Space for time substitution ---------------------------------------------

pred_time <- MCMCsummary(res_lm_faice, params = "mu_s_pred")[,1]

cor(pred_time, dat_sara$yvar)
cor(pred_time, log(dat_fled[28:66,3]))
cor(pred_time, log(dat_fled[28:66,2]))

ab1 <- (pred_time - mean(pred_time))/sd(pred_time)
ab2 <- (dat_sara$yvar- mean(dat_sara$yvar))/sd(dat_sara$yvar)
ab3 <- (log(dat_fled[28:66,3]) - mean(log(dat_fled[28:66,3])))/
  sd(log(dat_fled[28:66,3]))
#ab4 <- (log(dat_fled[28:66,2]) - mean(log(dat_fled[28:66,2])))/
  #sd(log(dat_fled[28:66,2]))

dat_time <- 
  data.frame(t = c(ab1, ab2, ab3), 
             year = rep(1979:2017, 3),
             group = rep(c("Colony Abundance\nPrediction (SFTS)", 
                           "Breeding Success", 
                           "Fledgling Abundance"), 
                         each = length(1979:2017)))

ggplot(data = dat_time) +
  geom_line(aes(x = year, y = t, col = group, linetype = group), size = 1.1,
            alpha = 0.8) +
  labs(y = "Standardized Time Series\n", x = "Years") +
  scale_color_manual(values = c("Colony Abundance\nPrediction (SFTS)" = "#FF800E", 
                                "Breeding Success" = "darkolivegreen",
                                "Fledgling Abundance" = "#006BA4")) +
  scale_linetype_manual(values = c("Colony Abundance\nPrediction (SFTS)" = 1, 
                                   "Breeding Success" = 2,
                                   "Fledgling Abundance" = 3)) +
  scale_x_continuous(breaks = seq(1980, 2018, by = 6)) +
  theme(panel.border = element_blank(),
        #axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")

ggsave("fig4.pdf", width = 150, height = 120, units = "mm", dpi = 600)
ggsave("fig4.jpeg", width = 150, height = 120, units = "mm", dpi = 600)

ggplot(data = filter(dat_time, group == "Breeding Success")) +
  geom_line(aes(x = year, y = t), size = 1.1, alpha = 0.8, 
            col = "darkolivegreen") +
  theme(panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        #axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank())
ggsave("ts1.pdf", width = 3, height = 2, units = "in") 

ggplot(data = filter(dat_time, group == "Fledgling Abundance")) +
  geom_line(aes(x = year, y = t), size = 1.1, alpha = 0.8, 
            col = "#006BA4") +
  theme(panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        #axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank())
ggsave("ts2.pdf", width = 3, height = 2, units = "in") 


# Random forest -----------------------------------------------------------

library(randomForest)
library(caret)
library(pdp)
library(permimp)

res_rf <- randomForest(y = dat_finn$y, x = dat_finn$X, ntree = 2000,
                       keep.forest = T, keep.inbag = T)
rfStats(res_rf)

var_imp <- permimp(res_rf)$values
var_imp <- var_imp[order(var_imp, decreasing = T)]

# Plot variable importance
plot_imp <- function(varimp) {
  theme_set(theme_bw())
  ggplot(mapping = aes(x = factor(names(varimp), 
                                  levels = names(varimp)),
                       y = varimp)) +
    geom_col(fill = "dark blue", alpha = 0.8) +
    labs(y = "Variable Importance") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 6, angle = 45, vjust = 0.8),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text.y = element_text(size = 8))
}

plot_imp(var_imp)
ggsave("fig_S2.jpeg", width = 6, height = 5, units = "in") 

# Plot partial dependence plots
par_N1 <- partial(res_rf, pred.var = c("aice_nonbreed"), 
                  rug = T, progress = "text")

par_N2 <- partial(res_rf, pred.var = c("fdice_rearing"), 
                  rug = T, progress = "text")

ggplot() +
  geom_line(data = par_N1, aes(x = aice_nonbreed, y = yhat), 
             col = "firebrick4", alpha = 0.9) +
  labs(y = "Colony Abundance (log)", 
       x = "SIC (Nonbreed)") +
  theme(panel.border = element_blank(),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.grid = element_blank())

ggsave("fig_S3.jpeg", width = 6, height = 5, units = "in") 

ggplot() +
  geom_line(data = par_N2, aes(x = fdice_rearing, y = yhat), 
            col = "firebrick4", alpha = 0.9) +
  labs(y = "Colony Abundance (log)", 
       x = "NOW") +
  theme(panel.border = element_blank(),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.grid = element_blank())

ggsave("fig_S4.jpeg", width = 6, height = 5, units = "in") 
