require(dplyr)
require(readr)
require(tsibble)
require(feasts)
require(lubridate)
require(ggplot2)
require(bsts)
require(tidyr)
require(RColorBrewer)


# ------ Data loading and cleaning ------

# Load data:
dat <- read_csv("./data/un_data.csv") %>% 
  setNames(., c("year", "month", "y")) %>% 
  dplyr::filter(!month %in% c("Total", "Unknown"))

# Cleaning and other preparations:
dat <- dat %>% 
  mutate(dm = paste(year, substr(month, 1, 3)) %>% yearmonth(.)) %>% 
  as_tsibble(., index = dm) %>% 
  fill_gaps() %>% 
  mutate(y = y /10000) %>% # смертность, выраженная в десятках тыс. человек
  dplyr::select(y, dm)

# Pre-covid data:
pre_covid <- dat %>% 
  filter(dat$dm < yearmonth("2020 Apr"))

# Training, validation and test sets:
n_pre_covid <- nrow(pre_covid)
train <- pre_covid[1:(n_pre_covid - 6), ]
valid <- pre_covid[(n_pre_covid - 5):n_pre_covid, ]
test <- dat %>% dplyr::filter(dat$dm >= yearmonth("2020 Apr"))


# ------ Exploratory analysis ------

# Figure 1:
ggplot(dat, aes(dm, y)) +
  geom_line() +
  geom_smooth(se = FALSE) +
  theme_minimal() +
  ylim(c(0, 1.5)) + 
  labs(x = "Time", y = "Mortality (tens of thousands)")


# Figure 2
train %>% gg_season(y, pal = brewer.pal(6, "OrRd")) + 
  ylim(c(0, NA)) +
  theme_minimal() + 
  labs(x = "Month", y = "Mortality (tens of thousands)")


# Figure 3:
train %>% 
  mutate(y = y) %>% 
  gg_lag(geom = "point", lags = 1:6,
         col = "black", alpha = 0.4) + 
  theme_minimal()


# ------ Models ------

# Model with trend, seasonal and autoregressive components:
y <- train$y

ss <- list()
ss <- AddLocalLinearTrend(ss, y)
ss <- AddSeasonal(ss, y, nseasons = 12)
ss <- AddAutoAr(ss, y, lag = 3)
m0 <- bsts(y, ss, niter = 2000, ping = 0, seed = 42)

# Goodness of fit metrics:
summary(m0)

# Figure 4:
par(mar = c(5.1, 4.1, 1, 1))
plot(m0, ylab = "Mortality (tens of thousands)",
     xlab = "Observation index")

# Figure 5:
plot(m0, y = "components", same.scale = TRUE)


# Model without the autoregressive component:
ss <- list()
ss <- AddLocalLinearTrend(ss, y)
ss <- AddSeasonal(ss, y, nseasons = 12)
m1 <- bsts(y, ss,
           niter = 2000, ping = 0, seed = 42)

# Goodness of fit metrics:
summary(m1)

# Figure 6:
models_to_compare <- list("m0" = m0, "m1" = m1)
CompareBstsModels(models_to_compare, colors = c("black", "red"))


# Quality of predictions on validation set:

# Utility function to calculate MAE:
mae <- function(observed, predicted){
  mean(abs(observed - predicted))
}

# MAE for models m0 and m1:
sapply(list("m0" = (predict(m0, horizon = 6) %>% .$median)*10000, 
            "m1" = (predict(m1, horizon = 6) %>% .$median)*10000), 
       mae, observed = valid$y * 10000) %>% round()


# Figure 7:
m0_pred <- predict(m0, horizon = 6)
point_predictions <- tibble(dm = valid$dm,
                            y = m0_pred$median,
                            ll95 = m0_pred$interval[1, ],
                            ul95 = m0_pred$interval[2, ]) %>% 
  as_tsibble(., index = dm)

bind_rows(train, valid) %>% 
  dplyr::filter(dm >= yearmonth("2015 Oct"),
                dm < yearmonth("2019 Oct")) %>% 
  ggplot(., aes(dm, y)) + geom_line() +
  geom_ribbon(data = point_predictions,
              aes(ymin = ll95, ymax = ul95),
              fill = "lightblue") +
  geom_line(data = point_predictions, col = "blue") +
  geom_point(data = valid, aes(dm, y), col = "red") +
  ylim(c(0, NA)) +
  theme_minimal() +
  labs(x = "Time", y = "Mortality (tens of thousands)")


# Fitting the final model:
y_full <- c(y, valid$y)

ss <- list()
ss <- AddLocalLinearTrend(ss, y_full)
ss <- AddSeasonal(ss, y_full, nseasons = 12)
ss <- AddAutoAr(ss, y_full, lag = 3)
m_final <- bsts(y_full, ss,
                niter = 2000, ping = 0, seed = 42)


# Predicting the next 3 months:
m_final_pred <- predict(m_final, horizon = 3)


# Figure 8:
point_predictions <- tibble(dm = test$dm,
                            y = m_final_pred$median,
                            ll95 = m_final_pred$interval[1, ],
                            ul95 = m_final_pred$interval[2, ]) %>% 
  as_tsibble(., index = dm)

bind_rows(train, valid) %>% 
  dplyr::filter(dm >= yearmonth("2015 Apr")) %>% 
  ggplot(., aes(dm, y)) + geom_line() +
  geom_ribbon(data = point_predictions,
              aes(ymin = ll95, ymax = ul95),
              fill = "lightblue") +
  geom_line(data = point_predictions, col = "blue") +
  geom_point(data = test, aes(dm, y), col = "red") +
  ylim(c(0, NA)) +
  theme_minimal() +
  labs(x = "Time", y = "Mortality (tens of thousands")


# ------ Excess deaths calculations ------

# For convinience, let's convert monthly mortality values
# to the original scale by multiplying by 10000:
m_dist <- (m_final_pred$distribution*10000) %>% round()

dim(m_dist)


# Absolute growth in mortality:
covid <- test$y

covid_y <- test$y * 10000
excess_mortality <- apply(m_dist, 1, 
                          FUN = function(x){covid_y - x}) %>% t()

excess_mortality <- excess_mortality %>% 
  as_tibble(.name_repair = "minimal") %>% 
  setNames(., c("apr", "may", "jun")) %>% 
  rowwise() %>% 
  mutate(total = sum(c(apr, may, jun)))

result_abs <- excess_mortality %>%
  dplyr::select(apr:total) %>% 
  pivot_longer(cols = apr:total) %>% 
  group_by(month = name) %>% 
  summarise(tibble(
    med = median(value),
    ll95 = quantile(value, p = 0.025),
    ul95 = quantile(value, p = 0.975)
  )) %>% 
  arrange(med)

# Figure 9:
result_abs %>% 
  dplyr::filter(month != "total") %>% 
  mutate(month = factor(month, ordered = TRUE,
                        levels = c("apr", "may", "jun"))) %>% 
  ggplot(., aes(month, med)) +
  geom_point(size = 4) +
  geom_hline(aes(yintercept = 0), linetype = 2, col = "gray60") +
  geom_errorbar(aes(ymin = ll95, ymax = ul95), width = 0.1) +
  theme_minimal() +
  labs(x = "Month of 2020", y = "Number of excess deaths")


# Percentage change in mortality:
covid_sum <- sum(covid_y)

excess_mortality_rel <- apply(m_dist, 1, 
                              FUN = function(x){
                                monthly <- (covid_y - x) / x * 100
                                x_sum <- sum(x)
                                total <- (covid_sum - x_sum) / x_sum * 100
                                c(monthly, total)
                              }) %>% t() %>% 
  as_tibble(.name_repair = "minimal") %>% 
  setNames(., c("apr", "may", "jun", "total"))

result_rel <- excess_mortality_rel %>%
  dplyr::select(apr:total) %>% 
  pivot_longer(cols = apr:total) %>% 
  group_by(month = name) %>% 
  summarise(tibble(
    med = median(value),
    ll95 = quantile(value, p = 0.025),
    ul95 = quantile(value, p = 0.975)
  )) %>% 
  mutate(month = factor(month, ordered = TRUE,
                        levels = c("apr", "may", "jun", "total"))) %>% 
  arrange(month)


# Figure 10:
result_rel %>% 
  dplyr::filter(month != "total") %>% 
  ggplot(., aes(month, med)) +
  geom_point(size = 4) +
  geom_hline(aes(yintercept = 0), linetype = 2, col = "gray60") +
  geom_errorbar(aes(ymin = ll95, ymax = ul95), width = 0.15) +
  theme_minimal() +
  labs(x = "Month of 2020", y = "Excess deaths (% increase)")


# Excess death rate per 100,000 people:

# Populations size:
population <- 9449323

# Most likely value:
(6731 * 100000) / population

# Lower limit:
(4189 * 100000) / population

# Upper limit:
(9335 * 100000) / population
