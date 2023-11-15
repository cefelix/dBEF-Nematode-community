#checking out hurdle models:

#https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/#normally-distributed-outcomes-with-zeros


#### libraries ####
library(tidyverse)       # ggplot, dplyr, and friends
library(brms)            # Bayesian modeling through Stan
library(emmeans)         # Calculate marginal effects in fancy ways
library(tidybayes)       # Manipulate Stan objects in a tidy way
library(broom)           # Convert model objects to data frames
library(broom.mixed)     # Convert brms model objects to data frames
library(scales)          # For formatting numbers with commas, percents, and dollars
library(patchwork)       # For combining plots
library(ggh4x)           # For nested facets in ggplot
library(ggtext)          # Use markdown and HTML in ggplot text
library(MetBrewer)       # Use pretty artistic colors
library(gapminder)       # Country-year panel data from the Gapminder project
library(palmerpenguins)  # Penguin data!

# Use the cmdstanr backend for Stan because it's faster and more modern than the
# default rstan You need to install the cmdstanr package first
# (https://mc-stan.org/cmdstanr/) and then run cmdstanr::install_cmdstan() to
# install cmdstan on your computer.
options(mc.cores = 4,
        brms.backend = "cmdstanr")


#### explore data ####
penguins <- palmerpenguins::penguins |> 
  drop_na(sex) |> 
  # Make a bunch of weight values 0
  mutate(prob_zero = ifelse(flipper_length_mm < 190, 0.3, 0.02),
         will_be_zero = rbinom(n(), 1, prob = prob_zero),
         body_mass_g = ifelse(will_be_zero, 0, body_mass_g)) |> 
  select(-prob_zero, -will_be_zero) |> 
  mutate(is_zero = body_mass_g == 0)

#how many zeros in data?
penguins |> 
  count(is_zero) |> 
  mutate(prop = n / sum(n) )

#check the distribution:
penguins |> 
  mutate(body_mass_g = ifelse(is_zero, -0.1, body_mass_g)) |> 
  ggplot(aes(x = body_mass_g)) +
  geom_histogram(aes(fill = is_zero), binwidth = 100,
                 boundary = 0, color = "white") +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c("darkorange", "purple")) +
  scale_x_continuous(labels = comma_format()) +
  labs(x = "Body mass (g)", y = "Count", fill = "Is zero?",
       title = "Distribution of penguin body mass") +
  theme(legend.position = "bottom")

# Set some global Stan options
CHAINS <- 4
ITER <- 2000
WARMUP <- 1000
BAYES_SEED <- 1234

# Use the Johnson color palette
clrs <- MetBrewer::met.brewer("Johnson")

# Tell bayesplot to use the Johnson palette (for things like pp_check())
bayesplot::color_scheme_set(c("grey30", clrs[2], clrs[1], clrs[3], clrs[5], clrs[4]))

# Custom ggplot theme to make pretty plots
# Get the font at https://fonts.google.com/specimen/Jost
theme_nice <- function() {
  theme_minimal(base_family = "Jost") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(family = "Jost", face = "bold"),
          axis.title = element_text(family = "Jost Medium"),
          strip.text = element_text(family = "Jost", face = "bold",
                                    size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA))
}


#### regular OLS model on non-exponential outcome ####

ggplot(penguins, aes(x = bill_length_mm, y = body_mass_g)) +
  geom_point(aes(color = species), size = 1.5) + 
  geom_smooth(method = "lm", color = "#F012BE") +
  geom_smooth(data = filter(penguins, body_mass_g != 0), method = "lm", color = "#0074D9") +
  scale_y_continuous(labels = label_comma()) +
  scale_color_manual(values = c("darkorange", "purple", "cyan4")) +
  labs(x = "Bill length (mm)", y = "Body mass (g)", color = "Species",
       title = 'OLS models <span style="color:#F012BE;">with</span> and <span style="color:#0074D9;">without</span> zeros') +
  theme_nice() +
  theme(legend.position = "bottom",
        plot.title = element_markdown())

#models with and without zeros:
model_mass_basic <- brm(
  bf(body_mass_g ~ bill_length_mm + species),
  data = penguins,
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  silent = 2
)
pp_check(model_mass_basic)

model_mass_basic_no_zero <- brm(
  bf(body_mass_g ~ bill_length_mm + species),
  data = filter(penguins, !is_zero),
  family = gaussian(),
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  silent = 2
)


#### hurdle lognormal model on a non-exponential outcome ####

penguins |> 
  mutate(body_mass_g = log1p(body_mass_g)) |> 
  mutate(body_mass_g = ifelse(is_zero, -0.01, body_mass_g)) |> 
  ggplot(aes(x = body_mass_g)) +
  geom_histogram(aes(fill = is_zero), binwidth = 0.1,
                 boundary = 0, color = "white") +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c("darkorange", "purple")) +
  labs(x = "Body mass (g)", y = "Count", fill = "Is zero?",
       title = "Distribution of logged penguin body mass") +
  theme_nice() +
  theme(legend.position = "bottom")

#check out non-transformed bodymasses:
penguins$body_mass_g %>% density %>% plot()

model_mass_hurdle_log <- brm(
  bf(body_mass_g ~ bill_length_mm + species,
     hu ~ flipper_length_mm),
  data = penguins,
  family = hurdle_lognormal(),
  cores = 3,
  chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
  silent = 2
)
