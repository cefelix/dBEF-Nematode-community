# a coding example to understand how to model zero inflated lognormal data

#from: https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/

# load data####
  
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
  
  set.seed(1234)
  gapminder <- gapminder::gapminder |> 
    filter(continent != "Oceania") |> 
    # Make a bunch of GDP values 0
    mutate(prob_zero = ifelse(lifeExp < 50, 0.3, 0.02),
           will_be_zero = rbinom(n(), 1, prob = prob_zero),
           gdpPercap = ifelse(will_be_zero, 0, gdpPercap)) |> 
    select(-prob_zero, -will_be_zero) |> 
    # Make a logged version of GDP per capita
    mutate(log_gdpPercap = log1p(gdpPercap)) |> 
    mutate(is_zero = gdpPercap == 0)


#explore data####
  gapminder |> 
    count(is_zero) |> 
    mutate(prop = n / sum(n))
  
  plot_dist_unlogged <- gapminder |> 
    mutate(gdpPercap = ifelse(is_zero, -0.1, gdpPercap)) |> 
    ggplot(aes(x = gdpPercap)) +
    geom_histogram(aes(fill = is_zero), binwidth = 5000, 
                   boundary = 0, color = "white") +
    geom_vline(xintercept = 0) + 
    scale_x_continuous(labels = label_dollar(scale_cut = cut_short_scale())) +
    scale_fill_manual(values = c(clrs[4], clrs[1]), 
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = "GDP per capita", y = "Count", fill = "Is zero?",
         subtitle = "Nice and exponentially shaped, with a bunch of zeros") +
    theme_nice() +
    theme(legend.position = "bottom")
  
  plot_dist_logged <- gapminder |> 
    mutate(log_gdpPercap = ifelse(is_zero, -0.1, log_gdpPercap)) |> 
    ggplot(aes(x = log_gdpPercap)) +
    geom_histogram(aes(fill = is_zero), binwidth = 0.5, 
                   boundary = 0, color = "white") +
    geom_vline(xintercept = 0) +
    scale_x_continuous(labels = label_math(e^.x)) +
    scale_fill_manual(values = c(clrs[4], clrs[1]), 
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = "GDP per capita", y = "Count", fill = "Is zero?",
         subtitle = "Nice and normally shaped, with a bunch of zeros;\nit's hard to interpret intuitively though") +
    theme_nice() +
    theme(legend.position = "bottom")
  
  (plot_dist_unlogged | plot_dist_logged) +
    plot_layout(guides = "collect") +
    plot_annotation(title = "GDP per capita, original vs. logged",
                    theme = theme(plot.title = element_text(family = "Jost", face = "bold"),
                                  legend.position = "bottom"))
  
# OLS using unlogged data####
  
  #a plot with OLS geom_smooth:
  ggplot(gapminder, aes(x = lifeExp, y = gdpPercap)) +
    geom_point(aes(color = continent), size = 1, alpha = 0.5) + 
    geom_smooth(method = "lm", color = "#0074D9") +
    scale_y_continuous(labels = label_dollar(scale_cut = cut_short_scale())) +
    scale_color_manual(values = clrs) +
    labs(x = "Life expectancy", y = "GDP per capita", color = NULL,
         title = "OLS model fit on unlogged GDP") +
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
    theme_nice() +
    theme(legend.position = "bottom")
  
  #the corresponding model:
  model_gdp_basic <- brm(
    bf(gdpPercap ~ lifeExp),
    data = gapminder,
    family = gaussian(),
    chains = CHAINS, iter = ITER, warmup = WARMUP, seed = BAYES_SEED,
    silent = 2
  )
  
  tidy(model_gdp_basic)
  