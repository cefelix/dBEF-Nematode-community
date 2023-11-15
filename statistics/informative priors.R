#lets try to add priors to narrow the posterior predictions:
get_prior(Fu_per100gLog ~ sowndiv*treatment + (1|block/plot),
          data = dBEF_nem21, family = "gaussian")
summary(m.Fu.12_nestRE_b)
priors <- c(
  prior(uniform(-100, 100), class="b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(student_t(3, 0, 2.5), class = "sigma"),
  prior(student_t(3, 5.6, 2.5), class = "intercept")
) 
