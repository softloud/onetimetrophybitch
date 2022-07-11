set.seed(42)
library(tidyverse)
library(viridis)

obs <- 100


simdat <- 
    tibble(
        x = rbeta(obs, 5, 2),
        discipline = "A"
    ) %>%   
    bind_rows(
        tibble(
            x = rbeta(obs, 2, 6),
            discipline = "B"
        ) 
    ) 

quants <-    
    simdat %>% 
    group_by(discipline) %>% 
    summarise(
        first = quantile(x, 0.25),
        second = quantile(x, 0.5),
        third = quantile(x, 0.75)
    )


plotdat <- 
    simdat %>% 
    left_join(quants) %>% 
    mutate(
        effort = case_when(
            x > third ~ "publishes, teaches, service",
            x > second ~ "publishes and teaches",
            x > first ~ "publishes",
            TRUE ~ "graduates"
        )
    )




plotdat  %>%
    ggplot(aes(x = discipline, y = x)) +
    geom_boxplot(
        alpha = 0.8
    ) +
    geom_jitter(aes(colour = effort), alpha = 0.7) +
    labs(
        title = "Discipline matters more than effort",
        subtitle = "When considering postgraduate study, consider the opportunities across domains of interest" %>% 
            str_wrap(60),
        x = "Made up disciplines (these are randomly generated data)" %>% 
            str_wrap(30),
        y = "Totally made up probability of obtaining a domain position" %>% 
            str_wrap(30),
        caption = "Think of each point as a person; an aspiring scholar with a passion. No matter how hard a scholar in Discipline B optimises within discipline, according to various advice provided by mentors (e.g., teach, publish, contribute to community), that scholar will at best achieve the opportunities available to the lowest quartile of Discipline A." %>% 
            str_wrap(80)
    ) +
    ylim(0, 1) +
    theme_minimal(
        base_size = 18,
        base_family = "serif"
    ) +
    theme(
        legend.direction = "horizontal",
        legend.position = "top",
        axis.text.y = element_blank()
    ) +
    scale_color_viridis(
        direction = -1,
        discrete = TRUE,
        guide =
    guide_legend(
        ncol = 2
    )
    ) #-> this_plot

# write_rds(this_plot, "_posts/2022-07-09-ikigai/discipline-vis.rds")

# scale_color_discrete(        guide =
#                                  guide_legend(
#                                      ncol = 2)
# ) +
# 
# 
# scale_color_brewer(palette = "Dark2",
# ) 
