# Reshape proportions and analyze changes associated with treatment
library(tidyverse)

load(file = "./output/C8/linseedProportions.RData")
cellProportions2 <- t(cellProportions) %>%
  as.data.frame(.) %>%
  rownames_to_column(.) %>% 
  as_tibble(.)
names(cellProportions2) <- str_replace_all(string = names(cellProportions2), pattern = " ", replacement = "")


# Diagnostic plot: distribution of cell type proportions
as.data.frame(cellProportions) %>%
  rownames_to_column() %>%
  gather(-1 ,key = 'ID', value = 'Proportion') %>%
  ggplot(., mapping = aes(x = Proportion, fill = rowname)) +
  geom_density(alpha = 0.3)


tidyprops <- gather(cellProportions2, key = "CellType", value = "prop", ... = -rowname) %>%
  mutate(.,
         subj = str_extract(rowname, pattern = "[0-9]{1,2}") %>%
           as.factor(),
         cond = str_replace_all(str_extract(rowname, "_(L|P)_"), "_", "") %>%
           factor(levels = c('P', 'L')),
         week = str_extract(rowname, pattern = "W[:digit:]") %>%
           factor()
  )

## Cell type proportions before and after treatment
ggplot(data = tidyprops, aes(x = cond, y = prop, group = subj)) +
  geom_line(aes(alpha=0.3)) +
  facet_wrap(.~CellType)


## Import responder vs non-responder information
metadata <- read_tsv(file = "../Data/Clinical_characteristics_Herlevstudy.txt") %>%
  select(., contains("Screening"), contains("Treatment"), contains("Group")) %>%
  gather(., key = "Visit", value = "Treatment", ... = contains("Treatment"))

metadata$Visit <- str_extract(string = metadata$Visit, pattern = "[0-9]") %>% 
  paste0("W", .) %>%
  as.factor()
metadata$Treatment <- str_extract(string = metadata$Treatment, pattern = "^[A,P]") %>%
  str_replace(., "A", "L") %>%
  as.factor()

names(metadata) <- c('subj', 'respond', 'week', 'cond')

# Merge with data (use treatment for visit 1/2 as check)
metaprops <- merge(tidyprops, metadata, all.x = T, all.y = T) %>% 
  na.exclude() %>%
  as_tibble()

ggplot(data = metaprops, aes(x = cond, y = prop, group = subj)) +
  geom_line(aes(alpha=0.3)) +
  facet_grid(respond~CellType)

# Check if initial cell proportions between groups differ
metaprops_parallel <- select(metaprops, -"week", -"rowname") %>%
  spread(., "cond", "prop") %>%
  group_by(., CellType)

# Test for changes in proportions between control and treatment across each cell type
general_responses <- summarise(metaprops_parallel, 
                               p.prop = wilcox.test(P, L, paired = T)$p.value,
                               ci.prop.min = wilcox.test(P, L, conf.int = T, paired = T)$conf.int[1],
                               ci.prop.max = wilcox.test(P, L, conf.int = T, paired = T)$conf.int[2],
                               p.student = t.test(P, L, paired = T)$p.value
)
general_responses$p.prop <- pmin(general_responses$p.prop * nrow(general_responses), 1)
general_responses$p.student <- pmin(general_responses$p.student * nrow(general_responses), 1)
general_responses

# Test for changes in proportions separating by response group
metaprops_parallel <- group_by(metaprops_parallel, respond, CellType)

grouped_responses <- summarise(metaprops_parallel, 
                               p.wilcox = wilcox.test(P, L, paired = T)$p.value,
                               ci.wilcox.min = wilcox.test(P, L, conf.int = T, paired = T)$conf.int[1],
                               ci.wilcox.max = wilcox.test(P, L, conf.int = T, paired = T)$conf.int[2],
                               p.student = t.test(P, L, paired = T)$p.value
)
grouped_responses$p.wilcox <- pmin(grouped_responses$p.wilcox * nrow(grouped_responses), 1)
grouped_responses$p.student <- pmin(grouped_responses$p.student * nrow(grouped_responses), 1)

grouped_responses

# Using GLMs
glm1 <- filter(metaprops, respond!="Outlier") %>%
  glm(data = ., formula = prop ~  CellType/(cond*respond), family = 'Gamma')
plot(glm1)
summary(glm1)

# Test for differences in initial cell proportions
PlaceboProps <- select(metaprops, -"week", -"rowname") %>%
  filter(., cond == "P") %>%
  spread(., respond, prop, ) %>%
  group_by(., CellType)

summarise(PlaceboProps, 
          p.prop = wilcox.test(Group1, Group2)$p.value
)

# Test whether responder population has greater changes in cell populations than non-responders

deltaProps <- mutate(metaprops_parallel, 
                     diff = P-L
) %>%
  spread(., respond, diff) %>%
  summarise(., 
            mean.group1 = mean(Group1, na.rm = T),
            mean.group2 = mean(Group2, na.rm = T),
            sd.group1 = sd(Group1, na.rm = T),
            sd.group2 = sd(Group2, na.rm = T),
            p.prop = wilcox.test(Group1, Group2)$p.value,
            ci.l.prop = wilcox.test(Group1, Group2, conf.int = T)$conf.int[1],
            ci.u.prop = wilcox.test(Group1, Group2, conf.int = T)$conf.int[2],
  )
deltaProps
