# =========================================================================================
#
# SST assessment Length figures
#
# =========================================================================================
# Pull length data
source("code/length_data_pull.r")

lls.sst.len %>% 
  write_csv(paste0(dat_path, "/goa_sst_lls_lengths", YEAR, ".csv"))
bts.sst.len %>% 
  write_csv(paste0(dat_path, "/goa_sst_bts_lengths", YEAR, ".csv"))
fsh.sst.len %>% 
  write_csv(paste0(dat_path, "/goa_sst_fishery_lengths", YEAR, ".csv"))

# Remove a known bad data point where year = 2017, area = W Yak Slope, Length = 29
# These are lengths from deep applied to a shallow strata (not real)
# Group the LL numbers by year/length
ll.len = lls.sst.len %>% 
  mutate(rpn = ifelse(year == 2017 & area_code == 38 & length == 29, NA, rpn)) %>% 
  group_by(year, length) %>% 
  summarize(freq = sum(rpn, na.rm = T))

ll.len$calc = ll.len$freq*ll.len$length

ll.means = ll.len %>% group_by(year) %>% 
  summarize(tot = sum(freq), l_calc = sum(calc))

ll.means$mean = ll.means$l_calc/ll.means$tot
ll.len = merge(ll.len, ll.means, by=c("year"))

ll.len %>%  
  ggplot(aes(x=length, y=after_stat(density), weighted.mean=freq)) +
  geom_histogram(alpha=0.25, binwidth=1, col="black") +
  facet_wrap(~year, ncol=3) +
  theme(legend.position = "top") +
  xlab("Length (cm)") +
  ylab("Proportion of LL survey RPNs") +
  geom_text(aes(x=18, y=0.075, label=year)) +
  geom_text(aes(x=18, y=0.055, label=paste0("(", format(round(mean, digits=1), nsmall = 1) , " cm)"))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.099)) +
  scale_x_continuous(breaks=seq(15,75,10), limits=c(12,78)) +
  theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=12), strip.background=element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text=element_blank()) 

ggsave(file=paste0("results/", YEAR, "/SST_LLS_Lengths.png"), height = 15, width = 12, dpi=600)

ggplot(ll.len, aes(year, mean)) + geom_point() + geom_line()

ll.len$survey = "LLS"

#####################
# Group the BTS numbers by year/length
len = bts.sst.len %>%
  filter(year < 2025) |> 
  mutate(length = length_mm / 10) %>%
  group_by(year, length) %>%
  summarize(freq = sum(population_count))

len$calc = len$freq*len$length

means = len %>% group_by(year) %>%
  summarize(tot = sum(freq), l_calc = sum(calc))

means$mean = means$l_calc/means$tot
len = merge(len, means, by=c("year"))

len %>%
        ggplot(aes(x=length, y=after_stat(density), weighted.mean=freq)) +
        # theme_linedraw() +
        geom_histogram(alpha=0.25, binwidth=1, col="black") +
        facet_wrap(~year, ncol=2) +
        theme(legend.position = "top") +
        xlab("Length (cm)") +
        ylab("Proportion of trawl survey population") +
        geom_text(aes(x=60, y=0.055, label=year)) +
        geom_text(aes(x=60, y=0.03, label=paste0("(", format(round(mean, digits=1), nsmall = 1) , " cm)"))) +
        scale_y_continuous(expand=c(0,0), limits=c(0,0.09)) +
        scale_x_continuous(expand=c(0,0), breaks=seq(0,110,10)) +
        theme_bw() +
        theme(axis.title=element_text(size=14), axis.text=element_text(size=12), strip.background=element_blank(),
                panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text=element_blank())

# Group the BTS numbers by year/length
len2 = bts.sst.len %>% 
  filter(year < 2025) %>%
  mutate(length = length_mm / 10,
         strata = ifelse(area_id %in% c(10:13, 110:112, 210, 310), 'WGOA (0-500 m)',
                         ifelse(area_id %in% c(20:35, 120:134, 220:232, 32, 320, 330), 'CGOA (0-500 m)',
                                ifelse(area_id %in% c(40:50, 140:151, 240:251, 340:351), 'EGOA (0-500 m)',
                                       ifelse(area_id == 410, 'WGOA (501-700 m)',
                                              ifelse(area_id == 510, 'WGOA (701-1000 m)',
                                                     ifelse(area_id %in% c(420, 430), 'CGOA (501-700 m)',
                                                            ifelse(area_id %in% c(520, 530), 'CGOA (701-1000 m)',
                                                                   ifelse(area_id %in% c(440, 450), 'EGOA (501-700 m)',
                                                                          ifelse(area_id %in% c(540, 550), 'EGOA (701-1000 m)', NA)))))))))) %>% 
  filter(!is.na(strata)) %>% 
  group_by(strata, length) %>% 
  summarize(freq = sum(population_count)) 

len2$calc = len2$freq*len2$length

means2 = len2 %>% group_by(strata) %>% 
  summarize(tot = sum(freq), l_calc = sum(calc))

means2$mean = means2$l_calc/means2$tot
len2 = merge(len2, means2, by=c("strata"))

len2 %>%  
  ggplot(aes(x=length, y=after_stat(density), weighted.mean=freq)) +
  # theme_linedraw() +
  geom_histogram(alpha=0.25, binwidth=1, col="black") +
  facet_wrap(~strata, ncol=3) +
  # facet_wrap(~year, ncol=2) +
  theme(legend.position = "top") +
  xlab("Length (cm)") +
  ylab("Proportion of trawl survey population") +
  geom_text(aes(x=60, y=0.055, label=strata)) +
  geom_text(aes(x=60, y=0.03, label=paste0("(", format(round(mean, digits=1), nsmall = 1) , " cm)"))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.115)) +
  scale_x_continuous(expand=c(0,0), breaks=seq(0,110,10)) +
  theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=12), strip.background=element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text=element_blank()) 

ggsave(file = paste0("results/", YEAR, "/SST_BTS_strata_Lengths.png"), height = 10, width = 7, dpi=600)

len$survey = "BTS"

# Look at a time series of the mean length by surveys
ggplot(ll.means, aes(year, mean))  + 
  geom_point(size=4) +
  geom_point(data=means, aes(x=year, y=mean), size=4, pch=21) +
  ylab("Mean length (cm)") +
  scale_x_continuous(breaks=seq(1990,2020,5)) +
  theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=12), strip.background=element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text=element_blank()) 

ggsave(file = paste0("results/", YEAR, "/SST_Time_Series_length_Comp.png"), height = 3.5, width = 6, dpi=600)

# IS there a relationship between mean of LL lengths and mean of BTS lengthS?
comb = merge(ll.means, means, by="year")

ggplot(comb, aes(mean.x, mean.y)) + 
  geom_point()

# NO relationship present

# Summary of entire survey datasets
all.bts = bts.sst.len %>% 
  filter(year < 2025) |> 
  mutate(length = length_mm / 10) %>% 
  group_by(length) %>% 
  summarize(freq = sum(population_count))

all.bts$calc = as.numeric(all.bts$freq)*as.numeric(all.bts$length)
all.bts$survey = "BTS"

# all.bts = all.bts[ , c(3,2,4,5)]

bts.mean = all.bts %>% summarize(tot = sum(freq), l_calc = sum(calc))

bts_mean = data.frame("Mean" = bts.mean$l_calc/bts.mean$tot)

all.bts %>%  
  ggplot(aes(x=length, y=after_stat(density), weighted.mean=freq)) +
  geom_histogram(alpha=0.25, binwidth=1, col="black") +
  theme(legend.position = "top") +
  xlab("Length (cm)") +
  ylab("Proportion of trawl survey population") +
  geom_text(aes(x=10, y=0.05, label=paste0("(", format(round(bts_mean, digits=1), nsmall = 1) , " cm)"))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.07)) +
  scale_x_continuous(breaks=seq(10,110,10)) +
  theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=12), strip.background=element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text=element_blank()) 

# ggsave(file="sst_all_bts_len.png", height = 4, width = 5, dpi=600)

all.ll = lls.sst.len %>% group_by(length) %>% 
  summarize(freq = sum(rpn, na.rm = TRUE))

all.ll$calc = as.numeric(all.ll$freq)*as.numeric(all.ll$length)
all.ll$survey = "LLS"

ll.mean = all.ll %>% 
  summarize(tot = sum(freq), l_calc = sum(calc))

ll_mean = data.frame("Mean" = ll.mean$l_calc/ll.mean$tot)

all.ll %>%  
  ggplot(aes(x=length, y=after_stat(density), weighted.mean=freq)) +
  geom_histogram(alpha=0.25, binwidth=1, col="black") +
  theme(legend.position = "top") +
  xlab("Length (cm)") +
  ylab("Proportion of LL survey RPNs") +
  geom_text(aes(x=10, y=0.06, label=paste0("(", format(round(ll_mean, digits=1), nsmall = 1) , " cm)"))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.085)) +
  scale_x_continuous(breaks=seq(10,110,10)) +
  theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=12), strip.background=element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text=element_blank()) 

# ggsave(file="sst_all_ll_len.png", height = 4, width = 5, dpi=600)

# Combine for plot
all = rbind(all.bts, all.ll)

ll_mean$survey = "LLS"
bts_mean$survey = "BTS"
all.mean = rbind(ll_mean, bts_mean)

ggplot() +
  geom_histogram(data=all, aes(x=length, y=after_stat(density), weighted.mean=freq, fill = survey),
                 position = 'identity', alpha=0.5, binwidth=1, col="black") +
  theme(legend.position = "top") +
  scale_fill_discrete(type = c('red', 'blue')) +
  xlab("Length (cm)") +
  ylab("Length composition by survey") +
  geom_text(data=all.mean, aes(x=c(70,70),  y=c(0.055, .045), label=paste0(survey, " mean length: ", format(round(Mean, digits=1), nsmall = 1) , " cm"))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.085)) +
  scale_x_continuous(expand=c(0,5), breaks=seq(5,85,10)) +
  # facet_wrap(~survey, ncol=1) +
  theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=12), strip.background=element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text=element_blank()) 

ggsave(file = paste0("results/", YEAR, "/all_len.png"), height = 4, width = 6, dpi=600)
# Alternative using ggridges

# Fishery lengths
lengths <- fsh.sst.len %>% 
  filter(gear == 1 | gear == 8) %>%
  filter(nmfs_area > 609, nmfs_area < 651, !nmfs_area == 649) %>% 
  mutate(Gear = ifelse(gear == 1, "Trawl", "Longline"))

mean_lengths <- lengths %>% 
  group_by(Gear) %>% 
  summarize(Mean = mean(length))

ggplot(lengths) + 
  geom_histogram(aes(x=length, y=after_stat(density), weighted.mean=frequency, fill = Gear),
                 position = 'identity', alpha=0.5, binwidth=1, col="black") +
  scale_fill_discrete(type = c('blue', 'red')) +
  xlab("Length (cm)") +
  ylab("Length composition by gear type") +
  geom_text(data=mean_lengths, aes(x=c(60,60),  y=c(0.055, .05), label=paste0(Gear, " mean length: ", format(round(Mean, digits=1), nsmall = 1) , " cm"))) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.075)) +
  scale_x_continuous(expand=c(0,5), breaks=seq(5,85,10)) +
  theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=12), strip.background=element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.text=element_blank()) 

    
ggsave(file = paste0("results/", YEAR, "/all_fish_len.png"), height = 4, width = 7, dpi=600)
