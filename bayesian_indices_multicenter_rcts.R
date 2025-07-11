rm(list=ls())
library(tidyverse)
library(docxtractr)
library(rstanarm)
library(bayestestR)
library(insight)
library(datawizard)
library(see)
options(scipen=999)

df <- read_docx("test.docx")

# Figure 1
# Illustration of the indices for study 12

rm(dummy,mod,perf)
dummy.rep <- data.frame()
dummy <- docx_extract_tbl(df, 12)[1:2,1:3]
dummy$Outcome.Positive <- as.numeric(dummy$Outcome.Positive)
dummy$Outcome.Negative <- as.numeric(dummy$Outcome.Negative)
fisher.test(dummy[1:2,2:3])$p.value
for (id in 1:as.integer(dummy[1,2])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Intervention",Outcome = 1))}
for (id in 1:as.integer(dummy[1,3])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Intervention",Outcome = 0))}
for (id in 1:as.integer(dummy[2,2])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Control",Outcome = 1))}
for (id in 1:as.integer(dummy[2,3])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Control",Outcome = 0))}

mod  <- stan_glm(Outcome ~ Treatment,data = dummy.rep, family = "binomial",chains = 4, iter = 20000, warmup = 5000,refresh = 0)
fig1A <- plot(pd(mod))+
  theme_classic()+
  scale_y_discrete(expand = c(0,0))+
  scale_x_continuous(breaks = seq(-.4,1.4,by=0.2))+
  theme(
    axis.text.y=element_blank(),
    legend.position = c(0.8,0.8),
    text = element_text(size=14)
  )+
  ylab("Posterior distribution")+
  xlab("Treatment effect (log odds ratio scale)")+
  ggsci::scale_color_jama()+
  ggsci::scale_fill_jama()
fig1B <- plot(rope(mod))+
  theme_classic()+
  scale_y_discrete(expand = c(0,0))+
  scale_x_continuous(breaks = seq(-.4,1.4,by=0.2))+
  theme(
    axis.text.y=element_blank(),
    legend.position = c(0.8,0.8),
    text = element_text(size=14)
  )+
  ylab("")+
  xlab("Treatment effect (log odds ratio scale)")+
  ggsci::scale_color_npg()+
  ggsci::scale_fill_npg()

ggpubr::ggarrange(
  fig1A,
  fig1B,
  labels = c("A","B"),
  nrow=1
)

ggsave("figure1.jpeg",dpi=600,width=9.5,height = 4.6)

perf <- describe_posterior(mod)

n.studies <- 56
df.indices <- data.frame()
for (mystudy in 1:n.studies){
  print(mystudy)
rm(dummy,mod,perf)
dummy.rep <- data.frame()
dummy <- docx_extract_tbl(df, mystudy)[1:2,1:3]
dummy$Outcome.Positive <- as.numeric(dummy$Outcome.Positive)
dummy$Outcome.Negative <- as.numeric(dummy$Outcome.Negative)
fisher.test(dummy[1:2,2:3])$p.value
for (id in 1:as.integer(dummy[1,2])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Intervention",Outcome = 1))}
for (id in 1:as.integer(dummy[1,3])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Intervention",Outcome = 0))}
for (id in 1:as.integer(dummy[2,2])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Control",Outcome = 1))}
for (id in 1:as.integer(dummy[2,3])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Control",Outcome = 0))}

mod  <- stan_glm(Outcome ~ Treatment,data = dummy.rep, family = "binomial",chains = 4, iter = 20000, warmup = 5000,refresh = 0)
perf <- describe_posterior(mod)
df.indices <- rbind(
  df.indices,
  data.frame(mystudy,N = sum(dummy[1:2,2:3]),p = fisher.test(dummy[1:2,2:3])$p.value,pd = perf$pd[2],rope = perf$ROPE_Percentage[2])
)
}


df.indices %>% 
  mutate(p.cat = factor(case_when(
    p<0.001~"p<0.001",
    p>=0.001 & p<0.05~"0.001≤p<0.05",
    p>=0.05 & p<0.1~"0.05≤p<0.1",
    p>=0.1~"p≥0.1"
  ),levels = c("p<0.001","0.001≤p<0.05","0.05≤p<0.1","p≥0.1"))) %>% 
  rename(`P value` = p.cat) %>% 
  ggplot(aes(y=pd,x=rope,size=`P value`,color=`P value`,label=mystudy))+
  geom_point(alpha=0.1,show.legend = T)+#(shape=1)+
  geom_text(show.legend = F)+
  ggsci::scale_color_nejm()+
  guides(color= guide_legend(override.aes = list(alpha = 1)), size=guide_legend())+
  theme_classic()+
  theme(
    text = element_text(size=14),
    legend.position = c(0.2,0.26)
    )+
  scale_x_continuous(labels = scales::percent,n.breaks = 10,limits = c(0,1))+
  scale_y_continuous(labels = scales::percent,n.breaks = 10,limits = c(0.5,1))+
  ylab("Probability of Direction (PD)")+
  xlab("Percentage in Region of Practical Equivalence (ROPE 95%)")

ggsave("figure2.jpeg",dpi=600,width=9,height = 6)

df.indices %>% 
  select(mystudy,p) %>% 
  summarise(
    mymedian = median(p),
    mylower  = quantile(p,c(0.25)),
    myupper  = quantile(p,c(0.75))
  ) %>% round(2)

df.indices %>% 
  select(mystudy,pd) %>% 
  summarise(
    mymedian = median(pd),
    mylower  = quantile(pd,c(0.25)),
    myupper  = quantile(pd,c(0.75))
  ) %>% round(2)

df.indices %>% 
  select(mystudy,rope) %>% 
  summarise(
    mymedian = median(rope),
    mylower  = quantile(rope,c(0.25)),
    myupper  = quantile(rope,c(0.75))
  ) %>% round(2)

