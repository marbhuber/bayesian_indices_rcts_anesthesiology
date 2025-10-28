#######################################################################################################
# Illustrating Bayesian indices of effect existence and practical significance in Anesthesiology trials
# Markus Huber (markus.huber@insel.ch)
# 16 Octo 2025
#######################################################################################################

# setup
setwd("~/Desktop/research/Bayes/anesthesiology_bayesian_indices/revisions")
rm(list=ls())
library(tidyverse)
library(docxtractr)
library(rstanarm)
library(bayestestR)
library(insight)
library(datawizard)
library(see)
library(fitdistrplus)
options(scipen=999)

set.seed(1234)
df <- read_docx("test.docx")

##########################################
##########################################
# Figure 1
# Illustration of the indices for study 12
##########################################
##########################################

rm(dummy,mod,perf)
dummy.rep              <- data.frame()
dummy                  <- docx_extract_tbl(df, 12)[1:2,1:3]
dummy$Outcome.Positive <- as.integer(dummy$Outcome.Positive)
dummy$Outcome.Negative <- as.integer(dummy$Outcome.Negative)

mytest = fisher.test(dummy[1:2,2:3])
# p-value
mytest$p.value # 0.057
# odds ratio
with(mytest,
     paste0(
       format(round(estimate,2),nsmall=2),
       " 95%-CI: ",
       format(round(conf.int[1],2),nsmall=2),
       " - ",
       format(round(conf.int[2],2),nsmall=2)))

# create data frame (1 is positive outcome)
for (id in 1:as.integer(dummy[1,2])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Intervention",Outcome = 1))}
for (id in 1:as.integer(dummy[1,3])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Intervention",Outcome = 0))}
for (id in 1:as.integer(dummy[2,2])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Control",Outcome = 1))}
for (id in 1:as.integer(dummy[2,3])){dummy.rep <- rbind(dummy.rep,data.frame(Treatment = "Control",Outcome = 0))}

# fit bayesian logistic regression model
mod  <- stan_glm(Outcome ~ Treatment,data = dummy.rep, family = "binomial",chains = 4, iter = 20000, warmup = 5000,refresh = 0)

# posterior (odds rather than log odds, but normal fit on log odds scale)
mod.posterior  <- get_parameters(mod)$TreatmentIntervention
fit.posterior  <- fitdist(mod.posterior, "norm") # mean: 0.4336, sd: 0.226
posterior.mean <- fit.posterior$estimate[1] %>% as.numeric()
posterior.sd   <- fit.posterior$estimate[2] %>% as.numeric()

or = seq(0,5,by=0.0001)

fit.posterior.or    <- fitdist(exp(mod.posterior), "lnorm") 
sample.posterior.or <- dlnorm(or,fit.posterior.or$estimate["meanlog"],fit.posterior.or$estimate["sdlog"])

# (1A): probability of direction

fig1A <- data.frame(or,ymin=0,ymax=sample.posterior.or) %>% 
  mutate(Direction = case_when(
    or<1~"Negative",
    TRUE~"Postive")) %>% 
  ggplot(aes(x=or,ymin=ymin,ymax=ymax,color=Direction,fill=Direction))+
  geom_ribbon(alpha=0.3)+
  theme_classic()+
  theme(
    axis.text.y=element_blank(),
    legend.position = c(0.8,0.8),
    text = element_text(size=14)
  )+
  xlab("Treatment effect (odds ratio)")+
  labs(title = "Probability of direction (pd)")+
  ylab("Posterior distribution")+
  scale_x_continuous(n.breaks = 10,limits = c(0,3.5))+
  ggsci::scale_color_jama()+
  ggsci::scale_fill_jama()+
  geom_vline(xintercept = 1,linetype="dashed",color="darkgrey")

# (1B): ROPE full

df.fig1B <- data.frame(or,ymin=0,ymax=sample.posterior.or) %>% 
  mutate(ROPE = case_when(
    or<exp(-0.18)~"Outside",
    (or>=exp(-0.18) & or<=exp(0.18)) ~"Inside",
     or>exp(0.18)~"Outside")) %>% 
  mutate(cat = case_when(
    or<exp(-0.18)~"1",
    (or>=exp(-0.18) & or<=exp(0.18)) ~"2",
    or>exp(0.18)~"3"))

fig1B <- ggplot()+
  geom_ribbon(data=filter(df.fig1B,cat=="1"),aes(x=or,ymin=ymin,ymax=ymax,color=ROPE,fill=ROPE),alpha=0.3)+
  geom_ribbon(data=filter(df.fig1B,cat=="2"),aes(x=or,ymin=ymin,ymax=ymax,color=ROPE,fill=ROPE),alpha=0.3)+
  geom_ribbon(data=filter(df.fig1B,cat=="3"),aes(x=or,ymin=ymin,ymax=ymax,color=ROPE,fill=ROPE),alpha=0.3)+
  theme_classic()+
  theme(
    axis.text.y=element_blank(),
    legend.position = c(0.8,0.8),
    text = element_text(size=14)
  )+
  xlab("Treatment effect (odds ratio)")+
  labs(title = "ROPE full")+
  ylab("Posterior distribution")+
  scale_x_continuous(n.breaks = 10,limits = c(0,3.5))+
  ggsci::scale_color_bmj()+
  ggsci::scale_fill_bmj()+
  geom_vline(xintercept = exp(-0.18),linetype="dashed",color="darkgrey")+
  geom_vline(xintercept = exp(+0.18),linetype="dashed",color="darkgrey")

ggpubr::ggarrange(
  fig1A,
  fig1B,
  labels = c("A","B"),
  nrow=1
)

ggsave("figure1_revised.jpeg",dpi=600,height=4,width = 8)
ggsave("figure1_revised.pdf",dpi=600,height=4,width = 8)

# credible interval
ci_eti  <- ci(mod, method = "ETI")
df.95cr <- data.frame(x=NA,y=1.4,xmin = ci_eti$CI_low[2],xmax = ci_eti$CI_high[2])

# result <- bayesfactor_parameters(mod, null = 0)
# plot(result)
# 
# +
#   theme_classic()+
#   scale_y_discrete(expand = c(0,0))+
#   theme(
#     axis.text.y=element_blank(),
#     legend.position = c(0.8,0.8),
#     text = element_text(size=14)
#   )+
#   xlab("Treatment effect (odds ratio)")+
#   labs(title = "Prior and posterior distribution")+
#   scale_x_continuous(limits = c(-0.7,1.5),n.breaks = 10)
# 
#   geom_vline(xintercept = -0.7,color="black")
# 
# fig1B <- plot(pd(mod),log_scale = FALSE)+
#   theme_classic()+
#   scale_y_discrete(expand = c(0,0))+
#   theme(
#     axis.text.y=element_blank(),
#     legend.position = c(0.8,0.8),
#     text = element_text(size=14)
#   )+
#   ylab("Posterior distribution")+
#   xlab("Treatment effect (odds ratio)")+
#   ggsci::scale_color_jama()+
#   ggsci::scale_fill_jama()+
#   # geom_vline(xintercept = -0.7,color="black")+
#   scale_x_continuous(trans='exp',limits = c(-0.7,1.5),breaks = c(-0.5,0.25,0,0.25,0.5,0.75,1,1.25))
# 
# fig1C <- plot(rope(mod,ci=1))+
#   theme_classic()+
#   scale_y_discrete(expand = c(0,0))+
#   scale_x_continuous(breaks = seq(-.4,1.4,by=0.2))+
#   theme(
#     axis.text.y=element_blank(),
#     legend.position = c(0.8,0.8),
#     text = element_text(size=14)
#   )+
#   ylab("Posterior distribution")+
#   xlab("Treatment effect (odds ratio)")+
#   ggsci::scale_color_npg()+
#   ggsci::scale_fill_npg()+
#   scale_x_continuous(trans='exp',limits = c(-0.7,1.5),breaks = c(-0.5,0.25,0,0.25,0.5,0.75,1,1.25))+
#   geom_vline(xintercept = -0.7,color="black")

perf <- describe_posterior(mod,ci=0.95,rope_ci = 1)
perf$CI_low[2] %>% exp() %>% round(2)
perf$CI_high[2] %>% exp() %>% round(2)

#############################
#
#############################

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
perf <- describe_posterior(mod,ci=0.95,rope_ci = 1)
df.indices <- rbind(
  df.indices,
  data.frame(mystudy,N = sum(dummy[1:2,2:3]),p = fisher.test(dummy[1:2,2:3])$p.value,pd = perf$pd[2],rope = perf$ROPE_Percentage[2])
)
}

# summary measures
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

# actual figure

df.indices.fig <- df.indices %>% 
  mutate(`Size (N)` = factor(case_when(
    (N>=0 & N<=500)~"0-500",
    (N>500 & N<=1000)~"501-1'000",
    (N>1000 & N<=10000)~"1'001-10'000",
    (N>10000)~">10'000"),
    levels = c("0-500","501-1'000","1'001-10'000",">10'000"))) %>% 
  mutate(
    p = case_when(
      p<0.0001~0.0001,
      TRUE~p
    )
  )
    
# A
fig2A <- df.indices.fig %>% 
  ggplot(aes(x=p,y=pd,color=`Size (N)`,label=mystudy))+
  geom_point(alpha=0.01,show.legend = T)+
  geom_text(show.legend = F)+
  ggsci::scale_color_nejm()+
  guides(color= guide_legend(override.aes = list(alpha = 1)), size=guide_legend())+
  theme_classic()+
  theme(
    text = element_text(size=14),
    legend.position = "bottom"
  )+
  scale_x_log10(breaks = c(0.0001,0.001,0.01,0.05,0.1,0.5,1),labels = c("0.0001","0.001","0.01","0.05","0.1","0.5","1"))+#(n.breaks = 10,limits = c(0.0000001,1))+
  scale_y_continuous(labels = scales::percent,n.breaks = 10,limits = c(0.5,1))+
  # labs(title = "Probability of Direction")+
  ylab("pd")+
  xlab("P-value")

fig2B <- df.indices.fig %>% 
  ggplot(aes(x=p,y=rope,color=`Size (N)`,label=mystudy))+
  geom_point(alpha=0.01,show.legend = T)+
  geom_text(show.legend = F)+
  ggsci::scale_color_nejm()+
  guides(color= guide_legend(override.aes = list(alpha = 1)), size=guide_legend())+
  theme_classic()+
  theme(
    text = element_text(size=14),
    legend.position = "bottom"
  )+
  scale_x_log10(breaks = c(0.0001,0.001,0.01,0.05,0.1,0.5,1),labels = c("0.0001","0.001","0.01","0.05","0.1","0.5","1"))+#(n.breaks = 10,limits = c(0.0000001,1))+
  scale_y_continuous(labels = scales::percent,n.breaks = 10,limits = c(0,1))+
  # labs(title = "Percentage in Region of Practical Equivalence")+
  ylab(bquote(ROPE[full]))+
  xlab("P-value")

fig2C <- df.indices.fig %>% 
  ggplot(aes(x=pd,y=rope,color=`Size (N)`,label=mystudy))+
  geom_point(alpha=0.01,show.legend = T)+
  geom_text(show.legend = F)+
  ggsci::scale_color_nejm()+
  guides(color= guide_legend(override.aes = list(alpha = 1)), size=guide_legend())+
  theme_classic()+
  theme(
    text = element_text(size=14),
    legend.position = "bottom"
  )+
  scale_y_continuous(labels = scales::percent,n.breaks = 10,limits = c(0,1))+
  scale_x_continuous(labels = scales::percent,n.breaks = 10,limits = c(0.5,1))+
  # labs(title = "Association")+
  ylab(bquote(ROPE[full]))+
  xlab("Probability of direction (pd)")

ggpubr::ggarrange(
  fig2A,
  fig2B,
  fig2C,
  labels = c("A","B","C"),
  nrow=3,
  common.legend = T
)

ggsave("figure2_revised.jpeg",dpi=600,height=9.5,width = 7)
ggsave("figure2_revised.pdf",dpi=600,height=9.5,width = 7)

# check P values of 12 and 15
df.indices$p[12] %>% round(3)
df.indices$p[15] %>% round(3)

