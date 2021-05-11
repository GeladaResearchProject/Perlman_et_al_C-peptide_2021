## Data
load(file="/Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021/Perlman_et_al_C-peptide_2021.RData")


## Libraries
library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)

library(ggplot2)
library(ggeffects)
library(visreg)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)
library(effectsize)
library(sjPlot)

#### ANALYSES #### -------------------------------------------------------------
## WORKING DIRECTORY # -----------------------------------
# setwd("//Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021")


# Data
# propFA <- read.csv("df_propFA.csv",header=T,sep=",", fill=T)
# feed <- read.csv("df_feed.csv",header=T,sep=",", fill=T)
# move <- read.csv("df_move.csv",header=T,sep=",", fill=T)
# cpep1 <- read.csv("df_cpep1.csv",header=T,sep=",", fill=T)
# cpep2 <- read.csv("df_cpep2.csv",header=T,sep=",", fill=T)
# df_fig1 <- read.csv("df_fig1.csv",header=T,sep=",", fill=T)
# df_weather1<-as_tibble(read.csv("//Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021/df_weather1.csv",header=TRUE,stringsAsFactors=FALSE))
# df_weather2<-as_tibble(read.csv("//Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021/df_weather2.csv",header=TRUE,stringsAsFactors=FALSE))
# df_season<-as_tibble(read.csv("//Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021/df_season.csv",header=TRUE,stringsAsFactors=FALSE))

# To save all objects to RData file (done already)
# save.image(file="Perlman_et_al_C-peptide_2021.RData") 

#### ANALYSES 1. SEASONAL VARIATION IN ABOVE-GROUND FEEDING --------------------

# Feeding represented the majority of the gelada diet.
mean(feed$PropFeed)
min(feed$PropFeed)
max(feed$PropFeed)

# Above-ground feeding accounted for the largest percent of the feeding budget.
propFA$PropFA <- propFA$FA/propFA$F.Total
df.FA <- propFA %>% na.omit() %>% group_by(month, season) %>% summarise(mean=mean(PropFA))
summary(df.FA)

# Above-ground feeding was positively associated with cumulative rainfall.
model1 <- glmer.nb(FA ~ scale(Rain30) + MinT + offset(log(F.Total)) + (1|ID) + (1|Unit), data=propFA)
summary(model1)

# 4 driest months:
df.dry <- df.FA %>% filter(season == "Dry")
summary(df.dry)

# 4 wettest months:
df.wet <- df.FA %>% filter(season == "Wet")
summary(df.wet)


#### ANALYSES 2. SEASONAL VARIATION IN FORAGING EFFORT AND ENERGETIC STATUS ----

# Time spent feeding is negatively associated with rainfall (not significant).
model2 <- glmer.nb(sum ~ Rain30 + MinT + offset(log(total.obs)) + (1|ID) + (1 |Unit), data=feed)
summary(model2)


# Time spent moving is negatively associated with rainfall (significant).
model3 <- glmer.nb(sum ~ Rain30 + MinT + offset(log(total.obs)) + (1|ID) + (1|Unit), data=move)
summary(model3)


# UCPs are negatively associated with rainfall (significant).
cpep1$Time <- strptime(cpep1$Time, format="%I:%M")
cpep1$Time <- as.POSIXct(cpep1$Time)
model4 <- glmer(UCP ~  scale(Time)*scale(Rain30) + MinT + (1|ID), family = Gamma(link = "log"), data=cpep1)
summary(model4)



# UCPs are negatively associated with feeding above ground (significant).
cpep2$Time <- strptime(cpep2$Time, format="%I:%M")
cpep2$Time <- as.POSIXct(cpep2$Time)
model5 <- lmer(UCP ~  scale(Time)*scale(PropFA) + MinT + (1|ID), data=cpep2)
summary(model5)

#### FIGURES #### -------------------------------------------------------------


##### FIGURE 1 | Rainfall and feeding behavior ####   

# Custom labels
str_pad_custom <- function(season){
  new_labels <- stringr::str_pad(season,15, "right")
  return(new_labels)}

# Create Fig 1
fig1 <- ggplot(subset(df_fig1, !is.na(season)), aes(x=season, y=sum/F.Total*100, fill=Activity3)) +
  geom_boxplot(position=position_dodge(0.5), width=0.4, colour="black", lwd = 0.5, outlier.shape = NA, key_glyph=rectangle_key_glyph(color="black", size=0.4, padding=margin(1,1,1,1)))

# Format Fig 1
fig1 <- fig1 +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_rect(size = 1, colour = "black")) + 
  theme(axis.text = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text( size = 10)) +
  theme(plot.margin = unit(c(2.25,1.25,2.25,1.25), "lines")) +
  theme(axis.ticks.x = element_blank(), axis.ticks.y = element_line(size =0.5), axis.ticks.length = unit(.1, "cm")) +
  scale_fill_manual(values=c("yellowgreen", "darkgreen"), labels = str_pad_custom) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color="black", size =10, vjust=1),
        axis.title.y = element_text(vjust = 1.25, angle = 90, size = 10), axis.text.y = element_text(color="black", size=8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.1), shape=16, size=0.5, colour="black") +
  ylab("Percent of feeding time") + scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.title = element_blank(), legend.text=element_text(size=8)) + theme(legend.key.size=unit(0.3, "cm")) +
  theme(legend.direction = "horizontal") + theme(legend.position = c(0.5, 1)) + theme(legend.box.margin=margin(c(0,0,15,0))) +
  theme(plot.background = element_blank(), legend.background = element_blank(),legend.box.background = element_blank(), legend.spacing.x = unit(0.1, 'cm')) +
  guides(fill = guide_legend(override.aes = list(shape = 0, size = 1, colour=NA))) +
  theme(axis.title = element_text(family = "Arial"))

fig1

ggsave(fig1, file="Fig1.png", width =8, height = 8, units = c("cm"))


dev.off()

##### FIGURE 2ab | Rainfall and activity budgets ####

# Models for visualizations
modelFEED <- glmer.nb(sum ~ Rain30 + MinT + (1|ID) + (1|Unit), data=feed)
modelMOVE <- glmer.nb(sum ~ Rain30 + MinT + (1|ID) + (1|Unit), data=move)

# Create Fig 2a
fig2a <- visreg(modelFEED,"Rain30",type="contrast",scale="linear", line=list(col=c("black")), ylim=c(-4,2), gg=T,
                points=list(pch=21, bg="olivedrab3", col="black"), xlab="Cumulative rainfall (mm)", ylab="Feeding time residuals",band=F) 

# Format Fig 2a
fig2a <- fig2a + theme_bw() + geom_point(pch=21, size = 1.5, bg="olivedrab3", col="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_rect(size = 1, colour = "black")) +
  theme(axis.ticks=element_line(colour="black", size = 0.5), axis.ticks.length = unit(.15, "cm")) +
  theme(plot.margin = unit(c(1.25,1.25,1.25,1.25), "lines")) +
  theme(axis.title.x = element_text(size = 12, vjust = -1), axis.text.x = element_text(color="black", size = 12), axis.title.y = element_text(vjust = 2.5, angle = 90, size = 12),
        axis.text.y = element_text(color="black", size = 12)) + 
  scale_y_continuous(limits=c(-4,2))

fig2a

# Create Fig 2b
fig2b <- visreg(modelMOVE,"Rain30",type="contrast",scale="linear", line=list(col=c("black")), ylim=c(-4,2), gg=T,
                points=list(pch=21, bg="darkgoldenrod1", col="black"), xlab="Cumulative rainfall (mm)", ylab="Moving time residuals",band=F) 

# Format Fig 2b
fig2b <- fig2b + theme_bw() + geom_point(pch=21, size = 1.5, bg="darkgoldenrod1", col="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_rect(size = 1, colour = "black")) +
  theme(axis.ticks=element_line(colour="black", size = 0.5), axis.ticks.length = unit(.15, "cm")) +
  theme(plot.margin = unit(c(1.25,1.25,1.25,1.25), "lines")) +
  theme(axis.title.x = element_text(size = 12, vjust = -1), axis.text.x = element_text(color="black", size = 12), axis.title.y = element_text(vjust = 2.5, angle = 90, size = 12),
        axis.text.y = element_text(color="black", size = 12)) + 
  scale_y_continuous(limits=c(-4,2))

fig2b

## Fig 2ab final
fig2ab <- ggarrange(fig2a, fig2b, ncol=2, nrow=1, legend = "top", labels = c("A", "B"), vjust = 1.5, hjust = -0.2, font.label = list(size = 16, color = "black", face = "bold"))
fig2ab


ggsave("fig2ab.png", width = 16, height = 8, units = "cm")


dev.off()


##### FIGURE 3 | UCP concentrations ####
cpep1$scaleTime <- scale(cpep1$Time)
cpep1$scaleRain30 <- scale(cpep1$Rain30)
model4 <- glmer(UCP ~  scaleTime*scaleRain30 + MinT + (1|ID), family = Gamma(link = "log"), data=cpep1)


# Fig 3a
set_theme(
  base = theme_bw(), axis.textcolor="black")

fig3a <- plot_model(model4, transform=NULL,type="std",  
                   show.values=TRUE, value.offset=0.3,
                   vline.color="gray",
                   order.terms = c( 2, 3, 4, 1),
                   show.p=FALSE,
                   title="", colors=c("slateblue4", "black"))
fig3a <- fig3a +
  font_size(axis_title.x = 12, labels.x=12, labels.y=12) +
  scale_x_discrete(labels = c('Time','Time*Rainfall','Temperature', 'Rainfall')) +
  theme(axis.title.x = element_text(vjust = -1, color="black")) +
  panel_border(size=1, color="black") +   theme(plot.margin = unit(c(1,1,1,1), "lines"))

fig3a

# Fig 3b
modelUCP <- glmer(UCP ~  Rain30 + MinT + (1 | ID), family = Gamma(link = "log"), data=cpep1)

fig3b <- visreg(modelUCP,"Rain30",type="contrast",scale="linear", line=list(col=c("black")), ylim=c(-3,3), gg=T,
                points=list(pch=21, bg="slateblue3", col="black"), xlab="Cumulative rainfall (mm)", ylab="Urinary C-peptide residuals",band=F) 

fig3b <- fig3b + theme_bw() + geom_point(pch=21, size = 1.25, bg="slateblue3", col="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_rect(size = 1, colour = "black")) +
  theme(axis.ticks=element_line(colour="black", size = 0.5), axis.ticks.length = unit(.15, "cm")) +
  theme(axis.title.x = element_text(size = 12, vjust = -1), axis.text.x = element_text(color="black", size = 12), axis.title.y = element_text(vjust = 2.5, angle = 90, size = 12),
        axis.text.y = element_text(color="black", size = 12)) + scale_x_continuous(limits=c(0,500)) +
  theme(plot.margin = unit(c(1,1,1,1.5), "lines"))


fig3b


## Fig 3ab final
fig3ab <- ggarrange(fig3a, fig3b, ncol=2, nrow=1, 
                    legend = "top", labels = c("A", "B"), 
                    vjust = 1.5, hjust = -0.2, 
                    font.label = list(size = 16, color = "black", face = "bold"),
                    widths = c(0.45, 0.4), heights = c(1.1,1.1))
fig3ab

ggsave("fig3ab.png", width = 20, height = 8, units = "cm")

dev.off()

##### FIGURE 4 | Seasonality of Rainfall, Temperature, Foraging effort, and UCP concentrations####


## Rainfall and Feeding-aboveground
df_weather1<-as_tibble(read.csv("//Users/Rachel/Desktop/CURRENT/Perlman_et_al_C-peptide_2021/df_weather1.csv",header=TRUE,stringsAsFactors=FALSE))
df_season<-as_tibble(read.csv("//Users/Rachel/Desktop/CURRENT/Perlman_et_al_C-peptide_2021/df_season.csv",header=TRUE,stringsAsFactors=FALSE))

df_weather1 <- na.omit(df_weather1)

df_weather1 <- df_weather1 %>% filter(panel != "Total takeovers")
df_weather1 <- df_weather1 %>% filter(panel != "Time spent feeding")
df_weather1 <- df_weather1 %>% filter(panel != "Time spent moving")
df_weather1 <- df_weather1 %>% filter(panel != "Urinary C-peptide")
df_weather1 <- df_weather1 %>% filter(panel != "FA")


fig_4a<-ggplot(data = df_weather1, mapping = aes(x = as.factor(x), y = y, fill=as.factor(panel))) +
  geom_rect(data = df_season, aes(xmin = x-0.5, xmax = x+0.5, ymin = -Inf, ymax = Inf, fill=panel), fill = c("#abd9e9","#abd9e9","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#abd9e9","#abd9e9","#abd9e9"), alpha = 0.5) +
  scale_x_discrete(labels = c("S","O","N","D","J","F","M","A","M","J","J","A")) +
  scale_y_continuous(sec.axis = sec_axis(~./5, name =""), limits = c(0,100), labels = scales::number_format(accuracy = 0.1)) + 
  theme(axis.title.y.left=element_text(vjust=0.37)) +  
  theme(axis.title.x=element_blank()) +
  xlab("Months") +
  ylab("% of total feeding time") + 
  theme(legend.title=element_blank(), legend.position = c(0.48,0.9), legend.background = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0)) +
  geom_text(aes(x = 6.25, y = 100, label = "DRY", vjust=-2), size=6) +
  geom_text(aes(x = 11, y = 100, label = "WET", vjust=-2), size=6) +
  geom_text(aes(x = 1.5, y = 100, label = "WET", vjust=-2), size=6) +
  coord_cartesian(clip = "off") +
  theme(axis.title = element_text(size = 20), axis.text=element_text(size=16)) +
  geom_bar(data = subset(df_weather1, panel %in% "Rainfall"), stat = "identity", aes(y=y*5, group = as.factor(panel), shape=as.factor(panel)), fill = "white", colour="black") +
  geom_errorbar(data = subset(df_weather1, panel %in% "Rainfall"), aes(ymin=y*5-se, ymax=y*5+se, group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  geom_point(data = subset(df_weather1, panel %in% "Feeding above-ground"), aes(group = as.factor(panel), shape=as.factor(panel)), colour = "black", size=3) +
  geom_line(data = subset(df_weather1, panel %in% "Feeding above-ground"), aes(group=panel), position=position_dodge(0.5), colour="black", linetype=2, size=1) +
  geom_errorbar(data = subset(df_weather1, panel %in% "Feeding above-ground"), aes(ymin=y-se, ymax=y+se, group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  scale_shape_manual(labels = c("Feeding above-ground", "Rainfall"), values = c(16,22), name="", guide="legend") +
  scale_fill_manual(labels = c("Feeding above-ground", "Rainfall"), values = c("black","white"), name = "", guide=FALSE) +
  theme(legend.key=element_rect(fill=NA)) +
  guides(shape = guide_legend(override.aes = list(fill = c("black","white"), colour = c("black","black"), size=c(3,6)), nrow=2)) +
  theme(plot.margin=unit(c(1.5,1,1,2.2),"cm")) +
  theme(legend.text=element_text(size=16)) +
  labs(tag = "Total daily rain (mm)") +
  theme(plot.tag.position = c(1,0.55), plot.tag=element_text(size=20,angle=270))

fig_4a


dev.off()


## Seasonal - Time spent feeding
df_weather1<-as_tibble(read.csv("//Users/Rachel/Desktop/CURRENT/Perlman_et_al_C-peptide_2021/df_weather1.csv",header=TRUE,stringsAsFactors=FALSE))

df_weather1 <- na.omit(df_weather1)

df_weather1 <- df_weather1 %>% filter(panel != "Total takeovers")
df_weather1 <- df_weather1 %>% filter(panel != "Time spent moving")
df_weather1 <- df_weather1 %>% filter(panel != "Urinary C-peptide")
df_weather1 <- df_weather1 %>% filter(panel != "FA")
df_weather1 <- df_weather1 %>% filter(panel != "Feeding above-ground")

df_weather1 <- df_weather1 %>% filter(panel != "Rainfall")

theme_update(plot.title = element_text(hjust = 0.5, vjust=1))

fig_4b<-ggplot(data = df_weather1, mapping = aes(x = as.factor(x), y = y, fill=as.factor(panel))) +
  geom_rect(data = df_season, aes(xmin = x-0.5, xmax = x+0.5, ymin = -Inf, ymax = Inf, fill=panel), fill = c("#abd9e9","#abd9e9","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#abd9e9","#abd9e9","#abd9e9"), alpha = 0.5) +
  scale_x_discrete(labels = c("S","O","N","D","J","F","M","A","M","J","J","A")) +
  scale_y_continuous(limits=c(55,77), breaks=c(55,60,65,70,75), labels = scales::number_format(accuracy = 0.1)) +
  theme(axis.title.y=element_text(vjust=2)) +
  theme(axis.title.x=element_blank()) +
  xlab("Months") +
  ylab("% of activity budget") +
  theme(axis.title = element_text(size = 20), axis.text=element_text(size=16)) +
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0)) +
  geom_point(data = subset(df_weather1, panel %in% "Time spent feeding"), aes(y=y*100, group = as.factor(panel), shape=as.factor(panel)), colour = "black", size=3) +
  geom_line(data = subset(df_weather1, panel %in% "Time spent feeding"), aes(y=y*100, group=panel), position=position_dodge(0.5), colour="black", linetype=2, size=1) +
  geom_errorbar(data = subset(df_weather1, panel %in% "Time spent feeding"), aes(ymin=(y-se)*100, ymax=(y+se)*100, group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  scale_shape_manual(labels ="Time spent feeding", values = 16, name="", guide="legend") +
  scale_fill_manual(labels ="Time spent feeding", values = "black", name = "", guide=FALSE) +
  theme(plot.margin=unit(c(0.2,2.5,1,2.5),"cm"))
fig_4b

dev.off()

## Seasonal - Time spent moving
df_weather1<-as_tibble(read.csv("//Users/Rachel/Desktop/CURRENT/Perlman_et_al_C-peptide_2021/df_weather1.csv",header=TRUE,stringsAsFactors=FALSE))

df_weather1 <- na.omit(df_weather1)

df_weather1 <- df_weather1 %>% filter(panel != "Total takeovers")
df_weather1 <- df_weather1 %>% filter(panel != "Time spent feeding")
df_weather1 <- df_weather1 %>% filter(panel != "Urinary C-peptide")
df_weather1 <- df_weather1 %>% filter(panel != "FA")
df_weather1 <- df_weather1 %>% filter(panel != "Feeding above-ground")
df_weather1 <- df_weather1 %>% filter(panel != "Rainfall")

fig_4c<-ggplot(data = df_weather1, mapping = aes(x = as.factor(x), y = y, fill=as.factor(panel))) +
  geom_rect(data = df_season, aes(xmin = x-0.5, xmax = x+0.5, ymin = -Inf, ymax = Inf, fill=panel), fill = c("#abd9e9","#abd9e9","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#abd9e9","#abd9e9","#abd9e9"), alpha = 0.5) +
  scale_x_discrete(labels = c("S","O","N","D","J","F","M","A","M","J","J","A")) +
  scale_y_continuous(limits=c(9,20), breaks=c(10,12,14,16,18,20), labels = scales::number_format(accuracy = 0.1)) +
  theme(axis.title.y=element_text(vjust=2)) +
  theme(axis.title.x=element_blank()) +
  xlab("Months") +
  ylab("% of activity budget") +
  theme(axis.title = element_text(size = 20), axis.text=element_text(size=16)) + 
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0)) +
  geom_point(data = subset(df_weather1, panel %in% "Time spent moving"), aes(y=y*100, group = as.factor(panel), shape=as.factor(panel)), colour = "black", size=3) +
  geom_line(data = subset(df_weather1, panel %in% "Time spent moving"), aes(y=y*100, group=panel), position=position_dodge(0.5), colour="black", linetype=2, size=1) +
  geom_errorbar(data = subset(df_weather1, panel %in% "Time spent moving"), aes(ymin=(y-se)*100, ymax=(y+se)*100, group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  scale_shape_manual(labels ="Time spent moving", values = 16, name="", guide="legend") +
  scale_fill_manual(labels ="Time spent moving", values = "black", name = "", guide=FALSE) +
  theme(plot.margin=unit(c(0.2,2.5,1,2.5),"cm"))
fig_4c

dev.off()

## Seasonal - C-peptide
df_weather1<-as_tibble(read.csv("//Users/Rachel/Desktop/CURRENT/Perlman_et_al_C-peptide_2021/df_weather1.csv",header=TRUE,stringsAsFactors=FALSE))

df_weather1 <- na.omit(df_weather1)

df_weather1 <- df_weather1 %>% filter(panel != "Total takeovers")
df_weather1 <- df_weather1 %>% filter(panel != "Time spent feeding")
df_weather1 <- df_weather1 %>% filter(panel != "Time spent moving")
df_weather1 <- df_weather1 %>% filter(panel != "FA")
df_weather1 <- df_weather1 %>% filter(panel != "Feeding above-ground")
df_weather1 <- df_weather1 %>% filter(panel != "Rainfall")

fig_4d<-ggplot(data = df_weather1, mapping = aes(x = as.factor(x), y = y, fill=as.factor(panel))) +
  geom_rect(data = df_season, aes(xmin = x-0.5, xmax = x+0.5, ymin = -Inf, ymax = Inf, fill=panel), fill = c("#abd9e9","#abd9e9","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#abd9e9","#abd9e9","#abd9e9"), alpha = 0.5) +
  scale_x_discrete(labels = c("S","O","N","D","J","F","M","A","M","J","J","A")) +
  scale_y_continuous(limits=c(0,3.5), labels = scales::number_format(accuracy = 0.1)) +
  theme(axis.title.y=element_text(vjust=3.8)) +
  theme(axis.title.x = element_text(vjust=-0.5)) +
  xlab("Months") +
  ylab("Urinary C-peptide") +
  theme(axis.title = element_text(size = 20), axis.text=element_text(size=16)) +
  theme(legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0)) +
  geom_point(data = subset(df_weather1, panel %in% "Urinary C-peptide"), aes(y=y, group = as.factor(panel), shape=as.factor(panel)), colour = "black", size=3) +
  geom_line(data = subset(df_weather1, panel %in% "Urinary C-peptide"), aes(y=y, group=panel), position=position_dodge(0.5), colour="black", linetype=2, size=1) +
  geom_errorbar(data = subset(df_weather1, panel %in% "Urinary C-peptide"), aes(ymin=y-se, ymax=y+se, group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  scale_shape_manual(labels ="Urinary C-peptide", values = 16, name="", guide="legend") +
  scale_fill_manual(labels ="Urinary C-peptide", values = "black", name = "", guide=FALSE) +
  theme(plot.margin=unit(c(0.2,2.5,0.5,2.78),"cm"))

fig_4d

dev.off()

# Fig 4 final
fig4 <- grid.arrange(arrangeGrob(fig_4a), 
                     arrangeGrob(fig_4b, top=textGrob("Time spent feeding", gp=gpar(fontsize=20), hjust=0.3)), 
                     arrangeGrob(fig_4c, top=textGrob("Time spent moving", gp=gpar(fontsize=20), hjust=0.3)), 
                     arrangeGrob(fig_4d, top=textGrob("Urinary C-peptide", gp=gpar(fontsize=20), hjust=0.3)), 
                     ncol=1,  heights = c(1/4, 1/4, 1/4, 1/4))


ggsave("//Users/Rachel/Desktop/Fig4.png", plot=fig4, width = 8, height = 20, units = "in")


dev.off()

#### SUPPLEMENTAL DATA #### -------------------------------------------------------------




# Fig S1 | Ecological seasonality at Simien Mountains National Park ####

# SMGRP data
df_season<-as_tibble(read.csv("//Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021/df_season.csv",header=TRUE,stringsAsFactors=FALSE))
df_weather2<-as_tibble(read.csv("//Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021/df_weather2.csv",header=TRUE,stringsAsFactors=FALSE))

df_weather2 <- df_weather2 %>% filter(panel != "Maximum Temperature")
df_weather2 <- df_weather2 %>% filter(panel != "Total takeovers")
df_weather2 <- df_weather2 %>% filter(panel != "Daily Rain")
df_weather2 <- df_weather2 %>% filter(panel != "Minimum Temperature Study Period")
df_weather2 <- df_weather2 %>% filter(panel != "Rainfall Study Period")

# Fig S1a
fig_S1a<-ggplot(data = df_weather2, mapping = aes(x = as.factor(x), y = y, fill=as.factor(panel))) +
  geom_rect(data = df_season, aes(xmin = x-0.5, xmax = x+0.5, ymin = -Inf, ymax = Inf, fill=panel), fill = c("#abd9e9","#abd9e9","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#abd9e9","#abd9e9","#abd9e9"), alpha = 0.5) +
  scale_x_discrete(labels = c("S","O","N","D","J","F","M","A","M","J","J","A")) +
  scale_y_continuous(sec.axis = sec_axis(~.*55, name ="", labels=NULL, breaks=c(0,100,200,300,400,500)), limits = c(0,10.3), labels = scales::number_format(accuracy = 0.1)) + 
  xlab("Months") +
  ggtitle("(A) SMGRP 2006-2020") +
  ylab("Minimum temperature (°C)") + 
  theme(axis.title = element_text(size = 22), axis.text=element_text(size=16)) +
  theme(axis.title.y.left=element_text(vjust=0.37), axis.title.x=element_text(hjust=1.2, vjust=0.5), axis.ticks.length.y.right=unit(0,"cm"))  +
  theme(legend.title=element_blank(), legend.position = c(0.5,1.05), legend.background = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5, size=20)) +
  geom_bar(data = subset(df_weather2, panel %in% "Cumulative rainfall"), stat = "identity", aes(y=y/55, group = as.factor(panel), shape=as.factor(panel)), fill = "white", colour="black") +
  geom_errorbar(data = subset(df_weather2, panel %in% "Cumulative rainfall"), aes(ymin=y/55-(se/55), ymax=y/55+(se/55), group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  
  geom_point(data = subset(df_weather2, panel %in% "Minimum temperature"), aes(group = as.factor(panel), shape=as.factor(panel)), colour = "black", fill="black", size=3) +
  geom_line(data = subset(df_weather2, panel %in% "Minimum temperature"), aes(group=panel), position=position_dodge(0.5), colour="black", linetype=2, size=1) +
  geom_errorbar(data = subset(df_weather2, panel %in% "Minimum temperature"), aes(ymin=y-se, ymax=y+se, group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  
  scale_shape_manual(labels = c("Minimum temperature", "Cumulative rainfall"), values = c(22,22), name="", guide="legend") +
  scale_fill_manual(labels = c("Minimum temperature", "Cumulative rainfall"), values = c("black","white"), name = "", guide=FALSE) +
  guides(shape = guide_legend(override.aes = list(fill = c("black","white"), colour = c("black","black"), size=4), nrow=1)) +
  theme(plot.margin=unit(c(1,0.2,1,1),"cm")) +
  theme(legend.text=element_text(size=22))
fig_S1a



# Study period
df_weather2<-as_tibble(read.csv("//Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021/df_weather2.csv",header=TRUE,stringsAsFactors=FALSE))

df_weather2 <- df_weather2 %>% filter(panel != "Maximum Temperature")
df_weather2 <- df_weather2 %>% filter(panel != "Minimum temperature")
df_weather2 <- df_weather2 %>% filter(panel != "Total takeovers")
df_weather2 <- df_weather2 %>% filter(panel != "Daily Rain")
df_weather2 <- df_weather2 %>% filter(panel != "Cumulative rainfall")


# Fig S1b
fig_S1b<-ggplot(data = df_weather2, mapping = aes(x = as.factor(x), y = y, fill=as.factor(panel))) +
  geom_rect(data = df_season, aes(xmin = x-0.5, xmax = x+0.5, ymin = -Inf, ymax = Inf, fill=panel), fill = c("#abd9e9","#abd9e9","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#fdae61","#abd9e9","#abd9e9","#abd9e9"), alpha = 0.5) +
  scale_x_discrete(labels = c("S","O","N","D","J","F","M","A","M","J","J","A")) +
  scale_y_continuous(sec.axis = sec_axis(~.*55, name ="Cumulative rainfall (mm)", breaks=c(0,100,200,300,400,500)), limits = c(0,10.3), labels = NULL) + 
  xlab("Months") +
  ggtitle("(B) Study Period 2017-2018") +
  theme(axis.title = element_text(size = 22, color="white"), axis.text=element_text(size=16)) +
  theme(legend.title=element_blank(), legend.position = c(0.5,1.05), legend.background = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5, size=20), axis.title.y.left=element_text(size=0), axis.ticks.length.y.left=unit(0,"cm")) +
  geom_bar(data = subset(df_weather2, panel %in% "Rainfall Study Period"), stat = "identity", aes(y=y/55, group = as.factor(panel), shape=as.factor(panel)), fill = "white", colour="black") +
  geom_errorbar(data = subset(df_weather2, panel %in% "Rainfall Study Period"), aes(ymin=y/55-(se/55), ymax=y/55+(se/55), group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  
  geom_point(data = subset(df_weather2, panel %in% "Minimum Temperature Study Period"), aes(group = as.factor(panel), shape=as.factor(panel)), colour = "black", fill="black", size=3) +
  geom_line(data = subset(df_weather2, panel %in% "Minimum Temperature Study Period"), aes(group=panel), position=position_dodge(0.5), colour="black", linetype=2, size=1) +
  geom_errorbar(data = subset(df_weather2, panel %in% "Minimum Temperature Study Period"), aes(ymin=y-se, ymax=y+se, group=panel), width=.1, position=position_dodge(0.5), colour="black") +
  
  scale_shape_manual(labels = c("Minimum Temperature", "Cumulative rainfall"), values = c(22,22), name="", guide="legend") +
  scale_fill_manual(labels = c("Minimum Temperature", "Cumulative rainfall"), values = c("black","white"), name = "", guide=FALSE) +
  guides(shape = guide_legend(override.aes = list(fill = c("black","white"), colour = c("black","black"), size=4), nrow=1)) +
  theme(plot.margin=unit(c(1,1,1,0.2),"cm")) +
  theme(legend.text=element_text(size=22)) +
  labs(tag = "Cumulative rainfall (mm)") +
  theme(plot.tag.position = c(1,0.55), plot.tag=element_text(size=20,angle=270))

fig_S1b


# Fig S1
figS1 <- ggarrange(fig_S1a, fig_S1b, ncol=2, nrow=1, 
                   legend = "top", labels = c("", ""), vjust = 1.5, hjust = -0.2, 
                   font.label = list(size = 16, color = "black", face = "bold"),
                   align="h", common.legend = TRUE)
figS1

ggsave(figS1, file="figs1.png", width =32, height = 16, units = c("cm"))


# Fig S2 | Urinary C-peptide parallelism ####
parallelism <- as_tibble(read.csv("//Users/Rachel/Desktop/Perlman_et_al_C-peptide_2021/df_cpep_parallelism.csv",header=TRUE,stringsAsFactors=FALSE))
parallelism$Binding <- parallelism$Binding*100

# ANCOVA
ancova <- aov(Binding ~ Concentration*Label, data=parallelism)
summary(ancova)

# Fig S2
parallelism$Label <- factor(parallelism$Label, levels = c("Standard", "Gelada"))

figS2 <- ggplot(parallelism, aes(Concentration, Binding, by=Label)) + geom_point(aes(shape=Label, size=Label), stroke=1) + 
  scale_x_continuous(trans='log10', name="Concentration") + 
  scale_y_continuous(name ="% Binding", limits = c(0,100)) +
  scale_shape_manual(values=c(16, 2)) + scale_size_manual(values=c(4,3)) +
  stat_smooth(method="lm", se=F, aes(linetype=Label), color="black") +
  theme_classic(base_line_size = 1) + 
  theme(plot.margin = unit(c(0.5,1,1,1), "cm")) + 
  theme(axis.title.x = element_text(vjust=-1, size=16), 
        axis.text.x = element_text(size=16, color="black"),
        axis.title.y = element_text(vjust = 1.5, angle = 90, size = 16), 
        axis.text.y = element_text(size=16, color="black")) +
  theme(axis.ticks.length = unit(0.1,"cm")) + 
  theme(legend.position="top", legend.title=element_blank(), legend.text=element_text(size=16))

figS2

ggsave(figS2, file="c-pep parallelism.png", width =16, height = 16, units = c("cm"))


# FigS3 | Distribution of UCP data ####

# Fig S3
figS3 <- ggplot(cpep2, aes(x=UCP)) + geom_histogram(color="black", fill="gray", binwidth=6, lwd=0.5, boundary=0)

figS3 <- figS3 +   theme_classic() + 
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
  theme(axis.title.x = element_text(size=10), 
        axis.text.x = element_text(size=8, color="black"),
        axis.title.y = element_text(vjust = 1.5, angle = 90, size = 10), 
        axis.text.y = element_text(size=8, color="black")) +
  scale_y_continuous(name ="Number of samples", breaks = seq(from = 0, to = 150, by = 25), limits=c(0,152.5), expand=c(0,0)) +
  scale_x_continuous(name ="Urinary C-peptide values", expand=c(0,0)) +
  theme(axis.ticks.length = unit(0.1,"cm"))

figS3

ggsave(figS3, file="c-pep distribution.png", width =8, height = 8, units = c("cm"))

# FigS4 | Variation in UCP across time of sample collection and season ####

# Fig S4
cpep1$Time <- strptime(cpep1$Time, format="%I:%M")
cpep1$Time <- as.POSIXct(cpep1$Time)

cpep1$season <- NA
cpep1$season[cpep1$Month == "Dec_2017"] <- "Dry season"
cpep1$season[cpep1$Month == "Jan_2018"] <- "Dry season"
cpep1$season[cpep1$Month == "Feb_2018"] <- "Dry season"
cpep1$season[cpep1$Month == "Mar_2018"] <- "Dry season"

cpep1$season[cpep1$Month == "May_2018"] <- "Wet season"
cpep1$season[cpep1$Month == "Jun_2018"] <- "Wet season"
cpep1$season[cpep1$Month == "Jul_2018"] <- "Wet season"
cpep1$season[cpep1$Month == "Aug_2018"] <- "Wet season"



figS4 <- ggplot(data=subset(cpep1, !is.na(season)), aes(x=Time, y=UCP, color=season, fill=season)) + 
  theme_minimal() + theme(legend.position="none") +
  geom_point(shape=21, color="black") + facet_wrap(~ season, ncol = 1)  + 
  theme(strip.background = element_rect(fill=NA, color=NA), strip.text=element_text(size=24)) +
  ylab("Urinary C-peptide") + 
  scale_colour_manual(values = c("black","black")) +
  scale_fill_manual(values=c("#fdae61","#abd9e9")) +
  geom_smooth(method="loess") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  theme(axis.title.x = element_text(vjust=-1, size=24), 
        axis.text.x = element_text(size=24, color="black"),
        axis.title.y = element_text(vjust = 1.5, angle = 90, size = 24), 
        axis.text.y = element_text(size=24, color="black"))


figS4

ggsave(figS4, file="c-pep time of sample collection.png", width =16, height = 24, units = c("cm"))

# FigS5 | UCP concentrations are negatively associated with mean minimum temperature ####

# Fig S5
modelUCP <- glmer(UCP ~  Rain30 + MinT + (1 | ID), family = Gamma(link = "log"), data=cpep1)

figS5 <- visreg(modelUCP,"MinT",type="contrast",scale="linear", line=list(col=c("black")), ylim=c(-3,3), gg=T,
                points=list(pch=21, bg="slateblue3", col="black"), xlab="Minimum temperature (mm)", ylab="Urinary C-peptide residuals",band=F) 

figS5 <- figS5 + theme_bw() + geom_point(pch=21, size = 1.25, bg="lightblue3", col="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_rect(size = 1, colour = "black")) +
  theme(axis.ticks=element_line(colour="black", size = 0.5), axis.ticks.length = unit(.15, "cm")) +
  theme(axis.title.x = element_text(size = 12, vjust = -1), axis.text.x = element_text(color="black", size = 12), axis.title.y = element_text(vjust = 2.5, angle = 90, size = 12),
        axis.text.y = element_text(color="black", size = 12)) + scale_x_continuous(limits=c(6,10.5)) +
  theme(plot.margin = unit(c(1,1,1,1.5), "lines"))


figS5

ggsave(figS5, file="FigS5.png", width =8, height = 8, units = c("cm"))

