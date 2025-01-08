
library(stringr)
library(ggplot2)
library(dplyr)
#Figure 2
GridYear<-read.csv("GridYearRevision.csv")
GridYear$Site<-paste(GridYear$Grid, GridYear$Year, GridYear$Season, sep=" ")
opp<-subset(GridYear, Species=="Opossum")
racc<-subset(GridYear, Species=="Raccoon")
opp<-opp %>% arrange(Site)
racc<-racc %>% arrange(Site)
all<-bind_cols(opp, racc)
names(all)[names(all) == "mean...6"] <- "OpossumDensity"
names(all)[names(all) == "mean...17"] <- "RaccoonDensity"
names(all)[names(all) == "Grid...2"] <- "Grid"
submeany <- aggregate(all$OpossumDensity, by = list(all$Grid), mean)
submeanx <- aggregate(all$RaccoonDensity, by = list(all$Grid), mean)
alltot<-bind_cols(submeanx, submeany)
names(alltot)[names(alltot) == "x...4"] <- "Opossum"
names(alltot)[names(alltot) == "x...2"] <- "Raccoon"
names(alltot)[names(alltot) == "Group.1...1"] <- "Grid"
alltot$Habitat<- str_sub(alltot$Grid, 1,1)
alltot$Habitat[alltot$Habitat == "B"] <- "Bottomland hardwood"
alltot$Habitat[alltot$Habitat == "P"] <- "Pine upland"
alltot$Habitat[alltot$Habitat == "R"] <- "Riparian forest"
alltot$Habitat[alltot$Habitat == "W"] <- "Isolated wetland"

vars <- c("#DF536B", "#61D04F","#2297e6", "#F5C710")
textPart1 <- "paste(italic(r), \" = 0.667\")"  
textPart2 <- "paste(italic(p), \" < 0.001\")" 
petal.lm<-lm(Opossum~Raccoon, alltot)

ggplot(alltot, aes(x = Raccoon, y = Opossum, color=Habitat, label=Grid))+geom_point(size=3)+
  scale_color_manual(values=vars)+
  geom_abline(slope = coef(petal.lm)[["Raccoon"]], 
              intercept = coef(petal.lm)[["(Intercept)"]], linewidth=0.8)+geom_abline(slope = coef(petal.lm)[["Raccoon"]]+0.613, 
                                                                                      intercept = coef(petal.lm)[["(Intercept)"]], linewidth=0.75, linetype="longdash")+
  geom_abline(slope = coef(petal.lm)[["Raccoon"]]-0.613, 
              intercept = coef(petal.lm)[["(Intercept)"]], linewidth=0.75,  linetype="longdash")+
  theme(panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'), legend.text=element_text(size=10), legend.title = element_text(size=12),
        axis.text=element_text(size=14, color="black"), legend.position = "bottom",
        axis.title = element_text(size=14, vjust=3.9))+
  theme(plot.margin=unit(c(0.9,0.9,0.5,0.5), "cm"), axis.ticks.length.x =unit(.25, "cm"))+
  labs(x = expression(paste("Estimated raccoon density ","(animals/km"^2 ,")")), 
       y= expression(paste("Estimated opossum density ", " (animals/km"^2 ,")"))) + 
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12), expand=c(0,0), limits = c(0,13))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), expand=c(0,0),limits = c(0,13))+
  annotate("text", x=10, y=1.8, label= paste(textPart1), parse=TRUE, size=5)+
  annotate("text", x=10, y=0.6, label= paste(textPart2), parse=TRUE, size=5)

#Figure 3
vars <- c("#DF536B", "#2297e6", "#F5C710")
Density<-read.csv("DensityRevision.csv")
ggplot(Density, aes(x=Habitat, y=mean, color=Species))+geom_point(size=3)+geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                                                                                        data=Density, width=0.2)+
  facet_grid(~Season, scales="free_x")+
  theme(strip.background=element_rect(color="black", fill="darkslategray3"), panel.border = element_rect(colour='black', fill='NA'), 
        panel.background = element_rect(fill='gray97'), panel.grid = element_line(colour = NA), panel.spacing = unit(1, "lines"),
        legend.key = element_rect(fill='white'), legend.position="bottom", 
        legend.text=element_text(size=12), legend.title = element_text(size=12),
        axis.text=element_text(size=14), axis.title.x = element_blank(),axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=14, vjust=3.9), strip.text = element_text(size=14, face="bold"))+
  theme(plot.margin=unit(c(0.3,0.3,0.5,0.5), "cm"), axis.ticks.length.x =unit(.25, "cm"))+
  scale_y_continuous(breaks=c(2,4,6,8,10), limits = c(0,10))+
  labs(y = expression(paste("Estimated density  ", "(animals/km "^2 ,")")))+scale_color_manual(values=vars)+
  scale_x_discrete(labels=c('Bottomland\nhardwood','Upland\npine', 'Riparian\nforest', 'Isolated\nwetlands'))

