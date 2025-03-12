library(dplyr)
library(tidyr)
library(glmmTMB)
library(ggplot2)
library(lubridate)


df.pa <- readRDS("Data/PoultHen/PA.rds")
st <- "PA"

  df1 <- df.pa
  df1 <- df1[df1$Unknown/df1$Total<0.25 & df1$PHratio<=16,] #filter too many unknown sex-age and too many poults/hen
  df1 <- df1[!(df1$Hens>=8 & df1$Poults==0),] #filter too many hens without poults
  
  df1$doy.scale <- scale(df1$doy)
  df1$doy.2 <- df1$doy.scale^2
  
  
###### Statistical analysis of data  

  df.all <- df1[df1$Year>2018,]
   df.all$Year <- as.factor(df.all$Year)
   df.all$Year <- as.factor(df.all$Year)
   df.all$MU <- as.factor(df.all$MU)
   
  ### Make each hen in a record a separate observation
  df1.all <- uncount(df.all,Hens)
  
##################################################################      
### Model proportion of hens with poults using glmmTMB

### Vector-valued model https://m-clark.github.io/posts/2020-03-01-random-categorical/

m.hwb.all <- glmmTMB(HWB ~ 1 + Year + doy.scale + doy.2 + (1|MU), dispformula = ~0, 
                        data = df1.all, family = binomial())
  summary(m.hwb.all)

################################################################  
### Predict proportion with poults (1-proportion without poults)
  aug31 <- (244-attr(df1$doy.scale,"scaled:center"))/attr(df1$doy.scale,"scaled:scale")
  aug31.2 <- ((244-attr(df1$doy.scale,"scaled:center"))/attr(df1$doy.scale,"scaled:scale"))^2
  mindoy <- format((min(df1$doy,na.rm=T)-attr(df1$doy.scale,"scaled:center"))/attr(df1$doy.scale,"scaled:scale"),digits=3)
  maxdoy <- format((max(df1$doy,na.rm=T)-attr(df1$doy.scale,"scaled:center"))/attr(df1$doy.scale,"scaled:scale"),digits=3)
  



##### ALL SIGHTINGS #####################
pred.pop <- expand.grid(Year = unique(df1.all$Year), doy.scale=aug31, doy.2=aug31.2, MU=unique(df1.all$MU))
final4.all <- cbind(pred.pop,predict(m.hwb.all, newdata=pred.pop, se.fit=TRUE, type="response"))
final4.all$fit <- final4.all$fit
final4.all$lcl95 <- (final4.all$fit - 1.96*final4.all$se.fit) 
final4.all$ucl95 <- (final4.all$fit + 1.96*final4.all$se.fit)
final4.all$class <- final4.all$MU


p4.all <-  ggplot(data=final4.all) +
  geom_errorbar(aes(x=Year, ymin=lcl95, ymax=ucl95, group=MU, color=class), width=0.2, position=position_dodge(width=0.2)) +
  geom_line(aes(x=Year, y=fit, group=MU, color=class), position=position_dodge(width=0.2)) +
  geom_point(aes(x=Year, y=fit, group=MU, color=class), position=position_dodge(width=0.2)) +
  theme_classic() +
  scale_y_continuous(breaks=seq(0,1, by=.05)) +
  labs(y="Proportion of hens with poults", x="Year", color="WMU group",
       title="Proportion of hens with poults ALL SIGHTINGS") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),
        legend.position="none")
p4.all

pred.pop4a <- expand.grid(Year = unique(df1.all$Year), doy.scale = seq(mindoy,maxdoy,by=0.2), MU = unique(df1.all$MU))
pred.pop4a$doy.2 <- pred.pop4a$doy.scale^2
final4a <- cbind(pred.pop4a,predict(m.hwb.all, newdata=pred.pop4a, se.fit=TRUE, type="response"))
final4a$fit <- final4a$fit
final4a$lcl95 <- (final4a$fit - 1.96*final4a$se.fit) 
final4a$ucl95 <- (final4a$fit + 1.96*final4a$se.fit)
final4a$class <- factor(final4a$MU, levels=c(1:10))

final4a$origin <- as.Date("2019-01-01",tz = "UTC") - days(1)
final4a$doy.date <- as.Date(final4a$doy.scale*attr(df1$doy.scale,"scaled:scale")+attr(df1$doy.scale,"scaled:center"), 
                            origin = final4a$origin, tz = "UTC") 
final4a$State <- st
write.csv(final4a[final4a$Year==2019,],paste0("./doy/",st,"doy_hwb.csv"),row.names = F)

p4a.all <-  ggplot(data=final4a[final4a$Year==2019,]) +
  #  geom_errorbar(aes(x=doy, ymin=lcl95, ymax=ucl95, group=MU, color=class), width=0.2, position=position_dodge(width=0.4)) +
  geom_line(aes(x=doy.date, y=fit, group=MU, color=class)) + #, position=position_dodge(width=0.4)) +
  geom_point(aes(x=doy.date, y=fit, group=MU, color=class)) + #, position=position_dodge(width=0.4)) +
  theme_classic() +
  scale_y_continuous(breaks=seq(0,1, by=.05)) +
  scale_x_date(date_breaks = "1 week", date_labels =  "%d-%b") +
  labs(y="Proportion of hens with poults", x="Year", color="WMU group",
       title="Effect of Day of Year on Hens with Poults ALL SIGHTINGS") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12), legend.title=element_text(size=14), legend.key.width = unit(2.5, "line"))
p4a.all






