


# load packages
library(here)
library(tidyverse)
library(colorspace)
library(purrr)

# import data

db <- read.table(paste(here(),"submission","Database_Ventotene_2014.txt",sep="/"),
                 sep="|",header = TRUE, stringsAsFactors = TRUE)

str(db)

#########################################################################
# SiteID:    name of the sampling sites
# Lat:     latitude, gps position of trap
# Long:    longitude, gps position of trap
# week:    week of collection
# Species: mosquito species trappes, Aedes (albopictus) / Culex (pipiens)
# value:   number of mosquito collected
#########################################################################



# importing island boundary 
library(sp)
library(rgdal)
# reference systems
crsll  = CRS("+proj=longlat +datum=WGS84") # CRS lat long point
crsp   = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
crsutm = CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs") # CRS procida

# import region Lazion coastline
lazio <- readOGR(dsn="dati/Linea_costa_WGS84-33_12.shp")
proj4string(lazio) <- crsutm

utm.trap <- spTransform(SpatialPoints(list(
  cbind(db$Long,db$Lat)), 
  proj4string=crsll),crsutm)


ggplot(db,aes(Long,Lat))+geom_point(col="white")+
  geom_text(aes(Long,Lat,label=SiteID))+theme_bw()


# reordening site ID along natural to urban gradient
db$IDsite <- factor(db$SiteID,
                    levels = unique(db$SiteID[order(utm.trap@coords[,1],decreasing = FALSE)]),
                    labels = 1:10)

table(db$IDsite,db$SiteID)


ggplot(db,aes(Long,Lat))+geom_point(col="white")+
  geom_text(aes(Long,Lat,label=IDsite))+theme_bw()

library(rgeos)
bbx <- readWKT("POLYGON((365000 4516000, 368000 4516000, 368000 4518600, 365000 4516000))") 
proj4string(bbx) <- crsutm

# cut island from region coastline file
vent <- gIntersection(lazio, bbx)

plot(lazio, xlim =c(365700,367500),ylim=c(4516000.0 ,4518200) )
plot(bbx,add=TRUE)
plot(vent)


# compute human presence nearby traps

# obtain openstreetmap data
# library(osmar)
# src <- osmsource_api(url = "https://api.openstreetmap.org/api/0.6/")
# 
# bb <- corner_bbox(13.3949,40.7828,13.4485,40.8060)
# ua <- get_osm(bb, source = src)
# 
# 
# ll_poly <- as_sp(ua, "polygons")
# spplot(ll_poly, c("version"))
# 
# # identify buildings
# bg_ids <- find(ua, way(tags(k %in% c("building"))))
# bg_ids <- find_down(ua, way(bg_ids))
# # buildings 
# bg     <- subset(ua, ids = bg_ids)
# build  <- as_sp(bg, "polygons")

library(osmdata)

#getbb("Ventotene") %>%
q <-   opq(bbox = 'Ventotene')  %>%
  add_osm_feature(key = 'building') %>%
  osmdata_sp()
q
sp::plot(q$osm_polygons)
build <- q$osm_polygons
proj4string(build)

save(build, file = "Building.RData")

# change reference system
utm.build <- spTransform(build,proj4string(vent))

plot(vent)
plot(utm.build,add=T)

# defi ne a radius around each trap
radius <- 250
area   <- radius*radius*pi
btrap  <- gBuffer(utm.trap,width=radius,byid=T)

plot(vent)
plot(utm.trap,add=T)
plot(btrap,add=T)


# transform spatialLine to spatialPolygon

library(sf)
vent2     <- st_as_sf(vent) 
vent.pol  <- st_polygonize(vent2)
ventotene <- as(vent.pol, "Spatial") # If you want sp
class(ventotene)

ventotene2 <- spTransform(ventotene,crsutm)#gUnaryUnion(ventotene)
proj4string(ventotene2)
proj4string(btrap)


# cut off area in the sea and compute % of buildings in the area

# example
btrapnosea <- gIntersection(ventotene2,btrap[1],byid=TRUE)
build.buff <- gIntersection(utm.build,btrapnosea,byid=TRUE)
plot(ventotene2);plot(btrapnosea,add=T);plot(build.buff,add=T)

db$build <- NA
for(i in 1:nrow(db)){
  btrapnosea  <- gIntersection(ventotene2,btrap[i],byid=TRUE)
  area        <- gArea(btrapnosea)
  build.buff  <- gIntersection(utm.build,btrapnosea,byid=TRUE)
  buildings   <- gArea(build.buff)
  db$build[i] <- 100*buildings/area
}

db %>%group_by(IDsite,build)%>%count()%>%arrange(build)


btrapnosea  <- gIntersection(ventotene2,btrap,byid=TRUE)

par(mar=c(0,0,0,0))
plot(vent,lwd=2)
plot(utm.build,add=T)
plot(utm.trap,add=T,col="red",pch=19,cex=1.5)
plot(btrapnosea,add=T,lwd=2,fill="red")

isl <- fortify(ventotene2)
bld <- fortify(utm.build,by=ID)
nos <- fortify(btrapnosea)
trp <- as.data.frame(utm.trap@coords)

f1a <-ggplot()+geom_path(data= isl,aes(long,lat))+
  geom_polygon(data= nos,aes(long,lat,group=group),col="black",fill="grey90",alpha=0.05)+
  geom_path(data= bld,aes(long,lat,group=group))+
  geom_point(data = trp, aes(X1,X2),col="red",size=2)+
  theme_bw() + theme(legend.position = "top",
                     panel.grid.major.x = element_blank(),
                     panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
                     plot.title = element_text(size = rel(1.5), 
                                               face = "bold", vjust = 1.5),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     strip.background =  element_rect(fill="white"),
                     strip.text.x = element_text(size=15,face="italic"))+
  coord_quickmap()
f1a
ggsave(f1a, filename = "fig1a.pdf",width = 5, height = 4,dpi=600)

# exploratory analysis ###########################

summary(db)


# descriptive statistics #########################

# samplings when
range(unique(db$week))
length(unique(db$week))

# samplings where
length(unique(db$SiteID))

# trap position


trap.pos  <-db %>% group_by(Long,Lat) %>% count() 

trap.dist  <- spDists(as.matrix(trap.pos[,1:2]),longlat = T)*1000 #km to m
trap.dist2 <- trap.dist[lower.tri(trap.dist)]
hist(trap.dist2,breaks = 5)
plot(trap.dist2)
median(trap.dist2);mean(trap.dist2);max(trap.dist2);min(trap.dist2)

which.min.pos <- function(x){ which.min(x[x>0])}
apply(trap.dist,2,which.min.pos)
mean(c(208.7258,291.9419,378.4361,317.5654,161.0367,236.2301))
median(c(208.7258,291.9419,378.4361,317.5654,161.0367,236.2301))



db %>% group_by(Species) %>% 
  summarise(Tot    = sum(value),
            Mean   = mean(value),
            std    = sd(value),
            Median = median(value),
            quantL = quantile(value,0.25),
            quantH = quantile(value,0.75))



library(MASS)
db %>% group_by(Species,IDsite,build) %>% 
  summarise(Tot    = sum(value)) %>%
  split(.$Species) %>% map(.,~glm.nb(Tot~build, data = .)) %>%
  map(summary)

# no missing values
table(db$week,db$SiteID)


# mosquito capture distribution

db %>% ggplot(.,aes(value))+
  facet_wrap(~Species)+geom_histogram()+  theme_bw()

db %>% ggplot(.,aes(value))+
  facet_wrap(Species~IDsite)+geom_histogram()+  theme_bw()

# % of zeroes

db %>% mutate(zero = case_when(value == 0~1,
                               TRUE~0)) %>%
  summarise(p_zero = mean(zero))


db %>% mutate(zero = case_when(value == 0~1,
                               TRUE~0)) %>%
  group_by(Species) %>%
  summarise(p_zero = mean(zero))


db %>% mutate(zero = case_when(value == 0~1,
                               TRUE~0)) %>%
  group_by(SiteID) %>%
  summarise(p_zero = mean(zero))


db %>% mutate(zero = case_when(value == 0~1,
                               TRUE~0)) %>%
  group_by(Species,SiteID) %>%
  summarise(p_zero = mean(zero))%>%
  ggplot(.,aes(x=SiteID,y=p_zero,fill=Species))+
  geom_bar(stat="identity",position = position_dodge())+theme_bw()+
  scale_fill_manual(values = diverge_hcl(2))

# overall mosquito dynamics

f1c <- db %>% group_by(Species,week) %>% 
  summarise(Tot    = sum(value)) %>%
  ggplot(., aes(x=week,y=Tot,col=Species))+
  geom_line()+ geom_point()+
  geom_vline(xintercept = 32,linetype = "dashed")+
  scale_color_manual(values = diverging_hcl(2)) +
  theme_bw()  +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5),
        strip.background =  element_rect(fill="white"),
        strip.text.x = element_text(size=15,face="italic")
  )
f1c 

ggsave(f1c, filename = "fig1c.pdf",width = 5, height = 4,dpi=600)

# a disinfestation was reported on week 32, which is consistent with the observed pattern



# mosquito dynamics (spatio/temporal)

db %>% ggplot(., aes(x=week,y=value,col=Species))+
  facet_wrap(~IDsite)+
  geom_point()+  theme_bw()


db %>% group_by(IDsite,Long,Lat)%>%
  summarise(Tot = sum(value))%>%
  ggplot(.)+
    geom_point(aes(x=Long,y=Lat,size=Tot))+  theme_bw()


db %>% filter(Species == "Aedes") %>% group_by(IDsite,Long,Lat)%>%
  summarise(Tot = sum(value))%>%
  ggplot(.)+
  geom_point(aes(x=Long,y=Lat,size=Tot))+  theme_bw()


db %>% filter(Species == "Culex") %>% group_by(IDsite,Long,Lat)%>%
  summarise(Tot = sum(value))%>%
  ggplot(.)+
  geom_point(aes(x=Long,y=Lat,size=Tot))+  theme_bw()

# first detection

db %>% group_by(Species,week) %>% summarise(det = sum(value>0)) %>%
  filter(week <25)

db %>% group_by(Species,week,IDsite) %>% summarise(det = sum(value>0)) %>%
  filter(week ==21)

db %>% group_by(Species,week,IDsite) %>% summarise(det = sum(value>0)) %>%
  filter(week ==18)
# 
# f1b <- db %>% group_by(Species,IDsite) %>% summarise(Tot  = sum(value),
#                                                      Mean = mean(value))%>%
# ggplot(., aes(x=IDsite,y=Mean))+ theme_bw()+
#   geom_bar(stat="identity",col="black",fill="grey60")+
#   facet_grid( ~Species)+
#   labs(x = "ID Site", 
#        y = "N° of trapped mosquito")  +
#   theme_bw()  +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
#         plot.title = element_text(size = rel(1.5), 
#                                   face = "bold", vjust = 1.5),
#         axis.title = element_text(face = "bold"),
#         axis.title.y = element_text(vjust= 1.8),
#         axis.title.x = element_text(vjust= -0.5),
#         strip.background =  element_rect(fill="white"),
#         strip.text.x = element_text(size=15,face="italic")
#   )
# f1b


f1b <-   ggplot(db, aes(IDsite))+  # facet_grid( ~Species)+
  stat_summary(aes(y = value), fun = "mean", geom = "bar",fill = "grey70",col="black")+
  stat_summary(aes(y = value),geom="errorbar",fun.data = "mean_cl_boot", 
               colour = "black",width=0.1)+
  facet_grid( ~Species)+
  labs(x = "ID Site", 
       y = "N° of trapped mosquito")  +
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5),
        strip.background =  element_rect(fill="white"),
        strip.text.x = element_text(size=15,face="italic")
  )
f1b

ggsave(f1b, filename = "fig1b.pdf",width = 10, height = 4.5,dpi=600)


library(ggridges)
db %>% group_by(Species,IDsite,week) %>% summarise(Tot = sum(value))%>%
  ggplot(., aes(x=IDsite,y=week,height= Tot,group =week))+ theme_bw()+
  geom_ridgeline(alpha=0.75)+
  facet_grid( ~Species)+
  labs(x = "ID Trap", 
       y = "N° of trapped mosquito")  +
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5),
        strip.background =  element_rect(fill="white"),
        strip.text.x = element_text(size=15,face="italic")
  )


db$Species_site <- factor(paste(db$Species,db$IDsite,sep="_"))

db %>% ggplot(.,aes(x=IDsite,y=value,fill=Species))+geom_boxplot()+
  labs(x = "ID Site", 
       y = "N° of trapped mosquito")  +
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5),
        strip.background =  element_rect(fill="white"),
        strip.text.x = element_text(size=15,face="italic")
  )+
  scale_fill_manual(values = diverge_hcl(2))+
  scale_color_manual(values = diverge_hcl(2))







#### model #######

library(mgcv)
# citation("mgcv")
# Wood, S.N. (2017) 
# Generalized Additive Models: An Introduction with R (2nd edition). 
# Chapman and Hall/CRC.

# poisson distribution

m1 <- gam(value ~ Species+build+Species:build+
            s(week,by=Species) , 
          family = poisson, data = db)

# overdispersion
sum(resid(m1,type="pearson")^2)/(m1$df.residual)

# negative binomial distribution
m2 <- gam(value ~ Species+build+Species:build+
            s(week,by=Species) , 
          family = nb(), data = db)


# overdispersion
sum(resid(m2,type="pearson")^2)/(m2$df.residual) #ok

# model validation / check assumptions
gam.check(m2)
qq.gam(m2)

plot(resid(m2,type="pearson"))
plot(fitted(m2),resid(m2,type="pearson"))
plot(db$IDsite, resid(m2,type="pearson"))
plot(db$Species, resid(m2,type="pearson"))
plot(db$Lat, resid(m2,type="pearson"))
plot(db$Long,resid(m2,type="pearson"))
plot(factor(db$week),resid(m2,type="pearson"))
plot(factor(db$week),db$value)
# quite ok

db$E <- resid(m2,type="pearson")
ggplot(db, aes(week,E,col=Species))+
  facet_wrap(~IDsite)+geom_point()+theme_bw()
# temporal autocorrelation

for(i in levels(db$Species_site)){
  acf(db$E[db$Species_site==i] )}
# quite ok

# spatial autocorrelation
library(gstat)
mydb <- data.frame(e = db$E, lt = db$Lat,ln =db$Long)
coordinates(mydb)<- c("ln","lt")
vario <- variogram(object = e~1,
                   data = mydb,
                   cressie = TRUE)
vario
plot(vario,pch=19)
# ok

# not perfect but acceptable

# model results
summary(m2)

# increase from natural to urban gradient (albopcitus)
nat <- exp(coef(m2)[1] + coef(m2)[3]*0)
urb <- exp(coef(m2)[1] + coef(m2)[3]*20)

urb/nat



# plotting results
par(mfrow=c(1,2))
plot(m2)


newdb <- expand.grid(week = seq(min(db$week),max(db$week),by=1),
                     Species = levels(db$Species),
                     build = seq(0,20,by=1))

dim(newdb) 

newdb$fit <- predict(m2,newdata = newdb,type="response")
newdb$hi  <- exp(predict(m2,newdata = newdb) + 1.96* predict(m2,newdata = newdb,se.fit = TRUE)$se.fit)
newdb$lo  <- exp(predict(m2,newdata = newdb) - 1.96* predict(m2,newdata = newdb,se.fit = TRUE)$se.fit)

range(db$week)

f2a <- newdb %>% filter(week%in% c(35))%>%
  ggplot(.,aes(build,fit,col=Species))+
  geom_ribbon(aes(x=build,ymin= lo,ymax=hi,fill=Species),col="white",alpha=0.2)+
  geom_line()+
  xlab("% Area covered by buildings")+
  ylab("Mean mosquito abundance")+theme(legend.position = "top")+
  scale_fill_manual(values = diverge_hcl(2))+
  scale_color_manual(values = diverge_hcl(2))+
  ylim(c(0,23))+ 
  theme_bw()  +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5),
        strip.background =  element_rect(fill="white"),
        strip.text.x = element_text(size=15,face="italic")
  )

f2a

f2b <- newdb %>% filter(build == 10)%>%
  ggplot(.,aes(week,fit,col=Species))+
  geom_ribbon(aes(x=week,ymin= lo,ymax=hi,fill=Species),col="white",alpha=0.2)+
  geom_line()+
  xlab("Week")+
  ylab("Mean mosquito abundance")+
  scale_fill_manual(values = diverge_hcl(2))+
  scale_color_manual(values = diverge_hcl(2))+
  ylim(c(0,23)) +
  theme_bw()  +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5),
        strip.background =  element_rect(fill="white"),
        strip.text.x = element_text(size=15,face="italic")
  )
f2b



library(gridExtra)
f2 <- grid.arrange(f2a,f2b,ncol=2)
ggsave(f2, filename = "fig2.pdf",width = 8, height = 4,dpi=600)



# difference between temporal smoothers
library(mgcViz)
set.seed(20200702)
m3 <- getViz(m2)
plotDiff(s1 = sm(m3, 1), s2 = sm(m3, 2)) + l_ciPoly() + 
  l_fitLine() + geom_hline(yintercept = 0, linetype = 2)




# probability of zero presence

# check how size is computed in gam
# x <- rnbinom(1000, mu = 5,size =0.5)
# y <- rnorm(1000, mean =0 ,sd =3)
# mm <- gam(x~s(y),family = nb())
# summary(mm)

newdb$zerop <- dnbinom(0,mu = newdb$fit,size = 0.435)


dbfig3a <- newdb %>% filter(Species == "Aedes",build %in% c(0,5,10,15,20))
dbfig3b <- newdb %>% filter(Species == "Culex",build %in% c(10))

f3 <-  ggplot()+
    geom_line(data= dbfig3a,aes(week,1-zerop,col=factor(build)),size=1.1)+
    geom_line(data= dbfig3b,aes(week,1-zerop),col="black",linetype="solid",size=1.1)+
    xlab("Week")+
  ylab("Probability of observing a mosquito")+
  scale_fill_manual(values = diverge_hcl(5))+
  scale_color_manual(values = diverge_hcl(5))+
  ylim(c(0,1))+
  theme_bw()  +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(vjust= 1.8),
        axis.title.x = element_text(vjust= -0.5),
        strip.background =  element_rect(fill="white"),
        strip.text.x = element_text(size=15,face="italic")
  )+
  guides(col =guide_legend(title="% Area covered by buildings"))
f3

ggsave(f3, filename = "fig3.pdf",width = 6, height = 5,dpi=600)

# 
# f2d <- newdb %>% filter(Species == "Culex",build %in% c(0,0.05,0.1,0.15,0.2))%>%
#   ggplot(.,aes(week,zerop,col=factor(100*build)))+
#   geom_line(size=1.1)+
#   xlab("Week")+
#   ylab(expression(paste("Probability of observing ", italic(" Cx. pipiens"))))+
#   scale_fill_manual(values = diverge_hcl(5))+
#   scale_color_manual(values = diverge_hcl(5))+
#   ylim(c(0,1))+
#   theme_bw()  +
#   theme(legend.position = "top",
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey60",linetype = "dashed"),
#         plot.title = element_text(size = rel(1.5), 
#                                   face = "bold", vjust = 1.5),
#         axis.title = element_text(face = "bold"),
#         axis.title.y = element_text(vjust= 1.8),
#         axis.title.x = element_text(vjust= -0.5),
#         strip.background =  element_rect(fill="white"),
#         strip.text.x = element_text(size=15,face="italic")
#   )+
#   guides(col =guide_legend(title="% Area covered by buildings"))
# 
# f2d



# radius sensitivity

# defi ne a radius around each trap
radius <- seq(50,500,by=50)
res.rad <- array(NA,dim = c(4,2,length(radius)))
for(j in 1:length(radius)){
area   <- radius*radius[j]*pi
btrap  <- gBuffer(utm.trap,width=radius[j],byid=T)

db$build <- NA
for(i in 1:nrow(db)){
  btrapnosea  <- gIntersection(ventotene2,btrap[i],byid=TRUE)
  area        <- gArea(btrapnosea)
  build.buff  <- gIntersection(utm.build,btrapnosea,byid=TRUE)
  if(length(build.buff) >0 ){
  buildings   <- gArea(build.buff)
  db$build[i] <- buildings/area}else{db$build[i] = 0}
}

mod <- gam(value ~ Species+build+Species:build+
            s(week,by=Species) , 
            family = nb(), data = db)

res.rad[,1,j] <-  summary(mod)$p.coef
res.rad[,2,j] <-  summary(mod)$p.pv
}

# parameters are significant?
apply(res.rad[,2,]<0.05,1,sum)
# intercept / Species / buildings / interaction
res.rad[,2,]<0.05
# buildings not significant when radius <= 150
# interaction not signifant whan radius <= 200


# parameters sign similar?
# intercept / Species / buildings / interaction
res.rad[,1,]
# buildings consistently positive when radius > 50
# interaction consistently negative when radius > 100






