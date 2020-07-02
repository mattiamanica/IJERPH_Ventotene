


# load packages
library(here)
library(tidyverse)

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
library(osmar)
src <- osmsource_api(url = "https://api.openstreetmap.org/api/0.6/")

bb <- corner_bbox(13.3949,40.7828,13.4485,40.8060)
ua <- get_osm(bb, source = src)


ll_poly <- as_sp(ua, "polygons")
spplot(ll_poly, c("version"))

# identify buildings
bg_ids <- find(ua, way(tags(k %in% c("building"))))
bg_ids <- find_down(ua, way(bg_ids))
# buildings 
bg     <- subset(ua, ids = bg_ids)
build  <- as_sp(bg, "polygons")

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
  db$build[i] <- buildings/area
}




# exploratory analysis ###########################

summary(db)

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

db %>% group_by(Species,week) %>% summarise(Tot = sum(value)) %>%
  ggplot(., aes(x=week,y=Tot,col=Species))+
  geom_line()+  theme_bw()

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



db %>% group_by(Species,IDsite) %>% summarise(Tot = sum(value))%>%
ggplot(., aes(x=IDsite,y=Tot))+ theme_bw()+
  geom_bar(stat="identity",col="black",fill="grey60")+
  facet_grid( ~Species)+
  labs(x = "ID Trap", 
       y = "N? of trapped mosquito")  +
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


library(ggridges)
db %>% group_by(Species,IDsite,week) %>% summarise(Tot = sum(value))%>%
  ggplot(., aes(x=IDsite,y=week,height= Tot,group =week))+ theme_bw()+
  geom_ridgeline(alpha=0.75)+
  facet_grid( ~Species)+
  labs(x = "ID Trap", 
       y = "N? of trapped mosquito")  +
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

















