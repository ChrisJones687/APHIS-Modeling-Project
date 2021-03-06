pest_vars <<- list(host1_rast = NULL,host1_score = NULL, host2_rast=NULL,host2_score=NULL,host3_rast=NULL,host3_score=NULL, host4_rast=NULL,host4_score=NULL,host5_rast=NULL,host5_score=NULL,
                   host6_rast=NULL,host6_score=NULL,host7_rast=NULL,host7_score=NULL,host8_rast=NULL,host8_score=NULL,host9_rast=NULL,host9_score=NULL,host10_rast=NULL,host10_score=NULL,
                   allTrees=NULL,initialPopulation=NULL, start=2000, end=2010, seasonality = 'NO', s1 = 1 , s2 = 12, sporeRate = 4.4, windQ =NULL, windDir=NULL, tempQ="NO", tempData=NULL, 
                   precipQ="NO", precipData=NULL, kernelType ='Cauchy', kappa = 2, number_of_hosts = 1, scale1 = 20.57, seed_n = 42, scale2 = NULL, gamma = 1,  time_step = "weeks")


pest_vars$host1_rast = raster("./layers/UMCA_den_100m.img")
pest_vars$host2_rast = raster("./layers/OAKS_den_100m.img")
pest_vars$host3_rast = raster("./layers/LIDE3_den_100m.img")
pest_vars$allTrees = raster("./layers/TPH_den_100m.img")
pest_vars$initialPopulation = raster ("./layers/init_2000_cnt.img")
pest_vars$start = 2000
pest_vars$end = 2010
pest_vars$seasonality = 'YES'
pest_vars$s1 = 1
pest_vars$s2 = 9
pest_vars$sporeRate =4.4 ## spore rate default of original attempt
pest_vars$windQ ='YES'
pest_vars$windDir = 'NE'
pest_vars$tempQ = "YES"
pest_vars$precipQ = "YES"
pest_vars$tempData = './layers/weather/weatherCoeff_2000_2014.nc'
pest_vars$precipData = pest_vars$tempData
pest_vars$scale2 =NULL
pest_vars$number_of_hosts = 2
pest_vars$host1_score = 10
pest_vars$host2_score = 0
pest_vars$host3_score = 5
pest_vars$host4_score = NULL
pest_vars$host5_score = NULL
pest_vars$host6_score = NULL
pest_vars$host7_score = NULL
pest_vars$host8_score = NULL
pest_vars$host9_score = NULL
pest_vars$host10_score = NULL

dataList <- do.call(pest, pest_vars)

#setwd("C:\\Users\\chris\\Dropbox\\Projects\\Code\\Aphis Modeling Project")
I_oaks_rast2 <- pest(host1,host2,allTrees,initialPopulation, start, end, SS, s1, s2, sporeRate, windQ, windDir, tempData)
dataList <- pest(host1_rast =host1_rast,host2_rast=host2_rast,host3_rast = host3_rast, allTrees= allTrees,initialPopulation= initialPopulation, start =start, end=end, seasonality=SS, s1=s1, s2=s2, sporeRate=sporeRate, windQ=windQ, windDir=windDir, tempData=tempData, precipData=precipData, number_of_hosts = number_of_hosts, host1_score = host1_score, host2_score = host2_score, host3_score = host3_score)


dataReturn <- dataList[[1]]
I_oaks_rast2 <- dataList[[2]]

make1 <- dataReturn[,1:3]
make2 <- dataReturn[,c(1,4:5)]
names(make1) <- c('Year','Area','Count')
names(make2) <- c('Year','Area','Count')
make1$Host <- 'Tanoak'
make2$Host <- 'Oaks'
make <- rbind(make1,make2)
dataForPlot <- make

leaflet("stateData", height = "600px") %>% addMiniMap(position = "topright")

title = "Model Output"
theme = theme_set(theme_classic())
theme = theme_update(text = element_text(family = "Arial"), legend.position="top", legend.title=element_blank(),legend.spacing=unit(-0.5,"lines"), plot.background = element_rect(fill = "#3F3E3E", colour = "#3F3E3E"), panel.background = element_rect(fill = "#3F3E3E", colour = "#3F3E3E"), legend.background = element_rect(fill = "#3F3E3E"))
theme = theme_update(axis.text = element_text(size = 12, color="white"), axis.ticks=element_blank(), plot.title = element_text(hjust = 0.5,colour="white", size =18), axis.line = element_line(colour="white"),axis.title=element_text(size=16, vjust=0,35,color="white"),legend.text=element_text(size=12,color="white"))
ggplot(dataForPlot, aes(x=Year, y=Area, color=factor(Host)))+geom_line(aes(Year, Area), size = 1.5)+scale_color_manual(values=c("#54ACC1", "#ADBD60"))+scale_fill_manual(values=c("#54ACC1", "#ADBD60"))+
  ggtitle(title)+theme(text = element_text(family = "sans"))+
  scale_x_continuous(name="Year", breaks=seq(start, end, 2))+
  scale_y_continuous(name="yName")+guides(col=guide_legend(ncol=3),shape=guide_legend(ncol = 1))


## test sample to make sure it works as intended
prop_S_host1 =.1
prop_S_host2 =.1
prop_S_host3 =.3
prop_S_host4 =.05
prop_S_host5 =.05
prop_S_host6 =.05
prop_S_host7 =.05
prop_S_host8 =.2
prop_S_host9 =0
prop_S_host10 =.1
sample(seq_len(10), 1, replace = FALSE, prob = c(prop_S_host1, prop_S_host2,prop_S_host3, prop_S_host4, prop_S_host5,prop_S_host6,prop_S_host7,prop_S_host8,prop_S_host9,prop_S_host10))


## test implementation of zip file
write.csv(dataReturn, file = "output/data.csv")
writeRaster(dataList[[2]], file = "output/host1.tif")
files2zip <- dir("output", full.names = TRUE)
zip(zipfile = "C:/Users/Chris/Desktop/testZip2.zip", files = files2zip)
do.call(file.remove, list(list.files("output", full.names = TRUE)))

## move from .nc to .tif
file.nc <- "./layers/weather/weatherCoeff_2000_2014.nc"
temp <- raster(file.nc, varname = "Ccoef")
file.tif <- 'layers/weatherCoeff_2000_2014.tif'
precip <- raster(file.nc, varname = "Mcoef")
ccf.array2 <- as.array(precip)



writeRaster(d, "C:/Users/Chris/Desktop/dodtest.tif", overwrite = TRUE, format = 'GTiff')

for (tt in tstep){
  cnt <- cnt + 1
  if (cnt %in% yearlyoutputlist){
    print(cnt)
  }
}
class(cnt)
