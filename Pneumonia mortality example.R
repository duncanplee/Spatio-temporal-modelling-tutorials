############################
#### R code for the tutorial
############################

set.seed(1)
##############################################
#### Section 3 - Data and exploratory analysis
##############################################
#### Read in the data
dat <- read.csv(file="EnglandLUAdata.csv")
head(dat, n=3)


#### Summarise the variables
summary(dat[ ,4:7])


#### Add the SMR
library(dplyr)
dat <- dat %>% mutate(SMR=dat$Y/dat$E)


#### Scatterplot
library(GGally)
ggpairs(dat, columns=6:8)


#### Boxplots for temporal trend
library(ggplot2)
ggplot(dat, aes(x = factor(Year), y = SMR)) +
    geom_boxplot(fill="red", alpha=0.7) + 
    scale_x_discrete(name = "Year", breaks=c(2002, 2005, 2008, 2011, 2014, 2017), labels=c("2002", "2005", "2008", "2011", "2014", "2017")) +
    scale_y_continuous(name = "SMR") + 
    theme(text=element_text(size=16), plot.title=element_text(size=18, face="bold")) 



#### Read in the spatial object
library(rgdal)
LA <- readOGR(dsn = "LocalAuthorities.shp")


#### Compute the average SMR for each year
by_LA <- group_by(dat, Code)
averageSMR <- summarize(by_LA, SMR = mean(SMR, na.rm=T))


#### Merge the data and shapefile
library(leaflet)
averageSMR.LA <- merge(x=LA, y=averageSMR, by.x="lad09cd", by.y="Code", all.x=FALSE)
averageSMR.LA.ll <- spTransform(averageSMR.LA, CRS("+proj=longlat +datum=WGS84 +no_defs"))
variable <- averageSMR.LA.ll@data$SMR
colours <- colorNumeric(palette = "YlOrBr", domain = variable, reverse=FALSE)
leaflet(data=averageSMR.LA.ll) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(variable), 
                color="",
                fillOpacity = 0.7, weight = 1, smoothFactor = 0.5,
                opacity = 1.0) %>%
    addLegend(pal = colours, values = variable, 
              opacity = 1, title="SMR") %>%
    addScaleBar(position="bottomleft")



#### Fit a simple Poisson log-linear model
model1 <- glm(formula=Y~offset(log(E)) + IMD + PM25, family="poisson", data=dat)
round(cbind(model1$coefficients, confint(model1)),4)


#### compute the residuals from this model
dat$residuals <- residuals(model1)
residuals2010 <- filter(dat, Year==2010)
residuals2010.LA <- merge(x=LA, y=residuals2010, by.x="lad09cd", by.y="Code", all.x=FALSE)


#### Construct the spatial objects
library(spdep)
W.nb <- poly2nb(residuals2010.LA, row.names = residuals2010.LA@data$lad09cd)
W <- nb2mat(W.nb, style = "B")
W.list <- nb2listw(W.nb, style = "B")


#### Conduct Moran's I test
moran.mc(x = residuals2010.LA$residuals, listw = W.list, nsim = 10000)



#####################################################################
#### Section 4 - Spatio-temporal modelling and convergence assessment
#####################################################################
#### Order the data according to the neighbourhood matrix and year
lookup <- data.frame(Code=residuals2010.LA@data$lad09cd, spatialorder=1:nrow(residuals2010.LA@data))
dat.temp <- merge(x=dat, y=lookup, by="Code")
dat.ordered <- arrange(dat.temp, Year, spatialorder)


#### Fit the model
library(CARBayesST)
chain1 <- ST.CARar(formula=Y~offset(log(E)) + PM25 + IMD, family="poisson", data=dat.ordered, W=W, burnin=200000, n.sample=2200000, thin=1000, verbose=FALSE)
chain2 <- ST.CARar(formula=Y~offset(log(E)) + PM25 + IMD, family="poisson", data=dat.ordered, W=W, burnin=200000, n.sample=2200000, thin=1000, verbose=FALSE)
chain3 <- ST.CARar(formula=Y~offset(log(E)) + PM25 + IMD, family="poisson", data=dat.ordered, W=W, burnin=200000, n.sample=2200000, thin=1000, verbose=FALSE)

#### Check convergence - traceplot
library(coda)
beta.samples <- mcmc.list(chain1$samples$beta, chain2$samples$beta, chain3$samples$beta)
plot(beta.samples)


#### Check convergence - Gelman-Rubin plot
gelman.diag(beta.samples)


#### Model summary
print(chain1)



##########################
#### Section 5 - Inference
##########################
#### Effects of covariates on disease risk
sd(dat.ordered$PM25)
sd(dat.ordered$IMD)
beta.samples.combined <- rbind(chain1$samples$beta, chain2$samples$beta, chain3$samples$beta)
round(quantile(exp(sd(dat.ordered$PM25) * beta.samples.combined[ ,2]), c(0.5, 0.025, 0.975)),3)
round(quantile(exp(sd(dat.ordered$IMD) * beta.samples.combined[ ,3]), c(0.5, 0.025, 0.975)),3)


#### Compute the risk distributions
fitted.samples.combined <- rbind(chain1$samples$fitted, chain2$samples$fitted, chain3$samples$fitted)
n.samples <- nrow(fitted.samples.combined)
n.all <- ncol(fitted.samples.combined)
risk.samples.combined <- fitted.samples.combined / matrix(rep(dat.ordered$E, n.samples), nrow=n.samples, ncol=n.all, byrow=TRUE) 


#### Compute the areal unit average risk for each year
N <- length(table(dat.ordered$Year))
risk.trends <- array(NA, c(n.samples, N))
for(i in 1:n.samples)
{
    risk.trends[i, ] <- tapply(risk.samples.combined[i, ], dat.ordered$Year, mean)
}


#### Plot the average risk trends
time.trends <- as.data.frame(t(apply(risk.trends, 2, quantile, c(0.5, 0.025, 0.975))))
time.trends <- time.trends %>% mutate(Year=names(table(dat.ordered$Year)))
colnames(time.trends)[1:3] <- c("Median","LCI", "UCI")

ggplot(time.trends, aes(x = factor(Year), y = Median, group=1)) +
    geom_line(col="red") + 
    geom_line(aes(x=factor(Year), y=LCI), col="red", lty=2) +
    geom_line(aes(x=factor(Year), y=UCI), col="red", lty=2) + 
    scale_x_discrete(name = "Year", breaks=c(2002, 2005, 2008, 2011, 2014, 2017), labels=c("2002", "2005", "2008", "2011", "2014", "2017")) +
    scale_y_continuous(name = "Risk") + 
    theme(text=element_text(size=16), plot.title=element_text(size=18, face="bold")) 


#### Spatial pattern in disease risk in the last year - mean and PEP
risk.samples.2010 <- risk.samples.combined[ ,dat.ordered$Year==2010]
risk.2010 <- apply(risk.samples.2010, 2, median)
pep.2010 <- apply(risk.samples.2010 > 1, 2, mean)


#### Map the results
residuals2010.LA$risk.2010 <- risk.2010
residuals2010.LA$pep.2010 <- pep.2010
residuals2010.LA.ll <- spTransform(residuals2010.LA, CRS("+proj=longlat +datum=WGS84 +no_defs"))



colours <- colorNumeric(palette = "YlOrBr", domain = residuals2010.LA.ll@data$risk.2010, reverse=FALSE)
leaflet(data=residuals2010.LA.ll) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(risk.2010), 
                color="",
                fillOpacity = 0.7, weight = 1, smoothFactor = 0.5,
                opacity = 1.0) %>%
    addLegend(pal = colours, values = risk.2010, 
              opacity = 1, title="Risk") %>%
    addScaleBar(position="bottomleft")


colours <- colorNumeric(palette = "YlOrBr", domain = residuals2010.LA.ll@data$pep.2010, reverse=FALSE)
leaflet(data=residuals2010.LA.ll) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(pep.2010), 
                color="",
                fillOpacity = 0.7, weight = 1, smoothFactor = 0.5,
                opacity = 1.0) %>%
    addLegend(pal = colours, values = pep.2010, 
              opacity = 1, title="PEP") %>%
    addScaleBar(position="bottomleft")



#### Compute the median risk for each area
risk.median <- apply(risk.samples.combined, 2, median)
inequality <- tapply(risk.median, dat.ordered$Year, IQR)
ggplot(data.frame(inequality, Year=names(inequality)), aes(x = factor(Year), y = inequality, group=1)) +
    geom_line(col="red") + 
    scale_x_discrete(name = "Year", breaks=c(2002, 2005, 2008, 2011, 2014, 2017), labels=c("2002", "2005", "2008", "2011", "2014", "2017")) +
    scale_y_continuous(name = "Inequality") + 
    theme(text=element_text(size=16), plot.title=element_text(size=18, face="bold")) 

