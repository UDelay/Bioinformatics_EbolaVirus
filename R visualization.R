##### Welcome to R visualization ####
# Get data:

#install.packages("rlang")
library(gapminder)

# Charge libraries:
library(ggplot2)
library(gganimate)
#install.packages("gifski")
#install.packages("av")

# Make a ggplot, but add frame=year: one image per year
ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, color = continent)) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  # gganimate specific bits:
  labs(title = 'Year: {frame_time}', x = 'GDP per capita', y = 'life expectancy') +
  transition_time(year) +
  ease_aes('linear')

# Save at gif:
anim_save("271-ggplot2-animated-gif-chart-with-gganimate1.gif")



#### World map ####


library(leaflet)

m <- leaflet() %>% 
  addTiles() %>% 
  setView( lng = 2.34, lat = 48.85, zoom = 5 ) %>% 
  addProviderTiles("NASAGIBS.ViirsEarthAtNight2012")
m


### chord diagram ###


# Load package

#devtools::install_github("mattflor/chorddiag")
library(chorddiag)

# Create dummy data
m <- matrix(c(11975,  5871, 8916, 2868,
              1951, 10048, 2060, 6171,
              8010, 16145, 8090, 8045,
              1013,   990,  940, 6907),
            byrow = TRUE,
            nrow = 4, ncol = 4)

# A vector of 4 colors for 4 groups
haircolors <- c("black", "blonde", "brown", "red")
dimnames(m) <- list(have = haircolors,
                    prefer = haircolors)
groupColors <- c("#000000", "#FFDD89", "#957244", "#F26223")

# Build the chord diagram:
p <- chorddiag(m, groupColors = groupColors, groupnamePadding = 20)
p



#### Hexagonal plots on United States map ####


library(viridis)
library(tidyverse)
library(geojsonio)
library(RColorBrewer)
library(rgdal)
library(broom)



# Download the Hexagones boundaries at geojson format here: https://team.carto.com/u/andrew/tables/andrew.us_states_hexgrid/public/map.

# Load this file. (Note: I stored in a folder called DATA)
spdf <- geojson_read("DATA/us_states_hexgrid.geojson.json",  what = "sp")

# Bit of reformating
spdf@data = spdf@data %>%
  mutate(google_name = gsub(" \\(United States\\)", "", google_name))



spdf@data = spdf@data %>% mutate(google_name = gsub(" \\(United States\\)", "", google_name))
spdf_fortified <- tidy(spdf, region = "google_name")

# Calculate the centroid of each hexagon to add the label:
library(rgeos)
centers <- cbind.data.frame(data.frame(gCentroid(spdf, byid=TRUE), id=spdf@data$iso3166_2))


# Now I can plot this shape easily as described before:
ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group), fill="skyblue", color="white") +
  geom_text(data=centers, aes(x=x, y=y, label=id)) +
  theme_void() +
  coord_map()

# Distribution of the marriage rate?
data <- read.table("https://raw.githubusercontent.com/holtzy/R-graph-gallery/master/DATA/State_mariage_rate.csv", header=T, sep=",", na.strings="---")

# Distribution of the marriage rate?

spdf_fortified <- spdf_fortified %>%
  left_join(. , data, by=c("id"="state")) 

# Make a first chloropleth map
ggplot() +
  geom_polygon(data = spdf_fortified, aes(fill =  y_2015, x = long, y = lat, group = group)) +
  scale_fill_gradient(trans = "log") +
  theme_void() +
  coord_map()


# Prepare binning
spdf_fortified$bin <- cut( spdf_fortified$y_2015 , breaks=c(seq(5,10), Inf), labels=c("5-6", "6-7", "7-8", "8-9", "9-10", "10+" ), include.lowest = TRUE )

# Prepare a color scale coming from the viridis color palette
my_palette <- rev(magma(8))[c(-1,-8)]

# plot
ggplot() +
  geom_polygon(data = spdf_fortified, aes(fill = bin, x = long, y = lat, group = group) , size=0, alpha=0.9) +
  geom_text(data=centers, aes(x=x, y=y, label=id), color="white", size=3, alpha=0.6) +
  theme_void() +
  scale_fill_manual( 
    values=my_palette, 
    name="Wedding per 1000 people in 2015", 
    guide = guide_legend( keyheight = unit(3, units = "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) 
  ) +
  ggtitle( "A map of marriage rates, state by state" ) +
  theme(
    legend.position = c(0.5, 0.9),
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 22, hjust=0.5, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm"))
  )


#### let's start from basic R plots #####


### scatter plots ###
dataSet  = data.frame(aa = c(1,2,3,4),
                      bb = c(2,1,3,4),
                      cc = c(10,8,4,1))

#plot(dataSet$aa, dataSet$bb)
#plot(dataSet$aa, dataSet$bb, xlab="AA", ylab = "BB")

plot(dataSet$aa, dataSet$bb, xlab="AA", ylab = "BB", col = c("red", "red", "blue", "blue"), pch=16, cex=5)



#### box plots ####

dataSet  = data.frame(aa = c(1,2,3,4),
                      bb = c(2,1,3,4),
                      cc = c(10,8,4,1))

boxplot(dataSet)




#### simple statistical tests ####

t.test(dataSet$aa, dataSet$cc)


fit = lm(bb ~ aa, dataSet)
summary(fit)

#plot(fit) ### diagnostic plots

plot(dataSet$bb, predict(fit), pch=16, cex=3)
abline(a=0,b=1, col="red", lty="dashed")



### heatmap.2 ###

require(gplots)


heatmapExampleData = "heatmapExample.RData"
load(heatmapExampleData) 

blueCol = "#1D2088"
redCol = "#E60012"

heatmap.2(heatmapExample, 
          trace="none",
          key=T,
          sepcolor = "black", 
          colsep = c(0,dim(heatmapExample)[2]),
          rowsep = c(0,dim(heatmapExample)[1]),
          cexRow = 1.2, cexCol = 1.2,
          col = colorRampPalette(c(blueCol, "white", "white", redCol))(16))





#### ggplots ####
require(ggplot2)
require(ggsci)
data(iris)
head(iris)

### scatter plots ###
ggplot(iris) + geom_point(aes(x=Sepal.Length, y=Sepal.Width, colour=Species, size=3)) + scale_colour_aaas()

### scatter plots with white background ###
ggplot(iris) + geom_point(aes(x=Sepal.Length, y=Sepal.Width, colour=Species, size=3)) + scale_colour_aaas() + theme_bw()




### box plots ###
ggplot(iris) + geom_boxplot(aes(x=Species, y=Sepal.Width))

### box + jitterplots ###
ggplot(iris) + geom_boxplot(aes(x=Species, y=Sepal.Width)) + geom_jitter(aes(x=Species, y=Sepal.Width))

### box + jitterplots + group colors + white background ###
ggplot(iris) + geom_boxplot(aes(x=Species, y=Sepal.Width, fill=Species)) + geom_jitter(aes(x=Species, y=Sepal.Width)) + theme_bw()


### lineplots + group colors + white background ###
ggplot(iris) + geom_line(aes(x=Sepal.Length, y=Sepal.Width, colour=Species)) + theme_bw()



### scatter plots + fitted lines (linear regressions) + white background ###
ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width, colour=Species)) + geom_point() + geom_smooth(formula = y ~ x, method = "lm", aes(fill=Species)) + theme_bw()



