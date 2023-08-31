
#install.packages("ggraph")
require(ggraph)
require(igraph)

g=make_graph(~A-B-C-A, D-E-F-D, A-F, A-H, H-G-J-K-L-M, D-O, E-P, Q-R-S, S-Q) 
plot(g)

g #graph object 

### drawing network graph ###
ggraph(g) + geom_edge_link() + geom_node_point()

### vertex = node #node name check
V(g)$name

### edge 
E(g)

###adjacency matrix
as_adjacency_matrix(g) #matrix make 

### edges
as_edgelist(g)
as_data_frame(g)



### extracting centrality ###

sort(degree(g))
sort(betweenness(g)) #A>H>...
sort(closeness(g))
sort(eigen_centrality(g)$vector)
sort(page.rank(g)$vector)


de = degree(g)
be = betweenness(g)
cl = closeness(g)

### clustering coefficients ### (connectivity high)
g.cluster=transitivity(g, "global")
g.cluster


l.cluster=transitivity(g, "local")
names(l.cluster) = V(g)$name
l.cluster
#sort(l.cluster) > class time 

## cliques(most dense connected), component
largest_cliques(g)
components(g)



### module (community) detection ###

eb=edge.betweenness.community(g) #edge betweeness
eb
modularity(eb)
plot(eb, g)


com=cluster_louvain(g) #louvain algorithm
modularity(com)
plot(com, g)
#plot(com, g, vertex.label="")


wt=walktrap.community(g) #random-walk algorithm
wt
modularity(wt)
plot(wt, g)
plot(wt, g, vertex.label="")



### network visualization ###
library(RColorBrewer)
colors=brewer.pal(length(com),'Accent') #make a color palette
V(g)$color=colors[membership(com)] #assign each vertex a color based on the community assignment

set.seed(2)
plot(g, vertex.label="", edge.width=E(g)$weight*5)

plot(g, layout=layout_with_fr(g))




#size by degree
plot(g, vertex.label="", vertex.color="gold", edge.color="slateblue", 
     vertex.size=de*5, edge.width=E(g)$weight*5)

#size by betweenness
plot(g, vertex.label="", vertex.color="gold", edge.color="slateblue", 
     vertex.size=be*0.5, edge.width=E(g)$weight*5)

#size by closeness
plot(g, vertex.label="", vertex.color="gold", edge.color="slateblue", 
     vertex.size=cl*1000, edge.width=E(g)$weight*5)





### example! of air transport network
#install.packages("igraphdata")
library(igraphdata) 
data(USairports) 
USairports


#airport information
airports=read.csv('https://raw.githubusercontent.com/jpatokal/openflights/master/data/airports.dat', header=F)
head(airports)

V(USairports)$lat=airports[match(V(USairports)$name, airports[,5]), 7] 
V(USairports)$long=airports[match(V(USairports)$name, airports[,5]), 8]

usair=as.undirected(simplify(USairports))
usair=delete.vertices(usair, which(is.na(V(usair)$lat)==TRUE))

#remove nodes in the Eastern and Southern Hemispheres (US territories). This will make the plot easier to see.
usair=delete.vertices(usair, which(V(usair)$lat<0))
usair=delete.vertices(usair, which(V(usair)$long>0))

#keep only the largest connected component of the network ("giant component"). this also makes the network easier to see.
decomp=decompose.graph(usair)
usair=decomp[[1]]

set.seed(3)
l=layout_with_fr(usair)
par(mar=c(1,1,1,1))
plot(usair, layout=l, vertex.label="", vertex.size=3)

par(mar=c(5,4,4,1))
plot(degree(usair), betweenness(usair), pch=19) #pch=19 uses filled circles
sort(betweenness(usair),T)[1:5]
#ANC      SEA      DEN      MSP      ORD 
#64157.39 22259.14 21397.72 19975.46 16739.42 
#ANC = Alaska #ANC degree low but betweenness high
#SEA = Seattle
#DEN = Denver 

longlat=matrix(c(V(usair)$long, V(usair)$lat), ncol=2) #set up layout matrix 
par(mar=c(1,1,1,1))
plot(usair, layout=longlat, vertex.label="", vertex.size=3)

#require(devtools)
#install_github("dyerlab/popgraph")
#BiocManager::install("ggmap")
library(popgraph)
library(ggmap)

location = c( mean(V(usair)$long), mean(V(usair)$lat))
map=get_stamenmap(location, bbox=c(left=-179, bottom=15, right=-65, top=75), zoom=2, source="stamen", maptype="terrain-background", filetype="png")

p=ggmap(map)
p = p + geom_edgeset(aes(x=long, y=lat), usair, colour=gray(0.1, 0.3), size=1) #plot will also include the network edges, in white
p = p + geom_nodeset(aes(x=long, y=lat), usair, size=1, colour="tomato") #plot will also include the network nodes, color-coded by region
p #plot color
