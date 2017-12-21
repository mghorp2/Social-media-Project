library(dplyr)
library(igraph)
library(randomForest)
library(caret)
library(cluster)

getwd()
# Save the data file to a location on your hard drive and specify the path here (Windows systems use forward slashes)
dir_path <-"D:/UIC/Fall 2017/Social Media Networking/Project"
setwd(dir_path)
# clear everything out of memory
rm(list=ls())  

# Load primary school data, contact data
infile_edges<-"Edges.csv"
infile_nodes<-"Nodes.csv"
infile_predictors="With Predictors.csv"

## Load package

edge_frame=read.csv(infile_edges, header = TRUE, sep = ",")
node_frame=read.csv(infile_nodes, header = TRUE, sep = ",")
predictor_var=read.csv(infile_predictors,header = TRUE, sep = ",")

g.class=graph.data.frame(edge_frame, directed = FALSE, vertices= node_frame)

# Edges
ecount(g.class)
## Vertices
vcount(g.class)
is.weighted(g.class)

V(g.class)$name
E(g.class)$weight
#V(g_primschool)[V(g_primschool)$classname=="1B"]

is.simple(g.class)
is.connected(g.class,mode="weak")
is.connected(g.class,mode="strong")

community.significance.test <- function(graph, vs, ...) {
  if (is.directed(graph)) stop("This method requires an undirected graph")
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  # Total degree among nodes in the vs list, minus the degree within the subgraph 
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}

stud.class <- get.vertex.attribute(g.class, "Identity")
stud.weight <- E(g.class)$weight
plot(g.class)

E(g.class)$weight 
g.class_simpl<-simplify(g.class)

### The above should default to the option below, to sum the existing edge weights ### when combining them
##g_acq_simpl<-simplify(g_acq,edge.attr.comb="sum" )

E(g.class_simpl)$weight 
# Will use the inverse of log weight for shortest path calculations
inv_weight<-1/log(E(g.class_simpl)$weight  + 1)
num_weight<-E(g.class_simpl)$weight 
length(inv_weight)
#get.adjacency(g.class_simpl)
diameter(g.class_simpl,directed = FALSE,weights = NA)
average.path.length(g.class_simpl,directed = FALSE)

B <-betweenness(g.class_simpl,v=V(g.class_simpl),directed = FALSE, nobigint= TRUE, normalized = TRUE, weights = inv_weight)
V(g.class_simpl)[order(B,decreasing = T)[1:4]]
edgebetweens<-edge.betweenness(g.class_simpl, e=E(g.class_simpl), directed = TRUE)
E(g.class_simpl)[order(B, decreasing=T)[1:4]]

C<-closeness(g.class_simpl,v=V(g.class_simpl), mode = "all",weights = inv_weight)
summary(C)
V(g.class_simpl)[order(C,decreasing = T)[1:4]]
hist(C, col="lightblue", xlab="Vertex Closeness", ylab="Frequency", main="Closeness Distribution")
V(g.class_simpl)[order(C,decreasing = T)[1:4]]

transitivity(g.class_simpl)
clusters(g.class_simpl)
reciprocity(g.class_simpl)



#making width of edges proportional to edge weight

E(g.class_simpl)$width <- 0.02*(num_weight)
plot(g.class_simpl)
# making the size of node proportional to its degree centrality
Degree_cen <- degree(g.class_simpl)
summary(Degree_cen)
print(Degree_cen)
degree_c<-data.frame(Degree_cen)
degree_c$nodes<-NA
degree_c$nodes<-rownames(degree_c)
degr <- degree_c %>% arrange(desc(Degree_cen))
d<-top_n(degr,3,degr$Degree_cen)

hist(Degree_cen, col="lightblue", xlab="Vertex Degree", ylab="Frequency", main="Degree Distribution")
hist(graph.strength(g.class_simpl), col="pink",xlab="Vertex Strength", ylab="Frequency", main="vertex strength distibution")

a.nn.deg <- graph.knn(g.class_simpl,V(g.class_simpl))$knn
plot(Degree_cen, a.nn.deg, log="xy", 
     col="goldenrod", xlab=c("Log Vertex Degree"),
     ylab=c("Log Average Neighbor Degree"), main = "Average neighbor degree versus vertex degree")


constraints_class <- round(constraint(g.class_simpl, nodes=V(g.class_simpl)), digits=4)
# Local clustering coefficients
clustering_class <- transitivity(g.class_simpl, type="local", vids=V(g.class_simpl)) 
clustering_class


node_frame<-data.frame(B, constraints_class, clustering_class,Degree_cen)
a_node<-aggregate(B ~ clustering_class, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Clustering", ylab="Average Betweenness of nodes")

par(mfrow=c(2, 1))
a_node<-aggregate(betweens_SAP ~ degree_sap, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Betweenness")
a_node<-aggregate(clustering_SAP ~ degree_sap, data=node_frame, mean)
plot(a_node, col="blue", log="xy", xlab="Degree", ylab="Average Clustering")

V(g.class_simpl)$size<- 4*sqrt(Degree_cen)
v.lable<-V(g.class_simpl)$size
V(g.class_simpl)$shape<-"circle"
V(g.class_simpl)$color<-4*sqrt(Degree_cen)
V(g.class_simpl)[37]$shape<-"square"
V(g.class_simpl)[37]$color<-"blue"
plot(g.class_simpl)


# Using kamada kawai layout
set.seed(33)
l1 <- layout.reingold.tilford(g.class_simpl,circular=T)
plot(g.class_simpl,layout=l1)
l2 <- layout.reingold.tilford(g.class_simpl)
plot(g.class_simpl,layout=l2)


table(sapply(cliques(g.class_simpl), length))

cliques(g.class_simpl)[sapply(cliques(g.class_simpl), length) == 4]



# community detection algorithms


community.significance.test <- function(graph, vs, ...) {
  if (is.directed(graph)) stop("This method requires an undirected graph")
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  # Total degree among nodes in the vs list, minus the degree within the subgraph 
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}

stud.class <- get.vertex.attribute(g.class_simpl, "classname")

class_fg <- fastgreedy.community(g.class_simpl,weights=inv_weight)
length(class_fg)
sizes(class_fg)
plot(class_fg,g.class_simpl)
c.m.f <- membership(class_fg)



class_wt <- walktrap.community(g.class_simpl,weights=inv_weight)
length(class_wt)
sizes(class_wt)
plot(class_wt,g.class_simpl)
c.m.w <- membership(class_wt)
c.m.w


class_eb <- edge.betweenness.community(g.class_simpl,weights=inv_weight)
length(class_eb)
sizes(class_eb)
plot(class_eb,g.class_simpl)
c.m.e <- membership(class_eb)
c.m.e
f.clusters<-data.frame(class_eb$membership)
f.names<-data.frame(class_eb$names)
edge_bet<-cbind(f.names,f.clusters)

combined_data<-merge(x = predictor_var, y = edge_bet, by.x = "Name",by.y="class_eb.names", all.x = TRUE)
head(combined_data,4)


# using Random forest for checking the variable importance

fit <- randomForest(x = combined_data[,2:5], y = as.factor(combined_data[,6]), ntree = 10000, proximity = TRUE, oob.prox = TRUE)
importance(fit)
varImpPlot(fit,main="variable importance")

# using Random forest for clustering

rf.fit <- randomForest(x = combined_data[,2:5], y = NULL, ntree = 10000, proximity = TRUE, oob.prox = TRUE)
hclust.rf <- hclust(as.dist(1-rf.fit$proximity), method = "ward.D2")
rf.cluster = cutree(hclust.rf, k=10)
combined_data$rf.clusters <- rf.cluster
head(combined_data,4)
write.csv(combined_data,file="clusters.csv")

par(mfrow=c(1,2))
plot(combined_data$class_eb.membership,combined_data$Name,col=combined_data$class_eb.membership,xlab = "Cluster Numbers",ylab= "Students",main="Clustering based on Girvan-Newman (Edge Betweeness) Algorithm")
plot(rf.cluster,combined_data$Name,col=rf.cluster,xlab = "Cluster Numbers",ylab= "Students",main="Clustering based on Random forest")

table(rf.cluster, combined_data$class_eb.membership)
