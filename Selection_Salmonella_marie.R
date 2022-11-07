# SiSAGE strain selection


# 0. Libraries  
library(cluster) 
library(fpc)
library(ggplot2)
library(dplyr)
library(dendextend)
library(UpSetR)

data<-read.csv("sisage.csv",sep=";")
data$CONTEXTE<-as.factor(data$CONTEXTE)
data$MATRICE1<-as.factor(data$MATRICE1)
data$REGION<-as.factor(data$REGION)
data$SECTEUR<-as.factor(data$SECTEUR)
data$Annee<-as.factor(data$Annee)
data$mois<-as.factor(data$mois)
data$LIEUPVT<-as.factor(data$LIEUPVT)


raw_data<-read.csv("Variant_SISAGE.csv",sep=";")

prepare_input <- function (data, col_select=c(1:length(data)), id =1, date=0,
                           m=FALSE, a=TRUE) {
  duplicate <- which(duplicated(data[,id]))
  if (length(duplicate)!=0) {
    stop ( "duplicate ID : " , paste(data[duplicate,id], collapse = ", "))
  } else {
    print("no duplicate ID")
  }
  rownames(data) <- data[,id]
  col_select<-col_select[! col_select %in% date]
  clean_data <- data[,col_select]
  if (date!=0){
    if(a){
      clean_data$ANNEE <- substr(data[,date],7,10)
    }
    if(m){
      clean_data$MOIS <- substr(data[,date],4,5)
    } 
  }
  clean_data[] <- lapply(clean_data, factor) 
  print(colnames(clean_data))
  str(clean_data)
  clean_data
}

prepared_data <- prepare_input(raw_data,col_select = c(7,10,13,14,17), 
                               id=1,date=17,m=T,a=F)

assess_gower <- function(data, col_select=c(1:length(data)), id =1, date=0,
                         m=FALSE, a=TRUE, weights = rep.int(1, p), graph = TRUE){
  clean_data <- prepare_input(data,col_select=col_select, id =id, date=date,
                              m=m, a=a)
  p=length(clean_data)
  gower_dist <- daisy(clean_data, metric = "gower", weights = weights)
  aggl.clust.c <- hclust(gower_dist, method = "complete")
  if (graph) {
    plot(aggl.clust.c,main = "Agglomerative, complete linkages")
  }
  return(list(clean_data=clean_data,gower_dist=gower_dist, aggl.clust.c=aggl.clust.c))
}

gower <- assess_gower(raw_data,col_select = c(7,10,13,14,17), 
                      id=1,date=17,m=T,a=F)

 # 1. Distance

 
 gower.dist <- daisy(clean_data, metric = c("gower"))

 # 2. Clustering
 
 #divisive.clust <- diana(as.matrix(gower.dist), diss = TRUE, keep.diss = TRUE)
 #par(cex=0.6)
 #plot(divisive.clust, main = "Divisive") 

 aggl.clust.c <- hclust(gower.dist, method = "complete")
 plot(aggl.clust.c,main = "Agglomerative, complete linkages")
 
 # 3. Assess clustering
 
 # Cluster stats comes out as list while it is more convenient to look at it as a table
 # This code below will produce a dataframe with observations in columns and variables in row
 # Not quite tidy data, which will require a tweak for plotting, but I prefer this view as an output here as I find it more comprehensive 

 cstats.table <- function(dist, tree, k,l) {
   n=k/l
   clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between",
                     "wb.ratio","dunn2","avg.silwidth")
   clust.size <- c("cluster.size")
   stats.names <- c()
   row.clust <- c()
   output.stats <- matrix(ncol = n, nrow = length(clust.assess))
   cluster.sizes <- matrix(ncol =n, nrow = n)
   for(i in 1:n){
       row.clust[i] <- paste("Cluster-", i*l, " size")
   }
    for(i in 1:n){
     print(i*l)
     stats.names[i] <- paste("Test", i*l-1)
     
     for(j in seq_along(clust.assess)){
       output.stats[j, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i*l))[clust.assess])[j]
       
     }
     
     for(d in 1:n) {
       cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i*l))[clust.size])[d*l]
       dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
       cluster.sizes[d, i]
       
     }
   }
   output.stats.df <- data.frame(output.stats)
   cluster.sizes <- data.frame(cluster.sizes)
   cluster.sizes[is.na(cluster.sizes)] <- 0
   rows.all <- c(clust.assess, row.clust)
   # rownames(output.stats.df) <- clust.assess
   output <- rbind(output.stats.df, cluster.sizes)
   colnames(output) <- stats.names[1:n]
   rownames(output) <- rows.all
   is.num <- sapply(output, is.numeric)
   output[is.num] <- lapply(output[is.num], round, 2)
   output
 }

 
 #stats.df.divisive <- cstats.table(gower.dist, divisive.clust, 15)
 #stats.df.divisive 
 
 
 # stats.df.aggl <-cstats.table(gower.dist, aggl.clust.c, 90) #complete linkages looks like the most balanced approach
 # stats.df.aggl
 
 
 data_fig_5<-data.frame(t(cstats.table(gower$gower_dist, gower$aggl.clust.c, 350,10)))
 write.table(data_fig_5,"data_fig_350_5mtd.csv", sep = ";", row.names = T)
 # Elbow
 
 # Agglomerative clustering,provides a more ambiguous picture
 ggplot(data = data_fig_5, 
        aes(x=cluster.number, y=within.cluster.ss)) + 
   geom_point()+
   geom_line()+
   ggtitle("Agglomerative clustering") +
   labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
   theme(plot.title = element_text(hjust = 0.5)) 
 
 # Silhouette
 
 ggplot(data = data_fig_5, 
        aes(x=cluster.number, y=avg.silwidth)) + 
   geom_point()+
   geom_line()+
   ggtitle("Agglomerative clustering") +
   labs(x = "Num.of clusters", y = "Average silhouette width") +
   theme(plot.title = element_text(hjust = 0.5)) 
 
 # 4. Visualisation
 dendro <- as.dendrogram(gower$aggl.clust.c)
 dendro.col <- dendro %>%
   set("branches_k_color", k = 60, value =rainbow(n=60)    ) %>%
   set("branches_lwd", 0.6) %>%
   set("labels_colors", 
       value = c("darkslategray")) %>% 
   set("labels_cex", 0.5)
 ggd1 <- as.ggdend(dendro.col)
 ggplot(ggd1, theme = theme_minimal()) +
   labs(x = "Num. observations", y = "Height", title = "Dendrogram, k = 60")
 
 # Radial plot looks less cluttered (and cooler)
 ggplot(ggd1, labels = T) + 
   scale_y_reverse(expand = c(0.2, 0)) +
   coord_polar(theta="x")
 
 
 
 clust.num <- cutree(gower$aggl.clust.c, k = 60) 
 data.cl <- cbind(prepared_data, clust.num)
write.table(file="sisage_60.csv",data.cl,sep = ";", row.names = TRUE)
 clusplot(data.cl, clust.num, 
          color=TRUE, shade=TRUE, labels=0, lines=0, 
          main = "Customer clusters (k=60)", 
          cex = 0.3)
 
 clust.num <- cutree(gower$aggl.clust.c, k = 100) 
 data.cl <- cbind(prepared_data, clust.num)
 write.table(file="sisage_100.csv",data.cl,sep = ";", row.names = TRUE)
 clusplot(data.cl, clust.num, 
          color=TRUE, shade=TRUE, labels=0, lines=0, 
          main = "Customer clusters (k=100)", 
          cex = 0.3)

 clust.num <- cutree(gower$aggl.clust.c, k = 310) 
 data.cl <- cbind(prepared_data, clust.num)
 write.table(file="sisage_310.csv",data.cl,sep = ";", row.names = TRUE)
 clusplot(data.cl, clust.num, 
          color=TRUE, shade=TRUE, labels=0, lines=0, 
          main = "Customer clusters (k=310)", 
          cex = 0.3)
 
 
 
 ###matrice upset
 
 col.names <- NULL
 metadata <- colnames(prepared_data)
 vect_all <- NULL
 for (i in 1:ncol(prepared_data)){
   vect <- unique(prepared_data[,i])
   vect_all <- c(vect_all, as.vector(vect))
   for (j in 1:length(vect)){
     col.names <- c(col.names,paste(metadata[i],"_",vect[j]))
   }
 }
 
 upset_data <- matrix(0,nrow(prepared_data),length(vect_all))
 rownames(upset_data) <- rownames(prepared_data)
 colnames(upset_data) <- col.names
 for (j in 1:ncol(prepared_data)){
   for (i in 1:nrow(prepared_data)){  
     for (k in 1:length(col.names)){
       if (is.na(prepared_data[i,j])==F & is.na(vect_all[k])==F){
         if(prepared_data[i,j]==vect_all[k]){
           upset_data[i,k] <- 1
         } else {
         }
       }
     }
   }
 }
 
 upset_data <- as.data.frame(upset_data)
 upset(upset_data, nsets = 10,)
 ###
 
 