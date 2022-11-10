# SiSAGE strain selection


# 0. Libraries  
library(cluster) 
library(fpc)
library(ggplot2)
library(dplyr)
library(dendextend)
library(UpSetR)

raw_data<-read.csv("Variant_SISAGE.csv",sep=";")

# 1. Distance and Clustering

## Fonction pr?paration tableau de donn?es :
# data = tableau de donn?es complet 
# col_select = colonnes ? selectionner pour le clustering (ex : c(7,10,13,14,17))
# id = colonne correspondant ? l'identifiant de la souche 
# date = colonne correspondant ? la date
# m et a = TRUE si volont? de cr?er les colonnes MOIS et ANNEE sinon FALSE

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
  clean_data[is.na(clean_data)] <- "non spécifié"
  clean_data[] <- lapply(clean_data, factor) 
  print(colnames(clean_data))
  str(clean_data)
  clean_data
}

prepared_data <- prepare_input(raw_data,col_select = c(7,10,13,14,17,24), 
                                id=1,date=17,m=TRUE,a=TRUE)
## Fonction assess_gower : pr?pare le tableau de donn?e, calcul la matrice de distance et fait le clustering
# reprend les entr?es de prepare_input
# weights = Poids des colonnes (ex:  c(1,2,1,1,1) pour donner un poids de 2 ? la 2e colonne)
# graph = TRUE pour avoir le dendrogramme, FAUX sinon

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

gower <- assess_gower(raw_data,col_select = c(7,10,13,14,17,24), 
                      id=1,date=17,m=TRUE,a=TRUE, graph = T)


 # 2. Assess clustering
 
 # Cluster stats comes out as list while it is more convenient to look at it as a table
 # This code below will produce a dataframe with observations in columns and variables in row
 # Not quite tidy data, which will require a tweak for plotting, but I prefer this view as an output here as I find it more comprehensive 

# dist = matrice de distance de gower (gower$gower_dist)
# tree = clustering r?alis? sur la matrice de distance de goxer (gower$aggl.clust.c)
# k = nombre de cluster ? tester 
# l = pas pour tester les cluster 
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

 
 data_fig<-data.frame(t(cstats.table(dist = gower$gower_dist, tree = gower$aggl.clust.c, k = 50,l = 10)))

 # Elbow
 
 # Agglomerative clustering,provides a more ambiguous picture
 ggplot(data = data_fig, 
        aes(x=cluster.number, y=within.cluster.ss)) + 
   geom_point()+
   geom_line()+
   ggtitle("Agglomerative clustering") +
   labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
   theme(plot.title = element_text(hjust = 0.5)) 
 
 # Silhouette
 
 ggplot(data = data_fig, 
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
 
 
 
 clust.num <- cutree(gower$aggl.clust.c, k = 20) 
 data.cl <- cbind(prepared_data, clust.num)
write.table(file="sisage_20_matrice2.csv",data.cl,sep = ";", row.names = TRUE)
 clusplot(data.cl, clust.num, 
          color=TRUE, shade=TRUE, labels=0, lines=0, 
          main = "Customer clusters (k=20)", 
          cex = 0.3)
 
 clust.num <- cutree(gower$aggl.clust.c, k = 150) 
 data.cl <- cbind(prepared_data, clust.num)
 write.table(file="sisage_150_matrice2.csv",data.cl,sep = ";", row.names = TRUE)
 clusplot(data.cl, clust.num, 
          color=TRUE, shade=TRUE, labels=0, lines=0, 
          main = "Customer clusters (k=150)", 
          cex = 0.3)

 clust.num <- cutree(gower$aggl.clust.c, k = 310) 
 data.cl <- cbind(prepared_data, clust.num)
 write.table(file="sisage_310.csv",data.cl,sep = ";", row.names = TRUE)
 clusplot(data.cl, clust.num, 
          color=TRUE, shade=TRUE, labels=0, lines=0, 
          main = "Customer clusters (k=310)", 
          cex = 0.3)
 
 
 
 ###matrice upset
prepare_upset<- function(prepared_data){
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
        if (is.na(prepared_data[i,j])==F){
          if(prepared_data[i,j]==vect_all[k]){
            upset_data[i,k] <- 1
          }
        }
      }
    }
  }
  upset_data <- as.data.frame(upset_data)
  upset_data
}
upset_data <- prepare_upset(prepared_data)
upset(upset_data, nsets = 55, order.by = "freq", mb.ratio = c(0.30,0.70))
 ###
 #https://veroniquetremblay.github.io/analyse_de_donnees_et_apprentissage_statistique_en_R/classification-non-supervisee.html
#sm() = Simple Matching Coefficient (SM)
 prox.sm <- sm(upset_data)
 clust.sm <- hclust(prox.sm, method = "complete")
plot(clust.sm,main = "Agglomerative, SM")
#Comparaison dendrogrram
dend_gower <- as.dendrogram(gower$aggl.clust.c)
dend_sm <- as.dendrogram(clust.sm)
dend_list <- dendlist(dend_gower, dend_sm)
cor.dendlist(dend_list)
dendlist(dend_gower, dend_sm) %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram()                       # Draw the two dendrograms



profil_data <-  prepared_data
profil_data$profils <- paste(prepared_data$CONTEXTE,"_",prepared_data$MATRICE1,
                             "_",prepared_data$REGION,"_",prepared_data$SECTEUR,
                             "_",prepared_data$LIEUPVT,"_",prepared_data$ANNEE,
                             "_",prepared_data$MOIS) 
profils <- table(profil_data$profils)
View(profils)