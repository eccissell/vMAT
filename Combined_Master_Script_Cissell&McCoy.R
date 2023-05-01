#Combined script for analyses presented in Cissell & McCoy: Top-heavy trophic structure within benthic viral dark matter

#SECTION 1 - VIRAL SEQUENCE CURATION
#Load libraries
library("tidyverse")
library(tidyr)
library(ggplot2)
library(conflicted)
library(forcats)
#Read in tsvs
setwd("/Users/cissell/Desktop/VIRSORTER2_Diel")
viral_tsv=dir(patter=".tsv")
n = length(viral_tsv)
viral_tsv_list=vector("list",n)
for(i in 1:n) viral_tsv_list[[i]] = read.table(viral_tsv[i],header=T)

#Merge by contig id
D1=merge(viral_tsv_list[[1]],viral_tsv_list[[2]],by="contig_id",all=T)
D2=merge(viral_tsv_list[[3]],viral_tsv_list[[4]],by="contig_id",all=T)
D3=merge(viral_tsv_list[[5]],viral_tsv_list[[6]],by="contig_id",all=T)
D4=merge(viral_tsv_list[[7]],viral_tsv_list[[8]],by="contig_id",all=T)
R1=merge(viral_tsv_list[[9]],viral_tsv_list[[10]],by="contig_id",all=T)
R2=merge(viral_tsv_list[[11]],viral_tsv_list[[12]],by="contig_id",all=T)
R3=merge(viral_tsv_list[[13]],viral_tsv_list[[14]],by="contig_id",all=T)
R4=merge(viral_tsv_list[[15]],viral_tsv_list[[16]],by="contig_id",all=T)

#Make single dataframe
viral_merged=bind_rows(D1,D2,D3,D4,R1,R2,R3,R4)

#Now apply manual curation
#Keep1:viral_gene>0
keep1=viral_merged%>%
  subset(viral_genes>0)
keep1_remain=viral_merged%>%
  subset(!viral_genes>0)
#Keep2_1:host_gene=0
keep2_1=keep1_remain%>%
  subset(host_genes==0)
keep2_1_remain=keep1_remain%>%
  subset(!host_genes==0)
#Keep2_2:score>=0.95
keep2_2=keep2_1_remain%>%
  subset(max_score>=0.95)
keep2_2_remain=keep2_1_remain%>%
  subset(!max_score>=0.95)
#Keep2_3:hallmark>2
keep2_3=keep2_2_remain%>%
  subset(hallmark>2)
keep2_3_remain=keep2_2_remain%>%
  subset(!hallmark>2)
#Manual:Host<=1 & length >=10000
keepman=keep2_3_remain%>%
  subset(host_genes<=1 & length>=10000)
#Man merge
manmerge=rbind(keepman,keep2_3)
#Now screen those contigs in manmerge manually
#remove
##D2_2973455||full
manmerge_cur=manmerge%>%
  subset(!contig_id=="D2_2973455||full")

#Now merge all the kept phages
viral_curated=bind_rows(keep1,keep2_1,keep2_2,manmerge_cur)

#Ensure deduped
viral_dedup=viral_curated[!duplicated(viral_curated),] #Fine

#Export combined as a tsv
write.table(viral_curated[,1], file="/Users/cissell/Desktop/VIRSORTER2_Diel/curated_ids.tsv", quote=F, sep='\t', col.names=NA)

#Export prophage list
viral_prophage=viral_curated%>%
  subset(provirus=="Yes")
write.table(viral_prophage[,1], file="/Users/cissell/Desktop/VIRSORTER2_Diel/prophages_list.tsv", quote=F, sep='\t', col.names=NA)

#Export genes<4 list
viral_vibrant_cutoff=viral_curated%>%
  subset(!total_genes>=4)

##Now we have depreplicated at 95% within mat
dedup_ids=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/IDS/cat_dedup_headers.tsv",header=T)

#Remove IDS after || in curated dataframe for merging
#Just do it in Excel becuase it is easier.
write.table(viral_curated, file="/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/curated.tsv", quote=F, sep='\t', col.names=NA)
#Now read back in
viral_curated_simp=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/curated.tsv",header=T)

#Now merge to subset by clustered at 95%
viral_clust_curate=merge(viral_curated_simp,dedup_ids,by="contig_id")

#Let's clean up lengths for pro vs non pro regions
#Start by splitting into 2 dataframes, pro vs non-pro
viral_clust_curate_pro=viral_clust_curate%>%
  subset(provirus=="Yes")
viral_clust_curate_non=viral_clust_curate%>%
  subset(provirus=="No")

#clean up columns from non
viral_clust_curate_non_clean=viral_clust_curate_non%>%
  select(-c(proviral_length,
            host_length,
            region_types,
            region_lengths,
            region_coords_bp,
            region_coords_genes,
            region_viral_genes,
            region_host_genes,
            length))
#clean up columns in pro and rename length column
viral_clust_curate_pro_clean=viral_clust_curate_pro%>%
  select(-c(host_length,
            region_types,
            region_lengths,
            region_coords_bp,
            region_coords_genes,
            region_viral_genes,
            region_host_genes,
            length,
            contig_length)) %>%
  rename(contig_length=proviral_length)
#recombine
viral_clust_curate_clean=rbind(viral_clust_curate_non_clean,viral_clust_curate_pro_clean)

###SECTION 2: VIRAL TAXONOMIC PREDICTIONS AND HIVE PLOT
##Create Hive Plot Visualizations for viral network
#Load libraries
library(igraph)
library(ggraph)
library(dplyr)
set.seed(999)

#Read in graph edges that I converted to a tsv
diel_viral_network_edge=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/VCONTACT/c1.txt",
                                   header=T)
#Simplify redundancy
diel_viral_network_simp=diel_viral_network_edge[!duplicated(data.frame(t(apply(diel_viral_network_edge[1:2], 1, sort)), diel_viral_network_edge$Weight)),]

#Read in graph vertices and metadata
diel_viral_network_vert=read.csv("/Users/cissell/Desktop/VIRSORTER2_Diel/VCONTACT/genome_by_genome_overview.csv")
#Add in classifications from PhaGCN for my phages (requires outputs from Manual_Viral_Ref_script)
diel_phagcn_hive=rename(diel_phagcn,vertex=contig_id)
diel_phagcn_hive=diel_phagcn_hive%>%
  select(-c(idx,score))
vertex_merge=merge(diel_viral_network_vert,diel_phagcn_hive,by="vertex",all.x=T)
str(vertex_merge)
#Coalesce to create new merged taxonomy
vertex_merge=vertex_merge%>%
  mutate(mergetax=coalesce(refFamily,prediction))


######
#Extract clusters among my samples and Ref for vContact taxonomy to add to PhaGcn tax
######
#First we need to extract the 2 column we care about to create a list of common VCs
vcon_tax_1=diel_viral_network_vert%>%
  select(c("VC","Source"))
#Now distinct
vcon_tax_2=vcon_tax_1%>%
  distinct(VC,Source)
str(vcon_tax_2)

#Now extract duplicated
vcon_tax_3=vcon_tax_2%>%
  subset(duplicated(VC)) %>%
  select(-Source)
vcon_tax_3$VC=as.vector(vcon_tax_3$VC)

#Export
write.table(vcon_tax_3, file="/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/vcon_tax_clusters.tsv", quote=F, sep='\t', col.names=NA)
#After cleaning up read back in
vcon_tax_3=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/vcon_tax_clusters.tsv",
                      header=T)

#Now extract those clusters from original and manually inspect
vcon_tax_4=merge(vertex_merge,vcon_tax_3,by="VC",all.y=T)

#Manually create excel sheet because there are so few
write.table(vcon_tax_4, file="/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/vcon_tax_clusters_merged.tsv", quote=F, sep='\t', col.names=NA)
#Read in that manually created sheet
vcon_tax_5=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/vcon_tax_clusters_merged_final.txt",header=T)
vcon_tax_5=vcon_tax_5%>%
  select(-c(VC,PhaGCN_TAX))

#and merge to create new column with vcon taxonomy
vertex_merge_and_vcon=merge(vertex_merge,vcon_tax_5,by="vertex",all.x=T)
######
#Now use clusters just among mats where at least one member has an assigned taxonomy to further resolve taxonomy
######
mat_tax_1=vertex_merge%>%
  select(c("VC","Source","prediction","vertex")) %>%
  subset(Source=="Mat")
mat_tax_1=rename(mat_tax_1,contig_id=vertex)
mat_tax_1=mat_tax_1[!mat_tax_1$VC=="",]
#add in PhaGCN score
mat_tax_2=merge(mat_tax_1,diel_phagcn,by="contig_id",all.x=T)
#subset out lower than 0.75 score, but first make NA score 2 so it doesnt get filterd
mat_tax_2$score=as.character(mat_tax_2$score)
mat_tax_2$score[is.na(mat_tax_2$score)]=2
mat_tax_3=mat_tax_2%>%
  subset(score >=0.75)
#Now distinct
mat_tax_4=mat_tax_3[duplicated(mat_tax_3$VC),]

#manually go through this (more of a bummer than the last for sure, but doable)
write.table(mat_tax_4, file="/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/mat_tax_clusters_merged.tsv", quote=F, sep='\t', col.names=NA)

#read back in
mat_clust_tax=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/mat_tax_clusters_merged_simp_final.txt", header = T)

#merge back in
vertex_merge_and_vcon_and_matclust=merge(vertex_merge_and_vcon,mat_clust_tax,
                                         by="vertex",all.x=T)

#Coalesce into mergetax column predictions from vcon and mat clusts
vertex_merge_and_vcon_and_matclust=vertex_merge_and_vcon_and_matclust%>%
  mutate(mergetax=coalesce(mergetax,VCON_TAX))
vertex_merge_and_vcon_and_matclust=vertex_merge_and_vcon_and_matclust%>%
  mutate(mergetax=coalesce(mergetax,Mat_clust_tax_predict))

#Replace NA with "unknown'
vertex_merge_and_vcon_and_matclust$mergetax=as.character(vertex_merge_and_vcon_and_matclust$mergetax)
vertex_merge_and_vcon_and_matclust$mergetax=vertex_merge_and_vcon_and_matclust$mergetax %>%
  replace_na('Unknown')
vertex_merge_and_vcon_and_matclust$mergetax=as.factor(vertex_merge_and_vcon_and_matclust$mergetax)

#now lets subset to only color code major 3 families, lumping everything else as "other"
#gsub isn't working right so just export and import
write.table(vertex_merge_and_vcon_and_matclust, file="/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/vertex_merge_and_vcon_and_matclust.tsv", quote=F, sep='\t', col.names=NA)

vertex_fams_reduced=read.csv("/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/vertex_fams_reduced.csv")

#now lets subset to only color code major 3 families, lumping everything else as "other"
#gsub isn't working right so just export and import
write.table(vertex_merge_and_vcon_and_matclust, file="/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/vertex_merge_and_vcon_and_matclust.tsv", quote=F, sep='\t', col.names=NA)

vertex_fams_reduced=read.csv("/Users/cissell/Desktop/VIRSORTER2_Diel/Curate/vertex_fams_reduced.csv")
diel_viral_network_simp[which(!diel_viral_network_simp$from %in% vertex_merge_and_vcon_and_matclust$vertex),]

#Create a iGraph object from subset frame
viral_net_subs=graph_from_data_frame(diel_viral_network_simp,directed=F,vertices=vertex_fams_reduced)
print(viral_net_subs, e=TRUE, v=TRUE)

write.table(diel_viral_network_simp, file="/Users/cissell/Desktop/diel_viral_network_simp.tsv", quote=F, sep='\t', col.names=NA)
write.table(vertex_fams_reduced, file="/Users/cissell/Desktop/vertex_fams_reduced.tsv", quote=F, sep='\t', col.names=NA)
diel_viral_network_simp = read.table("/Users/cissell/Desktop/diel_viral_network_simp.tsv",header=T)
vertex_fams_reduced = read.xls("/Users/cissell/Desktop/vertex_fams_reduced.xlsx")

#Let's add a column for Mat spatial identity
vertex_fams_reduced = vertex_fams_reduced %>%
  mutate(Mat=case_when(
    str_detect(vertex,"D1_") ~"D1",
    str_detect(vertex,"D2_") ~"D2",
    str_detect(vertex,"D3_") ~"D3",
    str_detect(vertex,"D4_") ~"D4"
  )) %>%
  mutate(Mat = ifelse(is.na(Mat),"Ref",Mat))

#Get degree for sorting
V(viral_net_subs)$degree = degree(viral_net_subs,mode='all')
E(viral_net_subs)$betweeness=estimate_edge_betweenness(viral_net_subs,cutoff=0)
median(E(viral_net_subs)$betweeness)
E(viral_net_subs)$betweeness=ifelse(E(viral_net_subs)$betweeness <median(E(viral_net_subs)$betweeness),"Below","Above")
mean(E(viral_net_subs)$Weight)
E(viral_net_subs)$Weight=ifelse(E(viral_net_subs)$Weight <median(E(viral_net_subs)$Weight),"Low","High")
#Plot for visualizing dark matter and spread across diversity of ref genomes
dark_matter_net=ggraph(viral_net_subs, 'hive',
                       axis=Source,
                       sort.by=degree) +
  geom_axis_hive()+
  geom_edge_hive(alpha=0.04) +
  geom_node_point(aes(colour=mergetax),alpha=0.9) +
  theme_classic() +
  scale_color_manual(values=c("#382954",
                              "#aee3d2",
                              "#fa82a7",
                              "#3497a8",
                              "#8b9e90"))+
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  labs(colour="Family")
dark_matter_net
ggsave("/Users/cissell/Desktop/dm_net.svg",dark_matter_net)

#Plot for visualizing 'unknown' tax gene associations
unknown_net=ggraph(viral_net_subs, 'hive',
                   axis=Source,
                   sort.by=mergetax) +
  geom_axis_hive()+
  geom_edge_hive(aes(color=betweeness),alpha=0.04) +
  geom_node_point(aes(colour=mergetax),alpha=0.9) +
  theme_classic() +
  scale_color_manual(values=c("#382954",
                              "#aee3d2",
                              "#fa82a7",
                              "#3497a8",
                              "#8b9e90"))+
  new_scale_color()+
  scale_edge_color_manual(values=c("#77a0a1",
                                   "#d4485e"))+
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  labs(colour="Family")
unknown_net
ggsave("/Users/cissell/Desktop/uk_net.svg",unknown_net)

#Classic network visualization

classic_net_source=ggraph(viral_net_subs,
       layout='stress') +
  geom_edge_link(show.legend=F,
                 alpha=0.7) +
  geom_node_point(aes(color=Source),
                  alpha=0.8) +
  theme_classic() +
  scale_color_manual(values=c("#382954",
                              "#aee3d2",
                              "#fa82a7",
                              "#3497a8",
                              "#8b9e90")) +
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  labs(colour="Source")
classic_net_source
ggsave("/Users/cissell/Desktop/classic_net_source.svg",classic_net_source)

classic_net_tax=ggraph(viral_net_subs,
                          layout='stress') +
  geom_edge_link(show.legend=F,
                 alpha=0.7) +
  geom_node_point(aes(color=mergetax),
                  alpha=0.8) +
  theme_classic() +
  scale_color_manual(values=c("#382954",
                              "#aee3d2",
                              "#fa82a7",
                              "#3497a8",
                              "#8b9e90")) +
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  labs(colour="Source")
ggsave("/Users/cissell/Desktop/classic_net_tax.svg",classic_net_tax)

classic_net_mat=ggraph(viral_net_subs,
                          layout='stress') +
  geom_edge_link(show.legend=F,
                 alpha=0.7) +
  geom_node_point(aes(color=Mat),
                  alpha=0.8) +
  theme_classic() +
  scale_color_manual(values=c("#382954",
                              "#aee3d2",
                              "#fa82a7",
                              "#3497a8",
                              "#8b9e90")) +
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  labs(colour="Mat")
classic_net_mat
ggsave("/Users/cissell/Desktop/classic_net_mat.svg",classic_net_mat)


#Get network properties
gorder(viral_net_subs)
#Density
edge_density(viral_net_subs,loops=F) #0.00051
#Diameter
diameter(viral_net_subs,
         directed=F,
         weights=NA) #25

#Mean Degree and Distribution
deg = degree(viral_net_subs,
             mode="all")
mean(deg) #7.85
hist(deg)
dev.print(png,fil="/Users/cissell/Documents/Manuscripts/Diel_Cycle_Phages_Only/FIGURES/deghist_fullref_mat.png",type="quartz",antialias="default",width=7,
          height=5,units="in",res=1300)
#Cumulative degree distribution
deg.dist = degree_distribution(viral_net_subs,
                               cumulative=T,
                               mode="all")
full_degdist=plot(x=0:max(deg),
     y=1-deg.dist,
     pch=19,
     cex=1.2,
     col="orange",
     xlab="Degree",
     ylab="Cumulative Frequency")
full_degdist
dev.print(png,fil="/Users/cissell/Documents/Manuscripts/Diel_Cycle_Phages_Only/FIGURES/degdist_fullref_mat.png",type="quartz",antialias="default",width=7,
          height=5,units="in",res=1300)
#Community Analysis
ceb = cluster_edge_betweenness(viral_net_subs)
#Assign communities as vertex properties
V(viral_net_subs)$community = membership(ceb)

algorithm(ceb) #edge_between
modularity(ceb) #0.942
length(ceb) #7754
mean(sizes(ceb))
crossing(ceb,viral_net_subs)
ggraph(ceb,layout="dendrogram", circular = T) +
  geom_edge_diagonal()
plot(viral_net_subs)
#assortativity
assortativity_nominal(viral_net_subs,
                      as.factor(V(viral_net_subs)$Source),
                      directed=F) #0.806
assortativity_nominal(viral_net_subs,
                      as.factor(V(viral_net_subs)$refFamily),
                      directed=F) #0.803
assortativity_nominal(viral_net_subs,
                      as.factor(V(viral_net_subs)$Mat),
                      directed=F) #0.418

#After review, more traditional
library(igraph)
library(ggraph)
library(dplyr)
set.seed(999)
library(gdata)
library(stringr)
library(ggnewscale)

#Get descriptive stats like from other network for this network
descript = vertex_fams_reduced %>%
  subset(VCStatus == "Clustered") %>%
  mutate(VC = as.factor(VC)) %>%#2011 unique VCs
  group_by(VC) %>%
  mutate(Source = as.factor(Source)) %>%
  summarise(count = n_distinct(Source)) %>%
  mutate(count = as.factor(count))
count(descript,count)
#1      1983
#2        28
28/2011 #0.0139
1983/2011 #0.986

#Now get total counts (vOTU level) of those that clustered
descript2 = vertex_fams_reduced %>%
  subset(VCStatus == "Clustered") %>%
  mutate(VC = as.factor(VC)) %>%#2011 unique VCs
  group_by(VC) %>%
  mutate(Source = as.factor(Source)) %>%
  summarise(count = n_distinct(Source)) %>%
  mutate(count = as.factor(count)) %>%
  subset(count == 2)

clustvotu = vertex_fams_reduced %>%
  subset(VCStatus == "Clustered") %>%
  subset(Source == "Mat") %>%
  mutate(VC = as.factor(VC)) %>%
  subset(VC %in% descript2$VC)
 #70 total vOTUs

#We need to get counts of vOTUs from vMAT that share edges (not clustering necessarily) with RefSeq nodes
#We can also plot here a subset network that only includes those edges between differenct sources for a more zoomed in view
CrossEdges = which(sapply(E(viral_net_subs), function(x) {
  length(unique(V(viral_net_subs)[ends(viral_net_subs, x)]$Source)) != 1 }))
SUB = subgraph.edges(viral_net_subs, CrossEdges)

sub_net_ref=ggraph(SUB,
                      layout='stress') +
  geom_edge_link(show.legend=F,
                 alpha=0.7) +
  geom_node_point(aes(color=refFamily),
                  alpha=0.8) +
  theme_classic() +
  scale_color_manual(values=c("#382954",
                              "#aee3d2",
                              "#fa82a7",
                              "#3497a8",
                              "#8b9e90")) +
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  labs(colour="RefSeq Family")
sub_net_ref
ggsave("/Users/cissell/Desktop/sub_net_ref.svg",sub_net_ref)

#Now get count of Mat vertices
vertex_attr(SUB)
gorder(SUB)
votussubcount=(V(SUB)$Source=="Mat")
votussubcount = as.data.frame(votussubcount)
votussubcount = votussubcount %>%
  mutate(votussubcount = as.factor(votussubcount))
count(votussubcount,votussubcount)
#707 vOTUs have a significant edge to a RefSeq virus

#Community analysis on subset network
ceb2 = cluster_edge_betweenness(SUB)
modularity(ceb2) #0.88
length(ceb2)
mean(sizes(ceb2))
max(sizes(ceb2))
crossing(ceb2,SUB)
dendPlot(ceb2,mode="hclust",
         cex=0.3)
V(SUB)$community = ceb2$membership

sub_net_ref=ggraph(SUB,
                   layout='stress') +
  geom_edge_link(show.legend=F,
                 alpha=0.7) +
  geom_node_point(aes(color=refFamily),
                  alpha=0.8) +
  theme_classic() +
  scale_color_manual(values=c("#382954",
                              "#aee3d2",
                              "#fa82a7",
                              "#3497a8",
                              "#8b9e90")) +
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  labs(colour="RefSeq Family")
sub_net_ref

#Edge density and diameter of subset network
edge_density(SUB,loops=F) #0.00518
#Diameter
diameter(SUB,
         directed=F,
         weights=NA) #19

#Mean Degree and Distribution
deg_sub = degree(SUB,
             mode="all")
mean(deg_sub) #7.23
hist(deg_sub)
dev.print(png,fil="/Users/cissell/Documents/Manuscripts/Diel_Cycle_Phages_Only/FIGURES/deghist_subref_mat.png",type="quartz",antialias="default",width=7,
          height=5,units="in",res=1300)
#Cumulative degree distribution
deg.dist_sub = degree_distribution(SUB,
                               cumulative=T,
                               mode="all")
sub_degdist=plot(x=0:max(deg_sub),
                  y=1-deg.dist_sub,
                  pch=19,
                  cex=1.2,
                  col="orange",
                  xlab="Degree",
                  ylab="Cumulative Frequency")
dev.print(png,fil="/Users/cissell/Documents/Manuscripts/Diel_Cycle_Phages_Only/FIGURES/degdist_subref_mat.png",type="quartz",antialias="default",width=7,
          height=5,units="in",res=1300)

#After review, more traditional
####Diel Hive Plot for environmental context with permafrost and GOV2.0
library(igraph)
library(ggraph)
library(dplyr)
library(tidyr)
set.seed(999)

#Read in graph edges
diel_enviro_network_edge=read.table("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\c1_ntw.txt", header=F)
str(diel_enviro_network_edge)
colnames(diel_enviro_network_edge)=c("from","to","Weight")

#Simplify Redundancy
diel_enviro_network_edge_simp=diel_enviro_network_edge[!duplicated(data.frame(t(apply(diel_enviro_network_edge[1:2],1,sort)),diel_enviro_network_edge$Weight)),]

#Read in graph vertices and metadata
diel_enviro_network_vert=read.csv("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\genome_by_genome_overview.csv")

#Add in classifications to my phages? Unsure if this is really necessary. Hold off for now and revisist later on.

##Create an iGraph object from edge list with weights
viral_env_net=graph_from_data_frame(diel_enviro_network_edge_simp,directed=F,vertices=diel_enviro_network_vert)
print(viral_env_net,e=T,v=T)
V(viral_env_net)$degree = degree(viral_env_net,mode='all')

#png("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\hive.png",
#    width=6,
#    height=6,
#    units="in",
#    res=1080)
ggraph(viral_env_net,'hive',
       axis=Source,
       sort.by=degree) +
  geom_axis_hive()+
  geom_edge_hive(alpha=.04) +
  geom_node_point(aes(colour=Source),
                  alpha=0.9) +
  theme_classic()+
  scale_color_manual(values=c("#9aacb8",
                              "#3c455c",
                              "#f27349"))+
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())
#dev.off()
ggsave("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\hive.png",
       device="png",
       width=12,
       height=12,
       units="in",
       dpi=600)

#Now let's do a more traditional network visualization following reviewer comments
#First, we will try an LGL algorithm
ggraph(viral_env_net,
       layout = 'lgl') +
  geom_edge_link(show.legend = F) +
  geom_node_point(aes(colour=Source),
                  alpha=0.9) +
  theme_classic()+
  scale_color_manual(values=c("#9aacb8",
                              "#3c455c",
                              "#f27349"))+
  theme(axis.ticks=element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())

#Also, we need to get other descriptive statistics about the network topology (general network properties) to include in the supp and somewhat in the main results text
#Density
edge_density(viral_env_net,loops=F) #0.00029

#Diameter
diameter(viral_env_net,directed=F,weights=NA) #34

#Mean degree and distribution
degree = degree(viral_env_net, mode="all")
mean(degree) #60.09
hist(degree)
ggsave("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\govdeghist.png",
       device="png",
       width=6,
       height=6,
       units="in",
       dpi=600)
#Cumulative degree distribution
deg.dist = degree_distribution(viral_env_net, cumulative=T,mode="all")
plot( x=0:max(degree), y=1-deg.dist, pch=19, cex=1.2, col="orange",

      xlab="Degree", ylab="Cumulative Frequency")
ggsave("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\govdegdist.png",
       device="png",
       width=6,
       height=6,
       units="in",
       dpi=600)
-#Community analysis
ceb = cluster_edge_betweenness(viral_env_net)

#Assortativity coefficient
assortativity_nominal(viral_env_net,
                      as.factor(V(viral_env_net)$Source),
                      directed=F) #0.45




##Now extract summary stats about the clustering patterns, and make bar plots
#Overall barplot of clustering vs not irrespective of source
#Simplify Overlap factor levels in Excel and read in
diel_enviro_overlap_simp=read.csv("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\genome_by_genome_overview_overlap_simp.csv")
diel_enviro_net_stats_general=diel_enviro_overlap_simp%>%
  group_by(VC.Status)%>%
  summarise(n=n())%>%
  mutate(freq=100*(n/sum(n))) %>%
  ungroup()%>%
  mutate(set="overall")

#Now do overall for just vmat
diel_enviro_net_stats_vmatgeneral=diel_enviro_overlap_simp%>%
  subset(Source=='vMAT')%>%
  group_by(VC.Status)%>%
  summarise(n=n())%>%
  mutate(freq=100*(n/sum(n))) %>%
  ungroup()%>%
  mutate(set="vmatoverall")

#Now do overall for just vFROST
diel_enviro_net_stats_vfrostgeneral=diel_enviro_overlap_simp%>%
  subset(Source=='Permafrost')%>%
  group_by(VC.Status)%>%
  summarise(n=n())%>%
  mutate(freq=100*(n/sum(n))) %>%
  ungroup()%>%
  mutate(set="vfrostoverall")

#Now do overall for just GOV2.0
diel_enviro_net_stats_govgeneral=diel_enviro_overlap_simp%>%
  subset(Source=='GOV2.0')%>%
  group_by(VC.Status)%>%
  summarise(n=n())%>%
  mutate(freq=100*(n/sum(n))) %>%
  ungroup()%>%
  mutate(set="govoverall")

#Now do clustering of vMAT with GOV and vfrost and other vMAT
#make dataframes with each cluster type, and then find the merges
#Mat only
diel_enviro_net_stats_vmatonly=diel_enviro_overlap_simp%>%
  subset(VC.Status=="Clustered")%>%
  subset(Source=='vMAT')%>%
  subset(!duplicated(VC.Subcluster))
#GOV only
diel_enviro_net_stats_govonly=diel_enviro_overlap_simp%>%
  subset(VC.Status=="Clustered")%>%
  subset(Source=='GOV2.0')%>%
  subset(!duplicated(VC.Subcluster))
#Perm only
diel_enviro_net_stats_permonly=diel_enviro_overlap_simp%>%
  subset(VC.Status=="Clustered")%>%
  subset(Source=='Permafrost')%>%
  subset(!duplicated(VC.Subcluster))
#Intersect
mat_sea_clustersect=merge(diel_enviro_net_stats_vmatonly,diel_enviro_net_stats_govonly,by="VC.Subcluster")#258
mat_perm_clustersect=merge(diel_enviro_net_stats_vmatonly,diel_enviro_net_stats_permonly,by="VC.Subcluster")#32
tri_clustersect=merge(mat_perm_clustersect,mat_sea_clustersect,by="VC.Subcluster")#11
#So counts are MAt-GOV: 258-11 = 247 ; 14.6
#MAT-PERM: 32 - 11 = 21 ; 1.24
#Tri = 11 ; 0.65
#MAT-MAT = 1689 - (247+21+11) = 1410 ; 83.5

#Following review, do this also for vFROST and GOV
# vFROSTIntersect
frost_gov_clustersect=merge(diel_enviro_net_stats_permonly,diel_enviro_net_stats_govonly,by="VC.Subcluster")#112
#So counts are frost-GOV: 112-11 = 101; 29.4
#frost-mat: 32 - 11 = 21 ; 6.1
#Tri = 11 ; 3.2
#Frost-Frost = 343 - (101+21+11) = 210 ; 61.2

# GOVintersect
#So counts are GOV-frsot: 112-11 = 101; 0.81
#gov-mat: 258-11 = 247 ; 1.97
#Tri = 11 ; 0.088
#Frost-Frost = 12517 - (101+247+11) = 12158 ; 97.1


#Get count of number of unique clusters for whole network
diel_enviro_net_stats_overallclust=diel_enviro_overlap_simp%>%
  subset(VC.Status=="Clustered")%>%
  subset(!duplicated(VC.Subcluster))
#14,158 overall unique subclusters

#Get count of number of unique clusters for vmat network
diel_enviro_net_stats_overallclust=diel_enviro_overlap_simp%>%
  subset(VC.Status=="Clustered")%>%
  subset(Source=="vMAT")%>%
  subset(!duplicated(VC.Subcluster))
#4,177 original subclusters with 1,689 unique subclusters

#Create and excel sheet that summarizes all of this and do massive barplot
library(gdata)
sum_bar_df=read.csv('C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\env_net_summary_df.csv')

env_bar=ggplot(data=sum_bar_df,aes(fill=color,
                                   y=value,
                                   x=position)) +
  geom_col(colour="black",
           width=0.75,
           position=position_dodge(0.75)) +
  theme_classic()+
  xlab("VC Status")+
  ylab("Percent Belonging (%)")+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values=c("black","black","black","black","black"))
env_bar
ggsave("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\env_net_stats.svg",
       device="svg",
       width=8.5,
       height=3,
       units="in",
       dpi=1500)



#Examine Relationship between viral length and clustering status
#Read in viral lengths
viral_lengths=read.table('C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\curated.tsv',header=T)
viral_lengths_sub=viral_lengths%>%
  select(c(Genome,contig_length))
viral_clusts_simp=diel_enviro_overlap_simp%>%
  select(c(Genome,VC.Status))
length_clust_merge=merge(viral_clusts_simp,viral_lengths_sub,by="Genome")

ggplot(data=length_clust_merge,aes(x=VC.Status,y=log10(contig_length)))+
  geom_boxplot()+
  theme_classic()+
  xlab("VC Status")+
  ylab("Contig Length (log[bp])")+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=14,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))
ggsave("C:\\Users\\ecissell\\Documents\\Diel_Viruses_2022\\length_clust.png",
       device="png",
       width=12,
       height=8,
       units="in",
       dpi=600)

library("DHARMa")
lengthglm1=glm(log10(contig_length)~VC.Status,data=length_clust_merge)
summary(lengthglm1)
check_gamma_model1=simulateResiduals(fittedModel = lengthglm1, n=500 )
plot(check_gamma_model1)
#Fails miserably after trying families and normalization, switch to pairwise wilcox
pairwise.wilcox.test(x=length_clust_merge$contig_length,g=length_clust_merge$VC.Status)


#Obtain number of vMAT - vMAT clusters that span spatially distinct mats
#All vMAT clusters 1689
vmat_cluster_analysis=diel_enviro_overlap_simp%>%
  subset(Source=='vMAT') %>%
  subset(VC.Status=="Clustered")
  separate(Genome, into=c("MatID","SeqID",sep = 2)) %>%
  group_by(VC) %>%
  summarise(count = n_distinct(MatID)) %>%
  mutate(count = as.factor(count)) %>%
  ungroup()
count(vmat_cluster_analysis,count)
# 1       361
# 2       971
# 3       260
# 4        97
361/1689 #0.214
971/1689 #0.575
260/1689 #0.154
97/1689 #0.057

#vMAT - vMAT only 1410
vmat_cluster_analysis=diel_enviro_overlap_simp%>%
  subset(Source=='vMAT') %>%
  subset(VC.Status=="Clustered") %>%
  filter(!VC %in% mat_sea_clustersect$VC.x) %>%
  filter(!VC %in% mat_perm_clustersect$VC.x) %>%
  filter(!VC %in% tri_clustersect$VC.x) %>%
  mutate(VC = as.factor(VC))
  separate(Genome, into=c("MatID","SeqID",sep = 2)) %>%
  group_by(VC) %>%
  summarise(count = n_distinct(MatID)) %>%
  mutate(count = as.factor(count)) %>%
  ungroup()
count(vmat_cluster_analysis,count)
# 1       158
# 2       925
# 3       239
# 4        88
158/1410 #0.112
925/1410 #0.65.6
239/1410 #0.17.0
88/1410 #0.0624
158+925+239+88 #1410 sanity check

##SECTION 3 - VIRAL ABUNDANCES BY FAMILY
#Now read in PhaGCN results for taxonomic annotation
diel_phagcn=read.csv("/Users/cissell/Desktop/VIRSORTER2_Diel/PhagCN_result/final_prediction.csv")
#merge
viral_clust_curate_tax=merge(viral_clust_curate_clean,diel_phagcn,by="contig_id",all.x=T)

#merge in vcon taxonomy (from other script)
vcon_tax_6=dplyr::rename(vcon_tax_5,contig_id=vertex)
viral_clust_curate_tax=merge(viral_clust_curate_tax,vcon_tax_6,by="contig_id",all.x=T)

#merge in mat cluster taxonomy (from other script)
mat_clust_tax2=dplyr::rename(mat_clust_tax,contig_id=vertex)
viral_clust_curate_tax=merge(viral_clust_curate_tax,mat_clust_tax2,by="contig_id",all.x=T)

#coalesce these predictions
viral_clust_curate_tax=viral_clust_curate_tax%>%
  mutate(mergetax=coalesce(prediction,VCON_TAX))
viral_clust_curate_tax=viral_clust_curate_tax%>%
  mutate(mergetax=coalesce(mergetax,Mat_clust_tax_predict))
#Make NA taxonomy "Unknown"
viral_clust_curate_tax$mergetax=as.character(viral_clust_curate_tax$mergetax)
viral_clust_curate_tax$mergetax=viral_clust_curate_tax$mergetax %>%
  replace_na('Unknown')
viral_clust_curate_tax$mergetax=as.factor(viral_clust_curate_tax$mergetax)
#How many families
str(viral_clust_curate_tax) #10 families

#Now let's plot lengths by family
diel_viral_length=ggplot(data=viral_clust_curate_tax,aes(x=mergetax,
                                                         y=log10(contig_length),
                                                         fill=mergetax)) +
  geom_violin(size=0,
              alpha=0.8) +
  geom_point(aes(fill=mergetax,color=mergetax),position=position_jitterdodge(),pch=21,
             alpha=0.4,
             size=1.5)+
  coord_flip()+
  ylab("Viral sequence length (bp)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  guides(fill="none")+
  guides(color="none")+
  scale_x_discrete(limits = rev(levels(viral_clust_curate_tax$mergetax)))+
  scale_fill_manual(values=c("#0c0304",
                             "#a6b1ec",
                             "#241729",
                             "#aee3d2",
                             "#357aa2",
                             "#413f81",
                             "#3db4ad",
                             "#382954",
                             "#fa82a7",
                             "#3497a8",
                             "#8b9e90"))+
  scale_color_manual(values=c("#0c0304",
                              "#a6b1ec",
                              "#241729",
                              "#aee3d2",
                              "#357aa2",
                              "#413f81",
                              "#3db4ad",
                              "#382954",
                              "#fa82a7",
                              "#3497a8",
                              "#8b9e90"))
diel_viral_length
ggsave("/Users/cissell/Desktop/viral_length.svg",diel_viral_length,width=6.5,height=4)

#Obtain mean lengths for the new table in revised manuscript
table1lengths = viral_clust_curate_tax %>%
  group_by(mergetax) %>%
  summarise(mean = mean(contig_length))
#Now let's do stacked bar plot of proportions by family from counts and from abundance
viral_clust_curate_tax_props=viral_clust_curate_tax%>%
  group_by(mergetax)%>%
  dplyr::summarise(n=n())%>%
  mutate(freq=100*(n/sum(n))) %>%
  ungroup()%>%
  mutate(set="vMAT")
diel_viral_props=ggplot(data=viral_clust_curate_tax_props,aes(x=set,y=freq,fill=mergetax))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance (Raw Count %)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line.y = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"),
        axis.line.x= element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  coord_cartesian(ylim=c(84,100))+
  scale_fill_manual(values=c("#0c0304",
                             "#a6b1ec",
                             "#241729",
                             "#aee3d2",
                             "#357aa2",
                             "#413f81",
                             "#3db4ad",
                             "#382954",
                             "#fa82a7",
                             "#3497a8",
                             "#8b9e90"))
diel_viral_props
ggsave("/Users/cissell/Desktop/viral_props.svg",diel_viral_props,width=4.5,height=5)

#Make a dataframe with just contig name and tax
viral_tax_simpcol=viral_clust_curate_tax%>%
  select(contig_id,mergetax)

#Now let's do stacked bar plot of proportions removing unknown
viral_clust_curate_tax_props_unkrem=viral_clust_curate_tax%>%
  subset(!mergetax=="Unknown")%>%
  group_by(mergetax)%>%
  dplyr::summarise(n=n())%>%
  mutate(freq=100*(n/sum(n))) %>%
  ungroup()%>%
  mutate(set="vMAT")
View(viral_clust_curate_tax_props_unkrem)
diel_viral_props_unk_rem=ggplot(data=viral_clust_curate_tax_props_unkrem,aes(x=set,y=freq,fill=mergetax))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance (Raw Count %)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=13,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line.y = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"),
        axis.line.x= element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_manual(values=c("#0c0304",
                             "#a6b1ec",
                             "#241729",
                             "#aee3d2",
                             "#357aa2",
                             "#413f81",
                             "#3db4ad",
                             "#382954",
                             "#fa82a7",
                             "#3497a8"))
diel_viral_props_unk_rem
ggsave("/Users/cissell/Desktop/viral_props_nounk.svg",diel_viral_props_unk_rem,width=4.5,height=5)

#Now stacked bar plot of RNA abundance to TPM (this will reqire some calculations, and rely on df established in other paired scripts)

#First read in the RNA viral counts that werent inncluded in the host phage pairings
rviruscount1=read.xls("/Users/cissell/Desktop/VIRSORTER2_Diel/RNA_PHAGE_MAPS/R1_phage_counts.xlsx")
rviruslength1=read.xls("/Users/cissell/Desktop/VIRSORTER2_Diel/RNA_PHAGE_MAPS/R1_phage_length.xlsx")
rviruscm1=merge(rviruscount1,rviruslength1,by="contig_id")
rviruscount2=read.xls("/Users/cissell/Desktop/VIRSORTER2_Diel/RNA_PHAGE_MAPS/R2_phage_counts.xlsx")
rviruslength2=read.xls("/Users/cissell/Desktop/VIRSORTER2_Diel/RNA_PHAGE_MAPS/R2_phage_length.xlsx")
rviruscm2=merge(rviruscount2,rviruslength2,by="contig_id")
rviruscount3=read.xls("/Users/cissell/Desktop/VIRSORTER2_Diel/RNA_PHAGE_MAPS/R3_phage_counts.xlsx")
rviruslength3=read.xls("/Users/cissell/Desktop/VIRSORTER2_Diel/RNA_PHAGE_MAPS/R3_phage_length.xlsx")
rviruscm3=merge(rviruscount3,rviruslength3,by="contig_id")
rviruscount4=read.xls("/Users/cissell/Desktop/VIRSORTER2_Diel/RNA_PHAGE_MAPS/R4_phage_counts.xlsx")
rviruslength4=read.xls("/Users/cissell/Desktop/VIRSORTER2_Diel/RNA_PHAGE_MAPS/R4_phage_length.xlsx")
rviruscm4=merge(rviruscount4,rviruslength4,by="contig_id")


#Now df it up
rna_phage_sum_bar1=rna_phagediel1%>%
  select(-c(Geneid,Start,End))%>%
  mutate(PLength=Length)%>%
  select(-Length)%>%
  rbind(rviruscm1)%>%
  group_by(contig_id)%>%
  dplyr::summarise(across(1:6,sum))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)%>%
  merge(viral_tax_simpcol,by="contig_id")

rna_phage_sum_bar2=rna_phagediel2%>%
  select(-c(Geneid,Start,End))%>%
  mutate(PLength=Length)%>%
  select(-Length)%>%
  rbind(rviruscm2)%>%
  group_by(contig_id)%>%
  dplyr::summarise(across(1:6,sum))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)%>%
  merge(viral_tax_simpcol,by="contig_id")

rna_phage_sum_bar3=rna_phagediel3%>%
  select(-c(Geneid,Start,End))%>%
  mutate(PLength=Length)%>%
  select(-Length)%>%
  rbind(rviruscm3)%>%
  group_by(contig_id)%>%
  dplyr::summarise(across(1:6,sum))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)%>%
  merge(viral_tax_simpcol,by="contig_id")

rna_phage_sum_bar4=rna_phagediel4%>%
  select(-c(Geneid,Start,End))%>%
  mutate(PLength=Length)%>%
  select(-Length)%>%
  rbind(rviruscm4)%>%
  group_by(contig_id)%>%
  dplyr::summarise(across(1:6,sum))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)%>%
  merge(viral_tax_simpcol,by="contig_id")

#Now calculate TPM in RNA for each distinct mat, then we will bring back together
#Make new df
rpk1bar=rna_phage_sum_bar1
str(rpk1bar)
#Put length in KB
rpk1bar$PLength=rpk1bar$PLength/1000
#Create RPK for Mat 1
rpk1bar$P1_1=rpk1bar$P1_1/rpk1bar$PLength
rpk1bar$P1_2=rpk1bar$P1_2/rpk1bar$PLength
rpk1bar$P1_3=rpk1bar$P1_3/rpk1bar$PLength
rpk1bar$P1_4=rpk1bar$P1_4/rpk1bar$PLength
rpk1bar$P1_5=rpk1bar$P1_5/rpk1bar$PLength
#create scaling factors
#Mat 1
sum11bar=(sum(rpk1bar$P1_1)/1000000)
sum12bar=(sum(rpk1bar$P1_2)/1000000)
sum13bar=(sum(rpk1bar$P1_3)/1000000)
sum14bar=(sum(rpk1bar$P1_4)/1000000)
sum15bar=(sum(rpk1bar$P1_5)/1000000)
#Create TPM for Mat 1
rpk1bar$P1_1=rpk1bar$P1_1/sum11bar
rpk1bar$P1_2=rpk1bar$P1_2/sum12bar
rpk1bar$P1_3=rpk1bar$P1_3/sum13bar
rpk1bar$P1_4=rpk1bar$P1_4/sum14bar
rpk1bar$P1_5=rpk1bar$P1_5/sum15bar

sum(rpk1bar$P1_1)

#Bring it back to a sensical name now
D1_count_tpmbar=rpk1bar

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Summarize within family
D1_family_countsbar=D1_count_tpmbar%>%
  dplyr::select(-c(contig_id,PLength))%>%
  group_by(mergetax)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D1_family_countsbar$P1_1)
sum(D1_family_countsbar$P1_2)
sum(D1_family_countsbar$P1_3)
sum(D1_family_countsbar$P1_4)
sum(D1_family_countsbar$P1_5)
#Nice
#Pivot longer for later merging and plotting
D1_count_vmr_longbar=D1_family_countsbar%>%
  dplyr::rename(P_1=P1_1,
                P_2=P1_2,
                P_3=P1_3,
                P_4=P1_4,
                P_5=P1_5)%>%
  pivot_longer(
    cols=-c(mergetax),
    names_to=c(".value","Time"),
    names_sep="_")%>%
  mutate(Mat="D1")
D1_count_vmr_longbar=as.data.frame(D1_count_vmr_longbar)
D1_count_vmr_longbar$Mat=as.factor(D1_count_vmr_longbar$Mat)
D1_count_vmr_longbar$Time=as.factor(D1_count_vmr_longbar$Time)
str(D1_count_vmr_longbar)

#D2
#Make new df
rpk2bar=rna_phage_sum_bar2
str(rpk2bar)
#Put length in KB
rpk2bar$PLength=rpk2bar$PLength/1000
#Create RPK for Mat 1
rpk2bar$P2_1=rpk2bar$P2_1/rpk2bar$PLength
rpk2bar$P2_2=rpk2bar$P2_2/rpk2bar$PLength
rpk2bar$P2_3=rpk2bar$P2_3/rpk2bar$PLength
rpk2bar$P2_4=rpk2bar$P2_4/rpk2bar$PLength
rpk2bar$P2_5=rpk2bar$P2_5/rpk2bar$PLength
#create scaling factors
#Mat 2
sum21bar=(sum(rpk2bar$P2_1)/1000000)
sum22bar=(sum(rpk2bar$P2_2)/1000000)
sum23bar=(sum(rpk2bar$P2_3)/1000000)
sum24bar=(sum(rpk2bar$P2_4)/1000000)
sum25bar=(sum(rpk2bar$P2_5)/1000000)
#Create TPM for Mat 2
rpk2bar$P2_1=rpk2bar$P2_1/sum21bar
rpk2bar$P2_2=rpk2bar$P2_2/sum22bar
rpk2bar$P2_3=rpk2bar$P2_3/sum23bar
rpk2bar$P2_4=rpk2bar$P2_4/sum24bar
rpk2bar$P2_5=rpk2bar$P2_5/sum25bar

sum(rpk2bar$P2_1)

#Bring it back to a sensical name now
D2_count_tpmbar=rpk2bar

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Summarize within family
D2_family_countsbar=D2_count_tpmbar%>%
  dplyr::select(-c(contig_id,PLength))%>%
  group_by(mergetax)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D2_family_countsbar$P2_1)
sum(D2_family_countsbar$P2_2)
sum(D2_family_countsbar$P2_3)
sum(D2_family_countsbar$P2_4)
sum(D2_family_countsbar$P2_5)
#Nice
#Pivot longer for later merging and plotting
D2_count_vmr_longbar=D2_family_countsbar%>%
  dplyr::rename(P_1=P2_1,
                P_2=P2_2,
                P_3=P2_3,
                P_4=P2_4,
                P_5=P2_5)%>%
  pivot_longer(
    cols=-c(mergetax),
    names_to=c(".value","Time"),
    names_sep="_")%>%
  mutate(Mat="D2")
D2_count_vmr_longbar=as.data.frame(D2_count_vmr_longbar)
D2_count_vmr_longbar$Mat=as.factor(D2_count_vmr_longbar$Mat)
D2_count_vmr_longbar$Time=as.factor(D2_count_vmr_longbar$Time)
str(D2_count_vmr_longbar)

#D3
#Make new df
rpk3bar=rna_phage_sum_bar3
str(rpk3bar)
#Put length in KB
rpk3bar$PLength=rpk3bar$PLength/1000
#Create RPK for Mat 1
rpk3bar$P3_1=rpk3bar$P3_1/rpk3bar$PLength
rpk3bar$P3_2=rpk3bar$P3_2/rpk3bar$PLength
rpk3bar$P3_3=rpk3bar$P3_3/rpk3bar$PLength
rpk3bar$P3_4=rpk3bar$P3_4/rpk3bar$PLength
rpk3bar$P3_5=rpk3bar$P3_5/rpk3bar$PLength
#create scaling factors
#Mat 3
sum31bar=(sum(rpk3bar$P3_1)/1000000)
sum32bar=(sum(rpk3bar$P3_2)/1000000)
sum33bar=(sum(rpk3bar$P3_3)/1000000)
sum34bar=(sum(rpk3bar$P3_4)/1000000)
sum35bar=(sum(rpk3bar$P3_5)/1000000)
#Create TPM for Mat 3
rpk3bar$P3_1=rpk3bar$P3_1/sum31bar
rpk3bar$P3_2=rpk3bar$P3_2/sum32bar
rpk3bar$P3_3=rpk3bar$P3_3/sum33bar
rpk3bar$P3_4=rpk3bar$P3_4/sum34bar
rpk3bar$P3_5=rpk3bar$P3_5/sum35bar

sum(rpk3bar$P3_1)

#Bring it back to a sensical name now
D3_count_tpmbar=rpk3bar

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Summarize within family
D3_family_countsbar=D3_count_tpmbar%>%
  dplyr::select(-c(contig_id,PLength))%>%
  group_by(mergetax)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D3_family_countsbar$P3_1)
sum(D3_family_countsbar$P3_2)
sum(D3_family_countsbar$P3_3)
sum(D3_family_countsbar$P3_4)
sum(D3_family_countsbar$P3_5)
#Nice
#Pivot longer for later merging and plotting
D3_count_vmr_longbar=D3_family_countsbar%>%
  dplyr::rename(P_1=P3_1,
                P_2=P3_2,
                P_3=P3_3,
                P_4=P3_4,
                P_5=P3_5)%>%
  pivot_longer(
    cols=-c(mergetax),
    names_to=c(".value","Time"),
    names_sep="_")%>%
  mutate(Mat="D3")
D3_count_vmr_longbar=as.data.frame(D3_count_vmr_longbar)
D3_count_vmr_longbar$Mat=as.factor(D3_count_vmr_longbar$Mat)
D3_count_vmr_longbar$Time=as.factor(D3_count_vmr_longbar$Time)
str(D3_count_vmr_longbar)

#D4
#Make new df
rpk4bar=rna_phage_sum_bar4
str(rpk4bar)
#Put length in KB
rpk4bar$PLength=rpk4bar$PLength/1000
#Create RPK for Mat 1
rpk4bar$P4_1=rpk4bar$P4_1/rpk4bar$PLength
rpk4bar$P4_2=rpk4bar$P4_2/rpk4bar$PLength
rpk4bar$P4_3=rpk4bar$P4_3/rpk4bar$PLength
rpk4bar$P4_4=rpk4bar$P4_4/rpk4bar$PLength
rpk4bar$P4_5=rpk4bar$P4_5/rpk4bar$PLength
#create scaling factors
#Mat 4
sum41bar=(sum(rpk4bar$P4_1)/1000000)
sum42bar=(sum(rpk4bar$P4_2)/1000000)
sum43bar=(sum(rpk4bar$P4_3)/1000000)
sum44bar=(sum(rpk4bar$P4_4)/1000000)
sum45bar=(sum(rpk4bar$P4_5)/1000000)
#Create TPM for Mat 4
rpk4bar$P4_1=rpk4bar$P4_1/sum41bar
rpk4bar$P4_2=rpk4bar$P4_2/sum42bar
rpk4bar$P4_3=rpk4bar$P4_3/sum43bar
rpk4bar$P4_4=rpk4bar$P4_4/sum44bar
rpk4bar$P4_5=rpk4bar$P4_5/sum45bar

sum(rpk4bar$P4_1)

#Bring it back to a sensical name now
D4_count_tpmbar=rpk4bar

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Summarize within family
D4_family_countsbar=D4_count_tpmbar%>%
  dplyr::select(-c(contig_id,PLength))%>%
  group_by(mergetax)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D4_family_countsbar$P4_1)
sum(D4_family_countsbar$P4_2)
sum(D4_family_countsbar$P4_3)
sum(D4_family_countsbar$P4_4)
sum(D4_family_countsbar$P4_5)
#Nice

#Pivot longer for later merging and plotting
D4_count_vmr_longbar=D4_family_countsbar%>%
  dplyr::rename(P_1=P4_1,
                P_2=P4_2,
                P_3=P4_3,
                P_4=P4_4,
                P_5=P4_5)%>%
  pivot_longer(
    cols=-c(mergetax),
    names_to=c(".value","Time"),
    names_sep="_")%>%
  mutate(Mat="D4")
D4_count_vmr_longbar=as.data.frame(D4_count_vmr_longbar)
D4_count_vmr_longbar$Mat=as.factor(D4_count_vmr_longbar$Mat)
D4_count_vmr_longbar$Time=as.factor(D4_count_vmr_longbar$Time)
str(D4_count_vmr_longbar)

#Merge all together now
Phage_host_count_mergebar=rbind(D1_count_vmr_longbar,D2_count_vmr_longbar,D3_count_vmr_longbar,D4_count_vmr_longbar)

#Add molecule column
Phage_host_count_mergebar=Phage_host_count_mergebar%>%
  mutate(Molecule="R")
#Now lets merge in the DNA to make one plot, then the other stacks will be little cutouts in that
D1_DNA_phage_bar=D1_Phage_counts%>%
  mutate(contig_id=Phage)%>%
  select(-Phage)%>%
  merge(viral_tax_simpcol,by="contig_id")

D2_DNA_phage_bar=D2_Phage_counts%>%
  mutate(contig_id=Phage)%>%
  select(-Phage)%>%
  merge(viral_tax_simpcol,by="contig_id")

D3_DNA_phage_bar=D3_Phage_counts%>%
  mutate(contig_id=Phage)%>%
  select(-Phage)%>%
  merge(viral_tax_simpcol,by="contig_id")

D4_DNA_phage_bar=D4_Phage_counts%>%
  mutate(contig_id=Phage)%>%
  select(-Phage)%>%
  merge(viral_tax_simpcol,by="contig_id")

#Now calculate CPM in DNA for each distinct mat, then we will bring back together
#Make new df
rpk1bard=D1_DNA_phage_bar
str(rpk1bard)
#Put length in KB
rpk1bard$PLength=rpk1bard$PLength/1000
#Create RPK for Mat 1
rpk1bard$P1_1=rpk1bard$P1_1/rpk1bard$PLength
rpk1bard$P1_3=rpk1bard$P1_3/rpk1bard$PLength
rpk1bard$P1_5=rpk1bard$P1_5/rpk1bard$PLength
#create scaling factors
#Mat 1
sum11bard=(sum(rpk1bard$P1_1)/1000000)
sum13bard=(sum(rpk1bard$P1_3)/1000000)
sum15bard=(sum(rpk1bard$P1_5)/1000000)
#Create TPM for Mat 1
rpk1bard$P1_1=rpk1bard$P1_1/sum11bard
rpk1bard$P1_3=rpk1bard$P1_3/sum13bard
rpk1bard$P1_5=rpk1bard$P1_5/sum15bard

sum(rpk1bard$P1_1)

#Bring it back to a sensical name now
D1_count_tpmbard=rpk1bard

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Summarize within family
D1_family_countsbard=D1_count_tpmbard%>%
  select(-c(contig_id,PLength))%>%
  group_by(mergetax)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D1_family_countsbard$P1_1)
sum(D1_family_countsbard$P1_3)
sum(D1_family_countsbard$P1_5)
#Nice
#Pivot longer for later merging and plotting
D1_count_vmr_longbard=D1_family_countsbard%>%
  dplyr::rename(P_1=P1_1,
                P_3=P1_3,
                P_5=P1_5)%>%
  pivot_longer(
    cols=-c(mergetax),
    names_to=c(".value","Time"),
    names_sep="_")%>%
  mutate(Mat="D1")
D1_count_vmr_longbard=as.data.frame(D1_count_vmr_longbard)
D1_count_vmr_longbard$Mat=as.factor(D1_count_vmr_longbard$Mat)
D1_count_vmr_longbard$Time=as.factor(D1_count_vmr_longbard$Time)
str(D1_count_vmr_longbard)

#D2
#Make new df
rpk2bard=D2_DNA_phage_bar
str(rpk2bard)
#Put length in KB
rpk2bard$PLength=rpk2bard$PLength/1000
#Create RPK for Mat 1
rpk2bard$P2_1=rpk2bard$P2_1/rpk2bard$PLength
rpk2bard$P2_3=rpk2bard$P2_3/rpk2bard$PLength
rpk2bard$P2_5=rpk2bard$P2_5/rpk2bard$PLength
#create scaling factors
#Mat 2
sum21bard=(sum(rpk2bard$P2_1)/1000000)
sum23bard=(sum(rpk2bard$P2_3)/1000000)
sum25bard=(sum(rpk2bard$P2_5)/1000000)
#Create TPM for Mat 2
rpk2bard$P2_1=rpk2bard$P2_1/sum21bard
rpk2bard$P2_3=rpk2bard$P2_3/sum23bard
rpk2bard$P2_5=rpk2bard$P2_5/sum25bard

sum(rpk2bard$P2_1)

#Bring it back to a sensical name now
D2_count_tpmbard=rpk2bard

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Summarize within family
D2_family_countsbard=D2_count_tpmbard%>%
  select(-c(contig_id,PLength))%>%
  group_by(mergetax)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D2_family_countsbard$P2_1)
sum(D2_family_countsbard$P2_3)
sum(D2_family_countsbard$P2_5)
#Nice
#Pivot longer for later merging and plotting
D2_count_vmr_longbard=D2_family_countsbard%>%
  dplyr::rename(P_1=P2_1,
                P_3=P2_3,
                P_5=P2_5)%>%
  pivot_longer(
    cols=-c(mergetax),
    names_to=c(".value","Time"),
    names_sep="_")%>%
  mutate(Mat="D2")
D2_count_vmr_longbard=as.data.frame(D2_count_vmr_longbard)
D2_count_vmr_longbard$Mat=as.factor(D2_count_vmr_longbard$Mat)
D2_count_vmr_longbard$Time=as.factor(D2_count_vmr_longbard$Time)
str(D2_count_vmr_longbard)

#D3
#Make new df
rpk3bard=D3_DNA_phage_bar
str(rpk3bard)
#Put length in KB
rpk3bard$PLength=rpk3bard$PLength/1000
#Create RPK for Mat 1
rpk3bard$P3_1=rpk3bard$P3_1/rpk3bard$PLength
rpk3bard$P3_3=rpk3bard$P3_3/rpk3bard$PLength
rpk3bard$P3_5=rpk3bard$P3_5/rpk3bard$PLength
#create scaling factors
#Mat 3
sum31bard=(sum(rpk3bard$P3_1)/1000000)
sum33bard=(sum(rpk3bard$P3_3)/1000000)
sum35bard=(sum(rpk3bard$P3_5)/1000000)
#Create TPM for Mat 3
rpk3bard$P3_1=rpk3bard$P3_1/sum31bard
rpk3bard$P3_3=rpk3bard$P3_3/sum33bard
rpk3bard$P3_5=rpk3bard$P3_5/sum35bard

sum(rpk3bard$P3_1)

#Bring it back to a sensical name now
D3_count_tpmbard=rpk3bard

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Summarize within family
D3_family_countsbard=D3_count_tpmbard%>%
  select(-c(contig_id,PLength))%>%
  group_by(mergetax)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D3_family_countsbard$P3_1)
sum(D3_family_countsbard$P3_3)
sum(D3_family_countsbard$P3_5)
#Nice
#Pivot longer for later merging and plotting
D3_count_vmr_longbard=D3_family_countsbard%>%
  dplyr::rename(P_1=P3_1,
                P_3=P3_3,
                P_5=P3_5)%>%
  pivot_longer(
    cols=-c(mergetax),
    names_to=c(".value","Time"),
    names_sep="_")%>%
  mutate(Mat="D3")
D3_count_vmr_longbard=as.data.frame(D3_count_vmr_longbard)
D3_count_vmr_longbard$Mat=as.factor(D3_count_vmr_longbard$Mat)
D3_count_vmr_longbard$Time=as.factor(D3_count_vmr_longbard$Time)
str(D3_count_vmr_longbard)

#D4
#Make new df
rpk4bard=D4_DNA_phage_bar
str(rpk4bard)
#Put length in KB
rpk4bard$PLength=rpk4bard$PLength/1000
#Create RPK for Mat 1
rpk4bard$P4_1=rpk4bard$P4_1/rpk4bard$PLength
rpk4bard$P4_3=rpk4bard$P4_3/rpk4bard$PLength
rpk4bard$P4_5=rpk4bard$P4_5/rpk4bard$PLength
#create scaling factors
#Mat 4
sum41bard=(sum(rpk4bard$P4_1)/1000000)
sum43bard=(sum(rpk4bard$P4_3)/1000000)
sum45bard=(sum(rpk4bard$P4_5)/1000000)
#Create TPM for Mat 4
rpk4bard$P4_1=rpk4bard$P4_1/sum41bard
rpk4bard$P4_3=rpk4bard$P4_3/sum43bard
rpk4bard$P4_5=rpk4bard$P4_5/sum45bard

sum(rpk4bard$P4_1)

#Bring it back to a sensical name now
D4_count_tpmbard=rpk4bard

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Summarize within family
D4_family_countsbard=D4_count_tpmbard%>%
  select(-c(contig_id,PLength))%>%
  group_by(mergetax)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D4_family_countsbard$P4_1)
sum(D4_family_countsbard$P4_3)
sum(D4_family_countsbard$P4_5)
#Nice

#Pivot longer for later merging and plotting
D4_count_vmr_longbard=D4_family_countsbard%>%
  dplyr::rename(P_1=P4_1,
                P_3=P4_3,
                P_5=P4_5)%>%
  pivot_longer(
    cols=-c(mergetax),
    names_to=c(".value","Time"),
    names_sep="_")%>%
  mutate(Mat="D4")
D4_count_vmr_longbard=as.data.frame(D4_count_vmr_longbard)
D4_count_vmr_longbard$Mat=as.factor(D4_count_vmr_longbard$Mat)
D4_count_vmr_longbard$Time=as.factor(D4_count_vmr_longbard$Time)
str(D4_count_vmr_longbard)

#Merge all DNA together now
Phage_host_count_mergebard=rbind(D1_count_vmr_longbard,D2_count_vmr_longbard,D3_count_vmr_longbard,D4_count_vmr_longbard)
Phage_host_count_mergebard=Phage_host_count_mergebard%>%
  mutate(Molecule="D")
Phage_D_R_bar=rbind(Phage_host_count_mergebard,Phage_host_count_mergebar)
Phage_D_R_bar$mergetax=as.factor(Phage_D_R_bar$mergetax)
str(Phage_D_R_bar)

#Plot
library(metR)
library(magrittr)
library(ggnewscale)
diel_viral_fam_tpm=ggplot(data=Phage_D_R_bar,aes(x=mergetax,
                                                 y=log10(P),
                                                 group=interaction(Molecule,mergetax))) +
  geom_boxplot(aes(fill=mergetax,
                   color=Molecule),
               size=0.5,
               alpha=0.8,
               outlier.shape=NA,
               position=position_dodge(-0.9))+
  coord_flip()+
  ylab("Viral Abundance (logCPM or logTPM)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  guides(fill="none")+
  scale_x_discrete(limits = rev(levels(Phage_D_R_bar$mergetax)))+
  scale_fill_manual(values=c("#0c0304",
                             "#a6b1ec",
                             "#241729",
                             "#aee3d2",
                             "#357aa2",
                             "#413f81",
                             "#3db4ad",
                             "#382954",
                             "#fa82a7",
                             "#3497a8",
                             "#8b9e90"))+
  scale_color_manual(values=c("#4a6274","#94acbf"))+
  new_scale_color()+
  geom_point(aes(color=Molecule),position=position_dodge(-0.9),
             pch=19,
             alpha=0.4,
             size=1)+
  scale_color_manual(values=c("#4a6274","#94acbf"))
diel_viral_fam_tpm
ggsave("/Users/cissell/Desktop/viral_abunds.svg",diel_viral_fam_tpm,width=7,height=3.55)

#Plot RNA by sampling time
diel_viral_fam_tpm_time=ggplot(data=subset(Phage_D_R_bar,Molecule=="R"),aes(x=mergetax,
                                                 y=log10(P),
                                                 group=interaction(mergetax,Time))) +
  geom_boxplot(aes(fill=Time),
               color="black",
               size=0.5,
               alpha=0.8,
               outlier.shape=NA,
               position=position_dodge(-0.9))+
  coord_flip()+
  ylab("Viral Abundance (logTPM)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_x_discrete(limits = rev(levels(Phage_D_R_bar$mergetax)))+
  scale_fill_manual(values=c("#90a295",
                             "#9c8cdb",
                             "#455765",
                             "#add8e6",
                             "#cd9fb2"))+
#  scale_color_manual(values=c("#4a6274","#94acbf"))+
  new_scale_color()+
  geom_point(aes(color=Molecule),position=position_dodge(-0.9),
             pch=19,
             alpha=0.4,
             size=1)+
  scale_color_manual(values=c("#4a6274","#94acbf"))
diel_viral_fam_tpm_time
ggsave("/Users/cissell/Desktop/viral_abunds.svg",diel_viral_fam_tpm,width=7,height=3.55)

#Plot DNA by sampling time
diel_viral_fam_tpm_time_d=ggplot(data=subset(Phage_D_R_bar,Molecule=="D"),aes(x=mergetax,
                                                                            y=log10(P),
                                                                            group=interaction(mergetax,Time))) +
  geom_boxplot(aes(fill=Time),
               color="black",
               size=0.5,
               alpha=0.8,
               outlier.shape=NA,
               position=position_dodge(-0.9))+
  coord_flip()+
  ylab("Viral Abundance (logCPM)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_x_discrete(limits = rev(levels(Phage_D_R_bar$mergetax)))+
  scale_fill_manual(values=c("#90a295",
 #                            "#9c8cdb",
                             "#455765",
 #                            "#add8e6",
                             "#cd9fb2"))+
  #  scale_color_manual(values=c("#4a6274","#94acbf"))+
  new_scale_color()+
  geom_point(aes(color=Molecule),position=position_dodge(-0.9),
             pch=19,
             alpha=0.4,
             size=1)+
  scale_color_manual(values=c("#4a6274","#94acbf"))
diel_viral_fam_tpm_time_d
ggsave("/Users/cissell/Desktop/dielviral_famabunddnatime.svg",diel_viral_fam_tpm_time_d,width=7,height=7)

#Plot RNA by sampling time
diel_viral_fam_tpm_time_r=ggplot(data=subset(Phage_D_R_bar,Molecule=="R"),aes(x=mergetax,
                                                                              y=log10(P),
                                                                              group=interaction(mergetax,Time))) +
  geom_boxplot(aes(fill=Time),
               color="black",
               size=0.5,
               alpha=0.8,
               outlier.shape=NA,
               position=position_dodge(-0.9))+
  coord_flip()+
  ylab("Viral Activity (logTPM)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_x_discrete(limits = rev(levels(Phage_D_R_bar$mergetax)))+
  scale_fill_manual(values=c("#90a295",
                             "#9c8cdb",
                             "#455765",
                             "#add8e6",
                             "#cd9fb2"))+
  #  scale_color_manual(values=c("#4a6274","#94acbf"))+
  new_scale_color()+
  geom_point(aes(color=Molecule),position=position_dodge(-0.9),
             pch=19,
             alpha=0.4,
             size=1)+
  scale_color_manual(values=c("#4a6274","#94acbf"))
diel_viral_fam_tpm_time_r
ggsave("/Users/cissell/Desktop/dielviral_famabundrnatime.svg",diel_viral_fam_tpm_time_r,width=7,height=7)

#Plot DNA and RNA by sampled Mat identity
diel_viral_fam_cpm_mat=ggplot(data=subset(Phage_D_R_bar,Molecule=="D"),aes(x=mergetax,
                                                 y=log10(P),
                                                 group=interaction(mergetax,Mat))) +
  geom_boxplot(aes(fill=Mat),
               color="black",
               size=0.5,
               alpha=0.8,
               outlier.shape=NA,
               position=position_dodge(-0.9))+
  coord_flip()+
  ylab("Viral Abundance (logCPM)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_x_discrete(limits = rev(levels(Phage_D_R_bar$mergetax)))+
  scale_fill_manual(values=c("#90a295",
                             "#9c8cdb",
                             "#455765",
                             "#add8e6",
                             "#cd9fb2"))+
  scale_color_manual(values=c("#4a6274","#94acbf"))+
  new_scale_color()+
  geom_point(aes(color=Molecule),position=position_dodge(-0.9),
             pch=19,
             alpha=0.4,
             size=1)+
  scale_color_manual(values=c("#4a6274","#94acbf"))
diel_viral_fam_cpm_mat
ggsave("/Users/cissell/Desktop/dielviral_famabunddnamat.svg",diel_viral_fam_cpm_mat,width=7,height=7)


diel_viral_fam_tpm_mat=ggplot(data=subset(Phage_D_R_bar,Molecule=="R"),aes(x=mergetax,
                                                                           y=log10(P),
                                                                           group=interaction(mergetax,Mat))) +
  geom_boxplot(aes(fill=Mat),
               color="black",
               size=0.5,
               alpha=0.8,
               outlier.shape=NA,
               position=position_dodge(-0.9))+
  coord_flip()+
  ylab("Viral Abundance (logTPM)")+
  theme_classic()+
  xlab("")+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  scale_x_discrete(limits = rev(levels(Phage_D_R_bar$mergetax)))+
  scale_fill_manual(values=c("#90a295",
                             "#9c8cdb",
                             "#455765",
                             "#add8e6",
                             "#cd9fb2"))+
  scale_color_manual(values=c("#4a6274","#94acbf"))+
  new_scale_color()+
  geom_point(aes(color=Molecule),position=position_dodge(-0.9),
             pch=19,
             alpha=0.4,
             size=1)+
  scale_color_manual(values=c("#4a6274","#94acbf"))
diel_viral_fam_tpm_mat
ggsave("/Users/cissell/Desktop/dielviral_famabundrnamat.svg",diel_viral_fam_tpm_mat,width=7,height=7)

##SECTION 4
###########
#Predicting phage-host linkages by combining multiple in-silico approaches
###########
library(dplyr)
library(gdata)
library(tidyverse)
library(conflicted)
library(car)
library(lme4)
conflict_prefer("filter", "dplyr")
conflict_prefer("print", "EnvStats")
conflict_prefer("mutate", "dplyr")
conflict_prefer("rename", "dplyr")
######
#Read in files
######
#Read in files for D1
D1_crispr=read.table("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/CRISPR/DR1_blast_spacer.txt",header=T)
D1_nt=read.table("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/BLAST/D1_blast_nt_id.txt",header=T)
D1_dstar=read.csv("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/PHIST/D1_predictions.csv")

#Read in files for D2
D2_crispr=read.table("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/CRISPR/DR2_blast_spacer.txt",header=T)
D2_nt=read.table("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/BLAST/D2_blast_nt_id.txt",header=T)
D2_dstar=read.csv("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/PHIST/D2_predictions.csv")

#Read in files for D3
D3_crispr=read.table("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/CRISPR/DR3_blast_spacer.txt",header=T)
D3_nt=read.table("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/BLAST/D3_blast_nt_id.txt",header=T)
D3_dstar=read.csv("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/PHIST/D3_predictions.csv")

#Read in files for D4
D4_crispr=read.table("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/CRISPR/DR4_blast_spacer.txt",header=T)
D4_nt=read.table("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/BLAST/D4_blast_nt_id.txt",header=T)
D4_dstar=read.csv("/Users/cissell/Desktop/VIRSORTER2_diel/Host_Phage/PHIST/D4_predictions.csv")
#######
######D1
#Filter spacer by mismatch ≤1 and e-value ≤1, then best bitscore
D1_crispr_filter=D1_crispr %>%
  filter(mismatch<=1 & e_value <=1) %>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  dplyr::mutate(Source="crispr")

#Order by phage ID
D1_crispr_filter=D1_crispr_filter[order(D1_crispr_filter$Phage),]

#Filter nt by ≥50 bitscore and ≤0.001 e value
D1_nt_filter=D1_nt%>%
  filter(bitscore>=50 & e_value <=0.001)

#Extract only highest bitscore for each unique phage ID
D1_nt_high=D1_nt_filter%>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="nt")

#Order by phage ID
D1_nt_high=D1_nt_high[order(D1_nt_high$Phage),]

#Remove blanks from Kmer dataframe
D1_dstar=na.omit(D1_dstar)
D2_dstar=na.omit(D2_dstar)
D3_dstar=na.omit(D3_dstar)
D4_dstar=na.omit(D4_dstar)
#Find median k-mer matches number to filter by
median(D1_dstar$common.kmers)
median(D2_dstar$common.kmers)
median(D3_dstar$common.kmers)
median(D4_dstar$common.kmers)
mean(c(41,25,23,53))#36

#Filter d2star by common.kmer>36
D1_dstar_filter=D1_dstar%>%
  filter(common.kmers>36) %>%
  mutate(Source="kmer")

#Now create smaller dataframes from each for merging
D1_crisp_merger=D1_crispr_filter[,c(2,3,14)]
D1_nt_merger=D1_nt_high[,c(2,3,14)]
D1_dstar_merger=D1_dstar_filter[,c(1,2,6)]

#Now merge
D1_ph_merge1=merge(D1_crisp_merger,D1_nt_merger, by="Phage",all=T)
D1_ph_merge_2=merge(D1_ph_merge1,D1_dstar_merger,by="Phage",all=T)

#Now create new column with single best prediction
D1_ph=D1_ph_merge_2
D1_ph=rename(D1_ph,MAG_cr=MAG.x)
D1_ph=rename(D1_ph,MAG_nt=MAG.y)
D1_ph$MAG_cr=as.character(D1_ph$MAG_cr)
D1_ph$MAG_nt=as.character(D1_ph$MAG_nt)
D1_ph$MAG=as.character(D1_ph$MAG)
D1_ph=D1_ph%>%
  mutate(best=ifelse(is.na(MAG_cr) & is.na(MAG_nt),MAG[],
                     ifelse(is.na(MAG_cr),MAG_nt[],MAG_cr[])))

#Get length of unique phages to make sure all persist to end
length(unique(D1_ph$Phage))#1212

##Eliminate redundancy based on conditions
#Create dataframe without duplicates for later combining
D1_ph_nd=D1_ph[!c(duplicated(D1_ph$Phage) | duplicated(D1_ph$Phage, fromLast=T)),]
#Check
duplicated(D1_ph_nd$Phage)
length(unique(D1_ph_nd$Phage))#1051
#Create dataframe with all duplicates
D1_ph_d=D1_ph[duplicated(D1_ph$Phage) | duplicated(D1_ph$Phage, fromLast=T),]
#Check
length(unique(D1_ph_d$Phage))#161, adds to 1212 so still good
#Focus dataframe on just meaningful columns
D1_ph_d=D1_ph_d[,-c(3,5:8)]
#Create dummy column and fill with NA
D1_ph_d$Con=NA
#Fill dummy column with Yes if crispr and nt match
D1_ph_d$Con[D1_ph_d$MAG_cr==D1_ph_d$MAG_nt] = "Yes"
#Fill dummy column with no if crisp and nt do not match
D1_ph_d$Con[D1_ph_d$MAG_cr!=D1_ph_d$MAG_nt] = "No"

#Extract Concencus dataframe for later combining
Concensus1=D1_ph_d%>%
  filter(Con=="Yes")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best)

#2 of 161 recovered
#These 2 cases had a duplicate 'No' that we need to remove from the frame for downstream
D1_ph_d_non=D1_ph_d%>%
  group_by(Phage)%>%
  filter(!Phage%in%Concensus1$Phage)


#Pipeline to remove those that had concencus and then if not bring down to a single column with redundant hits
D1_ph_dup=D1_ph_d_non%>%
  group_by(Phage)%>%
  filter(is.na(Con))%>%
  mutate(Con=ifelse(is.na(MAG_nt),MAG_cr[],MAG_nt[]))%>%
  select(-MAG_cr,-MAG_nt) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(D1_ph_dup$Phage)) #130 of 161 recovered, that means No should have 29

#Pipeline to extract duplicates supported with multiple identical assignments and make nonredundant for later combining and to create vector for subsetting other dataframe
D1_ph_d2=D1_ph_dup%>%
  group_by(Phage)%>%
  filter(c(duplicated(Con)) | duplicated(Con, fromLast=T))%>%
  filter(duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)%>%
  filter(!duplicated(Phage))
length(unique(D1_ph_d2$Phage)) #44

#Create vector for subsetting
r1=as.vector(D1_ph_d2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
D1_ph_nd_2=D1_ph_dup%>%
  group_by(Phage)%>%
  filter(!Phage%in%r1)%>%
  filter(!duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)
length(unique(D1_ph_nd_2$Phage)) #86 recovered, D1_ph_nd_2+D1_ph_d2=130 recovered

#Pipeline to extract thos assigned 'No', use CRISPR as host, and remove redundancy
#First extract those supported multiple times by Crispr
Noncensus1=D1_ph_d_non%>%
  filter(Con=="No")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(Noncensus1$Phage))

Noncensus1_2=Noncensus1%>%
  group_by(Phage)%>%
  filter(c(duplicated(best)) | duplicated(best, fromLast=T)) %>%
  filter(duplicated(Phage))%>%
  mutate(best=best)%>%
  select(-vid)%>%
  filter(!duplicated(Phage))
length(unique(Noncensus1_2$Phage)) #16

#Create another vector for subsetting
rnoncen1=as.vector(Noncensus1_2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
Noncensus1_3=Noncensus1%>%
  group_by(Phage)%>%
  filter(!Phage%in%rnoncen1)%>%
  filter(!duplicated(Phage))%>%
  select(-vid)
length(unique(Noncensus1_3$Phage)) #13

#Make original dataframe fit the columnns
D1_ph_nd_o=D1_ph_nd%>%
  select(Phage,best)
#Create final dataframe
D1_phage_host_final=bind_rows(D1_ph_nd_o,
                              D1_ph_d2,
                              D1_ph_nd_2,
                              Concensus1,
                              Noncensus1_2,
                              Noncensus1_3)
#Make hosts factor
D1_phage_host_final$best=as.factor(D1_phage_host_final$best)
#Rename best to Host
D1_phage_host_final=rename(D1_phage_host_final,Host=best)

#Add column with matID
D1_phage_host_final=D1_phage_host_final%>%
  mutate(Mat="D1")

#Ensure there are no duplicates
length(unique(D1_phage_host_final$Phage)) #1212
duplicated(D1_phage_host_final$Phage)




#######
#######
#######
#######
######D2
#Filter spacer by mismatch ≤1 and e-value ≤1, then best bitscore
D2_crispr_filter=D2_crispr %>%
  filter(mismatch<=1 & e_value <=1) %>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="crispr")

#Order by phage ID
D2_crispr_filter=D2_crispr_filter[order(D2_crispr_filter$Phage),]

#Filter nt by ≥50 bitscore and ≤0.001 e value
D2_nt_filter=D2_nt%>%
  filter(bitscore>=50 & e_value <=0.001)

#Extract only highest bitscore for each unique phage ID
D2_nt_high=D2_nt_filter%>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="nt")

#Order by phage ID
D2_nt_high=D2_nt_high[order(D2_nt_high$Phage),]

#Filter d2star by common.kmer>36
D2_dstar_filter=D2_dstar%>%
  filter(common.kmers>36) %>%
  mutate(Source="kmer")

#Now create smaller dataframes from each for merging
D2_crisp_merger=D2_crispr_filter[,c(2,3,14)]
D2_nt_merger=D2_nt_high[,c(2,3,14)]
D2_dstar_merger=D2_dstar_filter[,c(1,2,6)]

#Now merge
D2_ph_merge1=merge(D2_crisp_merger,D2_nt_merger, by="Phage",all=T)
D2_ph_merge_2=merge(D2_ph_merge1,D2_dstar_merger,by="Phage",all=T)

#Now create new column with single best prediction
D2_ph=D2_ph_merge_2
D2_ph=rename(D2_ph,MAG_cr=MAG.x)
D2_ph=rename(D2_ph,MAG_nt=MAG.y)
D2_ph$MAG_cr=as.character(D2_ph$MAG_cr)
D2_ph$MAG_nt=as.character(D2_ph$MAG_nt)
D2_ph$MAG=as.character(D2_ph$MAG)
D2_ph=D2_ph%>%
  mutate(best=ifelse(is.na(MAG_cr) & is.na(MAG_nt),MAG[],
                     ifelse(is.na(MAG_cr),MAG_nt[],MAG_cr[])))

#Get length of unique phages to make sure all persist to end
length(unique(D2_ph$Phage))#1398

##Eliminate redundancy based on conditions
#Create dataframe without duplicates for later combining
D2_ph_nd=D2_ph[!c(duplicated(D2_ph$Phage) | duplicated(D2_ph$Phage, fromLast=T)),]
#Check
length(unique(D2_ph_nd$Phage))#1285
#Create dataframe with all duplicates
D2_ph_d=D2_ph[duplicated(D2_ph$Phage) | duplicated(D2_ph$Phage, fromLast=T),]
#Check
length(unique(D2_ph_d$Phage))#113, adds to 1398 so still good
#Focus dataframe on just meaningful columns
D2_ph_d=D2_ph_d[,-c(3,5:8)]
#Create dummy column and fill with NA
D2_ph_d$Con=NA
#Fill dummy column with Yes if crispr and nt match
D2_ph_d$Con[D2_ph_d$MAG_cr==D2_ph_d$MAG_nt] = "Yes"
#Fill dummy column with no if crisp and nt do not match
D2_ph_d$Con[D2_ph_d$MAG_cr!=D2_ph_d$MAG_nt] = "No"

#Extract Concencus dataframe for later combining
Concensus2=D2_ph_d%>%
  filter(Con=="Yes")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best)

#1 of 113 recovered
#Make sure if cases had a duplicate 'No' that we need to remove from the frame for downstream
D2_ph_d_non=D2_ph_d%>%
  group_by(Phage)%>%
  filter(!Phage%in%Concensus2$Phage)


#Pipeline to remove those that had concencus and then if not bring down to a single column with redundant hits
D2_ph_dup=D2_ph_d_non%>%
  group_by(Phage)%>%
  filter(is.na(Con))%>%
  mutate(Con=ifelse(is.na(MAG_nt),MAG_cr[],MAG_nt[]))%>%
  select(-MAG_cr,-MAG_nt) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(D2_ph_dup$Phage)) #90 of 113 recovered, that means No should have 22

#Pipeline to extract duplicates supported with multiple identical assignments and make nonredundant for later combining and to create vector for subsetting other dataframe
D2_ph_d2=D2_ph_dup%>%
  group_by(Phage)%>%
  filter(c(duplicated(Con)) | duplicated(Con, fromLast=T))%>%
  filter(duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)%>%
  filter(!duplicated(Phage))
length(unique(D2_ph_d2$Phage)) #32

#Create vector for subsetting
r2=as.vector(D2_ph_d2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
D2_ph_nd_2=D2_ph_dup%>%
  group_by(Phage)%>%
  filter(!Phage%in%r2)%>%
  filter(!duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)
length(unique(D2_ph_nd_2$Phage)) #58 recovered, D2_ph_nd_2+D2_ph_d2=90 recovered

#Pipeline to extract thos assigned 'No', use CRISPR as host, and remove redundancy
#First extract those supported multiple times by Crispr
Noncensus2=D2_ph_d_non%>%
  filter(Con=="No")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(Noncensus2$Phage)) #22

Noncensus2_2=Noncensus2%>%
  group_by(Phage)%>%
  filter(c(duplicated(best)) | duplicated(best, fromLast=T)) %>%
  filter(duplicated(Phage))%>%
  mutate(best=best)%>%
  select(-vid)%>%
  filter(!duplicated(Phage))
length(unique(Noncensus2_2$Phage)) #12

#Create another vector for subsetting
rnoncen2=as.vector(Noncensus2_2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
Noncensus2_3=Noncensus2%>%
  group_by(Phage)%>%
  filter(!Phage%in%rnoncen2)%>%
  filter(!duplicated(Phage))%>%
  select(-vid)
length(unique(Noncensus2_3$Phage)) #10

#Make original dataframe fit the columnns
D2_ph_nd_o=D2_ph_nd%>%
  select(Phage,best)
#Create final dataframe
D2_phage_host_final=bind_rows(D2_ph_nd_o,
                              D2_ph_d2,
                              D2_ph_nd_2,
                              Concensus2,
                              Noncensus2_2,
                              Noncensus2_3)
#Make hosts factor
D2_phage_host_final$best=as.factor(D2_phage_host_final$best)
#Rename best to Host
D2_phage_host_final=rename(D2_phage_host_final,Host=best)

#Add column with matID
D2_phage_host_final=D2_phage_host_final%>%
  mutate(Mat="D2")

#Ensure there are no duplicates
length(unique(D2_phage_host_final$Phage)) #1398/1398
duplicated(D2_phage_host_final$Phage)



#######
#######
#######
#######
######D3
#Filter spacer by mismatch ≤1 and e-value ≤1, then best bitscore
D3_crispr_filter=D3_crispr %>%
  filter(mismatch<=1 & e_value <=1) %>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="crispr")

#Order by phage ID
D3_crispr_filter=D3_crispr_filter[order(D3_crispr_filter$Phage),]

#Filter nt by ≥50 bitscore and ≤0.001 e value
D3_nt_filter=D3_nt%>%
  filter(bitscore>=50 & e_value <=0.001)

#Extract only highest bitscore for each unique phage ID
D3_nt_high=D3_nt_filter%>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="nt")

#Order by phage ID
D3_nt_high=D3_nt_high[order(D3_nt_high$Phage),]

#Filter d2star by common.kmer>36
D3_dstar_filter=D3_dstar%>%
  filter(common.kmers>36) %>%
  mutate(Source="kmer")

#Now create smaller dataframes from each for merging
D3_crisp_merger=D3_crispr_filter[,c(2,3,14)]
D3_nt_merger=D3_nt_high[,c(2,3,14)]
D3_dstar_merger=D3_dstar_filter[,c(1,2,6)]

#Now merge
D3_ph_merge1=merge(D3_crisp_merger,D3_nt_merger, by="Phage",all=T)
D3_ph_merge_2=merge(D3_ph_merge1,D3_dstar_merger,by="Phage",all=T)

#Now create new column with single best prediction
D3_ph=D3_ph_merge_2
D3_ph=rename(D3_ph,MAG_cr=MAG.x)
D3_ph=rename(D3_ph,MAG_nt=MAG.y)
D3_ph$MAG_cr=as.character(D3_ph$MAG_cr)
D3_ph$MAG_nt=as.character(D3_ph$MAG_nt)
D3_ph$MAG=as.character(D3_ph$MAG)
D3_ph=D3_ph%>%
  mutate(best=ifelse(is.na(MAG_cr) & is.na(MAG_nt),MAG[],
                     ifelse(is.na(MAG_cr),MAG_nt[],MAG_cr[])))

#Get length of unique phages to make sure all persist to end
length(unique(D3_ph$Phage))#1724

##Eliminate redundancy based on conditions
#Create dataframe without duplicates for later combining
D3_ph_nd=D3_ph[!c(duplicated(D3_ph$Phage) | duplicated(D3_ph$Phage, fromLast=T)),]
#Check
length(unique(D3_ph_nd$Phage))#1554
#Create dataframe with all duplicates
D3_ph_d=D3_ph[duplicated(D3_ph$Phage) | duplicated(D3_ph$Phage, fromLast=T),]
#Check
length(unique(D3_ph_d$Phage))#170, adds to 1724 so still good
#Focus dataframe on just meaningful columns
D3_ph_d=D3_ph_d[,-c(3,5:8)]
#Create dummy column and fill with NA
D3_ph_d$Con=NA
#Fill dummy column with Yes if crispr and nt match
D3_ph_d$Con[D3_ph_d$MAG_cr==D3_ph_d$MAG_nt] = "Yes"
#Fill dummy column with no if crisp and nt do not match
D3_ph_d$Con[D3_ph_d$MAG_cr!=D3_ph_d$MAG_nt] = "No"

#Extract Concencus dataframe for later combining
Concensus3=D3_ph_d%>%
  filter(Con=="Yes")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best)

#3 of 170 recovered, but this one has a redundandt so clean and 2 recovered for counting
Concensus3=Concensus3%>%
  filter(!duplicated(Phage))
#Make sure if cases had a duplicate 'No' that we need to remove from the frame for downstream
D3_ph_d_non=D3_ph_d%>%
  group_by(Phage)%>%
  filter(!Phage%in%Concensus3$Phage)


#Pipeline to remove those that had concencus and then if not bring down to a single column with redundant hits
D3_ph_dup=D3_ph_d_non%>%
  group_by(Phage)%>%
  filter(is.na(Con))%>%
  mutate(Con=ifelse(is.na(MAG_nt),MAG_cr[],MAG_nt[]))%>%
  select(-MAG_cr,-MAG_nt) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(D3_ph_dup$Phage)) #129 of 170 recovered, that means No should have 39

#Pipeline to extract duplicates supported with multiple identical assignments and make nonredundant for later combining and to create vector for subsetting other dataframe
D3_ph_d2=D3_ph_dup%>%
  group_by(Phage)%>%
  filter(c(duplicated(Con)) | duplicated(Con, fromLast=T))%>%
  filter(duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)%>%
  filter(!duplicated(Phage))
length(unique(D3_ph_d2$Phage)) #41

#Create vector for subsetting
r3=as.vector(D3_ph_d2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
D3_ph_nd_2=D3_ph_dup%>%
  group_by(Phage)%>%
  filter(!Phage%in%r3)%>%
  filter(!duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)
length(unique(D3_ph_nd_2$Phage)) #88 recovered, D3_ph_nd_2+D3_ph_d2=129 recovered

#Pipeline to extract thos assigned 'No', use CRISPR as host, and remove redundancy
#First extract those supported multiple times by Crispr
Noncensus3=D3_ph_d_non%>%
  filter(Con=="No")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(Noncensus3$Phage)) #39

Noncensus3_2=Noncensus3%>%
  group_by(Phage)%>%
  filter(c(duplicated(best)) | duplicated(best, fromLast=T)) %>%
  filter(duplicated(Phage))%>%
  mutate(best=best)%>%
  select(-vid)%>%
  filter(!duplicated(Phage))
length(unique(Noncensus3_2$Phage)) #18

#Create another vector for subsetting
rnoncen3=as.vector(Noncensus3_2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
Noncensus3_3=Noncensus3%>%
  group_by(Phage)%>%
  filter(!Phage%in%rnoncen3)%>%
  filter(!duplicated(Phage))%>%
  select(-vid)
length(unique(Noncensus3_3$Phage)) #21

#Make original dataframe fit the columnns
D3_ph_nd_o=D3_ph_nd%>%
  select(Phage,best)
#Create final dataframe
D3_phage_host_final=bind_rows(D3_ph_nd_o,
                              D3_ph_d2,
                              D3_ph_nd_2,
                              Concensus3,
                              Noncensus3_2,
                              Noncensus3_3)
#Make hosts factor
D3_phage_host_final$best=as.factor(D3_phage_host_final$best)
#Rename best to Host
D3_phage_host_final=rename(D3_phage_host_final,Host=best)

#Add column with matID
D3_phage_host_final=D3_phage_host_final%>%
  mutate(Mat="D3")

#Ensure there are no duplicates
length(unique(D3_phage_host_final$Phage)) #1724/1724
duplicated(D3_phage_host_final$Phage)



#######
#######
#######
#######
######D4
#Filter spacer by mismatch ≤1 and e-value ≤1, then best bitscore
D4_crispr_filter=D4_crispr %>%
  filter(mismatch<=1 & e_value <=1) %>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="crispr")

#Order by phage ID
D4_crispr_filter=D4_crispr_filter[order(D4_crispr_filter$Phage),]

#Filter nt by ≥50 bitscore and ≤0.001 e value
D4_nt_filter=D4_nt%>%
  filter(bitscore>=50 & e_value <=0.001)

#Extract only highest bitscore for each unique phage ID
D4_nt_high=D4_nt_filter%>%
  group_by(Phage) %>%
  filter(bitscore==max(bitscore)) %>%
  ungroup() %>%
  mutate(Source="nt")

#Order by phage ID
D4_nt_high=D4_nt_high[order(D4_nt_high$Phage),]

#Filter d2star by common.kmer>36
D4_dstar_filter=D4_dstar%>%
  filter(common.kmers>36) %>%
  mutate(Source="kmer")

#Now create smaller dataframes from each for merging
D4_crisp_merger=D4_crispr_filter[,c(2,3,14)]
D4_nt_merger=D4_nt_high[,c(2,3,14)]
D4_dstar_merger=D4_dstar_filter[,c(1,2,6)]

#Now merge
D4_ph_merge1=merge(D4_crisp_merger,D4_nt_merger, by="Phage",all=T)
D4_ph_merge_2=merge(D4_ph_merge1,D4_dstar_merger,by="Phage",all=T)

#Now create new column with single best prediction
D4_ph=D4_ph_merge_2
D4_ph=rename(D4_ph,MAG_cr=MAG.x)
D4_ph=rename(D4_ph,MAG_nt=MAG.y)
D4_ph$MAG_cr=as.character(D4_ph$MAG_cr)
D4_ph$MAG_nt=as.character(D4_ph$MAG_nt)
D4_ph$MAG=as.character(D4_ph$MAG)
D4_ph=D4_ph%>%
  mutate(best=ifelse(is.na(MAG_cr) & is.na(MAG_nt),MAG[],
                     ifelse(is.na(MAG_cr),MAG_nt[],MAG_cr[])))

#Get length of unique phages to make sure all persist to end
length(unique(D4_ph$Phage))#1205

##Eliminate redundancy based on conditions
#Create dataframe without duplicates for later combining
D4_ph_nd=D4_ph[!c(duplicated(D4_ph$Phage) | duplicated(D4_ph$Phage, fromLast=T)),]
#Check
length(unique(D4_ph_nd$Phage))#1088
#Create dataframe with all duplicates
D4_ph_d=D4_ph[duplicated(D4_ph$Phage) | duplicated(D4_ph$Phage, fromLast=T),]
#Check
length(unique(D4_ph_d$Phage))#117, adds to 1205 so still good
#Focus dataframe on just meaningful columns
D4_ph_d=D4_ph_d[,-c(3,5:8)]
#Create dummy column and fill with NA
D4_ph_d$Con=NA
#Fill dummy column with Yes if crispr and nt match
D4_ph_d$Con[D4_ph_d$MAG_cr==D4_ph_d$MAG_nt] = "Yes"
#Fill dummy column with no if crisp and nt do not match
D4_ph_d$Con[D4_ph_d$MAG_cr!=D4_ph_d$MAG_nt] = "No"

#Extract Concencus dataframe for later combining
Concensus4=D4_ph_d%>%
  filter(Con=="Yes")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best)

#6 of 170 recovered, but this one has a redundandt so clean and 5 recovered for counting
Concensus4=Concensus4%>%
  filter(!duplicated(Phage))
#Make sure if cases had a duplicate 'No' that we need to remove from the frame for downstream
D4_ph_d_non=D4_ph_d%>%
  group_by(Phage)%>%
  filter(!Phage%in%Concensus4$Phage)


#Pipeline to remove those that had concencus and then if not bring down to a single column with redundant hits
D4_ph_dup=D4_ph_d_non%>%
  group_by(Phage)%>%
  filter(is.na(Con))%>%
  mutate(Con=ifelse(is.na(MAG_nt),MAG_cr[],MAG_nt[]))%>%
  select(-MAG_cr,-MAG_nt) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(D4_ph_dup$Phage)) #93 of 117 recovered, that means No should have 19

#Pipeline to extract duplicates supported with multiple identical assignments and make nonredundant for later combining and to create vector for subsetting other dataframe
D4_ph_d2=D4_ph_dup%>%
  group_by(Phage)%>%
  filter(c(duplicated(Con)) | duplicated(Con, fromLast=T))%>%
  filter(duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)%>%
  filter(!duplicated(Phage))
length(unique(D4_ph_d2$Phage)) #32

#Create vector for subsetting
r4=as.vector(D4_ph_d2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
D4_ph_nd_2=D4_ph_dup%>%
  group_by(Phage)%>%
  filter(!Phage%in%r4)%>%
  filter(!duplicated(Phage))%>%
  mutate(best=Con)%>%
  select(-Con,-vid)
length(unique(D4_ph_nd_2$Phage)) #61 recovered, D4_ph_nd_2+D4_ph_d2=93 recovered

#Pipeline to extract thos assigned 'No', use CRISPR as host, and remove redundancy
#First extract those supported multiple times by Crispr
Noncensus4=D4_ph_d_non%>%
  filter(Con=="No")%>%
  mutate(best=MAG_cr)%>%
  select(Phage,best) %>%
  mutate(ID=row_number(),
         V="v")%>%
  unite(vid,c(V,ID),sep="")
length(unique(Noncensus4$Phage)) #19

Noncensus4_2=Noncensus4%>%
  group_by(Phage)%>%
  filter(c(duplicated(best)) | duplicated(best, fromLast=T)) %>%
  filter(duplicated(Phage))%>%
  mutate(best=best)%>%
  select(-vid)%>%
  filter(!duplicated(Phage))
length(unique(Noncensus4_2$Phage)) #7

#Create another vector for subsetting
rnoncen4=as.vector(Noncensus4_2$Phage)

#Pipeline to remove those already in above dataframe and then pick the first of remaining duplicates
Noncensus4_3=Noncensus4%>%
  group_by(Phage)%>%
  filter(!Phage%in%rnoncen4)%>%
  filter(!duplicated(Phage))%>%
  select(-vid)
length(unique(Noncensus4_3$Phage)) #12

#Make original dataframe fit the columnns
D4_ph_nd_o=D4_ph_nd%>%
  select(Phage,best)
#Create final dataframe
D4_phage_host_final=bind_rows(D4_ph_nd_o,
                              D4_ph_d2,
                              D4_ph_nd_2,
                              Concensus4,
                              Noncensus4_2,
                              Noncensus4_3)
#Make hosts factor
D4_phage_host_final$best=as.factor(D4_phage_host_final$best)
#Rename best to Host
D4_phage_host_final=rename(D4_phage_host_final,Host=best)

#Add column with matID
D4_phage_host_final=D4_phage_host_final%>%
  mutate(Mat="D4")

#Ensure there are no duplicates
length(unique(D4_phage_host_final$Phage)) #1205/1205
duplicated(D4_phage_host_final$Phage)

##A total of 5539 host-phage pairings / 13,026 phages
5539/13026
#42.5% Assignment
##############################################
########################
#Now we have extracted the single best host per phage
########################
##############################################
########################
########################
##############################################
#Read in host taxonomy
D1_host_tax=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Host_Phage/Host_TAX/D1_gtdbtk.txt",header=T)
D2_host_tax=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Host_Phage/Host_TAX/D2_gtdbtk.txt",header=T)
D3_host_tax=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Host_Phage/Host_TAX/D3_gtdbtk.txt",header=T)
D4_host_tax=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Host_Phage/Host_TAX/D4_gtdbtk.txt",header=T)

#Read in host counts
D1_host_counts=read.xls("/Users/cissell/Desktop/D1_MAG_counts.xlsx")
D2_host_counts=read.xls("/Users/cissell/Desktop/D2_MAG_counts.xlsx")
D3_host_counts=read.xls("/Users/cissell/Desktop/D3_MAG_counts.xlsx")
D4_host_counts=read.xls("/Users/cissell/Desktop/D4_MAG_counts.xlsx")

#Clean up <phylum> from host count frames
D1_host_counts_clean=D1_host_counts%>%
  select(-Phylum)
D2_host_counts_clean=D2_host_counts%>%
  select(-Phylum)
D3_host_counts_clean=D3_host_counts%>%
  select(-Phylum)
D4_host_counts_clean=D4_host_counts%>%
  select(-Phylum)

#Read in Phage counts
D1_Phage_counts=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Host_Phage/ABUNDANCE/D1_phage_counts.tsv",header=T)
D2_Phage_counts=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Host_Phage/ABUNDANCE/D2_phage_counts.tsv",header=T)
D3_Phage_counts=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Host_Phage/ABUNDANCE/D3_phage_counts.tsv",header=T)
D4_Phage_counts=read.table("/Users/cissell/Desktop/VIRSORTER2_Diel/Host_Phage/ABUNDANCE/D4_phage_counts.tsv",header=T)


#Clean up phage taxonomy from other companion Rscript for merging
phage_tax_import=rename(viral_clust_curate_tax,Phage=contig_id)

###NOW LETS BUILD THE PER MAT FRAMES WITH PAIRS,TAX,COUNTS,and LENGTHS and calculate TPM
##D1
#Merge D1 with taxonomy of phage and host
D1_tax_only=merge(D1_phage_host_final,phage_tax_import,by="Phage",all.x=T)
D1_tax_only=D1_tax_only%>%
  select(Phage,Host,mergetax)
D1_tax_only=merge(D1_tax_only,D1_host_tax,by="Host",all.x=T)

#Add in phage counts
D1_tax_counts1=merge(D1_tax_only,D1_Phage_counts,by="Phage") #Remeber this will decrease becasue RNA phages will be removed (we want this, for now)
#Add in host counts
D1_tax_counts2=merge(D1_tax_counts1,D1_host_counts_clean,by="Host",all.x=T)

#Make new df
rpk1=D1_tax_counts2
str(rpk1)
#Put length in KB
rpk1$PLength=rpk1$PLength/1000
rpk1$HLength=rpk1$HLength/1000
#Create RPK for Mat 1
rpk1$P1_1=rpk1$P1_1/rpk1$PLength
rpk1$H1_1=rpk1$H1_1/rpk1$HLength
rpk1$P1_3=rpk1$P1_3/rpk1$PLength
rpk1$H1_3=rpk1$H1_3/rpk1$HLength
rpk1$P1_5=rpk1$P1_5/rpk1$PLength
rpk1$H1_5=rpk1$H1_5/rpk1$HLength
#Need to extract a new dataframe of only unique so each MAG only occurs once
hpk1=rpk1[!duplicated(rpk1$Host),]
#create scaling factors
#Mat 1
sum11=((sum(rpk1$P1_1)+sum(hpk1$H1_1))/1000000)
sum13=((sum(rpk1$P1_3)+sum(hpk1$H1_3))/1000000)
sum15=((sum(rpk1$P1_5)+sum(hpk1$H1_5))/1000000)
#Create TPM for Mat 1
rpk1$P1_1=rpk1$P1_1/sum11
rpk1$H1_1=rpk1$H1_1/sum11
rpk1$P1_3=rpk1$P1_3/sum13
rpk1$H1_3=rpk1$H1_3/sum13
rpk1$P1_5=rpk1$P1_5/sum15
rpk1$H1_5=rpk1$H1_5/sum15
#Wont make sense just yet becuase of duplicates

#Bring it back to a sensical name now
D1_count_tpm=rpk1

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
D1_host_lost=D1_count_tpm%>%
  select(c(Host,H1_1,H1_3,H1_5,Phylum,Class))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D1_count_tpm2=D1_count_tpm%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P1_1","P1_3","P1_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D1_count_merger=merge(D1_count_tpm2,D1_host_lost,by="Host",all=T)
D1_count_merger$Phylum=as.factor(D1_count_merger$Phylum)
D1_count_merger$Class=as.factor(D1_count_merger$Class)
str(D1_count_merger)
#Summarize within phyla
D1_phyla_counts=D1_count_merger%>%
  select(-c(Host,Class))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D1_phyla_counts$P1_1)+sum(D1_phyla_counts$H1_1)
sum(D1_phyla_counts$P1_3)+sum(D1_phyla_counts$H1_3)
sum(D1_phyla_counts$P1_5)+sum(D1_phyla_counts$H1_5)
#Nice

#Create each individual VMR and logVMR for each pair
D1_count_vmrs=D1_phyla_counts%>%
  mutate(VMR1=P1_1/H1_1)%>%
  mutate(logVMR1=log10(P1_1/H1_1))%>%
  mutate(VMR3=P1_3/H1_3)%>%
  mutate(logVMR3=log10(P1_3/H1_3))%>%
  mutate(VMR5=P1_5/H1_5)%>%
  mutate(logVMR5=log10(P1_5/H1_5))%>%
  rename(P_1=P1_1,
         P_3=P1_3,
         P_5=P1_5,
         H_1=H1_1,
         H_3=H1_3,
         H_5=H1_5,
         VMR_1=VMR1,
         VMR_3=VMR3,
         VMR_5=VMR5,
         logVMR_1=logVMR1,
         logVMR_3=logVMR3,
         logVMR_5=logVMR5)%>%
  mutate(Mat="D1")

#Pivot longer for later merging and plotting
D1_count_vmr_long=D1_count_vmrs%>%
  pivot_longer(
    cols=-c(Phylum,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D1_count_vmr_long=as.data.frame(D1_count_vmr_long)
D1_count_vmr_long$Mat=as.factor(D1_count_vmr_long$Mat)
D1_count_vmr_long$Time=as.factor(D1_count_vmr_long$Time)
str(D1_count_vmr_long)

#Summarize within class
D1_class_counts=D1_count_merger%>%
  select(-c(Host,Phylum))%>%
  group_by(Class)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D1_class_counts$P1_1)+sum(D1_class_counts$H1_1)
sum(D1_class_counts$P1_3)+sum(D1_class_counts$H1_3)
sum(D1_class_counts$P1_5)+sum(D1_class_counts$H1_5)
#Nice

#Create each individual VMR and logVMR for each pair
D1_class_count_vmrs=D1_class_counts%>%
  mutate(VMR1=P1_1/H1_1)%>%
  mutate(logVMR1=log10(P1_1/H1_1))%>%
  mutate(VMR3=P1_3/H1_3)%>%
  mutate(logVMR3=log10(P1_3/H1_3))%>%
  mutate(VMR5=P1_5/H1_5)%>%
  mutate(logVMR5=log10(P1_5/H1_5))%>%
  rename(P_1=P1_1,
         P_3=P1_3,
         P_5=P1_5,
         H_1=H1_1,
         H_3=H1_3,
         H_5=H1_5,
         VMR_1=VMR1,
         VMR_3=VMR3,
         VMR_5=VMR5,
         logVMR_1=logVMR1,
         logVMR_3=logVMR3,
         logVMR_5=logVMR5)%>%
  mutate(Mat="D1")

#Pivot longer for later merging and plotting
D1_class_count_vmr_long=D1_class_count_vmrs%>%
  pivot_longer(
    cols=-c(Class,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D1_class_count_vmr_long=as.data.frame(D1_class_count_vmr_long)
D1_class_count_vmr_long$Mat=as.factor(D1_class_count_vmr_long$Mat)
D1_class_count_vmr_long$Time=as.factor(D1_class_count_vmr_long$Time)
str(D1_class_count_vmr_long)

###################
##D2
###################
#Merge D2 with taxonomy of phage and host
D2_tax_only=merge(D2_phage_host_final,phage_tax_import,by="Phage",all.x=T)
D2_tax_only=D2_tax_only%>%
  select(Phage,Host,mergetax)
D2_tax_only=merge(D2_tax_only,D2_host_tax,by="Host",all.x=T)

#Add in phage counts
D2_tax_counts1=merge(D2_tax_only,D2_Phage_counts,by="Phage") #Remeber this will decrease becasue RNA phages will be removed (we want this, for now)
#Add in host counts
D2_tax_counts2=merge(D2_tax_counts1,D2_host_counts_clean,by="Host",all.x=T)

#Make new df
rpk2=D2_tax_counts2
str(rpk2)
#Put length in KB
rpk2$PLength=rpk2$PLength/1000
rpk2$HLength=rpk2$HLength/1000
#Create RPK for Mat 1
rpk2$P2_1=rpk2$P2_1/rpk2$PLength
rpk2$H2_1=rpk2$H2_1/rpk2$HLength
rpk2$P2_3=rpk2$P2_3/rpk2$PLength
rpk2$H2_3=rpk2$H2_3/rpk2$HLength
rpk2$P2_5=rpk2$P2_5/rpk2$PLength
rpk2$H2_5=rpk2$H2_5/rpk2$HLength
#Need to extract a new dataframe of only unique so each MAG only occurs once
hpk2=rpk2[!duplicated(rpk2$Host),]
#create scaling factors
#Mat 1
sum21=((sum(rpk2$P2_1)+sum(hpk2$H2_1))/1000000)
sum23=((sum(rpk2$P2_3)+sum(hpk2$H2_3))/1000000)
sum25=((sum(rpk2$P2_5)+sum(hpk2$H2_5))/1000000)
#Create TPM for Mat 1
rpk2$P2_1=rpk2$P2_1/sum21
rpk2$H2_1=rpk2$H2_1/sum21
rpk2$P2_3=rpk2$P2_3/sum23
rpk2$H2_3=rpk2$H2_3/sum23
rpk2$P2_5=rpk2$P2_5/sum25
rpk2$H2_5=rpk2$H2_5/sum25
#Wont make sense just yet becuase of duplicates

#Bring it back to a sensical name now
D2_count_tpm=rpk2

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
D2_host_lost=D2_count_tpm%>%
  select(c(Host,H2_1,H2_3,H2_5,Phylum,Class))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D2_count_tpm2=D2_count_tpm%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P2_1","P2_3","P2_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D2_count_merger=merge(D2_count_tpm2,D2_host_lost,by="Host",all=T)
D2_count_merger$Phylum=as.factor(D2_count_merger$Phylum)
D2_count_merger$Class=as.factor(D2_count_merger$Class)
str(D2_count_merger)
#Summarize within phyla
D2_phyla_counts=D2_count_merger%>%
  select(-c(Host,Class))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D2_phyla_counts$P2_1)+sum(D2_phyla_counts$H2_1)
sum(D2_phyla_counts$P2_3)+sum(D2_phyla_counts$H2_3)
sum(D2_phyla_counts$P2_5)+sum(D2_phyla_counts$H2_5)
#Nice

#Create each individual VMR and logVMR for each pair
D2_count_vmrs=D2_phyla_counts%>%
  mutate(VMR1=P2_1/H2_1)%>%
  mutate(logVMR1=log10(P2_1/H2_1))%>%
  mutate(VMR3=P2_3/H2_3)%>%
  mutate(logVMR3=log10(P2_3/H2_3))%>%
  mutate(VMR5=P2_5/H2_5)%>%
  mutate(logVMR5=log10(P2_5/H2_5))%>%
  rename(P_1=P2_1,
         P_3=P2_3,
         P_5=P2_5,
         H_1=H2_1,
         H_3=H2_3,
         H_5=H2_5,
         VMR_1=VMR1,
         VMR_3=VMR3,
         VMR_5=VMR5,
         logVMR_1=logVMR1,
         logVMR_3=logVMR3,
         logVMR_5=logVMR5)%>%
  mutate(Mat="D2")

#Pivot longer for later merging and plotting
D2_count_vmr_long=D2_count_vmrs%>%
  pivot_longer(
    cols=-c(Phylum,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D2_count_vmr_long=as.data.frame(D2_count_vmr_long)
D2_count_vmr_long$Mat=as.factor(D2_count_vmr_long$Mat)
D2_count_vmr_long$Time=as.factor(D2_count_vmr_long$Time)
str(D2_count_vmr_long)

#Summarize within class
D2_class_counts=D2_count_merger%>%
  select(-c(Host,Phylum))%>%
  group_by(Class)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D2_class_counts$P2_1)+sum(D2_class_counts$H2_1)
sum(D2_class_counts$P2_3)+sum(D2_class_counts$H2_3)
sum(D2_class_counts$P2_5)+sum(D2_class_counts$H2_5)
#Nice

#Create each individual VMR and logVMR for each pair
D2_class_count_vmrs=D2_class_counts%>%
  mutate(VMR1=P2_1/H2_1)%>%
  mutate(logVMR1=log10(P2_1/H2_1))%>%
  mutate(VMR3=P2_3/H2_3)%>%
  mutate(logVMR3=log10(P2_3/H2_3))%>%
  mutate(VMR5=P2_5/H2_5)%>%
  mutate(logVMR5=log10(P2_5/H2_5))%>%
  rename(P_1=P2_1,
         P_3=P2_3,
         P_5=P2_5,
         H_1=H2_1,
         H_3=H2_3,
         H_5=H2_5,
         VMR_1=VMR1,
         VMR_3=VMR3,
         VMR_5=VMR5,
         logVMR_1=logVMR1,
         logVMR_3=logVMR3,
         logVMR_5=logVMR5)%>%
  mutate(Mat="D2")

#Pivot longer for later merging and plotting
D2_class_count_vmr_long=D2_class_count_vmrs%>%
  pivot_longer(
    cols=-c(Class,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D2_class_count_vmr_long=as.data.frame(D2_class_count_vmr_long)
D2_class_count_vmr_long$Mat=as.factor(D2_class_count_vmr_long$Mat)
D2_class_count_vmr_long$Time=as.factor(D2_class_count_vmr_long$Time)
str(D2_class_count_vmr_long)
###################
##D3
###################
#Merge D3 with taxonomy of phage and host
D3_tax_only=merge(D3_phage_host_final,phage_tax_import,by="Phage",all.x=T)
D3_tax_only=D3_tax_only%>%
  select(Phage,Host,mergetax)
D3_tax_only=merge(D3_tax_only,D3_host_tax,by="Host",all.x=T)

#Add in phage counts
D3_tax_counts1=merge(D3_tax_only,D3_Phage_counts,by="Phage") #Remeber this will decrease becasue RNA phages will be removed (we want this, for now)
#Add in host counts
D3_tax_counts2=merge(D3_tax_counts1,D3_host_counts_clean,by="Host",all.x=T)

#Make new df
rpk3=D3_tax_counts2
str(rpk3)
#Put length in KB
rpk3$PLength=rpk3$PLength/1000
rpk3$HLength=rpk3$HLength/1000
#Create RPK for Mat 1
rpk3$P3_1=rpk3$P3_1/rpk3$PLength
rpk3$H3_1=rpk3$H3_1/rpk3$HLength
rpk3$P3_3=rpk3$P3_3/rpk3$PLength
rpk3$H3_3=rpk3$H3_3/rpk3$HLength
rpk3$P3_5=rpk3$P3_5/rpk3$PLength
rpk3$H3_5=rpk3$H3_5/rpk3$HLength
#Need to extract a new dataframe of only unique so each MAG only occurs once
hpk3=rpk3[!duplicated(rpk3$Host),]
#create scaling factors
#Mat 1
sum31=((sum(rpk3$P3_1)+sum(hpk3$H3_1))/1000000)
sum33=((sum(rpk3$P3_3)+sum(hpk3$H3_3))/1000000)
sum35=((sum(rpk3$P3_5)+sum(hpk3$H3_5))/1000000)
#Create TPM for Mat 1
rpk3$P3_1=rpk3$P3_1/sum31
rpk3$H3_1=rpk3$H3_1/sum31
rpk3$P3_3=rpk3$P3_3/sum33
rpk3$H3_3=rpk3$H3_3/sum33
rpk3$P3_5=rpk3$P3_5/sum35
rpk3$H3_5=rpk3$H3_5/sum35
#Wont make sense just yet becuase of duplicates

#Bring it back to a sensical name now
D3_count_tpm=rpk3

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
D3_host_lost=D3_count_tpm%>%
  select(c(Host,H3_1,H3_3,H3_5,Phylum,Class))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D3_count_tpm2=D3_count_tpm%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P3_1","P3_3","P3_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D3_count_merger=merge(D3_count_tpm2,D3_host_lost,by="Host",all=T)
D3_count_merger$Phylum=as.factor(D3_count_merger$Phylum)
D3_count_merger=D3_count_merger%>%
  mutate_at("Phylum",str_replace,"_F","")
str(D3_count_merger)
#Summarize within phyla
D3_phyla_counts=D3_count_merger%>%
  select(-c(Host,Class))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D3_phyla_counts$P3_1)+sum(D3_phyla_counts$H3_1)
sum(D3_phyla_counts$P3_3)+sum(D3_phyla_counts$H3_3)
sum(D3_phyla_counts$P3_5)+sum(D3_phyla_counts$H3_5)
#Nice

#Create each individual VMR and logVMR for each pair
D3_count_vmrs=D3_phyla_counts%>%
  mutate(VMR1=P3_1/H3_1)%>%
  mutate(logVMR1=log10(P3_1/H3_1))%>%
  mutate(VMR3=P3_3/H3_3)%>%
  mutate(logVMR3=log10(P3_3/H3_3))%>%
  mutate(VMR5=P3_5/H3_5)%>%
  mutate(logVMR5=log10(P3_5/H3_5))%>%
  rename(P_1=P3_1,
         P_3=P3_3,
         P_5=P3_5,
         H_1=H3_1,
         H_3=H3_3,
         H_5=H3_5,
         VMR_1=VMR1,
         VMR_3=VMR3,
         VMR_5=VMR5,
         logVMR_1=logVMR1,
         logVMR_3=logVMR3,
         logVMR_5=logVMR5)%>%
  mutate(Mat="D3")

#Pivot longer for later merging and plotting
D3_count_vmr_long=D3_count_vmrs%>%
  pivot_longer(
    cols=-c(Phylum,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D3_count_vmr_long=as.data.frame(D3_count_vmr_long)
D3_count_vmr_long$Mat=as.factor(D3_count_vmr_long$Mat)
D3_count_vmr_long$Time=as.factor(D3_count_vmr_long$Time)
str(D3_count_vmr_long)

#Summarize within class
D3_class_counts=D3_count_merger%>%
  select(-c(Host,Phylum))%>%
  group_by(Class)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D3_class_counts$P3_1)+sum(D3_class_counts$H3_1)
sum(D3_class_counts$P3_3)+sum(D3_class_counts$H3_3)
sum(D3_class_counts$P3_5)+sum(D3_class_counts$H3_5)
#Nice

#Create each individual VMR and logVMR for each pair
D3_class_count_vmrs=D3_class_counts%>%
  mutate(VMR1=P3_1/H3_1)%>%
  mutate(logVMR1=log10(P3_1/H3_1))%>%
  mutate(VMR3=P3_3/H3_3)%>%
  mutate(logVMR3=log10(P3_3/H3_3))%>%
  mutate(VMR5=P3_5/H3_5)%>%
  mutate(logVMR5=log10(P3_5/H3_5))%>%
  dplyr::rename(P_1=P3_1,
                P_3=P3_3,
                P_5=P3_5,
                H_1=H3_1,
                H_3=H3_3,
                H_5=H3_5,
                VMR_1=VMR1,
                VMR_3=VMR3,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_3=logVMR3,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D3")

#Pivot longer for later merging and plotting
D3_class_count_vmr_long=D3_class_count_vmrs%>%
  pivot_longer(
    cols=-c(Class,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D3_class_count_vmr_long=as.data.frame(D3_class_count_vmr_long)
D3_class_count_vmr_long$Mat=as.factor(D3_class_count_vmr_long$Mat)
D3_class_count_vmr_long$Time=as.factor(D3_class_count_vmr_long$Time)
str(D3_class_count_vmr_long)

###################
##D4
###################
#Merge D4 with taxonomy of phage and host
D4_tax_only=merge(D4_phage_host_final,phage_tax_import,by="Phage",all.x=T)
D4_tax_only=D4_tax_only%>%
  select(Phage,Host,mergetax)
D4_tax_only=merge(D4_tax_only,D4_host_tax,by="Host",all.x=T)

#Add in phage counts
D4_tax_counts1=merge(D4_tax_only,D4_Phage_counts,by="Phage") #Remeber this will decrease becasue RNA phages will be removed (we want this, for now)
#Add in host counts
D4_tax_counts2=merge(D4_tax_counts1,D4_host_counts_clean,by="Host",all.x=T)

#Make new df
rpk4=D4_tax_counts2
str(rpk4)
#Put length in KB
rpk4$PLength=rpk4$PLength/1000
rpk4$HLength=rpk4$HLength/1000
#Create RPK for Mat 1
rpk4$P4_1=rpk4$P4_1/rpk4$PLength
rpk4$H4_1=rpk4$H4_1/rpk4$HLength
rpk4$P4_3=rpk4$P4_3/rpk4$PLength
rpk4$H4_3=rpk4$H4_3/rpk4$HLength
rpk4$P4_5=rpk4$P4_5/rpk4$PLength
rpk4$H4_5=rpk4$H4_5/rpk4$HLength
#Need to extract a new dataframe of only unique so each MAG only occurs once
hpk4=rpk4[!duplicated(rpk4$Host),]
#create scaling factors
#Mat 1
sum41=((sum(rpk4$P4_1)+sum(hpk4$H4_1))/1000000)
sum43=((sum(rpk4$P4_3)+sum(hpk4$H4_3))/1000000)
sum45=((sum(rpk4$P4_5)+sum(hpk4$H4_5))/1000000)
#Create TPM for Mat 1
rpk4$P4_1=rpk4$P4_1/sum41
rpk4$H4_1=rpk4$H4_1/sum41
rpk4$P4_3=rpk4$P4_3/sum43
rpk4$H4_3=rpk4$H4_3/sum43
rpk4$P4_5=rpk4$P4_5/sum45
rpk4$H4_5=rpk4$H4_5/sum45
#Wont make sense just yet becuase of duplicates

#Bring it back to a sensical name now
D4_count_tpm=rpk4

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
D4_host_lost=D4_count_tpm%>%
  select(c(Host,H4_1,H4_3,H4_5,Phylum,Class))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D4_count_tpm2=D4_count_tpm%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P4_1","P4_3","P4_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D4_count_merger=merge(D4_count_tpm2,D4_host_lost,by="Host",all=T)
D4_count_merger$Phylum=as.factor(D4_count_merger$Phylum)
str(D4_count_merger)
#Summarize within phyla
D4_phyla_counts=D4_count_merger%>%
  select(-c(Host,Class))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D4_phyla_counts$P4_1)+sum(D4_phyla_counts$H4_1)
sum(D4_phyla_counts$P4_3)+sum(D4_phyla_counts$H4_3)
sum(D4_phyla_counts$P4_5)+sum(D4_phyla_counts$H4_5)
#Nice

#Create each individual VMR and logVMR for each pair
D4_count_vmrs=D4_phyla_counts%>%
  mutate(VMR1=P4_1/H4_1)%>%
  mutate(logVMR1=log10(P4_1/H4_1))%>%
  mutate(VMR3=P4_3/H4_3)%>%
  mutate(logVMR3=log10(P4_3/H4_3))%>%
  mutate(VMR5=P4_5/H4_5)%>%
  mutate(logVMR5=log10(P4_5/H4_5))%>%
  rename(P_1=P4_1,
         P_3=P4_3,
         P_5=P4_5,
         H_1=H4_1,
         H_3=H4_3,
         H_5=H4_5,
         VMR_1=VMR1,
         VMR_3=VMR3,
         VMR_5=VMR5,
         logVMR_1=logVMR1,
         logVMR_3=logVMR3,
         logVMR_5=logVMR5)%>%
  mutate(Mat="D4")

#Pivot longer for later merging and plotting
D4_count_vmr_long=D4_count_vmrs%>%
  pivot_longer(
    cols=-c(Phylum,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D4_count_vmr_long=as.data.frame(D4_count_vmr_long)
D4_count_vmr_long$Mat=as.factor(D4_count_vmr_long$Mat)
D4_count_vmr_long$Time=as.factor(D4_count_vmr_long$Time)
str(D4_count_vmr_long)

#Summarize within class
D4_class_counts=D4_count_merger%>%
  select(-c(Host,Phylum))%>%
  group_by(Class)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D4_class_counts$P4_1)+sum(D4_class_counts$H4_1)
sum(D4_class_counts$P4_3)+sum(D4_class_counts$H4_3)
sum(D4_class_counts$P4_5)+sum(D4_class_counts$H4_5)
#Nice

#Create each individual VMR and logVMR for each pair
D4_class_count_vmrs=D4_class_counts%>%
  mutate(VMR1=P4_1/H4_1)%>%
  mutate(logVMR1=log10(P4_1/H4_1))%>%
  mutate(VMR3=P4_3/H4_3)%>%
  mutate(logVMR3=log10(P4_3/H4_3))%>%
  mutate(VMR5=P4_5/H4_5)%>%
  mutate(logVMR5=log10(P4_5/H4_5))%>%
  dplyr::rename(P_1=P4_1,
                P_3=P4_3,
                P_5=P4_5,
                H_1=H4_1,
                H_3=H4_3,
                H_5=H4_5,
                VMR_1=VMR1,
                VMR_3=VMR3,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_3=logVMR3,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D4")

#Pivot longer for later merging and plotting
D4_class_count_vmr_long=D4_class_count_vmrs%>%
  pivot_longer(
    cols=-c(Class,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D4_class_count_vmr_long=as.data.frame(D4_class_count_vmr_long)
D4_class_count_vmr_long$Mat=as.factor(D4_class_count_vmr_long$Mat)
D4_class_count_vmr_long$Time=as.factor(D4_class_count_vmr_long$Time)
str(D4_class_count_vmr_long)

#########################
##All PHYLA
#########################
conflict_prefer("lmer", "lme4")
#Now combine all 4 mats into one dataframe
Phage_host_count_merge=rbind(D1_count_vmr_long,D2_count_vmr_long,D3_count_vmr_long,D4_count_vmr_long)
Phage_host_count_merge=Phage_host_count_merge%>%
  mutate_at("Phylum",str_replace,"_F","")
str(Phage_host_count_merge)
Phage_host_count_merge$Phylum=as.factor(Phage_host_count_merge$Phylum)
#Make log abundances a column for better liklihood calculation
Phage_host_count_merge=Phage_host_count_merge%>%
  mutate(logP=log10(P))%>%
  mutate(logH=log10(H))

#Plot overall
overall_p_h=ggplot(Phage_host_count_merge, aes(y=logP,x=logH))+
  geom_smooth(method=lm,se=T,color="black",
              aes(weight=wt))+
  geom_point(alpha=0.8, aes(color=Phylum))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral Abundance ("~log[10]~TPM~')'),
       x=expression("Host Abundance ("~log[10]~TPM~')'))+
  scale_color_viridis_d(direction=-1)+
  theme(legend.position = "none")
overall_p_h
ggsave("/Users/cissell/Desktop/overall_wls.svg",width=3.25,height=3.25,units="in",overall_p_h)

#Plot abundances per phyla
ggplot(Phage_host_count_merge, aes(y=log10(P),x=log10(H)))+
  geom_smooth(method=lm,se=F,
              aes(color=Phylum))+
  geom_point(alpha=0.8,
             aes(color=Phylum))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral Abundance ("~log[10]~TPM~')'),
       x=expression("Host Abundance ("~log[10]~TPM~')'))+
  scale_color_viridis_d(direction=-1)

#Plot VMR across time
ggplot(Phage_host_count_merge, aes(y=logVMR,x=Time))+
  geom_boxplot(aes(fill=Time),
               outlier.shape = NA,
               size=1,
               alpha=0.9)+
  scale_fill_manual(values=c("#90a295","#455765","#cd9fb2"),name="Time")+
  geom_point(aes(fill=Time),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMR ("~log[10]~')'),
       x="Sampling Time Point")
dev.print(png,fil="vmr_time_box.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

#Plot VMR across Mat
ggplot(Phage_host_count_merge, aes(y=logVMR,x=Mat))+
  geom_boxplot(aes(fill=Mat),
               outlier.shape = NA,
               size=1,
               alpha=0.9)+
  #  scale_fill_manual(values=c("#90a295","#455765","#cd9fb2"),name="Time")+
  geom_point(aes(fill=Time),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMR ("~log[10]~')'),
       x="Sampling Time Point")
dev.print(png,fil="vmr_time_box.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)
#plot logVMR across Phyla overall with pairwise significance
vmr_box=ggplot(Phage_host_count_merge, aes(y=logVMR,x=fct_rev(Phylum)))+
  geom_hline(yintercept=0,
             linetype="dashed",
             col="black",
             alpha=0.6)+
  geom_boxplot(aes(fill=Phylum),
               outlier.shape = NA,
               size=1,
               alpha=0.9)+
  geom_point(aes(fill=Phylum),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMR ("~log[10]~')'),
       x="Phylum")+
  scale_fill_viridis_d(direction=-1)+
  theme(legend.position="none")+
  coord_flip()
vmr_box
ggsave("/Users/cissell/Desktop/vmr_box.svg",width=6.5,height=4.5,units="in",vmr_box)

##########
#ALL CLASS
##########
#Combine all 4 mats into one dataframe
Phage_host_count_merge_class=rbind(D1_class_count_vmr_long,D2_class_count_vmr_long,D3_class_count_vmr_long,D4_class_count_vmr_long)
str(Phage_host_count_merge_class)
Phage_host_count_merge_class$Class=as.factor(Phage_host_count_merge_class$Class)
#Make log abundances a column for better liklihood calculation
Phage_host_count_merge_class=Phage_host_count_merge_class%>%
  mutate(logP=log10(P))%>%
  mutate(logH=log10(H))

#Spearman correlation
cor.test(Phage_host_count_merge_class$logH,Phage_host_count_merge_class$logP, method="spearman")
#S = 100814, p-value < 2.2e-16
#rho 0.8406632

#Plot overall
overall_p_h_class=ggplot(Phage_host_count_merge_class, aes(y=logP,x=logH))+
  geom_smooth(method=lm,se=T,color="black",
              aes(weight=wt_class))+
  geom_point(alpha=0.8, aes(color=Class))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral Abundance ("~log[10]~TPM~')'),
       x=expression("Host Abundance ("~log[10]~TPM~')'))+
  scale_color_viridis_d(direction=-1)+
  theme(legend.position = "none")
overall_p_h_class
ggsave("/Users/cissell/Desktop/overall_wls_class.svg",width=6.5,height=3.25,units="in",overall_p_h_class)

#Plot abundances per phyla
class_per=ggplot(Phage_host_count_merge_class, aes(y=log10(P),x=log10(H)))+
  geom_smooth(method=lm,se=F,
              aes(color=Class))+
  geom_point(alpha=0.8,
             aes(color=Class))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral Abundance ("~log[10]~TPM~')'),
       x=expression("Host Abundance ("~log[10]~TPM~')'))+
  scale_color_viridis_d(direction=-1) +
  theme(legend.position = "none")
class_per
ggsave("/Users/cissell/Desktop/overall_wls_class.svg",width=6.5,height=3.25,units="in",overall_p_h_class)
#Model logVMR across Phyla

clm3_class=lmer(data=Phage_host_count_merge_class,Phage_host_count_merge_class$logVMR~Phage_host_count_merge_class$Class+(1|Phage_host_count_merge_class$Mat))
summary(clm3_class)
#Test assumptions
library("DHARMa")
check_vmr_count_model3_class <- simulateResiduals(fittedModel = clm3_class, n = 999)
plot(check_vmr_count_model3_class)
#Check significance
Anova(clm3_class)

#plot logVMR across Class overall
library(forcats)
vmr_box_class=ggplot(Phage_host_count_merge_class, aes(y=logVMR,x=fct_rev(Class)))+
  geom_hline(yintercept=0,
             linetype="dashed",
             col="black",
             alpha=0.6)+
  geom_boxplot(aes(fill=Class),
               outlier.shape = NA,
               size=1,
               alpha=0.9)+
  geom_point(aes(fill=Class),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMR ("~log[10]~')'),
       x="Class")+
  scale_fill_viridis_d(direction=-1)+
  theme(legend.position="none")+
  coord_flip()
vmr_box_class
ggsave("/Users/cissell/Desktop/vmr_box_class.svg",width=6.5,height=5.5,units="in",vmr_box_class)

#####################################
#####Taxonomy of infection Chord Diagram
#####################################
#First subset tax tables down to mergetax and phylum
D1_tax_only_sub=D1_tax_only%>%
  select(-c(Host,Phage,Class))

D2_tax_only_sub=D2_tax_only%>%
  select(-c(Host,Phage,Class))

D3_tax_only_sub=D3_tax_only%>%
  select(-c(Host,Phage,Class))

D4_tax_only_sub=D4_tax_only%>%
  select(-c(Host,Phage,Class))

mergetax_chord_df=rbind(D1_tax_only_sub,D2_tax_only_sub,D3_tax_only_sub,D4_tax_only_sub)
str(mergetax_chord_df)
#Get stats on infectivity
chord_stats = mergetax_chord_unk_rem %>%
  group_by(Phylum) %>%
  summarise(count = n_distinct(mergetax))


mergetax_chord_unk_rem=mergetax_chord_df%>%
  subset(!mergetax=="Unknown")
str(mergetax_chord_unk_rem)

chord_net=graph_from_data_frame(mergetax_chord_unk_rem,directed=T)
chord_adjlist=as.matrix(get.adjacency(chord_net))
chord_adjlist=as.data.frame(chord_adjlist)
# I need a long format
chord_adjlist_long=chord_adjlist %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  mutate_at("rowname",str_replace,"_F","")%>%
  mutate_at("key",str_replace,"_F","")
chord_adjlist_long_clean=chord_adjlist_long%>%
  subset(!rowname=="Bacteroidota")%>%
  subset(!rowname=="Cyanobacteria")%>%
  subset(!rowname=="Desulfobacterota")%>%
  subset(!rowname=="Gemmatimonadota")%>%
  subset(!rowname=="J088")%>%
  subset(!rowname=="Patescibacteria")%>%
  subset(!rowname=="Planctomycetota")%>%
  subset(!rowname=="Proteobacteria")%>%
  subset(!rowname=="Verrucomicrobiota")%>%
  subset(!key=="Ackermannviridae")%>%
  subset(!key=="Autographiviridae")%>%
  subset(!key=="Demerecviridae")%>%
  subset(!key=="Drexlerviridae")%>%
  subset(!key=="Herelleviridae")%>%
  subset(!key=="Microviridae")%>%
  subset(!key=="Myoviridae")%>%
  subset(!key=="Podoviridae")%>%
  subset(!key=="Siphoviridae")
#Now begin building chord diagram to get an overall sense of where to split to build the master split plot
library(circlize)
library(viridis)
svg("/Users/cissell/Desktop/marked_chord.svg")
circos.clear()
circos.par(start.degree = 270, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))
groupColors=c(Ackermannviridae="#0c0304",
              Autographiviridae="#a6b1ec",
              Demerecviridae="#aee3d2",
              Drexlerviridae="#357aa2",
              Herelleviridae="#413f81",
              Microviridae="#3db4ad",
              Myoviridae="#382954",
              Podoviridae="#fa82a7",
              Siphoviridae="#3497a8",
              Bacteroidota="#a0da39",
              Cyanobacteria="#4ac16d",
              Desulfobacterota="#2db27d",
              Gemmatimonadota="#21918c",
              J088="#277f8e",
              Patescibacteria="#365c8d",
              Planctomycetota="#3f4788",
              Proteobacteria="#46327e",
              Verrucomicrobiota="#440154")
chordDiagram(chord_adjlist_long_clean,
             big.gap=50,
             transparency = 0.25,
             grid.col = groupColors,
             directional = 1,
             direction.type = c("arrows", "diffHeight"),
             diffHeight  = -0.04,
             annotationTrack = "grid",
             annotationTrackHeight = c(0.05, 0.1),
             link.arr.type = "big.arrow",
             link.sort = TRUE,
             link.largest.ontop = F,
             order=c("Ackermannviridae",
                     "Autographiviridae",
                     "Demerecviridae",
                     "Drexlerviridae",
                     "Herelleviridae",
                     "Microviridae",
                     "Myoviridae",
                     "Podoviridae",
                     "Siphoviridae",
                     "Bacteroidota",
                     "Cyanobacteria",
                     "Desulfobacterota",
                     "Gemmatimonadota",
                     "J088",
                     "Patescibacteria",
                     "Planctomycetota",
                     "Proteobacteria",
                     "Verrucomicrobiota"))
# Add text and axis initially to add as scale, then plot again and remove
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {

    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")

    # Add names to the sector.
    circos.text(
      x = mean(xlim),
      y = 3.2,
      labels = sector.index,
      facing = "bending",
      cex = 0.8
    )

    # Add graduation on axis
    circos.axis(
      h = "top",
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>200, yes = 50, no = 25)),
      major.tick.percentage = 0.5,
      labels.niceFacing = T)
  }
)
dev.off()
circos.clear()

#Now plot without gradation and labels
svg("/Users/cissell/Desktop/unmarked_chord.svg")
circos.clear()
circos.par(start.degree = 270, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))
groupColors=c(Ackermannviridae="#0c0304",
              Autographiviridae="#a6b1ec",
              Demerecviridae="#aee3d2",
              Drexlerviridae="#357aa2",
              Herelleviridae="#413f81",
              Microviridae="#3db4ad",
              Myoviridae="#382954",
              Podoviridae="#fa82a7",
              Siphoviridae="#3497a8",
              Bacteroidota="#a0da39",
              Cyanobacteria="#4ac16d",
              Desulfobacterota="#2db27d",
              Gemmatimonadota="#21918c",
              J088="#277f8e",
              Patescibacteria="#365c8d",
              Planctomycetota="#3f4788",
              Proteobacteria="#46327e",
              Verrucomicrobiota="#440154")
chordDiagram(chord_adjlist_long_clean,
             big.gap=50,
             transparency = 0.25,
             grid.col = groupColors,
             directional = 1,
             direction.type = c("arrows", "diffHeight"),
             diffHeight  = -0.04,
             annotationTrack = "grid",
             annotationTrackHeight = c(0.05, 0.1),
             link.arr.type = "big.arrow",
             link.sort = TRUE,
             link.largest.ontop = F,
             order=c("Ackermannviridae",
                     "Autographiviridae",
                     "Demerecviridae",
                     "Drexlerviridae",
                     "Herelleviridae",
                     "Microviridae",
                     "Myoviridae",
                     "Podoviridae",
                     "Siphoviridae",
                     "Bacteroidota",
                     "Cyanobacteria",
                     "Desulfobacterota",
                     "Gemmatimonadota",
                     "J088",
                     "Patescibacteria",
                     "Planctomycetota",
                     "Proteobacteria",
                     "Verrucomicrobiota"))
dev.off()



###########
#Diel Phage Activity
###########
library(dplyr)
library(gdata)
library(tidyverse)
library(conflicted)
library(car)
library(ggplot2)
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")
library(lme4)

##Read in counts
rna_phagediel1=read.table('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D1_RNA_phages_counts.tsv',header=T)
rna_phagediel2=read.table('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D2_RNA_phages_counts.tsv',header=T)
rna_phagediel3=read.table('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D3_RNA_phages_counts.tsv',header=T)
rna_phagediel4=read.table('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D4_RNA_phages_counts.tsv',header=T)

#Read in main Annotations
rna_phagekegg1=read.xls('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D1_main_annotations.xlsx')
rna_phagekegg2=read.xls('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D2_main_annotations.xlsx')
rna_phagekegg3=read.xls('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D3_main_annotations.xlsx')
rna_phagekegg4=read.xls('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D4_main_annotations.xlsx')
#Read in leftover annotations
rna_phagekegg11=read.xls('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D1_supp_annotations.xlsx')
rna_phagekegg22=read.xls('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D2_supp_annotations.xlsx')
rna_phagekegg33=read.xls('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D3_supp_annotations.xlsx')
rna_phagekegg44=read.xls('/Users/cissell/Desktop/VIRSORTER2_Diel/DRAM_AND_FEATURES/D4_supp_annotations.xlsx')

#Merge main and leftover
rna_phagekegg_merge1=rbind(rna_phagekegg1,rna_phagekegg11)
rna_phagekegg_merge2=rbind(rna_phagekegg2,rna_phagekegg22)
rna_phagekegg_merge3=rbind(rna_phagekegg3,rna_phagekegg33)
rna_phagekegg_merge4=rbind(rna_phagekegg4,rna_phagekegg44)

#SUBSET BY CURATED LIST
#Make blanks NA
rna_phagekegg_merge1=rna_phagekegg_merge1%>%
  mutate(across(everything(),~ifelse(.=="",NA,as.character(.))))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)
rna_phagekegg_merge2=rna_phagekegg_merge2%>%
  mutate(across(everything(),~ifelse(.=="",NA,as.character(.))))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)
rna_phagekegg_merge3=rna_phagekegg_merge3%>%
  mutate(across(everything(),~ifelse(.=="",NA,as.character(.))))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)
rna_phagekegg_merge4=rna_phagekegg_merge4%>%
  mutate(across(everything(),~ifelse(.=="",NA,as.character(.))))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)
#Make characters not factors
rna_phagekegg_merge1$kegg_hit=as.character(rna_phagekegg_merge1$kegg_hit)
rna_phagekegg_merge1$viral_hit=as.character(rna_phagekegg_merge1$viral_hit)
rna_phagekegg_merge1$pfam_hits=as.character(rna_phagekegg_merge1$pfam_hits)
rna_phagekegg_merge1$vogdb_description=as.character(rna_phagekegg_merge1$vogdb_description)

rna_phagekegg_merge2$kegg_hit=as.character(rna_phagekegg_merge2$kegg_hit)
rna_phagekegg_merge2$viral_hit=as.character(rna_phagekegg_merge2$viral_hit)
rna_phagekegg_merge2$pfam_hits=as.character(rna_phagekegg_merge2$pfam_hits)
rna_phagekegg_merge2$vogdb_description=as.character(rna_phagekegg_merge2$vogdb_description)

rna_phagekegg_merge3$kegg_hit=as.character(rna_phagekegg_merge3$kegg_hit)
rna_phagekegg_merge3$viral_hit=as.character(rna_phagekegg_merge3$viral_hit)
rna_phagekegg_merge3$pfam_hits=as.character(rna_phagekegg_merge3$pfam_hits)
rna_phagekegg_merge3$vogdb_description=as.character(rna_phagekegg_merge3$vogdb_description)

rna_phagekegg_merge4$kegg_hit=as.character(rna_phagekegg_merge4$kegg_hit)
rna_phagekegg_merge4$viral_hit=as.character(rna_phagekegg_merge4$viral_hit)
rna_phagekegg_merge4$pfam_hits=as.character(rna_phagekegg_merge4$pfam_hits)
rna_phagekegg_merge4$vogdb_description=as.character(rna_phagekegg_merge4$vogdb_description)

#Create unified annotation column by series of rules across database hits
rna_phagekegg_merge11=rna_phagekegg_merge1%>%
  mutate(annotate=ifelse(is.na(viral_hit) & is.na(pfam_hits) & is.na(vogdb_description),kegg_hit[],ifelse(is.na(viral_hit) & is.na (pfam_hits),vogdb_description[],ifelse(is.na(viral_hit),pfam_hits[],viral_hit[]))))

rna_phagekegg_merge22=rna_phagekegg_merge2%>%
  mutate(annotate=ifelse(is.na(viral_hit) & is.na(pfam_hits) & is.na(vogdb_description),kegg_hit[],ifelse(is.na(viral_hit) & is.na (pfam_hits),vogdb_description[],ifelse(is.na(viral_hit),pfam_hits[],viral_hit[]))))

rna_phagekegg_merge33=rna_phagekegg_merge3%>%
  mutate(annotate=ifelse(is.na(viral_hit) & is.na(pfam_hits) & is.na(vogdb_description),kegg_hit[],ifelse(is.na(viral_hit) & is.na (pfam_hits),vogdb_description[],ifelse(is.na(viral_hit),pfam_hits[],viral_hit[]))))

rna_phagekegg_merge44=rna_phagekegg_merge4%>%
  mutate(annotate=ifelse(is.na(viral_hit) & is.na(pfam_hits) & is.na(vogdb_description),kegg_hit[],ifelse(is.na(viral_hit) & is.na (pfam_hits),vogdb_description[],ifelse(is.na(viral_hit),pfam_hits[],viral_hit[]))))

#Now merge in counts
rna_phagekegg_merge_counts1=merge(rna_phagekegg_merge11,
                                  rna_phagediel1,by='Geneid',
                                  all=T)
rna_phagekegg_merge_counts2=merge(rna_phagekegg_merge22,
                                  rna_phagediel2,by='Geneid',
                                  all=T)
rna_phagekegg_merge_counts3=merge(rna_phagekegg_merge33,
                                  rna_phagediel3,by='Geneid',
                                  all=T)
rna_phagekegg_merge_counts4=merge(rna_phagekegg_merge44,
                                  rna_phagediel4,by='Geneid',
                                  all=T)

#######################
#Now merge in Phage Taxonomy and Host prediction + taxonomy, then subset to only what we need and summarize within viral family
rna_phagekegg_merge_counts1_tax=rna_phagekegg_merge_counts1%>%
  mutate(contig_id=contig_id.y)%>%
  select(-c(contig_id.x,contig_id.y))%>%
  merge(viral_tax_simpcol,by="contig_id")%>%
  select(-c(contig_id,
            Geneid,
            gene_position,
            start_position,
            end_position,
            kegg_id,
            kegg_hit,
            viral_id,
            viral_hit,
            pfam_hits,
            cazy_hits,
            vogdb_description,
            vogdb_hits,
            annotate,
            Start,
            End,
            Length))%>%
  group_by(mergetax)%>%
  summarise(across(1:5,sum))


rna_phagekegg_merge_counts2_tax=rna_phagekegg_merge_counts2%>%
  mutate(contig_id=contig_id.y)%>%
  select(-c(contig_id.x,contig_id.y))%>%
  merge(viral_tax_simpcol,by="contig_id")%>%
  select(-c(contig_id,
            Geneid,
            gene_position,
            start_position,
            end_position,
            kegg_id,
            kegg_hit,
            viral_id,
            viral_hit,
            pfam_hits,
            cazy_hits,
            vogdb_description,
            vogdb_hits,
            annotate,
            Start,
            End,
            Length))%>%
  group_by(mergetax)%>%
  summarise(across(1:5,sum))

rna_phagekegg_merge_counts3_tax=rna_phagekegg_merge_counts3%>%
  mutate(contig_id=contig_id.y)%>%
  select(-c(contig_id.x,contig_id.y))%>%
  merge(viral_tax_simpcol,by="contig_id")%>%
  select(-c(contig_id,
            Geneid,
            gene_position,
            start_position,
            end_position,
            kegg_id,
            kegg_hit,
            viral_id,
            viral_hit,
            pfam_hits,
            cazy_hits,
            vogdb_description,
            vogdb_hits,
            annotate,
            Start,
            End,
            Length))%>%
  group_by(mergetax)%>%
  summarise(across(1:5,sum))

rna_phagekegg_merge_counts4_tax=rna_phagekegg_merge_counts4%>%
  mutate(contig_id=contig_id.y)%>%
  select(-c(contig_id.x,contig_id.y))%>%
  merge(viral_tax_simpcol,by="contig_id")%>%
  select(-c(contig_id,
            Geneid,
            gene_position,
            start_position,
            end_position,
            kegg_id,
            kegg_hit,
            viral_id,
            viral_hit,
            pfam_hits,
            cazy_hits,
            vogdb_description,
            vogdb_hits,
            annotate,
            Start,
            End,
            Length))%>%
  group_by(mergetax)%>%
  summarise(across(1:5,sum))

#########################################################
####Do VMAR calculations now
###################################
#First we need to get the length / count sum across ORFs per genome, filter to the curated list, and add in host predictions
rna_phage_sum_orf1=rna_phagediel1%>%
  select(-c(Geneid,Start,End))%>%
  group_by(contig_id)%>%
  summarise(across(1:6,sum))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)%>%
  mutate(Phage=contig_id)%>%
  select(-contig_id)%>%
  merge(D1_phage_host_final,by="Phage")%>%
  mutate(PLength=Length)%>%
  select(-Length)

rna_phage_sum_orf2=rna_phagediel2%>%
  select(-c(Geneid,Start,End))%>%
  group_by(contig_id)%>%
  summarise(across(1:6,sum))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)%>%
  mutate(Phage=contig_id)%>%
  select(-contig_id)%>%
  merge(D2_phage_host_final,by="Phage")%>%
  mutate(PLength=Length)%>%
  select(-Length)

rna_phage_sum_orf3=rna_phagediel3%>%
  select(-c(Geneid,Start,End))%>%
  group_by(contig_id)%>%
  summarise(across(1:6,sum))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)%>%
  mutate(Phage=contig_id)%>%
  select(-contig_id)%>%
  merge(D3_phage_host_final,by="Phage")%>%
  mutate(PLength=Length)%>%
  select(-Length)

rna_phage_sum_orf4=rna_phagediel4%>%
  select(-c(Geneid,Start,End))%>%
  group_by(contig_id)%>%
  summarise(across(1:6,sum))%>%
  ungroup()%>%
  filter(contig_id%in%viral_clust_curate_tax$contig_id)%>%
  mutate(Phage=contig_id)%>%
  select(-contig_id)%>%
  merge(D4_phage_host_final,by="Phage")%>%
  mutate(PLength=Length)%>%
  select(-Length)

#Read in my modified Host RNA counts per gene region for ORF summing
host_orf_counts1=read.xls("/Users/cissell/Desktop/D1_RNA_counts.xlsx")
host_orf_counts2=read.xls("/Users/cissell/Desktop/D2_RNA_counts.xlsx")
host_orf_counts3=read.xls("/Users/cissell/Desktop/D3_RNA_counts.xlsx")
host_orf_counts4=read.xls("/Users/cissell/Desktop/D4_RNA_counts.xlsx")

#Now sum lengths and counts to Host, plus add in host taxonomy
host_orf_sumcounts1=host_orf_counts1%>%
  select(-c(contig))%>%
  group_by(Host)%>%
  dplyr::summarise(across(1:6,sum))%>%
  ungroup()%>%
  merge(D1_host_tax,by="Host")%>%
  mutate(HLength=Length)%>%
  select(-Length)

host_orf_sumcounts2=host_orf_counts2%>%
  select(-c(contig))%>%
  group_by(Host)%>%
  dplyr::summarise(across(1:6,sum))%>%
  ungroup()%>%
  merge(D2_host_tax,by="Host")%>%
  mutate(HLength=Length)%>%
  select(-Length)

host_orf_sumcounts3=host_orf_counts3%>%
  select(-c(contig))%>%
  group_by(Host)%>%
  dplyr::summarise(across(1:6,sum))%>%
  ungroup()%>%
  merge(D3_host_tax,by="Host")%>%
  mutate(HLength=Length)%>%
  select(-Length)

host_orf_sumcounts4=host_orf_counts4%>%
  select(-c(contig))%>%
  group_by(Host)%>%
  dplyr::summarise(across(1:6,sum))%>%
  ungroup()%>%
  merge(D4_host_tax,by="Host")%>%
  mutate(HLength=Length)%>%
  select(-Length)

##Merge together phage and host counts
VMAR_host_phage_count_merge_1=merge(rna_phage_sum_orf1,
                                    host_orf_sumcounts1,
                                    by="Host")

VMAR_host_phage_count_merge_2=merge(rna_phage_sum_orf2,
                                    host_orf_sumcounts2,
                                    by="Host")

VMAR_host_phage_count_merge_3=merge(rna_phage_sum_orf3,
                                    host_orf_sumcounts3,
                                    by="Host")

VMAR_host_phage_count_merge_4=merge(rna_phage_sum_orf4,
                                    host_orf_sumcounts4,
                                    by="Host")

#Now calculate TPM in RNA for each distinct mat, then we will bring back together
#Make new df
rpk1mar=VMAR_host_phage_count_merge_1
str(rpk1mar)
#Put length in KB
rpk1mar$PLength=rpk1mar$PLength/1000
rpk1mar$HLength=rpk1mar$HLength/1000
#Create RPK for Mat 1
rpk1mar$P1_1=rpk1mar$P1_1/rpk1mar$PLength
rpk1mar$R1_1=rpk1mar$R1_1/rpk1mar$HLength
rpk1mar$P1_2=rpk1mar$P1_2/rpk1mar$PLength
rpk1mar$R1_2=rpk1mar$R1_2/rpk1mar$HLength
rpk1mar$P1_3=rpk1mar$P1_3/rpk1mar$PLength
rpk1mar$R1_3=rpk1mar$R1_3/rpk1mar$HLength
rpk1mar$P1_4=rpk1mar$P1_4/rpk1mar$PLength
rpk1mar$R1_4=rpk1mar$R1_4/rpk1mar$HLength
rpk1mar$P1_5=rpk1mar$P1_5/rpk1mar$PLength
rpk1mar$R1_5=rpk1mar$R1_5/rpk1mar$HLength
#Need to extract a new dataframe of only unique so each MAG only occurs once
hpk1mar=rpk1mar[!duplicated(rpk1mar$Host),]
#create scaling factors
#Mat 1
sum11mar=((sum(rpk1mar$P1_1)+sum(hpk1mar$R1_1))/1000000)
sum12mar=((sum(rpk1mar$P1_2)+sum(hpk1mar$R1_2))/1000000)
sum13mar=((sum(rpk1mar$P1_3)+sum(hpk1mar$R1_3))/1000000)
sum14mar=((sum(rpk1mar$P1_4)+sum(hpk1mar$R1_4))/1000000)
sum15mar=((sum(rpk1mar$P1_5)+sum(hpk1mar$R1_5))/1000000)
#Create TPM for Mat 1
rpk1mar$P1_1=rpk1mar$P1_1/sum11mar
rpk1mar$R1_1=rpk1mar$R1_1/sum11mar
rpk1mar$P1_2=rpk1mar$P1_2/sum12mar
rpk1mar$R1_2=rpk1mar$R1_2/sum12mar
rpk1mar$P1_3=rpk1mar$P1_3/sum13mar
rpk1mar$R1_3=rpk1mar$R1_3/sum13mar
rpk1mar$P1_4=rpk1mar$P1_4/sum14mar
rpk1mar$R1_4=rpk1mar$R1_4/sum14mar
rpk1mar$P1_5=rpk1mar$P1_5/sum15mar
rpk1mar$R1_5=rpk1mar$R1_5/sum15mar
#Wont make sense just yet becuase of duplicates

#Bring it back to a sensical name now
D1_count_tpmmar=rpk1mar

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
D1_host_lostmar=D1_count_tpmmar%>%
  select(c(Host,R1_1,R1_2,R1_3,R1_4,R1_5,Phylum))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D1_count_tpm2mar=D1_count_tpmmar%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P1_1","P1_2","P1_3","P1_4","P1_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D1_count_mergermar=merge(D1_count_tpm2mar,D1_host_lostmar,by="Host",all=T)
D1_count_mergermar$Phylum=as.factor(D1_count_mergermar$Phylum)
str(D1_count_mergermar)
#Summarize within phyla
D1_phyla_countsmar=D1_count_mergermar%>%
  select(-c(Host))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D1_phyla_countsmar$P1_1)+sum(D1_phyla_countsmar$R1_1)
sum(D1_phyla_countsmar$P1_2)+sum(D1_phyla_countsmar$R1_2)
sum(D1_phyla_countsmar$P1_3)+sum(D1_phyla_countsmar$R1_3)
sum(D1_phyla_countsmar$P1_4)+sum(D1_phyla_countsmar$R1_4)
sum(D1_phyla_countsmar$P1_5)+sum(D1_phyla_countsmar$R1_5)
#Nice

#Create each individual VMR and logVMR for each pair
D1_count_vmrsmar=D1_phyla_countsmar%>%
  mutate(VMR1=P1_1/R1_1)%>%
  mutate(logVMR1=log10(P1_1/R1_1))%>%
  mutate(VMR2=P1_2/R1_2)%>%
  mutate(logVMR2=log10(P1_2/R1_2))%>%
  mutate(VMR3=P1_3/R1_3)%>%
  mutate(logVMR3=log10(P1_3/R1_3))%>%
  mutate(VMR4=P1_4/R1_4)%>%
  mutate(logVMR4=log10(P1_4/R1_4))%>%
  mutate(VMR5=P1_5/R1_5)%>%
  mutate(logVMR5=log10(P1_5/R1_5))%>%
  dplyr::rename(P_1=P1_1,
                P_2=P1_2,
                P_3=P1_3,
                P_4=P1_4,
                P_5=P1_5,
                R_1=R1_1,
                R_2=R1_2,
                R_3=R1_3,
                R_4=R1_4,
                R_5=R1_5,
                VMR_1=VMR1,
                VMR_2=VMR2,
                VMR_3=VMR3,
                VMR_4=VMR4,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_2=logVMR2,
                logVMR_3=logVMR3,
                logVMR_4=logVMR4,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D1")

#Pivot longer for later merging and plotting
D1_count_vmr_longmar=D1_count_vmrsmar%>%
  pivot_longer(
    cols=-c(Phylum,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D1_count_vmr_longmar=as.data.frame(D1_count_vmr_longmar)
D1_count_vmr_longmar$Mat=as.factor(D1_count_vmr_longmar$Mat)
D1_count_vmr_longmar$Time=as.factor(D1_count_vmr_longmar$Time)
str(D1_count_vmr_longmar)

###############
#D2
###############
#Make new df
rpk2mar=VMAR_host_phage_count_merge_2
str(rpk2mar)
#Put length in KB
rpk2mar$PLength=rpk2mar$PLength/1000
rpk2mar$HLength=rpk2mar$HLength/1000
#Create RPK for Mat 2
rpk2mar$P2_1=rpk2mar$P2_1/rpk2mar$PLength
rpk2mar$R2_1=rpk2mar$R2_1/rpk2mar$HLength
rpk2mar$P2_2=rpk2mar$P2_2/rpk2mar$PLength
rpk2mar$R2_2=rpk2mar$R2_2/rpk2mar$HLength
rpk2mar$P2_3=rpk2mar$P2_3/rpk2mar$PLength
rpk2mar$R2_3=rpk2mar$R2_3/rpk2mar$HLength
rpk2mar$P2_4=rpk2mar$P2_4/rpk2mar$PLength
rpk2mar$R2_4=rpk2mar$R2_4/rpk2mar$HLength
rpk2mar$P2_5=rpk2mar$P2_5/rpk2mar$PLength
rpk2mar$R2_5=rpk2mar$R2_5/rpk2mar$HLength
#Need to extract a new dataframe of only unique so each MAG only occurs once
hpk2mar=rpk2mar[!duplicated(rpk2mar$Host),]
#create scaling factors
#Mat 2
sum21mar=((sum(rpk2mar$P2_1)+sum(hpk2mar$R2_1))/1000000)
sum22mar=((sum(rpk2mar$P2_2)+sum(hpk2mar$R2_2))/1000000)
sum23mar=((sum(rpk2mar$P2_3)+sum(hpk2mar$R2_3))/1000000)
sum24mar=((sum(rpk2mar$P2_4)+sum(hpk2mar$R2_4))/1000000)
sum25mar=((sum(rpk2mar$P2_5)+sum(hpk2mar$R2_5))/1000000)
#Create TPM for Mat 2
rpk2mar$P2_1=rpk2mar$P2_1/sum21mar
rpk2mar$R2_1=rpk2mar$R2_1/sum21mar
rpk2mar$P2_2=rpk2mar$P2_2/sum22mar
rpk2mar$R2_2=rpk2mar$R2_2/sum22mar
rpk2mar$P2_3=rpk2mar$P2_3/sum23mar
rpk2mar$R2_3=rpk2mar$R2_3/sum23mar
rpk2mar$P2_4=rpk2mar$P2_4/sum24mar
rpk2mar$R2_4=rpk2mar$R2_4/sum24mar
rpk2mar$P2_5=rpk2mar$P2_5/sum25mar
rpk2mar$R2_5=rpk2mar$R2_5/sum25mar
#Wont make sense just yet becuase of duplicates

#Bring it back to a sensical name now
D2_count_tpmmar=rpk2mar

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
D2_host_lostmar=D2_count_tpmmar%>%
  select(c(Host,R2_1,R2_2,R2_3,R2_4,R2_5,Phylum))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D2_count_tpm2mar=D2_count_tpmmar%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P2_1","P2_2","P2_3","P2_4","P2_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D2_count_mergermar=merge(D2_count_tpm2mar,D2_host_lostmar,by="Host",all=T)
D2_count_mergermar$Phylum=as.factor(D2_count_mergermar$Phylum)
str(D2_count_mergermar)
#Summarize within phyla
D2_phyla_countsmar=D2_count_mergermar%>%
  select(-c(Host))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D2_phyla_countsmar$P2_1)+sum(D2_phyla_countsmar$R2_1)
sum(D2_phyla_countsmar$P2_2)+sum(D2_phyla_countsmar$R2_2)
sum(D2_phyla_countsmar$P2_3)+sum(D2_phyla_countsmar$R2_3)
sum(D2_phyla_countsmar$P2_4)+sum(D2_phyla_countsmar$R2_4)
sum(D2_phyla_countsmar$P2_5)+sum(D2_phyla_countsmar$R2_5)
#Nice

#Create each individual VMR and logVMR for each pair
D2_count_vmrsmar=D2_phyla_countsmar%>%
  mutate(VMR1=P2_1/R2_1)%>%
  mutate(logVMR1=log10(P2_1/R2_1))%>%
  mutate(VMR2=P2_2/R2_2)%>%
  mutate(logVMR2=log10(P2_2/R2_2))%>%
  mutate(VMR3=P2_3/R2_3)%>%
  mutate(logVMR3=log10(P2_3/R2_3))%>%
  mutate(VMR4=P2_4/R2_4)%>%
  mutate(logVMR4=log10(P2_4/R2_4))%>%
  mutate(VMR5=P2_5/R2_5)%>%
  mutate(logVMR5=log10(P2_5/R2_5))%>%
  dplyr::rename(P_1=P2_1,
                P_2=P2_2,
                P_3=P2_3,
                P_4=P2_4,
                P_5=P2_5,
                R_1=R2_1,
                R_2=R2_2,
                R_3=R2_3,
                R_4=R2_4,
                R_5=R2_5,
                VMR_1=VMR1,
                VMR_2=VMR2,
                VMR_3=VMR3,
                VMR_4=VMR4,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_2=logVMR2,
                logVMR_3=logVMR3,
                logVMR_4=logVMR4,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D2")

#Pivot longer for later merging and plotting
D2_count_vmr_longmar=D2_count_vmrsmar%>%
  pivot_longer(
    cols=-c(Phylum,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D2_count_vmr_longmar=as.data.frame(D2_count_vmr_longmar)
D2_count_vmr_longmar$Mat=as.factor(D2_count_vmr_longmar$Mat)
D2_count_vmr_longmar$Time=as.factor(D2_count_vmr_longmar$Time)
str(D2_count_vmr_longmar)

###############
#D3
###############
#Make new df
rpk3mar=VMAR_host_phage_count_merge_3
str(rpk3mar)
#Put length in KB
rpk3mar$PLength=rpk3mar$PLength/1000
rpk3mar$HLength=rpk3mar$HLength/1000
#Create RPK for Mat 2
rpk3mar$P3_1=rpk3mar$P3_1/rpk3mar$PLength
rpk3mar$R3_1=rpk3mar$R3_1/rpk3mar$HLength
rpk3mar$P3_2=rpk3mar$P3_2/rpk3mar$PLength
rpk3mar$R3_2=rpk3mar$R3_2/rpk3mar$HLength
rpk3mar$P3_3=rpk3mar$P3_3/rpk3mar$PLength
rpk3mar$R3_3=rpk3mar$R3_3/rpk3mar$HLength
rpk3mar$P3_4=rpk3mar$P3_4/rpk3mar$PLength
rpk3mar$R3_4=rpk3mar$R3_4/rpk3mar$HLength
rpk3mar$P3_5=rpk3mar$P3_5/rpk3mar$PLength
rpk3mar$R3_5=rpk3mar$R3_5/rpk3mar$HLength
#Need to extract a new dataframe of only unique so each MAG only occurs once
hpk3mar=rpk3mar[!duplicated(rpk3mar$Host),]
#create scaling factors
#Mat 3
sum31mar=((sum(rpk3mar$P3_1)+sum(hpk3mar$R3_1))/1000000)
sum32mar=((sum(rpk3mar$P3_2)+sum(hpk3mar$R3_2))/1000000)
sum33mar=((sum(rpk3mar$P3_3)+sum(hpk3mar$R3_3))/1000000)
sum34mar=((sum(rpk3mar$P3_4)+sum(hpk3mar$R3_4))/1000000)
sum35mar=((sum(rpk3mar$P3_5)+sum(hpk3mar$R3_5))/1000000)
#Create TPM for Mat 3
rpk3mar$P3_1=rpk3mar$P3_1/sum31mar
rpk3mar$R3_1=rpk3mar$R3_1/sum31mar
rpk3mar$P3_2=rpk3mar$P3_2/sum32mar
rpk3mar$R3_2=rpk3mar$R3_2/sum32mar
rpk3mar$P3_3=rpk3mar$P3_3/sum33mar
rpk3mar$R3_3=rpk3mar$R3_3/sum33mar
rpk3mar$P3_4=rpk3mar$P3_4/sum34mar
rpk3mar$R3_4=rpk3mar$R3_4/sum34mar
rpk3mar$P3_5=rpk3mar$P3_5/sum35mar
rpk3mar$R3_5=rpk3mar$R3_5/sum35mar
#Wont make sense just yet becuase of duplicates

#Bring it back to a sensical name now
D3_count_tpmmar=rpk3mar

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
D3_host_lostmar=D3_count_tpmmar%>%
  select(c(Host,R3_1,R3_2,R3_3,R3_4,R3_5,Phylum))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D3_count_tpm3mar=D3_count_tpmmar%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P3_1","P3_2","P3_3","P3_4","P3_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D3_count_mergermar=merge(D3_count_tpm3mar,D3_host_lostmar,by="Host",all=T)
D3_count_mergermar$Phylum=as.factor(D3_count_mergermar$Phylum)
D3_count_mergermar=D3_count_mergermar%>%
  mutate_at("Phylum",str_replace,"_F","")
str(D3_count_mergermar)
#Summarize within phyla
D3_phyla_countsmar=D3_count_mergermar%>%
  select(-c(Host))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D3_phyla_countsmar$P3_1)+sum(D3_phyla_countsmar$R3_1)
sum(D3_phyla_countsmar$P3_2)+sum(D3_phyla_countsmar$R3_2)
sum(D3_phyla_countsmar$P3_3)+sum(D3_phyla_countsmar$R3_3)
sum(D3_phyla_countsmar$P3_4)+sum(D3_phyla_countsmar$R3_4)
sum(D3_phyla_countsmar$P3_5)+sum(D3_phyla_countsmar$R3_5)
#Nice

#Create each individual VMR and logVMR for each pair
D3_count_vmrsmar=D3_phyla_countsmar%>%
  mutate(VMR1=P3_1/R3_1)%>%
  mutate(logVMR1=log10(P3_1/R3_1))%>%
  mutate(VMR2=P3_2/R3_2)%>%
  mutate(logVMR2=log10(P3_2/R3_2))%>%
  mutate(VMR3=P3_3/R3_3)%>%
  mutate(logVMR3=log10(P3_3/R3_3))%>%
  mutate(VMR4=P3_4/R3_4)%>%
  mutate(logVMR4=log10(P3_4/R3_4))%>%
  mutate(VMR5=P3_5/R3_5)%>%
  mutate(logVMR5=log10(P3_5/R3_5))%>%
  dplyr::rename(P_1=P3_1,
                P_2=P3_2,
                P_3=P3_3,
                P_4=P3_4,
                P_5=P3_5,
                R_1=R3_1,
                R_2=R3_2,
                R_3=R3_3,
                R_4=R3_4,
                R_5=R3_5,
                VMR_1=VMR1,
                VMR_2=VMR2,
                VMR_3=VMR3,
                VMR_4=VMR4,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_2=logVMR2,
                logVMR_3=logVMR3,
                logVMR_4=logVMR4,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D3")

#Pivot longer for later merging and plotting
D3_count_vmr_longmar=D3_count_vmrsmar%>%
  pivot_longer(
    cols=-c(Phylum,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D3_count_vmr_longmar=as.data.frame(D3_count_vmr_longmar)
D3_count_vmr_longmar$Mat=as.factor(D3_count_vmr_longmar$Mat)
D3_count_vmr_longmar$Time=as.factor(D3_count_vmr_longmar$Time)
str(D3_count_vmr_longmar)

###############
#D4
###############
#Make new df
rpk4mar=VMAR_host_phage_count_merge_4
str(rpk4mar)
#Put length in KB
rpk4mar$PLength=rpk4mar$PLength/1000
rpk4mar$HLength=rpk4mar$HLength/1000
#Create RPK for Mat 2
rpk4mar$P4_1=rpk4mar$P4_1/rpk4mar$PLength
rpk4mar$R4_1=rpk4mar$R4_1/rpk4mar$HLength
rpk4mar$P4_2=rpk4mar$P4_2/rpk4mar$PLength
rpk4mar$R4_2=rpk4mar$R4_2/rpk4mar$HLength
rpk4mar$P4_3=rpk4mar$P4_3/rpk4mar$PLength
rpk4mar$R4_3=rpk4mar$R4_3/rpk4mar$HLength
rpk4mar$P4_4=rpk4mar$P4_4/rpk4mar$PLength
rpk4mar$R4_4=rpk4mar$R4_4/rpk4mar$HLength
rpk4mar$P4_5=rpk4mar$P4_5/rpk4mar$PLength
rpk4mar$R4_5=rpk4mar$R4_5/rpk4mar$HLength
#Need to extract a new dataframe of only unique so each MAG only occurs once
hpk4mar=rpk4mar[!duplicated(rpk4mar$Host),]
#create scaling factors
#Mat 4
sum41mar=((sum(rpk4mar$P4_1)+sum(hpk4mar$R4_1))/1000000)
sum42mar=((sum(rpk4mar$P4_2)+sum(hpk4mar$R4_2))/1000000)
sum43mar=((sum(rpk4mar$P4_3)+sum(hpk4mar$R4_3))/1000000)
sum44mar=((sum(rpk4mar$P4_4)+sum(hpk4mar$R4_4))/1000000)
sum45mar=((sum(rpk4mar$P4_5)+sum(hpk4mar$R4_5))/1000000)
#Create TPM for Mat 4
rpk4mar$P4_1=rpk4mar$P4_1/sum41mar
rpk4mar$R4_1=rpk4mar$R4_1/sum41mar
rpk4mar$P4_2=rpk4mar$P4_2/sum42mar
rpk4mar$R4_2=rpk4mar$R4_2/sum42mar
rpk4mar$P4_3=rpk4mar$P4_3/sum43mar
rpk4mar$R4_3=rpk4mar$R4_3/sum43mar
rpk4mar$P4_4=rpk4mar$P4_4/sum44mar
rpk4mar$R4_4=rpk4mar$R4_4/sum44mar
rpk4mar$P4_5=rpk4mar$P4_5/sum45mar
rpk4mar$R4_5=rpk4mar$R4_5/sum45mar
#Wont make sense just yet becuase of duplicates

#Bring it back to a sensical name now
D4_count_tpmmar=rpk4mar

#Summarize within host
#Summarize within Host
conflict_prefer("select", "dplyr")
#Create dataframe to merge back after summarizing
D4_host_lostmar=D4_count_tpmmar%>%
  select(c(Host,R4_1,R4_2,R4_3,R4_4,R4_5,Phylum))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D4_count_tpm4mar=D4_count_tpmmar%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P4_1","P4_2","P4_3","P4_4","P4_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D4_count_mergermar=merge(D4_count_tpm4mar,D4_host_lostmar,by="Host",all=T)
D4_count_mergermar$Phylum=as.factor(D4_count_mergermar$Phylum)
str(D4_count_mergermar)
#Summarize within phyla
D4_phyla_countsmar=D4_count_mergermar%>%
  select(-c(Host))%>%
  group_by(Phylum)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D4_phyla_countsmar$P4_1)+sum(D4_phyla_countsmar$R4_1)
sum(D4_phyla_countsmar$P4_2)+sum(D4_phyla_countsmar$R4_2)
sum(D4_phyla_countsmar$P4_3)+sum(D4_phyla_countsmar$R4_3)
sum(D4_phyla_countsmar$P4_4)+sum(D4_phyla_countsmar$R4_4)
sum(D4_phyla_countsmar$P4_5)+sum(D4_phyla_countsmar$R4_5)
#Nice

#Create each individual VMR and logVMR for each pair
D4_count_vmrsmar=D4_phyla_countsmar%>%
  mutate(VMR1=P4_1/R4_1)%>%
  mutate(logVMR1=log10(P4_1/R4_1))%>%
  mutate(VMR2=P4_2/R4_2)%>%
  mutate(logVMR2=log10(P4_2/R4_2))%>%
  mutate(VMR3=P4_3/R4_3)%>%
  mutate(logVMR3=log10(P4_3/R4_3))%>%
  mutate(VMR4=P4_4/R4_4)%>%
  mutate(logVMR4=log10(P4_4/R4_4))%>%
  mutate(VMR5=P4_5/R4_5)%>%
  mutate(logVMR5=log10(P4_5/R4_5))%>%
  dplyr::rename(P_1=P4_1,
                P_2=P4_2,
                P_3=P4_3,
                P_4=P4_4,
                P_5=P4_5,
                R_1=R4_1,
                R_2=R4_2,
                R_3=R4_3,
                R_4=R4_4,
                R_5=R4_5,
                VMR_1=VMR1,
                VMR_2=VMR2,
                VMR_3=VMR3,
                VMR_4=VMR4,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_2=logVMR2,
                logVMR_3=logVMR3,
                logVMR_4=logVMR4,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D4")

#Pivot longer for later merging and plotting
D4_count_vmr_longmar=D4_count_vmrsmar%>%
  pivot_longer(
    cols=-c(Phylum,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D4_count_vmr_longmar=as.data.frame(D4_count_vmr_longmar)
D4_count_vmr_longmar$Mat=as.factor(D4_count_vmr_longmar$Mat)
D4_count_vmr_longmar$Time=as.factor(D4_count_vmr_longmar$Time)
str(D4_count_vmr_longmar)

#########################
##All
#########################
#Now combine all 4 mats into one dataframe
Phage_host_count_mergemar=rbind(D1_count_vmr_longmar,D2_count_vmr_longmar,D3_count_vmr_longmar,D4_count_vmr_longmar)
Phage_host_count_mergemar=Phage_host_count_mergemar%>%
  mutate_at("Phylum",str_replace,"_F","")
str(Phage_host_count_mergemar)
Phage_host_count_mergemar$Phylum=as.factor(Phage_host_count_mergemar$Phylum)
#Make log abundances a column for better liklihood calculation
Phage_host_count_mergemar=Phage_host_count_mergemar%>%
  mutate(logP=log10(P))%>%
  mutate(logH=log10(R))

#Plot overall
overall_p_hmar=ggplot(Phage_host_count_mergemar, aes(y=logP,x=logH))+
  geom_smooth(method=lm,se=T,color="black")+
  geom_point(alpha=0.8, aes(color=Phylum))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral RNA Abundance ("~log[10]~TPM~')'),
       x=expression("Host RNA Abundance ("~log[10]~TPM~')'))+
  scale_color_viridis_d(direction=-1)+
  theme(legend.position = "none")
overall_p_hmar
ggsave("/Users/cissell/Desktop/overall_wlsmar.svg",width=3.25,height=3.25,units="in",overall_p_hmar)

#Plot abundances per phyla
ggplot(Phage_host_count_mergemar, aes(y=log10(P),x=log10(R)))+
  geom_smooth(method=lm,se=F,
              aes(color=Phylum))+
  geom_point(alpha=0.8,
             aes(color=Phylum))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral Abundance ("~log[10]~TPM~')'),
       x=expression("Host Abundance ("~log[10]~TPM~')'))+
  scale_color_viridis_d(direction=-1)

#Plot VMAR across time
ggplot(Phage_host_count_mergemar, aes(y=logVMR,x=Time))+
  geom_boxplot(aes(fill=Phylum),
               outlier.shape = NA,
               size=1,
               alpha=0.7)+
  #  scale_fill_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2"),name="Time") +
  geom_point(aes(color=Phylum),position=position_jitterdodge(0.2),pch=19,
             alpha=0.9,
             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMAR ("~log[10]~')'),
       x="Sampling Time Point")+
  scale_fill_viridis_d(direction=-1)+
  scale_color_viridis_d(direction=-1)
dev.print(png,fil="vmr_time_box.png",type="quartz",antialias="default",width=6.5,
          height=6.5,units="in",res=1300)

#Plot VMR across time but this time wiht phylum as axis and time as groups
vmrar_box_all=ggplot(Phage_host_count_mergemar, aes(y=logVMR,x=fct_rev(Phylum)))+
  geom_boxplot(aes(fill=rev(Time)),
               outlier.shape = NA,
               size=1,
               alpha=0.7)+
  scale_fill_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6","#cd9fb2"),name="Time") +
  #  geom_point(aes(color=Phylum),position=position_jitterdodge(0.2),pch=19,
  #             alpha=0.9,
  #             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMAR ("~log[10]~')'),
       x="Sampling Time Point")+
  coord_flip()
#  scale_fill_viridis_d(direction=-1)+
#  scale_color_viridis_d(direction=-1)
vmrar_box_all
ggsave("/Users/cissell/Desktop/vmr_time_box.svg",width=6.5,height=8.5,units="in",vmrar_box_all)

##################
##Class Resolved
##################
#Create dataframe to merge back after summarizing
D1_host_lostmar_class=D1_count_tpmmar%>%
  select(c(Host,R1_1,R1_2,R1_3,R1_4,R1_5,Class))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D1_count_tpm2mar_class=D1_count_tpmmar%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P1_1","P1_2","P1_3","P1_4","P1_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D1_count_mergermar_class=merge(D1_count_tpm2mar_class,D1_host_lostmar_class,by="Host",all=T)
D1_count_mergermar_class$Class=as.factor(D1_count_mergermar_class$Class)
str(D1_count_mergermar_class)
#Summarize within class
D1_phyla_countsmar_class=D1_count_mergermar_class%>%
  select(-c(Host))%>%
  group_by(Class)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D1_phyla_countsmar_class$P1_1)+sum(D1_phyla_countsmar_class$R1_1)
sum(D1_phyla_countsmar_class$P1_2)+sum(D1_phyla_countsmar_class$R1_2)
sum(D1_phyla_countsmar_class$P1_3)+sum(D1_phyla_countsmar_class$R1_3)
sum(D1_phyla_countsmar_class$P1_4)+sum(D1_phyla_countsmar_class$R1_4)
sum(D1_phyla_countsmar_class$P1_5)+sum(D1_phyla_countsmar_class$R1_5)
#Nice

#Create each individual VMR and logVMR for each pair
D1_count_vmrsmar_class=D1_phyla_countsmar_class%>%
  mutate(VMR1=P1_1/R1_1)%>%
  mutate(logVMR1=log10(P1_1/R1_1))%>%
  mutate(VMR2=P1_2/R1_2)%>%
  mutate(logVMR2=log10(P1_2/R1_2))%>%
  mutate(VMR3=P1_3/R1_3)%>%
  mutate(logVMR3=log10(P1_3/R1_3))%>%
  mutate(VMR4=P1_4/R1_4)%>%
  mutate(logVMR4=log10(P1_4/R1_4))%>%
  mutate(VMR5=P1_5/R1_5)%>%
  mutate(logVMR5=log10(P1_5/R1_5))%>%
  dplyr::rename(P_1=P1_1,
                P_2=P1_2,
                P_3=P1_3,
                P_4=P1_4,
                P_5=P1_5,
                R_1=R1_1,
                R_2=R1_2,
                R_3=R1_3,
                R_4=R1_4,
                R_5=R1_5,
                VMR_1=VMR1,
                VMR_2=VMR2,
                VMR_3=VMR3,
                VMR_4=VMR4,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_2=logVMR2,
                logVMR_3=logVMR3,
                logVMR_4=logVMR4,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D1")

#Pivot longer for later merging and plotting
D1_count_vmr_longmar_class=D1_count_vmrsmar_class%>%
  pivot_longer(
    cols=-c(Class,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D1_count_vmr_longmar_class=as.data.frame(D1_count_vmr_longmar_class)
D1_count_vmr_longmar_class$Mat=as.factor(D1_count_vmr_longmar_class$Mat)
D1_count_vmr_longmar_class$Time=as.factor(D1_count_vmr_longmar_class$Time)
str(D1_count_vmr_longmar_class)

####D2
#Create dataframe to merge back after summarizing
D2_host_lostmar_class=D2_count_tpmmar%>%
  select(c(Host,R2_1,R2_2,R2_3,R2_4,R2_5,Class))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D2_count_tpm2mar_class=D2_count_tpmmar%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P2_1","P2_2","P2_3","P2_4","P2_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D2_count_mergermar_class=merge(D2_count_tpm2mar_class,D2_host_lostmar_class,by="Host",all=T)
D2_count_mergermar_class$Class=as.factor(D2_count_mergermar_class$Class)
str(D2_count_mergermar_class)
#Summarize within class
D2_phyla_countsmar_class=D2_count_mergermar_class%>%
  select(-c(Host))%>%
  group_by(Class)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D2_phyla_countsmar_class$P2_1)+sum(D2_phyla_countsmar_class$R2_1)
sum(D2_phyla_countsmar_class$P2_2)+sum(D2_phyla_countsmar_class$R2_2)
sum(D2_phyla_countsmar_class$P2_3)+sum(D2_phyla_countsmar_class$R2_3)
sum(D2_phyla_countsmar_class$P2_4)+sum(D2_phyla_countsmar_class$R2_4)
sum(D2_phyla_countsmar_class$P2_5)+sum(D2_phyla_countsmar_class$R2_5)
#Nice

#Create each individual VMR and logVMR for each pair
D2_count_vmrsmar_class=D2_phyla_countsmar_class%>%
  mutate(VMR1=P2_1/R2_1)%>%
  mutate(logVMR1=log10(P2_1/R2_1))%>%
  mutate(VMR2=P2_2/R2_2)%>%
  mutate(logVMR2=log10(P2_2/R2_2))%>%
  mutate(VMR3=P2_3/R2_3)%>%
  mutate(logVMR3=log10(P2_3/R2_3))%>%
  mutate(VMR4=P2_4/R2_4)%>%
  mutate(logVMR4=log10(P2_4/R2_4))%>%
  mutate(VMR5=P2_5/R2_5)%>%
  mutate(logVMR5=log10(P2_5/R2_5))%>%
  dplyr::rename(P_1=P2_1,
                P_2=P2_2,
                P_3=P2_3,
                P_4=P2_4,
                P_5=P2_5,
                R_1=R2_1,
                R_2=R2_2,
                R_3=R2_3,
                R_4=R2_4,
                R_5=R2_5,
                VMR_1=VMR1,
                VMR_2=VMR2,
                VMR_3=VMR3,
                VMR_4=VMR4,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_2=logVMR2,
                logVMR_3=logVMR3,
                logVMR_4=logVMR4,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D2")

#Pivot longer for later merging and plotting
D2_count_vmr_longmar_class=D2_count_vmrsmar_class%>%
  pivot_longer(
    cols=-c(Class,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D2_count_vmr_longmar_class=as.data.frame(D2_count_vmr_longmar_class)
D2_count_vmr_longmar_class$Mat=as.factor(D2_count_vmr_longmar_class$Mat)
D2_count_vmr_longmar_class$Time=as.factor(D2_count_vmr_longmar_class$Time)
str(D2_count_vmr_longmar_class)

#####R3
#Create dataframe to merge back after summarizing
D3_host_lostmar_class=D3_count_tpmmar%>%
  select(c(Host,R3_1,R3_2,R3_3,R3_4,R3_5,Class))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D3_count_tpm2mar_class=D3_count_tpmmar%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P3_1","P3_2","P3_3","P3_4","P3_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D3_count_mergermar_class=merge(D3_count_tpm2mar_class,D3_host_lostmar_class,by="Host",all=T)
D3_count_mergermar_class$Class=as.factor(D3_count_mergermar_class$Class)
str(D3_count_mergermar_class)
#Summarize within class
D3_phyla_countsmar_class=D3_count_mergermar_class%>%
  select(-c(Host))%>%
  group_by(Class)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D3_phyla_countsmar_class$P3_1)+sum(D3_phyla_countsmar_class$R3_1)
sum(D3_phyla_countsmar_class$P3_2)+sum(D3_phyla_countsmar_class$R3_2)
sum(D3_phyla_countsmar_class$P3_3)+sum(D3_phyla_countsmar_class$R3_3)
sum(D3_phyla_countsmar_class$P3_4)+sum(D3_phyla_countsmar_class$R3_4)
sum(D3_phyla_countsmar_class$P3_5)+sum(D3_phyla_countsmar_class$R3_5)
#Nice

#Create each individual VMR and logVMR for each pair
D3_count_vmrsmar_class=D3_phyla_countsmar_class%>%
  mutate(VMR1=P3_1/R3_1)%>%
  mutate(logVMR1=log10(P3_1/R3_1))%>%
  mutate(VMR2=P3_2/R3_2)%>%
  mutate(logVMR2=log10(P3_2/R3_2))%>%
  mutate(VMR3=P3_3/R3_3)%>%
  mutate(logVMR3=log10(P3_3/R3_3))%>%
  mutate(VMR4=P3_4/R3_4)%>%
  mutate(logVMR4=log10(P3_4/R3_4))%>%
  mutate(VMR5=P3_5/R3_5)%>%
  mutate(logVMR5=log10(P3_5/R3_5))%>%
  dplyr::rename(P_1=P3_1,
                P_2=P3_2,
                P_3=P3_3,
                P_4=P3_4,
                P_5=P3_5,
                R_1=R3_1,
                R_2=R3_2,
                R_3=R3_3,
                R_4=R3_4,
                R_5=R3_5,
                VMR_1=VMR1,
                VMR_2=VMR2,
                VMR_3=VMR3,
                VMR_4=VMR4,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_2=logVMR2,
                logVMR_3=logVMR3,
                logVMR_4=logVMR4,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D3")

#Pivot longer for later merging and plotting
D3_count_vmr_longmar_class=D3_count_vmrsmar_class%>%
  pivot_longer(
    cols=-c(Class,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D3_count_vmr_longmar_class=as.data.frame(D3_count_vmr_longmar_class)
D3_count_vmr_longmar_class$Mat=as.factor(D3_count_vmr_longmar_class$Mat)
D3_count_vmr_longmar_class$Time=as.factor(D3_count_vmr_longmar_class$Time)
str(D3_count_vmr_longmar_class)

####D4
#Create dataframe to merge back after summarizing
D4_host_lostmar_class=D4_count_tpmmar%>%
  select(c(Host,R4_1,R4_2,R4_3,R4_4,R4_5,Class))%>%
  group_by(Host)%>%
  filter(!duplicated(Host))
#Summarize within host
D4_count_tpm2mar_class=D4_count_tpmmar%>%
  select(-Phage,PLength,HLength)%>%
  group_by(Host)%>%
  summarise_at(c("P4_1","P4_2","P4_3","P4_4","P4_5"),sum,na.rm=T)%>%
  ungroup()
#Merge host tpm back in
D4_count_mergermar_class=merge(D4_count_tpm2mar_class,D4_host_lostmar_class,by="Host",all=T)
D4_count_mergermar_class$Class=as.factor(D4_count_mergermar_class$Class)
str(D4_count_mergermar_class)
#Summarize within class
D4_phyla_countsmar_class=D4_count_mergermar_class%>%
  select(-c(Host))%>%
  group_by(Class)%>%
  summarise_all(sum)
#Check for sense in TPM just to be sure
sum(D4_phyla_countsmar_class$P4_1)+sum(D4_phyla_countsmar_class$R4_1)
sum(D4_phyla_countsmar_class$P4_2)+sum(D4_phyla_countsmar_class$R4_2)
sum(D4_phyla_countsmar_class$P4_3)+sum(D4_phyla_countsmar_class$R4_3)
sum(D4_phyla_countsmar_class$P4_4)+sum(D4_phyla_countsmar_class$R4_4)
sum(D4_phyla_countsmar_class$P4_5)+sum(D4_phyla_countsmar_class$R4_5)
#Nice

#Create each individual VMR and logVMR for each pair
D4_count_vmrsmar_class=D4_phyla_countsmar_class%>%
  mutate(VMR1=P4_1/R4_1)%>%
  mutate(logVMR1=log10(P4_1/R4_1))%>%
  mutate(VMR2=P4_2/R4_2)%>%
  mutate(logVMR2=log10(P4_2/R4_2))%>%
  mutate(VMR3=P4_3/R4_3)%>%
  mutate(logVMR3=log10(P4_3/R4_3))%>%
  mutate(VMR4=P4_4/R4_4)%>%
  mutate(logVMR4=log10(P4_4/R4_4))%>%
  mutate(VMR5=P4_5/R4_5)%>%
  mutate(logVMR5=log10(P4_5/R4_5))%>%
  dplyr::rename(P_1=P4_1,
                P_2=P4_2,
                P_3=P4_3,
                P_4=P4_4,
                P_5=P4_5,
                R_1=R4_1,
                R_2=R4_2,
                R_3=R4_3,
                R_4=R4_4,
                R_5=R4_5,
                VMR_1=VMR1,
                VMR_2=VMR2,
                VMR_3=VMR3,
                VMR_4=VMR4,
                VMR_5=VMR5,
                logVMR_1=logVMR1,
                logVMR_2=logVMR2,
                logVMR_3=logVMR3,
                logVMR_4=logVMR4,
                logVMR_5=logVMR5)%>%
  mutate(Mat="D4")

#Pivot longer for later merging and plotting
D4_count_vmr_longmar_class=D4_count_vmrsmar_class%>%
  pivot_longer(
    cols=-c(Class,Mat),
    names_to=c(".value","Time"),
    names_sep="_")
D4_count_vmr_longmar_class=as.data.frame(D4_count_vmr_longmar_class)
D4_count_vmr_longmar_class$Mat=as.factor(D4_count_vmr_longmar_class$Mat)
D4_count_vmr_longmar_class$Time=as.factor(D4_count_vmr_longmar_class$Time)
str(D4_count_vmr_longmar_class)

#########################
##All CLASSS
#########################
#Now combine all 4 mats into one dataframe
Phage_host_count_mergemar_class=rbind(D1_count_vmr_longmar_class,D2_count_vmr_longmar_class,D3_count_vmr_longmar_class,D4_count_vmr_longmar_class)
str(Phage_host_count_mergemar_class)
Phage_host_count_mergemar_class$Class=as.factor(Phage_host_count_mergemar_class$Class)

#Make log abundances a column for better liklihood calculation
Phage_host_count_mergemar_class=Phage_host_count_mergemar_class%>%
  mutate(logP=log10(P))%>%
  mutate(logH=log10(R))

#Spearman correlation
cor.test(Phage_host_count_mergemar_class$logH,Phage_host_count_mergemar_class$logP, method="spearman")
#S = 504900, p-value < 2.2e-16
#  rho = 0.8276374

#Plot overall
overall_p_hmar_class=ggplot(Phage_host_count_mergemar_class, aes(y=logP,x=logH))+
  geom_smooth(method=lm,se=T,color="black")+
  geom_point(alpha=0.8, aes(color=Class))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral RNA Abundance ("~log[10]~TPM~')'),
       x=expression("Host RNA Abundance ("~log[10]~TPM~')'))+
  scale_color_viridis_d(direction=-1)+
  theme(legend.position = "none")
overall_p_hmar_class
ggsave("/Users/cissell/Desktop/overall_wlsmar_class.svg",width=6.5,height=3.25,units="in",overall_p_hmar_class)

#Plot abundances per class
ggplot(Phage_host_count_mergemar_class, aes(y=log10(P),x=log10(R)))+
  geom_smooth(method=lm,se=F,
              aes(color=Class))+
  geom_point(alpha=0.8,
             aes(color=Class))+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("Viral Abundance ("~log[10]~TPM~')'),
       x=expression("Host Abundance ("~log[10]~TPM~')'))+
  scale_color_viridis_d(direction=-1)

#Model logVMAR across Class
clm3mar_class=lmer(data=Phage_host_count_mergemar_class,Phage_host_count_mergemar_class$logVMR~Phage_host_count_mergemar_class$Class+Phage_host_count_mergemar_class$Mat$Mat+(1|Phage_host_count_mergemar_class$Time))
summary(clm3mar_class)
#Test assumptions
library("DHARMa")
check_vmr_count_model3mar_class <- simulateResiduals(fittedModel = clm3mar_class, n = 999)
plot(check_vmr_count_model3mar_class)

Anova(clm3mar_class)#Significant differences by phyla and mat in VMAR


#plot logVMR across Class overall with pairwise significance
vmr_boxmar_class=ggplot(Phage_host_count_mergemar_class, aes(y=logVMR,x=fct_rev(Class)))+
  geom_hline(yintercept=0,
             linetype="dashed",
             col="black",
             alpha=0.6)+
  geom_boxplot(aes(fill=Class),
               outlier.shape = NA,
               size=1,
               alpha=0.9)+
  geom_point(aes(fill=Class),position=position_jitterdodge(0.2),colour="black",pch=21,
             alpha=0.7,
             size=1.5)+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y=expression("VMAR ("~log[10]~')'),
       x="Class")+
  scale_fill_viridis_d(direction=-1)+
  theme(legend.position="none")+
  coord_flip()
vmr_boxmar_class
ggsave("/Users/cissell/Desktop/vmr_boxmar_class.svg",width=6.5,height=5.5,units="in",vmr_boxmar_class)

##############
####Temporal Decay Plot
##############
#We will do per mat and combine into one plot
#Mat 1
#First we need to summarize ORF length and mapping to get single length and RNA abund per vOTU
rna_phagediel1_brayprep=rna_phagediel1 %>%
  select(-c(Geneid,Start,End)) %>%
  group_by(contig_id) %>%
  mutate(Length=sum(Length)) %>%
  mutate(P1_1=sum(P1_1)) %>%
  mutate(P1_2=sum(P1_2)) %>%
  mutate(P1_3=sum(P1_3)) %>%
  mutate(P1_4=sum(P1_4)) %>%
  mutate(P1_5=sum(P1_5)) %>%
  ungroup() %>%
  distinct(contig_id, .keep_all=T)
#Make new df
rna_phagediel1_brayprep_new=rna_phagediel1_brayprep
str(rna_phagediel1_brayprep_new)
#Put length in KB
rna_phagediel1_brayprep_new$Length=rna_phagediel1_brayprep_new$Length/1000
#Create RPK
rna_phagediel1_brayprep_new$P1_1=rna_phagediel1_brayprep_new$P1_1/rna_phagediel1_brayprep_new$Length
rna_phagediel1_brayprep_new$P1_2=rna_phagediel1_brayprep_new$P1_2/rna_phagediel1_brayprep_new$Length
rna_phagediel1_brayprep_new$P1_3=rna_phagediel1_brayprep_new$P1_3/rna_phagediel1_brayprep_new$Length
rna_phagediel1_brayprep_new$P1_4=rna_phagediel1_brayprep_new$P1_4/rna_phagediel1_brayprep_new$Length
rna_phagediel1_brayprep_new$P1_5=rna_phagediel1_brayprep_new$P1_5/rna_phagediel1_brayprep_new$Length

#create scaling factors
sumP11po_bray=(sum(rna_phagediel1_brayprep_new$P1_1)/1000000)
sumP12po_bray=(sum(rna_phagediel1_brayprep_new$P1_2)/1000000)
sumP13po_bray=(sum(rna_phagediel1_brayprep_new$P1_3)/1000000)
sumP14po_bray=(sum(rna_phagediel1_brayprep_new$P1_4)/1000000)
sumP15po_bray=(sum(rna_phagediel1_brayprep_new$P1_5)/1000000)

#Create TPM
rna_phagediel1_brayprep_new$P1_1=rna_phagediel1_brayprep_new$P1_1/sumP11po_bray
rna_phagediel1_brayprep_new$P1_2=rna_phagediel1_brayprep_new$P1_2/sumP12po_bray
rna_phagediel1_brayprep_new$P1_3=rna_phagediel1_brayprep_new$P1_3/sumP13po_bray
rna_phagediel1_brayprep_new$P1_4=rna_phagediel1_brayprep_new$P1_4/sumP14po_bray
rna_phagediel1_brayprep_new$P1_5=rna_phagediel1_brayprep_new$P1_5/sumP15po_bray

#Check for sense
sum(rna_phagediel1_brayprep_new$P1_1) #Good
sum(rna_phagediel1_brayprep_new$P1_2) #Good
sum(rna_phagediel1_brayprep_new$P1_3) #Good
sum(rna_phagediel1_brayprep_new$P1_4) #Good
sum(rna_phagediel1_brayprep_new$P1_5) #Good

#Remove length and pivot wider
rna_phagediel1_brayprep_new1 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_1) %>%
  pivot_wider(names_from = contig_id, values_from =P1_1)
rna_phagediel1_brayprep_new2 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_2) %>%
  pivot_wider(names_from = contig_id, values_from =P1_2)
rna_phagediel1_brayprep_new3 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_3) %>%
  pivot_wider(names_from = contig_id, values_from =P1_3)
rna_phagediel1_brayprep_new4 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_4) %>%
  pivot_wider(names_from = contig_id, values_from =P1_4)
rna_phagediel1_brayprep_new5 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_5) %>%
  pivot_wider(names_from = contig_id, values_from =P1_5)
rna_phagediel1_brayprep_bind=rbind(rna_phagediel1_brayprep_new1,
                                   rna_phagediel1_brayprep_new2,
                                   rna_phagediel1_brayprep_new3,
                                   rna_phagediel1_brayprep_new4,
                                   rna_phagediel1_brayprep_new5)
rna_phagediel1_bray_dist_matrix=as.matrix(vegdist(rna_phagediel1_brayprep_bind,method="bray"))
View(rna_phagediel1_bray_dist_matrix)
rna_phagediel1_bray_dist_matrix=as.data.frame(rna_phagediel1_bray_dist_matrix)
#Add distance column for regression
#split and recombine
rna_phagediel1_bray_dist_matrix_1=rna_phagediel1_bray_dist_matrix%>%
  select(1) %>%
  mutate(distance=c(0,6,12,18,24))
rna_phagediel1_bray_dist_matrix_1=dplyr::rename(rna_phagediel1_bray_dist_matrix_1,Dist=1)

rna_phagediel1_bray_dist_matrix_2=rna_phagediel1_bray_dist_matrix%>%
  select(2) %>%
  mutate(distance=c(6,0,6,12,18))
rna_phagediel1_bray_dist_matrix_2=dplyr::rename(rna_phagediel1_bray_dist_matrix_2,Dist=1)

rna_phagediel1_bray_dist_matrix_3=rna_phagediel1_bray_dist_matrix%>%
  select(3) %>%
  mutate(distance=c(12,6,0,6,12))
rna_phagediel1_bray_dist_matrix_3=dplyr::rename(rna_phagediel1_bray_dist_matrix_3,Dist=1)

rna_phagediel1_bray_dist_matrix_4=rna_phagediel1_bray_dist_matrix%>%
  select(4) %>%
  mutate(distance=c(18,12,6,0,6))
rna_phagediel1_bray_dist_matrix_4=dplyr::rename(rna_phagediel1_bray_dist_matrix_4,Dist=1)

rna_phagediel1_bray_dist_matrix_5=rna_phagediel1_bray_dist_matrix%>%
  select(5) %>%
  mutate(distance=c(24,18,12,6,0))
rna_phagediel1_bray_dist_matrix_5=dplyr::rename(rna_phagediel1_bray_dist_matrix_5,Dist=1)

rna_phagediel1_bray_distances=rbind(rna_phagediel1_bray_dist_matrix_1,
                                    rna_phagediel1_bray_dist_matrix_2,
                                    rna_phagediel1_bray_dist_matrix_3,
                                    rna_phagediel1_bray_dist_matrix_4,
                                    rna_phagediel1_bray_dist_matrix_5)
rna_phagediel1_bray_distances=rna_phagediel1_bray_distances[!duplicated(rna_phagediel1_bray_distances$Dist),]

rna_phagediel1_bray_distances = rna_phagediel1_bray_distances%>%
  mutate(mat="D1")

##Mat 2
#First we need to summarize ORF length and mapping to get single length and RNA abund per vOTU
rna_phagediel2_brayprep=rna_phagediel2 %>%
  select(-c(Geneid,Start,End)) %>%
  group_by(contig_id) %>%
  mutate(Length=sum(Length)) %>%
  mutate(P2_1=sum(P2_1)) %>%
  mutate(P2_2=sum(P2_2)) %>%
  mutate(P2_3=sum(P2_3)) %>%
  mutate(P2_4=sum(P2_4)) %>%
  mutate(P2_5=sum(P2_5)) %>%
  ungroup() %>%
  distinct(contig_id, .keep_all=T)
#Make new df
rna_phagediel2_brayprep_new=rna_phagediel2_brayprep
str(rna_phagediel2_brayprep_new)
#Put length in KB
rna_phagediel2_brayprep_new$Length=rna_phagediel2_brayprep_new$Length/1000
#Create RPK
rna_phagediel2_brayprep_new$P2_1=rna_phagediel2_brayprep_new$P2_1/rna_phagediel2_brayprep_new$Length
rna_phagediel2_brayprep_new$P2_2=rna_phagediel2_brayprep_new$P2_2/rna_phagediel2_brayprep_new$Length
rna_phagediel2_brayprep_new$P2_3=rna_phagediel2_brayprep_new$P2_3/rna_phagediel2_brayprep_new$Length
rna_phagediel2_brayprep_new$P2_4=rna_phagediel2_brayprep_new$P2_4/rna_phagediel2_brayprep_new$Length
rna_phagediel2_brayprep_new$P2_5=rna_phagediel2_brayprep_new$P2_5/rna_phagediel2_brayprep_new$Length

#create scaling factors
sumP21po_bray=(sum(rna_phagediel2_brayprep_new$P2_1)/1000000)
sumP22po_bray=(sum(rna_phagediel2_brayprep_new$P2_2)/1000000)
sumP23po_bray=(sum(rna_phagediel2_brayprep_new$P2_3)/1000000)
sumP24po_bray=(sum(rna_phagediel2_brayprep_new$P2_4)/1000000)
sumP25po_bray=(sum(rna_phagediel2_brayprep_new$P2_5)/1000000)

#Create TPM
rna_phagediel2_brayprep_new$P2_1=rna_phagediel2_brayprep_new$P2_1/sumP21po_bray
rna_phagediel2_brayprep_new$P2_2=rna_phagediel2_brayprep_new$P2_2/sumP22po_bray
rna_phagediel2_brayprep_new$P2_3=rna_phagediel2_brayprep_new$P2_3/sumP23po_bray
rna_phagediel2_brayprep_new$P2_4=rna_phagediel2_brayprep_new$P2_4/sumP24po_bray
rna_phagediel2_brayprep_new$P2_5=rna_phagediel2_brayprep_new$P2_5/sumP25po_bray

#Check for sense
sum(rna_phagediel2_brayprep_new$P2_1) #Good
sum(rna_phagediel2_brayprep_new$P2_2) #Good
sum(rna_phagediel2_brayprep_new$P2_3) #Good
sum(rna_phagediel2_brayprep_new$P2_4) #Good
sum(rna_phagediel2_brayprep_new$P2_5) #Good

#Remove length and pivot wider
rna_phagediel2_brayprep_new1 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_1) %>%
  pivot_wider(names_from = contig_id, values_from =P2_1)
rna_phagediel2_brayprep_new2 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_2) %>%
  pivot_wider(names_from = contig_id, values_from =P2_2)
rna_phagediel2_brayprep_new3 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_3) %>%
  pivot_wider(names_from = contig_id, values_from =P2_3)
rna_phagediel2_brayprep_new4 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_4) %>%
  pivot_wider(names_from = contig_id, values_from =P2_4)
rna_phagediel2_brayprep_new5 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_5) %>%
  pivot_wider(names_from = contig_id, values_from =P2_5)
rna_phagediel2_brayprep_bind=rbind(rna_phagediel2_brayprep_new1,
                                   rna_phagediel2_brayprep_new2,
                                   rna_phagediel2_brayprep_new3,
                                   rna_phagediel2_brayprep_new4,
                                   rna_phagediel2_brayprep_new5)
rna_phagediel2_bray_dist_matrix=as.matrix(vegdist(rna_phagediel2_brayprep_bind,method="bray"))
View(rna_phagediel2_bray_dist_matrix)
rna_phagediel2_bray_dist_matrix=as.data.frame(rna_phagediel2_bray_dist_matrix)
#Add distance column for regression
#split and recombine
rna_phagediel2_bray_dist_matrix_1=rna_phagediel2_bray_dist_matrix%>%
  select(1) %>%
  mutate(distance=c(0,6,12,18,24))
rna_phagediel2_bray_dist_matrix_1=dplyr::rename(rna_phagediel2_bray_dist_matrix_1,Dist=1)

rna_phagediel2_bray_dist_matrix_2=rna_phagediel2_bray_dist_matrix%>%
  select(2) %>%
  mutate(distance=c(6,0,6,12,18))
rna_phagediel2_bray_dist_matrix_2=dplyr::rename(rna_phagediel2_bray_dist_matrix_2,Dist=1)

rna_phagediel2_bray_dist_matrix_3=rna_phagediel2_bray_dist_matrix%>%
  select(3) %>%
  mutate(distance=c(12,6,0,6,12))
rna_phagediel2_bray_dist_matrix_3=dplyr::rename(rna_phagediel2_bray_dist_matrix_3,Dist=1)

rna_phagediel2_bray_dist_matrix_4=rna_phagediel2_bray_dist_matrix%>%
  select(4) %>%
  mutate(distance=c(18,12,6,0,6))
rna_phagediel2_bray_dist_matrix_4=dplyr::rename(rna_phagediel2_bray_dist_matrix_4,Dist=1)

rna_phagediel2_bray_dist_matrix_5=rna_phagediel2_bray_dist_matrix%>%
  select(5) %>%
  mutate(distance=c(24,18,12,6,0))
rna_phagediel2_bray_dist_matrix_5=dplyr::rename(rna_phagediel2_bray_dist_matrix_5,Dist=1)

rna_phagediel2_bray_distances=rbind(rna_phagediel2_bray_dist_matrix_1,
                                    rna_phagediel2_bray_dist_matrix_2,
                                    rna_phagediel2_bray_dist_matrix_3,
                                    rna_phagediel2_bray_dist_matrix_4,
                                    rna_phagediel2_bray_dist_matrix_5)
rna_phagediel2_bray_distances=rna_phagediel2_bray_distances[!duplicated(rna_phagediel2_bray_distances$Dist),]

rna_phagediel2_bray_distances = rna_phagediel2_bray_distances%>%
  mutate(mat="D2")

##Mat 3
#First we need to summarize ORF length and mapping to get single length and RNA abund per vOTU
rna_phagediel3_brayprep=rna_phagediel3 %>%
  select(-c(Geneid,Start,End)) %>%
  group_by(contig_id) %>%
  mutate(Length=sum(Length)) %>%
  mutate(P3_1=sum(P3_1)) %>%
  mutate(P3_2=sum(P3_2)) %>%
  mutate(P3_3=sum(P3_3)) %>%
  mutate(P3_4=sum(P3_4)) %>%
  mutate(P3_5=sum(P3_5)) %>%
  ungroup() %>%
  distinct(contig_id, .keep_all=T)
#Make new df
rna_phagediel3_brayprep_new=rna_phagediel3_brayprep
str(rna_phagediel3_brayprep_new)
#Put length in KB
rna_phagediel3_brayprep_new$Length=rna_phagediel3_brayprep_new$Length/1000
#Create RPK
rna_phagediel3_brayprep_new$P3_1=rna_phagediel3_brayprep_new$P3_1/rna_phagediel3_brayprep_new$Length
rna_phagediel3_brayprep_new$P3_2=rna_phagediel3_brayprep_new$P3_2/rna_phagediel3_brayprep_new$Length
rna_phagediel3_brayprep_new$P3_3=rna_phagediel3_brayprep_new$P3_3/rna_phagediel3_brayprep_new$Length
rna_phagediel3_brayprep_new$P3_4=rna_phagediel3_brayprep_new$P3_4/rna_phagediel3_brayprep_new$Length
rna_phagediel3_brayprep_new$P3_5=rna_phagediel3_brayprep_new$P3_5/rna_phagediel3_brayprep_new$Length

#create scaling factors
sumP31po_bray=(sum(rna_phagediel3_brayprep_new$P3_1)/1000000)
sumP32po_bray=(sum(rna_phagediel3_brayprep_new$P3_2)/1000000)
sumP33po_bray=(sum(rna_phagediel3_brayprep_new$P3_3)/1000000)
sumP34po_bray=(sum(rna_phagediel3_brayprep_new$P3_4)/1000000)
sumP35po_bray=(sum(rna_phagediel3_brayprep_new$P3_5)/1000000)

#Create TPM
rna_phagediel3_brayprep_new$P3_1=rna_phagediel3_brayprep_new$P3_1/sumP31po_bray
rna_phagediel3_brayprep_new$P3_2=rna_phagediel3_brayprep_new$P3_2/sumP32po_bray
rna_phagediel3_brayprep_new$P3_3=rna_phagediel3_brayprep_new$P3_3/sumP33po_bray
rna_phagediel3_brayprep_new$P3_4=rna_phagediel3_brayprep_new$P3_4/sumP34po_bray
rna_phagediel3_brayprep_new$P3_5=rna_phagediel3_brayprep_new$P3_5/sumP35po_bray

#Check for sense
sum(rna_phagediel3_brayprep_new$P3_1) #Good
sum(rna_phagediel3_brayprep_new$P3_2) #Good
sum(rna_phagediel3_brayprep_new$P3_3) #Good
sum(rna_phagediel3_brayprep_new$P3_4) #Good
sum(rna_phagediel3_brayprep_new$P3_5) #Good

#Remove length and pivot wider
rna_phagediel3_brayprep_new1 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_1) %>%
  pivot_wider(names_from = contig_id, values_from =P3_1)
rna_phagediel3_brayprep_new2 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_2) %>%
  pivot_wider(names_from = contig_id, values_from =P3_2)
rna_phagediel3_brayprep_new3 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_3) %>%
  pivot_wider(names_from = contig_id, values_from =P3_3)
rna_phagediel3_brayprep_new4 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_4) %>%
  pivot_wider(names_from = contig_id, values_from =P3_4)
rna_phagediel3_brayprep_new5 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_5) %>%
  pivot_wider(names_from = contig_id, values_from =P3_5)
rna_phagediel3_brayprep_bind=rbind(rna_phagediel3_brayprep_new1,
                                   rna_phagediel3_brayprep_new2,
                                   rna_phagediel3_brayprep_new3,
                                   rna_phagediel3_brayprep_new4,
                                   rna_phagediel3_brayprep_new5)
rna_phagediel3_bray_dist_matrix=as.matrix(vegdist(rna_phagediel3_brayprep_bind,method="bray"))
View(rna_phagediel3_bray_dist_matrix)
rna_phagediel3_bray_dist_matrix=as.data.frame(rna_phagediel3_bray_dist_matrix)
#Add distance column for regression
#split and recombine
rna_phagediel3_bray_dist_matrix_1=rna_phagediel3_bray_dist_matrix%>%
  select(1) %>%
  mutate(distance=c(0,6,12,18,24))
rna_phagediel3_bray_dist_matrix_1=dplyr::rename(rna_phagediel3_bray_dist_matrix_1,Dist=1)

rna_phagediel3_bray_dist_matrix_2=rna_phagediel3_bray_dist_matrix%>%
  select(2) %>%
  mutate(distance=c(6,0,6,12,18))
rna_phagediel3_bray_dist_matrix_2=dplyr::rename(rna_phagediel3_bray_dist_matrix_2,Dist=1)

rna_phagediel3_bray_dist_matrix_3=rna_phagediel3_bray_dist_matrix%>%
  select(3) %>%
  mutate(distance=c(12,6,0,6,12))
rna_phagediel3_bray_dist_matrix_3=dplyr::rename(rna_phagediel3_bray_dist_matrix_3,Dist=1)

rna_phagediel3_bray_dist_matrix_4=rna_phagediel3_bray_dist_matrix%>%
  select(4) %>%
  mutate(distance=c(18,12,6,0,6))
rna_phagediel3_bray_dist_matrix_4=dplyr::rename(rna_phagediel3_bray_dist_matrix_4,Dist=1)

rna_phagediel3_bray_dist_matrix_5=rna_phagediel3_bray_dist_matrix%>%
  select(5) %>%
  mutate(distance=c(24,18,12,6,0))
rna_phagediel3_bray_dist_matrix_5=dplyr::rename(rna_phagediel3_bray_dist_matrix_5,Dist=1)

rna_phagediel3_bray_distances=rbind(rna_phagediel3_bray_dist_matrix_1,
                                    rna_phagediel3_bray_dist_matrix_2,
                                    rna_phagediel3_bray_dist_matrix_3,
                                    rna_phagediel3_bray_dist_matrix_4,
                                    rna_phagediel3_bray_dist_matrix_5)
rna_phagediel3_bray_distances=rna_phagediel3_bray_distances[!duplicated(rna_phagediel3_bray_distances$Dist),]

rna_phagediel3_bray_distances = rna_phagediel3_bray_distances%>%
  mutate(mat="D3")

##Mat 4
#First we need to summarize ORF length and mapping to get single length and RNA abund per vOTU
rna_phagediel4_brayprep=rna_phagediel4 %>%
  select(-c(Geneid,Start,End)) %>%
  group_by(contig_id) %>%
  mutate(Length=sum(Length)) %>%
  mutate(P4_1=sum(P4_1)) %>%
  mutate(P4_2=sum(P4_2)) %>%
  mutate(P4_3=sum(P4_3)) %>%
  mutate(P4_4=sum(P4_4)) %>%
  mutate(P4_5=sum(P4_5)) %>%
  ungroup() %>%
  distinct(contig_id, .keep_all=T)
#Make new df
rna_phagediel4_brayprep_new=rna_phagediel4_brayprep
str(rna_phagediel4_brayprep_new)
#Put length in KB
rna_phagediel4_brayprep_new$Length=rna_phagediel4_brayprep_new$Length/1000
#Create RPK
rna_phagediel4_brayprep_new$P4_1=rna_phagediel4_brayprep_new$P4_1/rna_phagediel4_brayprep_new$Length
rna_phagediel4_brayprep_new$P4_2=rna_phagediel4_brayprep_new$P4_2/rna_phagediel4_brayprep_new$Length
rna_phagediel4_brayprep_new$P4_3=rna_phagediel4_brayprep_new$P4_3/rna_phagediel4_brayprep_new$Length
rna_phagediel4_brayprep_new$P4_4=rna_phagediel4_brayprep_new$P4_4/rna_phagediel4_brayprep_new$Length
rna_phagediel4_brayprep_new$P4_5=rna_phagediel4_brayprep_new$P4_5/rna_phagediel4_brayprep_new$Length

#create scaling factors
sumP41po_bray=(sum(rna_phagediel4_brayprep_new$P4_1)/1000000)
sumP42po_bray=(sum(rna_phagediel4_brayprep_new$P4_2)/1000000)
sumP43po_bray=(sum(rna_phagediel4_brayprep_new$P4_3)/1000000)
sumP44po_bray=(sum(rna_phagediel4_brayprep_new$P4_4)/1000000)
sumP45po_bray=(sum(rna_phagediel4_brayprep_new$P4_5)/1000000)

#Create TPM
rna_phagediel4_brayprep_new$P4_1=rna_phagediel4_brayprep_new$P4_1/sumP41po_bray
rna_phagediel4_brayprep_new$P4_2=rna_phagediel4_brayprep_new$P4_2/sumP42po_bray
rna_phagediel4_brayprep_new$P4_3=rna_phagediel4_brayprep_new$P4_3/sumP43po_bray
rna_phagediel4_brayprep_new$P4_4=rna_phagediel4_brayprep_new$P4_4/sumP44po_bray
rna_phagediel4_brayprep_new$P4_5=rna_phagediel4_brayprep_new$P4_5/sumP45po_bray

#Check for sense
sum(rna_phagediel4_brayprep_new$P4_1) #Good
sum(rna_phagediel4_brayprep_new$P4_2) #Good
sum(rna_phagediel4_brayprep_new$P4_3) #Good
sum(rna_phagediel4_brayprep_new$P4_4) #Good
sum(rna_phagediel4_brayprep_new$P4_5) #Good

#Remove length and pivot wider
rna_phagediel4_brayprep_new1 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_1) %>%
  pivot_wider(names_from = contig_id, values_from =P4_1)
rna_phagediel4_brayprep_new2 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_2) %>%
  pivot_wider(names_from = contig_id, values_from =P4_2)
rna_phagediel4_brayprep_new3 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_3) %>%
  pivot_wider(names_from = contig_id, values_from =P4_3)
rna_phagediel4_brayprep_new4 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_4) %>%
  pivot_wider(names_from = contig_id, values_from =P4_4)
rna_phagediel4_brayprep_new5 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_5) %>%
  pivot_wider(names_from = contig_id, values_from =P4_5)
rna_phagediel4_brayprep_bind=rbind(rna_phagediel4_brayprep_new1,
                                   rna_phagediel4_brayprep_new2,
                                   rna_phagediel4_brayprep_new3,
                                   rna_phagediel4_brayprep_new4,
                                   rna_phagediel4_brayprep_new5)
rna_phagediel4_bray_dist_matrix=as.matrix(vegdist(rna_phagediel4_brayprep_bind,method="bray"))
View(rna_phagediel4_bray_dist_matrix)
rna_phagediel4_bray_dist_matrix=as.data.frame(rna_phagediel4_bray_dist_matrix)
#Add distance column for regression
#split and recombine
rna_phagediel4_bray_dist_matrix_1=rna_phagediel4_bray_dist_matrix%>%
  select(1) %>%
  mutate(distance=c(0,6,12,18,24))
rna_phagediel4_bray_dist_matrix_1=dplyr::rename(rna_phagediel4_bray_dist_matrix_1,Dist=1)

rna_phagediel4_bray_dist_matrix_2=rna_phagediel4_bray_dist_matrix%>%
  select(2) %>%
  mutate(distance=c(6,0,6,12,18))
rna_phagediel4_bray_dist_matrix_2=dplyr::rename(rna_phagediel4_bray_dist_matrix_2,Dist=1)

rna_phagediel4_bray_dist_matrix_3=rna_phagediel4_bray_dist_matrix%>%
  select(3) %>%
  mutate(distance=c(12,6,0,6,12))
rna_phagediel4_bray_dist_matrix_3=dplyr::rename(rna_phagediel4_bray_dist_matrix_3,Dist=1)

rna_phagediel4_bray_dist_matrix_4=rna_phagediel4_bray_dist_matrix%>%
  select(4) %>%
  mutate(distance=c(18,12,6,0,6))
rna_phagediel4_bray_dist_matrix_4=dplyr::rename(rna_phagediel4_bray_dist_matrix_4,Dist=1)

rna_phagediel4_bray_dist_matrix_5=rna_phagediel4_bray_dist_matrix%>%
  select(5) %>%
  mutate(distance=c(24,18,12,6,0))
rna_phagediel4_bray_dist_matrix_5=dplyr::rename(rna_phagediel4_bray_dist_matrix_5,Dist=1)

rna_phagediel4_bray_distances=rbind(rna_phagediel4_bray_dist_matrix_1,
                                    rna_phagediel4_bray_dist_matrix_2,
                                    rna_phagediel4_bray_dist_matrix_3,
                                    rna_phagediel4_bray_dist_matrix_4,
                                    rna_phagediel4_bray_dist_matrix_5)
rna_phagediel4_bray_distances=rna_phagediel4_bray_distances[!duplicated(rna_phagediel4_bray_distances$Dist),]

rna_phagediel4_bray_distances = rna_phagediel4_bray_distances%>%
  mutate(mat="D4")

#Merge into 1 bray df for all mats
phage_bray_all=rbind(rna_phagediel1_bray_distances,
                     rna_phagediel2_bray_distances,
                     rna_phagediel3_bray_distances,
                     rna_phagediel4_bray_distances)

#Remove zeros because they aren't really meaningful
phage_bray_all=phage_bray_all%>%
  subset(!Dist==0)
phage_bray_all$mat=as.factor(phage_bray_all$mat)
#quick exploratory plot
plot(phage_bray_all$Dist~phage_bray_all$distance)

#Fit with GAMM
library(mgcv)
bray_phage_gamm=gam(Dist~s(distance,k=3),method="REML",data=phage_bray_all)
summary(bray_phage_gamm)
anova(bray_phage_gamm)

#Not significant, plot the points only, no smoother
#Plot Bray across time
bray_all=ggplot(phage_bray_all, aes(y=Dist,x=distance))+
  geom_point(aes(color=mat),
             size=3,
             alpha=0.9,
             position=position_jitterdodge(0.2),
             pch=19)+
  scale_color_manual(values=c("#90a295","#9c8cdb","#455765","#add8e6"),name="Mat") +
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.title=element_text(size=13,color="black"),
        axis.line = element_line(size=1.1))+
  labs(y="Bray-Curtis dissimilarity",
       x="Temporal distance (Hrs)")
bray_all
ggsave("/Users/cissell/Desktop/bray_all.svg",width=6.5,height=6.5,units="in",bray_all)

#Not do evenness for every single mat for every single time point
#Sort by abundance from each sample in separate df
###Mat 1
#time 1
rna_phagedielevenness11 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_1) %>%
  arrange(desc(P1_1))
rna_phagedielevenness11$contig_id=factor(rna_phagedielevenness11$contig_id, levels = rna_phagedielevenness11$contig_id)
#time 2
rna_phagedielevenness12 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_2) %>%
  arrange(desc(P1_2))
rna_phagedielevenness12$contig_id=factor(rna_phagedielevenness12$contig_id, levels = rna_phagedielevenness12$contig_id)
#time 3
rna_phagedielevenness13 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_3) %>%
  arrange(desc(P1_3))
rna_phagedielevenness13$contig_id=factor(rna_phagedielevenness13$contig_id, levels = rna_phagedielevenness13$contig_id)
#time 4
rna_phagedielevenness14 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_4) %>%
  arrange(desc(P1_4))
rna_phagedielevenness14$contig_id=factor(rna_phagedielevenness14$contig_id, levels = rna_phagedielevenness14$contig_id)
#time 5
rna_phagedielevenness15 = rna_phagediel1_brayprep_new %>%
  select(contig_id,P1_5) %>%
  arrange(desc(P1_5))
rna_phagedielevenness15$contig_id=factor(rna_phagedielevenness15$contig_id, levels = rna_phagedielevenness15$contig_id)

###Mat 2
#time 1
rna_phagedielevenness21 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_1) %>%
  arrange(desc(P2_1))
rna_phagedielevenness21$contig_id=factor(rna_phagedielevenness21$contig_id, levels = rna_phagedielevenness21$contig_id)
#time 2
rna_phagedielevenness22 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_2) %>%
  arrange(desc(P2_2))
rna_phagedielevenness22$contig_id=factor(rna_phagedielevenness22$contig_id, levels = rna_phagedielevenness22$contig_id)
#time 3
rna_phagedielevenness23 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_3) %>%
  arrange(desc(P2_3))
rna_phagedielevenness23$contig_id=factor(rna_phagedielevenness23$contig_id, levels = rna_phagedielevenness23$contig_id)
#time 4
rna_phagedielevenness24 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_4) %>%
  arrange(desc(P2_4))
rna_phagedielevenness24$contig_id=factor(rna_phagedielevenness24$contig_id, levels = rna_phagedielevenness24$contig_id)
#time 5
rna_phagedielevenness25 = rna_phagediel2_brayprep_new %>%
  select(contig_id,P2_5) %>%
  arrange(desc(P2_5))
rna_phagedielevenness25$contig_id=factor(rna_phagedielevenness25$contig_id, levels = rna_phagedielevenness25$contig_id)

###Mat 3
#time 1
rna_phagedielevenness31 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_1) %>%
  arrange(desc(P3_1))
rna_phagedielevenness31$contig_id=factor(rna_phagedielevenness31$contig_id, levels = rna_phagedielevenness31$contig_id)
#time 2
rna_phagedielevenness32 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_2) %>%
  arrange(desc(P3_2))
rna_phagedielevenness32$contig_id=factor(rna_phagedielevenness32$contig_id, levels = rna_phagedielevenness32$contig_id)
#time 3
rna_phagedielevenness33 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_3) %>%
  arrange(desc(P3_3))
rna_phagedielevenness33$contig_id=factor(rna_phagedielevenness33$contig_id, levels = rna_phagedielevenness33$contig_id)
#time 4
rna_phagedielevenness34 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_4) %>%
  arrange(desc(P3_4))
rna_phagedielevenness34$contig_id=factor(rna_phagedielevenness34$contig_id, levels = rna_phagedielevenness34$contig_id)
#time 5
rna_phagedielevenness35 = rna_phagediel3_brayprep_new %>%
  select(contig_id,P3_5) %>%
  arrange(desc(P3_5))
rna_phagedielevenness35$contig_id=factor(rna_phagedielevenness35$contig_id, levels = rna_phagedielevenness35$contig_id)

###Mat 4
#time 1
rna_phagedielevenness41 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_1) %>%
  arrange(desc(P4_1))
rna_phagedielevenness41$contig_id=factor(rna_phagedielevenness41$contig_id, levels = rna_phagedielevenness41$contig_id)
#time 2
rna_phagedielevenness42 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_2) %>%
  arrange(desc(P4_2))
rna_phagedielevenness42$contig_id=factor(rna_phagedielevenness42$contig_id, levels = rna_phagedielevenness42$contig_id)
#time 3
rna_phagedielevenness43 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_3) %>%
  arrange(desc(P4_3))
rna_phagedielevenness43$contig_id=factor(rna_phagedielevenness43$contig_id, levels = rna_phagedielevenness43$contig_id)
#time 4
rna_phagedielevenness44 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_4) %>%
  arrange(desc(P4_4))
rna_phagedielevenness44$contig_id=factor(rna_phagedielevenness44$contig_id, levels = rna_phagedielevenness44$contig_id)
#time 5
rna_phagedielevenness45 = rna_phagediel4_brayprep_new %>%
  select(contig_id,P4_5) %>%
  arrange(desc(P4_5))
rna_phagedielevenness45$contig_id=factor(rna_phagedielevenness45$contig_id, levels = rna_phagedielevenness45$contig_id)

#Now plot the abundance 'curves'
rank_11=ggplot(rna_phagedielevenness11, aes(y=log10(1+P1_1),x=contig_id))+
  geom_bar(stat="identity",
           color="#90a295",
           fill = "#90a295")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_11
ggsave("/Users/cissell/Desktop/even_11.svg",rank_11,width=1.625,height=1.8)

rank_12=ggplot(rna_phagedielevenness12, aes(y=log10(1+P1_2),x=contig_id))+
  geom_bar(stat="identity",
           color="#9c8cdb",
           fill = "#9c8cdb")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_12
ggsave("/Users/cissell/Desktop/even_12.svg",rank_12,width=1.625,height=1.8)

rank_13=ggplot(rna_phagedielevenness13, aes(y=log10(1+P1_3),x=contig_id))+
  geom_bar(stat="identity",
           color="#455765",
           fill = "#455765")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_13
ggsave("/Users/cissell/Desktop/even_13.svg",rank_13,width=1.625,height=1.8)

rank_14=ggplot(rna_phagedielevenness14, aes(y=log10(1+P1_4),x=contig_id))+
  geom_bar(stat="identity",
           color="#add8e6",
           fill = "#add8e6")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_14
ggsave("/Users/cissell/Desktop/even_14.svg",rank_14,width=1.625,height=1.8)

rank_15=ggplot(rna_phagedielevenness15, aes(y=log10(1+P1_5),x=contig_id))+
  geom_bar(stat="identity",
           color="#cd9fb2",
           fill = "#cd9fb2")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_15
ggsave("/Users/cissell/Desktop/even_15.svg",rank_15,width=1.625,height=1.8)



rank_21=ggplot(rna_phagedielevenness21, aes(y=log10(1+P2_1),x=contig_id))+
  geom_bar(stat="identity",
           color="#90a295",
           fill = "#90a295")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_21
ggsave("/Users/cissell/Desktop/even_21.svg",rank_21,width=1.625,height=1.8)

rank_22=ggplot(rna_phagedielevenness22, aes(y=log10(1+P2_2),x=contig_id))+
  geom_bar(stat="identity",
           color="#9c8cdb",
           fill = "#9c8cdb")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_22
ggsave("/Users/cissell/Desktop/even_22.svg",rank_22,width=1.625,height=1.8)

rank_23=ggplot(rna_phagedielevenness23, aes(y=log10(1+P2_3),x=contig_id))+
  geom_bar(stat="identity",
           color="#455765",
           fill = "#455765")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_23
ggsave("/Users/cissell/Desktop/even_23.svg",rank_23,width=1.625,height=1.8)

rank_24=ggplot(rna_phagedielevenness24, aes(y=log10(1+P2_4),x=contig_id))+
  geom_bar(stat="identity",
           color="#add8e6",
           fill = "#add8e6")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_24
ggsave("/Users/cissell/Desktop/even_24.svg",rank_24,width=1.625,height=1.8)

rank_25=ggplot(rna_phagedielevenness25, aes(y=log10(1+P2_5),x=contig_id))+
  geom_bar(stat="identity",
           color="#cd9fb2",
           fill = "#cd9fb2")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_25
ggsave("/Users/cissell/Desktop/even_25.svg",rank_25,width=1.625,height=1.8)


rank_31=ggplot(rna_phagedielevenness31, aes(y=log10(1+P3_1),x=contig_id))+
  geom_bar(stat="identity",
           color="#90a295",
           fill = "#90a295")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_31
ggsave("/Users/cissell/Desktop/even_31.svg",rank_31,width=1.625,height=1.8)

rank_32=ggplot(rna_phagedielevenness32, aes(y=log10(1+P3_2),x=contig_id))+
  geom_bar(stat="identity",
           color="#9c8cdb",
           fill = "#9c8cdb")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_32
ggsave("/Users/cissell/Desktop/even_32.svg",rank_32,width=1.625,height=1.8)

rank_33=ggplot(rna_phagedielevenness33, aes(y=log10(1+P3_3),x=contig_id))+
  geom_bar(stat="identity",
           color="#455765",
           fill = "#455765")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_33
ggsave("/Users/cissell/Desktop/even_33.svg",rank_33,width=1.625,height=1.8)

rank_34=ggplot(rna_phagedielevenness34, aes(y=log10(1+P3_4),x=contig_id))+
  geom_bar(stat="identity",
           color="#add8e6",
           fill = "#add8e6")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_34
ggsave("/Users/cissell/Desktop/even_34.svg",rank_34,width=1.625,height=1.8)

rank_35=ggplot(rna_phagedielevenness35, aes(y=log10(1+P3_5),x=contig_id))+
  geom_bar(stat="identity",
           color="#cd9fb2",
           fill = "#cd9fb2")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_35
ggsave("/Users/cissell/Desktop/even_35.svg",rank_35,width=1.625,height=1.8)



rank_41=ggplot(rna_phagedielevenness41, aes(y=log10(1+P4_1),x=contig_id))+
  geom_bar(stat="identity",
           color="#90a295",
           fill = "#90a295")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_41
ggsave("/Users/cissell/Desktop/even_41.svg",rank_41,width=1.625,height=1.8)

rank_42=ggplot(rna_phagedielevenness42, aes(y=log10(1+P4_2),x=contig_id))+
  geom_bar(stat="identity",
           color="#9c8cdb",
           fill = "#9c8cdb")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_42
ggsave("/Users/cissell/Desktop/even_42.svg",rank_42,width=1.625,height=1.8)

rank_43=ggplot(rna_phagedielevenness43, aes(y=log10(1+P4_3),x=contig_id))+
  geom_bar(stat="identity",
           color="#455765",
           fill = "#455765")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_43
ggsave("/Users/cissell/Desktop/even_43.svg",rank_43,width=1.625,height=1.8)

rank_44=ggplot(rna_phagedielevenness44, aes(y=log10(1+P4_4),x=contig_id))+
  geom_bar(stat="identity",
           color="#add8e6",
           fill = "#add8e6")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_44
ggsave("/Users/cissell/Desktop/even_44.svg",rank_44,width=1.625,height=1.8)

rank_45=ggplot(rna_phagedielevenness45, aes(y=log10(1+P4_5),x=contig_id))+
  geom_bar(stat="identity",
           color="#cd9fb2",
           fill = "#cd9fb2")+
  theme_classic()+
  theme(axis.text = element_text(size=12,color="black"),
        axis.line = element_line(size=1.1),
        axis.text.x = element_blank(),
        axis.title = element_blank())+
  theme(legend.position="none")
rank_45
ggsave("/Users/cissell/Desktop/even_45.svg",rank_45,width=1.625,height=1.8)
