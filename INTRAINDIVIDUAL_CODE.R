# INTRA-INDIVIDUAL PROJECT CX
# Author: Isabelle Laforest-Lapointe (2016)

# Can work with two different OTU table (the first is the complete one and the second contains
# only the OTUs that were represented at least by 20 sequences in the dataset to go trhough analyses faster)
setwd('/data/users/isabelle/Temperate_phyllo_2014/Chimeras/uclust/')
otus.biom.uclust.nc<-read_biom("uclust_non_chimeric_otu_table.biom")
otus.biom.uclust.nc.select<-read_biom("uclust_non_chimeric_otu_table_more_than_20")

# CONVERT BIOM TABLE TO COMMUNITY MATRIX

otus.comm.uclust.nc<-t(as.matrix(biom_data(otus.biom.uclust.nc)))
otus.comm.uclust.nc<-as.data.frame(otus.comm.uclust.nc)
dim(otus.comm.uclust.nc)

otus.comm.uclust.nc.select<-t(as.matrix(biom_data(otus.biom.uclust.nc.select)))
otus.comm.uclust.nc.select<-as.data.frame(otus.comm.uclust.nc.select)
dim(otus.comm.uclust.nc.select)

# CHANGING SAMPLE NAMES THAT ARE ERRONEOUS
which(row.names(otus.comm.uclust.nc.select)=="AS9ACSASutton4")
which(row.names(otus.comm.uclust.nc.select)=="AS10ACSASutton5")
which(row.names(otus.comm.uclust.nc.select)=="AS11ACSASutton6")
row.names(otus.comm.uclust.nc.select)[175]<-"AS9ACSASutton1"
row.names(otus.comm.uclust.nc.select)[114]<-"AS10ACSASutton2"
row.names(otus.comm.uclust.nc.select)[76]<-"AS11ACSASutton3"

which(row.names(otus.comm.uclust.nc)=="AS9ACSASutton4")
which(row.names(otus.comm.uclust.nc)=="AS10ACSASutton5")
which(row.names(otus.comm.uclust.nc)=="AS11ACSASutton6")
row.names(otus.comm.uclust.nc)[175]<-"AS9ACSASutton1"
row.names(otus.comm.uclust.nc)[114]<-"AS10ACSASutton2"
row.names(otus.comm.uclust.nc)[76]<-"AS11ACSASutton3"

# IMPORT TAXONOMY
# Note that we have to do some editing “by hand” of the taxonomy assignments file
# set tab, ; and _ as separators (can do this in a spreadsheet editor)
otus.taxonomy <- read.delim("/data/users/isabelle/Temperate_phyllo_2014/R_analyses/tax_blast_non_chimeric_otu_uclust.txt", row.names=1)
head(otus.taxonomy)
summary(otus.taxonomy)
otus.taxonomy[-c(which(is.na(otus.taxonomy$ABUNDANCE)==TRUE)),]->otus.taxonomy

## HOW TO TAKE OUT ALL CHLOROPLASTS FROM OUR DATA

dim(otus.taxonomy) # 45580 8
c<-which(otus.taxonomy$CLASS=="Chloroplast")
otus.taxonomy<-otus.taxonomy[-c,]
dim(otus.taxonomy) # 45297 8
dim(otus.comm.uclust.nc) #219 45580
dim(otus.comm.uclust.nc.select) #219 4784
otus.comm.uclust.nc<-otus.comm.uclust.nc[,row.names(otus.taxonomy)]
otus.taxonomy.select<-otus.taxonomy[colnames(otus.comm.uclust.nc.select),]
dim(otus.comm.uclust.nc.select) #219 4784
dim(otus.comm.uclust.nc) # 219 45272

# IMPORT METADATA
metadata_intraindividual<-read.csv("/data/users/isabelle/Temperate_phyllo_2014/Chimeras/uclust/metadata_intraindividual.csv", row.names=1)
metadata_interindividual<-read.csv("/data/users/isabelle/Temperate_phyllo_2014/Chimeras/uclust/metadata_interindividual.csv", row.names=1)
metadata_for_dataset3<-read.csv("/data/users/isabelle/Temperate_phyllo_2014/Chimeras/uclust/metadata_for_dataset3.csv", row.names=1)
metadata_for_nmds<-read.csv("/data/users/isabelle/Temperate_phyllo_2014/Chimeras/uclust/metadata_for_nmds.csv", row.names=1)

dim(metadata_intraindividual) #30  3
dim(metadata_interindividual) #30  4
dim(metadata_for_dataset3) #55  4
dim(metadata_for_nmds) #60 4

## MAKE SURE ALL SAMPLES MATCH IN DIFFERENT FILES BEFORE PROCEEDING

# Dataset 1 Intra-individual
comm.sub.intra <- otus.comm.uclust.nc[rownames(metadata_intraindividual),]
comm.sub.s.intra <- otus.comm.uclust.nc.select[rownames(metadata_intraindividual),]
comm.sub.intra<-comm.sub.intra[,-c(which(apply(comm.sub.intra,2,sum)<=4))]
comm.sub.s.intra<-comm.sub.s.intra[,-c(which(apply(comm.sub.s.intra,2,sum)<=4))]

# Dataset 2 Inter-individual and Inter-specific
comm.sub.inter <- otus.comm.uclust.nc[rownames(metadata_interindividual),]
comm.sub.s.inter <- otus.comm.uclust.nc.select[rownames(metadata_interindividual),]
comm.sub.inter<-comm.sub.inter[,-c(which(apply(comm.sub.inter,2,sum)<=4))]
comm.sub.s.inter<-comm.sub.s.inter[,-c(which(apply(comm.sub.s.inter,2,sum)<=4))]

# Dataset3 Intra-individual, Inter-individual and Inter-specific
comm.sub.dataset3 <- otus.comm.uclust.nc[rownames(metadata_for_dataset3),]
comm.sub.s.dataset3 <- otus.comm.uclust.nc.select[rownames(metadata_for_dataset3),]
comm.sub.s.dataset3<-comm.sub.s.dataset3[,-c(which(apply(comm.sub.s.dataset3,2,sum)<=4))]
comm.sub.dataset3<-comm.sub.dataset3[,-c(which(apply(comm.sub.dataset3,2,sum)<=4))]

# Dataset for NMDS (figure 3)
comm.sub.nmds <- otus.comm.uclust.nc[rownames(metadata_for_nmds)[1:55],]
comm.sub.nmds <- rbind(comm.sub.nmds,comm.sub.nmds[1:5,])
comm.sub.s.nmds <- otus.comm.uclust.nc.select[rownames(metadata_for_nmds)[1:55],]
comm.sub.s.nmds <- rbind(comm.sub.s.nmds,comm.sub.s.nmds[1:5,])
comm.sub.nmds<-comm.sub.nmds[,-c(which(apply(comm.sub.nmds,2,sum)<=4))]
comm.sub.s.nmds<-comm.sub.s.nmds[,-c(which(apply(comm.sub.s.nmds,2,sum)<=4))]

dim(comm.sub.intra) #30 3039
dim(comm.sub.s.intra) #30 2501
dim(comm.sub.inter) #30 3422
dim(comm.sub.s.inter) #30 2969
dim(comm.sub.dataset3) #55 5005
dim(comm.sub.s.dataset3) #55 3531
dim(comm.sub.nmds) #60 5160
dim(comm.sub.s.nmds) #60 3553

comm.taxo <- otus.taxonomy[colnames(comm.sub.dataset3),]
dim(comm.taxo) #5005    8
comm.taxo.s <- otus.taxonomy[colnames(comm.sub.s.dataset3),]
dim(comm.taxo.s) #3531    8

# INFORMATION ABOUT SEQUENCES IN SAMPLE

min(apply(comm.sub.intra,1,sum)) #6199
min(apply(comm.sub.s.intra,1,sum)) #6168
min(apply(comm.sub.inter,1,sum)) #8600
min(apply(comm.sub.s.inter,1,sum)) #8590
min(apply(comm.sub.dataset3,1,sum)) #6256
min(apply(comm.sub.s.dataset3,1,sum)) #6202
min(apply(comm.sub.nmds,1,sum)) #6259
min(apply(comm.sub.s.nmds,1,sum)) #6202

## RAREFACTION OF THE DATA

comm.subH3.rare5000.intra <- rrarefy(comm.sub.intra, sample=5000)
comm.subH3.s.rare5000.intra <- rrarefy(comm.sub.s.intra, sample=5000)
comm.subH3.rare6000.inter<- rrarefy(comm.sub.inter, sample=6000)
comm.subH3.s.rare6000.inter <- rrarefy(comm.sub.s.inter, sample=6000)
comm.subH3.rare5000.dataset3 <- rrarefy(comm.sub.dataset3, sample=5000)
comm.subH3.s.rare5000.dataset3 <- rrarefy(comm.sub.s.dataset3, sample=5000)
comm.subH3.rare5000.nmds <- rrarefy(comm.sub.nmds, sample=5000)
comm.subH3.s.rare5000.nmds <- rrarefy(comm.sub.s.nmds, sample=5000)

## ORDINATION ON THE RAREFIED DATA

comm.subH3.rare5000.intra.mds <- metaMDS(comm.subH3.rare5000.intra)
comm.subH3.s.rare5000.intra.mds <- metaMDS(comm.subH3.s.rare5000.intra)
comm.subH3.rare6000.inter.mds <- metaMDS(comm.subH3.rare6000.inter)
comm.subH3.s.rare6000.inter.mds <- metaMDS(comm.subH3.s.rare6000.inter)
comm.subH3.rare5000.dataset3.mds <- metaMDS(comm.subH3.rare5000.dataset3)
comm.subH3.s.rare5000.dataset3.mds <- metaMDS(comm.subH3.s.rare5000.dataset3)
comm.subH3.rare5000.nmds.mds <- metaMDS(comm.subH3.rare5000.nmds)
comm.subH3.s.rare5000.nmds.mds <- metaMDS(comm.subH3.s.rare5000.nmds)

##################################################################################################
########################################### PERMANOVA ############################################
##################################################################################################

# ADONIS: TEST ANOVA-LIKE INTERACTIONS AMONG VARIABLES
# PERMANOVA ON BRAY-CURTIS DISSIMILARITIES

# Dataset 1 Intra-individual
model_5Kintra<-adonis(comm.subH3.rare5000.intra ~ IND_ID+Canopy_Location, data=metadata_intraindividual)
model_5Kintra$aov.tab # 64%

# Dataset 2 Inter-individual and Inter-specific
model_6Kinter<-adonis(comm.subH3.rare6000.inter ~ Species, data=metadata_interindividual, strata=metadata_interindividual$Site)
model_6Kinter$aov.tab # 47%

# Dataset3 Intra-individual, Inter-individual and Inter-specific
model_5Kdataset3<-adonis(comm.subH3.rare5000.dataset3 ~ Species/IND_ID+Canopy_Location, data=metadata_for_dataset3)
model_5Kdataset3$aov.tab # 48%,6% and 31%

# FIGURE 3 (intra-individual and host species)

sc_other_ind<-scores(comm.subH3.rare5000.nmds.mds)[1:30,]
test <- qplot(scores(comm.subH3.rare5000.nmds.mds)[1:30,1], scores(comm.subH3.rare5000.nmds.mds)[1:30,2],
              colour=Species, geom="point", xlab="NMDS 1", ylab="NMDS 2",
              data=metadata_for_nmds[1:30,]) + 
  stat_ellipse(level=0.95, alpha = 1/2, geom="polygon", aes(fill=Species)) + 
  theme_bw() +
  theme(axis.title.x = element_text(face="bold",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.text.x  = element_text(size=14,color="black"),
        axis.text.y  = element_text(size=14,color="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.position="right") + 
  annotate("rect", xmin = sc_other_ind[1,1]-0.05, xmax = sc_other_ind[1,1]+0.05, ymin = sc_other_ind[1,2]-0.05, ymax = sc_other_ind[1,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[2,1]-0.05, xmax = sc_other_ind[2,1]+0.05, ymin = sc_other_ind[2,2]-0.05, ymax = sc_other_ind[2,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[3,1]-0.05, xmax = sc_other_ind[3,1]+0.05, ymin = sc_other_ind[3,2]-0.05, ymax = sc_other_ind[3,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[4,1]-0.05, xmax = sc_other_ind[4,1]+0.05, ymin = sc_other_ind[4,2]-0.05, ymax = sc_other_ind[4,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[5,1]-0.05, xmax = sc_other_ind[5,1]+0.05, ymin = sc_other_ind[5,2]-0.05, ymax = sc_other_ind[5,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[6,1]-0.05, xmax = sc_other_ind[6,1]+0.05, ymin = sc_other_ind[6,2]-0.05, ymax = sc_other_ind[6,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[7,1]-0.05, xmax = sc_other_ind[7,1]+0.05, ymin = sc_other_ind[7,2]-0.05, ymax = sc_other_ind[7,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[8,1]-0.05, xmax = sc_other_ind[8,1]+0.05, ymin = sc_other_ind[8,2]-0.05, ymax = sc_other_ind[8,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[9,1]-0.05, xmax = sc_other_ind[9,1]+0.05, ymin = sc_other_ind[9,2]-0.05, ymax = sc_other_ind[9,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[10,1]-0.05, xmax = sc_other_ind[10,1]+0.05, ymin = sc_other_ind[10,2]-0.05, ymax = sc_other_ind[10,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[11,1]-0.05, xmax = sc_other_ind[11,1]+0.05, ymin = sc_other_ind[11,2]-0.05, ymax = sc_other_ind[11,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[12,1]-0.05, xmax = sc_other_ind[12,1]+0.05, ymin = sc_other_ind[12,2]-0.05, ymax = sc_other_ind[12,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[13,1]-0.05, xmax = sc_other_ind[13,1]+0.05, ymin = sc_other_ind[13,2]-0.05, ymax = sc_other_ind[13,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[14,1]-0.05, xmax = sc_other_ind[14,1]+0.05, ymin = sc_other_ind[14,2]-0.05, ymax = sc_other_ind[14,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[15,1]-0.05, xmax = sc_other_ind[15,1]+0.05, ymin = sc_other_ind[15,2]-0.05, ymax = sc_other_ind[15,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[16,1]-0.05, xmax = sc_other_ind[16,1]+0.05, ymin = sc_other_ind[16,2]-0.05, ymax = sc_other_ind[16,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[17,1]-0.05, xmax = sc_other_ind[17,1]+0.05, ymin = sc_other_ind[17,2]-0.05, ymax = sc_other_ind[17,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[18,1]-0.05, xmax = sc_other_ind[18,1]+0.05, ymin = sc_other_ind[18,2]-0.05, ymax = sc_other_ind[18,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[19,1]-0.05, xmax = sc_other_ind[19,1]+0.05, ymin = sc_other_ind[19,2]-0.05, ymax = sc_other_ind[19,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[20,1]-0.05, xmax = sc_other_ind[20,1]+0.05, ymin = sc_other_ind[20,2]-0.05, ymax = sc_other_ind[20,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[21,1]-0.05, xmax = sc_other_ind[21,1]+0.05, ymin = sc_other_ind[21,2]-0.05, ymax = sc_other_ind[21,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[22,1]-0.05, xmax = sc_other_ind[22,1]+0.05, ymin = sc_other_ind[22,2]-0.05, ymax = sc_other_ind[22,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[23,1]-0.05, xmax = sc_other_ind[23,1]+0.05, ymin = sc_other_ind[23,2]-0.05, ymax = sc_other_ind[23,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[24,1]-0.05, xmax = sc_other_ind[24,1]+0.05, ymin = sc_other_ind[24,2]-0.05, ymax = sc_other_ind[24,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[25,1]-0.05, xmax = sc_other_ind[25,1]+0.05, ymin = sc_other_ind[25,2]-0.05, ymax = sc_other_ind[25,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[26,1]-0.05, xmax = sc_other_ind[26,1]+0.05, ymin = sc_other_ind[26,2]-0.05, ymax = sc_other_ind[26,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[27,1]-0.05, xmax = sc_other_ind[27,1]+0.05, ymin = sc_other_ind[27,2]-0.05, ymax = sc_other_ind[27,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[28,1]-0.05, xmax = sc_other_ind[28,1]+0.05, ymin = sc_other_ind[28,2]-0.05, ymax = sc_other_ind[28,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[29,1]-0.05, xmax = sc_other_ind[29,1]+0.05, ymin = sc_other_ind[29,2]-0.05, ymax = sc_other_ind[29,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[30,1]-0.05, xmax = sc_other_ind[30,1]+0.05, ymin = sc_other_ind[30,2]-0.05, ymax = sc_other_ind[30,2]+0.05, alpha = .2, colour="black")+
  geom_point(aes(scores(comm.subH3.rare5000.nmds.mds)[1:30,1], scores(comm.subH3.rare5000.nmds.mds)[1:30,2],
                 shape=Canopy_Location, fill=Species, size=4),metadata_for_nmds[1:30,]) +
  geom_point(aes(scores(comm.subH3.rare5000.nmds.mds)[31:60,1], scores(comm.subH3.rare5000.nmds.mds)[31:60,2],
                 shape=Canopy_Location, fill=Species, size=4),metadata_for_nmds[31:60,]) + 
  scale_fill_manual(name="Species",values = c("#2ca25f", "#de2d26", "#feb24c", "#756bb1", "#2c7bb6")) + 
  scale_color_manual(values = c("#2ca25f", "#de2d26", "#feb24c", "#756bb1", "#2c7bb6")) + 
  scale_shape_manual(name="Canopy Location",values=c(8,24,22,23,25,11),labels=c("Bottom","East","North","West","South","Top")) +
  scale_size(guide = "none") +
  guides(color=FALSE)
test

sc_other_ind<-scores(comm.subH3.s.rare5000.nmds.mds)[1:30,]
test <- qplot(scores(comm.subH3.s.rare5000.nmds.mds)[31:60,1], scores(comm.subH3.s.rare5000.nmds.mds)[31:60,2],
              colour=Species, geom="point", xlab="NMDS 1", ylab="NMDS 2",
              data=metadata_for_nmds[31:60,]) + 
  stat_ellipse(level=0.95, alpha = 1/2, geom="polygon", aes(fill=Species)) + 
  theme_bw() +
  theme(axis.title.x = element_text(face="bold",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.text.x  = element_text(size=14,color="black"),
        axis.text.y  = element_text(size=14,color="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size=14),
        legend.position="right") + 
  annotate("rect", xmin = sc_other_ind[1,1]-0.05, xmax = sc_other_ind[1,1]+0.05, ymin = sc_other_ind[1,2]-0.05, ymax = sc_other_ind[1,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[2,1]-0.05, xmax = sc_other_ind[2,1]+0.05, ymin = sc_other_ind[2,2]-0.05, ymax = sc_other_ind[2,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[3,1]-0.05, xmax = sc_other_ind[3,1]+0.05, ymin = sc_other_ind[3,2]-0.05, ymax = sc_other_ind[3,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[4,1]-0.05, xmax = sc_other_ind[4,1]+0.05, ymin = sc_other_ind[4,2]-0.05, ymax = sc_other_ind[4,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[5,1]-0.05, xmax = sc_other_ind[5,1]+0.05, ymin = sc_other_ind[5,2]-0.05, ymax = sc_other_ind[5,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[6,1]-0.05, xmax = sc_other_ind[6,1]+0.05, ymin = sc_other_ind[6,2]-0.05, ymax = sc_other_ind[6,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[7,1]-0.05, xmax = sc_other_ind[7,1]+0.05, ymin = sc_other_ind[7,2]-0.05, ymax = sc_other_ind[7,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[8,1]-0.05, xmax = sc_other_ind[8,1]+0.05, ymin = sc_other_ind[8,2]-0.05, ymax = sc_other_ind[8,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[9,1]-0.05, xmax = sc_other_ind[9,1]+0.05, ymin = sc_other_ind[9,2]-0.05, ymax = sc_other_ind[9,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[10,1]-0.05, xmax = sc_other_ind[10,1]+0.05, ymin = sc_other_ind[10,2]-0.05, ymax = sc_other_ind[10,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[11,1]-0.05, xmax = sc_other_ind[11,1]+0.05, ymin = sc_other_ind[11,2]-0.05, ymax = sc_other_ind[11,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[12,1]-0.05, xmax = sc_other_ind[12,1]+0.05, ymin = sc_other_ind[12,2]-0.05, ymax = sc_other_ind[12,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[13,1]-0.05, xmax = sc_other_ind[13,1]+0.05, ymin = sc_other_ind[13,2]-0.05, ymax = sc_other_ind[13,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[14,1]-0.05, xmax = sc_other_ind[14,1]+0.05, ymin = sc_other_ind[14,2]-0.05, ymax = sc_other_ind[14,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[15,1]-0.05, xmax = sc_other_ind[15,1]+0.05, ymin = sc_other_ind[15,2]-0.05, ymax = sc_other_ind[15,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[16,1]-0.05, xmax = sc_other_ind[16,1]+0.05, ymin = sc_other_ind[16,2]-0.05, ymax = sc_other_ind[16,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[17,1]-0.05, xmax = sc_other_ind[17,1]+0.05, ymin = sc_other_ind[17,2]-0.05, ymax = sc_other_ind[17,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[18,1]-0.05, xmax = sc_other_ind[18,1]+0.05, ymin = sc_other_ind[18,2]-0.05, ymax = sc_other_ind[18,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[19,1]-0.05, xmax = sc_other_ind[19,1]+0.05, ymin = sc_other_ind[19,2]-0.05, ymax = sc_other_ind[19,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[20,1]-0.05, xmax = sc_other_ind[20,1]+0.05, ymin = sc_other_ind[20,2]-0.05, ymax = sc_other_ind[20,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[21,1]-0.05, xmax = sc_other_ind[21,1]+0.05, ymin = sc_other_ind[21,2]-0.05, ymax = sc_other_ind[21,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[22,1]-0.05, xmax = sc_other_ind[22,1]+0.05, ymin = sc_other_ind[22,2]-0.05, ymax = sc_other_ind[22,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[23,1]-0.05, xmax = sc_other_ind[23,1]+0.05, ymin = sc_other_ind[23,2]-0.05, ymax = sc_other_ind[23,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[24,1]-0.05, xmax = sc_other_ind[24,1]+0.05, ymin = sc_other_ind[24,2]-0.05, ymax = sc_other_ind[24,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[25,1]-0.05, xmax = sc_other_ind[25,1]+0.05, ymin = sc_other_ind[25,2]-0.05, ymax = sc_other_ind[25,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[26,1]-0.05, xmax = sc_other_ind[26,1]+0.05, ymin = sc_other_ind[26,2]-0.05, ymax = sc_other_ind[26,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[27,1]-0.05, xmax = sc_other_ind[27,1]+0.05, ymin = sc_other_ind[27,2]-0.05, ymax = sc_other_ind[27,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[28,1]-0.05, xmax = sc_other_ind[28,1]+0.05, ymin = sc_other_ind[28,2]-0.05, ymax = sc_other_ind[28,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[29,1]-0.05, xmax = sc_other_ind[29,1]+0.05, ymin = sc_other_ind[29,2]-0.05, ymax = sc_other_ind[29,2]+0.05, alpha = .2, colour="black")+
  annotate("rect", xmin = sc_other_ind[30,1]-0.05, xmax = sc_other_ind[30,1]+0.05, ymin = sc_other_ind[30,2]-0.05, ymax = sc_other_ind[30,2]+0.05, alpha = .2, colour="black")+
  geom_point(aes(scores(comm.subH3.s.rare5000.nmds.mds)[31:60,1], scores(comm.subH3.s.rare5000.nmds.mds)[31:60,2],
                 shape=Canopy_Location, fill=Species, size=4),metadata_for_nmds[31:60,]) +
  geom_point(aes(scores(comm.subH3.s.rare5000.nmds.mds)[1:30,1], scores(comm.subH3.s.rare5000.nmds.mds)[1:30,2],
                 shape=Canopy_Location, fill=Species, size=4),metadata_for_nmds[1:30,]) + 
  scale_fill_manual(name="Species",values = c("#2ca25f", "#de2d26", "#feb24c", "#756bb1", "#2c7bb6")) + 
  scale_color_manual(values = c("#2ca25f", "#de2d26", "#feb24c", "#756bb1", "#2c7bb6")) + 
  scale_shape_manual(name="Canopy Location",values=c(8,24,22,23,25,11),labels=c("Bottom","East","North","West","South","Top")) +
  scale_size(guide = "none") +
  guides(color=FALSE)
test

##################################################################################################
########################################### ALPHA-DIV #########################################
##################################################################################################


x<-diversity(comm.sub.s.intra, index="shannon", MARGIN=1, base=exp(1))
metadata_intraindividual$alphadiv<-x
summary(fm1<-aov(alphadiv~Canopy_Location+Species, data=metadata_intraindividual))
post_oc<-TukeyHSD(fm1, "Species", ordered=TRUE)
post_oc
plot(post_oc)
plot(alphadiv~Canopy_Location+Species, data=metadata_intraindividual)
post_oc2<-TukeyHSD(fm1, "Canopy_Location", ordered=TRUE)
post_oc2
plot(post_oc2)

# Mean per species
mean_byspecies<-aggregate(alphadiv~Species, data=metadata_intraindividual,FUN=mean)
sd_byspecies<-aggregate(alphadiv~Species, data=metadata_intraindividual,FUN=sd)
#computation of the standard error of the mean
sem<-sd_byspecies$alphadiv/sqrt(dim(metadata_intraindividual)[1])
mean_byspecies$sem<-sem
mean_byspecies

##################################################################################################
############################################ BETADISPER ##########################################
##################################################################################################

# BETADISPER test of distance to centroids in different species
# Test of homogeneity of variance across groups

# Intra-individual (no longer included as plot)
# Yields a single significative difference between PIGL and BEPA

dis <- vegdist(comm.subH3.rare5000.intra,method="bray")
mod2 <- betadisper(dis, metadata_intraindividual$Species,type ="centroid") ## warnings
mod2
perm_mod2<-permutest(mod2, control = permControl(nperm = 999))
perm_mod2
anova(mod2)
boxplot(mod2)
TukeyHSD(mod2)
plot(TukeyHSD(mod2))

plot(mod2,main = "")
points(mod2$vectors[which(metadata_intraindividual$Species=="ABBA")[1:6],1:2], pch=8, col="#2ca25f", lwd=2, bg="#2ca25f", cex=2)
points(mod2$vectors[which(metadata_intraindividual$Species=="PIGL")[1:6],1:2], pch=8, col="#2c7bb6", lwd=2, bg="#2c7bb6", cex=2)
points(mod2$vectors[which(metadata_intraindividual$Species=="BEPA")[1:6],1:2], pch=8, col="#756bb1", lwd=2, bg="#756bb1", cex=2)
points(mod2$vectors[which(metadata_intraindividual$Species=="ACSA")[1:6],1:2], pch=8, col="#feb24c", lwd=2, bg="#feb24c", cex=2)
points(mod2$vectors[which(metadata_intraindividual$Species=="ACRU")[1:6],1:2], pch=8, col="#de2d26", lwd=2, bg="#de2d26", cex=2)

# INTERSPECIFIC
# Significative difference between PIGL and ACSA and ACRU

dis <- vegdist(comm.subH3.rare6000.inter,method="bray")
mod2 <- betadisper(dis, metadata_interindividual$Species,type ="centroid") ## warnings
mod2
perm_mod2<-permutest(mod2, control = permControl(nperm = 999))
perm_mod2
anova(mod2)
boxplot(mod2)
TukeyHSD(mod2)
plot(TukeyHSD(mod2))

# TEST OF INTRA_INDIVIDUAL VARIATION ACROSS HOST SPECIES
metadata_for_nmds[which(metadata_for_nmds$Species=="ABBA"),]->abba_test
abba_comm<-comm.subH3.rare5000.nmds[which(metadata_for_nmds$Species=="ABBA"),]
abba_test$test[1:6]<-"intra"
abba_test$test[7:12]<-"inter"
dis <- vegdist(abba_comm,method="bray")
mod2 <- betadisper(dis, abba_test$test,type ="centroid") ## warnings
mod2
perm_mod2<-permutest(mod2, control = permControl(nperm = 999))
anova(mod2)
boxplot(mod2)
TukeyHSD(mod2)
plot(TukeyHSD(mod2))

metadata_for_nmds[which(metadata_for_nmds$Species=="ACRU"),]->acru_test
acru_comm<-comm.subH3.rare5000.nmds[which(metadata_for_nmds$Species=="ACRU"),]
acru_test$test[1:6]<-"intra"
acru_test$test[7:12]<-"inter"
dis <- vegdist(acru_comm,method="bray")
mod2 <- betadisper(dis, acru_test$test,type ="centroid") ## warnings
mod2
perm_mod2<-permutest(mod2, control = permControl(nperm = 999))
anova(mod2)
boxplot(mod2)
TukeyHSD(mod2)
plot(TukeyHSD(mod2))

metadata_for_nmds[which(metadata_for_nmds$Species=="ACSA"),]->acsa_test
acsa_comm<-comm.subH3.rare5000.nmds[which(metadata_for_nmds$Species=="ACSA"),]
acsa_test$test[1:6]<-"intra"
acsa_test$test[7:12]<-"inter"
dis <- vegdist(acsa_comm,method="bray")
mod2 <- betadisper(dis, acsa_test$test,type ="centroid") ## warnings
mod2
perm_mod2<-permutest(mod2, control = permControl(nperm = 999))
anova(mod2)
boxplot(mod2)
TukeyHSD(mod2)
plot(TukeyHSD(mod2))

metadata_for_nmds[which(metadata_for_nmds$Species=="BEPA"),]->bepa_test
bepa_comm<-comm.subH3.rare5000.nmds[which(metadata_for_nmds$Species=="BEPA"),]
bepa_test$test[1:6]<-"intra"
bepa_test$test[7:12]<-"inter"
dis <- vegdist(bepa_comm,method="bray")
mod2 <- betadisper(dis, bepa_test$test,type ="centroid") ## warnings
mod2
perm_mod2<-permutest(mod2, control = permControl(nperm = 999))
anova(mod2)
boxplot(mod2)
TukeyHSD(mod2)
plot(TukeyHSD(mod2))

#Average distance to centroid:
#  inter  intra 
#0.2799 0.1928 


metadata_for_nmds[which(metadata_for_nmds$Species=="PIGL"),]->pigl_test
pigl_comm<-comm.subH3.rare5000.nmds.mds[which(metadata_for_nmds$Species=="PIGL"),]
pigl_test$test[1:6]<-"intra"
pigl_test$test[7:12]<-"inter"
dis <- vegdist(pigl_comm,method="bray")
mod2 <- betadisper(dis, pigl_test$test,type ="centroid") ## warnings
mod2
perm_mod2<-permutest(mod2, control = permControl(nperm = 999))
anova(mod2)
boxplot(mod2)
TukeyHSD(mod2)
plot(TukeyHSD(mod2))

# FIGURE 4
dis <- vegdist(comm.subH3.rare5000.nmds,method="bray")
metadata_for_nmds$test[c(1,6:10)]<-"abba_intra"
metadata_for_nmds$test[c(2,11:15)]<-"acsa_intra"
metadata_for_nmds$test[c(3,16:20)]<-"acru_intra"
metadata_for_nmds$test[c(4,21:25)]<-"bepa_intra"
metadata_for_nmds$test[c(5,26:30)]<-"pigl_intra"
metadata_for_nmds$test[c(31:32,53:56)]<-"abba_inter"
metadata_for_nmds$test[c(33:34,41:43,57)]<-"acsa_inter"
metadata_for_nmds$test[c(35:36,44:46,58)]<-"acru_inter"
metadata_for_nmds$test[c(37:38,47:49,59)]<-"bepa_inter"
metadata_for_nmds$test[c(39:40,50:52,60)]<-"pigl_inter"
mod2 <- betadisper(dis, metadata_for_nmds$test,type ="centroid") ## warnings
bplot<-as.data.frame(cbind(as.character(mod2$group),as.numeric(as.character(mod2$distances))))
colnames(bplot)<-c("group","distances")
bplot$distances<-as.numeric(as.character(bplot$distances))
bplot$group<-factor(bplot$group, levels(bplot$group)[c(2,1,4,3,6,5,8,7,10,9)])
p<-ggplot(bplot,aes(group, distances))
p+geom_boxplot(aes(fill=group)) +
  scale_fill_manual(values = c("#2ca25f","#2ca25f","#de2d26","#de2d26","#feb24c",
                               "#feb24c","#756bb1","#756bb1","#2c7bb6","#2c7bb6"))+
  guides(fill=FALSE)+
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14))+
  scale_y_continuous(name="Distance to centroid")

##################################################################################################
############################################ Partial RDA #########################################
##################################################################################################

# Partialing out species/individual effect to highlight residual effect of canopy location
rda.comm.sub<-rda(comm.subH3.s.rare5000.intra~Condition(Species)+Canopy_Location,metadata_intraindividual)
plot(rda.comm.sub)
rda.comm.sub

##################################################################################################
######################################### CORE MICROBIOME ########################################
##################################################################################################

z<-decostand(comm.sub.dataset3, method="pa")
order<-apply(z,1,sum)

# MEAN OTU RICHNESS PER TREE

mean(order)
se<-sd(order)/sqrt(length(order));se

# % OTUs on single TREES and on ALL TREES

OTUs<-apply(z,2,sum)
perc_OTU_only_on_one_sample<-length(which(OTUs==1))/length(OTUs);perc_OTU_only_on_one_sample
perc_on_all_trees<-length(which(OTUs==40))/length(OTUs);perc_on_all_trees

comm.taxo$occ<-OTUs/dim(comm.sub.dataset3)[1]
sum_seq_per_otus<-apply(comm.sub.dataset3,2,sum)
comm.taxo$seq<-sum_seq_per_otus
core<-tail(comm.taxo[order(comm.taxo$occ),],n=42);core
sum(core$seq)/sum(sum_seq_per_otus)

sum(core$seq[which(core$CLASS=="Alphaproteobacteria")])/sum(core$seq)
sum(core$seq[which(core$CLASS=="Betaproteobacteria")])/sum(core$seq)
sum(core$seq[which(core$CLASS=="Gammaproteobacteria")])/sum(core$seq)
sum(core$seq[which(core$CLASS=="Cytophagia")])/sum(core$seq)
sum(core$seq[which(core$CLASS=="Acidobacteriia")])/sum(core$seq)
sum(core$seq[which(core$CLASS=="Actinobacteria")])/sum(core$seq)
sum(core$seq[which(core$ORDER=="Rhizobiales")])/sum(core$seq)
sum(core$seq[which(core$FAMILY=="Methylocystaceae")])/sum(core$seq[which(core$ORDER=="Rhizobiales")])

core$prop<-core$seq/sum(core$seq)
write.csv(core,"core_20june2016.csv")
##################################################################################################
############################################# BARPLOTS ###########################################
################################################################################################## 

# CREATING CLASS/ORDER/GENUS SUMMARY PER SAMPLE
sample_class<-matrix(nrow=1,ncol=7)
sample_class<-as.data.frame(sample_class)
colnames(sample_class)<-c("Group.1","x","prop","sample_id","species","position","ind_id")
sample_order<-matrix(nrow=1,ncol=7)
sample_order<-as.data.frame(sample_order)
colnames(sample_order)<-c("Group.1","x","prop","sample_id","species","position","ind_id")
sample_genus<-matrix(nrow=1,ncol=7)
sample_genus<-as.data.frame(sample_genus)
colnames(sample_genus)<-c("Group.1","x","prop","sample_id","species","position","ind_id")

for (i in 1:dim(metadata_for_dataset3)[1]){
  x<-comm.sub.dataset3[i,]
  y<-comm.taxo[colnames(x),]
  sum_seq_x<-apply(x,2,sum)
  y$seq<-sum_seq_x
  agg_class<-aggregate(y$seq,by=list(y$CLASS),FUN=sum)
  agg_order<-aggregate(y$seq,by=list(y$ORDER),FUN=sum)
  agg_genus<-aggregate(y$seq,by=list(y$GENERA),FUN=sum)
  agg_class$prop<-100*(agg_class$x/sum(agg_class$x))
  agg_order$prop<-100*(agg_order$x/sum(agg_order$x))
  agg_genus$prop<-100*(agg_genus$x/sum(agg_genus$x))
  agg_class<-agg_class[-which(agg_class$prop<1),]
  agg_class$sample_id<-row.names(metadata_for_dataset3)[i]
  agg_order$sample_id<-row.names(metadata_for_dataset3)[i]
  agg_genus$sample_id<-row.names(metadata_for_dataset3)[i]
  agg_class$species<-metadata_for_dataset3[i,1]
  agg_order$species<-metadata_for_dataset3[i,1]
  agg_genus$species<-metadata_for_dataset3[i,1]
  agg_class$position<-metadata_for_dataset3[i,2]
  agg_order$position<-metadata_for_dataset3[i,2]
  agg_genus$position<-metadata_for_dataset3[i,2]
  agg_class$ind_id<-metadata_for_dataset3[i,3]
  agg_order$ind_id<-metadata_for_dataset3[i,3]
  agg_genus$ind_id<-metadata_for_dataset3[i,3]
  sample_class<-rbind(sample_class,agg_class)
  sample_order<-rbind(sample_order,agg_order)
  sample_genus<-rbind(sample_genus,agg_genus)
}

sample_class$Group.1<-as.factor(sample_class$Group.1)
sample_class$sample_id<-as.factor(sample_class$sample_id)
sample_class[-1,]->sample_class
sample_class<-sample_class[,-2]
summary(sample_class)
sample_order$Group.1<-as.factor(sample_order$Group.1)
sample_order$sample_id<-as.factor(sample_order$sample_id)
sample_order[-1,]->sample_order
sample_order<-sample_order[,-2]
summary(sample_order)
sample_genus$Group.1<-as.factor(sample_genus$Group.1)
sample_genus$sample_id<-as.factor(sample_genus$sample_id)
sample_genus[-1,]->sample_genus
sample_genus<-sample_genus[,-2]
summary(sample_genus)

library(tidyr)
test_class<-spread(sample_class,Group.1,prop)
summary(test_class)
test_order<-spread(sample_order,Group.1,prop)
test_genus<-spread(sample_genus,Group.1,prop)

test_class$species<-as.factor(test_class$species)
test_class$position<-as.factor(test_class$position)
test_class$ind_id<-as.factor(test_class$ind_id)

test_class[is.na(test_class)]<-0
test_order[is.na(test_order)]<-0
test_genus[is.na(test_genus)]<-0

for (i in 5:21){
  test_class[which(test_class[,i]<1),i]<-0
}

test_class$Other<-NA

for (i in 1:dim(test_class)[1]){
  test_class$Other[i]<-100-sum(test_class[i,5:21])
}
summary(test_class)
dim(test_class)
apply(test_class[,5:22],1,sum)
apply(test_class[,5:22],2,sum)

summary(test_class)
gathercols<-colnames(test_class)[5:dim(test_class)[2]]
CLASS<-"CLASS"
PROPORTION<-"PROPORTION"
long_stat_class<-gather_(test_class,CLASS,PROPORTION,gathercols)
mean_class<-aggregate(PROPORTION~CLASS,data=long_stat_class,mean)
order(mean_class$PROPORTION,decreasing=TRUE)
long_stat_class$CLASS<-as.factor(long_stat_class$CLASS)
levels(long_stat_class$CLASS)
long_stat_class$CLASS<-factor(long_stat_class$CLASS,levels(long_stat_class$CLASS)[order(mean_class$PROPORTION,decreasing=TRUE)])
long_stat_class2<-matrix(nrow=dim(long_stat_class)[1],ncol=dim(long_stat_class)[2])
i<-1
c<-1
for (i in 1:length(levels(long_stat_class$CLASS))){
  temp<-which(long_stat_class$CLASS==levels(long_stat_class$CLASS)[i])
  f<-c+length(temp)-1
  long_stat_class2[c(c:f),]<-as.matrix(long_stat_class[temp,])
  c<-c+length(temp)
}

as.data.frame(long_stat_class2)->long_stat_class
summary(long_stat_class)
colnames(long_stat_class)<-c("sample_id","Species","Canopy_Location","ind_id","CLASS","PROPORTION")
long_stat_class$PROPORTION<-as.numeric(as.character(long_stat_class$PROPORTION))
long_stat_class$CLASS<-factor(long_stat_class$CLASS,levels(long_stat_class$CLASS)[order(mean_class$PROPORTION,decreasing=TRUE)])
levels(long_stat_class$CLASS)

aggregate(long_stat_class$PROPORTION,by=list(long_stat_class$Species,long_stat_class$CLASS),FUN=mean)

# FIGURE 1
q<-ggplot(long_stat_class,aes(Species,I(PROPORTION/8),fill=CLASS))+
  geom_bar(stat="identity")+
  xlab("")+ylab("PROPORTION OF TOTAL COMMUNITY (%)")+
  ggtitle("")+
  theme(axis.title.y = element_text(face="bold",size=20),
        axis.text.x  = element_text(size=24,color="black"),
        axis.text.y  = element_text(size=20,color="black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=18,face="bold"))
q + scale_fill_manual(name="BACTERIAL CLASS",values=c("#67001f","#b2182b","#d6604d","#f4a582",
                                                      "#fddbc7","#f7f7f7","#d1e5f0","black", "#92c5de", "#4393c3",
                                                      '#2166ac','#053061',"blue","#5e4fa2","#3288bd",
                                                      "#66c2a5","#2ca25f","#e6f598"))


# FIGURE 2 a)
# ABBA
long_stat_class_abba1<-long_stat_class[which(long_stat_class$ind_id=="ABBAGatineau1"),]

q<-ggplot(long_stat_class_abba1,aes(Canopy_Location,PROPORTION,fill=CLASS))+
  geom_bar(stat="identity")+
  xlab("")+ylab("PROPORTION OF TOTAL COMMUNITY (%)")+
  ggtitle("")+
  theme(axis.title.y = element_text(face="bold",size=20),
        axis.text.x  = element_text(size=24,color="black"),
        axis.text.y  = element_text(size=20,color="black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=18,face="bold"))
q + scale_fill_manual(name="BACTERIAL CLASS",values=c("#67001f","#b2182b","#d6604d","#f4a582",
                                                      "#fddbc7","#f7f7f7","#d1e5f0","black", "#92c5de", "#4393c3",
                                                      '#2166ac','#053061',"blue","#5e4fa2","#3288bd",
                                                      "#66c2a5","#2ca25f","#e6f598"))+
  scale_x_discrete(labels = c("bottom" = "B","E" = "E",
                              "N" = "N","O" = "W","top" = "T"))

# FIGURE 2 b)
# PIGL
long_stat_class_pigl1<-long_stat_class[which(long_stat_class$ind_id=="PIGLGatineau1"),]

q<-ggplot(long_stat_class_pigl1,aes(Canopy_Location,PROPORTION,fill=CLASS))+
  geom_bar(stat="identity")+
  xlab("")+ylab("PROPORTION OF TOTAL COMMUNITY (%)")+
  ggtitle("")+
  theme(axis.title.y = element_text(face="bold",size=20),
        axis.text.x  = element_text(size=24,color="black"),
        axis.text.y  = element_text(size=20,color="black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=18,face="bold"))
q + scale_fill_manual(name="BACTERIAL CLASS",values=c("#67001f","#b2182b","#d6604d","#f4a582",
                                                      "#fddbc7","#f7f7f7","#d1e5f0","black", "#92c5de", "#4393c3",
                                                      '#2166ac','#053061',"blue","#5e4fa2","#3288bd",
                                                      "#66c2a5","#2ca25f","#e6f598"))+
  scale_x_discrete(labels = c("bottom" = "B","E" = "E",
                              "N" = "N","O" = "W","top" = "T"))

# FIGURE 2 c)
# ACRU
long_stat_class_acru1<-long_stat_class[which(long_stat_class$ind_id=="ACRUGatineau1"),]

q<-ggplot(long_stat_class_acru1,aes(Canopy_Location,PROPORTION,fill=CLASS))+
  geom_bar(stat="identity")+
  xlab("")+ylab("PROPORTION OF TOTAL COMMUNITY (%)")+
  ggtitle("")+
  theme(axis.title.y = element_text(face="bold",size=20),
        axis.text.x  = element_text(size=24,color="black"),
        axis.text.y  = element_text(size=20,color="black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=18,face="bold"))
q + scale_fill_manual(name="BACTERIAL CLASS",values=c("#67001f","#b2182b","#d6604d","#f4a582",
                                                      "#fddbc7","#f7f7f7","#d1e5f0","black", "#92c5de", "#4393c3",
                                                      '#2166ac','#053061',"blue","#5e4fa2","#3288bd",
                                                      "#66c2a5","#2ca25f","#e6f598"))+
  scale_x_discrete(labels = c("bottom" = "B","E" = "E",
                              "N" = "N","O" = "W","top" = "T"))

# FIGURE 2 d)
# ACSA
long_stat_class_acsa1<-long_stat_class[which(long_stat_class$ind_id=="ACSAGatineau1"),]

q<-ggplot(long_stat_class_acsa1,aes(Canopy_Location,PROPORTION,fill=CLASS))+
  geom_bar(stat="identity")+
  xlab("")+ylab("PROPORTION OF TOTAL COMMUNITY (%)")+
  ggtitle("")+
  theme(axis.title.y = element_text(face="bold",size=20),
        axis.text.x  = element_text(size=24,color="black"),
        axis.text.y  = element_text(size=20,color="black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=18,face="bold"))
q + scale_fill_manual(name="BACTERIAL CLASS",values=c("#67001f","#b2182b","#d6604d","#f4a582",
                                                      "#fddbc7","#f7f7f7","#d1e5f0","black", "#92c5de", "#4393c3",
                                                      '#2166ac','#053061',"blue","#5e4fa2","#3288bd",
                                                      "#66c2a5","#2ca25f","#e6f598"))+
  scale_x_discrete(labels = c("bottom" = "B","E" = "E",
                              "N" = "N","O" = "W","top" = "T"))

# FIGURE 2 e)
# BEPA
long_stat_class_bepa1<-long_stat_class[which(long_stat_class$ind_id=="BEPAGatineau1"),]

q<-ggplot(long_stat_class_bepa1,aes(Canopy_Location,PROPORTION,fill=CLASS))+
  geom_bar(stat="identity")+
  xlab("")+ylab("PROPORTION OF TOTAL COMMUNITY (%)")+
  ggtitle("")+
  theme(axis.title.y = element_text(face="bold",size=20),
        axis.text.x  = element_text(size=24,color="black"),
        axis.text.y  = element_text(size=20,color="black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=18,face="bold"))
q + scale_fill_manual(name="BACTERIAL CLASS",values=c("#67001f","#b2182b","#d6604d","#f4a582",
                                                      "#fddbc7","#f7f7f7","#d1e5f0","black", "#92c5de", "#4393c3",
                                                      '#2166ac','#053061',"blue","#5e4fa2","#3288bd",
                                                      "#66c2a5","#2ca25f","#e6f598"))+
  scale_x_discrete(labels = c("bottom" = "B","E" = "E",
                              "N" = "N","O" = "W","top" = "T"))

##################################################################################################
############################################## END ###############################################
##################################################################################################

