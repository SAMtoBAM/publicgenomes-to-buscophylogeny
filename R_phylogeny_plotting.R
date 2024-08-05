#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

install.packages("remotes")
remotes::install_github("djw533/micro.gen.extra")

library(ape)
library(phytools)
library(ggtree)
library(aplot)
library(treeio)
library(stringr)
library(ggbreak)
library(ggpubr)
library(ggplot2)
library(forcats)
library(dplyr)
library(micro.gen.extra)

raw_tree <- ape::read.tree("~/cluster1/projects/Penicillium/phylogenetics/busco_aa_mafft_trimal.concat.veryfasttree_doubleprecision_treefile")

ggtree(raw_tree)

##this step takes the full phylogeny, first scales down the size of the roqueforti region then collapses it into a triangle which represents the length of the branches min and max at the different triangle tips
##it also readds the tip labels in order to overwrite the previously added tips which gets messed up with the collapsing and re-scaling
#scaleClade(p1, 698, .1) %>% collapse(698, 'mixed', fill="darkgreen") + geom_tiplab(size=5, as_ylab = TRUE)

##additionally can add a right sided label to the node too
#p1+
#  geom_hilight(node=698, type = "gradient", gradient.direction = 'rt',alpha = .5, to.bottom=T, fill="steelblue", extend=1)+
#  geom_cladelab(node=698, label="roqueforti", angle=0, fontsize=4, vjust=.5)+
#  geom_hilight(node=545, type = "gradient", gradient.direction = 'rt',alpha = .5, to.bottom=T, fill="darkgreen", extend=1)+
#  geom_cladelab(node=545, label="Camemberti", angle=0, fontsize=4, vjust=.5)


###TO ROTATE A NODE (MUST USE THE GGTREE SPECIFIC TAG DUE TO CONFLICTS)
#ggtree::rotate(p1, 644)

#######extracting all of the nodes for each species

##rerooting tree based on the edge before fumigatus
nodes=grep("fumigatus", raw_tree$tip.label)
root=MRCA(raw_tree, nodes)
rooted_tree=ape::root.phylo(raw_tree, node = root-1)
rooted_tree=TreeTools::Preorder(rooted_tree)

##just change the species for one strain of citrinum (need to do this at the base of all the data too)
#rooted_tree$tiplabel=gsub('Psteckii.P2648', 'Pcitrinum.P2648', rooted_tree$tiplabel)

##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
species= unique(do.call('rbind', strsplit(as.character(rooted_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
sp_clades=as.data.frame(c())


for (sp in species) {
  nodes=grep(paste(sp), rooted_tree$tip.label)
  clade=MRCA(rooted_tree, nodes)
  output=print(paste(sp,clade))
  sp_clades=rbind(sp_clades, output)
}

##rename header temporarily
colnames(sp_clades) = "temp"

##split the column into two seperate columns with appropriate headers
sp_clades2=sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##remove the unnamed, rubens, crustosum and simplicissimum
sp_clades2=subset(sp_clades2, species != "Prubens" & species != "Pherquei" & species != "Punnamed" )
##rename chrysogenum to chrysogenum/rubens 
sp_clades2$species=gsub('Pchrysogenum', 'Pchrysogenum/Prubens', sp_clades2$species)
sp_clades2$species=gsub('Pmalachiteum', 'Pmalachiteum/Pherquei', sp_clades2$species)
#sp_clades2$species=gsub('Psolitum', 'Psolitum/Pcrustosum', sp_clades2$species)
#sp_clades2$species=gsub('Pjanthinellum', 'Pjanthinellum/Psimplicissimum', sp_clades2$species)


##extract each column as a list MIGHT NOT NEED THIS
species2=sp_clades2[,1]
clades=sp_clades2[,2]

##force a strict order for the plotting to maintain it
sp_clades2$species <- factor(sp_clades2$species, levels = unique(sp_clades2$species))
sp_clades2$node <- factor(sp_clades2$node, levels = unique(sp_clades2$node))
sp_clades3=data.frame(node=as.numeric(paste(clades)), species=species2)

##extract each column as a list MIGHT NOT NEED THIS
species2=sp_clades3[,1]
clades=sp_clades3[,2]

##now plot with higlighting values on top of rerooted tree by using the clades and species
##added a line around each species defined region
##currently have the tips/genome labelled to just to help
ggtree(rooted_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node), colour="black", linetype="dotted", size=0.25 ,alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_tiplab(size=2, as_ylab = TRUE)+
  scale_fill_viridis_b()+
  geom_rootedge(rootedge = 0.01, linewidth=0.2)+
  geom_treescale(x=0.1, y=150, width=0.05, color='black')

##hilighting the species using cladelab too
ggtree(rooted_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_rootedge(rootedge = 0.01, linewidth=0.2)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )+
  geom_treescale(x=0.1, y=150, width=0.05, color='black', offset = 5)

##highlighting the clades and placing the genome name at the tip at a reasonable size
ggtree(rooted_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node), colour="black", linetype="dotted", size=0.25 ,alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_tiplab(size=.25)+
  scale_fill_viridis_b()+
  geom_rootedge(rootedge = 0.01, linewidth=0.2)+
  geom_treescale(x=0.1, y=150, width=0.05, color='black')

###now can add the starship information
##reading in DUF3435 data
captains=read.csv(file="cluster2/projects/Penicillium/starfish/allgenomes.genome_strain.TEMP_STARFISH_RESULTS.PLUS_PHYLOGENY_COMMENTS.PLUS_MYBs.ONLY_IN_PHYLO.tsv", header=T, sep='\t')

###need to install aplot and then use this to align a plot of the cpatins againt the tree
ggtree(rooted_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  scale_fill_viridis_b()+
  geom_rootedge(rootedge = 0.01, linewidth=0.2)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="grey" )+
  geom_treescale(x=0.1, y=150, width=0.05, color='black', offset = 5)
##alt tree with bigger text and no different colours for the boxed and a box edge colour
ggtree(rooted_tree, size=1) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node),alpha = 0.1, to.bottom=T, extend=0.01, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.01, linewidth=1)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 3, barcolour="grey" )+
  geom_treescale(x=0.1, y=150, width=0.05, color='black', offset = 5)


##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=rooted_tree$tip.label
captains=captains[grepl(paste(genomes, collapse="|"), captains$genome2),]

###the colour scheme for the isolation catagories=
## orange = clinical
## green = environment
## yellow = env-clinical
## purple = env-saline (only in aspergillus)
## light-blue = env-food
## blue = food-production
## black = unknown
## grey = NA

##plot tree to align with the starships
g=ggtree(rooted_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_rootedge(rootedge = 0.01, linewidth=0.2)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )+
  geom_treescale(x=0.1, y=150, width=0.05, color='black', offset = 5)
##plotting the starships
p1=ggplot(captains, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2))+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))
p1 %>% insert_left(g, width=5)




##regerenate the tree, removing the aspergillus branches after using them to root
aspsp=c("Afumigatus.E142_GCA_949125165.1", "Afumigatus.A1160_ASM2422042v1", "Afumigatus.A1163_ASM15014v1", "Afumigatus.Af293_ASM265v1", "Afumigatus.Afir964_ASM2875220v1", "Afumigatus.C6_GCA_949125545.1", "Afumigatus.C87_GCA_949125185.1", "Aawamori.IFM58123_GCA_003850985.1", "Aniger.ATCC13157_ATCC13157_v1", "Aniger.KJC3_ASM2978390v1", "Aniger.H915-1_ASM174190v1", "Aniger.WU-2020_GCA_024862975.1", "Aniger.KYF3_ASM2978392v1", "Aniger.JA-B-2022_ASM2958203v1", "Aluchuensis.IFO4308_GCA_016861625.1", "Aluchuensis.RIB2601_GCA_016865315.1", "Atubingensis.C2-2_ASM1061485v1","Atubingensis.WU-2223L_ASM1334032v1","Aterreus.ATCC20542_ASM1680841v1","Aterreus.M6925_ASM983442v1","Aflavus.A5P1_ASM2958205v1","Aflavus.CA14_ASM1478422v2","Aflavus.AF13_ASM1411748v1","Aflavus.NRRL3357-2_ASM901741v1","Aflavus.NRRL3357_ASM1411746v1","Aoryzae.RIB40_ASM18445v3","Aoryzae.KSS2_ASM803225v1","Aoryzae.SU-16_ASM985666v1","Aoryzae.BCC7051_ASM200794v1","Aoryzae.KBP3_ASM803205v1","Aparasiticus.MRI410_ASM2850576v1","Asojae.SMF134_ASM827498v1")
rooted_trim_tree=drop.tip(rooted_tree, aspsp)

##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
species= unique(do.call('rbind', strsplit(as.character(rooted_trim_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
sp_clades=as.data.frame(c())


for (sp in species) {
  nodes=grep(paste(sp), rooted_trim_tree$tip.label)
  clade=MRCA(rooted_trim_tree, nodes)
  output=print(paste(sp,clade))
  sp_clades=rbind(sp_clades, output)
}

##rename header temporarily
colnames(sp_clades) = "temp"

##split the column into two seperate columns with appropriate headers
sp_clades2=sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##remove the unnamed, rubens, crustosum and simplicissimum
sp_clades2=subset(sp_clades2, species != "Prubens" & species != "Pherquei" & species != "Punnamed" )
##rename chrysogenum to chrysogenum/rubens 
sp_clades2$species=gsub('Pchrysogenum', 'Pchrysogenum/Prubens', sp_clades2$species)
sp_clades2$species=gsub('Pmalachiteum', 'Pmalachiteum/Pherquei', sp_clades2$species)
#sp_clades2$species=gsub('Psolitum', 'Psolitum/Pcrustosum', sp_clades2$species)
#sp_clades2$species=gsub('Pjanthinellum', 'Pjanthinellum/Psimplicissimum', sp_clades2$species)


##extract each column as a list MIGHT NOT NEED THIS
species2=sp_clades2[,1]
clades=sp_clades2[,2]

##force a strict order for the plotting to maintain it
sp_clades2$species <- factor(sp_clades2$species, levels = unique(sp_clades2$species))
sp_clades2$node <- factor(sp_clades2$node, levels = unique(sp_clades2$node))
sp_clades3=data.frame(node=as.numeric(paste(clades)), species=species2)

##extract each column as a list MIGHT NOT NEED THIS
species2=sp_clades3[,1]
clades=sp_clades3[,2]

ggtree(rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )

g=ggtree(rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )


genomes=rooted_trim_tree$tip.label
captains2=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]
p1=ggplot(captains2, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))
p1 %>% insert_left(g, width=2)


##plotting all together just split by isolation orgigins
ggplot(subset(captains2, isolation_simple != "" & isolation_simple != "unknown"), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")


####one last time, generating trees, 
###but this time using the captain information file to remove a bunch of genomes
###removing genomes that are not able to be made public for roqueforti, solitum and crustosum
toremove=subset(captains, species == "Proqueforti" & supp == "" | species == "Psolitum" & supp == "" | species == "Pcrustosum" & supp == "" )$genome2
captains2_public=subset(captains2, species == "Proqueforti" & supp != "" | species == "Psolitum" & supp != "" | species == "Pcrustosum" & supp != ""  | species != "Proqueforti" & species != "Psolitum" &  species != "Pcrustosum" & supp == "" )

rooted_trim_tree_public=drop.tip(rooted_trim_tree, toremove, trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(rooted_trim_tree))

species= unique(do.call('rbind', strsplit(as.character(rooted_trim_tree_public$tip.label),'.',fixed=TRUE))[,1] )
sp_clades=as.data.frame(c())
for (sp in species) {
  nodes=grep(paste(sp), rooted_trim_tree_public$tip.label)
  clade=MRCA(rooted_trim_tree_public, nodes)
  output=print(paste(sp,clade))
  sp_clades=rbind(sp_clades, output)
}
colnames(sp_clades) = "temp"
sp_clades2=sp_clades %>% tidyr::separate(temp, c('species', 'node'))
sp_clades2=subset(sp_clades2, species != "Prubens" & species != "Pherquei" & species != "Punnamed" )
sp_clades2$species=gsub('Pchrysogenum', 'Pchrysogenum/Prubens', sp_clades2$species)
sp_clades2$species=gsub('Pmalachiteum', 'Pmalachiteum/Pherquei', sp_clades2$species)
species2=sp_clades2[,1]
clades=sp_clades2[,2]
sp_clades2$species <- factor(sp_clades2$species, levels = unique(sp_clades2$species))
sp_clades2$node <- factor(sp_clades2$node, levels = unique(sp_clades2$node))
sp_clades3=data.frame(node=as.numeric(paste(clades)), species=species2)
species2=sp_clades3[,1]
clades=sp_clades3[,2]

ggtree(rooted_trim_tree_public, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )
###saved as A4 pdf ('phylogeny_complete.penicillium_public')

##same as above but with the node support values for all branches added (need to manually remove those inside a species branching)
ggtree(rooted_trim_tree_public, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )+
  geom_text2(aes(label=label, subset=!isTip), hjust=-.2, size=0.5)
###saved as A4 pdf ('phylogeny_complete.penicillium_public.support')

##plot the same phylogeny but with only boxes around the species levels and with the genome name next to the tip
ggtree(rooted_trim_tree_public, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node), colour="black", linetype="dotted", size=0.25 ,alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_tiplab(size=.25)+
  scale_fill_viridis_b()+
  geom_rootedge(rootedge = 0.01, linewidth=0.2)+
  geom_treescale(x=0.1, y=150, width=0.05, color='black')
###saved as A4 pdf ('phylogeny_complete.penicillium_public.genome_lab')

g=ggtree(rooted_trim_tree_public, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )

p1=ggplot(captains2_public, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))
p1 %>% insert_left(g, width=2)
###saved as A4 pdf ('Penicillium_phylogeny_public.captains_isolation')


##CAN ALSO SAVE THE ABOVE TREE FOR ANOTHER PLOT WITH ASPERGILLUS ALONGSIDE
PENICILLIUMsupport=ggtree(rooted_trim_tree_public, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )+
  geom_text2(aes(label=label, subset=!isTip), hjust=-.2, size=0.5)
PENICILLIUM=ggtree(rooted_trim_tree_public, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )

##alternative plot with a smaller bar graph and reducing the x axis to remove the long (and false) branches

p1=ggplot(captains2_public, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  ylim(0,50)
p1 %>% insert_left(g, width=5)

##and saving as another i,age to join with aspergillus
PENICILLIUM2support=p1 %>% insert_left(PENICILLIUMsupport, width=5)
PENICILLIUM2=p1 %>% insert_left(PENICILLIUM, width=5)



##plotting all together just split by isolation origins
##removing unknown, NA and clinical (only had 1 clinical sample)
##also removing the outliers (that look to be truly errors) of janthinelleum
ggplot(subset(captains2_public, isolation_simple != "" & isolation_simple != "unknown" & isolation_simple != "clinical" & captains < 50), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")
###saved as svg 300x600 ('Penicillium_no_phylogeny.captains_isolation.svg')

##can randomly sample a subset of sets within the 'isolation_simple' group
captains2_public_random=do.call(rbind,replicate(1000, subset(captains2_public, isolation_simple != "" & isolation_simple != "unknown" & isolation_simple != "clinical" & isolation_simple != "environment-clinical" & isolation_simple != "environment-food" ) %>% group_by(isolation_simple) %>% slice_sample(n=25)  %>% summarise(mean_captains=mean(captains)), simplify=FALSE))
##plot this
ggplot(subset(captains2_public_random, isolation_simple != "" & isolation_simple != "unknown" & isolation_simple != "clinical"), aes(x=isolation_simple, y=mean_captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 45)+
  labs(x="")
###saved as svg 250x450 ('Penicillium_public_randomisation.captains_isolation')


#########CASEI/BIF/CAM ###############
###will want to take subsets to align alongside other data
##can use the fasciculata1 genome graph group to begin
nodes=grep("caseifulvum|biforme|camemberti", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 2)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 2)+
  ggplot2::xlim(-.0005, 0.0052)
##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)
##saved as A5


##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Ppalitans', 'Pfuscoglaucum', 'Pcaseifulvum', 'Pcamemberti', 'Pbiforme'))
#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")
  
##split by isolation environment WILL NOT USE THIS
ggplot(captains_subset, aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  theme_pubr(x.text.angle = 65)+
  facet_grid(species ~ .)



###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.0001, y=31, width=0.001, color='black', offset = 0.5)+
  ggplot2::xlim(-.0005, 0.0052)


g2=scaleClade(g, 204, .1) %>%
  scaleClade(296, .1) %>%
  scaleClade(286, .1) %>%
  scaleClade(156, .1) %>%
  scaleClade(177, .1) %>%
  collapse(204, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(296, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(286, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(177, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(156, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=204, label="P.biforme", color='black', fontsize=3)+
  geom_cladelabel(node=296, label="P.camemberti", color='black', fontsize=3)+
  geom_cladelabel(node=286, label="P.caseifulvum", color='black', fontsize=3)+
  geom_cladelabel(node=177, label="P.fuscoglaucum", color='black', fontsize=3)+
  geom_cladelabel(node=156, label="P.palitans", color='black', fontsize=3)+
  geom_rootedge(rootedge = -0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")
  
multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)



######### SOLITUM ###############
###will want to take subsets to align alongside other data
##can use the fasciculata1 genome graph group to begin
nodes=grep("solitum", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 1)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.005, y=45, width=0.001, color='black', offset = 1)+
  ggplot2::xlim(-.0005, 0.02)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)
##saved as A5

##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Psolitum','Pcrustosum'))
#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation environment 
ggplot(captains_subset, aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  theme_pubr(x.text.angle = 65)+
  facet_grid(species ~ .)
##for the figure save as 

###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.005, y=19, width=0.001, color='black', offset = .15)+
  ggplot2::xlim(-.0005, 0.02)


g2=scaleClade(g, 62, .1) %>%
  scaleClade(89, .1) %>%
  collapse(89, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(62, 'mixed', fill="firebrick", alpha=0.75)  +
  geom_cladelabel(node=89, label="P.crustosum", color='black', fontsize=3)+
  geom_cladelabel(node=62, label="P.solitum", color='black', fontsize=3)+
  geom_rootedge(rootedge = -0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

##for the same plot to add to the side
##same plot as above splitting the isolation origins of solitum
##except changing the x axis label angle and removing the x axis title
ggplot(subset(captains_subset, species=="Psolitum" & isolation_simple != "NA"), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  theme_pubr(x.text.angle = 35)+
  labs(x="")
##saved as svg 200x350 captains_boxplot.isolation.sol_only

#######REDOING THE WHOLE SOLITUM/CRUSTOSUM SET BUT THIS TIME WITH THE PUBLIC REDUCED DATASET




#########ROQUEFORTI ###############
###will want to take subsets to align alongside other data
##can use the fasciculata1 genome graph group to begin
nodes=grep("roqueforti", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 2)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
  ggplot2::xlim(-.0005, 0.012)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)
##saved as A5

##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Pcarneum', 'Ppsychrosexualis', 'Proqueforti'))
#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation environment
ggplot(captains_subset, aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")+
  facet_grid(~fct_rev(species) ~ . , )


###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.0006, y=35, width=0.001, color='black', offset = 1)+
  ggplot2::xlim(-.0005, 0.012)



g2=scaleClade(g, 258, .1) %>%
  scaleClade(508, .1) %>%
  scaleClade(509, .1) %>%
  collapse(258, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(508, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(509, 'mixed', fill="firebrick", alpha=0.75)  +
  geom_cladelabel(node=258, label="P.roqueforti", color='black', fontsize=3)+
  geom_cladelabel(node=508, label="P.psychrosexualis", color='black', fontsize=3)+
  geom_cladelabel(node=509, label="P.carneum", color='black', fontsize=3)+
  geom_rootedge(rootedge = -0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

##for the same plot to add to the side
##same plot as above splitting the isolation origins of solitum
##except changing the x axis label angle and removing the x axis title
captains_subset$cluster=factor(captains_subset$cluster, levels=c('Termignon','non-Roquefort' , 'contaminants', 'silage3', 'silage2', 'silage', 'wood', 'Roquefort'))
ggplot(subset(captains_subset, cluster != ""), aes(x=cluster, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  labs(x="")+
  coord_flip()
##saved as svg 200x350 captains_boxplot.isolation.roq_only??

###############
###above it is not clear how roqueforti's structure impacts the starship distribution so we can re subset for just roqueforti and split by cluster
nodes=grep("roqueforti", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 0)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.0001, y=200, width=0.001, color='black', offset = 5)+
  ggplot2::xlim(-.0005, 0.0035)

genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=cluster, group=genome2))+
  coord_flip()+
  theme_tree2()+
  labs(fill="cluster")+
  scale_fill_lancet()
p1 %>% insert_left(g)



##now plot the clusters grouped
captains_subset$cluster=factor(captains_subset$cluster, levels=c('Termignon','non-Roquefort' , 'contaminants', 'silage3', 'silage2', 'silage', 'wood', 'Roquefort'))
ggplot(subset(captains_subset, species== "Proqueforti" & cluster != "undefined" & cluster != ""), aes(x=cluster, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  labs(x="")+
  coord_flip()

###try to label, collapse etc a tree using the clusters of roquefort populations
##defined the clusters
subset_cluster=unique(do.call('rbind', strsplit(as.character(subset(captains_subset, cluster != "undefined")$cluster),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_cluster=as.data.frame(c())
for (cl in subset_cluster) {
  nodes <- c()
  genome=subset(captains_subset, cluster == cl)[,3]
  print(paste(genome))
  for(gen in genome) {
    node=grep(paste(gen), subset_tree$tip.label)
    nodes=append(nodes, node)
  }
  print(nodes)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(cl,clade))
  subset_sp_cluster=rbind(subset_sp_cluster, output)
}


colnames(subset_sp_cluster) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_cluster2=subset_sp_cluster %>% tidyr::separate(sep=" ", temp, c('cluster', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_cluster2=subset_sp_cluster2[,1]
subset_clades=subset_sp_cluster2[,2]
##force a strict order for the plotting to maintain it
subset_sp_cluster2$cluster <- factor(subset_sp_cluster2$cluster, levels = unique(subset_sp_cluster2$cluster))
subset_sp_cluster2$node <- factor(subset_sp_cluster2$node, levels = unique(subset_sp_cluster2$node))
subset_sp_cluster3=data.frame(node=as.numeric(paste(subset_clades)), cluster=subset_cluster2)
##extract each column as a list MIGHT NOT NEED THIS
subset_cluster2=subset_sp_cluster3[,2]
subset_clades=subset_sp_cluster3[,1]

##with the cluster labels
ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_cluster3 , mapping=aes(node=node, label=cluster), fontsize = 2, barcolour="grey")+
  geom_treescale(x=0.0001, y=200, width=0.001, color='black', offset = 5)+
  ggplot2::xlim(-.0005, 0.0035)

##now without cluster labels and collapsing for a simplified version to place next to the captain count
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.0000001, y=80, width=0.001, color='black', offset = 2)+
  ggplot2::xlim(-.0005, 0.0035)

g2=scaleClade(g, 389, .1) %>%
  scaleClade(416, .1) %>%
  scaleClade(255, .1) %>%
  scaleClade(328, .1) %>%
  scaleClade(320, .1) %>%
  scaleClade(339, .1) %>%
  scaleClade(478, .1) %>%
  collapse(389, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(416, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(255, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(328, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(320, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(339, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(478, 'mixed', fill="firebrick", alpha=0.75)  +
  geom_cladelabel(node=389, label="wood", color='black', fontsize=3)+
  geom_cladelabel(node=416, label="Roquefort", color='black', fontsize=3)+
  geom_cladelabel(node=255, label="non-Roquefort", color='black', fontsize=3)+
  geom_cladelabel(node=328, label="contaminants", color='black', fontsize=3)+
  geom_cladelabel(node=320, label="Termignon", color='black', fontsize=3)+
  geom_cladelabel(node=339, label="silage2", color='black', fontsize=3)+
  geom_cladelabel(node=478, label="silage", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
captains_subset$cluster=factor(captains_subset$cluster, levels=c('Termignon','non-Roquefort' , 'contaminants', 'silage2', 'silage', 'wood', 'Roquefort'))
p2=ggplot(subset(captains_subset, species== "Proqueforti" & cluster != "undefined" & cluster != ""), aes(x=cluster, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc(label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  labs(x="")+
  coord_flip()


multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

#######REDOING THE WHOLE ROQUEFORTI SET BUT THIS TIME WITH THE PUBLIC REDUCED DATASET
###will want to take subsets to align alongside other data
##can use the fasciculata1 genome graph group to begin
nodes=grep("roqueforti", rooted_trim_tree_public$tip.label)
##this can add a grouping label to the full tree
rooted_trim_tree_public = groupOTU(rooted_trim_tree_public, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_trim_tree_public, nodes)
subset_tree = tree_subset(rooted_trim_tree_public, clade, levels_back = 2)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.0025, y=30, width=0.001, color='black', offset = 1)+
  ggplot2::xlim(-.0005, 0.011)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains2_public[grepl(paste(genomes, collapse="|"), captains2_public$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)
##saved as A5 ('phylogeny_public.captains_barplot.A5.roq_psy_car')

##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Pcarneum', 'Ppsychrosexualis', 'Proqueforti'))
#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation environment
ggplot(captains_subset, aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")+
  facet_grid(~fct_rev(species) ~ . , )


###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.003, y=16, width=0.001, color='black', offset = .1)+
  ggplot2::xlim(-.0005, 0.011)

g2=scaleClade(g, 99, .1) %>%
  scaleClade(190, .1) %>%
  scaleClade(191, .1) %>%
  collapse(99, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(190, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(191, 'mixed', fill="firebrick", alpha=0.75)  +
  geom_cladelabel(node=99, label="P.roqueforti", color='black', fontsize=3)+
  geom_cladelabel(node=190, label="P.psychrosexualis", color='black', fontsize=3)+
  geom_cladelabel(node=191, label="P.carneum", color='black', fontsize=3)+
  geom_rootedge(rootedge = -0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide) ('phylogeny_public_reduced.captains_boxplot.THIN.roq_psych_carn')

##for the same plot to add to the side
##same plot as above splitting the isolation origins of solitum
##except changing the x axis label angle and removing the x axis title
captains_subset$cluster=factor(captains_subset$cluster, levels=c('Termignon','non-Roquefort' , 'contaminants', 'silage3', 'silage2', 'silage', 'wood', 'Roquefort'))
ggplot(subset(captains_subset, cluster != ""), aes(x=cluster, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  labs(x="")+
  coord_flip()

###############
###above it is not clear how roqueforti's structure impacts the starship distribution so we can re subset for just roqueforti and split by cluster
nodes=grep("roqueforti", rooted_trim_tree_public$tip.label)
##this can add a grouping label to the full tree
rooted_trim_tree_public = groupOTU(rooted_trim_tree_public, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_trim_tree_public, nodes)
subset_tree = tree_subset(rooted_trim_tree_public, clade, levels_back = 0)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.0001, y=50, width=0.001, color='black', offset = 1)+
  ggplot2::xlim(-.0005, 0.003)

genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=cluster, group=genome2))+
  coord_flip()+
  theme_tree2()+
  labs(fill="cluster")+
  scale_fill_lancet()
p1 %>% insert_left(g)
###saved as A5 ('phylogeny_public.captains_barplot.A5.roq_clusters')


##now plot the clusters grouped
captains_subset$cluster=factor(captains_subset$cluster, levels=c('Termignon','non-Roquefort' , 'contaminants', 'silage3', 'silage2', 'silage', 'Roquefort' , 'wood'))
ggplot(subset(captains_subset, species== "Proqueforti" & cluster != "undefined" & cluster != ""), aes(x=cluster, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  labs(x="")+
  coord_flip()

###try to label, collapse etc a tree using the clusters of roquefort populations
##defined the clusters
subset_cluster=unique(do.call('rbind', strsplit(as.character(subset(captains_subset, cluster != "undefined")$cluster),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_cluster=as.data.frame(c())
for (cl in subset_cluster) {
  nodes <- c()
  genome=subset(captains_subset, cluster == cl)[,4]
  print(paste(genome))
  for(gen in genome) {
    node=grep(paste(gen), subset_tree$tip.label)
    nodes=append(nodes, node)
  }
  print(nodes)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(cl,clade))
  subset_sp_cluster=rbind(subset_sp_cluster, output)
}


colnames(subset_sp_cluster) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_cluster2=subset_sp_cluster %>% tidyr::separate(sep=" ", temp, c('cluster', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_cluster2=subset_sp_cluster2[,1]
subset_clades=subset_sp_cluster2[,2]
##force a strict order for the plotting to maintain it
subset_sp_cluster2$cluster <- factor(subset_sp_cluster2$cluster, levels = unique(subset_sp_cluster2$cluster))
subset_sp_cluster2$node <- factor(subset_sp_cluster2$node, levels = unique(subset_sp_cluster2$node))
subset_sp_cluster3=data.frame(node=as.numeric(paste(subset_clades)), cluster=subset_cluster2)
##extract each column as a list MIGHT NOT NEED THIS
subset_cluster2=subset_sp_cluster3[,2]
subset_clades=subset_sp_cluster3[,1]

##with the cluster labels
ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_cluster3 , mapping=aes(node=node, label=cluster), fontsize = 2, barcolour="grey")+
  geom_treescale(x=0.0001, y=75, width=0.001, color='black', offset = 1)+
  ggplot2::xlim(-.0005, 0.0035)

##now without cluster labels and collapsing for a simplified version to place next to the captain count
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.000000005, y=35, width=0.001, color='black', offset = 0.25)+
  ggplot2::xlim(-.0005, 0.0035)

g2=scaleClade(g, 143, .1) %>%
  scaleClade(168, .1) %>%
  scaleClade(96, .1) %>%
  scaleClade(124, .1) %>%
  scaleClade(118, .1) %>%
  scaleClade(132, .1) %>%
  scaleClade(181, .1) %>%
  collapse(143, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(168, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(96, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(124, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(118, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(132, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(181, 'mixed', fill="firebrick", alpha=0.75)  +
  geom_cladelabel(node=143, label="wood", color='black', fontsize=3)+
  geom_cladelabel(node=168, label="Roquefort", color='black', fontsize=3)+
  geom_cladelabel(node=96, label="non-Roquefort", color='black', fontsize=3)+
  geom_cladelabel(node=124, label="contaminants", color='black', fontsize=3)+
  geom_cladelabel(node=118, label="Termignon", color='black', fontsize=3)+
  geom_cladelabel(node=132, label="silage2", color='black', fontsize=3)+
  geom_cladelabel(node=181, label="silage", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
captains_subset$cluster=factor(captains_subset$cluster, levels=c('Termignon','non-Roquefort' , 'contaminants', 'silage2', 'silage', 'Roquefort' , 'wood'))
p2=ggplot(subset(captains_subset, species== "Proqueforti" & cluster != "undefined" & cluster != ""), aes(x=cluster, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc(label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  labs(x="")+
  coord_flip()


multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide) ('phylogeny_public_reduced.captains_boxplot.THIN.roq_only')






#########NORDICUM ###############
###will want to take subsets to align alongside other data
nodes=grep("nordicum", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 2)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.0001, y=15, width=0.001, color='black', offset =.2)+
  ggplot2::xlim(-.0005, 0.03)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)


##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Pnordicum', 'Pverrucosum', 'Paurantiogriseum', 'Pfreii', 'Punnamed', 'Pviridicatum', 'Ppolonicum'))
#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation environment DON'T NEED FOR NORDICUM
#ggplot(captains_subset, aes(x=isolation_simple, y=captains))+
#  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
#  scale_colour_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
#  theme_pubr()+
#  coord_flip()+
#  labs(x="")+
#  facet_grid(~fct_rev(species) ~ . , )




######### SALAMII ###############
###will want to take subsets to align alongside other data
nodes=grep("salami", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 1)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.01, y=20, width=0.001, color='black', offset =.2)+
  ggplot2::xlim(-.0005, 0.05)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73", "#0072B2"))+
  labs(fill="isolation")
p1 %>% insert_left(g)
##saved as A5

##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Polsonii', 'Psalamii'))
#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73", "#0072B2"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation environment ONLY TAKING SALAMII
ggplot(subset(captains_subset, species == "Psalamii"), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73", "#0072B2"))+
  theme_pubr(x.text.angle = 65)+
  labs(x="")

###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.01, y=9.5, width=0.001, color='black', offset =.1)+
  ggplot2::xlim(-.0005, 0.055)


g2=scaleClade(g, 37, .1) %>%
  scaleClade(45, .1) %>%
  collapse(37, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(45, 'mixed', fill="firebrick", alpha=0.75)  +
  geom_cladelabel(node=37, label="P.olsonii", color='black', fontsize=3)+
  geom_cladelabel(node=45, label="P.salamii", color='black', fontsize=3)+
  geom_rootedge(rootedge = -0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73", "#0072B2"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

##for the same plot to add to the side
##same plot as above splitting the isolation origins of solitum
##except changing the x axis label angle and removing the x axis title
captains_subset$species=factor(captains_subset$species, levels=c('Polsonii', 'Psalamii'))
ggplot(subset(captains_subset, species == "Psalamii" & isolation_simple != ""), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73", "#0072B2"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")
##saved as svg 200x350 captains_boxplot.isolation.salamii_only


######### NALGIOVENSE AND CHRYOSGENUM/RUBENS AND DESERTORIUM ###############
###will want to take subsets to align alongside other data
nodes=grep("nalgiovense", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 1)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )

##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))

subset_sp_clades2=subset(subset_sp_clades2, species != "Prubens" )
##rename chrysogenum to chrysogenum/rubens 
subset_sp_clades2$species=gsub('Pchrysogenum', 'Pchrysogenum/Prubens', subset_sp_clades2$species)
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]

subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]

captains_subset$species=gsub('Pchrysogenum', 'Pchrysogenum/XXXXX', captains_subset$species)
captains_subset$species=gsub('Prubens', 'Pchrysogenum/XXXXX', captains_subset$species)
captains_subset$species=gsub('XXXXX', 'Prubens', captains_subset$species)

###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.005, y=20, width=0.001, color='black', offset =.2)+
  ggplot2::xlim(-.0005, 0.03)
##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)


##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Pdesertorum', 'Pflavigenum', 'Pmononematosum', 'Pchrysogenum/Prubens', 'Pnalgiovense'))

#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation environment ONLY TAKING NALGIOVENSE
ggplot(subset(captains_subset, species == "Pnalgiovense"), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 65)+
  labs(x="")


###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.005, y=12, width=0.001, color='black', offset =.2)+
  ggplot2::xlim(-.0005, 0.025)



g2=scaleClade(g, 64, .1) %>%
  scaleClade(97, .1) %>%
  collapse(64, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(97, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=64, label="P.nalgiovense", color='black', fontsize=3)+
  geom_cladelabel(node=97, label="P.chrysogenum/P.rubens", color='black', fontsize=3)+
  geom_cladelabel(node=33, label="P.desertorum", color='black', fontsize=3)+
  geom_cladelabel(node=61, label="P.flavigenum", color='black', fontsize=3)+
  geom_cladelabel(node=62, label="P.mononematosum", color='black', fontsize=3)+
  geom_rootedge(rootedge = -0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

##for the same plot to add to the side
##same plot as above splitting the isolation origins of solitum
##except changing the x axis label angle and removing the x axis title
##split by isolation environment ONLY TAKING NALGIOVENSE
ggplot(subset(captains_subset, species == "Pnalgiovense" & isolation_simple != ""), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")



######### Janthinellum, MYSTERY ###############
###will want to take subsets to align alongside other data
nodes=grep("janthinellum", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 1)
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.01, y=10, width=0.001, color='black', offset =.2)+
  ggplot2::xlim(-.0005, 0.09)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73","#F0E442" ))+
  labs(fill="isolation")
p1 %>% insert_left(g)


##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Pjanthinellum', 'Prolfsii', 'Pochrochloron', 'Psubrubescens', 'Pdaleae', 'Pcataractarum', 'Pbrasilianum'))
#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442" ))+
  theme_pubr()+
  coord_flip()+
  labs(x="")


###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.01, y=5, width=0.001, color='black', offset =.05)+
  ggplot2::xlim(-.0005, 0.09)


g2=scaleClade(g, 28, .1) %>%
  scaleClade(25, .1) %>%
  scaleClade(17, .1) %>%
  collapse(28, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(25, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(17, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=28, label="P.brasilianum", color='black', fontsize=3)+
  geom_cladelabel(node=25, label="P.subrubescens", color='black', fontsize=3)+
  geom_cladelabel(node=17, label="P.janthinellum", color='black', fontsize=3)+
  geom_cladelabel(node=12, label="P.cataractarum", color='black', fontsize=3)+
  geom_cladelabel(node=11, label="P.daleae", color='black', fontsize=3)+
  geom_cladelabel(node=8, label="P.ochrochloron", color='black', fontsize=3)+
  geom_cladelabel(node=7, label="P.rolfsii", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","#F0E442" ))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

######### CANESCENS, MYSTERY ###############
###will want to take subsets to align alongside other data
nodes=grep("canescens", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 2)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))
subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]
###now plot all together
##need to always adjust the fixed values such as geom_hilight extend, treescale position and size and geom_rootedge length
#g=ggtree(subset_tree, size=0.2) +
#  guides(fill = "none")+
#  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extendto=1, colour="black", linetype="dotted", size=0.25)+
#  scale_fill_viridis_b()+
#  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
#  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey" )+
#  geom_treescale(x=0.0001, y=100, width=0.001, color='black', offset = 5)+
#  ggplot2::xlim(-.0005, 0.0052)
##without colours highlighting the clades
g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.01, y=10, width=0.001, color='black', offset =.2)+
  ggplot2::xlim(-.0005, 0.07)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00","#009E73","#56B4E9", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)
##saved as A5 pdf

##now plotting grouped captain counts
##first order according to phylogeny
captains_subset$species=factor(captains_subset$species, levels=c('Pnucicola', 'Pantarcticum', 'Parizonense', 'Pcanescens'))
#ggplot(captains_subset, aes(x=species, y=captains))+
#  geom_boxplot(outlier.shape = "", width=0.3)+
#  geom_jitter(width=0.075, alpha=0.5, aes(colour=isolation_simple))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  theme_pubr(x.text.angle = 65)
ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00","#009E73","#56B4E9", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")


###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.01, y=5.5, width=0.001, color='black', offset =.1)+
  ggplot2::xlim(-.0005, 0.08)


g2=scaleClade(g, 16, .1) %>%
  scaleClade(23, .1) %>%
  scaleClade(25, .1) %>%
  collapse(16, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(23, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(25, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=16, label="P.canescens", color='black', fontsize=3)+
  geom_cladelabel(node=23, label="P.arizonense", color='black', fontsize=3)+
  geom_cladelabel(node=11, label="P.nucicola", color='black', fontsize=3)+
  geom_cladelabel(node=25, label="P.antarcticum", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(captains_subset, aes(x=species, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00","#009E73","#56B4E9", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)


###########MALACHITEUM/HERQUIE AND SCLEROTIORUM, MYSTERY #####
###will want to take subsets to align alongside other data
nodes=grep("malachiteum", rooted_tree$tip.label)
##this can add a grouping label to the full tree
rooted_tree = groupOTU(rooted_tree, nodes)
##now can subset based on the node of the MRCA
##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
clade = MRCA(rooted_tree, nodes)
subset_tree = tree_subset(rooted_tree, clade, levels_back = 2)
##rotate the node of biformer and casei/camemberti
##now redo all the steps for the hilighting of species etc
##get a list of all the species by taking everything before the first dot (they have all been named properly for this sort of purpose)
subset_species=unique(do.call('rbind', strsplit(as.character(subset_tree$tip.label),'.',fixed=TRUE))[,1] )
##create empty dataframe
subset_sp_clades=as.data.frame(c())
for (sp in subset_species) {
  nodes=grep(paste(sp), subset_tree$tip.label)
  clade=MRCA(subset_tree, nodes)
  output=print(paste(sp,clade))
  subset_sp_clades=rbind(subset_sp_clades, output)
}
##rename header temporarily
colnames(subset_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
subset_sp_clades2=subset_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]
##force a strict order for the plotting to maintain it
subset_sp_clades2$species <- factor(subset_sp_clades2$species, levels = unique(subset_sp_clades2$species))
subset_sp_clades2$node <- factor(subset_sp_clades2$node, levels = unique(subset_sp_clades2$node))

subset_sp_clades2=subset(subset_sp_clades2, species != "Pherquei" )
##rename Pmalachiteum to Pmalachiteum/Pherquei 
subset_sp_clades2$species=gsub('Pmalachiteum', 'Pmalachiteum/Pherquei', subset_sp_clades2$species)
subset_species2=subset_sp_clades2[,1]
subset_clades=subset_sp_clades2[,2]

subset_sp_clades3=data.frame(node=as.numeric(paste(subset_clades)), species=subset_species2)
##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]

captains_subset$species=gsub('Pmalachiteum', 'Pmalachiteum/XXXXX', captains_subset$species)
captains_subset$species=gsub('Pherquei', 'Pmalachiteum/XXXXX', captains_subset$species)
captains_subset$species=gsub('XXXXX', 'Pherquei', captains_subset$species)


##extract each column as a list MIGHT NOT NEED THIS
subset_species2=subset_sp_clades3[,2]
subset_clades=subset_sp_clades3[,1]


g=ggtree(subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.04, y=8, width=0.001, color='black', offset = .1)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
genomes=subset_tree$tip.label
captains_subset=captains[grepl(paste(genomes, collapse="|"), captains$genome2), ]

p1=ggplot(captains_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)

##ADD LIKLIHOOD OF BRANCHING 



################ ASPERGILLUS #################

###generating a tree for all the aspergillus genomes
asp_raw_tree <- ape::read.tree("~/cluster2/projects/Penicillium/phylogenetics/ncbi_aspergillus/phylogeny/busco_aa_mafft_trimal.concat.veryfasttree_doubleprecision_treefile")
##plot the raw tree
ggtree(asp_raw_tree)

##rerooting tree
nodes=grep("roqueforti", asp_raw_tree$tip.label)
root=MRCA(asp_raw_tree, nodes)
asp_rooted_tree=ape::root.phylo(asp_raw_tree, node = root-1)
asp_rooted_tree=TreeTools::Preorder(asp_rooted_tree)

##reading in DUF3435 data
asp_captains=read.csv("cluster2/projects/ncbi_penicillium_aspergillus/genome_datasets/all_aspergillus/all_aspergillus.genome_strain.TEMP_STARFISH_RESULTS.PLUS_PHYLOGENY_COMMENTS.PLUS_MYBs.tsv", sep="\t" , header=T)
##remove the rows which contain genomes which were removed due to busco issues
asp_captains2=asp_captains[- grep("emoved", asp_captains$phylogeny_notes),]

###need to install aplot and then use this to align a plot of the cpatins againt the tree
g=ggtree(asp_rooted_tree)
asp_captains2$genome2=paste(asp_captains2$species, asp_captains2$genome, sep=".")
p1=ggplot(asp_captains2, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2))+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73", "purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme(legend.position='none')+
  geom_text(aes(label=species, y=captains+2), size=0.5)
p1 %>% insert_left(g, width=0.5)


##can specify certain clades to subset
##manual mods for old tree without renaming to look at sojae/parasiticus and flavus/oryzae

#asp_raw_tree$tip.label=gsub('Aspergillus.sp.GCA002894705', 'Aspergillus.parasiticus.GCA002894705', asp_raw_tree$tip.label)
#asp_captains2$genome2=gsub('Aspergillus.sp.GCA002894705', 'Aspergillus.parasiticus.GCA002894705', asp_captains2$genome2)


###niger and it'a naming convention is very complicated WILL HELP LATER
###because of this I will recalssify some strains
###primarily niger strains reclassified within the tubingensis and welwitshiae clades
###below will just rename them in our subset phylogeny and captains data
tubingensis_renamed=c('23625475','25769395','25769005','25769155','25769305','23625375','23625455','25769325','23625415','23625315','23625435','25768885','25769035','23625395','25769265','25768975','25769175','25769125','23509725','25769085','25769245','01515345','13618955','25769415','26284195','25769065','23625355','25769365','25769215')
welwit_renamed=c('25769515','23618435','27923945','15586215','27923805','23973445','27923985','27924005','15586235','27923745','25769485','19843515','27923785')
for (geno in tubingensis_renamed) {
  asp_rooted_tree$tip.label=gsub(paste("Aspergillus.niger.GCA0", geno, sep=""), paste("Aspergillus.tubingensis.GCA0", geno, sep=""), asp_rooted_tree$tip.label)
  asp_captains2$genome2=gsub(paste("Aspergillus.niger.GCA0", geno, sep=""), paste("Aspergillus.tubingensis.GCA0", geno, sep=""), asp_captains2$genome2)
}
for (geno in welwit_renamed) {
  asp_rooted_tree$tip.label=gsub(paste("Aspergillus.niger.GCA0", geno, sep=""), paste("Aspergillus.welwitschiae.GCA0", geno, sep=""), asp_rooted_tree$tip.label)
  asp_captains2$genome2=gsub(paste("Aspergillus.niger.GCA0", geno, sep=""), paste("Aspergillus.welwitschiae.GCA0", geno, sep=""), asp_captains2$genome2)
}

###can also just plot the captains all together split by isolation, therfore ignoring the phylogeny effects
##this helps for clinical samples which are often lowly sampled
ggplot(subset(asp_captains2, isolation_simple != "" & isolation_simple != "unknown" & isolation_simple != "environment-food"), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442", "#0072B2"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")
##saved as svg 300x600 "Aspergillus_no_phylogeny.captains_isolation.svg"

##simple plot in the end but removing the Penicillium outgroup
pensp=c("Pcaseifulvum.ESE00019","P.roqueforti.LCP06133","Proqueforti.LCP06133" ,"Proqueforti.LCP06136","Pfuscoglaucum.ESE00090" ,"Pbiforme.LCP05531","Psolitum.LCP06249")
asp_rooted_trim_tree=drop.tip(asp_rooted_tree, pensp)

ggtree(asp_rooted_trim_tree, size=0.1)+
  geom_tiplab(size=0.3)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)

asp_species= unique(do.call('rbind', strsplit(as.character(asp_rooted_trim_tree$tip.label),'.',fixed=TRUE))[,2] )
##remove the ".sp" category
asp_species=asp_species[asp_species %in% "sp" == FALSE]   
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_rooted_trim_tree$tip.label)
  clade=MRCA(asp_rooted_trim_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]
asp_sp_clades2$species <- sub("^", "A", asp_sp_clades2$species )
asp_sp_clades3$species <- sub("^", "A", asp_sp_clades3$species )

##plot entire phylogeny with boxes around the species groups and legible genome names at tips
ggtree(asp_rooted_trim_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node, fill=node), colour="black", linetype="dotted", size=0.25 ,alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_tiplab(size=.2)+
  scale_fill_viridis_b()+
  geom_rootedge(rootedge = 0.01, linewidth=0.2)+
  geom_treescale(x=0.1, y=150, width=0.05, color='black')
###saved as A4 pdf 'phylogeny_complete.aspergillus.genome_lab.pdf'

##tree labelling species groups 
ggtree(asp_rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )
###saved as A4 pdf 'phylogeny_complete.aspergillus.pdf'

##same as above but with node support
ggtree(asp_rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )+
  geom_text2(aes(label=label, subset=!isTip), hjust=-.2, size=0.5)
###saved as A4 pdf 'phylogeny_complete.aspergillus.support.pdf'

##plot nice tree with barplot associated
g=ggtree(asp_rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )

genomes=asp_rooted_trim_tree$tip.label
asp_captains3=asp_captains2[grepl(paste(genomes, collapse="|"), asp_captains2$genome2), ]
p1=ggplot(asp_captains3, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))
p1 %>% insert_left(g, width=2)
##save as A4 pdf "Aspergillus_phylogeny.captains_isolation.pdf"

##CAN ALSO SAVE THE ABOVE TREE FOR ANOTHER PLOT WITH ASPERGILLUS ALONGSIDE
ASPERGILLUSsupport=ggtree(asp_rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )+
  geom_text2(aes(label=label, subset=!isTip), hjust=-.2, size=0.5)
ASPERGILLUS=ggtree(asp_rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" )

##alternative plot with a smaller bar graph
p1 %>% insert_left(g, width=5)

##and saving as another i,age to join with aspergillus
ASPERGILLUS2support=p1 %>% insert_left(ASPERGILLUSsupport, width=5)
ASPERGILLUS2=p1 %>% insert_left(ASPERGILLUS, width=5)

#####NOW PLOT ASPERGILLUS AND PENICILLIUM TREES TOGETHER 
ggarrange(PENICILLIUM, ASPERGILLUS)
##now the same but with the bar plots alongside
##first need to convert the aplot images (used to make sure the barplot aligned with the phylogeny tips) to a grob
PENICILLIUM2grob=as.grob(PENICILLIUM2)
PENICILLIUM2supportgrob=as.grob(PENICILLIUM2support)
ASPERGILLUS2grob=as.grob(ASPERGILLUS2)
ASPERGILLUS2supportgrob=as.grob(ASPERGILLUS2support)
##now plot with and without support (can manually combine later only keeping nodes for species branches and larger)
ggarrange(PENICILLIUM2grob, ASPERGILLUS2grob)
#saved as A4 pdf "Aspergillus_phylogeny_AND_Penicillium_phylogeny_public.captains_isolation.pdf"
ggarrange(PENICILLIUM2supportgrob, ASPERGILLUS2supportgrob)
#saved as A4 pdf "Aspergillus_phylogeny_AND_Penicillium_phylogeny_public.captains_isolation.support.pdf"

###############################################
##now select for the species of interest
##first doing sojae and have to incoperate some outpgroups too
asp_nodes=grep("Aspergillus.sojae", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 2)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]

##now plot with captains
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.0002, y=14, width=0.001, color='black', offset = 0.25)+
  ggplot2::xlim(-.0005, 0.0035)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]

p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)


##now plotting grouped captain counts
##first order according to phylogeny
asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
asp_captains2_subset$species2=factor(asp_captains2_subset$species2, levels=c('A.parasiticus', 'A.sojae'))
ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

###trying to generate a cleaner image for presentation by collapsing nodes and only showing the boxplots
##this tree will have no labels initially, only manually done by by after collapsing etc
g=ggtree(asp_subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.00005, y=7.7, width=0.001, color='black', offset = 0.05)

g2=scaleClade(g, 24, .1) %>%
  scaleClade(30, .1) %>%
  scaleClade(19, .1) %>%
  collapse(24, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(30, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(19, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=24, label="A.sojae", color='black', fontsize=3)+
  geom_cladelabel(node=30, label="A.parasiticus", color='black', fontsize=3)+
  geom_cladelabel(node=19, label="A.parasiticus", color='black', fontsize=3)+
  geom_rootedge(rootedge = -0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)



###########################################
###can also run this for oryzae species
asp_nodes=grep("Aspergillus.oryzae", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 8)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]


##now plot with captains
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.0002, y=14, width=0.001, color='black', offset = 0.25)+
  ggplot2::xlim(-.0005, 0.011)

##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]

p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)

##plot the boxplots for each species coloured by the isolation origin
asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
asp_captains2_subset$species2=factor(asp_captains2_subset$species2, levels=c('A.texensis', 'A.minisclerotigenes', 'A.flavus', 'A.oryzae'))
ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation, using only oryzae 
##also remove env-clin and env-food with only one sample each
ggplot(subset(asp_captains2_subset,  species2 == "A.oryzae" & isolation_simple != "" & isolation_simple != "unknown" & isolation_simple != "environment-clinical" & isolation_simple != "environment-food"  ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73", "#0072B2"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")
##saved as svg 200x350 captains_boxplot.isolation.oryz_only


  
###generate the reduced clean phylogeny next to the boxplots for species
g=ggtree(asp_subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.001, y=55, width=0.001, color='black', offset = 0.2)

g2=scaleClade(g, 459, .1) %>%
  scaleClade(336, .1) %>%
  scaleClade(442, .1) %>%
  scaleClade(580, .1) %>%
  scaleClade(626, .1) %>%
  scaleClade(648, .1) %>%
  scaleClade(647, .1) %>%
  scaleClade(345, .1) %>%
  collapse(459, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(336, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(442, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(580, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(626, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(648, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(647, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(345, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=459, label="A.oryzae", color='black', fontsize=3)+
  geom_cladelabel(node=336, label="A.minisclerotigenes", color='black', fontsize=3)+
  geom_cladelabel(node=442, label="A.flavus", color='black', fontsize=3)+
  geom_cladelabel(node=580, label="A.flavus", color='black', fontsize=3)+
  geom_cladelabel(node=626, label="A.flavus", color='black', fontsize=3)+
  geom_cladelabel(node=648, label="A.flavus", color='black', fontsize=3)+
  geom_cladelabel(node=647, label="A.flavus", color='black', fontsize=3)+
  geom_cladelabel(node=345, label="A.flavus", color='black', fontsize=3)+
  geom_cladelabel(node=334, label="A.flavus", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)


###########################################
###can also run this for terreus species
asp_nodes=grep("Aspergillus.terreus", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 5)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]


g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% subset(asp_sp_clades3, species != "sp")$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=subset(asp_sp_clades3, species != "sp")  , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.08, y=13, width=0.001, color='black', offset = 0.25)+
  ggplot2::xlim(-.0005, 0.3)+
  geom_cladelabel(node=3, label="A.sp.GCA019721355", color='black', fontsize=2)+
  geom_cladelabel(node=4, label="A.sp.GCA001931935", color='black', fontsize=2)


##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]

p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","purple","#F0E442"))+
  labs(fill="isolation")
p1 %>% insert_left(g)

##plot the boxplots for each species coloured by the isolation origin
asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
asp_captains2_subset$species2=factor(asp_captains2_subset$species2, levels=c('A.neotritici','A.taichungensis','A.campestris','A.candidus','A.nanangensis','A.olivimuriae','A.sp','A.terreus'))
ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation
##just terreus
ggplot(subset(asp_captains2_subset, species2 == "A.terreus" & isolation_simple != "" & isolation_simple != "unknown" ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")


###generate the reduced clean phylogeny next to the boxplots for species
g=ggtree(asp_subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.1, y=6, width=0.001, color='black', offset = 0.2)+
  ggplot2::xlim(-.005, 0.3)

g2=scaleClade(g, 46, .1) %>%
  scaleClade(31, .1) %>%
  collapse(46, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(31, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=46, label="A.neotritici", color='black', fontsize=3)+
  geom_cladelabel(node=31, label="A.terreus", color='black', fontsize=3)+
  geom_cladelabel(node=1, label="A.nanangensis", color='black', fontsize=3)+
  geom_cladelabel(node=2, label="A.olivimuriae", color='black', fontsize=3)+
  geom_cladelabel(node=23, label="A.taichungensis", color='black', fontsize=3)+
  geom_cladelabel(node=24, label="A.campestris", color='black', fontsize=3)+
  geom_cladelabel(node=25, label="A.candidus", color='black', fontsize=3)+
  geom_cladelabel(node=3, label="A.sp", color='black', fontsize=3)+
  geom_cladelabel(node=4, label="A.sp", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.005, linewidth=0.5)

##same plot as from above
p2=ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)

##For terreus this is only one actual clinical strain left after filtering
##the other two high ones are from saline environments, but are actually the same strain (one a mutant of the other)
## the env-clin is from dalian/human feces as usual


###########################################
###can also run this for niger species
asp_nodes=grep("Aspergillus.niger", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 1)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]


##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]




##now plot with captains
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.001, y=60, width=0.001, color='black', offset = 0.5)+
  ggplot2::xlim(-.0005, 0.008)

p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)

##plot the boxplots for each species coloured by the isolation origin
asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
asp_captains2_subset$species2=factor(asp_captains2_subset$species2, levels=c('A.welwitschiae', 'A.niger'))
ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation
ggplot(subset(asp_captains2_subset, isolation_simple != "" & isolation_simple != "unknown" ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 35)+
  facet_grid(species2 ~ .)+
  labs(x="")

##just plot A.niger and remove the env-clin stupid samples
ggplot(subset(asp_captains2_subset, isolation_simple != "" & isolation_simple != "unknown" & isolation_simple != "environment-clinical"  & species2 == "A.niger" ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")
##saved as svg 250x450 captains_boxplot.isolation.niger_only

###generate the reduced clean phylogeny next to the boxplots for species
g=ggtree(asp_subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.0005, y=15, width=0.001, color='black', offset = 0.2)

g2=scaleClade(g, 93, .1) %>%
  scaleClade(164, .1) %>%
  collapse(93, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(164, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=93, label="A.niger", color='black', fontsize=3)+
  geom_cladelabel(node=164, label="A.welwitchiae", color='black', fontsize=3)+
  geom_rootedge(rootedge = -0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

##splitting up all the strains by isolation_simple for niger, welwitshiae and tubingensis which are vey mixed
##split by isolation
ggplot(subset(asp_captains2_subset, isolation_simple != "" & isolation_simple != "unknown"  ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 35)+
  facet_grid(species2 ~ .)+
  labs(x="")

###########################################
###can also run this for niger species
asp_nodes=grep("Aspergillus.luchuensis", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 2)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]


##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]




##now plot with captains
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.002, y=30, width=0.001, color='black', offset = 0.5)+
  ggplot2::xlim(-.0005, 0.018)

p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","#F0E442", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)

##plot the boxplots for each species coloured by the isolation origin
asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
asp_captains2_subset$species2=factor(asp_captains2_subset$species2, levels=c('A.eucalypticola', 'A.luchuensis', 'A.tubingensis'))
ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation
#ggplot(subset(asp_captains2_subset, isolation_simple != "" & isolation_simple != "unknown" ), aes(x=isolation_simple, y=captains))+
#  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
#  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
#  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
#  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442", "#0072B2", "#000000"))+
#  theme_pubr(x.text.angle = 35)+
#  facet_grid(species2 ~ .)+
#  labs(x="")

###generate the reduced clean phylogeny next to the boxplots for species
g=ggtree(asp_subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.002, y=10, width=0.001, color='black', offset = 0.5)
  
g2=scaleClade(g, 68, .1) %>%
  scaleClade(117, .1) %>%
  collapse(68, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(117, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=68, label="A.tubingensis", color='black', fontsize=3)+
  geom_cladelabel(node=117, label="A.luchuensis", color='black', fontsize=3)+
  geom_cladelabel(node=66, label="A.eucalypticola", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)




###########################################
###can also run this for fumigatus species
asp_nodes=grep("Aspergillus.fumigatus", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 1)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]


##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]




##now plot with captains
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.008, y=50, width=0.001, color='black', offset = 0.5)+
  ggplot2::xlim(-.0005, 0.03)


p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)

##same as above but zooming in on fumigatus with a bit of a xlim hack
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.008, y=50, width=0.001, color='black', offset = 0.5)+
  ggplot2::xlim(0.0225, 0.0275)


p1=ggplot(subset(asp_captains2_subset, species != "Aspergillus.oerlinghausenensis"), aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)
##saved as A5


##plot the boxplots for each species coloured by the isolation origin
asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
asp_captains2_subset$species2=factor(asp_captains2_subset$species2, levels=c('A.oerlinghausenensis', 'A.fumigatus'))
ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation
ggplot(subset(asp_captains2_subset, species2 == "A.fumigatus" & isolation_simple != "" & isolation_simple != "unknown"   ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 35)+
  facet_grid(species2 ~ .)+
  labs(x="")

###generate the reduced clean phylogeny next to the boxplots for species
g=ggtree(asp_subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.005, y=30, width=0.001, color='black', offset = 0.2)

g2=scaleClade(g, 339, .1) %>%
  collapse(339, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=339, label="A.fumigatus", color='black', fontsize=3)+
  geom_cladelabel(node=124, label="A.oerlinghausenensis", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

##only compare environment and clinical samples and the env-cli samples are basically all clonal in a single clade
ggplot(subset(asp_captains2_subset, species2 == "A.fumigatus" & isolation_simple != "environment-clinical" & isolation_simple != "" & isolation_simple != "unknown" & species2 != "A.minisclerotigenes"  ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")


####there is alos a group that produced many genomes within a couple of projects that were all treated the same
###therefore this controls for between project differences
###a seperate file with thee genomes has been manually generated
fum_cap_subset=read.csv(file="cluster2/projects/ncbi_penicillium_aspergillus/genome_datasets/all_aspergillus/fumigatus_captains.two_bioprojects.tsv", header=T, sep='\t')
ggplot(fum_cap_subset, aes(x=environment, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=environment), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")
##saved as svg 250x450 'captains_boxplot.isolation.fumi_only.svg'



###########################################
###can also run this for niger species
asp_nodes=grep("Aspergillus.udagawae", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 1)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]


##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]


##now plot with captains
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.005, y=4, width=0.001, color='black', offset = 0.1)+
  ggplot2::xlim(-.0009, 0.024)


p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73", "purple","#F0E442"))+
  labs(fill="isolation")
p1 %>% insert_left(g)

##plot the boxplots for each species coloured by the isolation origin
asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
asp_captains2_subset$species2=factor(asp_captains2_subset$species2, levels=c('A.pseudoviridinutans', 'A.felis', 'A.udagawae'))
ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

##split by isolation
ggplot(subset(asp_captains2_subset, species2 == "A.udagawae" & isolation_simple != "" & isolation_simple != "unknown" ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442"))+
  theme_pubr(x.text.angle = 35)+
  facet_grid(species2 ~ .)+
  labs(x="")

###generate the reduced clean phylogeny next to the boxplots for species
g=ggtree(asp_subset_tree, size=0.5) +
  guides(fill = "none")+
  geom_treescale(x=0.01, y=4, width=0.001, color='black', offset = 0.1)+
  ggplot2::xlim(-.0005, 0.025)
  

g2=scaleClade(g, 12, .1) %>%
  scaleClade(18, .1) %>%
  collapse(12, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(18, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=12, label="A.udagawae", color='black', fontsize=3)+
  geom_cladelabel(node=18, label="A.felis", color='black', fontsize=3)+
  geom_cladelabel(node=7, label="A.pseudoviridinutans", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
p2=ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)
##saved as pdf, 3x11.69 inches, landscape (so short and wide)

##no need to split by isolation, already very clear



##### CRISTATUS/CHEVALIERI/MONTEVIDENSIS
asp_nodes=grep("cristatus", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 4)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]


##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]


##now plot with captains
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.1, y=7, width=0.01, color='black', offset = 0.1)+
  ggplot2::xlim(-.0005, 0.3)


p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)

###generate the reduced clean phylogeny next to the boxplots for species
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_treescale(x=0.1, y=10, width=0.01, color='black', offset = 0.1)+
  ggplot2::xlim(-.0005, 0.3)


g2=scaleClade(g, 35, .1) %>%
  scaleClade(38, .1) %>%
  scaleClade(42, .1) %>%
  scaleClade(45, .1) %>%
  scaleClade(52, .1) %>%
  scaleClade(60, .1) %>%
  collapse(35, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(38, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(42, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(45, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(52, 'mixed', fill="firebrick", alpha=0.75) %>%
  collapse(60, 'mixed', fill="firebrick", alpha=0.75) +
  geom_cladelabel(node=35, label="A.ruber", color='black', fontsize=3)+
  geom_cladelabel(node=38, label="A.brunneus", color='black', fontsize=3)+
  geom_cladelabel(node=42, label="A.chevalieri", color='black', fontsize=3)+
  geom_cladelabel(node=45, label="A.cristatus", color='black', fontsize=3)+
  geom_cladelabel(node=52, label="A.montevidensis", color='black', fontsize=3)+
  geom_cladelabel(node=60, label="A.wentii", color='black', fontsize=3)+
  geom_cladelabel(node=4, label="A.glaucus", color='black', fontsize=3)+
  geom_cladelabel(node=19, label="A.intermedius", color='black', fontsize=3)+
  geom_rootedge(rootedge = 0.0005, linewidth=0.5)

##same plot as from above
asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
asp_captains2_subset$species2=factor(asp_captains2_subset$species2, levels=c('A.wentii', 'A.ruber', 'A.glaucus', 'A.brunneus', 'A.intermedius','A.montevidensis','A.chevalieri','A.cristatus'))
p2=ggplot(asp_captains2_subset, aes(x=species2, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#009E73","purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme_pubr()+
  coord_flip()+
  labs(x="")

multiplot(g2, p2, ncol = 2)

## CHEVALIERI/CRISTATUS/clade including amstelodami, montevidensis, intermedius and many un-identified sp.
##    the clade is interesting, compared to the outgroup brunneus and ruber they have ALOT of starships
##    however ruber and brunneus have a lot of starships too relative to other strains
##      all the brunneus strains are space-craft associated and ruber are env-clin and env-saline...
##    the further outgroup wentii have MUCH less than all and therefore suggest they all beyong wentii have MANY
##    TO NOTE; montevidensis is has also been found in food-production environments commonly https://www.jstage.jst.go.jp/article/jgam/advpub/0/advpub_2019.09.003/_pdf
##    TO NOTE: montevidensis and chevaleri also commonly found in clinical samples : https://academic.oup.com/mmy/article/56/5/541/4372452  and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5573881/



#######sydowii/versicolor
asp_nodes=grep("Aspergillus.versicolor", asp_rooted_tree$tip.label)
asp_rooted_tree = groupOTU(asp_rooted_tree, asp_nodes)

##give the nodes based on the grep above using species names
##can also replace the variable 'nodes' with just a node number 
asp_clade = MRCA(asp_rooted_tree, asp_nodes)
asp_subset_tree = tree_subset(asp_rooted_tree, asp_clade, levels_back = 2)

##labelling the tree with the species extracted
asp_species= unique(do.call('rbind', strsplit(as.character(asp_subset_tree$tip.label),'.',fixed=TRUE))[,2] )
##create empty dataframe
asp_sp_clades=as.data.frame(c())
for (sp in asp_species) {
  nodes=grep(paste("Aspergillus",sp,sep="."), asp_subset_tree$tip.label)
  clade=MRCA(asp_subset_tree, nodes)
  output=print(paste(sp,clade))
  asp_sp_clades=rbind(asp_sp_clades, output)
}
##rename header temporarily
colnames(asp_sp_clades) = "temp"
##split the column into two seperate columns with appropriate headers
asp_sp_clades2=asp_sp_clades %>% tidyr::separate(temp, c('species', 'node'))
##force a strict order for the plotting to maintain it
asp_sp_clades2$species <- factor(asp_sp_clades2$species, levels = unique(asp_sp_clades2$species))
asp_sp_clades2$node <- factor(asp_sp_clades2$node, levels = unique(asp_sp_clades2$node))
##extract each column as a list MIGHT NOT NEED THIS
asp_species2=asp_sp_clades2[,1]
asp_clades=asp_sp_clades2[,2]
asp_sp_clades3=data.frame(node=as.numeric(paste(asp_clades)), species=asp_species2)

asp_sp_clades2$species2=paste("Aspergillus",asp_sp_clades2$species, sep=".")
asp_species3=asp_sp_clades2[,3]

##now plot with captains
g=ggtree(asp_subset_tree, size=0.2) +
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node),alpha = 0.1, to.bottom=T, linetype="dotted", colour="firebrick")+
  geom_rootedge(rootedge = 0.0005, linewidth=0.2)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 2, barcolour="grey", angle=270, vjust=-0.5, hjust=0.5 )+
  geom_treescale(x=0.01, y=25, width=0.001, color='black', offset = 0.25)+
  ggplot2::xlim(-.0005, 0.04)
##extract name of genomes in tree
##then subset captains data frame for only those genomes
asp_genomes=asp_subset_tree$tip.label
asp_captains2_subset=asp_captains2[grepl(paste(asp_genomes, collapse="|"), asp_captains2$genome2), ]

p1=ggplot(asp_captains2_subset, aes(genome2, captains))+
  geom_col(aes(fill=isolation_simple, group=genome2), show.legend = F)+
  coord_flip()+
  theme_tree2()+
  scale_fill_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#000000"))+
  labs(fill="isolation")
p1 %>% insert_left(g)


asp_captains2_subset$species2=factor(do.call('rbind', strsplit(as.character(asp_captains2_subset$genome2),'.',fixed=TRUE))[,2], asp_species2)
asp_captains2_subset$species2=factor(paste("A", asp_captains2_subset$species2, sep = "."))
##split by isolation but combine all the species together due to low counts
##probably will not use
ggplot(subset(asp_captains2_subset, isolation_simple != "" & isolation_simple != "unknown" ), aes(x=isolation_simple, y=captains))+
  gghalves::geom_half_boxplot(side = "r", outlier.shape = "", width=0.4, position=position_nudge(x = .1))+
  ggpubr::geom_pwc( label='p.adj.signif', p.adjust.method = "bonferroni", vjust = .5, hide.ns = T)+
  geom_jitter(width=0.05, alpha=0.6, aes(colour=isolation_simple), show.legend = FALSE)+
  scale_colour_manual(values = c("#E69F00", "#009E73","purple","#F0E442","#000000"))+
  theme_pubr(x.text.angle = 35)+
  labs(x="")

####sydowii, creber and versicolor are made up of 56/75 of the Dalian Medical human feces samples.
###very hard to conclude anything using this so no more analysis
###the versicolor strain with the largest number of starships is env-clinical but not from dalian
###two sydowii strains with high starships are isolated from the roots of Phytolacca americana 
### This plant has high heavy metal tolerance and lives in extreme metal concentration soils
##the large number of sydowii strains from healthy fecal samples have high diversity


####not much to say for nidulans/unguis, some clinical samples with high starship count but very few and few samples in total
###calidoustus and close species are all environment or NA
###for uvarum, japonicus, fijensis has very few samples in total








####the phylogeny can also be used to indicate the diversity in gene clusters present across the two genera
clusters=read.csv(file="~/cluster2/projects/Penicillium/genome_graphs/cluster_analyses/population_BLAST/combined.blastn_all_95.species_strain.tsv", sep='\t', header=T)

##calculation proportion of each isolation type per cluster but only those that contain all the 'core' genes
clusters2=subset(clusters, core/count==1) %>%
  group_by(cluster, isolation) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
##plot it
ggplot(data=clusters2)+geom_col(aes(x=cluster, y=freq, fill=isolation), position = position_stack())+
  theme_pubr(x.text.angle = 35)+
  scale_fill_manual(values = c("#E69F00","#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))



###now we can try plot a heatmap next to the phylogeny
##first change the name of the genome column to match the 
clusters$genome2=clusters$genome
##get only those that have the full cluster
clusters3=subset(clusters, core==count)

##need to now fill the absences with NAs
##first get list of all genomes in tree and put into a dataframe column
genomelist=rooted_trim_tree_public$tip.label
genomelist2=asp_rooted_trim_tree$tip.label
genomespen=data.frame(unlist(genomelist))
names(genomespen)[1] <- "genome2"
genomesasp=data.frame(unlist(genomelist2))
names(genomesasp)[1] <- "genome2"
genomes=rbind(genomespen,genomesasp)
names(genomes)[1] <- "genome2"
genomes1=genomes
genomes1$cluster="lactose"
genomes2=genomes
genomes2$cluster="dityrosine"
genomes3=genomes
genomes3$cluster="salt"
genomes4=genomes
genomes4$cluster="arsenic"
genomes5=genomes
genomes5$cluster="ethanol"
genomes=rbind(genomes1,genomes2,genomes3,genomes4,genomes5)
##now merge this file with our cluster data suing the genome2 column name
clusters4=merge(clusters3,genomes,by = c("genome2","cluster"), all = TRUE)


##set the heatmap as one too
heatmap=ggplot(clusters4, aes(x=cluster, y=genome2))+
  geom_tile(aes(fill=isolation))+theme_pubr(x.text.angle = 25,  border  = FALSE, base_size = 5, legend = "none")+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y = element_blank())+
  xlab(NULL)+
  ylab(NULL)+
  scale_fill_manual(values = c("#E69F00","#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"),  na.value = "white")

##plot the trees
ASPERGILLUS=ggtree(asp_rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=asp_sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange" , fontface = "italic")
PENICILLIUM=ggtree(rooted_trim_tree_public, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades2$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades3 , mapping=aes(node=node, label=species), fontsize = 1, barcolour="orange", fontface = "italic")

##plot together with the tree variables already being used
heatmap %>% insert_left(PENICILLIUM, width = 5)
heatmap %>% insert_left(ASPERGILLUS, width = 5)


##however this can be very cluttered especially with all the species names etc
##therefore can we reduce the number of species names highlighted on the phylogeny to just those that contain at least one cluster
sp_clades4=filter(sp_clades3, species %in% clusters3$species)
asp_sp_clades4=filter(asp_sp_clades3, species %in% clusters3$species)
##rerun the trees
ASPERGILLUS=ggtree(asp_rooted_trim_tree, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% asp_sp_clades4$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=asp_sp_clades4 , mapping=aes(node=node, label=species), fontsize = 1.25, barcolour="orange" , fontface = "italic")
PENICILLIUM=ggtree(rooted_trim_tree_public, size=0.1)+
  geom_rootedge(rootedge = 0.02, linewidth=0.2)+
  geom_treescale(x=0.15, y=200, width=0.01, color='black', offset = 5)+
  guides(fill = "none")+
  geom_hilight(mapping=aes(subset=node %in% sp_clades4$node, fill=node),alpha = 0.5, to.bottom=T, extend=0.01)+
  geom_cladelab(data=sp_clades4 , mapping=aes(node=node, label=species), fontsize = 1.25, barcolour="orange", fontface = "italic")

##replot with the heatmap
PENCLUST=heatmap %>% insert_left(PENICILLIUM, width = 5)
ASPCLUST=heatmap %>% insert_left(ASPERGILLUS, width = 5)

##now convert the two phylogeny-heatmaps to variables and then to grobs
##first need to convert the aplot images (used to make sure the barplot aligned with the phylogeny tips) to a grob
PENCLUSTgrob=as.grob(PENCLUST)
ASPCLUSTgrob=as.grob(ASPCLUST)
##now plot together
ggarrange(PENCLUSTgrob, ASPCLUSTgrob)
#saved as 8.27x5 pdf landscape (A4 width made square) "phylogeny_cluster.8x5.pdf"

###recalculate proportions of isolation types using only the public data again
public=c(asp_rooted_trim_tree$tip.label , rooted_trim_tree_public$tip.label)
clusters5=filter(clusters3, genome2 %in% public)
##calculation proportion of each isolation type per cluster but only those that contain all the 'core' genes
clusters6=clusters5 %>%
  group_by(cluster, isolation) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
##plot it
ggplot(data=clusters6)+geom_col(aes(x=cluster, y=freq, fill=isolation), position = position_stack())+
  theme_pubr(x.text.angle = 35)+
  scale_fill_manual(values = c("#E69F00","#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))

##need to now generate a stat relative to the number of each type in total
##can combine the two captains files, then calculate it based on the isolation column
all_captains=full_join(captains2_public, asp_captains3)
prop_isolation=all_captains %>%
  group_by(isolation_simple) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
##plot it

ggplot(data=prop_isolation)+geom_col(aes(x="All genomes", y=freq, fill=isolation_simple), position = position_stack())+
  theme_pubr(x.text.angle = 35)+
  scale_fill_manual(values = c("#E69F00","#009E73", "purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  xlab(label = element_blank())

##plot both side by side for comparison
prop=ggplot(data=clusters6)+geom_col(aes(x=cluster, y=freq, fill=isolation), position = position_stack())+
  theme_pubr(x.text.angle = 35)+
  scale_fill_manual(values = c("#E69F00","#009E73","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  theme(legend.position = "none")+
  xlab(label = element_blank())
prop2=ggplot(data=prop_isolation)+geom_col(aes(x="All genomes", y=freq, fill=isolation_simple), position = position_stack())+
  theme_pubr(x.text.angle = 35)+
  scale_fill_manual(values = c("#E69F00","#009E73", "purple","#F0E442","#56B4E9", "#0072B2", "#000000"))+
  xlab(label = element_blank())
ggarrange(prop, prop2, common.legend = T, align = "h", widths = c(2,1), legend="right")
##saved as A5 pdf "clusters_proportion_isolation.pdf"



















######to show horizontal transfer we can test for whether similarities in starship genes are more similar than should be expected at the whole genome level
##to do this diamond blastp was run on all vs all for captains and a set of 50 random BUSCOs
## for the captains, only matches above 95% coverage, both ref and query, and 95% identity were kept
##then for each match, the BUSCOs were used to estimate the whole genomes average amino acid identity

##now we can read in the summary file and plot this
HGT=read.csv(file="~/cluster2/projects/Penicillium/mash_distances/captains.diamond.summary4.tsv", header=T, sep='\t')

##get the difference of the identities
HGT$difference=HGT$captain_identity-HGT$genome_AAI

##plot the points, highlithing those that have a difference greater than 5%
ggplot2::ggplot()+
  geom_point(data=HGT, aes(x=genome_AAI, y=captain_identity), alpha=0.5)+
  theme_pubr()+
  labs(x="Average Amino Acid Identity (AAI)", y="Pairwise Captain Identity")+
  geom_point(data=subset(HGT, difference > 5), aes(x=genome_AAI, y=captain_identity), alpha=0.5, colour="red")
##or colouring by difference
ggplot2::ggplot()+
  geom_point(data=HGT, aes(x=genome_AAI, y=captain_identity, colour=difference))+
  theme_pubr()+
  labs(x="AAI (%)", y="Pairwise Captain Identity (%)", colour="Captain identity-AAI (%)")+
  scale_colour_material("orange")
##save as svg 1200x300 'captain_and_genome_AAI' TO DO!!

##highlighting specific geno;es
ggplot2::ggplot()+
  geom_point(data=HGT, aes(x=genome_AAI, y=captain_identity), alpha=0.5)+
  theme_pubr()+
  labs(x="Average Amino Acid identity (AAI)", y="Pairwise Captain identity")+
  geom_point(data=subset(HGT, genome_ref == "Pantarcticum.IBT31339_ASM2897420v1" & genome_query == "Ppalitans.F6-4S-1A-F_GCA_019191055.1"), aes(x=genome_AAI, y=captain_identity), colour="red")



ggplot2::ggplot()+geom_jitter(data=HGT, aes(x=difference, y=""), alpha=0.5, width = 0, height=0.5)+theme_pubr()+labs(x="Pairwise Captain identity - AAI", y="")






##the pairs with the largest difference between whole genome AAI and captain identity were further evaluated
##the whole genomes were aligned and only well aligned regions were extracted and realigned
##now we can use the minimap2 paf alignment to generate some plots
##need the library "SVbyEyeé (has many dependencies that need BioCManager installation (will tell you the list if it doesn't install))
#devtools::install_github("daewoooo/SVbyEye", branch="master")
library(SVbyEye)

##first comparison
comp1=readPaf(paf.file="cluster2/projects/Penicillium/mash_distances/comp1.minimap2.contig_ROI.paf", include.paf.tags =T, restrict.paf.tags = "cg")
plotMiro(paf.table=comp1, color.by = 'identity', perc.identity.breaks = c(90, 95, 99))
plotAVA(paf.table=comp1, color.by = 'identity', perc.identity.breaks = c(90, 95, 99), binsize = 1000)
##add annotation of captain using a captain bed file
comp1plot=plotMiro(paf.table=comp1, color.by = 'identity', perc.identity.breaks = c(90, 95, 99))
comp1_tyr<- read.table("cluster2/projects/Penicillium/mash_distances/comp1.tyr.bed", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
comp1_tyr$seqnames=comp1_tyr$contig
comp1_tyr$strand=comp1_tyr$sense
comp1_tyr.gr <- GenomicRanges::makeGRangesFromDataFrame(comp1_tyr, ignore.strand = FALSE)
addAnnotation(ggplot.obj = comp1plot, annot.gr = comp1_tyr.gr, coordinate.space = 'self', y.label.id = 'seqnames', annotation.level = 0)



##2nd
comp2=readPaf(paf.file="cluster2/projects/Penicillium/mash_distances/comp2.minimap2.contig_ROI.paf", include.paf.tags =T, restrict.paf.tags = "cg")
plotMiro(paf.table=comp2, color.by = 'identity', perc.identity.breaks = c(90, 95, 99))
plotAVA(paf.table=comp2, color.by = 'identity', perc.identity.breaks = c(90, 95, 99), binsize = 1000)
##add annotation of captain using a captain bed file
comp2plot=plotMiro(paf.table=comp2, color.by = 'identity', perc.identity.breaks = c(90, 95, 99))
comp2_tyr<- read.table("cluster2/projects/Penicillium/mash_distances/comp2.tyr.bed", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
comp2_tyr$seqnames=comp2_tyr$contig
comp2_tyr$strand=comp2_tyr$sense
comp2_tyr.gr <- GenomicRanges::makeGRangesFromDataFrame(comp2_tyr, ignore.strand = FALSE)
addAnnotation(ggplot.obj = comp2plot, annot.gr = comp2_tyr.gr, coordinate.space = 'self', y.label.id = 'seqnames', annotation.level = 0)
