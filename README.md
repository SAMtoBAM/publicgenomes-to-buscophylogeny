# publicgenomes-to-buscophylogeny
Downloading sets of public genomes, analysing BUSCOs and generating trees for such genomes.   
This pipeline is basically just a simple way to genrate a gene-based phylogeny for a set of genomes already on ncbi  
This pipeline was used for O'Donnell et al. 2023 (DOI)  


###Example below is using the genus Aspergillus which contains ~1200 genomes as of 2023 (particularly flavus, oryzae and fumigatus)

####PUT AUTOMATIC DOWNLOADING STEP  
####RUNNING THIS ON ENTIRE GENUS GENOMES DOWNLOADED FROM NCBI  
####AFTER DOWNLOADING ALL GENOMES EXCLUDING THOSE THAT ARE MARKED AS CONCERNING THE RAW GENOME FILES ARE RENAMED TO KEEP JUST THE GCF/GCA NUMBER AND THE UNIQUE CONTIG/SCAFFOLD NUMBER GIVEN BY NCBI/SRA  
###JUST MOVE INTO THE 'DATA' REGION WITH ALL THE GENOMES AND RUN THE BELOW  

    mkdir ../genomes_renamed  
    ls **/*.fna | while read genome
    do
    assembly=$( echo $genome | awk -F "/" '{print $NF }' | awk -F "." '{print $1}' | sed 's/GCF_/GCF/g' | sed 's/GCA_/GCA/g'  )
    cat $genome | awk -v assembly="$assembly" '{if($0 ~ ">" ) {print ">"assembly"_"$1} else {print}}' | sed 's/_>/_/g' | sed 's/NW_/NW/g' | sed 's/NC_/NC/g' | gzip > ../data_renamed/${assembly}.fa.gz
    done  

####ALL THE COMPLETE DATASET FILES WERE THEN NAMED AS 'all_${genus}' AND PLACED IN A FOLDER CALLLED 'genome_datasets'  
####ACTUALLY IN ~/projects/ncbi_aspergillus_penicillium/starfish/  
###LASTLY TOOK THE FILE **/ncbi_dataset/data/data_summary.tsv AND TRANSFORMED IT TO A SIMPLER FORMAT WITH THE GCA TAG AND THE STRAIN NAME.  
###here we will just use the entire aspergillus dataset and run BUSCO and a quick veryfasttree analysis  
###THEN NEED TO CHECK IF THE NAMING IS ALL CORRECT!  

    ###first run busco
    genus="aspergillus"
    cd ~/projects/Penicillium/phylogenetics/
    mkdir ncbi_${genus} && cd ncbi_${genus}
    mkdir genomes
    mkdir busco_analyses
    #######MORE COMMENTS FOR CREATING THE ENVIRONMENT (including adding ,muscle, trimal and mafft)
    conda activate busco
    ##needed to install muscle into the busco conda env and expecially the v5 version due to changes in the output formats
    #mamba install -c bioconda muscle

    ##get all the genomes except for those outside of the penicillium and aspergillus genus such as fusarium
    ##the aspergillus genomes are kept as outgroups (only four of them)
    ##for mucor use the busco library '-l fungi_odb10'
    ls ../../../ncbi_penicillium_aspergillus/genome_datasets/all_${genus}/ncbi_dataset/data_renamed/* | while read genome
    do
    genome2=$( echo $genome | sed "s/.fa.gz//g" | awk -F "/" '{print $NF}'  )
    species=$( cat ../../../ncbi_penicillium_aspergillus/genome_datasets/all_${genus}/all_${genus}.genome_strain.tsv | awk -v genome2="$genome2" '{if($3 ~ genome2"." ) print $1}' )
    if [ ! -d "busco_analyses/${species}.${genome2}" ]
    then
    cp $genome genomes/${species}.${genome2}.fa.gz
    zcat genomes/${species}.${genome2}.fa.gz > ${species}.${genome2}.fa
    busco -m genome -i ${species}.${genome2}.fa -o busco_analyses/${species}.${genome2} -l eurotiales_odb10 -c ${threads}
    mv busco_analyses/${species}.${genome2}/run_eurotiales_odb10/* busco_analyses/${species}.${genome2}/
    rm ${species}.${genome2}.fa
    fi
    done

    ##also add a subset of Penicillium genomes for rooting
    for strain in ESE00019 ESE00090 LCP06249 LCP05531 LCP06133 LCP06136
    do
    cp -r ../../references/busco_analyses/*.${strain} ./busco_analyses/
    done


    ##now we need to check if some genomes have low BUSCO scores and therefore will be problematic
    ##we will put a threshold on 80% BUSCO score and remove genomes below this
    ls busco_analyses/ | while read folder
    do
    cat busco_analyses/$folder/short_summary.txt | grep "Complete and single" | awk -v folder="$folder" '{print folder"\t"$1"\t"$1/4191}'
    done | awk '{if($3 < .8) print $1}' | while read genome
    do
    rm -r busco_analyses/$genome
    rm genomes/${genome}.fa.gz
    done


    mkdir phylogeny
    cd phylogeny
    cat ../busco_analyses/**/full_table.tsv | grep -v "^#" | awk '$2=="Complete" {print $1}' | sort | uniq -c | awk '{print $2"\t"$1}' > busco_complete.tsv
    total=$( ls ../busco_analyses/ | wc -l )
    cat busco_complete.tsv | awk -v total="$total" '{if($2 > (total*0.985)) print $1}' > busco_complete.in_99p.tsv
    ##using this threshold we have 3258/4190 (~75%)
    total2=$( cat busco_complete.in_99p.tsv | wc -l  )
    ls ../busco_analyses/ | while read genome
    do
    cat busco_complete.in_99p.tsv | while read protein
    do
    ls ../busco_analyses/$genome/busco_sequences/single_copy_busco_sequences/ | grep "^${protein}.faa"
    done | wc -l | awk -v genome="$genome" -v total2="$total2" '{ print genome"\t"$1"\t"$1/total2}'
    done > busco_complete.in_99p.per_genome.tsv


    ##now we take this list, extract the protein name and set up a directory where we place all the sequences for the genomes where they were found
    ##during the files need to be renamed using a genome/strain identifier
    mkdir busco_aa
    ##get protein
    cat busco_complete.in_99p.tsv | while read protein
    do
    ##get list of genomes which contain a single complete copy
    rm  busco_aa/${protein}.concatenated.faa
    ls ../busco_analyses/**/busco_sequences/single_copy_busco_sequences/${protein}.faa | while read path
    do
    ##get a new genome name to add as a prefix to the protein file
    genome=$( echo $path | awk -F "/"  '{print $3}' )
    cat ../busco_analyses/${genome}/busco_sequences/single_copy_busco_sequences/${protein}.faa | awk -v genome="$genome" '{if($0 ~ ">") {print ">"genome} else {print}}' >> busco_aa/${protein}.concatenated.faa
    done
    done

    threads="10"
    ##we can also just use the proteins
    ##therefore need to feed these concatenated protein multi-fasta files individualy to MAFFT for multi-sequence alignment
    ##then we need to trim the alingments
    ##the busco environment has them installed otherwise use: mamba create -n mafft_trimal mafft trimal
    conda activate busco
    ##using default options here
    mkdir busco_aa_mafft
    mkdir busco_aa_mafft_trimal
    ls busco_aa/* | while read file
    do
    name=$( echo $file | awk -F "/" '{print $NF}' | sed 's/.faa//g' )
    mafft --auto --thread ${threads} $file > busco_aa_mafft/$name.mafft
    trimal -in busco_aa_mafft/$name.mafft -out busco_aa_mafft_trimal/$name.mafft.trimal -automated1
    done
    conda deactivate

    #######MORE COMMENTS FOR CREATING THE ENVIRONMENT
    conda activate veryfasttree
    ../../catfasta2phyml/catfasta2phyml.pl -c busco_aa_mafft_trimal/* > busco_aa_mafft_trimal.concat.phyml
    ../../catfasta2phyml/catfasta2phyml.pl -c -f busco_aa_mafft_trimal/* > busco_aa_mafft_trimal.concat.fa

    ##now can run veryfasttree
    ######MORE COMMENTS FOR THE OPTIONS
    ######NEED TO TRY ADDING MORE THREADS -threads 20
    VeryFastTree -double-precision -threads 5 busco_aa_mafft_trimal.concat.fa > busco_aa_mafft_trimal.concat.veryfasttree_treefile 


