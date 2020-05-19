### Description for usage of the scripts

\# create full length viruses from contigs

`for f in \*/contigs.fa; do echo ${f%%/contigs.fa}; python create_full_length_virus.py ${f} SARS_CoV2.fa` 
`${f%%contigs.fa}/quast/contigs_reports/nucmer_output/contigs.all_snps`
`${f%%contigs.fa}/quast/contigs_reports/nucmer_output/contigs.coords ${f%%contigs.fa}/quast/report.tsv ${f%\*/contigs.fa}.bam` 
`all_strains/; done &> genome_replace_new.log`

\# add root >NC\_045512.2|NC\_045512.2|2020-01-31 with sequence

\# run virusprepro (-> GitHub)

`conda activate virustracker`

`bash ~/Documents/VirusTracker/phylogeo_discrete/preprocessing/viralPrePro.sh strains_cds.fa "acr"`

`~/Documents/VirusTracker/phylogeo_discrete/preprocessing/software/ "NC_045512.2"`

\# copy samples.tsv file and create tree (optional)

`Rscript visualize_tree.Rscript RAxML_bestTree.acr_raxml.phy SARS-CoV-2-tree.pdf`

\# Root the tree by putting NC in extra bracket (manually, replace acr with strains in name)  
\# run raxML -f A  

`raxmlHPC -f A -t RAxML_bestTree.strains_raxml.phy -s acr_cds_mapped_c.aln -m GTRCAT -n "strains_raxml.phy"`

\# add labels (edge lengths)

`python3 add_labels.py RAxML_bestTree.strains_raxml.phy RAxML_nodeLabelledRootedTree.strains_raxml.phy`

\# take strains produced by this (RAxML_marginalAncestralStates.strains_raxml.phy)  
\# add all (leaf) strains to original mapping (acr_cds_mapped_c.aln -> strains_cds.aln)

`cat acr_cds_mapped_c.aln > strains_cds.aln`

`while read -r line; do if [[ $line == \>* ]]; then echo -e "\n$line"; else echo -ne ${line^^}; fi; done < strains_cds.aln > strains_cds.aln1`

`echo -e "\n" >> strains_cds.aln1`

`while read -r line; do if [[ $line == \>* ]]; then echo -ne "$line "; else echo ${line}; fi; done < strains_cds.aln1 > strains_cds.aln`

`rm -f strains_cds.aln1`

\# remove leading ">"

`sed -i '/^>/s/^.//' strains_cds.aln`

`cat RAxML_marginalAncestralStates.strains_raxml.phy >> strains_cds.aln`

\# replace ? by -

`sed -i 's:?:-:g' strains_cds.aln`

\# add "label location" to first line  
\# run count_links.R script

`Rscript count_links.R RAxML_nodeLabelledRootedTree.strains_raxml.phy strains_cds.aln > changes.tsv`

`sed -i 's:\[1\] ::g; s:"::g; s:$`::g; s:`::g' changes.tsv`

\#remove all empty changes
