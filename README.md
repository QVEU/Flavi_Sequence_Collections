# Zika Research Consortium (ZiRCon) Flavivirus Structural Protein Expression Library 

## Library Generation
All full-length coding sequences were retrieved from the Viral Pathogen Resource (ViPR), now part of the Bacterial and Viral Bioinformatics Resource Center (BV-BRC), for all available full-length flavivirus genome sequences (as of February 2021[1]). Sequences were codon aligned in DECIPHER 2.0[2] and the structural protein (CprME) regions were extracted from this alignment. CprME sequences were then aligned and clustered using SeqCluster.py (https://github.com/QVEU/QVEU_Code/blob/main/Alignment/SeqCluster_v2.py), an in-house script that uses the BIRCH clustering algorithm [3], implemented in scikit-learn [4], to choose a specified number of sequences with balanced representation based on pairwise Levenshtein distances [5]. Clustering was performed and representative sequences collected for each of three sets: all flaviviruses, ZIKV, and all DENV serotypes and a non-redundant combination of these collections were combined with curated sequences representing lab and vaccine strains to make up the final library. For completeness, all available non-redundant ZIKV sequences available in the database were included in the collection. A total of 372 CprME sequences with stop codons added, were ordered for synthesis from Twist Biosciences: 173 ZIKV, 121 DENV and 78 representing other flaviviruses. Each ORF was cloned into the NotI and EcoRI sites of a pcDNA3-1(-) vector. The input data, alignments, and information about the sequences are available at: https://github.com/QVEU/Flavi_Sequence_Collections/tree/main.

## References
1.	Olson, R. D. et al. Introducing the Bacterial and Viral Bioinformatics Resource Center (BV-BRC): a resource combining PATRIC, IRD and ViPR. Nucleic Acids Res. 51, D678–D689 (2023).
2.	Wright, E. Using DECIPHER v2.0 to analyze big biological sequence data in R. R J. 8, 352 (2016).
3.	Zhang, T., Ramakrishnan, R. & Livny, M. BIRCH. in Proceedings of the 1996 ACM SIGMOD international conference on Management of data - SIGMOD ’96 (ACM Press, New York, New York, USA, 1996). doi:10.1145/233269.233324.
4.	Pedregosa, F. et al. Scikit-learn: Machine Learning in Python. J. Mach. Learn. Res. 12, 2825–2830 (2011).
5.	Levenshtein, V. I. Binary Codes Capable of Correcting Deletions, Insertions and Reversals. Soviet Physics Doklady 10, 707 (1966).

![image](https://github.com/QVEU/Flavi_Sequence_Collections/assets/10180619/8b5d7173-fd32-4164-b8d5-744933d7247e)

## Contact
Patrick Dolan, NIAID-NIH
Patrick.Dolan@nih.gov
