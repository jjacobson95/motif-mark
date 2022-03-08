# Motif-mark  
  
Create Motif Mark png images using Cairo.  
  
  
Command line arguments:   
-f   =  fasta file(s)    
-m   =  motif file(s)   
  
Example Usage:  
  
**motif-mark-oop.py -f Figure_1.fasta Figure_2.fasta -m Fig_1_motifs.txt Fig_2_motifs.txt**  
  
This example will generate two motif mark figures.    
  
  
Limitations:   
- Up to 18 different motifs may be included, 10 or less is recommended.   
- Any number of genes per file is accepted.  
- Any number of figures may be generated.  
- Number of Fasta files must equal the number of motif files.  
  
  

