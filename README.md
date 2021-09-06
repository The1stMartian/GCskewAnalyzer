### C. Breuer (Formerly C. Merrikh) 9/2021

### Usage:
The GC skew analysis script is written in Python3 and has no dependencies. 
To execute the program, enter:<br>
> python GCS_v6.py

![GC Skew Analysis Script](https://github.com/The1stMartian/GCskew/blob/main/Data/cmdLine.png)

### Function:
The script calculates GC skew values using multiple calculations: <br><br>
     By gene: all nucleotides ("whole gene") or only nucleotides in codon positions 1,2 or 3<br>
     By sliding window: across the full genome using a sliding window of adjustable length<br>

### Outputs (6):
         
1) A flat file with the original gene annotation data plus four columns, one for each of the various whole-gene GC skew values (all nucleotides, or only nucleotides in codon position 1, codon position 2, or codon position 3).
2) Wiggle file of GC skew values calculated over a Sliding window across the whole genome
3) Wiggle file by Whole gene:  All nucleotides
4) Wiggle file by Whole gene:  Codon position 1
5) Wiggle file by Whole gene:  Codon position 2
6) Wiggle file by Whole gene:  Codon position 3

Example output visualized in [MochiView](http://www.johnsonlab.ucsf.edu/mochi/)
![GC skew data visualization](https://github.com/The1stMartian/GCskew/blob/main/Data/Mochi1.png)

### Notes:
a) For whole gene GC skew values, the value is relative to the <i>leading strand</i> of the replication fork.<br><br>
b) This is where ter information is used - nucleotide 1 of the fasta sequence is assumed to be part of the origin, and ter is given
by the user. Note: ori and ter coordinates can be looked up on the DoriC website: http://tubic.org/doric/public/index.php<br><br>
(In cases where ori is <i>not</i> the first nucleotide, or very close, the genome sequence needs to be manually rotated.)<br><br>
c) As a result, whole-gene GC skew values post-ter will have the opposite sign (-/+) as the value indicated by the sliding GC skew value which is <i>not</i> inverted post-ter.<br><br>
d) Discussion of script function and hyperparameter adjustment can be found at the bottom, and changed in the analyzeGCskew() function<br>

### Development/References:

I developed this script as part of my investigations into the evolutionary history of gene inversions:

     Gene inversion potentiates bacterial evolvability and virulence
        Christopher N Merrikh, Houra Merrikh
        Nature Communications, Nov. 2018.
        https://pubmed.ncbi.nlm.nih.gov/30405125/

I also used the updated the script with the codon-specific calculations for a response to a matters arising paper that challenged our initial findings:<br><br>
[Challenge](https://www.biorxiv.org/content/10.1101/2020.01.14.906818v1)<br>
[Our Response](https://www.biorxiv.org/content/10.1101/2020.05.26.117366v2)<br>

### Discussion
Our analysis of the codon-specific GC skew values suggests that the replication fork and gene inversions, together, explain the vast majority of GC skew values in bacterial genomes. As a result, one can interpret negative GC skew values as an indication that a gene is in an atypical orientation as the result of an inversion event. 

Rationale: Codon position 3 (CP3) nucleotides are highly mutable because they can often be changed without altering the encoded amino acid. Meanwhile, CP1 nucleotides are under strong selection because mutations usually result in a change in the encoded amino acid. Therefore, we did a CP1 vs. CP3 based analysis. We hypothesized that if the  replication fork is the primary driver of the GC skew, CP3-based GC skew values should generally be positive with respect to the leading strand of the replication fork due to their general mutability. Indeed that is what our results showed. This argues against the hypothesis that sequence context and other forms of selection are major drivers of the GC skew. 

Additionally, we hypothesized that all genes are capable of gaining a positive GC skew after millions of years of evolution, irrespective of the codon position in question. This was supported by our observation that even CP1 nucleotides, which are under strong selection, show the overwhelming pattern of positive GC skew values. By extension, negative GC skew values (over whole gene regions) should be caused exclusively (more or less) by recombination events that placed these positive-GC skew genes on the opposite strand relative to their typical orientation that was retained throughout millions of years of evolutionary history. Mathematically, this inverts the GC skew value (x * -1) because it exchanges the number of Gs and Cs in the top strand. If our hypothesis is correct, then the negative GC skew values of inverted genes should be well retained in CP1 nucleotides (due to selection on protein function) and less well retained in CP3 nucleotides (which can mutate without effect). Also, CP3-based GC skew values of inverted genes should always be roughly equal or <i>higher</i> than the skew of CP1 nucleotides. This prediction is based on our fundamental assumption (supported with a variety of data) that normal DNA replication is the primary driver of the GC skew, driving it in a positive direction, even after taking into account other influences. Thus, the initially negative GC skew of inverted genes should be increasingly positive, an effect that should be far more potent in CP3 nucleotides than CP1 nucleotides. This is exactly what we observed. 

The implications of this effect are significant as it implies that nearly all bacterial genes have been encoded primarily on the leading strand of each chromosome arm throughout nearly all of evolutionary history. Furthermore, it is only through fairly recent (in evolutionary terms at least - i.e. the past 100M years) inversion events that lagging strand genes have been produced. It demonstrates that negative GC skew regions can be reasonably interpreted as an indication that the region is currently in an inverted configuration. 

Interestingly, we also identified evidence that in some cases periods of lagging strand encoding can accelerate the evolution of genes due to head-on replication-transcription conflicts. This could help bacteria evolve and adapt to new environments. This includes evolving new antibiotic resistance or immune evasion properties. As such, this our results have clinical relevance.  

I encourage others to test out these patterns in their species of interest. 