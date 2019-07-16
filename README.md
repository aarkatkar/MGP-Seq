# MGP-Seq
A python application for predicting DNA, RNA, and amino acid sequences for genes in *Mus musculus.*

## Getting Started

### Prerequisites
+ Download and install Python (<a href="https://www.python.org/downloads/">https://www.python.org/downloads/</a>) to your computer. MGP-Seq is only compatible with Python versions 3.X. 
+ Download *MGP-Seq.py* to any location on your computer

#### Additional Prequisites for Python 3.6+ AND Mac OSX
You will have to upgrade *certifi.* There are two options for doing this:
1. Open the terminal and enter ```pip install --upgrade certifi```
2. Run ```Applications/Python 3.6/Install Certificates.command``` either through the terminal or by browsing to the file and double-clicking it


### Running *MGP-Seq.py*
Python scripts can be run using the steps described <a href="http://www.cs.bu.edu/courses/cs108/guides/runpython.html">here.</a>

## Predicting Sequences
After completing the following steps, click *Impute* to output the predicted sequence(s). An additional C57BL_6J reference is output each time.

### Gene Sequences
1. Search for your desired gene in the search box at the top of the window and hit Enter
2. Select the output directory. All downloads/exports will direct here
3. Download/Select the Reference Gene, SNPs, and Indels through their respective drop-downs
4. Select the desired strains in the panel on the right
5. If you want to exclude SNPs/Indels that are marked as "Low Confidence," check *High Confidence Only*

<p align="center">
  <img src="./resources/Gene.JPG" width="60%">
</p>

### mRNA Sequences
Complete all of the steps required to predict a gene sequence, in addition to the following:

6. Download/Select the desired Reference Transcript (Alignment file)
7. Specify an exon seed length for the alignment algorithm. The length should always be a positive integer and is optimally between 10-14. The default length = 12 will work for most cases. If an alignment fails, nudging the seed length to 11 or 13 might fix the problem

<p align="center">
  <img src="./resources/mRNA.JPG" width="60%">
</p>

### Protein Sequences
Complete all of the steps required to predict an mRNA sequence, in addition to the following:

8. Select which ORF will be translated in the output. *Minimize 5'-UTR* selects the ORF that starts at the very first start codon and ends at the following stop codon. *Maximize ORF* examines the ORF beginning with each potential start codon and translates the longest one. 

<p align="center">
  <img src="./resources/Protein.JPG" width="60%">
</p>

## Working with Local Files
You might want to use files stored locally on your computer rather than files queried from a database. The process is largely the same, but here are some helpful tips:

+ Under each drop-down menu, the "\*Open File Picker\*" option allows you to select local files instead of queried ones. 
+ When the you hit *Impute*, whatever text is written in the search box is attached to the end of output filenames. By typing the name of your gene into the search box prior to imputation, you can still allow your outputs to include the gene in their filenames

## Specifications for homemade files
### SNP/Indel Files
+ Homemade SNP/Indel files should follow the following format:

| Chr | Position  | Gene  | dbSNP       | Ref | PWK_PhJ |
| --- | --------  | ----  | ----------- | --- | ------  |
| 12  | 35498379  | Ahr   | rs241381828 | T   | TA      |
| 12  | 35502300  | Ahr   | rs255095624 | AAC | A       |

The first line indicates that at position 35498379, an insertion replaces T in the reference sequence with TA in the PWK_PhJ strain. The second line indicates that at position 35502300, a deletion replaces ACC with A in the same fashion. <i>Position</i> refers to the position in the chromosome of the first base in <i>Ref</i>. The file should be formatted as a Comma Separated Variable file with UTF-8 encoding. The *Gene Location* is subtracted from the Position column to determine the 0-based index in the gene that will be imputed. For example, if the SNP file says that a SNP occurs in position 101, and 100 is entered into the *Gene Location* box, the SNP will occur at the 2nd nucleotide of the reference gene (or the 1st 0-based index).

### Reference Gene/Reference Transcript
+ The file may contain a FASTA header beginning with a ">" symbol. Newline characters are acceptable. 
+ MGP-Seq is only designed to work with the base symbols A, C, G, and T. Other IUPAC nucleotide symbols are not supported, including U (which must be changed to T for the code to work).
+ MGP-Seq can parse Reference Gene files in a specific format to automatically identify the starting position of the gene. Homemade Reference Gene files might not follow this format and <b>you will have to type in the the position of the first base into the *Gene Location* box</b>. Alternatively, you can edit your file to follow the header format by specifying the gene sequence range after a colon. An example is given below.
```
>AHR:35497979-35534989
ATTTGGTGGTTAGTTTATTATCATATATCACTGTTTACATTTCCTTTGATATTTAAGGACATTATAAATCATACTT
TAACAAAAATCTATGCATGGCATATGTAAGCAAGTAACAGTTTAAAGGACAGCAAGGAGCTATGTAGACCTAAGAT
CTATTGATACTCAAACTAGCTTGTGCATAGTAACAATTGAAAAGTTGATACACTGGCATATTGGTTAAGTGTTCAG
GTGCTGTAAGAAACACCTAGAAAAGAAGAACTTCAATCAGACATACACAAGAGTCTTAGAATTTCATCTAGAAGTT
AATCCCTCAGTCTAGTGGTGCTCTAAGTGGTCTTTCCATTTGGCACTTGGCAATACAGACAAACAAATAGGCA...
```

