# Varlock

## introduction
The human genome is sequence of more than 3 billion DNA bases: adenine, cytozine, thymine and guanine. In technical therms it is a word over an alphabet of {A, C, T, G}. Any two individuals share almost all of the genome and each genome is unique to person.

Modern technologies provides us with efficient tools for reading DNA bases. Genome analysis starts by sequencing DNA sample, which is biochemical process of reading DNA fragments base by base. These fragments are called reads and they need to be aligned (mapped) to reference genome for biological intepretation. Aligned reads can be stored in BAM file format which stands for Binary Alignment Map. Aligned reads reveal where person's genome differs from the reference genome. Single difference is called genomic variant or simply variant. There exists separate file format derived from BAM for storing variants - VCF (Variant Call Format).

Most common type of variant is SNV - Single Nucleotide Variant, a variation in a single DNA base that occurs at a specific position in the genome. Different variants found at single position of the gene represent different forms of this gene, which are called alleles.

Presence of specific alleles in person's genome can tell about his health, disease predisposition, population of origin, physical traits, etc. These properties make the genomic data valuable source of research and diagnosis in future health care. Ongoing research in molecular biology makes genomic data extraction and further interpretation more affordable and complex. While access to meaningful data can provide better understanding of the studied phenomenon, it also provides space for potentional abuse because alleles can be also used for person identification.
 
  In most variant studies researchers need access only short regions of specific genes but genome studies are not limited only to variants. Some of them analyses genome as whole focusing on coverage, that means comparing numbers of mapped reads to specific regions of the genome, which can be useful e.g. for chromosome aneuploidy detection. 

  Due to mentioned properties of genomic data, BAM files need to be stored in secure form while preserving general information such as coverage. Nevertheless the authority managing genomic data must be able to provide access to unaltered region of BAM on valid request from the researcher. These contradicting requirements, when considering that common size of BAM is in order of gigabytes, are not easy to fullfil. To overcome this challenge we developed Varlock - a command line interface utility for securing BAM files in manner that enables distributing access to specific BAM regions among users.

## methods
Varlock is command line interface utility for depersonalization of BAM file based on allele counts and storing differences between original and depersonalized BAM in encrypted form. It uses RSA asymetric encryption so users must provide their own private / public key pairs. Varlock introduces two new binary file types for internal purposes. The first is VAC - Variant Allele Count, which can be thinked of as compact VCF file, the second is BDIFF - BAM difference file.

### VAC
Compact binary file format for storing SNV allele counts. VAC file does not contain header. Each record represents one SNV from VCF file and contains it's absolute position in the genome and four numbers representing counts of four different alleles (A, C, T, G) for that position. Varlock has command for converting VCF file to VAC file. 

### BDIFF
SNV differences between original and depersonalized BAM file are stored in binary BDIFF file format. This format specifies a header consisting of checksum of original BAM and range that BDIFF covers as two absolute positions on the genome. BDIFF record represents difference at one position in the genome. Because BAM can contain multiple alleles for one position, the difference is a mapping from depersonalized alleles to original ones. Record consists of absolute position on the genome and id of allele mapping, which is permutation of DNA bases. 

### encrypt
Process involving BAM depersonalization we simply call encryption. User must supply original BAM file, VAC file and his private key for DIFF encryption.


### decrypt

### reencrypt
