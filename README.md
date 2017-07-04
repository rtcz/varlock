# Varlock

## Introduction
The human genome is a sequence of more than 3 billion DNA bases: adenine, cytosine, thymine, and guanine. In technical terms, it is a word over an alphabet {A, C, T, G}. Any two human individuals share over 98% of all genome, but at the same time each genome is unique to each person.

Modern technologies provide efficient tools for determining the precise order of nucleotides within a DNA molecule - _DNA sequencing_. Usually, genome analysis starts by sequencing fragments from a DNA sample. These fragments are called _reads_ and they need to be aligned (mapped) to a reference genome for biological interpretation. Aligned reads are typically stored in a (Binary) Sequence Alignment Map file format (_SAM_ or _BAM_). These alignments reveal differences between the studied and the reference genome. A single difference is called a _genomic variant_ or simply a _variant_. A separate file format called _Variant Call Format_ (_VCF_) was derived from SAM for storing genomic variants.

The most common type of variant is a _Single Nucleotide Variant_ (_SNV_), a variation in a single DNA base that occurs at a specific position in the genome. Different variants found at a single position of a gene represent different forms of this gene called _alleles_.

Presence of specific alleles in a person's genome can tell about his health, predisposition to specific diseases, population of origin, physical traits, etc. These properties make the genomic data very valuable for research and diagnosis in a future health care. Genomic data extraction and interpretation is more affordable and complex every year as research in these fields progresses. While access to meaningful data is necessary for better understanding of human biology, it also provides space for potential abuse, since alleles can be very effectively used for person identification. 
 
Most variant studies require access only to short regions of specific genes and some analyses do not even require information about genomic variants at all. For example, chromosome aneuploidy detection needs only information about number of reads per chromosome or chromosome parts - information that cannot be abused for person identification. 

Due to this potential threat to privacy, BAM files need to be stored in a secure form while preserving general information such as position of reads. These depersonalized BAM files can be then made available to general public without fear of personal data leakage. Moreover, the authority managing genomic data must be able to provide access to unaltered regions of BAM on a valid request from an authorized user. These requirements are not easy to fulfill effectively, since the common size of BAM is in order of gigabytes. To overcome this challenge we developed Varlock - a command line interface utility for securing BAM files in a manner that enables distributing access to specific BAM regions among users.

## Methods
Varlock is a command line interface utility for depersonalization of BAM files by mutating bases at SNV positions from provided VCF file. The differences between the original and the depersonalized BAM are stored in encrypted form, so only authorized users can recover the original BAM. Varlock uses RSA asymmetric encryption so users must provide their own private/public key pairs. Two new binary file types are introduced for internal purposes: the first is _Variant Allele Count_ (_VAC_), which is essentially a compact VCF file. The second is _BAM difference_ (_BDIFF_), where the differences between the original and depersonalized BAM are stored. Varlock provides commands for encryption, decryption, granting access, and conversion VCF file to VAC file. 

### Variant Allele Count (VAC) file format
VAC is a compact binary file format for storing SNV allele counts. VAC file does not contain header. Each record represents one SNV from a VCF file and contains its absolute position in the genome and four numbers representing counts of different alleles (A, C, T, G) for that position. 

### BAM Difference (BDIFF) file format 
SNV differences between original and depersonalized BAM file are stored in a binary BDIFF file format. The format specifies a header consisting of the checksum of depersonalized BAM and an effective range that BDIFF covers as two absolute positions in the genome. A single BDIFF record represents a difference at one position in the genome. Since BAM can contain multiple alleles for one position, the difference is a map from depersonalized alleles to original ones. A record consists of an absolute position on the genome and a permutation of DNA bases at that position.  

### Encrypt
Process involving BAM depersonalization is simply called _encryption_. User must supply the original BAM file, VAC (or VCF) file, and his public key for encryption. Both VAC and BAM files are iterated at the same time to produce a new depersonalized BAM. Every time a genomic position of VAC record overlaps with  one or more alignments from BAM file, a permutation for this alignment column is created. 

The permutation depersonalizes original bases and it is created in an iterative manner. At each step a base chosen randomly from a multinomial distribution of allele counts defined by a VAC record is mapped to the most abundant base from the alignment column. Only remaining bases are considered in next steps. This effectively makes this permutation a representative of the SNV distribution in population.  

Bases in alignment column at position of SNV are permuted with the corresponding permutation. When the generated permutation is not identity, a new record is written to BDIFF file with the position of SNV and generated permutation. After all alignments or VAC records are processed, depersonalization algorithm has ended and remaining alignments are written. Checksum of original BAM file and VAC file is added to the header of the depersonalized BAM for validation and file management purposes. At the same time, a checksum of depersonalized BAM is added to the header of BDIFF file. BDIFF file format requires also a genomic range, the full range is used for BAM files encrypted for the first time to allow access to all data.

Last step of encryption is the actual encryption of BDIFF file. Public key can not be used for the encryption directly because the RSA encryption has a limited input size. Instead, the AES encryption with a randomly generated key is used to encrypt the whole file. The AES key is encrypted with the private key and stored as a part of encrypted BDIFF file. 
 
 After the encryption, the original BAM file can be deleted and the access to data is then restricted to the owner of the private key paired with the public key used for encryption. 

### Decrypt
Decrypt is the inverse command of encrypt. Whole depersonalized BAM or its parts can be restored to their original form when a BDIFF file and a corresponding private key is provided. The BAM is restored in the range defined by BDIFF or only part of this range is restored if requested by the user. The algorithm for restoring BAM first decrypts the AES key, which is then used to decrypt the rest of the BDIFF file. The original BAM is then restored similarly as it is depersonalized, only the permutation is not created randomly but taken from the corresponding BDIFF record as an inverse permutation. 

### Grant access (reencrypt)
Owner who encrypted the original BAM can use a public key of other user to grant him access to the BAM or its part. Owner must provide his private key to decrypt BDIFF file and the public key for the second encryption. The depersonalized BAM file is then sliced and only the desired part is decrypted and reencrypted. Only the authorized user can then use this new reencrypted BDIFF file (together with the depersonalized BAM) to decrypt the parts of the original BAM. This is used to provide the access only to a part of the original BAM to an authorized user, while preserving security of the rest of the file. 

## Discussion
A secure pseudo-random number generator or a true random number generator must be used in multinomial drawing and to generate AES keys to ensure security. Other security risk lies in an incomplete VCF file, since the rate of depersonalization is given by the provided VCF file, i.e., only the variants from the VCF file will be masked. The user is thus responsible for providing a VCF file covering all individual related SNVs. For purpose of depersonalization, a VCF with large population sample is needed, where the allele counts represent general distribution in population.

## Conclusion
Varlock is a tool suitable for storing BAM files in a depersonalized form with the option to restore their original content. The authority managing BAM files can distribute access to specific genomic regions of original BAM file among users, while other parts of the original BAM file stay inaccessible. Now, only the sequence is depersonalized, however, we intend to extend this concept in the future to other BAM file format fields such as CIGAR string, base quality, and optional fields.

# Development
## Testing


run all tests
python3 -m unittest discover

run single test
python3 -m unittest tests.<FILE_NAME>.<CLASS_NAME>.<METHOD_NAME>

## Examples
run example
python3 -m examples.bam_mutator

profile example
python3 -m cProfile -s tottime examples/bam_mutator.py | head -n 140


## Samtools
samtools view -H examples/resources/sample.mut.bam
samtools view -h resources/out.mut.bam >> out.mut.sam




