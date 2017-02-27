# Varlock

## introduction
The human genome is sequence of more than 3 billion DNA bases: adenine, cytozine, thymine and guanine. In technical therms it is a word over an alphabet of {A, C, T, G}. Any two individuals share almost all of the genome and each genome is unique to person.

Modern technologies provides us with efficient tools for reading DNA bases. Genome analysis starts by sequencing DNA sample, which is biochemical process of reading DNA fragments base by base. These fragments are called reads and they need to be aligned (mapped) to reference genome for biological intepretation. Aligned reads can be stored in BAM file format which stands for Binary Alignment Map. Aligned reads reveal where person's genome differs from the reference genome. Single difference is called genomic variant or simply variant. There exists separate file format derived from BAM for storing variants - VCF (Variant Call Format).

Most common type of variant is SNV - Single Nucleotide Variant, a variation in a single DNA base that occurs at a specific position in the genome. Different variants found at single position of the gene represent different forms of this gene, which are called alleles.

Presence of specific alleles in person's genome can tell about his health, disease predisposition, population of origin, physical traits, etc. These properties make the genomic data valuable source of research and diagnosis in future health care. Ongoing research in molecular biology makes genomic data extraction and further interpretation more affordable and complex. While access to meaningful data can provide better understanding of the studied phenomenon, it also provides space for potentional abuse because alleles can be also used for person identification.
 
  In most variant studies researchers need access only short regions of specific genes but genome studies are not limited only to variants. Some of them analyses genome as whole focusing on coverage, that means comparing numbers of mapped reads to specific regions of the genome, which can be useful e.g. for chromosome aneuploidy detection. 

  Due to mentioned properties of genomic data, BAM files need to be stored in secure form while preserving general information such as coverage. Nevertheless the authority managing genomic data must be able to provide access to unaltered region of BAM on valid request from the researcher. These contradicting requirements, when considering that common size of BAM is in order of gigabytes, are not easy to fullfil. To overcome this challenge we developed Varlock - a command line interface utility for securing BAM files in manner that enables distributing access to specific BAM regions among users.

## methods
Varlock is command line interface utility for depersonalization of BAM file by mutating bases at SNV positions and storing the differences between original and depersonalized BAM in encrypted form. It uses RSA asymetric encryption so users must provide their own private / public key pairs. Varlock introduces two new binary file types for internal purposes. The first is VAC - Variant Allele Count, which can be thinked of as compact VCF file, the second is BDIFF - BAM difference file. Varlock provides encryption command for depersonalization, its counterpart - decryption, reencryption and command to convert VCF file to VAC file.

### VAC
Compact binary file format for storing SNV allele counts. VAC file does not contain header. Each record represents one SNV from VCF file and contains it's absolute position in the genome and four numbers representing counts of four different alleles (A, C, T, G) for that position. Varlock provides command for converting VCF file to VAC file.

### BDIFF
SNV differences between original and depersonalized BAM file are stored in binary BDIFF file format. This format specifies a header consisting of checksum of original BAM and effective range that BDIFF covers as two absolute positions on the genome. BDIFF record represents difference at one position in the genome. Because BAM can contain multiple alleles for one position, the difference is a map from depersonalized alleles to original ones. Record consists of absolute position on the genome and id of DNA bases permutation. The permutation represent values of ordered mutation map with keys: A, C, T, G. 

### encrypt
Process involving BAM depersonalization we simply call encryption. User must supply the original BAM file, VAC file and his public key for BDIFF encryption. Both VAC and BAM files are iterated at the same time to produce a new mutated BAM. When genomic position of VAC record overlaps with alignments from BAM file, mutation map for this alignment column is created. 

Mutation map maps original bases to new ones and is created in iterative manner. At each step a base is drawed randomly from multinomial distribution of allele counts defined by VAC record. This is the resulting mutated base. At the same time most abundant base from the alignment column is taken as an original base. In the end of the step mapping from original to mutated base is added to mutation map. In next step only remaining bases from allele counts and alignment column respectively are considered. 

Bases in alignment column at position of SNV are mutated with mutation map. When non-synonymous mutation occurs new record is written to BDIFF file as id of bases permutation of mutation map values
After processing of all alignments or VAC records mutation algorithm is ended and remaining alignments are written. Checksums of original BAM file and VAC file are added to mutated BAM's header for validation in decrypt process and better file management. Likewise checksum of mutated BAM is added to BDIFF's header. Because BDIFF format also defines that genomic range must be present in header, full range is used.

Last step of encryption command is actual encryption of BDIFF file. Public key can not be used for encryption directly because of RSA encryption has limited input size. Instead AES encryption with randomly generated key is used to encrypt the whole file. Subsequently AES key is encrypted with the private key and stored as part of encrypted BDIFF file. 

### decrypt
Decrypt is reverse command for encrypt. Whole mutated BAM or it's subrange can be restored to it's original form with BDIFF file and correct private key provided. By default BAM is restored in range defined by BDIFF, which is full range when BDIFF is first created after encryption. User can also specify valid genomic range used to restore original BAM. Algorithm for restoring BAM is the same as for mutating except mutation map is not created by random but bases are taken from BDIFF record as mutation map keys and A, C, T, G are assigned to them as their values. Thus mutation map for decryption is inverted mutation map for encryption.

### reencrypt
Owner who encrypted BAM with his public key can use user's public key to grant him access to original BAM or it's subregion. Owner must provide his private key to decrypt BDIFF file and desired public key for second encryption. Intended use case is to provide genomic range so overlapping portion of BDIFF is sliced and encrypted. BAM file must be also provided for validation and genomic position conversion. When new reencrypted BDIFF is created, user can use it in decrypt command with his private key.

## discussion
Secure pseudo random number generator is used to generate AES keys and for multinomial distribution to increase security. Rate of depersonalization is given by VAC file generated from VCF file thus the user is responsible for providing VCF covering all individual related SNV's. User can provide any VCF, but for purpose of depersonalization VCF with large population sample is assumed so allele counts represent general state in population.

## conclusion
We created a tool suitable for storing BAM files in depersonalized form with option to restore their original content. Authority managing BAM files can distribute access to specific genomic regions of original BAM file among users by using reencrypt command. In the future we intend to focus on depersonalizing another BAM file format fields such as new format of CIGAR string, base quality and optional fields.
