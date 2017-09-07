# Varlock

## Introduction
The human genome is a sequence of more than 3 billion DNA bases: adenine, thymine, guanine and cytosine. In technical terms, it is a word over an alphabet {A, T, G, C}. Any two human individuals share over 99.9% of all genome, but at the same time each genome is unique to each person.

Modern technologies provide efficient tools for determining the precise order of nucleotides within a DNA molecule - _DNA sequencing_. Usually, genome analysis starts by sequencing fragments from a DNA sample. These fragments are called _reads_ and they need to be aligned (mapped) to a reference genome for biological interpretation. Aligned reads are typically stored in a (Binary) Sequence Alignment Map file format (_SAM_ or _BAM_). These alignments reveal differences between the studied and the reference genome. A single difference is called a _genomic variant_ or simply a _variant_. A separate file format called _Variant Call Format_ (_VCF_) was derived from SAM for storing genomic variants.

The most common variant is a _Single Nucleotide Variant_ (_SNV_), a variation in a single DNA base that occurs at a specific position in the genome. INDELs - insertions and deletions of DNA bases are also common variants. Different variants found at a single position of a gene represent different forms of this gene called _alleles_.

Presence of specific alleles in a person's genome can tell about his health, predisposition to specific diseases, population of origin, physical traits, etc. These properties make the genomic data very valuable for research and diagnosis in a future health care. Genomic data extraction and interpretation is more affordable and complex every year as research in these fields progresses. While access to meaningful data is necessary for better understanding of human biology, it also provides space for potential abuse, since alleles can be very effectively used for person identification. 
 
Most variant studies require access only to short regions of specific genes and some analyses do not even require information about genomic variants at all. For example, chromosome aneuploidy detection needs only information about number of reads per chromosome or chromosome parts - information that cannot be abused for person identification. 

Due to this potential threat to privacy, BAM files need to be stored in a secure form while preserving general information such as position of reads. These depersonalized BAM files can be then made available to general public without fear of personal data leakage. Moreover, the authority managing genomic data must be able to provide access to unaltered regions of BAM on a valid request from an authorized user. These requirements are not easy to fulfill effectively, since the common size of BAM is in order of gigabytes. To overcome this challenge we developed Varlock - a command line interface utility for securing BAM files in a manner that enables distributing access to specific BAM regions among users.

## Methods
Varlock is a command line interface utility for depersonalization of BAM files by altering SNVs and INDELs based on supplied VCF file. This method works with mapped alignments thus unmapped alignments are handled separately using stream cipher encryption.

The differences between the original and the depersonalized BAM are stored in encrypted form, so only authorized users can recover the original BAM. Varlock uses RSA asymmetric encryption, with pairs of public and private keys for encryption and decryption, respectively.

Two new binary file types are introduced for internal purposes: the first is _Variant Allele Count_ (_VAC_), which is essentially a compact VCF file that represents standard allele frequencies of a population. The second is _BAM difference_ (_BDIFF_), where the differences between the original and depersonalized BAM are compactly stored. Varlock provides commands for encryption, decryption, granting access, and conversion VCF file to VAC file. 

### Variant Allele Count (VAC) file format
VAC is a compact binary file format for storing SNV and INDEL allele frequencies among the population. VAC file contains header with number of stored SNVs and INDELs. Each record represents one variation in a VCF file and contains its absolute position in the genome and observed allele frequencies in the population. SNV record contains four numbers as frequencies of A,T,G,C bases respectively. INDEL record features a map of frequencies and coresponding DNA sequences. SNV and INDEL records are stored separately in two blocks and are sorted by genomic position.

### BAM Difference (BDIFF) file format 
SNV and INDEL differences between original and depersonalized BAM file are stored in a binary BDIFF file format. The file header contains number of SNVs and INDELs stored in the file and can contain arbitrary data. File index stored after the header of the file provides fast seeking of genomic positions. A single BDIFF record represents a difference at a single position in the genome. Since BAM file usually contains multiple alignments covering one genomic position (where they can differ) and individual can have different alleles of the same gene, the BAM difference is stored as a map representing bijection from set of original alleles to set of depersonalized ones and vice versa. In case of SNV the map is a bijection from set of four DNA letters to permuted set of itself. In case of INDEL the map is a bijection from set of known INDEL sequences at this genomic position to permuted set of itself.

BDIFF file act as the secret key in context of varlock encryption, thus it is actually never stored on disk in plain form. All the writing and reading is done in memory using binary stream.


### Encrypt
Process involving BAM depersonalization is simply called _encryption_. User must provide the original BAM file, VAC file, and his public key for encryption. Both VAC and BAM files are iterated at the same time to write a new depersonalized BAM together with BDIFF tracking all the changes. Every time a genomic position of VAC record intersects with one or more alignments from BAM file a pseudo random permutation of DNA base letters (in case of SNV) or INDEL sequences is generated. Those alignments that do not interfere with VAC records are written to depersonalized BAM as they are. Encryption always involves the whole content of supplied BAM, meaning that all alignments from original BAM are present in depersonalized BAM either in altered or unaltered.

The permutation is generated in an iterative manner. At each step of a SNV case a base is drawn by random from a multinomial distribution of allele counts defined by a VAC record. The base is then mapped to the most abundant base from the alignment column and only remaining bases are considered in next steps. If base is not found in alignment column it has zero frequency. Unknown base 'N' is always mapped to itself. This effectively makes this permutation a representative of the SNV distribution in population.

In case that VAC record represents INDEL, permutation set must be found first. Each alignment intersecting with a VAC record is searched for longest match by one of INDEL sequences supplied by VAC record. If match is found matching sequence is added to permutation set, otherwise alignment is omitted from encryption. Permutation generation is same as with SNV case, only instead of bases sequences are permuted now.

The permutation is used to alter intersecting alignments in the place of variant. Original bases or sequences are mapped to permuted ones. For each alignment in SNV case one base is replaced, while in INDEL case it is the longest matching sequence described before. When the generated permutation is not an identity, a new record is written to BDIFF file with the position of SNV and generated permutation. After all alignments or VAC records are processed, remaining alignments are written to the BAM file. Checksum of the original BAM file and checksum of the VAC file is added to the header of the depersonalized BAM for later validation and file management purposes. At the same time, a checksum of depersonalized BAM is added to the header of BDIFF file.

The process described above handles mapped alignments only because it is dependent on alignment position. Unmapped alignments are handled separately. Their sequences are encrypted completely using stream cipher encryption which has the ability to produce cipher with the same length as input message. At first secret key is randomly generated for all unmapped alignments. Later this key is stored within header of BDIFF file. When unmapped alignment is found, the secret key together with alignment's query template name is hashed by SHA-512 algorithm producing 512 bits long hash. Alignment's sequence is then encrypted using stream cipher encryption where obtained hash act as key stream. Each 2 bits of the key stream are used to encrypt one DNA base from the input sequence, thus key size is enough to uniquely encrypt sequence of 256 bases. If the sequence is longer, the key stream is repeated. Each DNA base is encrypted combining 2 bits code of the BASE and next 2 bits of the key stream using XOR operation, where the result is bit code of encrypted base. The original secret key is stored within header of produced BDIFF file.

Besides genomic range of alignments contained in BDIFF file there must be also information about effective range of BDIFF which is an exact range that depersonalization covers. It is intended for later use in reencryption operation. Range of depersonalization can be bigger or slightly smaller than the range of alignments in BAM file, denoted by start of first alignment and end of last alignment. Start and end positions of an effective range are stored in BDIFF header. Maximum / full effective range is written when depersonalizing BAM - encryption operation always covers the whole BAM content.

After the finished depersonalized BAM file is closed, it is checksumed for validation purposes and the checksum is stored in BDIFF header.

Last step of encryption is the actual encryption of BDIFF file. Public key can not be used for the encryption directly because the RSA encryption has a limited input size. Instead, the AES encryption with a randomly generated key is used to encrypt the whole file. The AES key is then encrypted with the private key and stored before the actual encrypted BDIFF file. In this manner an access to the content of original BAM is restricted to an owner of the private key paired with the public key used for encryption. Since genomic information may be fully reconstructed from depersonalized BAM file and associated (encrypted) BDIFF file, the original, unsecured BAM file may be deleted after the whole genome encryption.

To increase security and integrity of the workflow the plain BDIFF file is signed with supplied private key. The signature is subsequently 
 stored before the actual encrypted BDIFF file (next to encrypted AES key).
 

### Decrypt
Decrypt is the inverse operation to encrypt where private parts of BAM file get restored to original form. The user must provide depersonalized BAM and encrypted BDIFF file along with the RSA private key whose public conterpart was used in BDIFF encryption. Whole depersonalized BAM or its part can be restored to their original form when a BDIFF file and a corresponding private key is provided. The BAM file can be restored in the effective range defined by BDIFF or in a subrange of it provided by the user who can also decide how to handle unmapped alignments. They can be included in restored BAM or not. It is also possible to restore unmapped alignments only and ommit all mapped ones.

The algorithm first reads encrypted AES key and file signature both stored before the actual encrypted BDIFF file. AES key gets decrypted with owner's private key and BDIFF file is subsequently decrypted with it. The decrypted file is verified with owner's public key against the signature to prove its origin.

The original BAM is then restored in a same way it was depersonalized, only the permutations are not generated but obtained from the corresponding BDIFF records as an inverse permutations. Permuted base or sequences are mapped to original ones.

When the effective range of decryption is supplied, only alignments that are within or intersecting with this range are written to output BAM file. Moreover depersonalization covers exactly this range, so boundary alignments may be only partially restored. 

Implemented stream cipher encryption method is symetric thus unmapped alignments are decrypted using the very same method. Secret key for decryption is obtained from BDIFF header, where it was previously stored.


### Grant access (reencrypt)
Owner of the private key which was used to encrypt BDIFF file can grant access to restore the original data from depersonalized BAM file. BDIFF file is first decrypted by owner's private key and then encrypted by public key of another user (assuming he has the private counterpart), so only he can decrypt it. Owner can choose to grant only partial access to original data by providing subrange of BDIFF's effective range. In this case all BDIFF records outside of this range are omitted so this range becomes the new effective range and is stored in BDIFF header. During the proces decrypted BDIFF gets verified with owner's public key and new BDIFF is newly signed with his private key.

Granting access also involves depersonalized BAM because it contains the needed information to resolve the effective range and it is checksumed against the value stored in BDIFF header for validation, although no reading of alignments occur at this operation.

 The process may be repeated with different combinations of genomic ranges and public keys, producing separate access rights for individual users.

### Create VAC
User must convert his VCF file to VAC file to use it in encryption. Since VCF records are dynamicaly typed, variant type must be infered from alleles in the record. Only SNV's and INDEL's are converted to VAC file. We define SNV as record which contains at least one reference and one allele while each of them must be one of DNA bases. INDEL must also have at least one reference and one alternative allele while each allele can consist only of DNA bases. Moreover INDEL must have at least one allele with more than one base - otherwise it would be SNV. Besides VCF user must also supply BAM file as it contains information neccessary to resolve genomic indices. If generated VAC file meets user's assumption that it contains a all exploitable variants it can be widely reused for encryption of different samples.

## Discussion
Biggest security question lies in a comprehensiveness of supplied VCF file as it is assumed that it contains all SNV's and INDEL's that could be exploited in indentification of the person. Since that the robustness of this method heavily depends on the content of provided VCF file, i.e., only the variants from the VCF are depersonalized. The user is thus responsible for providing a VCF file covering all individual related SNVs. For purpose of depersonalization, a VCF with large population sample is needed, where the allele counts represent general distribution in population.
 
Permutation of bases and sequences in encryption operation should be representative of the variant distribution among the general population. When generating this permutation multinomial distribution is used to randomly draw alternative variant, but the reference variant is chosen in deterministic way (the most common one). It is matter of question wheter the reference variant should be chosen randomly using multinomial distribution of reference variant frequencies to be a better representative of general random variant..
 
 A secure pseudo-random number generator or a true random number generator must be used in multinomial drawing and to generate AES keys to ensure security.
 
 Size of an encrypted BAM file should not change significantly. Size of a BDIFF file is dependent on the number of covered variations from VAC file in an encrypted genomic range.

## Conclusion
Varlock is a tool suitable for storing BAM files in a depersonalized form with the option to restore their original content. The authority managing BAM files can distribute access to specific genomic regions of original BAM file among users, while other parts of the original BAM file stay inaccessible. Now, only the sequence is depersonalized, however, we intend to extend this concept in the future to other BAM file format fields such as CIGAR string, base quality, and optional fields of the alignment section.

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


## Samtools help
samtools view -H examples/resources/sample.mut.bam

samtools view -h resources/out.mut.bam >> out.mut.sam

less chr22.vcf.gz | egrep -v "^#" | cut -f 4,5 | awk 'length($1) > 2 || length($2) > 2'





