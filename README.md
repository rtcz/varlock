# Varlock

## Introduction
The human genome is a sequence of more than 3 billion DNA bases: adenine, thymine, guanine and cytosine. In technical terms, it is a word over an alphabet {A, T, G, C}. Any two human individuals share over 99.9% of all genome, but at the same time each genome is unique to each person.

Modern technologies provide efficient tools for determining the precise order of nucleotides within a DNA molecule - _DNA sequencing_. Usually, genome analysis starts by sequencing fragments from a DNA sample. These fragments are called _reads_ and they need to be aligned (mapped) to a reference genome for biological interpretation. Aligned reads or _alignments_ are typically stored in a (Binary) Sequence Alignment Map file format (_BAM_ or _SAM_). These alignments reveal differences between the studied and the reference genome. A single difference is called a _genomic variant_ or simply a _variant_. A separate file format called _Variant Call Format_ (_VCF_) was derived from SAM for storing genomic variants.

The most common variant is a _Single Nucleotide Variant_ (_SNV_), a variation in a single DNA base that occurs at a specific position in the genome. INDELs - insertions and deletions of DNA bases are also common variants. Different variants found at a single position of a gene represent different forms of this gene called _alleles_.

Presence of specific alleles in a person's genome can tell about his health, predisposition to specific diseases, population of origin, physical traits, etc. These properties make the genomic data very valuable for research and diagnosis in a future health care. Genomic data extraction and interpretation is more affordable and complex every year as research in these fields progresses. While access to meaningful data is necessary for better understanding of human biology, it also provides space for potential abuse, since alleles can be very effectively used for person identification. 
 
Most variant studies require access only to short regions of specific genes and some analyses do not even require information about genomic variants at all. For example, chromosome aneuploidy detection needs only information about number of reads per chromosome or chromosome parts - information that cannot be abused for person identification. 

Due to this potential threat to privacy, BAM files need to be stored in a secure form while preserving general information such as position of reads. These depersonalized BAM files can be then made available to general public without fear of personal data leakage. Moreover, the authority managing genomic data must be able to provide access to unaltered regions of BAM on a valid request from an authorized user. These requirements are not easy to fulfill effectively, since the common size of BAM is in order of gigabytes. To overcome this challenge we developed Varlock - a command line interface utility for securing BAM files in a manner that enables distributing access to specific BAM regions among users.

## Methods
Varlock is a command line interface utility for depersonalization of BAM files by altering SNVs and INDELs based on supplied VCF file. This method works with mapped alignments only, thus unmapped alignments are handled separately using stream cipher encryption.

The differences between the original and the depersonalized BAM are stored in encrypted form, so only authorized users can recover the original BAM. Varlock uses RSA asymmetric encryption, with pairs of public and private keys for encryption and decryption, respectively.

Two new binary file types are introduced for internal purposes: the first is _Variant Allele Count_ (_VAC_), which is essentially a compact VCF file that represents standard allele frequencies of a population. The second is _BAM difference_ (_BDIFF_), where the differences between the original and a depersonalized BAM are compactly stored. Varlock provides commands for encryption, decryption, granting access, and conversion VCF file to VAC file. 

### Variant Allele Count (VAC) file format
VAC is a compact binary file format for storing SNV and INDEL allele frequencies among the population. VAC file contains a header with a number of stored SNVs and INDELs. Each record represents one variation in a VCF file and contains its absolute position in the genome and observed allele frequencies in the population. Each SNV record contains four numbers as frequencies of A,T,G,C bases respectively. INDEL record features a map of frequencies and corresponding DNA sequences. SNV and INDEL records are stored separately in two blocks and are sorted by their genomic position.

### BAM Difference (BDIFF) file format 
SNV and INDEL differences between original and depersonalized BAM file are stored in a binary BDIFF file format. The BDIFF file format consists of a header, a file index, and a number of BDIFF records. The file header contains number of SNVs and INDELs stored in the file and arbitrary metadata. Secondly, the file index provides fast seeking of genomic positions. A single BDIFF record represents a difference at a single position in the genome. Since BAM file usually contains multiple alignments covering one genomic position (where they can differ) and individual can have different alleles of the same gene, the BAM difference is stored as a map representing bijection from set of original alleles to set of depersonalized ones and vice versa. In case of an SNV, the map is a bijection from the set of four DNA letters to its permutation. In case of an INDEL, the map is a bijection from a set of known INDEL sequences at this genomic position to its permutation.

BDIFF file acts as the secret key in the context of Varlock encryption, thus it is actually never stored on a medium in plain form. BDIFF file is encrypted symmetrically with the AES encryption using a randomly generated key. Only the AES key is then encrypted asymmetrically by the public key of a user and stored in the BDIFF header. 

### Encrypt
Process involving BAM depersonalization is simply called _encryption_. User must provide the original BAM file, VAC file, and his public key for encryption. Both VAC and BAM files are iterated at the same time to write a new depersonalized BAM together with BDIFF tracking all the changes. Every time a genomic position of VAC record intersects with one or more alignments from the BAM file, a pseudo random permutation of DNA base letters (in case of SNV) or INDEL sequences is generated. Alignments that do not overlap with any VAC records are written to the depersonalized BAM unchanged. Encryption always involves the whole content of the supplied BAM, meaning that all alignments from it are present in the depersonalized BAM either altered or unaltered.

#### Encryption of mapped alignments

The permutation required to encrypt each SNV record is generated in an iterative manner. At each step, a base is drawn randomly from a multinomial distribution of allele counts defined by the VAC record. The base is then mapped to the most abundant base from the alignment column of the VAC record and only remaining bases are considered in next steps. The base has zero frequency if it is not found in the alignment column. Unknown base 'N' is always mapped to itself. This effectively makes this permutation a representative of the SNV distribution in a population as defined by the VAC file.

In a case that a VAC record represents an INDEL, the encryption process is similar, but on a set defined by the VAC record. The permutation generation is the same as in the SNV case, but the resulting permutation has a variable length now. Since the supplied BAM file usually does not contain all of the possible INDELs on the position, only the relevant parts of the permutation are stored for storage efficiency.

The generated permutation is used to alter intersecting alignments in the position of a variant of the original BAM to create the depersonalized BAM. When the generated permutation is not an identity, a new record is written to BDIFF file with the position of SNV and the generated permutation. After all VAC records are processed, remaining alignments are written to the depersonalized BAM file unchanged. 

#### Encryption of unmapped reads

Unmapped reads are encrypted completely using stream cipher encryption which has the ability to produce cipher with the same length as the input. At first, a single secret key is randomly generated for all unmapped reads. This key is stored in the BDIFF file header. When an unmapped read is found, its template name and the secret key are hashed by SHA-512 algorithm producing 512 bits long hash. The hash is then used to encrypt the sequence of the read. Each 2 bits of the hash are used to encrypt one DNA base from the input sequence using a simple XOR operation. Thus, the key size is enough to uniquely encrypt sequence of 256 bases. If the sequence is longer, the key is repeated. Unknown bases ('N') are skipped in encryption.   

#### Finalization

Checksum of the original BAM file and checksum of the VAC file ia added to the header of the depersonalized BAM for validation and file management purposes. Afterwards, a checksum of the depersonalized BAM is added to the header of BDIFF file.

Due to the fact that positions without a VAC record are not encrypted or stored in the BDIFF file, BDIFF file header also contains the exact range that the depersonalization covers - _effective range_. It is clear this way, exactly which part of genome is encrypted. In the standard encryption operation, the effective range is maximal, i.e. the whole genome. An user that has access to a BDIFF file with a certain effective range can produce a 'subrange' BDIFF file for other user with a smaller (or equal) range. This process is called _reencryption_ and it is explained later. 

The actual encryption of BDIFF file is explained in the BDIFF file format specification. In this manner, an access to the content of the original BAM is restricted to the owner of the private key paired with the public key used for encryption. Since genomic information may be fully reconstructed from depersonalized BAM file and the associated (encrypted) BDIFF file, the original BAM file may (and should) be deleted after the whole genome encryption.

Furthermore, the plain BDIFF file is signed with supplied private key to increase security and integrity of the workflow. The signature is stored in the start of the encrypted BDIFF file.
 
### Decrypt

Decrypt is the inverse operation to encrypt where private parts of a depersonalized BAM file are restored to the original form. The user must provide a depersonalized BAM and encrypted BDIFF file along with the RSA private key whose public counterpart was used in the BDIFF encryption. Whole depersonalized BAM or its part can be restored to their original form when a BDIFF file and a corresponding private key is provided. Unmapped reads are handled separately and the user can decide whether to include them into decryption or not (provided he has granted access to them).

The algorithm first reads the encrypted AES key and the file signature both stored before the actual encrypted BDIFF file. The AES key is decrypted with the owner's private key and the BDIFF file is subsequently decrypted with it. The decrypted file is verified with owner's public key against the signature to prove its origin.

The original BAM is then restored in a same way it was depersonalized, only the permutations are not generated but obtained from the corresponding BDIFF records as inverse permutations. Permuted bases or sequences are mapped to original ones.

When the effective range of decryption is supplied, only alignments that are within or intersecting with this range are written to the output BAM file. Moreover, the depersonalization process covers exactly this range, so boundary alignments may be only partially restored. 

Implemented stream cipher encryption method is symmetric thus unmapped alignments are decrypted using the very same method. Secret key for decryption is obtained from BDIFF header, where it was previously stored.

### Grant access (reencrypt)

Owner of the private key which was used to encrypt a BDIFF file can grant access to restore the original data from a depersonalized BAM file. BDIFF file is first decrypted by the owner's private key and then encrypted by a public key of another user, so only he can decrypt it using his private key. Owner can choose to grant only partial access to original data by providing a subrange of BDIFF's effective range. In this case all BDIFF records outside of the new range are omitted and this range becomes the effective range of the new BDIFF file. Decrypted BDIFF file gets verified with owner's public key and the new BDIFF is signed with his private key.

The depersonalized BAM file is also needed for reencryption, since it contains the chromosome lengths needed to resolve the effective range. Besides that, the BAM file's checksum is compared to the checksum stored in the BDIFF file header. This ensures that the BDIFF file belongs to the BAM file and that the BAM file was not modified.

 The process can be repeated with different combinations of genomic ranges and public keys, producing separate access rights for individual users.

### VAC file creation

User must convert his VCF file to VAC file to use it in encryption. Since VCF records are dynamically typed, the variant type (SNV/INDEL/other) must be inferred from alleles in the record. Only SNVs and INDELs are converted to VAC file. Each record, where each of reference and alternative alleles is only one DNA base, is categorized as SNV record. Records that either reference or alternative alleles contain longer sequences are categorized as INDELs. VAC creation process furthermore needs information about the genomic index (chromosome lengths of reference genome), included in each BAM file. Thus a target BAM file is needed for VAC creation. The generated VAC file then can be used to encrypt multiple BAM files, provided that they used the same reference genome. 

## Discussion

The quality of depersonalisation is directly linked with the quality and quantity of variants in the provided VCF file, since variants not present in the VCF file are not encrypted. Therefore, the biggest security question lies in a comprehensiveness of supplied VCF file, since it is assumed that the VCF file contains all variants that could be exploited in the person identification. For the purpose of depersonalization, a VCF with large population sample is needed, where the allele counts represent general distribution in population. In the Varlock concept, the user is responsible to provide a good quality VCF file for his population of interest.
 
Permutations of bases and sequences in encryption operation are generated by a random draw from the alternative variants according to multinomial distribution of their abundances. Therefore, they should represent the variant distribution among the general population. However, the reference variant is chosen in a deterministic way (the most common one from the remaining set). It is a matter of question, whether the reference variant should be chosen randomly in the same way as alternative variant. Probability properties of permutation generation with both (and possibly more) approaches should be evaluated to agree upon a final decision.  

 A secure pseudo-random number generator or a true random number generator must be used to generate cryptographic keys and in the multinomial drawing to ensure security.
 
 The size of the encrypted BAM file should not change significantly, since the symmetric cipher is used. However, the size of a BDIFF file is dependent on the number of covered variations from the VCF file, thus an attacker can guess the number of changes made, but this is not considered a security threat.
 
### Random mutations

The encryption process currently works only on variants in the VCF file. The VCF file and the depersonalized BAM file is considered public, therefore everybody can reconstruct the parts of the original BAM file that do not have a corresponding variant in the VCF file. This introduces a vulnerability when there exist variants outside the VCF, that can be used for person identification. This vulnerability can be partially negated with the introduction of _randomly mutated bases or INDELs_. 

The concept is pretty simple - after the depersonalization of the BAM file, a set of genomic positions is picked and a random "de novo" variant is introduced on each position. These variants are written in the BDIFF file as the normal variants. The decryption then works as usual.

The number of the new variants should be high enough to disallow attacks in-between the variants from the VCF file. On the other hand, it should be the lowest possible, since the size of the BDIFF file and time cost of all operations linearly increases with this number.       

### Cigar strings and sequencing quality

The BAM file does not contain only sequence alignments, but other data, that should be taken care of in the depersonalization process. For example, cigar string is a compressed representation of an alignment, that can change when the alignment sequence changes. Therefore, cigar strings should be adjusted to corresponding new alignments, otherwise it would be easy (in some cases) to guess where a change was made. 

The mapping quality is other type of data that is available for every alignment in the BAM file. While it is not needed to do anything with the qualities in case of a SNV, in case of an INDEL the length of the sequence changes, so it is necessary to adjust the length of the mapping quality in this case. When the length of the depersonalized BAM is greater, thus we introduce an insertion, we can provide any qualities as long as these qualities will not betray the position of the change. However, when we introduce a deletion, the number of qualities fall and we need to store the remaining ones into the BDIFF file for each of the alignments. This brings another component into the BDIFF file, since until now, it contained only data about particular positions or variants and not about single alignments. This can potentially make problems in case of high coverage BAM files.           

It is possible, that also other (meta-)data in the BAM files are able to reveal the introduced change and the algorithm should be adjusted to include them.    

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





