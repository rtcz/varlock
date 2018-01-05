# Personal genome anonymization

## Introduction
Human genome is a linear sequence of more than 3 billion DNA nucleotides, called bases. In technical terms, it is a word over an alphabet {A, T, G, C}, representing 4 nucleotide types: adenine, thymine, guanine, and cytosine. Any two human individuals share over 99.9% of all genome, but the remaining difference is sufficient to ensure uniqueness of each person's genome.

Modern technologies provide efficient tools for determining the precise order of nucleotides within a DNA molecule, with a laboratory process called _DNA sequencing_. Usually, genome analysis starts by sequencing fragments from a DNA sample. Sequenced fragments, named _reads_, represents only a small portion of genomic information, generally around few hundred subsequent bases. They need to be compared with a representative human genome sequence, _reference genome_ for a biological interpretation. 

Comparison analysis starts with a process called _alignment_ or _mapping_ that assigns the most similar region of the reference to a read. Reads without a sufficient similarity match are considered _unmapped_. Aligned reads (_alignments_) and unmapped reads are stored in a (Binary) Sequence Alignment Map file format (_BAM_ or _SAM_). These alignments reveal differences between the studied and the reference genome. 

A single difference between genomes is called a _genomic variant_ or simply a _variant_. The most common variant is a _Single Nucleotide Variant_ (_SNV_), a variation in a single DNA base that occurs at a specific position in the genome. Other common type of variants is INDEL - insertion or deletion of one or more bases. Different variants found at a single position of a gene represent different forms of this gene called _alleles_. Genomic variants are retrieved from aligned reads using _variant calling_ tools and they are then stored in _Variant Call Format_ (_VCF_) files. 

Presence of specific alleles in a person's genome can tell about his health, predisposition to specific diseases, population of origin, physical traits, etc. These properties make the genomic data very valuable for research and diagnosis in health care. Genomic data extraction and interpretation is more affordable and complex every year as research in these fields progresses. While access to meaningful data is necessary for better understanding of human biology, it also provides a space for potential abuse, since alleles can be very effectively used for identification of a person and his phenotypic traits. 
 
Most variant studies require limited access to short regions of specific genes only, such as diagnosis of specific disorder with a known set of causal genes. In addition, some analyses do not even require information about genomic variants at all. For example, chromosome aneuploidy detection needs only information about number of reads aligned to individual chromosomes that cannot be abused for person identification. 

Due to this potential threat to privacy, BAM files need to be stored in a secure form while preserving general information about alignments. These depersonalized BAM files can be then made available to the general public without fear of personal data leakage. Moreover, the authority managing genomic data must be able to provide access to unaltered regions of BAM on a valid request from an authorized user. These requirements are not easy to fulfill effectively, since the common size of BAM is in order of gigabytes. 

To overcome this challenge, we developed Varlock - a command line interface utility for securing BAM files that enables distributing access privileges to specific BAM regions among users. Varlock also includes separate command line interface utility called Anonymizer for irreversible depersonalization of FASTQ files with paired-end reads.


## Varlock
Varlock is a command line interface utility for reversible depersonalization of BAM file by altering contained SNVs and INDELs based on a supplied VCF file. This method works with mapped alignments only and unmapped alignments are handled separately using stream cipher encryption.

The differences between the original and the depersonalized BAM are stored in encrypted form, so only authorized users can recover their permitted parts of the original BAM. Varlock uses RSA asymmetric encryption with pairs of public and private keys for encryption and decryption, respectively.

Two new binary file types are introduced for internal purposes: the first is a population specific _Variant Allele Count_ (_VAC_) file, which is essentially a compact VCF file that represents standard allele frequencies of a population. The second is a BAM file specific _BAM difference_ (_BDIFF_) file, where the differences between the original and a depersonalized BAM are compactly stored. 

Varlock provides commands for encryption, decryption, granting access, and conversion from VCF files to VAC files. 

### Variant Allele Count (VAC) file format
VAC is a compact binary file format for storing SNV and INDEL allele frequencies among the population. VAC file contains a header with a number of stored SNVs and INDELs. Each record represents one variation in a VCF file and contains its absolute position in the genome and observed allele frequencies in the population. Each SNV record contains four numbers as frequencies of A,T,G,C bases respectively. INDEL record features a map of frequencies and corresponding DNA sequences. SNV and INDEL records are stored separately in two blocks and are sorted by their genomic position.

### BAM Difference (BDIFF) file format 
SNV and INDEL differences between original and depersonalized BAM file are stored in a binary BDIFF file format. The BDIFF file format consists of a header, a file index, and a number of BDIFF records. The file header contains number of SNVs and INDELs stored in the file and metadata required for decryption process. Secondly, the file index enables fast seeking of genomic positions. 

A single BDIFF record represents a difference at a single position in the genome. Since BAM file usually contains multiple alignments covering one genomic position (where they can differ) and individual can have different alleles of the same gene, the BAM difference is stored as a map representing bijection from set of original alleles to set of depersonalized ones and vice versa. In case of an SNV, the map is a bijection from the set of four DNA letters to its permutation. In case of an INDEL, the map is a bijection from a set of known INDEL sequences at this genomic position to its permutation.

BDIFF file acts as the secret key in the context of Varlock encryption, thus it is actually never stored on a medium in plain form. BDIFF file is encrypted symmetrically with the AES encryption using a randomly generated key. Only the AES key is then encrypted asymmetrically by the user's public key and stored in the BDIFF header. 

### Encryption
Depersonalization of a BAM file is called _encryption_. User must provide the original BAM file, a VAC file of the closest population, and his public key for encryption. Both VAC and BAM files are iterated at the same time to write a new depersonalized BAM together with BDIFF tracking all the changes. Every time a genomic position of a VAC record intersects with one or more alignments from the BAM file, a pseudo random permutation of DNA base letters (in the case of SNV) or INDEL sequences is generated. Alignments that do not overlap with any VAC records are written to the depersonalized BAM unchanged. Encryption always involves the whole content of the supplied BAM, meaning that all alignments from it are present in the depersonalized BAM either altered or unaltered.

#### Encryption of mapped alignments

The permutation required to encrypt each SNV record is generated in an iterative manner. At each step, a base is drawn randomly from a multinomial distribution of allele counts defined by the VAC record. The base is then mapped to the most abundant base from the alignment column of the VAC record and only remaining bases are considered in next steps. The base has zero frequency if it is not found in the alignment column. Unknown base 'N' is always mapped to itself. This effectively makes this permutation a representative of the SNV distribution in a population as defined by the VAC file.

In a case that a VAC record represents an INDEL, the encryption process is similar, but on a set defined by the VAC record. The permutation generation is the same as in the SNV case, but the resulting permutation has a variable length now. Since the supplied BAM file usually does not contain all of the possible INDELs on the position, only the relevant parts of the permutation are stored for storage efficiency.

The generated permutation is used to alter intersecting alignments in the position of a variant of the original BAM to create the depersonalized BAM. When the generated permutation is not an identity, a new record is written to BDIFF file with the position of SNV and the generated permutation. After all VAC records are processed, remaining alignments are written to the depersonalized BAM file unchanged. 

#### Encryption of unmapped reads

Unmapped reads are encrypted completely using stream cipher encryption which has the ability to produce cipher with the same length as the input. At first, a single secret key is randomly generated for all unmapped reads. This key is stored in the BDIFF file header. When an unmapped read is found, its template name and the secret key are hashed by SHA-512 algorithm producing 512 bits long hash. The hash is then used to encrypt the sequence of the read. Each 2 bits of the hash are used to encrypt one DNA base from the input sequence using a simple XOR operation. Thus, the key size is enough to uniquely encrypt sequence of 256 bases. If the sequence is longer, the key is repeated. Unknown bases ('N') are skipped in encryption.   

#### Finalization of encryption

Checksum of the original BAM file and checksum of the VAC file is added to the header of the depersonalized BAM for validation and file management purposes. Afterwards, a checksum of the depersonalized BAM is added to the header of BDIFF file.

Due to the fact that positions without a VAC record are not encrypted or stored in the BDIFF file, BDIFF file header also contains the exact range that the depersonalization covers - _effective range_. It is clear this way, exactly which part of genome is encrypted. In the standard encryption operation, the effective range is maximal, i.e. the whole genome. An user that has access to a BDIFF file with a certain effective range can produce a 'subrange' BDIFF file for other user with a smaller (or equal) range. This process is called _reencryption_ and it is explained later. 

The actual encryption of BDIFF file is explained in the BDIFF file format specification. In this manner, an access to the content of the original BAM is restricted to the owner of the private key paired with the public key used for encryption. Since genomic information may be fully reconstructed from depersonalized BAM file and the associated (encrypted) BDIFF file, the original BAM file may (and should) be deleted after the whole genome encryption.

Furthermore, the plain BDIFF file is signed with a supplied private key to validate its authenticity and integrity. The signature is stored at the beginning of the encrypted BDIFF file and it is verified in the decryption process using the associated public key.
 
### Decryption

Decryption is the inverse operation to encryption where secured parts of a depersonalized BAM file are restored to the original form. A user must provide depersonalized BAM and an associated, encrypted BDIFF file along with the RSA private key whose public counterpart was used in the BDIFF encryption. Whole depersonalized BAM or its part can be restored to its original form according to user privileges specified in the BDIFF file. Unmapped reads are handled separately and the user can decide whether to include them into decryption or not (provided he has granted access to them).

The algorithm first reads the encrypted AES key and the file signature both stored before the actual encrypted BDIFF file. The AES key is decrypted with the owner's private key and the BDIFF file is subsequently decrypted with it. The decrypted file is verified with owner's public key against the signature to prove its origin.

The original BAM is then restored in a same way it was depersonalized, only the permutations are not generated but obtained from the corresponding BDIFF records as inverse permutations. Permuted bases or sequences are mapped to original ones.

When the effective range of decryption is supplied, only alignments that are within or intersecting with this range are written to the output BAM file. Moreover, the depersonalization process covers exactly this range, so boundary alignments may be only partially restored. 

Unmapped alignments are decrypted with symmetric stream cipher decryption using secret key obtained from the BDIFF header.

### Grant access (reencrypt)

Owner of the private key which was used to encrypt a BDIFF file can grant access to restore the original data from a depersonalized BAM file. BDIFF file is first decrypted by the owner's private key and then encrypted by a public key of another user, so only he can decrypt it using his private key. Owner can choose to grant only partial access to original data by providing a subrange of BDIFF's effective range. In this case all BDIFF records outside of the new range are omitted and this range becomes the effective range of the new BDIFF file. Decrypted BDIFF file gets verified with owner's public key and the new BDIFF is signed with his private key.

The depersonalized BAM file is also needed for reencryption, since it contains the chromosome lengths needed to resolve the effective range. Besides that, the BAM file's checksum is compared to the checksum stored in the BDIFF file header. This ensures that the BDIFF file belongs to the BAM file and that the BAM file was not modified.

 The process can be repeated with different combinations of genomic ranges and public keys, producing separate access rights for individual users.

### VAC file creation

VAC file is derived from a VCF file with the variant frequencies of the closest population. Since VCF records are dynamically typed, the variant type (SNV/INDEL/other) must be inferred from alleles in the record. Only SNVs and INDELs are converted to VAC file. Each record, where each of reference and alternative alleles is only one DNA base, is categorized as SNV record. Records that either reference or alternative alleles contain longer sequences are categorized as INDELs. VAC creation process furthermore needs information about the genomic index (chromosome lengths of reference genome), included in each BAM file. Thus a target BAM file is needed for VAC creation. The generated VAC file then can be used to encrypt all BAM files for the population, provided that they used the same reference genome. 


## Anonymizer
Along with the Varlock, we provide a method for irreversible depersonalization of paired-end reads, called Anonymizer. Original sequences stored in FASTQ files are replaced with corresponding sequences from the reference genome. This way, sequence variation from the reference is removed, while mapping positions remain the same. Consistent mapping is important for several sequence analyses, such as non-invasive prenatal testing and copy number variation detection. The method is implemented as a command line tool that accepts a pair of FASTQ files and a reference genome together with its index.

Anonymizer maps paired-end reads twice. The first mapping generates a primary SAM file with alignments of reads to the reference. Reads from the SAM file are then converted to the FASTQ format and stored in the primary pair of FASTQ files. If a read was mapped, its sequence is replaced with the corresponding part of reference sequence, otherwise the original sequence is kept. This way, the personal sequence variation is eliminated.  
 
The second mapping aligns reads from the primary pair of FASTQ files and creates a secondary SAM file. The secondary pair of FASTQ files is then produced by using both primary and secondary SAM files together with the original pair of FASTQ files in the following manner: alignments from primary and secondary SAM files are iterated in parallel and mapped read positions are compared between them. If the mapped positions between primary and secondary reads match, they are considered as consistently mapped. In this case, the sequence pair from secondary SAM file (originating from the reference sequence) is used to write a secondary FASTQ pair. When a sequence pair is not consistently mapped, a sequence pair from the original pair of FASTQ files is used as output. The pair of secondary FASTQ files is the final output and it is considered as depersonalized.
 
Anonymizer uses internally "Bowtie2" mapper to map reads, albeit it can be substituted by any mapping tool that produce alignments according to the SAM file specification.

## Discussion

The quality of depersonalisation is directly linked with the quality and quantity of variants in the provided VCF file, since variants not present in the VCF file are not encrypted. Therefore, the biggest security question lies in a comprehensiveness of supplied VCF file, since it is assumed that the VCF file contains all variants that could be exploited in the person identification. For the purpose of depersonalization, a VCF with large population sample is needed, where the allele counts represent general distribution in population. If the frequencies are not available, the concept of the encryption may be extended with a random mutations of positions that are not specified in the VCF file.  
 
Permutations of bases and sequences in encryption operation are generated by a random draw from the alternative variants according to multinomial distribution of their abundances. Therefore, they should represent the variant distribution among the general population. However, the reference variant is chosen in a deterministic way (the most common one from the remaining set). It is a matter of question, whether the reference variant should be chosen randomly in the same way as alternative variant. Probability properties of permutation generation with both (and possibly more) approaches should be evaluated to agree upon a final decision.  

 A secure pseudo-random number generator or a true random number generator must be used to generate cryptographic keys and in the multinomial drawing to ensure security.
 
 The size of the encrypted BAM file should not change significantly, since the symmetric cipher is used. However, the size of a BDIFF file is dependent on the number of covered variations from the VCF file, thus an attacker can guess the number of changes made. This potential threat may be overcome by including additional special records, that would change the size of the file, but otherwise do not change the depersonalised BAM file.  
 
### Random mutations

The encryption process currently works only on variants in the VCF file. The VCF file and the depersonalized BAM file is considered public, therefore everybody can reconstruct the parts of the original BAM file that do not have a corresponding variant in the VCF file. This introduces a vulnerability when there exist variants outside the VCF, that can be used for person identification. This vulnerability can be partially negated with the introduction of _randomly mutated bases or INDELs_. 

The concept is pretty simple - after the depersonalization of the BAM file, a set of genomic positions is picked and a random "de novo" variant is introduced on each position. These variants are written in the BDIFF file as the normal variants. The decryption then works as usual.

The number of the new variants should be high enough to disallow attacks in-between the variants from the VCF file. On the other hand, it should be the lowest possible, since the size of the BDIFF file and time cost of all operations linearly increases with this number.       

### Cigar strings and sequencing quality

BAM files do not contain only sequence alignments, but other data, that should be taken care of in the depersonalization process. For example, cigar string is a compressed representation of an alignment, that can change when the alignment sequence changes. Therefore, cigar strings should be adjusted to corresponding new alignments, otherwise it would be easy (in some cases) to guess where a change was made. 

Every read in the BAM file contains its _sequencing quality_ or _base quality_. While it is not needed to do anything with the qualities in the case of a SNV, in the case of an INDEL the length of the sequence changes, so it is necessary to adjust the length of the mapping quality in this case. When the length of the depersonalized BAM is greater, thus we introduce an insertion, we can provide any qualities as long as these qualities will correspond to surrounding qualities of the changed position. However, when we introduce a deletion, the number of qualities decrease and we need to store the remaining ones into the BDIFF file for each of the alignments. This brings another component into the BDIFF file, since until now, it contained only data about particular positions or variants and not about single alignments. This can potentially make problems in the case of high coverage BAM files.           

It is possible, that also other (meta-)data in the BAM files are able to reveal the introduced change and the algorithm should be adjusted to include them.    

## Conclusion

We developed methods for reversible and irreversible removal of personal sequence variation from nucleotide sequences. Methods are demonstrated on reads obtained from the NGS sequencing, but may be easily adapted for any other sequence retrieval method, such as Third-generation sequencing. 

Anonymizer irreversibly replaces personal genomic identifiers with a corresponding sequence of the reference genome. The two-mapping process ensures consistency of mapping locations after depersonalisation.  

Varlock is a tool suitable for storing BAM files in a depersonalized form with the option to restore their original content. The authority managing BAM files can distribute access to specific genomic regions of original BAM file among users, while other parts of the original BAM file stay inaccessible. Now, only the sequence is depersonalized, however, we intend to extend this concept in the future to other BAM file format fields such as CIGAR string, base quality, and optional fields of the alignment section.

# Development
## Testing

run all tests:
<code>python3 -m unittest discover </code>

run single test:
<code>python3 -m unittest tests.<FILE_NAME>.<CLASS_NAME>.<METHOD_NAME></code>

## Examples
run example:
<code>python3 -m examples.bam_mutator</code>

profile example:
<code>python3 -m cProfile -s tottime examples/bam_mutator.py | head -n 140</code>

python3 varlock.py encrypt --key resources/jozko --pub_key resources/jozko.pub --bam resources/input.bam --vac resources/input.vac --out_bam resources/out.mut.bam --out_diff resources/out.diff.enc -v -p password
python3 varlock.py decrypt --key resources/jozko --bam resources/out.mut.bam --diff resources/out.diff.enc --out_bam out.bam
python3 varlock.py reencrypt -d resources/jozko -e resources/jozko.pub -b resources/out.mut.bam -s resources/input.diff -o resources/output.diff -v -p password
python3 varlock.py vac --bam examples/resources/sample.bam --vcf examples/resources/sample.vcf.gz --vac examples/resources/sample.vac


python3 varlock.py encrypt --key tests/resources/varlocker/admin --pub_key tests/resources/varlocker/admin.pub --bam tests/resources/varlocker/encrypt/input.bam --vac tests/resources/varlocker/encrypt/input.vac --out_bam tests/resources/varlocker/encrypt/out.mut.bam --out_diff tests/resources/varlocker/encrypt/out.diff.enc -p password -s seed


## Samtools help
<code>samtools view -H examples/resources/sample.mut.bam</code>

<code>samtools view -h resources/out.mut.bam >> out.mut.sam</code>

<code>less chr22.vcf.gz | egrep -v "^#" | cut -f 4,5 | awk 'length($1) > 2 || length($2) > 2'</code>





