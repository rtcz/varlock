# Varlock: the tool for pseudonymization of sequenced genome

## Varlock
Varlock is a command line interface utility for reversible allele maskig contained within a BAM file. The differences between the original and the masked BAM are stored in encrypted form, so only authorized users can recover their permitted parts of the original BAM. Varlock uses RSA asymmetric encryption with pairs of public and private keys for encryption and decryption, respectively.

Two new binary file types are introduced in the tool: the first is a population specific _Variant Allele Count_ (_VAC_) file, which is essentially a compact VCF file that represents standard allele frequencies of a population. The second is a BAM file specific _BAM difference_ (_BDIFF_) file, where the differences between the original and a depersonalized BAM are compactly stored. 

Varlock provides commands for encryption, decryption, granting access, and conversion from VCF files to VAC files. 

### Variant Allele Count (VAC) file format
VAC is a compact binary file format for storing SNV and INDEL allele frequencies among the population. VAC file contains a header with a number of stored SNVs and INDELs. Each record represents one variation in a VCF file and contains its absolute position in the genome and observed allele frequencies in the population. Each SNV record contains four numbers as frequencies of A,T,G,C bases respectively. INDEL record features a map of frequencies and corresponding DNA sequences. SNV and INDEL records are stored separately in two blocks and are sorted by their genomic position.

### BAM Difference (BDIFF) file format 
SNV and INDEL differences between original and depersonalized BAM file are stored in a binary BDIFF file format. The BDIFF file format consists of a header, a file index, and a number of BDIFF records. The file header contains number of SNVs and INDELs stored in the file and metadata required for decryption process. Secondly, the file index enables fast seeking of genomic positions. 

A single BDIFF record represents a difference at a single position in the genome. Since BAM file usually contains multiple alignments covering one genomic position (where they can differ) and individual can have different alleles of the same gene, the BAM difference is stored as a map representing bijection from set of original alleles to set of depersonalized ones and vice versa. In case of an SNV, the map is a bijection from the set of four DNA letters to its permutation. In case of an INDEL, the map is a bijection from a set of known INDEL sequences at this genomic position to its permutation.

BDIFF file acts as the secret key in the context of Varlock encryption, thus it is actually never stored on a medium in plain form. BDIFF file is encrypted symmetrically with the AES encryption using a randomly generated key. Only the AES key is then encrypted asymmetrically by the user's public key and stored in the BDIFF header. 

### Masking
User must provide the original BAM file, a VAC file of the closest population, and his public key for encryption. Both VAC and BAM files are iterated at the same time to write a new masked BAM together with BDIFF tracking all the changes. Every time a genomic position of a VAC record intersects with one or more alignments from the BAM file, a pseudo random permutation of DNA base letters (in the case of SNV) or INDEL sequences is generated. Alignments that do not overlap with any VAC records are written to the masked BAM unchanged. The masking always involves the whole content of the supplied BAM, meaning that all former alignments are present in the masked BAM either altered or unaltered.

# Development
## Testing

run all tests:
<code>python3 -m unittest discover </code>

run single test:
<code>python3 -m unittest tests.<FILE_NAME>.<CLASS_NAME>.<METHOD_NAME></code>

## Examples
python3 varlock.py encrypt --key resources/jozko --pub_key resources/jozko.pub --bam resources/input.bam --vac resources/input.vac --out_bam resources/out.mut.bam --out_diff resources/out.diff.enc -v -p password
python3 varlock.py decrypt --key resources/jozko --bam resources/out.mut.bam --diff resources/out.diff.enc --out_bam out.bam
python3 varlock.py reencrypt -d resources/jozko -e resources/jozko.pub -b resources/out.mut.bam -s resources/input.diff -o resources/output.diff -v -p password
python3 varlock.py vac --bam examples/resources/sample.bam --vcf examples/resources/sample.vcf.gz --vac examples/resources/sample.vac
python3 varlock.py encrypt --key tests/resources/varlocker/admin --pub_key tests/resources/varlocker/admin.pub --bam tests/resources/varlocker/encrypt/input.bam --vac tests/resources/varlocker/encrypt/input.vac --out_bam tests/resources/varlocker/encrypt/out.mut.bam --out_diff tests/resources/varlocker/encrypt/out.diff.enc -p password -s seed


