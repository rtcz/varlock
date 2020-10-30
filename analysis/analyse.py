import io
import os
import sys
from os.path import dirname, abspath

import matplotlib_venn as venn
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

# enable absolute import from this package
sys.path.append(dirname(dirname(abspath(__file__))))

from src.vac import Vac, pysam
from src.fasta_index import FastaIndex

sns.set_style('ticks')


def vcf2df(filename: str) -> pd.DataFrame:
    with open(filename, 'r') as file:
        lines = [line for line in file if not line.startswith('##')]

    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={
            '#CHROM': str,
            'POS': int,
            'ID': str,
            'REF': str,
            'ALT': str,
            'QUAL': float,
            'FILTER': str,
            'INFO': str
        },
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def vac2df(filename: str) -> pd.DataFrame:
    new_filename = os.path.splitext(filename)[0] + '.txt'
    Vac.vac2text(filename, new_filename)

    return pd.read_csv(new_filename, sep='\t', skiprows=2, index_col=0, names=['ref_id', 'counts'])


def countstr2altfreq(string: str, ref_id: int) -> float:
    counts = list(map(int, string.split(',')))
    return 1 - counts[ref_id] / sum(counts)


def depth(format: str, genotype: str) -> int:
    data = dict(zip(format.split(':'), genotype.split(':')))
    return int(data['DP'])


def altfreq(format: str, genotype: str) -> float:
    data = dict(zip(format.split(':'), genotype.split(':')))
    return sum(map(float, data['VAF'].split(',')))


class Analyser:
    def __init__(self, personal_vcf: str, masked_vcf: str, vac: str, out_dir: str):
        assert os.path.isdir(out_dir)
        self._out_dir = out_dir

        self._index = FastaIndex.from_vcf(pysam.VariantFile(personal_vcf))

        self._personal_df = vcf2df(personal_vcf)
        self._personal_df = self._personal_df[self._personal_df['FILTER'] == 'PASS']
        self._personal_df.index = self._personal_df.apply(
            lambda row: self._index.pos2index(row['CHROM'], row['POS'] - 1), axis=1
        )
        self._personal_df['depth'] = self._personal_df.apply(
            lambda row: depth(row['FORMAT'], row[self._personal_df.columns.get_loc("FORMAT") + 1]), axis=1
        )
        self._personal_df['alt_freq'] = self._personal_df.apply(
            lambda row: altfreq(row['FORMAT'], row[self._personal_df.columns.get_loc("FORMAT") + 1]), axis=1
        )
        self._personal_df = self._personal_df[(self._personal_df['depth'] > 30) & (self._personal_df['QUAL'] > 30)]

        self._masked_df = vcf2df(masked_vcf)
        self._masked_df = self._masked_df[self._masked_df['FILTER'] == 'PASS']
        self._masked_df.index = self._masked_df.apply(
            lambda row: self._index.pos2index(row['CHROM'], row['POS'] - 1), axis=1
        )
        self._masked_df['depth'] = self._masked_df.apply(
            lambda row: depth(row['FORMAT'], row[self._masked_df.columns.get_loc("FORMAT") + 1]), axis=1
        )
        self._masked_df['alt_freq'] = self._masked_df.apply(
            lambda row: altfreq(row['FORMAT'], row[self._masked_df.columns.get_loc("FORMAT") + 1]), axis=1
        )
        self._masked_df = self._masked_df[(self._masked_df['depth'] > 30) & (self._masked_df['QUAL'] > 30)]

        self._population_df = vac2df(vac)
        self._population_df['alt_freq'] = self._population_df.apply(
            lambda row: countstr2altfreq(row['counts'], int(row['ref_id'])), axis=1)

        self._personal_ids = set(self._personal_df.index)
        self._masked_ids = set(self._masked_df.index)
        self._public_ids = set(self._population_df.index)

        self._masked_set = self._public_ids \
            .intersection(self._personal_ids) \
            .difference(self._masked_ids)

        self._not_masked_set = self._public_ids \
            .intersection(self._personal_ids) \
            .intersection(self._masked_ids)

        self._introduced_set = self._public_ids \
            .intersection(self._masked_ids) \
            .difference(self._personal_ids)

        self._not_covered_set = self._personal_ids \
            .intersection(self._masked_ids) \
            .difference(self._public_ids)

        self._not_found_set = self._public_ids \
            .difference(self._personal_ids) \
            .difference(self._masked_ids)

        self._not_covered_hidden_set = self._personal_ids \
            .difference(self._public_ids) \
            .difference(self._masked_ids)

        self._not_covered_revealed_set = self._masked_ids \
            .difference(self._public_ids) \
            .difference(self._personal_ids)

    def sets(self):
        union = len(self._public_ids.union(self._personal_ids).union(self._masked_ids))
        print('A ∪ B ∪ C %d' % union)
        print('A (population variants) %d' % len(self._public_ids))
        print('B (personal variants) %d' % len(self._personal_ids))
        print('C (masked variants) %d' % len(self._masked_ids))
        print()

        value = len(self._public_ids.intersection(self._personal_ids).intersection(self._masked_ids))
        print('A ∩ B ∩ C %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._public_ids.intersection(self._personal_ids).difference(self._masked_ids))
        print('A ∩ B - C %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._personal_ids.intersection(self._masked_ids).difference(self._public_ids))
        print('B ∩ C - A %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._public_ids.intersection(self._masked_ids).difference(self._personal_ids))
        print('A ∩ C - B %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._public_ids.difference(self._masked_ids).difference(self._personal_ids))
        print('A - C - B %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._personal_ids.difference(self._public_ids).difference(self._masked_ids))
        print('B - A - C %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._masked_ids.difference(self._public_ids).difference(self._personal_ids))
        print('C - A - B %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._public_ids.intersection(self._personal_ids))
        print('A ∩ B %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._personal_ids.intersection(self._masked_ids))
        print('B ∩ C %d (%.2f%%)' % (value, value / union * 100))

        value = len(self._public_ids.intersection(self._masked_ids))
        print('A ∩ C %d (%.2f%%)' % (value, value / union * 100))

        plt.figure()
        venn.venn3_unweighted(
            [self._public_ids, self._personal_ids, self._masked_ids],
            set_labels=('A', 'B', 'C')
        )
        plt.savefig(os.path.join(self._out_dir, 'sets.png'))

    def vcfs(self):
        step_size = 0.1

        plt.figure()
        plt.title('Allele frequency by VCF')

        plt.hist(
            [
                self._population_df['alt_freq'].values,
                self._personal_df['alt_freq'].values,
                self._masked_df['alt_freq'].values
            ],
            bins=np.arange(0, 1.0 + step_size, step_size),
            label=['population VCF', 'personal VCF', 'masked VCF'],
            density=False
        )
        plt.yscale('log')
        plt.legend()
        plt.xlabel('allele frequency')
        plt.ylabel('occurrence')
        plt.xticks(np.arange(0, 1.0 + step_size, step_size))
        plt.savefig(os.path.join(self._out_dir, 'vcfs.png'), bbox_inches='tight', dpi=300)

    @DeprecationWarning
    def analyse_vcf(self):
        step_size = 0.1

        plt.figure()
        plt.title('Public allele frequency')

        plt.hist(
            [
                self._population_df[
                    (self._population_df[1] <= 1.0)
                    & (self._population_df[0].isin(self._personal_ids))
                    ]
                [1].values,
            ],
            bins=np.arange(0, 1.0 + step_size, step_size),
            label=['personal VCF'],
            density=False
        )
        plt.legend()
        plt.xlabel('allele frequency')
        plt.ylabel('occurrence')
        plt.xticks(np.arange(0, 1.0 + step_size, step_size))
        plt.savefig(self.OUT_DIR + 'vcf.png', bbox_inches='tight')

    def pop(self):
        step_size = 0.1

        plt.figure()
        plt.title('Allele Frequency by Category')

        result = plt.hist(
            [
                self._population_df[self._population_df.index.isin(self._masked_set)]['alt_freq'].values,
                self._population_df[self._population_df.index.isin(self._introduced_set)]['alt_freq'].values,
                self._population_df[self._population_df.index.isin(self._not_masked_set)]['alt_freq'].values,
            ],
            bins=np.arange(0, 1.0 + step_size, step_size),
            label=['masked', 'introduced', 'not masked'],
            density=False
        )
        # plt.yscale('log')
        plt.legend()
        plt.xlabel('population allele frequency')
        plt.ylabel('occurrence')
        plt.xticks(np.arange(0, 1.0 + step_size, step_size))
        plt.savefig(os.path.join(self._out_dir, 'categorized_af.png'), bbox_inches='tight', dpi=300)

        # vs
        masked_y = result[0][0]
        introduced_y = result[0][1]
        not_masked_y = result[0][2]

        plt.figure()
        plt.title('Masked vs Not Masked')
        width = step_size / 2
        bar1 = plt.bar(result[1][:-1] + (step_size / 2), not_masked_y / (not_masked_y + masked_y), width=width)
        bar2 = plt.bar(result[1][:-1] + (step_size / 2), masked_y / (not_masked_y + masked_y), width=width,
                       bottom=not_masked_y / (not_masked_y + masked_y))
        plt.legend((bar1, bar2), ['not masked', 'masked'])
        plt.xlabel('population allele frequency')
        plt.ylabel('ratio')
        plt.ylim((0, 1))
        plt.savefig(os.path.join(self._out_dir, 'masked_notmasked.png'), bbox_inches='tight', dpi=300)

        plt.figure()
        plt.title('Not Masked vs Introduced')
        width = step_size / 2
        bar1 = plt.bar(result[1][:-1] + (step_size / 2), not_masked_y / (not_masked_y + introduced_y), width=width)
        bar2 = plt.bar(result[1][:-1] + (step_size / 2), introduced_y / (not_masked_y + introduced_y), width=width,
                       bottom=not_masked_y / (not_masked_y + introduced_y))
        plt.legend((bar1, bar2), ['not masked', 'introduced'])
        plt.xlabel('population allele frequency')
        plt.ylabel('ratio')
        plt.ylim((0, 1))
        plt.savefig(os.path.join(self._out_dir, 'notmasked_introduced.png'), bbox_inches='tight', dpi=300)

        plt.figure()
        plt.title('Masked vs Introduced')
        width = step_size / 2
        bar1 = plt.bar(result[1][:-1] + (step_size / 2), masked_y / (masked_y + introduced_y), width=width)
        bar2 = plt.bar(result[1][:-1] + (step_size / 2), introduced_y / (masked_y + introduced_y), width=width,
                       bottom=masked_y / (masked_y + introduced_y))
        plt.legend((bar1, bar2), ['masked', 'introduced'])
        plt.xlabel('population allele frequency')
        plt.ylabel('ratio')
        plt.ylim((0, 1))
        plt.savefig(os.path.join(self._out_dir, 'masked_introduced.png'), bbox_inches='tight', dpi=300)

    def not_masked(self):
        # remove duplicates
        not_masked_personal_df = self._personal_df[self._personal_df.index.isin(self._not_masked_set)]

        not_masked_masked_df = self._masked_df[self._masked_df.index.isin(self._not_masked_set)]

        # inner join on position column
        joined_df = not_masked_personal_df.join(
            not_masked_masked_df,
            how='inner',
            lsuffix='_personal',
            rsuffix='_masked'
        )

        print(f'not_masked total len: {len(joined_df)}')
        joined_df = joined_df[joined_df['ALT_personal'] == joined_df['ALT_masked']]
        print(f'not_masked same alt len: {len(joined_df)}')

        # altered: 13 (0.29%) from 4416

        altered_df = joined_df[joined_df['alt_freq_personal'] != joined_df['alt_freq_masked']]  # type: pd.DataFrame
        not_altered_df = joined_df[joined_df['alt_freq_personal'] == joined_df['alt_freq_masked']]  # type: pd.DataFrame

        altered_set = set(altered_df.index)
        not_altered_set = set(joined_df.index).difference(altered_set)

        # all except one of altered differ only in AF
        masked_percentage = len(altered_set) / len(joined_df) * 100

        print('altered: %d (%.2f%%) from %d' % (len(altered_set), masked_percentage, len(joined_df)))
        # altered: 1463 (33.23%) from 4403

        # print(altered_df.iloc[3].values)
        # test_df = self._varlock_df[
        #     self._varlock_df['pos'].isin(altered_set)
        #     # & (self._varlock_df['hetero2homo'] == 1)
        # ]
        # print(test_df)
        # print(self._not_masked_set.difference(set(test_df['pos'].values)))
        # print(altered_set.difference(set(test_df['pos'].values)))
        # exit(0)

        step_size = 0.1

        altered_population_df = self._population_df[self._population_df.index.isin(altered_set)]
        not_altered_population_df = self._population_df[self._population_df.index.isin(not_altered_set)]

        plt.figure()
        plt.title('not masked position')
        result = plt.hist(
            [
                altered_population_df['alt_freq'].values,
                not_altered_population_df['alt_freq'].values,
            ],
            bins=np.arange(0, 1.0 + step_size, step_size),
            label=['masked allele', 'not masked allele'],
            density=False
        )
        # plt.yscale('log')
        plt.legend()
        plt.xlabel('population allele frequency')
        plt.ylabel('occurrence')
        plt.xticks(np.arange(0, 1.0 + step_size, step_size))
        plt.savefig(os.path.join(self._out_dir, 'not_masked_a.png'), bbox_inches='tight', dpi=300)

        altered_y = result[0][0]
        not_altered_y = result[0][1]

        plt.figure()
        plt.title('Altered vs Not Altered')
        width = step_size / 2
        bar1 = plt.bar(result[1][:-1] + (step_size / 2), altered_y / (altered_y + not_altered_y), width=width)
        bar2 = plt.bar(result[1][:-1] + (step_size / 2), not_altered_y / (altered_y + not_altered_y), width=width,
                       bottom=altered_y / (altered_y + not_altered_y))
        plt.legend((bar1, bar2), ['altered', 'not altered'])
        plt.xlabel('population allele frequency')
        plt.ylabel('ratio')
        plt.ylim((0, 1))
        plt.savefig(os.path.join(self._out_dir, 'not_masked_b.png'), bbox_inches='tight', dpi=300)

    @DeprecationWarning
    def analyse_rare(self):
        rare_not_masked_df = self._population_df[
            (self._population_df[0].isin(self._not_masked_set))
            & (self._population_df[1] < 0.1)
            ]
        rare_not_masked_set = set(rare_not_masked_df[0].values)

        # print('rare (<0.1) not_masked len %d' % len(rare_not_masked_set))

        rare_df = self._varlock_df[
            (self._varlock_df['pos'].isin(rare_not_masked_set))
            # & (self._varlock_df['is_mut'] == 0)
        ]

        # missing_set = rare_not_masked_set.difference(set(rare_df['pos'].values))
        # print(missing_set)
        print(len(rare_df))
        print(rare_df)

        # step_size = 0.01
        #
        # plt.figure()
        # plt.title('Rare Not Masked')
        #
        # plt.hist(
        #     [
        #         vardict_not_masked_df[1].values,
        #     ],
        #     bins=np.arange(0, 0.1 + step_size, step_size),
        #     label=['value'],
        #     density=False
        # )
        # # plt.yscale('log')
        # plt.legend()
        # plt.xlabel('allele frequency')
        # plt.ylabel('occurrence')
        # plt.xticks(np.arange(0, 0.1 + step_size, step_size))

    @DeprecationWarning
    def analyse_varlock(self):
        homo2homo_df = self._varlock_df[
            (self._varlock_df['homo2homo'] == 1)
            & (self._varlock_df['is_mut'] == 1)
            ]
        homo2hetero_df = self._varlock_df[
            (self._varlock_df['homo2hetero'] == 1)
            & (self._varlock_df['is_mut'] == 1)
            ]
        hetero2homo_df = self._varlock_df[
            (self._varlock_df['hetero2homo'] == 1)
            & (self._varlock_df['is_mut'] == 1)
            ]
        hetero2hetero_df = self._varlock_df[
            (self._varlock_df['hetero2hetero'] == 1)
            & (self._varlock_df['is_mut'] == 1)
            ]

        print('variants total %d' % len(self._varlock_df))
        print('homo2homo %d' % len(homo2homo_df))
        print('homo2hetero %d' % len(homo2hetero_df))
        print('hetero2homo %d' % len(hetero2homo_df))
        print('hetero2hetero %d' % len(hetero2hetero_df))

    @DeprecationWarning
    def analyse_not_covered(self):
        not_covered_df = self._personal_df[self._personal_df[3].isin(self._not_covered_set)]
        print(not_covered_df)
        print()

    @DeprecationWarning
    def analyse_not_covered_hidden(self):
        hidden_df = self._personal_df[self._personal_df[3].isin(self._not_covered_hidden_set)]
        print(hidden_df)
        print()

    @DeprecationWarning
    def analyse_not_covered_revealed(self):
        hidden_df = self._masked_df[self._masked_df[3].isin(self._not_covered_revealed_set)]
        print(len(hidden_df))
        print(hidden_df.values)

    @staticmethod
    def show():
        plt.show()


if __name__ == '__main__':
    # Analyser.extract_vcf_pos()

    analyser = Analyser(
        personal_vcf='/data/projects/exome/variant/grch38_decoy_alt-one/original/one_ZinaBednarikova.vcf',
        masked_vcf='/data/projects/varlock/variant/grch38_decoy_alt-one/original_masked_gnomad3nfe_onepanel/'
                   'one_ZinaBednarikova.deepvariant.vcf',
        vac='/data/projects/varlock/vac/gnomad3nfe_pass_sorted_multisnp_onepanel.vac',
        out_dir='/data/projects/varlock/analysis/gnomad3nfe_pass_deep'
    )
    analyser.sets()
    analyser.vcfs()
    analyser.pop()
    analyser.not_masked()
