import pandas as pd
import numpy as np
import wget
import os
from functools import reduce
import argparse
import argcomplete


def split_options():
    parser = argparse.ArgumentParser(
        description="Kegg pathway predictor",
        prefix_chars="--")

    parser.add_argument("-s", "--samples", type=str,
                        help="ITS1 classification counts",
                        action="store", required=True,
                        default=None)
    parser.add_argument("-l", "--list", type=str,
                        help="KEGG organism list",
                        action="store", required=True,
                        default=None)
    parser.add_argument("-r", "--resume",
                        help="Resume after kegg download",
                        action="store_true", required=False,
                        default=False)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def kegg2genus(genus_comp_path, kegg_organism_list, resume=False):
    '''
    This function interrogates KEGG database API, downloading pathways of
    species belonging to genera found in samples. Then, for a reliable genus
    representation of KEGG pathways, only pathways in common among species
    of the same genus were reported.
    The result of this function is a DataFrame reporting KEGG pathways as rows
    and genera as columns and filled by 0/1 as boolean.
    '''
    a = pd.read_csv(kegg_organism_list, sep='\t', header=None)
    fungi = a[[True if 'Fungi' in el else False for el in a[3]]]
    protista = a[[True if 'Protists' in el else False for el in a[3]]]
    organisms = pd.concat([fungi, protista], axis=0)
    organisms.columns = ['accession', 'abbreviation', 'name', 'path']

    organisms['genus'] = organisms.name.apply(lambda x: x.split(' ')[0])
    samples = pd.read_csv(genus_comp_path, index_col=0)

    genera = organisms.genus.unique()

    outdir = 'KEGG_selected_genera'
    if resume is False:
        c = 0
        cn = 0
        for g in samples.genus:
            if g in genera:
                c += 1
            else:
                cn += 1
                print(g)
        print('Genera in KEGG are %s and not in KEGG are %s' % (c, cn))
        os.mkdir(outdir)
        selgens = organisms[organisms.genus.isin(samples.genus)]
        for acc in selgens.abbreviation:
            genus = selgens[selgens.abbreviation == acc].genus.iloc[0]
            try:
                os.mkdir(os.path.join(outdir,genus))
            except FileExistsError:
                pass
            wget.download('https://rest.kegg.jp/list/pathway/%s' % acc, out='%s/%s/%s.kegg' % (
                outdir,
                genus,
                selgens[selgens.abbreviation == acc].name.iloc[0].replace(' ', '_')))

    alldf = list()
    for genus in [ el for el in os.listdir(outdir) if os.path.isdir(os.path.join(outdir,el)) ]:
        long = os.listdir(os.path.join(outdir, genus))
        if len(long) == 1:
            new = pd.read_csv(os.path.join('./KEGG_selected_genera/', genus, long[0]), sep='\t', header=None)
            new['npath'] = new[1].apply(lambda x: x.split(' - ')[0])
            new = new[['npath']]
            new.drop_duplicates(inplace=True)
            new[['npath']].to_csv(os.path.join('./KEGG_selected_genera/', "%s.csv" % genus))
            new[genus] = 1
            alldf.append(new)
        else:
            dfl = list()
            for elem in long:
                new = pd.read_csv(os.path.join('./KEGG_selected_genera/', genus, elem), sep='\t', header=None)
                new['npath'] = new[1].apply(lambda x: x.split(' - ')[0])
                dfl.append(new[['npath']])
            reduced = reduce(lambda x, y: x.merge(y, on='npath'), dfl)
            reduced.drop_duplicates(inplace=True)
            reduced.to_csv(os.path.join('./KEGG_selected_genera/', "%s.csv" % genus))
            reduced[genus] = 1
            alldf.append(reduced)

    final = reduce(lambda x, y: x.merge(y, on='npath', how='outer'), alldf)
    final.replace(np.NaN, 0, inplace=True)
    npath = final['npath']
    final = final[final.columns[1:]]
    final.mask(final > 0, 1, inplace=True)
    final['percent'] = final.sum(axis=1) / final.shape[1] * 100
    final['npath'] = npath
    final.to_csv('./KEGG4samples_bool.csv', quoting=False)
    # p4s = pd.DataFrame(index=range(0, 137))
    # p4s['npath'] = final['npath']
    # for p in [el for el in samples.columns if el.endswith('_count')]:
    #     finalT = final.T
    #     finalT['idx'] = list(finalT.index)
    #     p4s[p] = finalT[finalT.idx.isin(samples[samples[p] != 0]['genus'].unique())][range(0, 137)].sum(axis=0)
    #     iterp4s = (p4s[p4s.columns[1:-1]] > 0).iterrows()
    #     p4s['percentage'] = [el[1].value_counts()[True] / 60 * 100 for el in iterp4s]
    return final, samples


def gen_rel_abund4kegg2samples(genus4kegg, samples):
    '''
    To obtain a DataFrame with rows representing KEGG paths,
    columns representing samples and filled by the sum of
    relative abundance of genera supporting each path in each sample.
    '''
    # fix: bisogna normalizzare rispetto al totale di copertura del campione
    # ovvero la somma totale di reads per i generi coperti dall'analisi kegg

    newdf = pd.DataFrame(columns=['npath', *samples.columns[1:-1]])
    newdf.npath = genus4kegg.npath
    c = 0
    all_supporting = [el for el in genus4kegg.columns if el not in ('percent', 'npath')]
    for p in genus4kegg['npath']:
        for s in samples.columns[1:-1]:
            sub = genus4kegg[genus4kegg['npath'] == p].T
            # print(sub)
            cols = [el for el in sub[sub[c] != 0].index if el not in ('percent', 'npath')]
            # print(cols)
            a = sum(samples[samples['genus'].isin(cols)][s])
            b = sum(samples[samples['genus'].isin(all_supporting)][s])
            newdf.loc[newdf.npath == p, s] = a / b * 100
            # print(newdf)
            # print(p,s,len(cols),samples[samples['genus'].isin(cols)][s].shape)
        c += 1
    newdf.to_csv('./KEGG4samples_rel_freq.csv',quoting=False)


if __name__ == '__main__':
    options = split_options()
    k2g_table, gcp_table = kegg2genus(
        genus_comp_path=options.samples,
        kegg_organism_list=options.list,
        resume=options.resume
    )
    gen_rel_abund4kegg2samples(
        genus4kegg=k2g_table,
        samples=gcp_table
    )