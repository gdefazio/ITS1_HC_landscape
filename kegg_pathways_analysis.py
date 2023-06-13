import pandas as pd
import numpy as np
import wget
import os
from functools import reduce


def kegg2genus(genus_comp_path):
    '''
    This function interrogates KEGG database API, downloading pathways of
    species belonging to genera found in samples. Then, for a reliable genus
    representation of KEGG pathways, only pathways in common among species
    of the same genus were reported.
    The result of this function is a DataFrame reporting KEGG pathways as rows
    and genera as columns and filled by 0/1 as boolean.
    '''
    a = pd.read_csv('./KEGG_organisms_list.tsv', '\t', header=None)
    fungi = a[[True if 'Fungi' in el else False for el in a[3]]]
    protista = a[[True if 'Protists' in el else False for el in a[3]]]
    organisms = pd.concat([fungi, protista], axis=0)
    organisms.columns = ['accession', 'abbreviation', 'name', 'path']

    organisms['genus'] = organisms.name.apply(lambda x: x.split(' ')[0])
    samples = pd.read_csv(genus_comp_path, index_col=0)

    genera = organisms.genus.unique()

    c = 0
    cn = 0
    for g in samples.genus:
        if g in genera:
            c += 1
        else:
            cn += 1
            print(g)
    print('Genera in KEGG are %s and not in KEGG are %s' % (c, cn))
    selgens = organisms[organisms.genus.isin(samples.genus)]
    for acc in selgens.abbreviation:
        genus = selgens[selgens.abbreviation == acc].genus.iloc[0]
        try:
            os.mkdir(genus)
        except FileExistsError:
            pass
        wget.download('https://rest.kegg.jp/list/pathway/%s' % acc, out='%s/%s.kegg' % (
            genus, selgens[selgens.abbreviation == acc].name.iloc[0].replace(' ', '_')))

    alldf = list()
    for genus in os.listdir('./KEGG_selected_genera/'):
        long = os.listdir(os.path.join('./KEGG_selected_genera/', genus))
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
    final['percent'] = final.sum(axis=1) / 60 * 100
    final.to_csv('./4genusKEGGpathway.csv')
    return final




    p4s = pd.DataFrame(index=range(0, 137))
    p4s['npath'] = final['npath']
    for p in [el for el in samples.columns if el.endswith('_count')]:
        finalT = final.T
        finalT['idx'] = list(finalT.index)
        p4s[p] = finalT[finalT.idx.isin(samples[samples[p] != 0]['genus'].unique())][range(0, 137)].sum(axis=0)
        iterp4s = (p4s[p4s.columns[1:-1]] > 0).iterrows()
        p4s['percentage'] = [el[1].value_counts()[True] / 60 * 100 for el in iterp4s]

if __name__ == '__main__':

    gcp = '../samples/BioMas_raref_genus_collapsed.csv'
    kegg2genus(genus_comp_path=gcp)




    ##### somma delle frequenze relative dei generi che esprimono determinati KEGG ####

    genus4kegg = pd.read_csv('./4genusKEGGpathway.csv', index_col=0)
    samples = pd.read_csv('../samples/BioMas_raref_genus_collapsed.csv', index_col=0)
    newdf = pd.DataFrame(columns = ['npath', *samples.columns[1:-1]])
    newdf.npath = genus4kegg.npath
    c = 0
    for p in genus4kegg['npath']:
        for s in samples.columns[1:-1]:
            sub = genus4kegg[genus4kegg['npath'] == p].T
            cols = list(sub[sub[c] != 0].index)[1:-1]
            newdf.loc[newdf.npath == p, s] = sum(samples[samples['genus'].isin(cols)][s])/100
            #print(p,s,len(cols),samples[samples['genus'].isin(cols)][s].shape)
        c += 1
