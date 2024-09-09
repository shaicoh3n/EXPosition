import pandas as pd
data_path = r""
from pyfaidx import Fasta
# import os
# chrs = Fasta(os.path.join(data_path,'Chromosome','hg38.fa'))
chrs = Fasta(r"./Data/Chromosome/hg38.fa")

def reverse_complement(s: str, complement: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}) -> str:
    '''Performs reverse-complement of a sequence. Default is a DNA sequence.'''
    s_rev = s[::-1]
    lower = [b.islower() for b in list(s_rev)]
    bases = [complement.get(base, base) for base in list(s_rev.upper())]
    rev_compl = ''.join([b.lower() if l else b for l, b in zip(lower, bases)])
    return rev_compl

def read_from_chr(chr_num: str, strand: str, start: int,stop: int) -> str:
    start_zero_based = start-1
    stop_zero_based = stop-1
    seq = chrs['chr' + chr_num][start_zero_based:stop_zero_based].seq.upper()
    indices = list(range(start,stop))
    if strand=='-':
        seq = reverse_complement(seq)
        indices.reverse()
    return seq, indices


# if 1:
#     from pyfaidx import Fasta
#     import os
#     data_path = r"C:\Users\shaic\PycharmProjects\GUI\Data"
#     chrs = Fasta(os.path.join(data_path,'Chromosome','hg38.fa'))
#     chr_y = Fasta(os.path.join(data_path,'Chromosome','Homo_sapiens.GRCh38.dna_sm.chromosome.Y.fa'))
#
# def get_exp_mut_dict(single_site):
#     # debug = 0
#     # if debug:
#     #     ind = single_site['serial']
#     #     if single_site['strand'] == '+':
#     #         single_site['cut_coord'] = int(single_site['start'] + 17)
#     #     elif single_site['strand'] == '-':
#     #         single_site['cut_coord'] = int(single_site['start'] + 2)
#     #     muts_df = pd.read_csv(f'../mut_df_copy_{ind}.csv')
#     # else:
#     #     muts_df = single_site['muts'].copy()
#     muts_df = single_site['muts'].copy()
#     if single_site['strand'] == '+':
#         single_site['cut_coord'] = int(single_site['start'] + 17)
#     elif single_site['strand'] == '-':
#         single_site['cut_coord'] = int(single_site['start'] + 2)
#     cut_coord_zero_base = int(single_site['cut_coord'])
#     cut_coord_one_base = cut_coord_zero_base + 1
#     INS_ins_pos_zero_base = [cut_coord_zero_base - 1 if single_site['strand'] == '+' else cut_coord_zero_base][
#         0]  # cut_coord is the 18th nt in this setup - we want 17th
#     muts_df['del_start_one_base'] = [-1]*muts_df.shape[0]
#     muts_df['del_nts'] = ['']*muts_df.shape[0]
#     muts_df['INS_ins_pos_zero_base'] = [-1]*muts_df.shape[0]
#     muts_df['INS_ins_pos_one_base'] = [-1]*muts_df.shape[0]
#     if 'I' not in muts_df['mut_type'].values:
#         muts_df['ins_nt'] = ['']*muts_df.shape[0]
#
#     for i in muts_df.index:
#         if muts_df.loc[i, 'mut_type'] == 'D':
#             curr_del_start_zero_base = int(muts_df.loc[i, 'del_start'])
#             # elif strand=='-': #before changing in neg cases from start+2-del_start
#             #     curr_del_start_zero_base = int(muts_df.loc[i,'del_start']-muts_df.loc[i,'del_len'])+1
#             muts_df.loc[i, 'del_start_one_base'] = curr_del_start_zero_base + 1
#             curr_del_stop_zero_base = int(curr_del_start_zero_base + muts_df.loc[i, 'del_len'])  # + 1??
#             muts_df.loc[i, 'del_nts'], _ = read_from_chr(single_site['chr'].replace('chr', ''), '+', curr_del_start_zero_base,
#                                                       curr_del_stop_zero_base)
#         elif muts_df.loc[i, 'mut_type'] == 'I':
#             if single_site['strand'] == '-':
#                 muts_df.loc[i, 'ins_nt'] = reverse_complement(muts_df.loc[i, 'ins_nt'])
#             muts_df.loc[i, 'INS_ins_pos_zero_base'] = INS_ins_pos_zero_base
#             muts_df.loc[i, 'INS_ins_pos_one_base'] = muts_df.loc[i, 'INS_ins_pos_zero_base'] + 1
#
#     exp_mut_dict = {'site_chr': single_site['chr'], 'cut_site_pos': cut_coord_one_base,
#                     'site_strand': single_site['strand'], 'seq': single_site['seq'],
#                     'DEL_names': muts_df.loc[muts_df['mut_type'] == 'D', 'indels'].values,
#                     'DEL_prob': single_site['eff'] * muts_df.loc[muts_df['mut_type'] == 'D', 'frequency'].values,
#                     'DEL_del_start': [int(x) for x in
#                                       muts_df.loc[muts_df['mut_type'] == 'D', 'del_start_one_base'].values],
#                     'DEL_del_nts': muts_df.loc[muts_df['mut_type'] == 'D', 'del_nts'].values,
#                     'INS_names': muts_df.loc[muts_df['mut_type'] == 'I', 'indels'].values,
#                     'INS_prob': single_site['eff'] * muts_df.loc[muts_df['mut_type'] == 'I', 'frequency'].values,
#                     'INS_ins_mut': muts_df.loc[muts_df['mut_type'] == 'I', 'ins_nt'].values,
#                     'INS_ins_pos': [int(x) for x in
#                                     muts_df.loc[muts_df['mut_type'] == 'I', 'INS_ins_pos_one_base'].values],
#                     }
#     return exp_mut_dict




# def find_interest_genes(exp_mut_dict, genes_mod_df, feat_type):
#     # get interest genes for _xpresso and then for SpliceAI+TITER and combine the results.
#     if feat_type in ['transcription']:
#         # _Xpresso genes
#         xp_INS_affected_genes = find_affected_genes_model_mut(exp_mut_dict, 'XP', 'INS', genes_mod_df)
#         xp_DEL_affected_genes = find_affected_genes_model_mut(exp_mut_dict, 'XP', 'DEL', genes_mod_df)
#         # xp_SUB_affected_genes = find_affected_genes_model_mut(exp_mut_dict, 'XP', 'SUB', genes_mod_df)
#         xp_SUB_affected_genes = pd.DataFrame()
#         xp_affected_genes = pd.concat([xp_INS_affected_genes, xp_DEL_affected_genes, xp_SUB_affected_genes])
#     else:
#         xp_affected_genes = pd.DataFrame()
#
#     if feat_type in ['splicing', 'initiation','initiation_splicing']:
#         # SAI+TITER genes
#         sai_INS_affected_genes = find_affected_genes_model_mut(exp_mut_dict, 'ST', 'INS', genes_mod_df)
#         sai_DEL_affected_genes = find_affected_genes_model_mut(exp_mut_dict, 'ST', 'DEL', genes_mod_df)
#         # sai_SUB_affected_genes = find_affected_genes_model_mut(exp_mut_dict, 'ST', 'SUB', genes_mod_df)
#         sai_SUB_affected_genes = pd.DataFrame()
#         sai_affected_genes = pd.concat([sai_INS_affected_genes, sai_DEL_affected_genes, sai_SUB_affected_genes])
#     else:
#         sai_affected_genes = pd.DataFrame()
#
#     # have only one occurrence of each gene (even if it was affected by more than one model)
#     affected_genes = pd.concat([xp_affected_genes, sai_affected_genes], ignore_index=True)
#     if affected_genes.shape[0]>0:
#         for unq_gene in pd.unique(affected_genes.gene):
#             curr_gene_df_idxs = affected_genes[affected_genes['gene'] == unq_gene].index
#             if len(curr_gene_df_idxs) > 1:
#                 affected_genes.loc[curr_gene_df_idxs, 'model_mut'] = ','.join(
#                     affected_genes.loc[curr_gene_df_idxs, 'model_mut'].values)
#                 affected_genes = affected_genes.drop(curr_gene_df_idxs[1:])
#     return affected_genes

#
# def find_affected_genes_model_mut(exp_mut_dict, model_name, mut_type, genes_mod_df):
#     '''
#     genes "end" coord for neg strand genes are their start.
#     this function assumes all positions are one based
#     '''
#     chr_num = exp_mut_dict['site_chr']
#     probability = exp_mut_dict[mut_type + '_prob']
#     affected_genes = pd.DataFrame()
#     for i in range(len(probability)):
#         curr_prob = probability[i]
#         if curr_prob > 0:
#             if model_name == 'XP':
#                 upstream_margin = 7000  # from TSS
#                 downstream_margin = 3499  # from TSS
#                 if mut_type == 'DEL' or mut_type == 'SUB':
#                     # print(exp_mut_dict['DEL_names'][i])
#                     mut_range = [exp_mut_dict[mut_type + '_del_start'][i],
#                                  exp_mut_dict[mut_type + '_del_start'][i] + len(
#                                      exp_mut_dict[mut_type + '_del_nts'][i]) - 1]
#                     neg_upstream_margin=upstream_margin
#                     pos_downstream_margin=downstream_margin
#                     pos_gene_start_pos = mut_range[0]
#                     pos_gene_end_pos = mut_range[1]
#                     neg_gene_start_pos = mut_range[0]
#                     neg_gene_end_pos = mut_range[1]
#                 elif mut_type == 'INS':
#                     # print(exp_mut_dict['INS_names'][i])
#                     neg_upstream_margin = upstream_margin-1
#                     pos_downstream_margin = downstream_margin-1
#                     pos_gene_start_pos = exp_mut_dict['INS_ins_pos'][i]
#                     pos_gene_end_pos = exp_mut_dict['INS_ins_pos'][i]
#                     neg_gene_start_pos = exp_mut_dict['INS_ins_pos'][i]
#                     neg_gene_end_pos = exp_mut_dict['INS_ins_pos'][i]
#                 pos_affected_genes = genes_mod_df[(genes_mod_df['chrm'] == chr_num) & (genes_mod_df['strand'] == '+')
#                                                   & (
#                                                           ((genes_mod_df[
#                                                                 'start'] - pos_gene_end_pos <= upstream_margin) & (
#                                                                    pos_gene_end_pos <= genes_mod_df['start']))
#                                                           | ((pos_gene_start_pos - genes_mod_df[
#                                                       'start'] <= pos_downstream_margin) & (
#                                                                      genes_mod_df['start'] <= pos_gene_start_pos))
#                                                   )]
#                 neg_affected_genes = genes_mod_df[(genes_mod_df['chrm'] == chr_num) & (genes_mod_df['strand'] == '-')
#                                                   & (
#                                                           ((genes_mod_df[
#                                                                 'end'] - neg_gene_end_pos <= downstream_margin) & (
#                                                                    neg_gene_end_pos <= genes_mod_df['end']))
#                                                           | ((neg_gene_start_pos - genes_mod_df[
#                                                       'end'] <= neg_upstream_margin) & (
#                                                                      genes_mod_df['end'] <= neg_gene_start_pos))
#                                                   )]
#                 curr_affected_genes = pd.concat([pos_affected_genes, neg_affected_genes])
#             elif model_name == 'ST':
#                 upstream_margin = 0
#                 downstream_margin = 0
#                 if mut_type == 'INS':
#                     # print(exp_mut_dict['INS_names'][i])
#                     start_pos = exp_mut_dict['INS_ins_pos'][i]
#                     end_pos = exp_mut_dict['INS_ins_pos'][i]
#                     downstream_margin+=1
#                 else:
#                     # print(exp_mut_dict['DEL_names'][i])
#                     start_pos = exp_mut_dict[mut_type + '_del_start'][i]
#                     end_pos = start_pos + len(exp_mut_dict[mut_type + '_del_nts'][i]) - 1
#                 curr_affected_genes = genes_mod_df[(genes_mod_df['chrm'] == chr_num) &
#                                                    (((start_pos - genes_mod_df['start'] >= upstream_margin) &
#                                                      (genes_mod_df['end'] - start_pos >= downstream_margin)) |
#                                                     ((end_pos - genes_mod_df['start'] >= upstream_margin) &
#                                                      (genes_mod_df['end'] - end_pos >= downstream_margin))
#                                                     )].copy(deep=True)
#         else:
#             curr_affected_genes = pd.DataFrame()
#         curr_affected_genes['model_mut'] = model_name + '_' + mut_type + '_' + str(i)
#         # print(curr_affected_genes)
#         affected_genes = pd.concat([affected_genes, curr_affected_genes])
#     return affected_genes




# # A1BG		-	chr19		58345178	58353492
# # A2ML1		+	chr12		8822621	    8887001

# print(find_interest_genes(ST_exp_gene_neg, genes, 'splicing'))
# b=1

# genes = pd.read_csv(r"C:\Users\shaic\PycharmProjects\GitHub\EXPosition_6_11_22\code\Data\all_genome_genes.csv")

#
# # A1BG		-	chr19		58345178	58353492
# ST_exp_gene_neg = {'site_chr': 'chr19', 'cut_site_pos': 9000000, 'site_strand': '-',
#                 'DEL_names': ['DEL 58353493','DEL 58353493_plus_1',
#                               'DEL 58353492','DEL 58353491_3',
#                               'DEL 58345179','DEL 58345178',
#                               'DEL 58345177','DEL 58345178_1','DEL 58345175_4',
#                               ],
#                 'DEL_prob': [1]*9,
#                 'DEL_del_start': [58353493, 58353493,
#                                   58353492, 58353491,
#                                   58345179, 58345178,
#                                   58345177, 58345178,
#                                   58345175],
#                 'DEL_del_nts': ['A','AA','A','AAA','A','A','A','AA','AAAA'],
#                 'INS_names': ['INS 58353493','INS 58353493_plus_1',
#                               'INS 58353492','INS 58353491_3',
#                               'INS 58345179','INS 58345178',
#                               'INS 58345177','INS 58345178_1','INS 58345175_3',
#                               ],
#                 'INS_prob': [1]*9,
#                 'INS_ins_mut': ['A','AA','A','AAA','A','A','A','AA','AAA'],
#                 'INS_ins_pos': [58353493, 58353493,
#                                   58353492, 58353491,
#                                   58345179, 58345178,
#                                   58345177, 58345178,
#                                   58345175]
#                 }

#
# # # A2ML1		+	chr12		8822621	    8887001
# ST_exp_gene_pos = {'site_chr': 'chr12', 'cut_site_pos': 9000000, 'site_strand': '+',
#                 'DEL_names': ['DEL 8822620_2nts','DEL 8822620','DEL 8822621',
#                               'DEL 8822622','DEL 8822721',
#                               'DEL 8887000','DEL 8887001',
#                               'DEL 8887002','DEL 8887001_2nts'
#                               ],
#                 'DEL_prob': [1]*9,
#                 'DEL_del_start': [8822621-1,8822621-1,8822621,
#                                 8822621+1,8822621+100,
#                                 8887001-1,8887001,
#                                 8887001+1,8887001
#                               ],
#                 'DEL_del_nts': ['AA','A','AA','A','AA','A','A','A','AA'],
#                 'INS_names': ['INS 8822620','INS 8822621',
#                               'INS 8822622','INS 8822721',
#                               'INS 8887000','INS 8887001',
#                               'INS 8887002'
#                               ],
#                 'INS_prob': [1]*7,
#                 'INS_ins_mut': ['A']*7,
#                 'INS_ins_pos': [8822621-1,8822621,
#                                 8822621+1,8822621+100,
#                                 8887001-1,8887001,
#                                 8887001+1
#                               ]
#                 }

# # # A2ML1		+	chr12		8822621	    8887001
# XP_exp_gene_pos = {'site_chr': 'chr12', 'cut_site_pos': 9000000, 'site_strand': '+',
#                 'DEL_names': ['DEL 8822621_min_7001','DEL 8822621_min_7001_0',
#                               'DEL 8822621_min_7000','DEL 8822621_min_7000_6999',
#                               'DEL 8822621_min_3000','DEL 8822621',
#                               'DEL 8822621_plus_3000','DEL 8822621_plus_3498','DEL 8822621_plus_3499',
#                               'DEL 8822621_plus_3498_3500','DEL 8822621_plus_3500'
#                               ],
#                 'DEL_prob': [1]*11,
#                 'DEL_del_start': [8822621-7001,8822621-7001,
#                               8822621-7000,8822621-7000,
#                               8822621-3000,8822621-0,
#                               8822621+3000,8822621+3498,8822621+3499,
#                               8822621+3498,8822621+3500
#                               ],
#                 'DEL_del_nts': ['A','AA','A','AA','A','A','A','A','A','AAA','A'],
#                 'INS_names': ['INS 8822621_min_7001','INS 8822621_min_7001_2nts',
#                               'INS 8822621_min_7000','INS 8822621_min_7000_2nts',
#                               'INS 8822621_min_3000','INS 8822621',
#                               'INS 8822621_plus_3000','INS 8822621_plus_3498','INS 8822621_plus_3499',
#                               'INS 8822621_plus_3498_3nts','INS 8822621_plus_3500'
#                               ],
#                 'INS_prob': [1]*11,
#                 'INS_ins_mut': ['A']*11,
#                 'INS_ins_pos': [8822621-7001,8822621-7001,
#                               8822621-7000,8822621-7000,
#                               8822621-3000,8822621-0,
#                               8822621+3000,8822621+3498,8822621+3499,
#                               8822621+3498,8822621+3500
#                               ]
#                 }
#
# # A1BG		-	chr19		58345178	58353492
# XP_exp_gene_neg = {'site_chr': 'chr19', 'cut_site_pos': 9000000, 'site_strand': '-',
#                 'DEL_names': ['DEL 58353492_plus_7001','DEL 58353492_plus_7001_0',
#                               'DEL 58353492_plus_7000','DEL 58353492_plus_7000_6999',
#                               'DEL 58353492_plus_3000','DEL 58353492',
#                               'DEL 58353492_min_3000','DEL 58353492_min_3498','DEL 58353492_min_3499',
#                               'DEL 58353492_min_3498_3500','DEL 58353492_min_3500'
#                               ],
#                 'DEL_prob': [1]*11,
#                 'DEL_del_start': [58353492+7001,58353492+7000,
#                               58353492+7000,58353492+6999,
#                               58353492+3000,58353492+0,
#                               58353492-3000,58353492-3498,58353492-3499,
#                               58353492-3500,58353492-3500
#                               ],
#                 'DEL_del_nts': ['A','AA','A','AA','A','A','A','A','A','AAA','A'],
#                                 'INS_names': ['INS 58353492_plus_7001','INS 58353492_plus_7000_2nts',
#                               'INS 58353492_plus_7000','INS 58353492_plus_6999_2nts',
#                               'INS 58353492_plus_3000','INS 58353492',
#                               'INS 58353492_min_3000','INS 58353492_min_3498','INS 58353492_min_3499',
#                               'INS 58353492_min_3500_3nts','INS 58353492_min_3500'
#                               ],
#                 'INS_prob': [1]*11,
#                 'INS_ins_mut': ['A']*11,
#                 'INS_ins_pos': [58353492+7001,58353492+7000,
#                               58353492+7000,58353492+6999,
#                               58353492+3000,58353492+0,
#                               58353492-3000,58353492-3498,58353492-3499,
#                               58353492-3500,58353492-3500
#                               ]
#                 }