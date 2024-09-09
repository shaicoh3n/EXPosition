from time import time
import os
import numpy as np
import pandas as pd
import sys
import re
# import tensorflow as tf
# from tensorflow import keras
# from mimetypes import guess_type
# from optparse import OptionParser
from keras.models import load_model
# from Bio import SeqIO
# from functools import partial

from Models.Oncosplice.genomic.variant_utils import generate_mut_variant
from EXPosition_utils import read_from_chr, reverse_complement

xp_model_path = os.path.join('.','Models','Xpresso','humanMedian_trainepoch.11-0.426.h5')
model = load_model(xp_model_path)

def pred_transcription_for_mut_gene(interest_gene_df, exp_mut_dict, XP_muts):
    left_margin = 7000  # from TSS
    right_margin = 3499  # form TSS
    gene, transcripts, transcripts_start, transcripts_end, gene_strand = interest_gene_df[['gene','transcripts', 'transcripts_start', 'transcripts_end', 'strand']]
    transcripts = transcripts.split(',')
    transcripts_start = [int(x) for x in transcripts_start.split('/')]
    transcripts_end = [int(x) for x in transcripts_end.split('/')]
    max_transcript_score = -1
    max_tss_no_p_mut_scores = {x: 0 for x in XP_muts}
    checked_tsss = []
    max_transcript_name = 'None'
    if exp_mut_dict['ensembl_id'] != 'No Interest Transcript':
        interest_transcript_index = transcripts.index(exp_mut_dict['ensembl_id'])
        transcripts = [transcripts[interest_transcript_index]]
        transcripts_start = [transcripts_start[interest_transcript_index]]
        transcripts_end = [transcripts_end[interest_transcript_index]]
    for i_transc in range(len(transcripts)):
        curr_transcript = transcripts[i_transc]
        transc_start = transcripts_start[i_transc]
        transc_end = transcripts_end[i_transc]
        rev = gene_strand == '-'
        tss_pos, l_margin, r_margin, range_start, range_end = \
            xpresso_get_tss_margin_range_info(transc_start, transc_end, rev, left_margin, right_margin)
        if tss_pos in checked_tsss: #we already checked that tss so no need to rerun xpresso with the same mutations
            continue
        checked_tsss.append(tss_pos)
        fasta_headers, fasta_seqs, all_same = \
            xpresso_generate_gene_seqs(gene, gene_strand, tss_pos,
                                       exp_mut_dict, XP_muts)
        orig_tss_mut_scores = {x: [] for x in XP_muts}
        orig_tss_no_p_mut_scores = {x: 0 for x in XP_muts}
        if not all_same:
            xp_pred_df = xp_predict_seqs(fasta_headers, fasta_seqs)
            xp_pred_df['score'] = [10 ** float(x) for x in xp_pred_df['RAW_SCORE'].to_list()]
            ref_score = xp_pred_df.loc[xp_pred_df['ID'] == gene + ':Ref_0', 'score'].values[0]
            orig_tss_mut_scores = xp_get_score_for_mut_type('INS', xp_pred_df, orig_tss_mut_scores)
            orig_tss_mut_scores = xp_get_score_for_mut_type('DEL', xp_pred_df, orig_tss_mut_scores)
            orig_tss_mut_scores = xp_get_score_for_mut_type('SUB', xp_pred_df, orig_tss_mut_scores)
        else: #all mutations resulted in the same sequence i.e the reference seq
            continue

        orig_abs_score = 0
        for mut_name in orig_tss_mut_scores.keys():
            raw_mut_name = mut_name[3:6]
            i_prob = int(mut_name[-1])
            orig_abs_score += exp_mut_dict[raw_mut_name+ '_prob'][i_prob] * abs(orig_tss_mut_scores[mut_name][0]-ref_score) / ref_score
            orig_tss_no_p_mut_scores[mut_name] = abs(orig_tss_mut_scores[mut_name][0]-ref_score) / ref_score
        if orig_abs_score>max_transcript_score:
            max_transcript_score = orig_abs_score
            max_tss_no_p_mut_scores = orig_tss_no_p_mut_scores
            max_transcript_name = curr_transcript

    if max_transcript_score==-1: #this should happen since if we got to this function, some transcript was affected by some mutation
        max_transcript_score=0
    #normalize with clinvar max result
    # orig_abs_score = min(1,orig_abs_score/3.55)
    return max_transcript_score, max_tss_no_p_mut_scores, max_transcript_name

def xp_get_score_for_mut_type(mut_type,xp_pred_df,mut_scores):
    ids = xp_pred_df['ID'].values
    scores = xp_pred_df['score'].values
    mut_comp = re.compile(f'{mut_type}\s\d')

    #1st idx is location in xpred_df, 2nd item is the name of the mutation ,3rd idx is the location of the probability in the exp_mut_dict[mut_type+'prob']
    mut_xp_pred_df_idxs = [[i_xp_pred, mut_comp.findall(ids[i_xp_pred])[0] ,int(mut_comp.findall(ids[i_xp_pred])[0][4:])] for i_xp_pred in range(len(ids)) if mut_comp.findall(ids[i_xp_pred])]
    if len(mut_xp_pred_df_idxs)==0:
        b=1
    for i_xp_pred,mut_name,i_prob in mut_xp_pred_df_idxs:
        mut_scores['XP_' + '_'.join(mut_name.split(' '))].append(scores[i_xp_pred])
        # mut_score += exp_mut_dict[mut_type+'_prob'][i_prob] * curr_score
    return mut_scores

def xp_one_hot(seq):
    num_seqs = len(seq)
    seq_len = len(seq[0])
    seqindex = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    seq_vec = np.zeros((num_seqs,seq_len,4), dtype='bool')
    for i in range(num_seqs):
        thisseq = seq[i]
        for j in range(seq_len):
            try:
                seq_vec[i,j,seqindex[thisseq[j]]] = 1
            except:
                pass
    return seq_vec

def xp_predict_seqs(fasta_headers,fasta_seqs):
    i, bs, names, predictions, sequences = 0, 2, [], [], [] #reduce batch size bs if takes too much memory
    halflifedata = np.zeros((bs,6), dtype='float32')
    for sequence,name in zip(fasta_seqs,fasta_headers):
        sequences.append(sequence)
        names.append(name)
        i += 1
        if (len(sequence) != 10500):
            sys.exit( "Error in sequence %s, length is not equal to the required 10,500 nts. Please fix or pad with Ns if necessary." % name )
        if i % bs == 0:
            seq = xp_one_hot(sequences)
            predictions.extend(model.predict([seq, halflifedata], batch_size=bs,verbose=0).flatten().tolist())
            sequences = []

    remain = i % bs
    if remain > 0:
        halflifedata = np.zeros((remain,6), dtype='float32')
        seq = xp_one_hot(sequences)
        predictions.extend( model.predict([seq, halflifedata], batch_size=remain).flatten().tolist() )

    df = pd.DataFrame(np.column_stack((names, predictions)), columns=['ID','RAW_SCORE'])
    return df

def xpresso_get_tss_margin_range_info(transc_start: int, transc_end: int, rev: bool, left_margin: int, right_margin: int ) -> tuple:
        '''
        Returns the left and right margin (based on rev) and tarnscript range.
        '''
        if rev:
            tss_pos = transc_end
            l_margin, r_margin = right_margin, left_margin
        else:
            tss_pos = transc_start
            l_margin, r_margin = left_margin, right_margin

        range_start, range_end = tss_pos-l_margin, tss_pos+r_margin
        return tss_pos, l_margin, r_margin, range_start, range_end

def xpresso_generate_gene_seqs(gene: str, gene_strand: str, tss_pos: int,
                               exp_mut_dict: dict, curr_models_muts: list):
    '''
    This function gets a dataframe of mutation of a SINGLE gene and generates the
    corresponding Fasta file for Xpresso.
    '''
    l_margin=7000
    r_margin=3499
    fasta_seqs, fasta_headers = [], []
    all_same = True
    #generate start+end positions
    start_position_one_base = np.zeros(len(curr_models_muts),dtype=int)
    end_position_one_base = np.zeros(len(curr_models_muts),dtype=int)
    del_nts = []
    ins_nts = []
    mut_type = []
    for i_model_mut in range(len(curr_models_muts)):
        curr_model, curr_mut, curr_mut_idx = curr_models_muts[i_model_mut].split('_')
        curr_mut_idx = int(curr_mut_idx)
        if curr_mut=='SUB':
            start_position_one_base[i_model_mut] = exp_mut_dict['SUB_del_start'][curr_mut_idx]
            end_position_one_base[i_model_mut] = start_position_one_base[i_model_mut] + len(exp_mut_dict['SUB_del_nts'][curr_mut_idx]) - 1
            del_nts.append(exp_mut_dict['SUB_del_nts'][curr_mut_idx])
            ins_nts.append(exp_mut_dict['SUB_ins_mut'][curr_mut_idx])
            mut_type.append(exp_mut_dict['SUB_names'][curr_mut_idx] + ' ' + curr_mut + ' ' + str(curr_mut_idx))
        if curr_mut=='DEL':
            start_position_one_base[i_model_mut] = exp_mut_dict['DEL_del_start'][curr_mut_idx]
            end_position_one_base[i_model_mut] = start_position_one_base[i_model_mut] + len(exp_mut_dict['DEL_del_nts'][curr_mut_idx]) - 1
            del_nts.append(exp_mut_dict['DEL_del_nts'][curr_mut_idx])
            ins_nts.append('')
            mut_type.append(exp_mut_dict['DEL_names'][curr_mut_idx] + ' ' + curr_mut + ' ' + str(curr_mut_idx))
        if curr_mut=='INS':
            start_position_one_base[i_model_mut] = exp_mut_dict['INS_ins_pos'][curr_mut_idx]
            end_position_one_base[i_model_mut] = start_position_one_base[i_model_mut]
            del_nts.append('')
            ins_nts.append(exp_mut_dict['INS_ins_mut'][curr_mut_idx])
            mut_type.append(exp_mut_dict['INS_names'][curr_mut_idx] + ' ' + curr_mut + ' ' + str(curr_mut_idx))

    min_ok_dist_to_tss = 40 #if a gene was mutated up to 40 nts from its tss, then we try all combinations (for now, -10 to +10)
    #relative to gene strand
    if gene_strand=='+':
        mut_loc_rel_tss = [['downstream to tss' if x >= tss_pos else 'upstream to tss'][0] for x in end_position_one_base]
    else:
        mut_loc_rel_tss = [['downstream to tss' if x <= tss_pos else 'upstream to tss'][0] for x in end_position_one_base]
    # mut_close_to_tss = np.logical_or(np.abs(end_position_one_base - tss_pos) <= min_ok_dist_to_tss,
    #                                  np.abs(start_position_one_base - tss_pos) <= min_ok_dist_to_tss)
    num_times_check_around_del_tss = 80
    max_margin = num_times_check_around_del_tss + max(end_position_one_base-start_position_one_base)
    # max_margin=0

    r_margin += max_margin
    l_margin += max_margin

    extract_start = int(tss_pos - [l_margin if gene_strand == '+' else r_margin][0])
    extract_end = extract_start+10500+2*max_margin

    ref_p_seq, ref_p_indices = read_from_chr(exp_mut_dict['site_chr'].replace('chr',''), '+', extract_start, extract_end)
    # ref_n_seq = reverse_complement(ref_p_seq)
    # ref_n_indices = ref_p_indices
    # ref_n_indices.reverse()

    ref_seq = ref_p_seq
    if gene_strand=='-':
        ref_seq = reverse_complement(ref_p_seq)

    ref_seq_10500 = ref_seq[max_margin:-max_margin].upper()
    # ref_10500_indices = ref_p_indices[max_margin:-max_margin]
    # if gene_strand=='-':
    #     ref_10500_indices.reverse()

    fasta_seqs.append(ref_seq_10500)
    fasta_headers.append(f"{gene}:Ref_0")

    # from pyfaidx import Fasta
    # pos = Fasta(r"C:\Users\shaic\Desktop\AGRN_gene_pos_strand.fna")
    # pos_orig = pos['NC_000001.11:1020120-1056116'][:].seq.upper()
    # neg = Fasta(r"C:\Users\shaic\Desktop\A1BG_gene_neg_strand.fna")
    # neg_orig = neg['NC_000019.10:c58353492-58345183'][:].seq.upper()

    # mutated sequences. Each mutation generates a mutated sequence
    for i in range(len(start_position_one_base)):
        if mut_type[i][0] == 'D':
            curr_mut_type = 'DEL'
        elif mut_type[i][0] == 'I':
            curr_mut_type = 'INS'
        elif mut_type[i][0] == 'S':
            curr_mut_type = 'SUB'
        mut_indices = ref_p_indices
        mut_seq = ref_p_seq
        mut_seq, mut_indices, mut_affected_gene, _ = generate_mut_variant(seq=mut_seq, indices=mut_indices,
                                                          start_pos=start_position_one_base[i],
                                                          ref=del_nts[i], mut=ins_nts[i],
                                                          var_type=curr_mut_type)
        if not mut_affected_gene:
            fasta_seqs.append(ref_seq_10500)
            fasta_headers.append(f"{gene}:{start_position_one_base[i]}:{end_position_one_base[i]}:{del_nts[i]}:{ins_nts[i]}:{mut_type[i]}:Mut no effect:Mut")
            continue

        if gene_strand == '-':
            mut_seq = reverse_complement(mut_seq)
            mut_indices.reverse()

        delta_xs = [0]
        # if mut_close_to_tss[i]:
        #     delta_xs = [x for x in range(-num_times_check_around_del_tss//2, num_times_check_around_del_tss//2 + 1)]
        # else:
        #     delta_xs = [0]

        if mut_loc_rel_tss[i]=='downstream to tss':
            trunc_w_muts = True
        else:
            trunc_w_muts = False

        for dx in delta_xs:
            if curr_mut_type == 'INS':
                if trunc_w_muts: # (+,down)
                    mut_seq = mut_seq[max_margin + dx:]
                else: # (+,up)
                    mut_seq = mut_seq[len(ins_nts[i]) + max_margin + dx:]
            elif curr_mut_type == 'DEL':
                if trunc_w_muts: #(+,down)
                    mut_seq = mut_seq[max_margin+dx:]
                else: # (+,up)
                    mut_seq = mut_seq[max_margin-len(del_nts[i]) + dx:]
            elif curr_mut_type =='SUB':
                if trunc_w_muts:
                    mut_seq = mut_seq[max_margin+dx:]
                else:
                    mut_seq = mut_seq[(len(ins_nts[i])+max_margin-len(del_nts[i])) + dx:]

            mut_seq = mut_seq[:10500]
            if len(mut_seq)!=10500:
                print('seq isn''t 10500 nts long! for mutations within {} to {}'.format(start_position_one_base[i],end_position_one_base[i]))
            if all_same and mut_seq != ref_seq_10500:
                all_same = False
            # if gene_strand=='+':
            #     orig = pos_orig
            # else:
            #     orig = neg_orig
            # try:
            #     print(i)
            #     print(exp_mut_dict['target_seq'])
            #     print(f'site strand {exp_mut_dict["site_strand"]}')
            #     print(mut_type[i])
            #     print(mut_loc_rel_tss[i])
            #     print(f'gene strand {gene_strand}')
            #     q=mut_seq.index(orig[:30])
            #     if q==7000:
            #         print('tss pos GOOD')
            #     else:
            #         print(f'BAD q is {q}') #not 7000
            # except:
            #     print('BAD tss not found')
            # try:
            #     if mut_type[i][0]=='I':
            #         debug_ins_nts = ins_nts[i]
            #         if exp_mut_dict['site_strand']=='-':
            #             debug_ins_nts = reverse_complement(debug_ins_nts)
            #         target_seq = exp_mut_dict['target_seq'][:17] + debug_ins_nts + exp_mut_dict['target_seq'][17:]
            #     elif mut_type[i][0]=='D':
            #         import re
            #         del_len, del_start, _ = [int(x) for x in re.findall("-?\d+", mut_type[i])]
            #         target_seq = exp_mut_dict['target_seq'][:17+del_start] + exp_mut_dict['target_seq'][17+del_start+del_len:]
            #
            #     if gene_strand!=exp_mut_dict['site_strand']:
            #         target_seq = reverse_complement(target_seq)
            #     print(f'GOOD target seq found in mut: {mut_seq.index(target_seq)}')
            #     if mut_seq.upper()==ref_seq_10500:
            #         print('BAD mut == ref')
            # except:
            #     print('BAD target seq not found')
            fasta_seqs.append(mut_seq.upper())
            # <gene>:<start>:<end>:<Reference_Allele>:<Tumor_Seq_Allele2>:<Variant_Type>:<Variant_Classification>:Mut
            fasta_headers.append(f"{gene}:{start_position_one_base[i]}:{end_position_one_base[i]}:{del_nts[i]}:{ins_nts[i]}:{mut_type[i]}:{dx}:Mut")

    return fasta_headers, fasta_seqs, all_same

# def add_max_score_to_res_df(ref_score, curr_gene_scores, var_type, curr_site_gene_res_df):
#     if len(curr_gene_scores) != 0:
#         max_mut_score = max(curr_gene_scores)
#         max_mut_delta = (max_mut_score - ref_score) / ref_score
#         curr_site_gene_res_df[var_type] = [f"{100 * max_mut_delta:.4f}%"]
#     return curr_site_gene_res_df

#
# data_path = r"/tamir1/shaicohen1_share/Data/"
# #data_path = r"C:\Users\shaic\PycharmProjects\GitHub\EXPosition-GUI\Data"
# # import pickle
# # with open('shaked_sites_transcription.pkl', 'rb') as f:
# #     sites = pickle.load(f)
# sites = pd.read_csv(r"C:\Users\shaic\PycharmProjects\GitHub\EXPosition_6_11_22\code\new_format_sites_2.csv")
# out_sites = pred_transcription(sites, data_path,1,'transcription')
# n=6