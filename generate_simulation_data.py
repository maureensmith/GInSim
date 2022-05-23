#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 18:426:09 2021

@authors: Maureen Smith, Maria Trofimova
"""

import sys
import argparse
import os
import random
import math
import copy
import numpy as np

from sequence_evolution import sequenceEvol
from output_writer import writer



# initialize argument parser
parser = argparse.ArgumentParser(description='Simulate sequence evolution.')

# Required arguments
parser.add_argument('-file_prefix', required=True,
                    help='file_prefix for fasta file to write the simulated sequences to')

parser.add_argument('-output', required=True,
                    help='output path were the simulated fasta and table files are saved')

parser.add_argument('-L', type=int, required=True,
                    help='sequence length')

parser.add_argument('-N_init', type=int, required=True,
                    help='number of initial sequences')

parser.add_argument('-p_mut', type=float, required=True,
                    help='mutation rate')

parser.add_argument('-t_final', type=int, required=True,
                    help='number of time steps')

parser.add_argument('-p_death', type=float, nargs=1, required=False,
                    help='death rate')


# Optional arguments
parser.add_argument('-p_rep', type=float,
                    help='(optional) replication rate')

parser.add_argument('-p_rep2', type=float,
                    help='(optional) second replication rate')

# Switch
parser.add_argument('--switch_orig', action='store_true', required=False,
					help='(optional) In the case of introductions and two replication rates: '
                         'If true: switch to secon replication rate at the same time point as the original outbreak. '
                         'If false: switch is at the half of the respective time frame.')

parser.add_argument('-fitness', type=int,
					help='(optional) Draw random fitness values for mutated positions for the given amount.')


parser.add_argument('-init_seq', nargs='*',
                    help='(optional) Initial sequence. If not given, a random sequence is created.')

parser.add_argument('-sub_rel', type=float, nargs='*', required=False,
                    help='(optional) list of relative subsampling of the complete sequence set')

parser.add_argument('-sub_abs', type=int, nargs='*', required=False,
                    help='(optional) list of absolute amount of subsampled sequences per day')

parser.add_argument('-intros', type=int, nargs='*', required=False,
                    help='(optional) list of number of introductions')

parser.add_argument('--biasedSampling', action='store_true', required=False,
					help='(optional) When giving relative subsampling, give similar sequences a higher probability to be sampled.')

parser.add_argument('-p_dismiss', type=float,
                    help='(optional) Probability for each shredded part to be dismissed.')

parser.add_argument('-p_shred', type=float,
                    help='(optional) Give probability for each position to be cut.')


args = parser.parse_args()

print("*"*100)
print("Running simulation of sequence evolution with argruments\n")
for a in vars(args):
    print(a, getattr(args, a))
    #print(' {} {}'.format(a, getattr(args, a) or ''))
print("*"*100)


evol_sim = sequenceEvol(length=args.L,
                            p_repl=args.p_rep,
                            p_repl2=args.p_rep2,
                            p_death=args.p_death,
                            p_mut=args.p_mut,
                            N_init=args.N_init,
                            t_start=0,
                            t_final=args.t_final,
                            t_switch=math.floor(args.t_final / 2),
                            fitness=args.fitness
                        )

# run simulation
time_trajectory, newMut_per_t, delSpecies_per_t, fitness = evol_sim.evolve_poi()
#time_trajectory, newMut_per_t, delSpecies_per_t = evol_sim.evolve_gil()

wr = writer(outputpath=args.output, file_prefix=args.file_prefix)

# write initial sequence as reference
wr.write_reference(evol_sim.init_seq)

# always run with 0 introductions
num_intros = [0]
if args.intros is not None:
    num_intros = num_intros + args.intros

print("Intros: ", num_intros)
ts = range(args.t_final)

for ni in num_intros:
    # number of introduction per generation
    num_intros = np.zeros(args.t_final)
    num_intro_seq = np.zeros(args.t_final)
    file_suffix_intro = ""

    # needs to be copied, otherwise the original one is overwritten
    trajectory_withIntroduction = copy.deepcopy(time_trajectory)
    #print("Numer of introductions added: " + str(ni))
    if ni != 0:
        file_suffix_intro="_intros_"+str(ni)
        for t_start in random.choices(range(1, args.t_final), k=ni):
            # random number of copies
            #n_intro = random.randrange(11)
            n_intro = random.randrange(5)+1

            # instead of random sequence, use initial sequence with 3 percent mutation (was 10 before)
            seq = list(evol_sim.init_seq)
            mutPositions  = random.sample(range(args.L), round(args.L*0.03))
            for i in mutPositions:
                seq[i] = sequenceEvol.mutateBase(seq[i])
            evol_sim_intro = sequenceEvol(length=args.L,
                                    p_repl=args.p_rep,
                                    p_repl2=args.p_rep2,
                                    p_mut=args.p_mut,
                                    N_init=n_intro,
                                    t_start=t_start,
                                    #t_final=t_start, # for no replication
                                    t_final=args.t_final,
                                    init_seq="".join(seq)
                                    )

            trajectory_intro, newMut_per_t_intro, delSpecies_per_t_intro = evol_sim_intro.evolve_poi()
            for t in range(t_start, evol_sim_intro.t_final + 1):
                # if not empty
                if len(trajectory_intro) >= (t-t_start)+1 and trajectory_intro[t-t_start]:
                    # add evolving introduction species to initial outbreak
                    for s, ns in trajectory_intro[t-t_start].items():
                        trajectory_withIntroduction[t][s] = trajectory_withIntroduction[t].get(s, 0) + ns
                        # number of new introduced sequences at time t
                        num_intro_seq[t] += ns

            # number of new introduced sequence type at time t
            num_intros[t_start] += 1

    # header=hCoV-19/Italy/LAZ-INMI1-isl/2020|EPI_ISL_410545|2020-01-29
    header_prefix = ">NS|"+file_suffix_intro
    file_suffix = "NS"+file_suffix_intro
    # write fasta with all sequences
    df_NS = wr.write_fasta(file_suffix=file_suffix,
                           header_prefix=header_prefix,
                           species_dict=trajectory_withIntroduction,
                           initialSequence=evol_sim.init_seq)

    df_NS["true_N"] = df_NS["sampled_N"]
    # add intro column
    df_NS["numIntros"] = num_intros
    df_NS["numIntrosSeqs"] = num_intro_seq
    df_NS["delSpecies"] = delSpecies_per_t
    df_NS["newMut"] = newMut_per_t
    df_NS["meanFitness"] = fitness["meanFitness"]
    df_NS["minFitness"] = fitness["minFitness"]
    df_NS["maxFitness"] = fitness["maxFitness"]

    wr.write_table(table=df_NS, file_suffix=file_suffix)
    wr.write_config_yaml(file_suffix=file_suffix)


    #write fasta with shredded reads
    if args.p_shred is not None:
        print("---  Shred and subsample sequence set with cut probability " + str(args.p_shred) + " ---")
        header_prefix = header_prefix = ">NS|shredded|" + str(args.p_shred) + "|" + file_suffix_intro
        file_suffix = "NS_shredded_" + str(args.p_shred) + file_suffix_intro

        df_NS_shred = wr.write_fasta(file_suffix=file_suffix,
                               header_prefix=header_prefix,
                               species_dict=trajectory_withIntroduction,
                               initialSequence=evol_sim.init_seq,
                               p_shred=args.p_shred,
                               p_dismiss=args.p_dismiss)

        df_NS_shred["true_N"] = df_NS["sampled_N"]
        # add intro column
        df_NS_shred["numIntros"] = num_intros
        df_NS_shred["numIntrosSeqs"] = num_intro_seq
        df_NS_shred["delSpecies"] = delSpecies_per_t
        df_NS["newMut"] = newMut_per_t
        df_NS_shred["meanFitness"] = fitness["meanFitness"]
        df_NS_shred["minFitness"] = fitness["minFitness"]
        df_NS_shred["maxFitness"] = fitness["maxFitness"]

        wr.write_table(table=df_NS_shred, file_suffix=file_suffix)
        wr.write_config_yaml(file_suffix=file_suffix)

    # write fasta with all absolute subsample
    if args.sub_abs is not None:
        for s_abs in args.sub_abs:
            print("---  Subsample sequence set taking " + str(s_abs) + " ---")
            header_prefix = header_prefix=">WS|" + str(s_abs) + "|" + file_suffix_intro
            file_suffix = "WSABS_"+ str(s_abs) + file_suffix_intro

            subsampled_time_trajectory = []
            for t in ts:
                num_seq = sum(trajectory_withIntroduction[t].values())
                time_trajectory_sub = {}
                # take abolsute subsample or N(t) if less
                if(num_seq > s_abs):
                    if args.biasedSampling and t > 0:
                        # sequences from the step before
                        biasSequences = subsampled_time_trajectory[t-1]
                        minDiffPerSeq = np.array(len(trajectory_withIntroduction[t])*[args.L])
                        for i, seq in enumerate(trajectory_withIntroduction[t].keys()):
                            for bSeq in biasSequences.keys():
                                #take the minimal difference to one of the sequences of the step before
                                minDiffPerSeq[i] = min(minDiffPerSeq[i],sum(nucl1 != nucl2 for nucl1, nucl2 in zip(bSeq, seq)))
                        p_sample = 1/np.power(2, minDiffPerSeq)
                        # p / sum_i(p_i * n_i)
                        p_sample = p_sample/sum(p_sample*list(trajectory_withIntroduction[t].values()))
                        # create sublists of N_i times item i, and flatten the list
                        allSequences=[]
                        list(map(allSequences.extend, [[it[0]]*it[1] for it in list(trajectory_withIntroduction[t].items())]))

                        allProbabilities=[]
                        list(map(allProbabilities.extend, [[p_sample[i]]*count for i, count in enumerate(trajectory_withIntroduction[t].values())]))

                        seq_subset = np.random.choice(allSequences, p=allProbabilities, replace=False, size=s_abs)
                    else:
                        #sampling the particular sequences for each time point without replacemenet
                        seq_subset = random.sample(list(trajectory_withIntroduction[t].keys()), counts=trajectory_withIntroduction[t].values(), k=s_abs)


                    # sampling the particular sequences for each time point without replacemenet
                    #seq_subset = random.sample(list(trajectory_withIntroduction[t].keys()), counts=trajectory_withIntroduction[t].values(), k=s_abs)
                    # counts sampled sequences
                    for s in seq_subset:
                        time_trajectory_sub[s] = time_trajectory_sub.get(s, 0) + 1
                else:
                    time_trajectory_sub = trajectory_withIntroduction[t]
                subsampled_time_trajectory.append(time_trajectory_sub)

            df_WS_abs = wr.write_fasta(file_suffix=file_suffix, header_prefix=header_prefix, species_dict=subsampled_time_trajectory,initialSequence=evol_sim.init_seq)
            df_WS_abs["true_N"] = df_NS["true_N"]
            df_WS_abs["numIntros"] = num_intros
            df_WS_abs["numIntrosSeqs"] = num_intro_seq
            df_WS_abs["delSpecies"] = delSpecies_per_t
            df_WS_abs["newMut"] = newMut_per_t
            df_WS_abs["meanFitness"] = fitness["meanFitness"]
            df_WS_abs["minFitness"] = fitness["minFitness"]
            df_WS_abs["maxFitness"] = fitness["maxFitness"]

            wr.write_table(table=df_WS_abs, file_suffix=file_suffix)
            wr.write_config_yaml(file_suffix=file_suffix)

    # write fasta with all relative subsample
    if args.sub_rel is not None:
        for s_rel in args.sub_rel:
            print("--- Subsample sequence set taking " + str(s_rel) + " ---")
            header_prefix = header_prefix=">WS|" + str(s_rel) + "|" + file_suffix_intro
            file_suffix = "WSREL_" + str(s_rel) + file_suffix_intro

            # sinoidal
            if s_rel == 0:
                s_sin =[0.5*(1 + np.sin(x/10)) for x in ts]
                rel_sizes=[int(x) for x in round(np.multiply(df_NS["true_N"], s_sin))]
                subsample_size = round(sum(rel_sizes))
                time_subset = random.sample(ts, k=subsample_size, counts=rel_sizes)
            # start with |s_rel| and change to full sample set at half of the time
            elif s_rel < 0:
                s_half = [1] * args.t_final
                s_half[: math.floor(args.t_final / 2)] = [abs(s_rel)] * math.floor(args.t_final / 2)
                rel_sizes=[int(x) for x in round(np.multiply(df_NS["true_N"], s_half))]
                subsample_size = round(sum(rel_sizes))
                time_subset = random.sample(ts, k=subsample_size, counts=rel_sizes)
            else:
                # total sample set size
                subsample_size = round(sum(df_NS["true_N"]) * s_rel)
                # sampling without replacement from which time point the sequences are coming (weighted by the number seqs)
                time_subset = random.sample(ts, k=subsample_size, counts=df_NS["true_N"])
            subsampled_time_trajectory = []
            for t in ts:
                if args.biasedSampling and t > 0:
                    # sequences from the step before
                    biasSequences = subsampled_time_trajectory[t-1]
                    minDiffPerSeq = np.array(len(trajectory_withIntroduction[t])*[args.L])
                    for i, seq in enumerate(trajectory_withIntroduction[t].keys()):
                        for bSeq in biasSequences.keys():
                            #take the minimal difference to one of the sequences of the step before
                            minDiffPerSeq[i] = min(minDiffPerSeq[i],sum(nucl1 != nucl2 for nucl1, nucl2 in zip(bSeq, seq)))
                    p_sample = 1/np.power(2, minDiffPerSeq)
                    # p / sum_i(p_i * n_i)
                    p_sample = p_sample/sum(p_sample*list(trajectory_withIntroduction[t].values()))
                    # create sublists of N_i times item i, and flatten the list
                    allSequences=[]
                    list(map(allSequences.extend, [[it[0]]*it[1] for it in list(trajectory_withIntroduction[t].items())]))

                    allProbabilities=[]
                    list(map(allProbabilities.extend, [[p_sample[i]]*count for i, count in enumerate(trajectory_withIntroduction[t].values())]))

                    seq_subset = np.random.choice(allSequences, p=allProbabilities, replace=False, size=time_subset.count(t))
                else:
                    #sampling the particular sequences for each time point without replacemenet
                    seq_subset = random.sample(list(trajectory_withIntroduction[t].keys()), counts=trajectory_withIntroduction[t].values(), k=time_subset.count(t))

                # counts sampled sequences
                time_trajectory_sub = {}
                for s in seq_subset:
                    time_trajectory_sub[s] = time_trajectory_sub.get(s, 0) + 1
                subsampled_time_trajectory.append(time_trajectory_sub)

            df_WS_rel = wr.write_fasta(file_suffix=file_suffix, header_prefix=header_prefix, species_dict=subsampled_time_trajectory, initialSequence=evol_sim.init_seq)

            df_WS_rel["true_N"] = df_NS["true_N"]
            df_WS_rel["numIntros"] = num_intros
            df_WS_rel["numIntrosSeqs"] = num_intro_seq
            df_WS_rel["delSpecies"] = delSpecies_per_t
            df_WS_rel["newMut"] = newMut_per_t
            df_WS_rel["meanFitness"] = fitness["meanFitness"]
            df_WS_rel["minFitness"] = fitness["minFitness"]
            df_WS_rel["maxFitness"] = fitness["maxFitness"]

            wr.write_table(table=df_WS_rel, file_suffix=file_suffix)
            wr.write_config_yaml(file_suffix=file_suffix)
