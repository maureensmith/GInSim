#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import random
import math
from copy import deepcopy


class sequenceEvol:
	'''
	Each instance of this class represens one evolutionary trajectory run. Starting with a sequence set of equal sequences,
	growing and mutating over time.
	'''

	def __init__(self, length, p_mut, N_init, t_start, t_final, p_repl=None, p_repl2=None, rand_repl=None, p_death=None, t_switch=None, init_seq=None, fitness=None):
		self.L = length
		self.p_repl = p_repl
		self.p_repl2 = p_repl2
		self.rand_repl = rand_repl
		self.p_death = p_death
		self.p_mut = p_mut
		self.N = N_init
		self.t_start = t_start
		self.t_final = t_final
		self.t_switch = t_switch
		self.init_seq = self._initSeq(init_seq)
		self.fitness_per_pos = self._drawFitness(length, fitness)

	def _drawFitness(self, length, f=None):
		if f is None:
			return length*[1.0]
		else:
			fitness = length*[1.0]
			# only first half of sequence has an effect if mutated (should have same effect as if i draw affecting indices)
			#fitness[0:round((length/2)-1)] = np.random.lognormal(mean=0.0, sigma=1.0, size=round(length/2))
			#fitness[0:f-1] = np.random.lognormal(mean=0.0, sigma=np.sqrt(0.2/f), size=f)
			fitness[0:f-1] = np.random.lognormal(mean=0.0, sigma=np.sqrt(0.2/f), size=f)

			# for i in range(round((length/2)-1)):
			# 	r = 0
			# 	#while(r < 1.0 or r > 2.0):
			# 	while(r > 1.0):
			# 		r = np.random.lognormal(mean=0.0, sigma=0.2, size=1)[0]
			# 	fitness[i] = r
			print("fitness: ", fitness)
			return fitness

	def _initSeq(self, init_seq=None):
		'''
		Create initial sequence as array of each nucleotide. If no sequence is given, create a random one.
		'''
		if init_seq is None:
			alphabet = 'ACGT'
			init_seq = ''.join(random.SystemRandom().choice(alphabet)for _ in range(self.L))
		return (init_seq)

	def evolve_poi(self):
		'''
		Mutate initial sequence set over the course of t generations
		'''
		print("--- Run simulation ---")

		# start time
		t = self.t_start+1

		# the current replication rate of the time step
		curr_p_repl = self.p_repl

		curr_N = self.N

		#list of dictionaries per time point, containing each sequence species an their amount
		species_dict_per_t = []
		species_dict_per_t.append({self.init_seq: curr_N})



		delSpecies_per_t = np.zeros(self.t_final)
		delSpecies = 0
		#delSpecies_per_t.append(deepcopy(delSpecies))

		newMut_per_t = np.zeros(self.t_final)
		newMut = 0

		meanFitness_per_t = np.ones(self.t_final)
		minFitness_per_t = np.ones(self.t_final)
		maxFitness_per_t = np.ones(self.t_final)


		new_set = [self.init_seq]
		while t < self.t_final:
			#print("p_repl2 = ", self.p_repl2)
			print("Current time t = ", t)

			curr_p_repl = self._get_repl(t)
			print("Current time repl = ", curr_p_repl)
			# if self.rand_repl is None and self.p_repl is None:
			# 	# TODO check, which one was publication??
			# 	curr_p_repl = math.sin(t * 0.11) / 15 + 1.03
			# 	curr_p_repl = math.sin(t * 0.11) / 25 + 1.01
			# 	#curr_p_repl = math.sin(t * 0.1) / 20 + 1.05
			# elif(self.p_repl2 is not None):
			# 	# Switch the replication rate somewhere in the middle if that is the correct mode (p_repl2!=0)
			# 	if (t >= self.t_switch) and (self.p_repl2 > 0):
			# 		curr_p_repl = self.p_repl2
			# 	#if  repl2 is negative, switch only in one timestep and switch back to repl
			# 	elif(self.p_repl2 < 0):
			# 		if(t == self.t_switch):
			# 			curr_p_repl=abs(self.p_repl2)
			# 		else:
			# 			curr_p_repl=self.p_repl

			#print("Current replication rate = ", curr_p_repl)

			# draw number of sequences of the next generation
			curr_N = np.random.poisson(curr_p_repl * curr_N)
			print("Current N = ", curr_N)

			#print("Number of new sequences: ", curr_N)

			# collect sequence species for t
			species_dict = {}

			if curr_N > 0:
				# calculate fitness per species
				fitness_per_species = [np.prod([self.fitness_per_pos[i] for i in range(len(self.init_seq)) if self.init_seq[i] != seq[i]]) for seq in new_set]
				meanFitness_per_t[t]= np.mean(fitness_per_species)
				minFitness_per_t[t] = min(fitness_per_species)
				maxFitness_per_t[t] = max(fitness_per_species)

				#print(fitness_per_species)
				new_set = random.choices(new_set, k=curr_N, weights=fitness_per_species)
				#new_set = random.choices(new_set, k=curr_N)


				# draw number of mutations in new generation
				# if number of new mutation sizes is larger than sites (rare/not likely) take number of sites
				# (= 2 mutations on one site, ignore back mutation)
				num_mut_sites = min(np.random.poisson(self.p_mut * curr_N * self.L), curr_N * self.L)
				newMut_per_t[t] = num_mut_sites
				#t("Number of mutating sites: ",num_mut_sites)

				if num_mut_sites > 0:
					#for i in mut_sites:
					for _ in range(num_mut_sites):
						i = random.randrange(curr_N * self.L)
						seq_ind = math.floor(i / self.L)
						seq_pos = i % self.L
						seq = list(new_set[seq_ind])
						# seq[seq_pos] = self._mutateBase(new_set[seq_ind][seq_pos],initSeq[seq_pos])
						seq[seq_pos] = self.mutateBase(new_set[seq_ind][seq_pos])
						s = "".join(seq)
						new_set[seq_ind] = s

				# count species
				for s in new_set:
					species_dict[s] = species_dict.get(s, 0) + 1

			# add the sequences to the map
			species_dict_per_t.append(species_dict)
			t += 1

		return species_dict_per_t, newMut_per_t, delSpecies_per_t, \
			{"meanFitness":meanFitness_per_t, "minFitness":minFitness_per_t, "maxFitness":maxFitness_per_t}

	def evolve_gil(self):
		'''
		Birth-death-process and mutation of the given population
		'''
		print("--- Run simulation with gillespie---")
		# start time
		t = self.t_start

		# collect all sequences for a certain "day" buckets of this size
		t_bucket=1

		# the current replication rate of the time step
		curr_p_repl = self._get_repl(t)

		curr_N = self.N

		#death rate
		delta=0.5
		alpha=curr_p_repl*delta

		#birth rate (for now birth and death is for all species equivalent)
		#alpha=curr_p_repl-delta

		#list of dictionaries per time point, containing each sequence species an their amount
		species_dict_per_t = []
		species_dict = {self.init_seq: curr_N}
		species_dict_per_t.append(deepcopy(species_dict))

		delSpecies_per_t = np.zeros(self.t_final)
		delSpecies = 0
		#delSpecies_per_t.append(deepcopy(delSpecies))

		newMut_per_t = np.zeros(self.t_final)
		newMut = 0
		#newMut_per_t.append(deepcopy(newMut))
		meanFitness_per_t = np.ones(self.t_final)
		minFitness_per_t = np.ones(self.t_final)
		maxFitness_per_t = np.ones(self.t_final)


		curr_t = t+t_bucket
		while curr_t < self.t_final:
			# new day new sequences, keep species dict, as it is
			if(t>=curr_t):
				# add the sequences to the list
				print("NewMut: ",newMut)
				print("delSpecies: ",delSpecies)
				print("Current t:",curr_t)
				print("Number of previous sequences: ", sum(species_dict_per_t[curr_t-1].values()))
				print("Curr repl: ", curr_p_repl, " ", self._get_repl(curr_t-1))
				print("Number of new sequences: ", sum(species_dict.values()))
				print("Number of expected sequences: ", sum(species_dict_per_t[curr_t-1].values())*self._get_repl(curr_t-1))
				delSpecies = 0
				newMut = 0
				newMut_per_t[curr_t] = newMut
				delSpecies_per_t[curr_t] = delSpecies

				fitness_per_species = [np.prod([self.fitness_per_pos[i] for i in range(len(self.init_seq)) if self.init_seq[i] != seq[i]]) for seq in species_dict.keys()]
				meanFitness_per_t[curr_t]= np.mean(fitness_per_species*species_dict.values())
				minFitness_per_t[curr_t] = min(fitness_per_species)
				maxFitness_per_t[curr_t] = max(fitness_per_species)
				species_dict_per_t.append(deepcopy(species_dict))
				curr_t += t_bucket

			if(curr_N != 0):
				num_species = len(species_dict)

				# have the reactions for each species (can be adjusted to fitness)
				#react = num_species*[alpha, beta]

				# stoichiometric matrix,
				#update_matrix = [[0 for x in range(num_species)] for x in range(num_species*2)]
				#actually list is enough, because there is only one species involved
				# update_vector = [0 for x in range(num_species*2)]
				# for i in range(num_species):
				# 	# birth rate
				# 	#update_matrix[2*i][i] = 1
				# 	update_vector[2*i] = 1
				# 	#death rate
				# 	#update_matrix[2*i+1][i] = -1
				# 	update_vector[2*i+1] = -1
				update_vector = [1, -1]

				u1 = random.random()
				u2 = random.random()

				fitness_per_species = [np.prod([self.fitness_per_pos[i] for i in range(len(self.init_seq)) if self.init_seq[i] != seq[i]]) for seq in species_dict.keys()]
				#print(fitness_per_species)

				# reaction propensities here: N_i*a and N_i*d
				#a_m = np.repeat(list(species_dict.values()), 2)*(num_species*[alpha, delta])
				a_m = [alpha*curr_N, delta*curr_N]
				#print("Repl :", curr_p_repl, " a ", alpha, " b ", delta)
				#print(num_species)
				#print(fitness_per_species)
				#print(np.sum(a_m))
				#print(np.sum(a_blub))

				#time to next reaction TODO needs to be adjusted! is for constant rates
				tau = (1/sum(a_m)*np.log(1/u1))
				t += tau
				# get replication rate for new time point
				curr_p_repl = self._get_repl(t)

				#adjust alpha accordingly
				alpha=curr_p_repl*delta
				a_m = [alpha*curr_N, delta*curr_N]
				#print("Repl :", curr_p_repl, " a ", alpha, " b ", delta)
				#for i in range(num_species):
				#	a_m[2*i] *= fitness_per_species[i]/np.sum(fitness_per_species)

				# choose reaction
				r_idx = next(x for x, val in enumerate(np.cumsum([i/sum(a_m) for i in a_m])) if val > u2)
				#print(r_idx)
				# always two reactions per species
				#spec_idx = r_idx//2
				#spec_seq = list(species_dict.keys())[spec_idx]

				spec_prob = np.array([i/curr_N for i in species_dict.values()])
				#print(self.fitness_per_pos)
				#print(fitness_per_species)
				#spec_prob = np.array(species_dict.values())/curr_N

				if(r_idx==0):
					weighted_N = list(species_dict.values())*np.array(fitness_per_species)
					#print(weighted_N)
					spec_prob = weighted_N/sum(weighted_N)

				spec_seq = np.random.choice(list(species_dict.keys()), p=spec_prob)


				species_dict[spec_seq] = species_dict[spec_seq] + update_vector[r_idx]
				# delete species, if the last one died
				if species_dict[spec_seq] < 1:
					del species_dict[spec_seq]
					delSpecies += 1

				curr_N = sum(species_dict.values())

				if curr_N > 0:
					# draw number of mutations in new fraction of generation
					# if number of new mutation sizes is larger than sites (rare/not likely) take number of sites
					# (= 2 mutations on one site, ignore back mutation)
					nm = np.random.poisson(self.p_mut * curr_N * self.L * tau)
					num_mut_sites = min(nm, curr_N * self.L)
					#t("Number of mutating sites: ",num_mut_sites)
					newMut += num_mut_sites
					if num_mut_sites > 0:
						#print("MUTATION!!!")
						#print(blub, " mut ", self.p_mut, " N ", curr_N, " L ", self.L, " tau ", tau)
						#for i in mut_sites:
						for _ in range(num_mut_sites):

							i = random.randrange(curr_N * self.L)
							seq_ind = math.floor(i / self.L)
							seq_pos = i % self.L

							#seq = next(x for x, val in enumerate(np.cumsum(species_dict.values())) if x > seq_ind)
							cum_N = 0
							# choose sequence to mutate (over cumsum)
							for seq in species_dict.keys():
								cum_N += species_dict[seq]
								if cum_N > seq_ind:
									break

							s_mut = list(seq)
							# seq[seq_pos] = self._mutateBase(new_set[seq_ind][seq_pos],initSeq[seq_pos])
							s_mut[seq_pos] = self.mutateBase(s_mut[seq_pos])
							mut_seq = "".join(s_mut)

							# create or increment new species
							species_dict[mut_seq] = species_dict.get(mut_seq, 0) + 1

							# decrement mutated sequence, or delete if it was the last one
							if species_dict[seq] > 1:
								species_dict[seq] = species_dict[seq] - 1
							else:
								del species_dict[seq]
								delSpecies += 1
			else:
				t = curr_t


		return species_dict_per_t, newMut_per_t, delSpecies_per_t,\
			{"meanFitness":meanFitness_per_t, "minFitness":minFitness_per_t, "maxFitness":maxFitness_per_t}

	def _get_repl(self, t):
		curr_p_repl = self.p_repl
		#determine replication rate at time t
		if self.p_repl is None:
			# TODO check, which one was publication??
			curr_p_repl = math.sin(t * 0.11) / 15 + 1.03 #original
			curr_p_repl = math.sin(t * 0.11) / 25 + 1.01
			curr_p_repl = math.sin(t * 0.11) / 15 + 1.02
			#curr_p_repl = math.sin(t * 0.1) / 20 + 1.05
			#curr_p_repl = math.sin(t * 0.11) / 15 + 1.05
		if(self.p_repl2 is not None):
			# Switch the replication rate somewhere in the middle if that is the correct mode (p_repl2!=0)
			if (t >= self.t_switch) and (self.p_repl2 > 0):
				curr_p_repl = self.p_repl2
			#if  repl2 is negative, switch only in one timestep and switch back to repl
			elif(self.p_repl2 < 0):
				if(t == self.t_switch):
					curr_p_repl=abs(self.p_repl2)
				else:
					curr_p_repl=self.p_repl
		return curr_p_repl


	@staticmethod
	def mutateBase(base):
		'''
		Mutate individual nucleotide according
		to Kimura model, fixed transition probability
		'''

		stay = 0
		alpha = 0.2
		beta = 0.2
		gamma = 0.6

		transitions = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
		transversions_alpha = {'A': 'C', 'C': 'A', 'T': 'G', 'G': 'T'}
		transversions_beta = {'A': 'G', 'G': 'A', 'T': 'C', 'C': 'T'}

		r = np.random.uniform(0, 1)

		if stay <= r < (stay + gamma):
			return transitions[base]
		elif (stay + gamma) <= r < (stay + gamma + alpha):
			return transversions_alpha[base]
		elif (stay + gamma + alpha) <= r < (stay + gamma + alpha + beta):
			return transversions_beta[base]
		return base
