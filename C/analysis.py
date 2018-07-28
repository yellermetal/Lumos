#! /usr/bin/python2.7

'''
/*
*	Solstice, Eclipse & Lumos Algorithms Analysis
*	Author: Ariel Livshits 2018
*/
'''

import numpy as np
import shelve
from subprocess import call
	
def reconfiguration_penalty_sweep(params):

	completion_times = {'Solstice':[], 'Eclipse':[], 'Lumos':[]}
	decomp_lengths = {'Solstice':[], 'Eclipse':[], 'Lumos':[]}
	running_times = {'Solstice':[], 'Eclipse':[], 'Lumos':[]}

	trials = params['trials']
	delta = params['delta']
	matrix_dim = params['matrix_dim']
	small_m  = params['small_m']
	large_m  = params['large_m']
	small_coeff  = params['small_coeff']
	large_coeff  = params['large_coeff']
    
	alg_num = 3

	for size in matrix_dim:
		call(["./decomp", str(trials), str(size), str(delta), str(small_m), str(large_m), str(small_coeff[0]), str(small_coeff[1]), str(large_coeff[0]), str(large_coeff[1])])
		log = open("./logfile.txt","r")
		for num in range(alg_num):
			line = log.readline().split(" ")
			completion_times[line[0]].append((line[2],line[5]))
			decomp_lengths[line[0]].append((line[1],line[4]))
			running_times[line[0]].append((line[3],line[6]))
		log.close()

	data = shelve.open('SMG.db')

	data['completion_times'] = completion_times
	data['decomp_lengths'] = decomp_lengths
	data['running_times'] = running_times
	data['params'] = params

	data.close()

def main():
	
	params = {}
	params['trials'] = 100
	params['delta'] = 25
	params['matrix_dim'] = [96, 160, 384]
	params['small_m'] = 48
	params['large_m'] = 16
	params['small_coeff'] = [1,16]
	params['large_coeff'] = [16,100]
   
	reconfiguration_penalty_sweep(params)


if __name__ == "__main__":
	main()