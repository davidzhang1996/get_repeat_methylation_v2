import sys
import subprocess
import time 
import os.path
#from __future__ import print_function

start_time=time.time()

#1
def main(): 
	#print ("Sorting Sam2Bed file: %s seconds") %(time.time()-start_time)

	#f=open("adrenal_gland_intersect.bed", "w")
	#subprocess.call (["sort", "-k1,1", "-k2,2n", "adrenal_gland.bed"], stdout=f)
	#f.close()

	#print ("Intersecting start/stop with Sam2Bed File: %s seconds") %(time.time()-start_time)
	#f=open ("adrenal_gland_intersect_TAR1-1.bed", "w")
	#subprocess.call(["../../../bin/bedtools/bedtools2-master/bin/bedtools", "intersect", "-a", "adrenal_gland_sorted.bed", "-b", "rmsk_start_stop_L1MC5a_3end-1.bed", "-wo"], stdout=f)
	#f.close()

	print ("Inputting repeat copy reference cpgs: %s seconds") %(time.time()-start_time)
	repeat_copy_cpg_ref=get_copy_cpg_ref() #1a

	print ("Inputting repeat copy consensus cpgs: %s seconds") %(time.time()-start_time)
	repeat_copy_cpg_cons=get_copy_cpg_cons() #1b

	print ("Inputting consensus sequence c: %s seconds") %(time.time()-start_time)
	repeat_consensus_c=get_consensus_c() #1c


	print ("Processing intersected reads: %s seconds") %(time.time()-start_time)
	repeat_consensus_c=get_intersected_reads(repeat_consensus_c, repeat_copy_cpg_ref, repeat_copy_cpg_cons) #1d


	#testing(repeat_consensus_c) #1e
	#for item in repeat_consensus_c:
	#	print item 
	#	print repeat_consensus_c[item]
	print ("Writing output to Wig file: %s seconds") %(time.time()-start_time)
	format_output(repeat_consensus_c) #1f


#1a
def get_copy_cpg_ref(): 
	cpg_ref_dict={}
	with open ("rmsk_c_ref_chrm22.bed", "r") as rmsk_cpg_ref: 
		name, forward_list, reverse_list= "" , [], []

		for line in rmsk_cpg_ref:
			line_list=line.split()
			cpg_ref_dict, name, forward_list, reverse_list = get_cpg(line_list, cpg_ref_dict, name, forward_list, reverse_list) #1ai
		cpg_ref_dict[name]= [forward_list, reverse_list]

	return cpg_ref_dict

#1b
def get_copy_cpg_cons(): 
	cpg_cons_dict={}
	with open ("rmsk_c_cons_chrm22.bed", "r") as rmsk_cpg_cons: 
		name, forward_list, reverse_list= "", [], []

		for line in rmsk_cpg_cons: 
			line_list=line.split()
			cpg_cons_dict, name, forward_list, reverse_list = get_cpg (line_list, cpg_cons_dict, name, forward_list, reverse_list) #1bi

		cpg_cons_dict[name]= [forward_list, reverse_list]


	return cpg_cons_dict

#1c
def get_consensus_c(): 
	cons_c_dict={}
	with open ("repeat_consensus_seqs.bed", "r") as rmsk_cpg_ref: 
		name, forward_list, reverse_list, forward_methyl, reverse_methyl, forward_total, reverse_total= "" , [], [], [], [] ,[] ,[]

		for line in rmsk_cpg_ref:
			line_list=line.split()
			cons_c_dict, name, forward_list, reverse_list = get_cpg(line_list, cons_c_dict, name, forward_list, reverse_list) #1ci


		for item in cons_c_dict:
			forward_list=cons_c_dict[item][0]
			reverse_list=cons_c_dict[item][1]
			for number in xrange(len(cons_c_dict[item][0])):
				forward_methyl.append(0)
				forward_total.append(0)

			for number in xrange(len(cons_c_dict[item][1])): 
				reverse_methyl.append(0)
				reverse_total.append(0)

			cons_c_dict[item]= [[forward_list, forward_methyl, forward_total], [reverse_list, reverse_methyl, reverse_total]]
			forward_methyl, forward_total, reverse_methyl, reverse_total= [], [], [], []

		#print "Cons_c_dict", cons_c_dict
	return cons_c_dict

#1d
def get_intersected_reads(repeat_consensus_c_dictionary, repeat_copy_cpg_ref_dictionary, repeat_copy_cpg_cons_dictionary): 
	repeat_consensus_c, repeat_copy_cpg_ref, repeat_copy_cpg_cons = repeat_consensus_c_dictionary, repeat_copy_cpg_ref_dictionary, repeat_copy_cpg_cons_dictionary 
	with open ("adrenal_gland_intersect_chrm22.bed", "r") as adrenal_gland_intersect:
		read_number=1
		for line in adrenal_gland_intersect: 
			print ("Read Added: %s seconds" %(time.time()-start_time))
			print "Read Number: ", read_number
			line_list=line.split() 
			subfamily, rmsk_strand, sam_strand=line_list[9], line_list[11], line_list[5]
			methyl_positions, total_positions= analyze_overlap(repeat_copy_cpg_ref, line_list, subfamily, rmsk_strand, sam_strand) #1di


			subfamily_consensus=subfamily[:subfamily.index(".")]
			if ((sam_strand=="+" and rmsk_strand=="+") or (sam_strand=="-" and rmsk_strand=="-")): 
				n=0
				consensus_positions_list, copy_positions_list=repeat_consensus_c[subfamily_consensus][n][0], repeat_copy_cpg_cons[subfamily][0]
			elif ((sam_strand=="+" and rmsk_strand=="-") or (sam_strand=="-" and rmsk_strand=="+")):
				n=1
				consensus_positions_list, copy_positions_list=repeat_consensus_c[subfamily_consensus][n][0], repeat_copy_cpg_cons[subfamily][1]

			consensus_index = find_start_position(consensus_positions_list, copy_positions_list[0], len(consensus_positions_list), 0) #1dii

			repeat_consensus_c=update_consensus_seqs(repeat_consensus_c, consensus_positions_list, copy_positions_list, methyl_positions, total_positions, subfamily_consensus, rmsk_strand, sam_strand, consensus_index, n) #1diii

			read_number+=1

		#print "Forward Methyl", repeat_consensus_c[subfamily_consensus][0][1]
		#print "Forward Total", repeat_consensus_c[subfamily_consensus][0][2]	

	return repeat_consensus_c

#1e
def testing(repeat_consensus_c_dict): 
	repeat_consensus_c=repeat_consensus_c_dict
	subfamily=repeat_consensus_c.keys()[0]

	list_position=repeat_consensus_c[subfamily][0][0].index(393)
	length_positions_forward=len(repeat_consensus_c[subfamily][0][0])
	length_methyl_forward=len(repeat_consensus_c[subfamily][0][1])

#1f
def format_output(repeat_consensus_c_dict): 
	output_file=open("output", "w")
	output_file.close()

	with open("output", "a") as output_file:
		repeat_consensus_c=repeat_consensus_c_dict 
		#print "Repeat_consensus_c", repeat_consensus_c
		for item in repeat_consensus_c: 
			output_file.write("VariableStep" + "\t"+"subfam="+item+"\n")
			for strand in xrange(2):
				for bp in xrange(len(repeat_consensus_c[item][strand][1])):
					if (repeat_consensus_c[item][strand][2][bp]==0):
						methyl_score=0
					else: 
						methyl_score=float(repeat_consensus_c[item][strand][1][bp])/repeat_consensus_c[item][strand][2][bp]

					if(strand==0):
						output_file.write(str(repeat_consensus_c[item][strand][0][bp])+"\t"+str(methyl_score)+"\t" +"+"+"\n")

					if (strand==1):
						output_file.write(str(repeat_consensus_c[item][strand][0][bp])+"\t"+str(methyl_score)+"\t" + "-"+"\n")

#1ai, 1bi, 1ci
def get_cpg(line_list_list, cpg_dict_dictionary, name_string, forward_list_list, reverse_list_list):
	line_list, cpg_dict, name, forward_list, reverse_list= line_list_list, cpg_dict_dictionary, name_string, forward_list_list, reverse_list_list 
	position, strand=line_list[1], line_list[5]

	if (line_list[3]!=name):
		name=line_list[3]
		forward_list, reverse_list=[],[]

	if (strand=="+"):
		forward_list.append(int(position))

	elif (strand=="-"):
		reverse_list.append(int(position))

	cpg_dict[name]=[forward_list,reverse_list]

	
	return cpg_dict, name, forward_list, reverse_list

#1di
#Input:	The Bed file that contains the locations of the c's in the reference strand of the repeat copy
#		The line list of the repeat copy positions intersected with the Sam file 
#		The subfamily of the repeat copy 

#Output:A methyl list which contains the number of methylated cytosines in the sam reads for each cytosine in the consensus sequence. 
#		A total list that contains the total number of cytosines and thymines in the sam reads for each cytosine in the consensus sequence.

#Algorithm Idea- 1: The  
def analyze_overlap(repeat_copy_cpg_ref_dictionary, line_list_list, subfamily_string, rmsk_strand, sam_strand): 
	repeat_copy_cpg_ref, line_list, subfamily, rmsk_strand, sam_strand = repeat_copy_cpg_ref_dictionary, line_list_list, subfamily_string, rmsk_strand, sam_strand

	rmsk_chrm_start, rmsk_chrm_stop = int(line_list[7]), int(line_list[8])
	sam_chrm_start, sam_chrm_stop, sam_read = int(line_list[1]), int(line_list[2]), line_list[3]
	length_sam_read=len(sam_read)

	forward_strand=repeat_copy_cpg_ref[subfamily][0]
	reverse_strand=repeat_copy_cpg_ref[subfamily][1]

	methyl_positions=[]
	total_positions=[]

	if (sam_strand=="+"):
		ref_position_list=forward_strand
		bp_unmethyl, bp_methyl= "C", "T"
	if (sam_strand=="-"):
		ref_position_list=reverse_strand
		bp_unmethyl, bp_methyl="G", "A"

	rmsk_reference_position=ref_position_list[0] #Starting position of the first c in the repeat copy relative to the chromosome
	rmsk_relative2reference_list=[] #Contains the positions of the c's relative to the starting position of the rmsk_reference_position 


	sam_rel_start= sam_chrm_start-rmsk_reference_position #Readjusts the starting position of the sam read so it's relative to the start of the repeat copy
	sam_rel_stop= length_sam_read+sam_rel_start-1 #Calculates the stopping position of the sam read relative to the start of the repeat copy

	for cytosine_position in xrange(len(ref_position_list)): 
		rmsk_relative2reference_list.append(ref_position_list[cytosine_position]-rmsk_reference_position) #Appends the starting positions of the c's relative to the reference position 

	for rmsk_rel_position in xrange(len(rmsk_relative2reference_list)): 
		sam_bp_position= rmsk_relative2reference_list[rmsk_rel_position]-sam_rel_start #The position along the sam read that corresbonds to the cytosines in the rmsk repeat copy
																					  #It can be negative or positive, depending on how the sam read aligns with the repeat copy 

		if (sam_bp_position<0 or sam_bp_position>length_sam_read-1): #Checks to make sure that the sam_bp_position refers to a legitimate bp on the sam read
			total_positions.append(0)
			methyl_positions.append(0)

		else: 
			if (sam_read[sam_bp_position]==bp_unmethyl or sam_read[sam_bp_position]==bp_methyl):
				total_positions.append(1)
				if (sam_read[sam_bp_position]==bp_methyl): 
					methyl_positions.append(1)
				else: 
					methyl_positions.append(0)
			else: 
				total_positions.append(0)
				methyl_positions.append(0)

	return methyl_positions, total_positions 

#1dii
def find_start_position(consensus_positions_list_list, chrm_pos_int, length_list_int, position_int): 
	consensus_positions_list, chrm_pos, length_list, position= consensus_positions_list_list, chrm_pos_int, length_list_int, position_int
	if (length_list>0):
		halfWay=length_list/2 

		midListChrmPos=consensus_positions_list[halfWay]
		leftList, rightList= consensus_positions_list[:halfWay], consensus_positions_list[halfWay+1:]
		lengthLeft, lengthRight= len(leftList), len(rightList)

		if (chrm_pos>midListChrmPos):
			position=halfWay+find_start_position(rightList, chrm_pos, lengthRight, position)+1
		elif (chrm_pos < midListChrmPos):
			position=find_start_position(leftList, chrm_pos, lengthLeft, position)
		elif (chrm_pos==midListChrmPos):
			return halfWay 

	return position 

#1diii
def update_consensus_seqs(repeat_consensus_c_dict, consensus_positions_list_list, copy_positions_list_list, methyl_positions_int, total_positions_int, subfamily_consensus_string, rmsk_strand_string, sam_strand_string, consensus_index_int, n_int): 

	repeat_consensus_c, consensus_positions_list, copy_positions_list, methyl_positions, total_positions, subfamily_consensus, rmsk_strand, sam_strand, consensus_index, n = repeat_consensus_c_dict, consensus_positions_list_list, copy_positions_list_list, methyl_positions_int, total_positions_int, subfamily_consensus_string, rmsk_strand_string, sam_strand_string, consensus_index_int, n_int

	copy_index, consensus_length, copy_length = 0, len(consensus_positions_list), len(copy_positions_list)

	line_count=1

	#print "test"
	while (copy_index<copy_length and consensus_index<consensus_length and copy_index>=0 and consensus_index>=0): 
		#print "test"
		while (consensus_index<consensus_length and copy_index<copy_length and copy_index>=0 and consensus_index>=0): 
			#print "test"
			#print "Consensus_positions_list[consensus_index]: ", consensus_positions_list[consensus_index]
			#print "copy_positions_list[copy_index]: ", copy_positions_list[copy_index]
			if (consensus_positions_list[consensus_index]==copy_positions_list[copy_index]): 
				#print "test"
				#print "methyl_positions[copy_index]: ", methyl_positions[copy_index]
				#print "total_positions[copy_index]: ", total_positions[copy_index]
				repeat_consensus_c[subfamily_consensus][n][1][consensus_index]+=methyl_positions[copy_index]
				repeat_consensus_c[subfamily_consensus][n][2][consensus_index]+=total_positions[copy_index]
				copy_index+=1
				if (rmsk_strand=="+"):
					consensus_index+=1
				elif (rmsk_strand=="-"):
					consensus_index-=1
				break
			
			else: 
				if (rmsk_strand=="+"):
					consensus_index+=1
				elif (rmsk_strand=="-"): 
					consensus_index-=1
	#print "Methyl: ", repeat_consensus_c[subfamily_consensus][n][1]
	#print "Total: ", repeat_consensus_c[subfamily_consensus][n][2]
	return repeat_consensus_c


main()





