#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	construct_nat_pair(args.input,args.output,args.complementary_length)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def check_exon_overlap_between_mlt(mlt1_exon_number,mlt1_exon_start,mlt1_exon_end,mlt2_exon_number,mlt2_exon_start,mlt2_exon_end):
	overlap_flag = "no"
	mlt1_exon_start_list = [int(i) for i in mlt1_exon_start.split(",")[:-1]]
	mlt1_exon_end_list = [int(i) for i in mlt1_exon_end.split(",")[:-1]]
	mlt2_exon_start_list = [int(i) for i in mlt2_exon_start.split(",")[:-1]]
	mlt2_exon_end_list = [int(i) for i in mlt2_exon_end.split(",")[:-1]]
	for i in range(0,mlt1_exon_number):
		for j in range(0,mlt2_exon_number):
			if mlt1_exon_end_list[i] >= mlt2_exon_start_list[j] and mlt2_exon_end_list[j] >= mlt1_exon_start_list[i]:
				overlap_flag = "yes"
				break
	return overlap_flag

def construct_nat_pair(input_gpd,output_gpd,complementary_length):
	dic_iso_info = {}
	chr_list = []
	dic_chr_strand_info = {}
	for line in input_gpd: # sort -k3,3 -k4,4 -k5,5n -k6,6n
		gene_id,isoform_id,chr,strand,tss,tts,full_length_LR_count,LR_count,exon_number,exon_start,exon_end = line.strip().split("\t")[:11]
		dic_iso_info[isoform_id] = line.strip()
		if chr not in dic_chr_strand_info.keys():
			dic_chr_strand_info[chr] = {}
			chr_list.append(chr)
			dic_chr_strand_info[chr]["+"] = {}
			dic_chr_strand_info[chr]["-"] = {}
			dic_chr_strand_info[chr]["+"]["isoform_list"] = []
			dic_chr_strand_info[chr]["-"]["isoform_list"] = []
			dic_chr_strand_info[chr][strand][isoform_id] = [int(tss),int(tts),int(exon_number),exon_start,exon_end]
			dic_chr_strand_info[chr][strand]["isoform_list"].append(isoform_id)
		else:
			dic_chr_strand_info[chr][strand][isoform_id] = [int(tss),int(tts),int(exon_number),exon_start,exon_end]
			dic_chr_strand_info[chr][strand]["isoform_list"].append(isoform_id)
#			print dic_chr_strand_info[chr][strand]["isoform_list"]

	dic_iso_nat_marker = {}
	for chr in chr_list:
		if dic_chr_strand_info[chr]["+"]["isoform_list"] == [] or dic_chr_strand_info[chr]["-"]["isoform_list"] == []: continue # no sense-antisense pair for this chromosome
		for iso_plus in dic_chr_strand_info[chr]["+"]["isoform_list"]:
			p_tss,p_tts,p_en,p_es,p_ee = dic_chr_strand_info[chr]["+"][iso_plus]
			for iso_minus in dic_chr_strand_info[chr]["-"]["isoform_list"]:
				m_tss,m_tts,m_en,m_es,m_ee = dic_chr_strand_info[chr]["-"][iso_minus]
				if p_tts <= m_tss: 
					break # no overlap
				elif p_tss >= m_tts:
					continue
				else:
					if (min([m_tts,p_tts]) - max([m_tss,p_tss])) < complementary_length: 
#						print iso_plus+"\t"+iso_minus
						continue # check the complementary length
					overlap_flag = check_exon_overlap_between_mlt(p_en,p_es,p_ee,m_en,m_es,m_ee)
					if p_tss > m_tss and p_tts > m_tts:
						if overlap_flag == "yes":
							if iso_plus not in dic_iso_nat_marker.keys():
								dic_iso_nat_marker[iso_plus] = {}
								dic_iso_nat_marker[iso_plus][iso_minus] = "HTH"
							else:
								dic_iso_nat_marker[iso_plus][iso_minus] = "HTH"
						else:
							if iso_plus not in dic_iso_nat_marker.keys():
								dic_iso_nat_marker[iso_plus] = {}
								dic_iso_nat_marker[iso_plus][iso_minus] = "INT"
							else:
								dic_iso_nat_marker[iso_plus][iso_minus] = "INT"
					elif m_tss > p_tss and m_tts > p_tts:
						if overlap_flag == "yes":
							if iso_plus not in dic_iso_nat_marker.keys():
								dic_iso_nat_marker[iso_plus] = {}
								dic_iso_nat_marker[iso_plus][iso_minus] = "TTT"
							else:
								dic_iso_nat_marker[iso_plus][iso_minus] = "TTT"
						else:
							if iso_plus not in dic_iso_nat_marker.keys():
								dic_iso_nat_marker[iso_plus] = {}
								dic_iso_nat_marker[iso_plus][iso_minus] = "INT"
							else:
								dic_iso_nat_marker[iso_plus][iso_minus] = "INT"
					else:
						if iso_plus not in dic_iso_nat_marker.keys():
							dic_iso_nat_marker[iso_plus] = {}
							dic_iso_nat_marker[iso_plus][iso_minus] = "EMB"
						else:
							dic_iso_nat_marker[iso_plus][iso_minus] = "EMB"
		for iso_minus in dic_chr_strand_info[chr]["-"]["isoform_list"]:
			m_tss,m_tts,m_en,m_es,m_ee = dic_chr_strand_info[chr]["-"][iso_minus]
			for iso_plus in dic_chr_strand_info[chr]["+"]["isoform_list"]:
				p_tss,p_tts,p_en,p_es,p_ee = dic_chr_strand_info[chr]["+"][iso_plus]
				if m_tts <= p_tss:
					break # no overlap
				elif m_tss >= p_tts:
					continue
				else:
					if (min([m_tts,p_tts]) - max([m_tss,p_tss])) < complementary_length:
#						print iso_plus+"\t"+iso_minus
						continue # check the complementary length
					overlap_flag = check_exon_overlap_between_mlt(m_en,m_es,m_ee,p_en,p_es,p_ee)
					if m_tss > p_tss and m_tts > p_tts:
						if overlap_flag == "yes":
							if iso_minus not in dic_iso_nat_marker.keys():
								dic_iso_nat_marker[iso_minus] = {}
								dic_iso_nat_marker[iso_minus][iso_plus] = "TTT"
							else:
								dic_iso_nat_marker[iso_minus][iso_plus] = "TTT"
						else:
							if iso_minus not in dic_iso_nat_marker.keys():
								dic_iso_nat_marker[iso_minus] = {}
								dic_iso_nat_marker[iso_minus][iso_plus] = "INT"
							else:
								dic_iso_nat_marker[iso_minus][iso_plus] = "INT"
					elif p_tss > m_tss and p_tts > m_tts:
						if overlap_flag == "yes":
							if iso_minus not in dic_iso_nat_marker.keys():
								dic_iso_nat_marker[iso_minus] = {}
								dic_iso_nat_marker[iso_minus][iso_plus] = "HTH"
							else:
								dic_iso_nat_marker[iso_minus][iso_plus] = "HTH"
						else:
							if iso_minus not in dic_iso_nat_marker.keys():
								dic_iso_nat_marker[iso_minus] = {}
								dic_iso_nat_marker[iso_minus][iso_plus] = "INT"
							else:
								dic_iso_nat_marker[iso_minus][iso_plus] = "INT"
					else:
						if iso_minus not in dic_iso_nat_marker.keys():
							dic_iso_nat_marker[iso_minus] = {}
							dic_iso_nat_marker[iso_minus][iso_plus] = "EMB"
						else:
							dic_iso_nat_marker[iso_minus][iso_plus] = "EMB"

	for sense_iso in dic_iso_info.keys():
		if sense_iso in dic_iso_nat_marker.keys():
			dic_nat_type = {"HTH":[],"TTT":[],"EMB":[],"INT":[]}
			for antisense_iso in dic_iso_nat_marker[sense_iso].keys():
				dic_nat_type[dic_iso_nat_marker[sense_iso][antisense_iso]].append(antisense_iso)
			for nat_type in dic_nat_type.keys():
				if dic_nat_type[nat_type] == []:
					dic_nat_type[nat_type] = ["*"]
			print >>output_gpd, "\t".join([dic_iso_info[sense_iso],",".join(dic_nat_type["HTH"]),",".join(dic_nat_type["TTT"]),",".join(dic_nat_type["EMB"]),",".join(dic_nat_type["INT"])])
		else:
			print >>output_gpd, "\t".join([dic_iso_info[sense_iso],"*","*","*","*"])

	input_gpd.close()
	output_gpd.close()

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS (+)
6. TTS (+)
7. number of support full-length long reads
8. number of support total long reads
9. exon count
10. exon start set
11. exon end set
12. For novel isoform, derived genic locus
13. For novel isoform, overlap percentage with derived genic locus
14. For novel singleton isoform, if it is located at the last exon of any known isoform. If yes, isoform ID otherwise '-'
15. For novel singleton isoform, the overlap percentage with the the last exon
16. For novel multi-exon isoform, number of splice sites are detected by anno and/or short reads; and the total number of splice sites
17. For novel multi-exon isoform, if the multi-exon isoform is the subset (based on splice junction combination) of known multi-exon isoform, isoform ID if yes otherwise '-'
18. For novel isoform, maximal length of polyA track in defined region
19. For novel isoform, maximal percentage of nucleotide A in defined region
20. For NAT analysis, NAT pair set with the type "Head to Head (HTH, 5'to5')"
21. For NAT analysis, NAT pair set with the type "Tail to Tail (TTT, 3'to3')"
22. For NAT analysis, NAT pair set with the type "Embedded (EMB, one transcript contained entirely within the other transcript)"
23. For NAT analysis, NAT pair set with the type "Intronic (INT, two transcripts overlap only within introns)"'''

	parser = argparse.ArgumentParser(description="Construct natural antisense transcript (NAT) pairs",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file, this file should be 'sort -k3,3 -k4,4 -k5,5n -k6,6n'")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd file")
	parser.add_argument('-l','--complementary_length',type=int,default=50,help="Minimal complementary length between sense and antisense transcript pairs (bp)")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
