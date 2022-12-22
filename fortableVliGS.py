#uses the TAb3 from the cascade tables to create a prioritization table
#deletions and duplications can be inputed together
#insertions, inversions, translocations, and complex rearrangements can be inputed together
#the table from tab3 must be pasted into a txt file, comma separated. The output is, as well, a tab-separated txt file
#To run:
#
#python3 fortableVliGS.py infile del/dup/"" repeats_BD segmentalDuplications_BD GeneHancer_Cluster_interactions Loops_BD TADs_BD position_effect_BD observedvsexpected_list old/new conserved_regions_list HPOterms(comma separated) DEL/DUP/INV_DGRC_ID ""/INS_DGRC_ID caselist caseID
#
#infile is the tab3 input file; del/dup/"" indicates the type of input, being "" for all SVs;
#repeats_BD is a list of repetitive regions from UCSC genome browser, segmentalDuplications_BD is a list of segmental duplications from UCSC genome browser 
#GeneHancer_Cluster_interactions is the list of GeneHancer cluster of interactions
#Loops_BD is the list of loops of the chosen cell line
#TAD_BD is a list of TADs of the chosen cell line
#pos_effect_BD is a list of genes associated to position effect
#observedvsexpected_list s a list of genes and respective oe according to Gnomad
#old/new indicates the state of the tab3 input file: if tab3 was already including the merged BD is "new", otherwise is "old"
#conserved_regions_list is a list of conserved regions of the chosen species comparison, bin and identity
#HPOterms are the HPO terms associated to the case in analysis. They must be inputed separed by commas with no spaces between them.
#
#example for SVs:
#python3.6 fortableVliGS.py test_Svs.txt "" ~/for_cascat_tables/repeats ~/for_cascat_tables/segdups ~/for_cascat_tables/GeneHancer_cluster_interactions ~/for_cascat_tables/all_loops ~/software/bp2omim/data/hg38_tads/H1-ESC_Dixon_2015-raw_TADs.txt ~/for_cascat_tables/pos_effect.txt ~/software/bp2omim/data/oe_alt old ~/HCNE/HCE_conserved_mouse_90pc_50bin HP:0000798,HP:0012207 INV_DGRC_BD INS_DGRC_BD case_list DGRC0005
#
############Dependencies:
from sys import argv
import marrvel_input
import treat_genes_for_table
import cross_deletions_with_dgv_results_Vfor_new
import getcl
import subprocess
import bitly_api
import classes_for_cpx
import urllib
import requests
import time

def read_table(infile):
	"""reads the trans/inv/ins/cpx input file
	creates a dictionary with the type of alteration and the
	cluster ID as key, and a list of all the complete lines, 
	as is, as elements of a list. return the dictionary"""
	f=open(infile)
	dic={}
	a=1
	temp=""
	for i in f:
		if len(i)>5:
			line=i.split("\t")
			if temp=="":
				dic[line[0]+"_"+str(a)]=[i.strip("\n")]
				temp=line[1]+"_"+getch(line[4])+"_"+getch(line[5])
			elif temp==line[1]+"_"+getch(line[4])+"_"+getch(line[5]):
				dic[line[0]+"_"+str(a)].append(i.strip("\n"))
			else:
				a+=1
				dic[line[0]+"_"+str(a)]=[i.strip("\n")]
				temp=line[1]+"_"+getch(line[4])+"_"+getch(line[5])
	f.close()
	return dic

def read_dgrc_bd_cpx(is_ins, infile, caselist):
	f=open(infile)
	dd={}
	for i in f:
		line=i.split("\t")
		dd[0]=[line[1], getch(line[3]), int(line[4].replace(",","")), int(line[5].replace(",","")), line[10].split(",")]
	f.close()
	caselists=[]
	f=open(caselist)
	for i in f:
		caselists.append(i.strip())
	return dd, caselists


def do_it(dic, tt, reps, segdup, genehancer, loops, tadfile, poesf, oe, old_new, cncr, terms, DGRC_BD, DGRC_BD_ins, DGRC_BD_cpx, caselistt, caseid, gwas, gtex):
	"""directs the work. receives the dictionary with the elements of the input files created by
	read_table or read_for_del_dup, and a set of auxiliary files, and calls the functions for the
	respective SV/CNV type."""
	if tt!="":
		data_dic, caselist=read_dgrc_bd(False, DGRC_BD, caselistt)
	else:
		inv_data_dic, caselist=read_dgrc_bd(False, DGRC_BD, caselistt)
		cpx_data_dic, caselist=read_dgrc_bd_cpx(False, DGRC_BD_cpx, caselistt)
		ins_data_dic, caselist=read_dgrc_bd(True, DGRC_BD_ins, caselistt)
	for key,value in dic.items():
		if tt=="dup":
			read_del_dup(key, value, "/home/ficheiros_referencia/ACMG_59_Genes", "/home/dgrc/software/bp2omim/data/DDG2P_21_2_2019.csv", "/home/dgrc/software/bp2omim/data/clingen2", "/home/dgrc/software/bp2omim/data/panel_list", tt, reps, segdup, genehancer, loops, tadfile, poesf, oe, "/home/ficheiros_referencia/Gtext_highest_expressed_transcripts", "/home/ficheiros_referencia/databases/merged_database/merged_dup", "/home/ficheiros_referencia/databases/genes_syndrome_del_dup", "/home/ficheiros_referencia/clinical_exome", "/home/ficheiros_referencia/infertility", old_new, cncr, terms, data_dic, caseid, caselist, gwas, gtex)
		elif tt=="del":
			read_del_dup(key, value, "/home/ficheiros_referencia/ACMG_59_Genes", "/home/dgrc/software/bp2omim/data/DDG2P_21_2_2019.csv", "/home/dgrc/software/bp2omim/data/clingen2", "/home/dgrc/software/bp2omim/data/panel_list", tt, reps, segdup, genehancer, loops, tadfile, poesf, oe, "/home/ficheiros_referencia/Gtext_highest_expressed_transcripts", "/home/ficheiros_referencia/databases/merged_database/merged_del", "/home/ficheiros_referencia/databases/genes_syndrome_del_dup", "/home/ficheiros_referencia/clinical_exome", "/home/ficheiros_referencia/infertility", old_new, cncr, terms, data_dic, caseid, caselist, gwas, gtex)
		elif "ins" in key:
			read_clust_ins(value, "/home/ficheiros_referencia/ACMG_59_Genes", "/home/dgrc/software/bp2omim/data/DDG2P_21_2_2019.csv", "/home/dgrc/software/bp2omim/data/clingen2", "/home/dgrc/software/bp2omim/data/panel_list", reps, segdup, genehancer, loops, tadfile, poesf, oe, "/home/ficheiros_referencia/Gtext_highest_expressed_transcripts", "/home/ficheiros_referencia/databases/merged_database/merged_ins", "/home/ficheiros_referencia/databases/genes_syndrome_del_dup", "/home/ficheiros_referencia/clinical_exome","/home/ficheiros_referencia/infertility", old_new, cncr, terms, ins_data_dic, caseid, caselist, gwas, gtex)
		elif key.startswith("inv_"):
			read_clust_inv(value, "/home/ficheiros_referencia/ACMG_59_Genes", "/home/dgrc/software/bp2omim/data/DDG2P_21_2_2019.csv", "/home/dgrc/software/bp2omim/data/clingen2", "/home/dgrc/software/bp2omim/data/panel_list", reps, segdup, genehancer, loops, tadfile, poesf, oe, "/home/ficheiros_referencia/Gtext_highest_expressed_transcripts", "/home/ficheiros_referencia/databases/merged_database/merged_inv", "/home/ficheiros_referencia/databases/genes_syndrome_del_dup", "/home/ficheiros_referencia/clinical_exome","/home/ficheiros_referencia/infertility", old_new, cncr, terms, inv_data_dic, caseid, caselist, gwas, gtex)
		elif "trans" in key:
			read_clust_trans(value, "/home/ficheiros_referencia/ACMG_59_Genes", "/home/dgrc/software/bp2omim/data/DDG2P_21_2_2019.csv", "/home/dgrc/software/bp2omim/data/clingen2", "/home/dgrc/software/bp2omim/data/panel_list", reps, segdup, genehancer, loops, tadfile, poesf, oe, "/home/ficheiros_referencia/Gtext_highest_expressed_transcripts", "/home/ficheiros_referencia/databases/genes_syndrome_del_dup", "/home/ficheiros_referencia/clinical_exome","/home/ficheiros_referencia/infertility", old_new, cncr, terms, gwas, gtex)
		else:
			read_clust_cpx(key.split("_")[0], value, "/home/ficheiros_referencia/ACMG_59_Genes", "/home/dgrc/software/bp2omim/data/DDG2P_21_2_2019.csv", "/home/dgrc/software/bp2omim/data/clingen2", "/home/dgrc/software/bp2omim/data/panel_list", reps, segdup, genehancer, loops, tadfile, poesf, oe, "/home/ficheiros_referencia/Gtext_highest_expressed_transcripts", "/home/ficheiros_referencia/databases/merged_database/merged_cpx", "/home/ficheiros_referencia/databases/genes_syndrome_del_dup", "/home/ficheiros_referencia/clinical_exome","/home/ficheiros_referencia/infertility", old_new, cncr, terms, cpx_data_dic, caseid, caselist, gwas, gtex)

def read_bds(infile, is_cov):
	"""reads the merged BD for overlap.
	Returns a dictionary with the chromosome as key, and a list by entry of said chromosome.
	The entry is [start, end, ID, freq].
	As a fail safe, in the end, the function checks if every possible chromosome is in the BD.
	If it is not, adds the key and an empty value."""
	cdd=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
	f=open(infile)
	dic={}
	for i in f:
		line=i.split("\t")
		if is_cov==True:
			tot, freq=manage_freq(line[3], line[4])
		else:
			tot=line[3]
			freq=line[4]
		if line[0] not in dic:
			dic[line[0]]=[[float(line[1]), float(line[2]), tot, freq.strip()]]
		else:
			dic[line[0]].append([float(line[1]), float(line[2]), tot, freq.strip()])
	for el in cdd:#fail safe for absent chromosomes
		if el not in dic:
			dic[el]=[[]]
	f.close()
	return dic

def make_overlap(reg1, bd, ispath):
	"""makes the genomic overlap between the SV/CNV region and the BD entry.
	For it to be considered overlap, both regions must overlap the other by at
	least 70%; if pathogenic, if the BD CNV is contained on our CNV, it also retrives as
	overlap. If there is overlap, the ID of the BD SV/CNV and respective frequency is retrieved as a list
	otherwise na, na is retrived."""
	sz1=reg1[1]-reg1[0]
	result=["NF","NF"]#returned if there is no overlap
	if sz1>0:
		if len(bd[0])!=0:
			for regbd in bd:
				szbd=regbd[1]-regbd[0]
				if szbd>0:
					ovl=get_Overlap(reg1, [regbd[0], regbd[1]])
					if (ovl/sz1>=0.7 and ovl/szbd>=0.7) or (ovl/szbd==1 and ispath==True) or ((ovl/sz1==1 or ovl/szbd==1) and ispath=="ins"):
						result=[regbd[2], regbd[3]]#return if overlap; [BD CNV/SV ID, freq]
						break
	return result

def make_cat_freq(freq, DGRC):
	"""according to the frequency of the alteration
	with overlap, returns the freq category to which
	the SV/CNV belongs"""
	if freq=="na" or freq=="NF":
		if DGRC!="NF" and DGRC!="na":
			return DGRC
		else:
			return "1 - novel"
	elif float(freq[:-2])<=0.01:
		if DGRC=="NF":
			return "2 - MAF<0.01%"
		else:
			return "3 - MAF<0.1%"
	elif float(freq[:-2])<=0.1:
		return "3 - MAF<0.1%"
	elif float(freq[:-2])<=0.5:
		return "4 - MAF<0.5%"
	elif float(freq[:-2])<1.0:
		return "5 - MAF<1%"
	elif float(freq[:-2])>=1.0:
		return "6 - MAF>1%"

def read_syndromes(infile):
	"""read the syndromes table. Retrieves a 
	dictionary with	dic[gene name]=syndromes"""
	f=open(infile)
	dic={}
	for i in f:
		line=i.split("\t")
		dic[line[0]]=line[-1].strip()
	f.close()
	return dic


def read_del_dup(key, v, acmgg, ddg2p, cling, panel, tt, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, merged_bd, synd, clin, inf, old_new, cncr, terms, data_dic, case, caselist, gwas, gtexx):#############################
	"""writes the output with all the information regarding dels/dups."""
	syndromes=read_syndromes(synd)#reads the syndromes BD
	clinical_exome=read_syndromes(clin)# reads the clinical exome BD
	infertility=read_syndromes(inf)# reads the infertility BD
	acmg=read_acmg(acmgg)#reads the actionable genes BD
	out=open(tt+"tab4", "a")#opens the output file
	temp=[]
	dd2p=get_dd2p(ddg2p)#reads the ddg2p database
	clg=clingen(cling)# reads the clingen database
	line=v.split("\t")
	if old_new=="old":#since we are in a transition step, old is for the tab3 before using the merged BDs
		dgrc=selecta(data_dic, line[13], case, caselist)
		gg="\t".join(line[20:46])
		merg=read_bds(merged_bd, False)
		path=get_freq("\t".join(line[36:46]))
		if "cov" in key:
			mm=make_overlap([int(line[7]), int(line[8])], merg[seps(line[6], 0)], False)
		else:
			mm=make_overlap([int(line[3]), int(line[4])], merg[seps(line[2], 0)], False)
		merged="\t".join(mm)
		cla=make_cat_freq(mm[1], dgrc)
	#else:
	#	path=get_freq("\t".join(line[32:42]))
	#	dgrc=selecta(data_dic, line[13], case, caselist)
	#	gg="\t".join(line[16:42])
	#	merged="\t".join(line[14:16])
	#	cla=make_cat_freq(line[15], dgrc)
	if "cov" in key:
		alt="Coverage"
		key="NF"
		#"\t"+line[9]+"\t"+merged+"\t"+path+"\tConfirmed?\t"+cla+"\t"+dgrc+"\t"+cl+"\t"+firstA+"\t"+dell+"\t"+firstB+"\t"+secondA+"\t"+secondB+"\t"+alt+"\t"+pq_nomenclature(getch(line[6]), "", "", "", str(int(float(line[7]))), str(int(float(line[8]))), "", "", "", "", "", "del", "")+"\t"+key+"\tNF\t"+gg.strip()+"\n")
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link("hg38", getch(seps(line[6], 0)), str(int(line[7])-1000), str(int(line[8])+1000), [line[7], line[8]])+'";"'+'[GRCh38] chr'+getch(seps(line[6], 0))+':g.'+'{:,}'.format(int(line[7]))+'-'+'{:,}'.format(int(line[8]))+'")'
		GenArch, cl, firstA, secondA, firstB, secondB, dell, high_ph, pval_ph =get_genes_del_dup(line, acmg, clg, dd2p, panel, getch(line[6]), getch(line[6]), line[7], line[8],reps, segdup, genehancer, loops, tadfile, poesf, tt, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, gwas, gtexx)
		out.write("seq[GRCh38] "+tt+"("+seps(line[6], 0)+")("+seps(line[6], 1)+")\t"+regt+"\t"+merged+"\t"+path+"\t"+dgrc+"\t"+line[9]+"\t"+GenArch+"\tUser revision\t"+cla+"\t"+cl+"\timpact rating\t"+high_ph+"\t"+pval_ph+"\t"+firstA+"\t"+dell+"\t"+firstB+"\t"+secondA+"\t"+secondB+"\t"+pq_nomenclature(getch(line[6]), "", "", "", str(int(float(line[7]))), str(int(float(line[8]))), "", "", "", "", "", "del", "")+"\t"+alt+"\tNA\tNA\n")
	elif len(line[6])>1:
		alt="Cluster/Coverage"
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link("hg38", getch(seps(line[2], 0)), str(int(line[3])-1000), str(int(line[4])+1000), [line[3], line[4]])+'";"'+'[GRCh38] chr'+getch(seps(line[2], 0))+':g.'+'{:,}'.format(int(line[3]))+'-'+'{:,}'.format(int(line[4]))+'")'
		GenArch, cl, firstA, secondA, firstB, secondB, dell, high_ph, pval_ph =get_genes_del_dup(line, acmg, clg, dd2p, panel, getch(line[2]), getch(line[2]), line[3], line[4],reps, segdup, genehancer, loops, tadfile, poesf, tt, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, gwas, gtexx)
#		out.write("seq[GRCh38] "+tt+"("+seps(line[2], 0)+")("+seps(line[2], 1)+")\tchr"+getch(seps(line[2], 0))+":g."+'{:,}'.format(int(float(line[3])))+"-"+'{:,}'.format(int(float(line[4])))+"\t"+line[5]+"\t"+merged+"\t"+path+"\tConfirmed?\t"+cla+"\t"+dgrc+"\t"+cl+"\t"+firstA+"\t"+dell+"\t"+firstB+"\t"+secondA+"\t"+secondB+"\t"+alt+"\t"+pq_nomenclature(getch(line[2]), "", "", "", str(int(float(line[3]))), str(int(float(line[4]))), "", "", "", "", "", "del", "")+"\t"+key+"\t"+line[1]+"\t"+gg.strip()+"\n")
		print("seq[GRCh38] "+tt+"("+seps(line[2], 0), seps(line[2], 1), regt, merged, path, dgrc, line[5], GenArch,"\tUser revision\t",cla,cl,"\timpact rating\t", high_ph, pval_ph, firstA, dell, firstB, secondA, secondB, pq_nomenclature(getch(line[2]), "", "", "", str(int(float(line[3]))), str(int(float(line[4]))), "", "", "", "", "", "del", ""),alt,key,line[1])
		out.write("seq[GRCh38] "+tt+"("+seps(line[2], 0)+")("+seps(line[2], 1)+")\t"+regt+"\t"+merged+"\t"+path+"\t"+dgrc+"\t"+line[5]+"\t"+GenArch+"\tUser revision\t"+cla+"\t"+cl+"\timpact rating\t"+high_ph+"\t"+pval_ph+"\t"+firstA+"\t"+dell+"\t"+firstB+"\t"+secondA+"\t"+secondB+"\t"+pq_nomenclature(getch(line[2]), "", "", "", str(int(float(line[3]))), str(int(float(line[4]))), "", "", "", "", "", "del", "")+"\t"+alt+"\t"+key+"\t"+line[1]+"\n")
	elif line[1]!="0":
		alt="Cluster"
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link("hg38", getch(seps(line[2], 0)), str(int(line[3])-1000), str(int(line[4])+1000), [line[3], line[4]])+'";"'+'[GRCh38] chr'+getch(seps(line[2], 0))+':g.'+'{:,}'.format(int(line[3]))+'-'+'{:,}'.format(int(line[4]))+'")'
		GenArch, cl, firstA, secondA, firstB, secondB, dell, high_ph, pval_ph=get_genes_del_dup(line, acmg, clg, dd2p, panel, getch(line[2]), getch(line[2]), line[3], line[4],reps, segdup, genehancer, loops, tadfile, poesf, tt, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, gwas, gtexx)
		#print(cla,"seq[GRCh38] ","("+seps(line[2], 0)+")("+seps(line[2], 1)+")","chr"+getch(seps(line[2], 0))+":g."+str(int(float(line[3])))+"-"+str(int(float(line[4]))),merged,path,line[5],"\tConfirmed?\t",cla,dgrc,cl,firstA,dell,firstB,secondA,secondB,alt,pq_nomenclature(getch(line[2]), "", "", "", str(int(float(line[3]))), str(int(float(line[4]))), "", "", "", "", "", "del", ""),key,line[1],gg.strip())
		#out.write("seq[GRCh38] "+tt+"("+seps(line[2], 0)+")("+seps(line[2], 1)+")\t"+getch(seps(line[2], 0))+":g."+'{:,}'.format(int(float(line[3])))+"-"+'{:,}'.format(int(float(line[4])))+"\t"+line[5]+"\t"+merged+"\t"+path+"\tConfirmed?\t"+cla+"\t"+dgrc+"\t"+cl+"\t"+firstA+"\t"+dell+"\t"+firstB+"\t"+secondA+"\t"+secondB+"\t"+alt+"\t"+pq_nomenclature(getch(line[2]), "", "", "", str(int(float(line[3]))), str(int(float(line[4]))), "", "", "", "", "", "del", "")+"\t"+key+"\t"+line[1]+"\t"+gg.strip()+"\n")
		out.write("seq[GRCh38] "+tt+"("+seps(line[2], 0)+")("+seps(line[2], 1)+")\t"+regt+"\t"+merged+"\t"+path+"\t"+dgrc+"\t"+line[5]+"\t"+GenArch+"\tUser revision\t"+cla+"\t"+cl+"\timpact rating\t"+high_ph+"\t"+pval_ph+"\t"+firstA+"\t"+dell+"\t"+firstB+"\t"+secondA+"\t"+secondB+"\t"+pq_nomenclature(getch(line[2]), "", "", "", str(int(float(line[3]))), str(int(float(line[4]))), "", "", "", "", "", "del", "")+"\t"+alt+"\t"+key+"\t"+line[1]+"\n")
	out.close()

"""
def make_tiny_link(full_link):
	#BITLY_ACCESS_TOKEN ="f2a30283ebe1f5e7b9585b6b4b7a30fa9478fd04"
	#BITLY_ACCESS_TOKEN="c94473e8e7c5cda7431e44db490ce02b44c20daf"
	cutly_token="e60db944e869898b5731c1078e53e22016efb"
	#para bitlu
	#access = bitly_api.Connection(access_token = BITLY_ACCESS_TOKEN)
	#short_url = access.shorten(full_link)
	#return short_url['url']
	#para cutly
	r = requests.get('http://cutt.ly/api/api.php?key={}&short={}'.format(cutly_token, full_link))
	return(r.text)
"""
"""
def make_tiny_link(full_link):
	#BITLY_ACCESS_TOKEN ="f2a30283ebe1f5e7b9585b6b4b7a30fa9478fd04"
	#BITLY_ACCESS_TOKEN="c94473e8e7c5cda7431e44db490ce02b44c20daf"
	#para bitlu
	#access = bitly_api.Connection(access_token = BITLY_ACCESS_TOKEN)
	#short_url = access.shorten(full_link)
	#return short_url['url']
	#para cutly
	cutly_token="e60db944e869898b5731c1078e53e22016efb"
	url = urllib.parse.quote(full_link)
	r = requests.get('http://cutt.ly/api/api.php?key={}&short={}'.format(cutly_token, url))
	time.sleep(10)
	aa=r.json()
	return(aa["url"]["shortLink"])

"""
def make_tiny_link(full_link):
	return "mocklink"
	
def make_ucsc_link(genome_version, chrr, start, end, bp):
	pars=["ct_HescTADs_7758","ct_ChromatinLoopsinHESc_8728"]
	if isinstance(bp, list)==False:
		bp=[bp, int(float(bp))+5]
	if genome_version=="hg38":#com hi-c de hesc
		if pars[1]=="":
			return make_tiny_link("https://genome.ucsc.edu/s/dgrc/SVInterpreter?hgTracks?db=hg38&position=chr"+chrr+"%3A"+str(start)+"%2D"+str(end)+"&hgsid=1170104349_awgJTWcBgaaTBSf49ESJR0iv7pyr&hideTracks=1&h1hescInsitu=full&knownGene=pack&omimGene2=pack&unipFullSeq=pack&gtexGeneV8=pack&geneHancerClusteredInteractionsDoubleElite=pack&clinGenHaplo=pack&clinGenTriplo=pack&clinGenGeneDisease=pack&gwasCatalog=dense&ignoreCookie=1&"+pars[0]+"=full&highlight=hg38.chr"+chrr+":"+str(bp[0])+"-"+str(int(bp[1]))+"#00FF00")
		else:
			return make_tiny_link("https://genome.ucsc.edu/s/dgrc/SVInterpreter?hgTracks?db=hg38&position=chr"+chrr+"%3A"+str(start)+"%2D"+str(end)+"&hgsid=1170104349_awgJTWcBgaaTBSf49ESJR0iv7pyr&hideTracks=1&h1hescInsitu=full&knownGene=pack&omimGene2=pack&unipFullSeq=pack&gtexGeneV8=pack&geneHancerClusteredInteractionsDoubleElite=pack&clinGenHaplo=pack&clinGenTriplo=pack&clinGenGeneDisease=pack&gwasCatalog=dense&ignoreCookie=1&"+pars[0]+"=full&"+pars[1]+"=full&highlight=hg38.chr"+chrr+":"+str(bp[0])+"-"+str(int(bp[1]))+"#00FF00")
	else:#com hic de lcl
		if pars[1]=="":
			return make_tiny_link("https://genome.ucsc.edu/s/dgrc/hg19_SVInterpreter?hgTracks?db=hg19&position=chr"+chrr+"%3A"+str(start)+"%2D"+str(end)+"&hgsid=1166713929_q7NVByZy7wTT2a8BiJMTLvelqXSl&hideTracks=1&gm12878Insitu=full&knownGene=pack&omimGene2=pack&unipFullSeq=pack&gtexGeneV8=pack&geneHancerClusteredInteractionsDoubleElite=pack&clinGenHaplo=pack&clinGenTriplo=pack&clinGenGeneDisease=pack&gwasCatalog=dense&ignoreCookie=1&"+pars[0]+"=full&highlight=hg19."+chrr+":"+str(bp[0])+"-"+str(int(bp[1]))+"#00FF00")
		else:
			return make_tiny_link("https://genome.ucsc.edu/s/dgrc/hg19_SVInterpreter?hgTracks?db=hg19&position=chr"+chrr+"%3A"+str(start)+"%2D"+str(end)+"&hgsid=1166713929_q7NVByZy7wTT2a8BiJMTLvelqXSl&hideTracks=1&gm12878Insitu=full&knownGene=pack&omimGene2=pack&unipFullSeq=pack&gtexGeneV8=pack&geneHancerClusteredInteractionsDoubleElite=pack&clinGenHaplo=pack&clinGenTriplo=pack&clinGenGeneDisease=pack&gwasCatalog=dense&ignoreCookie=1&"+pars[0]+"=full&"+pars[1]+"=full&highlight=hg19."+chrr+":"+str(bp[0])+"-"+str(int(bp[1]))+"#00FF00")

####################################aqui!
	
def dividebystrand(v,st):
	a=0
	plus=[]
	minus=[]
	while a<len(v):
		if st[a]=="+":
			plus.append(v[a])
		else:
			minus.append(v[a])
		a+=1
	if len(plus)==1:
		plus.append(plus[0]+1)
	if len(minus)==1:
		minus.append(minus[0]+1)
	return sorted(plus), sorted(minus)

def get_cpx_regs(ttype,v):
	regc="NA"
	pqC="NA"
	temp=[]
	for i in v:
		line=i.split("\t")
		if temp==[]:
			temp.append(line[1])#ID cluster 0
			temp.append(line[4])#chrA 1
			temp.append(line[5])#chrB 2
			temp.append([int(line[6])])#posA 3
			temp.append([int(line[7])])#posB 4
			temp.append([line[12]])#strandA 5
			temp.append([line[13]])#strandB 6
		else: 
			temp[0]==line[1]
			temp[3].append(int(line[6]))
			temp[4].append(int(line[7]))
			temp[5].append(line[12])
			temp[6].append(line[13])
	plusA, minusA=dividebystrand(temp[3], temp[5])
	plusB, minusB=dividebystrand(temp[4], temp[6])
	name="seq[GRCh38] "+ttype+"("+seps(temp[1], 0)+")("+seps(temp[1], 1)+";"+seps(temp[2], 1)+")"
	if ttype=="dupinvdup":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(max(plusB)+1025), [min(minusA), max(plusB)+25])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(max(plusB)+25)+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(max(plusA)+1025), [min(minusA), max(plusA)+25])+'";"'+'[GRCh38] dup_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(max(plusA)+25)+'")'
		pqA=str(min(minusA))+"-"+str(max(plusA)+25)
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(max(plusB)+1025), [min(minusA), max(plusB)+25])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(max(plusB)+25)+'")'
		pqB=str(min(minusA))+"-"+str(max(plusB)+25)
		regc='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusB)-1000), str(max(plusB)+1025), [min(minusB), max(plusB)+25])+'";"'+'[GRCh38] dup_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusB))+'-'+'{:,}'.format(max(plusB)+25)+'")'
		pqC=str(min(minusB))+"-"+str(max(plusB)+25)
		sz=str(max(plusB)+25-min(minusA))
		rr=[min(minusA), max(plusB)+25]
	if ttype=="dupinvdel":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(min(minusB)+1000), [min(minusA), min(minusB)])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(min(minusB))+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(max(plusA)+1025), [min(minusA), max(plusA)+25])+'";"'+'[GRCh38] dup_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(max(plusA)+25)+'")'
		pqA=str(min(minusA))+"-"+str(max(plusA)+25)
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(max(plusB)+1025), [min(minusA), max(plusB)+25])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(max(plusB)+25)+'")'
		pqB=str(min(minusA))+"-"+str(max(plusB)+25)
		regc='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusB)-975), str(min(minusB)+1000), [max(plusB)+25, min(minusB)])+'";"'+'[GRCh38] del_chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusB)+25)+'-'+'{:,}'.format(min(minusB))+'")'
		pqC=str(max(plusB)+25)+"-"+str(min(minusB))
		sz=str(min(minusB)-min(minusA))
		rr=[min(minusA), min(minusB)]
	if ttype=="delinvdup":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(max(plusB)+1025), [max(plusA)+25, max(plusB)+25])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(max(plusB)+25)+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(min(minusA)+1000), [max(plusA)+25, min(minusA)])+'";"'+'[GRCh38] del_chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(min(minusA))+'")'
		pqA=str(max(plusA)+25)+"-"+str(min(minusA))
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(max(plusB)+1025), [min(minusA), max(plusB)+25])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(max(plusB)+25)+'")'
		pqB=str(min(minusA))+"-"+str(max(plusB)+25)
		regc='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusB)-1000), str(max(plusB)+1025), [min(minusB), max(plusB)+25])+'";"'+'[GRCh38] dup_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusB))+'-'+'{:,}'.format(max(plusB)+25)+'")'
		pqC=str(min(minusB))+"-"+str(max(plusB)+25)
		sz=str(max(plusB)+25-max(plusA)+25)
		rr=[max(plusA)+25, max(plusB)+25]
	if ttype=="delinvdel":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(min(minusB)+1000), [max(plusA)+25, min(minusB)])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(min(minusB))+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(min(minusA)+1000), [max(plusA)+25, min(minusA)])+'";"'+'[GRCh38] del_chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(min(minusA))+'")'
		pqA=str(max(plusA)+25)+"-"+str(min(minusA))
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA+plusB)-1000), str(max(plusB+minusA)+1000), [min(minusA+plusB), max(plusB+minusA)])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA+plusB))+'-'+'{:,}'.format(max(plusB+minusA))+'")'
		pqB=str(min(minusA))+"-"+str(max(plusB))
		regc='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusB)-975), str(min(minusB)+1000), [max(plusB)+25, min(minusB)])+'";"'+'[GRCh38] del_chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusB)+25)+'-'+'{:,}'.format(min(minusB))+'")'
		pqC=str(max(plusB))+"-"+str(min(minusB))
		sz=str(min(minusB)-max(plusA)+25)
		rr=[max(plusA)+25, min(minusB)]
	if ttype=="delinv":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(min(minusB)+1000), [max(plusA)+25, min(minusB)])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(min(minusB))+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(min(minusA)+1000), [max(plusA)+25, min(minusA)])+'";"'+'[GRCh38] del_chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(min(minusA))+'")'
		pqA=str(max(plusA)+25)+"-"+str(min(minusA))
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(min(minusB)+1000), [min(minusA), min(minusB)])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(min(minusB))+'")'
		pqB=str(min(minusA))+"-"+str(min(minusB))
		sz=str(min(minusB)-max(plusA)+25)
		rr=[max(plusA)+25, min(minusB)]
	if ttype=="invdel":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(min(minusB)+1000), [max(plusA)+25, min(minusB)])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(min(minusB))+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(max(plusB)+1025), [max(plusA)+25, max(plusB)+25])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(max(plusB)+25)+'")'
		pqA=str(max(plusA)+25)+"-"+str(max(plusB)+25)
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusB)-975), str(min(minusB)+1000), [max(plusB)+25, min(minusB)])+'";"'+'[GRCh38] del_chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusB)+25)+'-'+'{:,}'.format(min(minusB))+'")'
		pqB=str(min(minusB))+"-"+str(max(plusB)+25)
		sz=str(min(minusB)-max(plusA)+25)
		rr=[max(plusA)+25, min(minusB)]
	if ttype=="dupinv":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(min(minusB)+1000), [min(minusA), min(minusB)])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(min(minusB))+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(max(plusA)+1025), [min(minusA), max(plusA)+25])+'";"'+'[GRCh38] dup_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(max(plusA)+25)+'")'
		pqA=str(min(minusA))+"-"+str(max(plusA)+25)
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(min(minusB)+1000), [min(minusA), min(minusB)])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA))+'-'+'{:,}'.format(min(minusB))+'")'
		pqB=str(min(minusA))+"-"+str(min(minusB))
		sz=str(min(minusB)-min(minusA))
		rr=[min(minusA), min(minusB)]
	if ttype=="invdup":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(max(plusB)+1025), [max(plusA)+25, max(plusB)+25])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(max(plusB)+25)+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-975), str(max(plusB)+1025), [max(plusA)+25, max(plusB)+25])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA)+25)+'-'+'{:,}'.format(max(plusB)+25)+'")'
		pqA=str(max(plusA)+25)+"-"+str(max(plusB)+25)
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusB)-1000), str(max(plusB)+1025), [min(minusB), max(plusB)+25])+'";"'+'[GRCh38] dup_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusB))+'-'+'{:,}'.format(max(plusB)+25)+'")'
		pqB=str(min(minusB))+"-"+str(max(plusB)+25)
		sz=str(max(plusB)+25-max(plusA)+25)
		rr=[max(plusA)+25, max(plusB)+25]
	if ttype=="invtdup":
		regt='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-1000), str(max(plusB+plusA)+1025), [min(minusA), max(plusB+plusA)+25])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(minusA))+'-'+'{:,}'.format(max(plusB+plusA)+25)+'")'
		rega='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(minusA)-2000), str(min(minusA)+1000), [min(minusA)-1000, min(minusA)])+'";"'+'[GRCh38] inv_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(minusA)-1000)+'-'+'{:,}'.format(min(minusA))+'")'
		pqA=str(min(minusA)-1000)+"-"+str(min(minusA))
		regb='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA+plusB)-1000), str(min(minusB)+1025), [max(plusA+plusB), max(minusB)])+'";"'+'[GRCh38] dup_chr'+getch(temp[1])+':g.'+'{:,}'.format(min(plusA+plusB))+'-'+'{:,}'.format(max(minusB))+'")'
		pqB=str(max(plusA+plusB))+"-"+str(min(minusB))
		sz=str(min(minusB)-min(minusA))
		rr=[min(minusA), max(plusB+plusA)+25]
	#return getch(temp[1]), name, regt, rega, pqA, pqB, regb, regc, pqC, sz, rr, temp[0], str(len(temp[3]), temp[-2], temp[-1])
	return getch(temp[1]), name, regt, rega, pqA, pqB, regb, regc, pqC, sz, rr, temp[0], str(len(temp[3]))

# old_new, cncr, terms, data_dic, case, caselist):
def read_clust_cpx(ttype, v, acmgg, ddg2p, cling, panel, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, merged_bd, synd, clin, inf, old_new, cncr, terms, data_dic, case, caselist, gwas, gtexx):
	ch, name, regt, regA, pqA, pqB, regB, regC, pqC, sz, rr, cll, pr=get_cpx_regs(ttype,v)
	syndromes=read_syndromes(synd)
	clinical_exome=read_syndromes(clin)
	infertility=read_syndromes(inf)
	acmg=read_acmg(acmgg)
	dd2p=get_dd2p(ddg2p)
	clg=clingen(cling)
	out=open("tab4", "a")
	if old_new=="old":
		merg=read_bds(merged_bd, False)
		mm=make_overlap(rr, merg[ch], False)
	else:
		mm=line[28:30]
	dgrc=selecta_cpx(ttype, data_dic,rr,ch, case, caselist)#[DGRC]=[invdup, chr, start, end, [case1, case2...]]
	merged="\t".join(mm)
	cla=make_cat_freq(mm[1], dgrc)
	if regC=="NA":
		GenArch, cl, firstA, secondA, firstB, secondB, high_ph, pval_ph =get_genes_with_C2(acmg, clg, ddg2p, panel, ch, pqA, pqB, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, ttype, False, gwas, gtexx)
		#out.write("\t".join(["seq[GRCh38] "+ttype+"("+seps(temp[1], 0)+")("+seps(temp[1], 1)+";"+seps(temp[2], 1)+")", "chr"+getch(temp[1])+":"+pqtot,"[GRCh38] chr"+getch(temp[1])+":g."+pqAA+"\t[GRCh38] chr"+getch(temp[2])+":g."+pqBB, sz, merged, "NF","confirmed",cla, "NF", cl, firstA, firstB, secondA, secondB, "NA","NA",temp[0],str(len(temp[3])),"NF", "\n"]))
		out.write("\t".join([name, regt, regA, regB, regC, merged,dgrc, sz,"NA", GenArch, "",cla, cl,"",high_ph, pval_ph, firstA, firstB, "NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF", secondA, secondB,"NF","NF","NF","NF","NF","NF", "NA","NA",cll,pr, "\n"]))
	else:
		###para reg C tambem
		GenArch, cl, firstA, secondA, firstB, secondB, firstC, secondC, high_ph, pval_ph =get_genes_with_C3(acmg, clg, ddg2p, panel, ch, pqA, pqB, pqC, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, ttype, False, gwas, gtexx)
		out.write("\t".join([name, regt, regA, regB, regC, merged,dgrc, sz,"NA", GenArch, "",cla, cl,"", high_ph, pval_ph, firstA, firstB, firstC, secondA, secondB, secondC, "NA","NA",cll,pr, "\n"]))
	out.close()


def selecta_cpx(typee, data_dic,rr, ch, case, caselist):
	freq="NF"
	for key,value in data_dic.items():
		if value[0]==typee:
			if value[1]==ch:
				ovl=get_Overlap(rr, [value[2], value[3]]) 
				if (ovl/sz1>=0.7 and ovl/szbd>=0.7):
					if case in value[-1]:
						freq=key+" ("+str(len(value[-1])-1)+"/"+str(len(caselist)-1)+")"
					else:
						freq=key+" ("+str(len(value[-1]))+"/"+str(len(caselist))+")"
	return freq


def get_genes_with_C2(acmg, clingenn, ddg2p, pannel, chra, bp1, bp2, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex_file, synd, clinicalexome, infertility, cncr, terms, tt, is_ins, gwas, gtexx):
	dic_gtex=read_gtex(gtex_file)
	Conserved_reg=read_cncr(cncr)
	tadA_genes, tadA4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chra, int(bp1.split("-")[0]))
	tadB_genes, tadB4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chra, int(bp2.split("-")[0]))
	gA, biotypeA=treat_genes_for_table.getdisrp(chra, int(bp1.split("-")[0]), int(bp1.split("-")[1]))
	ddA=get_oe_for_tad(gA, oe)
	gB, biotypeB=treat_genes_for_table.getdisrp(chra, int(bp2.split("-")[0]), int(bp2.split("-")[1]))
	ddB=get_oe_for_tad(gB, oe)
	ta=get_oe_for_tad(tadA_genes, oe)
	tb=get_oe_for_tad(tadB_genes, oe)
	noncodA=[["NF","NF","NF"]]
	noncodB=[["NF","NF","NF"]]
	cpsA=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	cpsB=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	phenscoreA=""
	phenscoreB=""
	if ddA==["NF"]:
		firstA=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondA=repeats_segdups(reps, chra, bp1)+"\t"+repeats_segdups(segdup, chra, bp1)+"\t"+GeneHancer(tadA_genes, genehancer, chra, bp1)+"\t"+getloops(loops, chra, bp1)+"\t"+pos_effect(tadA4_genes, poesf)+"\t"+"& ".join(ta)
	if ddB==["NF"]:
		firstB=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondB=repeats_segdups(reps, chra, bp2)+"\t"+repeats_segdups(segdup, chra, bp2)+"\t"+GeneHancer(tadB_genes, genehancer, chra, bp2)+"\t"+getloops(loops, chra, bp2)+"\t"+pos_effect(tadB4_genes, poesf)+"\t"+"& ".join(tb)
	GenArch, secondA, secondB=edit_second(secondA, secondB, "")
	if ddA==["NF"] and ddB==["NF"]:
		cl, firstA, secondA, firstB, secondB=get_other_ops(chra, chra, bp1, bp2, "", Conserved_reg, firstA, secondA, firstB, secondB, "", "")
		return GenArch, cl, firstA, secondA, firstB, secondB, "NF", "NF"
	if ddA!=["NF"]:
		genesA=str(len(set(ddA)))
		j1A, j2A, j3A, j4A, j5A, j6A=[['NF'] for i in range(6)]
		for el in ddA:
			if el.split(" ")[0] in acmg:
				if "NF" in j1A: j1A.remove("NF")
				j1A.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2A: j2A.remove("NF")
				j2A.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3A: j3A.remove("NF")
				j3A.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4A: j4A.remove("NF")
				j4A.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5A: j5A.remove("NF")
				j5A.append(gwas[el.split(" ")[0]])#gwas
			if el.split(" ")[0] in gtexx:
				if "NF" in j6A: j6A.remove("NF")
				j6A.append(gtexx[el.split(" ")[0]])#gtex
		freqA, genecardA=get_oe(ddA)
		cpsA, noncodA, phenscoreA=get_phen_ins(ddA, clingenn, ddg2p, pannel, chra, str(bp1), "", dic_gtex, terms)
		firstA=[genesA,"& ".join(ddA),",".join(biotypeA), "",freqA, genecardA, cpsA[1],cpsA[2],phenscoreA,cpsA[3],cpsA[4],cpsA[5],cpsA[7],cpsA[0],cpsA[8], "& ".join(j5A), "& ".join(j6A), "& ".join(j3A), "& ".join(j4A),"& ".join(j2A),"& ".join(j1A),cpsA[6]]
	if ddB!=["NF"]:
		genesB=str(len(set(ddB)))
		j1B, j2B, j3B, j4B, j5B, j6B=[['NF'] for i in range(6)]
		for el in ddB:
			if el.split(" ")[0] in acmg:
				if "NF" in j1B: j1B.remove("NF")
				j1B.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2B: j2B.remove("NF")
				j2B.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3B: j3B.remove("NF")
				j3B.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4B: j4B.remove("NF")
				j4B.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5B: j5B.remove("NF")
				j5B.append(gwas[el.split(" ")[0]])#gwas
			if el.split(" ")[0] in gtexx:
				if "NF" in j6B: j6B.remove("NF")
				j6B.append(gtexx[el.split(" ")[0]])#gte			
		freqB, genecardB=get_oe(ddB)
		cpsB, noncodB, phenscoreB=get_phen_ins(ddB, clingenn, ddg2p, pannel, chra, str(bp2), "", dic_gtex, terms)
		firstB=[genesB,"& ".join(ddB), ",".join(biotypeB), "",freqB, genecardB, cpsB[1],cpsB[2],phenscoreB,cpsB[3],cpsB[4],cpsB[5],cpsB[7],cpsB[0],cpsB[8],"& ".join(j5B),"& ".join(j6B),  "& ".join(j3B), "& ".join(j4B), "& ".join(j2B), "& ".join(j1B), cpsB[6]]
	#print("slslslslsl", tt, ddA, classes_for_cpx.get_gene_class(ddA, cpsA, noncodA), classes_for_cpx.get_gene_class(ddB, cpsB, noncodB), chra, bp1.split("-"), bp2.split("-"))###calcula em que grau de prioritizacao é que aquilo fica
	cl, afectA, afectB= classes_for_cpx.interpret_to_class(tt, classes_for_cpx.get_gene_class(ddA, cpsA, noncodA), classes_for_cpx.get_gene_class(ddB, cpsB, noncodB), "", [chra, int(bp1.split("-")[0]), int(bp1.split("-")[1]), chra, int(bp2.split("-")[0]), int(bp2.split("-")[1])])###calcula em que grau de prioritizacao é que aquilo fica
	if firstA[0]!="0":
		#firstA[2]=afectA
		firstA[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chra+'%3A'+bp1+'";"'+afectA+'")'
	if firstB[0]!="0":
		#firstB[2]=afectB
		firstB[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chra+'%3A'+bp2+'";"'+afectB+'")'
	high_ph, pval_ph=get_highest_ph(phenscoreA, phenscoreB, "")
	if cl=="NF":
		cl, firstA, secondA, firstB, secondB=get_other_ops(chra, chra, bp1, bp2, "", Conserved_reg, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondA, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondB, "", "")
		return GenArch, cl, firstA, secondA, firstB, secondB, high_ph, pval_ph
	else:
		return GenArch, cl, "\t".join(firstA), secondA, "\t".join(firstB), secondB, high_ph, pval_ph

def get_highest_ph(phenscoreA, phenscoreB, phenscoreC):
	#print(phenscoreA, phenscoreB, phenscoreC)
	aa=phenscoreA.split("& ")+phenscoreB.split("& ")+phenscoreC.split("& ")
	#print("aa", aa)
	high="NF"
	pvv="NF"
	hh=[]
	pv=[]
	for el in aa:
		if "PhenSSc" in el:
			gg=el.split(" ")
			hh.append(float(gg[1]))
			pv.append(float(gg[4]))
	if len(hh)>0:
		high=max(hh)
		ind=hh.index(high)
		pvv=pv[ind]
	return str(high), str(pvv)


def get_genes_with_C3(acmg, clingenn, ddg2p, pannel, chra, bp1, bp2, bp3, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex_file, synd, clinicalexome, infertility, cncr, terms, tt, is_ins, gwas, gtexx):
	dic_gtex=read_gtex(gtex_file)
	Conserved_reg=read_cncr(cncr)
	tadA_genes, tadA4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chra, int(bp1.split("-")[0]))
	tadB_genes, tadB4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chra, int(bp2.split("-")[0]))
	tadC_genes, tadC4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chra, int(bp3.split("-")[0]))
	gA, biotypeA=treat_genes_for_table.getdisrp(chra, int(bp1.split("-")[0]), int(bp1.split("-")[1]))
	ddA=get_oe_for_tad(gA, oe)
	gB, biotypeB=treat_genes_for_table.getdisrp(chra, int(bp2.split("-")[0]), int(bp2.split("-")[1]))
	ddB=get_oe_for_tad(gB, oe)
	gC, biotypeC=treat_genes_for_table.getdisrp(chra, int(bp3.split("-")[0]), int(bp3.split("-")[1]))
	ddC=get_oe_for_tad(gC, oe)
	ta=get_oe_for_tad(tadA_genes, oe)
	tb=get_oe_for_tad(tadB_genes, oe)
	tc=get_oe_for_tad(tadC_genes, oe)
	noncodA=[["NF","NF","NF"]]
	noncodB=[["NF","NF","NF"]]
	noncodC=[["NF","NF","NF"]]
	cpsA=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	cpsB=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	cpsC=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	phenscoreA=""
	phenscoreB=""
	phenscoreC=""
	if ddA==["NF"]:
		firstA=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondA=repeats_segdups(reps, chra, bp1)+"\t"+repeats_segdups(segdup, chra, bp1)+"\t"+GeneHancer(tadA_genes, genehancer, chra, bp1)+"\t"+getloops(loops, chra, bp1)+"\t"+pos_effect(tadA4_genes, poesf)+"\t"+"& ".join(ta)
	if ddB==["NF"]:
		firstB=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondB=repeats_segdups(reps, chra, bp2)+"\t"+repeats_segdups(segdup, chra, bp2)+"\t"+GeneHancer(tadB_genes, genehancer, chra, bp2)+"\t"+getloops(loops, chra, bp2)+"\t"+pos_effect(tadB4_genes, poesf)+"\t"+"& ".join(tb)
	if ddC==["NF"]:
		firstC=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondC=repeats_segdups(reps, chra, bp3)+"\t"+repeats_segdups(segdup, chra, bp3)+"\t"+GeneHancer(tadC_genes, genehancer, chra, bp3)+"\t"+getloops(loops, chra, bp3)+"\t"+pos_effect(tadC4_genes, poesf)+"\t"+"& ".join(tc)
	GenArch, secondA, secondB, secondC =edit_second(secondA, secondB, secondC)
	if ddA==["NF"] and ddB==["NF"] and ddC==["NF"]:
		cl, firstA, secondA, firstB, secondB, firstC, secondC, =get_other_ops(chra, chra, bp1, bp2,bp3, Conserved_reg, firstA, secondA, firstB, secondB,  firstC, secondC)
		return GenArch, cl, firstA, secondA, firstB, secondB, firstC, secondC, "NF", "NF"
		#return GenArch, get_other_ops(chra, chra, bp1, bp2, bp3, Conserved_reg, firstA, secondA, firstB, secondB, firstC, secondC)
	if ddA!=["NF"]:
		genesA=str(len(set(ddA)))
		j1A, j2A, j3A, j4A, j5A, j6A=[['NF'] for i in range(6)]
		for el in ddA:
			if el.split(" ")[0] in acmg:
				if "NF" in j1A: j1A.remove("NF")
				j1A.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2A: j2A.remove("NF")
				j2A.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3A: j3A.remove("NF")
				j3A.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4A: j4A.remove("NF")
				j4A.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5A: j5A.remove("NF")
				j5A.append(gwas[el.split(" ")[0]])#gwas
			if el.split(" ")[0] in gtexx:
				if "NF" in j6A: j6A.remove("NF")
				j6A.append(gtexx[el.split(" ")[0]])#gte	
		freqA, genecardA=get_oe(ddA)
		cpsA, noncodA, phenscoreA=get_phen_ins(ddA, clingenn, ddg2p, pannel, chra, str(bp1), "", dic_gtex, terms)
		firstA=[genesA,"& ".join(ddA), ",".join(biotypeA), "",freqA, genecardA, cpsA[1],cpsA[2],phenscoreA,cpsA[3],cpsA[4],cpsA[5],cpsA[7],cpsA[0],cpsA[8],"& ".join(j5A),"& ".join(j6A),"& ".join(j3A),"& ".join(j4A),"& ".join(j2A),"& ".join(j1A),cpsA[6]]
	if ddB!=["NF"]:
		genesB=str(len(set(ddB)))
		j1B, j2B, j3B, j4B, j5B, j6B=[['NF'] for i in range(6)]
		for el in ddB:
			if el.split(" ")[0] in acmg:
				if "NF" in j1B: j1B.remove("NF")
				j1B.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2B: j2B.remove("NF")
				j2B.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3B: j3B.remove("NF")
				j3B.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4B: j4B.remove("NF")
				j4B.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5B: j5B.remove("NF")
				j5B.append(gwas[el.split(" ")[0]])#gwas
			if el.split(" ")[0] in gtexx:
				if "NF" in j6B: j6B.remove("NF")
				j6B.append(gtexx[el.split(" ")[0]])#gte	
		freqB, genecardB=get_oe(ddB)
		cpsB, noncodB, phenscoreB=get_phen_ins(ddB, clingenn, ddg2p, pannel, chra, str(bp2), "", dic_gtex, terms)
		firstB=[genesB,"& ".join(ddB),",".join(biotypeB), "",freqB, genecardB, cpsB[1],cpsB[2],phenscoreB,cpsB[3],cpsB[4],cpsB[5],cpsB[7],cpsB[0],cpsB[8],"& ".join(j5B),"& ".join(j6B),"& ".join(j3B),"& ".join(j4B),"& ".join(j2B),"& ".join(j1B),cpsB[6]]
	if ddC!=["NF"]:
		genesC=str(len(set(ddC)))
		j1C, j2C, j3C, j4C, j5C, j6C=[['NF'] for i in range(6)]
		for el in ddC:
			if el.split(" ")[0] in acmg:
				if "NF" in j1C: j1C.remove("NF")
				j1C.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2C: j2C.remove("NF")
				j2C.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3C: j3C.remove("NF")
				j3C.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4C: j4C.remove("NF")
				j4C.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5C: j5C.remove("NF")
				j5C.append(el.split("(")[0])#gwas
			if el.split(" ")[0] in gtexx:
				if "NF" in j6C: j6C.remove("NF")
				j6C.append(el.split("(")[0])#gte				
		freqC, genecardC=get_oe(ddC)
		cpsC, noncodC, phenscoreC=get_phen_ins(ddC, clingenn, ddg2p, pannel, chra, str(bp3), "", dic_gtex, terms)
		firstC=[genesC,"& ".join(ddC), ",".join(biotypeC), "",freqC, genecardC, cpsC[1],cpsC[2],phenscoreC,cpsC[3],cpsC[4],cpsC[5],cpsC[7],cpsC[0],cpsC[8],"& ".join(j5C),"& ".join(j6C),"& ".join(j3C),"& ".join(j4C),"& ".join(j2C),"& ".join(j1C),cpsC[6]]
		#print("get geneclass", "ddA", ddA, "cpsA", cpsA, "noncodA",noncodA, "ddB",ddB, "cpsB", cpsB, "noncodB",noncodB, "0", 0, "is_ins", is_ins)
	cl, afectA, afectB, afectC= classes_for_cpx.interpret_to_class(tt, classes_for_cpx.get_gene_class(ddA, cpsA, noncodA), classes_for_cpx.get_gene_class(ddB, cpsB, noncodB), classes_for_cpx.get_gene_class(ddC, cpsC, noncodC), [chra, int(bp1.split("-")[0]), int(bp1.split("-")[1]), int(bp2.split("-")[0]), int(bp2.split("-")[1]), int(bp3.split("-")[0]), int(bp3.split("-")[1])])###calcula em que grau de prioritizacao é que aquilo fica
	if firstA[0]!="0":
		#firstA[2]=afectA
		firstA[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chra+'%3A'+bp1+'";"'+afectA+'")'
	if firstB[0]!="0":
		#firstB[2]=afectB
		firstB[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chra+'%3A'+bp2+'";"'+afectB+'")'
	if firstC[0]!="0":
		#firstC[2]=afectC
		firstC[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chra+'%3A'+bp3+'";"'+afectC+'")'
	high_ph, pval_ph=get_highest_ph(phenscoreA, phenscoreB, phenscoreC)
	if cl=="NF":
		cl, firstA, secondA, firstB, secondB, firstC, secondC, =get_other_ops(chra, chra, bp1, bp2,bp3, "", Conserved_reg, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondA, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondB,  ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondC)
		return GenArch, cl, firstA, secondA, firstB, secondB, firstC, secondC, high_ph, pval_ph
		#return GenArch, get_other_ops(chra, chrb, bp1, bp2, bp3, Conserved_reg, firstA, secondA, firstB, secondB, firstC, secondC)
	else:
		return GenArch, cl, "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(firstC), secondC, high_ph, pval_ph
		
		
def read_clust_inv(v, acmgg, ddg2p, cling, panel, reps, segdup, genehancer, loops, tadfile, poesf,oe, gtex, merged_bd, synd, clin, inf, old_new, cncr, terms, data_dic, case, caselist, gwas, gtexx):
	syndromes=read_syndromes(synd)
	clinical_exome=read_syndromes(clin)
	infertility=read_syndromes(inf)
	acmg=read_acmg(acmgg)
	temp=[]
	dd2p=get_dd2p(ddg2p)
	clg=clingen(cling)
	out=open("tab4", "a")
	for i in v:
		line=i.split("\t")
		if temp==[]:
			temp.append(line[1])#clusterID 0
			temp.append(line[4])#chrAband 1
			temp.append(line[5])#chrBband 2
			temp.append([int(line[6])])#cordA 3
			temp.append([int(line[7])])#cordB 4
			temp.append([line[12]])#strandA 5
			temp.append([line[13]])#strandB 6
		else:
			temp[3].append(int(line[6]))#cordA 3
			temp[4].append(int(line[7]))#cordB 4
			temp[5].append(line[12])#strandA 5
			temp[6].append(line[13])#strandB 6
	plusA, minusA=dividebystrand(temp[3], temp[5])
	plusB, minusB=dividebystrand(temp[4], temp[6])
	szA=str(min(minusA)-max(plusA))
	pqA=str(max(plusA))+"-"+str(min(minusA))
	pqAA='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-1000), str(min(minusA)+1000), [max(plusA), min(minusA)])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA))+'-'+'{:,}'.format(min(minusA))+'")'
	szB=str(min(minusB)-max(plusB))
	pqB=str(max(plusB))+"-"+str(min(minusB))
	pqBB='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusB)-1000), str(min(minusB)+1000), [max(plusB), min(minusB)])+'";"'+"[GRCh38] chr"+getch(temp[2])+":g."+'{:,}'.format(max(plusB))+"-"+'{:,}'.format(min(minusB))+'")'
	sz=str(min(minusB)-max(plusA))
	pqtot='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-1000), str(min(minusB)+1000), [max(plusA), min(minusB)])+'";"'+'chr'+getch(temp[1])+':'+'{:,}'.format(max(plusA))+'-'+'{:,}'.format(min(minusB))+'")'
	if old_new=="old":
		merg=read_bds(merged_bd, False)
		mm=make_overlap([max(plusA), min(minusB)], merg[seps(temp[1], 0)], False)
		dgrc=selecta(data_dic, line[42], case, caselist)
	else:
		mm=line[24:26]
		dgrc=selecta(data_dic, line[30], case, caselist)
	merged="\t".join(mm)
	cla=make_cat_freq(mm[1], dgrc)
	GenArch, cl, firstA, secondA, firstB, secondB, high_ph, pval_ph =get_genes(acmg, clg, ddg2p, panel, getch(temp[1]), getch(temp[1]), pqA, pqB, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, "inv", False, gwas, gtexx)
	out.write("\t".join(["seq[GRCh38] inv("+seps(temp[1], 0)+")("+seps(temp[1], 1)+";"+seps(temp[2], 1)+")", pqtot, pqAA, pqBB, "NA", merged, dgrc, sz,"NA", GenArch, "confirmed", cla, cl,"impact rating", high_ph, pval_ph, firstA, firstB, "NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF", secondA, secondB,"NF","NF","NF","NF","NF","NF", pq_nomenclature(seps(temp[2], 0), seps(temp[1], 0), seps(temp[1], 1), seps(temp[2], 1), pqA, pqB, False, False, False, False, False, "inv", ""), "NF", temp[0], str(len(temp[3])),"\n"]))
	#out.write("\t".join(["seq[GRCh38] inv("+seps(temp[1], 0)+")("+seps(temp[1], 1)+";"+seps(temp[2], 1)+")", "chr"+getch(temp[1])+":"+pqtot,"[GRCh38] chr"+getch(temp[1])+":g."+pqAA+"\t[GRCh38] chr"+getch(temp[2])+":g."+pqBB,  sz,merged,"NF", "confirmed", cla, dgrc, cl, firstA, firstB, secondA, secondB, pq_nomenclature(seps(temp[2], 0), seps(temp[1], 0), seps(temp[1], 1), seps(temp[2], 1), pqA, pqB, False, False, False, False, False, "inv", ""), "NF", temp[0], str(len(temp[3])),"\n"]))
	out.close()
	
#cl, firstA, secondA, firstB, secondB, firstC, secondC =get_genes_with_C3(acmg, clg, ddg2p, panel, ch, pqA, pqB, pqC, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, ttype, False)
#out.write("\t".join([name, regt, regA, regB, regC, merged,dgrc, sz,"NF", "", "",cla, cl,"","high phenoscore","pvalue", firstA, firstB, firstC, secondA, secondB, secondC, "NA","NA",cl,pr,"NF", "\n"]))
	

def read_clust_trans(v, acmgg, ddg2p, cling, panel,reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, synd, clin, inf, old_new, cncr, terms, gwas, gtexx):
	syndromes=read_syndromes(synd)
	clinical_exome=read_syndromes(clin)
	infertility=read_syndromes(inf)
	acmg=read_acmg(acmgg)
	temp=[]
	dd2p=get_dd2p(ddg2p)
	clg=clingen(cling)
	out=open("tab4", "a")
	for i in v:
		line=i.split("\t")
		if temp==[]:
			temp.append(line[1]) #cluster ID 0
			temp.append(line[4]) #chrA 1
			temp.append(line[5]) #chrB 2
			temp.append([int(line[6])]) #cordA 3
			temp.append([int(line[7])]) #cordB 4
			temp.append([line[12]]) #strand A 5
			temp.append([line[13]]) #strand B 6
		else:
			temp[3].append(int(line[6]))#cordA 3
			temp[4].append(int(line[7]))#cordB 4
			temp[5].append(line[12])#strandA 5
			temp[6].append(line[13])#stransB 6
	plusA, minusA=dividebystrand(temp[3], temp[5])
	plusB, minusB=dividebystrand(temp[4], temp[6])
	sza=str(min(minusA)-max(plusA))
	szb=str(min(minusB)-max(plusB))
	pqA=str(max(plusA))+"-"+str(min(minusA))
	pqB=str(max(plusB))+"-"+str(min(minusB))
	pqAA='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-1000), str(min(minusA)+1000), [max(plusA), min(minusA)])+'";"'+'[GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA))+'-'+'{:,}'.format(min(minusA))+'")'
	pqBB='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[2]), str(max(plusB)-1000), str(min(minusB)+1000), [max(plusB), min(minusB)])+'";"'+"[GRCh38] chr"+getch(temp[2])+":g."+'{:,}'.format(max(plusB))+"-"+'{:,}'.format(min(minusB))+'")'
	cla=make_cat_freq("NF", "NF")
	GenArch, cl, firstA, secondA, firstB, secondB, high_ph, pval_ph =get_genes(acmg, clg, ddg2p, panel, getch(temp[1]), getch(temp[2]), pqA, pqB, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, "trans", False, gwas, gtexx)
	#out.write("\t".join(["seq[GRCh38] t("+seps(temp[1], 0)+";"+seps(temp[2], 0)+")("+seps(temp[1], 1)+";"+seps(temp[2], 1)+")","NA", "[GRCh38] chr"+getch(temp[1])+":g."+pqAA+"\t[GRCh38] chr"+getch(temp[2])+":g."+pqBB, sza, szb, "NF","NF",  "confirmed", cla, "NF", cl, firstA, firstB, secondA, secondB,pq_nomenclature(seps(temp[1], 0), seps(temp[2], 0), seps(temp[1], 1), seps(temp[2], 1), pqA, pqB, False, False, False, False, False, "trans", seps(temp[1], 0)), pq_nomenclature(seps(temp[1], 0), seps(temp[2], 0), seps(temp[1], 1), seps(temp[2], 1), pqA, pqB, False, False, False, False, False, "trans", seps(temp[2], 0)), temp[0], str(len(temp[3])), "\n"]))
	out.write("\t".join(["seq[GRCh38] t("+seps(temp[1], 0)+";"+seps(temp[2], 0)+")("+seps(temp[1], 1)+";"+seps(temp[2], 1)+")","NA", pqAA, pqBB,"NA", "NF","NF","NF",  sza, szb, GenArch, "confirmed", cla, cl, "impact rating", high_ph, pval_ph, firstA, firstB, "NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF", secondA, secondB, "NF","NF","NF","NF","NF","NF", pq_nomenclature(seps(temp[1], 0), seps(temp[2], 0), seps(temp[1], 1), seps(temp[2], 1), pqA, pqB, False, False, False, False, False, "trans", seps(temp[1], 0)), pq_nomenclature(seps(temp[1], 0), seps(temp[2], 0), seps(temp[1], 1), seps(temp[2], 1), pqA, pqB, False, False, False, False, False, "trans", seps(temp[2], 0)), temp[0], str(len(temp[3])), "\n"]))

	
def read_clust_ins(v, acmgg, ddg2p, cling, panel, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, merged_bd, synd, clin, inf, old_new, cncr, terms, data_dic, case, caselist, gwas, gtexx):
	syndromes=read_syndromes(synd)
	clinical_exome=read_syndromes(clin)
	infertility=read_syndromes(inf)
	acmg=read_acmg(acmgg)
	temp=[]
	dd2p=get_dd2p(ddg2p)
	clg=clingen(cling)
	out=open("tab4", "a")
	for i in v:
		line=i.split("\t")
		if temp==[]:
			temp.append(line[1])#cluster id 0
			temp.append(line[4])#chrA 1
			temp.append(line[5])#chrB 2
			temp.append([int(line[6])])#cordA 3
			temp.append([int(line[7])])#cordB 4
			temp.append([line[12]])#strandA 5
			temp.append([line[13]])#strandB 6
		else:
			temp[3].append(int(line[6]))
			temp[4].append(int(line[7]))
			temp[5].append(line[12])
			temp[6].append(line[13])
	plusA, minusA=dividebystrand(temp[3], temp[5])
	plusB, minusB=dividebystrand(temp[4], temp[6])
	a1="-/-"
	a2="+/+"
	if temp[5][0]!=temp[6][0]:
		a1="-/+"
		a2="+/-"
	if get_Overlap(plusA,minusA)>0 or max(minusA)<min(plusA):
		pqA=str(max(plusB))+"-"+str(min(minusB))
		pqAA='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[2]), str(max(plusB)-1000), str(min(minusB)+1000), [max(plusB), min(minusB)])+'";"'+"Sink [GRCh38] chr"+getch(temp[2])+":g."+'{:,}'.format(max(plusB))+"-"+'{:,}'.format(min(minusB))+'")'
		pqB=str(min(plusA+minusA)-25)+"-"+str(max(plusA+minusA)+25)
		pqBB='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(min(plusA+minusA)-1025), str(max(plusA+minusA)+1025), [min(plusA+minusA)-25, max(plusA+minusA)+25])+'";"'+'Source [GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(min(plusA+minusA)-25)+'-'+'{:,}'.format(max(plusA+minusA)+25)+'")'
		if old_new=="old":
			merg=read_bds(merged_bd, False)
			mm=make_overlap([max(plusB), min(minusB)], merg[seps(temp[2], 0)], "ins")
			if len(line)>36:
				dgrc=selecta(data_dic, line[43], case, caselist)
			else:
				dgrc=selecta(data_dic, line[30], case, caselist)
		else:
			if len(line)>36:
				mm=line[26:28]
				dgrc=selecta(data_dic, line[31], case, caselist)
			else:
				mm=line[24:26]
				dgrc=selecta(data_dic, line[26], case, caselist)
		merged="\t".join(mm)
		cla=make_cat_freq(mm[1], dgrc)
		GenArch, cl, firstA, secondA, firstB, secondB, high_ph, pval_ph =get_genes(acmg, clg, ddg2p, panel, getch(temp[2]), getch(temp[1]), pqA, pqB, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, "ins", True, gwas, gtexx)
		#print("seq[GRCh38] ins("+seps(temp[2], 0)+";"+seps(temp[1], 0)+")("+seps(temp[2], 1)+";"+seps(temp[1], 1)+")", "Source [GRCh38] chr"+getch(temp[1])+":g."+pqB+"\tSink [GRCh38] chr"+getch(temp[2])+":g."+pqA, merged, str(max(plusA+minusA)-min(plusA+minusA)), str(min(minusB)-max(plusB)), "confirmed", cla, dgrc, cl, firstB, firstA, secondB, secondA,  pq_nomenclature(seps(temp[2], 0), seps(temp[1], 0), seps(temp[2], 1), seps(temp[1], 1), pqA, pqB, False, False, False, False, False, "ins", "alal"), pq_nomenclature(seps(temp[2], 0), seps(temp[1], 0), seps(temp[2], 1), seps(temp[1], 1), pqA, pqB, False, False, False, False, False, "ins", seps(temp[2], 0)), temp[0], str(len(temp[3])))
		out.write("\t".join(["seq[GRCh38] ins("+seps(temp[2], 0)+";"+seps(temp[1], 0)+")("+seps(temp[2], 1)+";"+seps(temp[1], 1)+")", pqBB, pqBB, pqAA, "NA", merged, dgrc, str((max(plusA+minusA)+25)-(min(plusA+minusA)-25)), str(min(minusB)-max(plusB)), GenArch, "confirmed", cla, cl, "impact rating", high_ph, pval_ph, firstB, firstA, "NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF", secondB, secondA,"NF","NF","NF","NF","NF","NF", pq_nomenclature(seps(temp[2], 0), seps(temp[1], 0), seps(temp[2], 1), seps(temp[1], 1), pqA, pqB, False, False, False, False, False, "ins", "alal"), pq_nomenclature(seps(temp[2], 0), seps(temp[1], 0), seps(temp[2], 1), seps(temp[1], 1), pqA, pqB, False, False, False, False, False, "ins", seps(temp[2], 0)), temp[0], str(len(temp[3])), "\n"]))
	else:
		if old_new=="old":
			merg=read_bds(merged_bd, False)
			mm=make_overlap([max(plusA), min(minusA)], merg[seps(temp[1], 0)], "ins")
			if len(line)>36:
				dgrc=selecta(data_dic, line[43], case, caselist)
			else:
				dgrc=selecta(data_dic, line[30], case, caselist)
		else:
			if len(line)>36:
				mm=line[26:28]
				dgrc=selecta(data_dic, line[31], case, caselist)
			else:
				mm=line[24:26]
				dgrc=selecta(data_dic, line[26], case, caselist)
		merged="\t".join(mm)
		cla=make_cat_freq(mm[1], dgrc)
		pqA=str(max(plusA))+"-"+str(min(minusA))
		pqAA='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[1]), str(max(plusA)-1000), str(min(minusA)+1000), [max(plusA), min(minusA)])+'";"'+'Sink [GRCh38] chr'+getch(temp[1])+':g.'+'{:,}'.format(max(plusA))+'-'+'{:,}'.format(min(minusA))+'")'
		pqB=str(min(plusB+minusB)-25)+"-"+str(max(plusB+minusB)+25)
		pqBB='=HIPERLIGAÇÃO("'+make_ucsc_link('hg38', getch(temp[2]), str(min(plusB+minusB)-1025), str(max(plusB+minusB)+1025), [min(plusB+minusB)-25, max(plusB+minusB)+25])+'";"'+"Source [GRCh38] chr"+getch(temp[2])+":g."+'{:,}'.format(min(plusB+minusB)-25)+"-"+'{:,}'.format(max(plusB+minusB)+25)+'")'
		GenArch, cl, firstA, secondA, firstB, secondB, high_ph, pval_ph =get_genes(acmg, clg, ddg2p, panel, getch(temp[1]), getch(temp[2]), pqA, pqB, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex, syndromes, clinical_exome, infertility, cncr, terms, "ins", True, gwas, gtexx)
		out.write("\t".join(["seq[GRCh38] ins("+seps(temp[1], 0)+";"+seps(temp[2], 0)+")("+seps(temp[1], 1)+";"+seps(temp[2], 1)+")",pqBB, pqBB, pqAA, "NA", merged, dgrc, str((max(plusB+minusB)+25)-(min(plusB+minusB)-25)), str(min(minusA)-max(plusA)), GenArch, "confirmed", cla, cl,"impact rating", high_ph, pval_ph,firstB, firstA, "NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF", secondB, secondA,"NF","NF","NF","NF","NF","NF", pq_nomenclature(seps(temp[1], 0), seps(temp[2], 0), seps(temp[1], 1), seps(temp[2], 1), pqA, pqB, False, False, False, False, False, "ins", "alala"), pq_nomenclature(seps(temp[1], 0), seps(temp[2], 0), seps(temp[1], 1), seps(temp[2], 1), pqA, pqB, False, False, False, False, False, "ins", seps(temp[1], 0)),temp[0], str(len(temp[3])), "\n"]))
	out.close()

def read_dgrc_bd(is_ins, infile, caselist):
	f=open(infile)
	dd={}
	for i in f:
		line=i.split("\t")
		if is_ins==True:
			dd[line[0]]=line[15]
		else:
			dd[line[0]]=line[10]
	f.close()
	caselists=[]
	f=open(caselist)
	for i in f:
		caselists.append(i.strip())
	return dd, caselists
	
def selecta(data_dic, DGRC, case, caselist): 
	freq=""
	#print("first", DGRC, case, caselist)
	if DGRC in data_dic:
		cases=data_dic[DGRC].split(",")
		if case in cases and len(cases)>1:
			freq=" ("+str(len(cases)-1)+"/"+str(len(caselist)-1)+")"
		elif case in cases and len(cases)==1:
			DGRC="NF"
		elif case not in cases and case in caselist:
			freq=" ("+str(len(cases))+"/"+str(len(caselist)-1)+")"		
		else:
			freq=" ("+str(len(cases))+"/"+str(len(caselist))+")"
	else:
		DGRC="NF"
	#print("second", DGRC, freq)
	return DGRC+freq

	
def pq_nomenclature(chrA, chrB, bandA, bandB, pqAA, pqBB, isspec, isdelA, isdupA, isdelB, isdupB, typee, der):
	"""makes the breakpoint nomenclature according to ISCN2020"""
	dic_chr_ref={"1":"NC_000001.11", "2":"NC_000002.12", "3":"NC_000003.12", "4":"NC_000004.12", "5":"NC_000005.10", "6":"NC_000006.12", "7":"NC_000007.14", "8":"NC_000008.11", "9":"NC_000009.12", "10":"NC_000010.11", "11":"NC_000011.10", "12":"NC_000012.12", "13":"NC_000013.11", "14":"NC_000014.9", "15":"NC_000015.10", "16":"NC_000016.10", "17":"NC_000017.11", "18":"NC_000018.10", "19":"NC_000019.10", "20":"NC_000020.11", "21":"NC_000021.9", "22":"NC_000022.11", "X":"NC_000023.11", "Y":"NC_000024.10"}
	pqA=pqAA.strip()
	pqB=pqBB.strip()
	singA=pqA
	singB=pqB
	if typee=="ins":#chrA is the sink cromossome and chrB the source cromossome, in insertion
		if der==chrA:
			return dic_chr_ref[chrA.replace("chr","")]+":g.{"+put_comma(pqA)+"}ins["+dic_chr_ref[chrB.replace("chr","")]+":g.{"+put_comma(pqB)+"}]"
		else:
			return dic_chr_ref[chrB.replace("chr","")]+":g.{"+put_comma(pqB)+"}del"
	elif isspec==True:
		return chrA+bandA+": g.[{"+put_comma(pqAA)+"-"+put_comma(pqBB)+"}]"
	aa=pqA.split("-")
	bb=pqB.split("-")
	if typee=="del":
		return dic_chr_ref[chrA.replace("chr","")]+": g.{"+put_comma(pqA)+"_"+put_comma(pqB)+"}del"
	elif typee=="dup":
		return dic_chr_ref[chrA.replace("chr","")]+": g.{"+put_comma(pqA)+"_"+put_comma(pqB)+"}dup"
	elif typee=="inv":
		return dic_chr_ref[chrA.replace("chr","")]+":g.[{"+put_comma(pqA)+"_"+put_comma(pqB)+"}inv]"
	elif typee=="trans":
		if der==chrA:
			if "p" in bandA and "p" in bandB:
				cc="der("+chrA+") "+dic_chr_ref[chrA.replace("chr","")]+":g.pter_{"+put_comma(aa[0])+"}delins["+dic_chr_ref[chrB.replace("chr","")]+":g.pter_{"+put_comma(bb[-1])+"}]"
			if "p" in bandA and "q" in bandB:
				cc="der("+chrA+") "+dic_chr_ref[chrA.replace("chr","")]+":g.pter_{"+put_comma(aa[0])+"}delins["+dic_chr_ref[chrB.replace("chr","")]+":g.{"+put_comma(bb[-1])+"}_qterinv]"
			if "q" in bandA and "p" in bandB:
				cc="der("+chrA+") "+dic_chr_ref[chrA.replace("chr","")]+":g.{"+put_comma(aa[0])+"}_qterdelins["+dic_chr_ref[chrB.replace("chr","")]+":g.{"+put_comma(bb[-1])+"}_pterinv]"
			if "q" in bandA and "q" in bandB:
				cc="der("+chrA+") "+dic_chr_ref[chrA.replace("chr","")]+":g.{"+put_comma(aa[0])+"}_qterdelins["+dic_chr_ref[chrB.replace("chr","")]+":g.{"+put_comma(bb[-1])+"}_qter]"
		if der==chrB:
			if "p" in bandA and "p" in bandB:
				cc="der("+chrB+") "+dic_chr_ref[chrB.replace("chr","")]+":g.pter_{"+put_comma(bb[0])+"}delins["+dic_chr_ref[chrA.replace("chr","")]+":g.pter_{"+put_comma(aa[-1])+"}]"
			if "p" in bandA and "q" in bandB:
				cc="der("+chrB+") "+dic_chr_ref[chrB.replace("chr","")]+":g.{"+put_comma(bb[0])+"}_qterdelins["+dic_chr_ref[chrA.replace("chr","")]+"g.{"+put_comma(aa[-1])+"}_pterinv]"
			if "q" in bandA and "p" in bandB:
				cc="der("+chrB+") "+dic_chr_ref[chrB.replace("chr","")]+":g.pter_{"+put_comma(bb[0])+"}delins["+dic_chr_ref[chrA.replace("chr","")]+"g.{"+put_comma(aa[-1])+"}_qterinv]"
			if "q" in bandA and "q" in bandB:
				cc="der("+chrB+") "+dic_chr_ref[chrB.replace("chr","")]+":g.{"+put_comma(bb[0])+"}_qterdelins["+dic_chr_ref[chrA.replace("chr","")]+"g.{"+put_comma(aa[-1])+"}_qter]"
		return cc

def put_comma(v):#f'{value:,}
	if "-" in v:
		y=v.replace(",","")
		aa=y.split("-")
		return f'{int(float(aa[0])):,}'+"_"+f'{int(float(aa[1])):,}'
	elif "_" in v and "del" in v:
		y=v.replace(",","")
		u=y.replace("del","")
		aa=u.split("_")
		return f'{int(float(aa[0])):,}'+"_"+f'{int(float(aa[1])):,}'+"del"
	else:
		y=v.replace(",","")
		return f'{int(float(y)):,}'

def seps(f, v):
	if "q" in f:
		bb=f.split("q")[v]
		if v==1:
			aa="q"+bb
		else:
			aa=bb
	if "p" in f:
		bb=f.split("p")[v]
		if v==1:
			aa="p"+bb
		else:
			aa=bb
	return aa

def get_freq(gg):
	d=0.0
	de="NF"
	is_it=False
	for el in gg.split("\t"):
		if el!="na" and "%" in el:
			is_it=True
			if float(el.replace("%",""))>d:
				d=float(el.replace("%",""))
		elif el!="na" and "%" not in el:
			dgrc=el.strip()
	if is_it==True:
		de=dgrc
	else:
		de="NF"
	return de

def getch(c):
	aa=c.split("p")
	if len(aa)==1:
		aa=c.split("q")
	return aa[0]

def see_dd(dd):
	l=[]
	for el in dd:
		if "(" in el:
			l.append(el)
	return l

def edit_second(secondA, secondB, secondC):
	sep_secondA=secondA.split("\t")
	sep_secondB=secondB.split("\t")
	aa=0
	gA="NF"#genome architecture item
	if secondC!="":
		sep_secondC=secondC.split("\t")
		while(aa<2):
			if len(set(sep_secondA[aa].split(";")) & set(sep_secondB[aa].split(";")) & set(sep_secondC[aa].split(";")))>0:
				l=set(sep_secondA[aa].split(";")) & set(sep_secondB[aa].split(";")) & set(sep_secondC[aa].split(";"))
				l=list(map(lambda x: x.replace('NF', ''), l))
				l= list(filter(None, l))
				if len(l)>0:
					if gA=="NF":
						gA=l[0]
					sep_secondA[aa]="*"+sep_secondA[aa]
					sep_secondB[aa]="*"+sep_secondB[aa]
					sep_secondC[aa]="*"+sep_secondC[aa]
			elif len(set(sep_secondA[aa].split(";")) & set(sep_secondB[aa].split(";")))>0:
				l=set(sep_secondA[aa].split(";")) & set(sep_secondB[aa].split(";"))
				l=list(map(lambda x: x.replace('NF', ''), l))
				l= list(filter(None, l))
				if len(l)>0:
					if gA=="NF":
						gA=l[0]
					sep_secondA[aa]="*"+sep_secondA[aa]
					sep_secondB[aa]="*"+sep_secondB[aa]
			elif len(set(sep_secondA[aa].split(";")) & set(sep_secondC[aa].split(";")))>0:
				l=set(sep_secondA[aa].split(";")) & set(sep_secondC[aa].split(";"))
				l=list(map(lambda x: x.replace('NF', ''), l))
				l= list(filter(None, l))
				if len(l)>0:
					if gA=="NF":
						gA=l[0]
					sep_secondA[aa]="*"+sep_secondA[aa]
					sep_secondC[aa]="*"+sep_secondC[aa]
			elif len(set(sep_secondC[aa].split(";")) & set(sep_secondB[aa].split(";")))>0:
				l=set(sep_secondC[aa].split(";")) & set(sep_secondB[aa].split(";"))
				l=list(map(lambda x: x.replace('NF', ''), l))
				l= list(filter(None, l))
				if len(l)>0:
					if gA=="NF":
						gA=l[0]
					sep_secondC[aa]="*"+sep_secondC[aa]
					sep_secondB[aa]="*"+sep_secondB[aa]
			aa+=1
		return gA, "\t".join(sep_secondA), "\t".join(sep_secondB), "\t".join(sep_secondC)
	else:
		while(aa<2):
			if len(set(sep_secondA[aa].split(";")).intersection(set(sep_secondB[aa].split(";"))))>0:
				l=set(sep_secondA[aa].split(";")) & set(sep_secondB[aa].split(";"))
				l=list(map(lambda x: x.replace('NF', ''), l))
				l= list(filter(None, l))
				if len(l)>0:
					if gA=="NF":
						gA=l[0]
					sep_secondA[aa]="*"+sep_secondA[aa]
					sep_secondB[aa]="*"+sep_secondB[aa]
			aa+=1
		return gA, "\t".join(sep_secondA), "\t".join(sep_secondB)


def get_other_ops(chra, chrb, bp1, bp2, bp3, Conserved_reg, firstA, secondA, firstB, secondB, firstC, secondC):
	isfA, ccA=treat_genes_for_table.check_2_5kb(chra, int(bp1.split("-")[0]), int(bp1.split("-")[1]))
	isfB, ccB=treat_genes_for_table.check_2_5kb(chrb, int(bp2.split("-")[0]), int(bp2.split("-")[1]))
	cc=[]
	#print("lalala",ccA, ccB, isfA, isfB)
	if ccA!="":
		cc.append(ccA)
	if ccB!="":
		cc.append(ccB)
	if bp3!="":
		isfC, ccC=treat_genes_for_table.check_2_5kb(chra, int(bp3.split("-")[0]), int(bp3.split("-")[1]))
		if ccC!="":
			cc.append(ccC)
		if ccA=="" and ccB=="" and ccC=="":
			if compare_cncr(Conserved_reg, chra, [int(bp1.split("-")[0]), int(bp1.split("-")[1])])==False and compare_cncr(Conserved_reg, chrb, [int(bp2.split("-")[0]), int(bp2.split("-")[1])])==False and compare_cncr(Conserved_reg, chra, [int(bp3.split("-")[0]), int(bp3.split("-")[1])])==False:
				return "10 - Intergenic variant", "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(firstC), secondC
			else:
				return "9 - Conserved Region", "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(firstC), secondC
		elif isfA==True or isfB==True or isfC==True:
			return "7 - Promoter region: "+"& ".join(cc), "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(firstC), secondC
		else:
			return "8 - 2.5kb region: "+"& ".join(cc), "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(firstC), secondC
	else:
		if ccA=="" and ccB=="":
			if compare_cncr(Conserved_reg, chra, [int(bp1.split("-")[0]), int(bp1.split("-")[1])])==False and compare_cncr(Conserved_reg, chrb, [int(bp2.split("-")[0]), int(bp2.split("-")[1])])==False:
				return "10 - Intergenic variant", "\t".join(firstA), secondA, "\t".join(firstB), secondB
			else:
				return "9 - Conserved Region", "\t".join(firstA), secondA, "\t".join(firstB), secondB
		elif isfA==True or isfB==True:
			print(firstA)
			return "7 - Promoter region: " + "& ".join(cc), "\t".join(firstA), secondA, "\t".join(firstB), secondB
		else:
			return "8 - 2.5kb region: "+"& ".join(cc), "\t".join(firstA), secondA, "\t".join(firstB), secondB

def get_genes(acmg, clingenn, ddg2p, pannel, chra, chrb, bp1, bp2, reps, segdup, genehancer, loops, tadfile, poesf, oe, gtex_file, synd, clinicalexome, infertility, cncr, terms, tt, is_ins, gwas, gtexx):
	dic_gtex=read_gtex(gtex_file)
	Conserved_reg=read_cncr(cncr)
	tadA_genes, tadA4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chra, int(bp1.split("-")[0]))
	tadB_genes, tadB4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chrb, int(bp2.split("-")[0]))
	gA, biotypeA=treat_genes_for_table.getdisrp(chra, int(bp1.split("-")[0]), int(bp1.split("-")[1]))
	ddA=get_oe_for_tad(gA, oe)
	gB, biotypeB=treat_genes_for_table.getdisrp(chrb, int(bp2.split("-")[0]), int(bp2.split("-")[1]))
	ddB=get_oe_for_tad(gB, oe)
	ta=get_oe_for_tad(tadA_genes, oe)
	tb=get_oe_for_tad(tadB_genes, oe)	
	noncodA=[["NF","NF","NF"]]
	noncodB=[["NF","NF","NF"]]
	cpsA=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	cpsB=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	phenscoreA=""
	phenscoreB=""
	if ddA==["NF"]:
		firstA=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondA=repeats_segdups(reps, chra, bp1)+"\t"+repeats_segdups(segdup, chra, bp1)+"\t"+GeneHancer(tadA_genes, genehancer, chra, bp1)+"\t"+getloops(loops, chra, bp1)+"\t"+pos_effect(tadA4_genes, poesf)+"\t"+"& ".join(ta)
	if ddB==["NF"]:
		firstB=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondB=repeats_segdups(reps, chrb, bp2)+"\t"+repeats_segdups(segdup, chrb, bp2)+"\t"+GeneHancer(tadB_genes, genehancer, chrb, bp2)+"\t"+getloops(loops, chrb, bp2)+"\t"+pos_effect(tadB4_genes, poesf)+"\t"+"& ".join(tb)
	GenArch, secondA, secondB=edit_second(secondA, secondB, "")
	if ddA==["NF"] and ddB==["NF"]:
		cl, firstA, secondA, firstB, secondB = get_other_ops(chra, chrb, bp1, bp2, "", Conserved_reg, firstA, secondA, firstB, secondB, "", "")
		return GenArch, cl, firstA, secondA, firstB, secondB, "NF", "NF"
	if ddA!=["NF"]:
		genesA=str(len(set(ddA)))
		j1A, j2A, j3A, j4A, j5A, j6A=[['NF'] for i in range(6)]
		for el in ddA:
			if el.split(" ")[0] in acmg:
				if "NF" in j1A: j1A.remove("NF")
				j1A.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2A: j2A.remove("NF")
				j2A.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3A: j3A.remove("NF")
				j3A.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4A: j4A.remove("NF")
				j4A.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5A: j5A.remove("NF")
				j5A.append(gwas[el.split(" ")[0]])#gwas
			if el.split(" ")[0] in gtexx:
				if "NF" in j6A: j6A.remove("NF")
				j6A.append(gtexx[el.split(" ")[0]])#gtex
		freqA, genecardA=get_oe(ddA)
		cpsA, noncodA, phenscoreA=get_phen(ddA, clingenn, ddg2p, pannel, chra, str(bp1), "", dic_gtex, terms)
		print(genesA, ddA, biotypeA, freqA, genecardA, cpsA[1],cpsA[2],phenscoreA,cpsA[3],cpsA[4],cpsA[5],cpsA[7],cpsA[0],cpsA[8], j5A, j6A, j3A, j4A,j2A,j1A,cpsA[6])
		firstA=[genesA,"& ".join(ddA),",".join(biotypeA), "",freqA, genecardA, cpsA[1],cpsA[2],phenscoreA,cpsA[3],cpsA[4],cpsA[5],cpsA[7],cpsA[0],cpsA[8], "& ".join(j5A), "& ".join(j6A), "& ".join(j3A),"& ".join(j4A),"& ".join(j2A),"& ".join(j1A),cpsA[6]]
	if ddB!=["NF"]:
		genesB=str(len(set(ddB)))
		j1B, j2B, j3B, j4B, j5B, j6B=[['NF'] for i in range(6)]
		for el in ddB:
			if el.split(" ")[0] in acmg:
				if "NF" in j1B: j1B.remove("NF")
				j1B.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2B: j2B.remove("NF")
				j2B.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3B: j3B.remove("NF")
				j3B.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4B: j4B.remove("NF")
				j4B.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5B: j5B.remove("NF")
				j5B.append(gwas[el.split(" ")[0]])#gwas
			if el.split(" ")[0] in gtexx:
				if "NF" in j6B: j6B.remove("NF")
				j6B.append(gtexx[el.split(" ")[0]])#gtex		
		freqB, genecardB=get_oe(ddB)
		if tt!="ins":
			cpsB, noncodB, phenscoreB=get_phen(ddB, clingenn, ddg2p, pannel, chrb, str(bp2), "", dic_gtex, terms)
		else:
			cpsB, noncodB, phenscoreB=get_phen_ins(ddB, clingenn, ddg2p, pannel, chrb, str(bp2), "", dic_gtex, terms)
		firstB=[genesB,"& ".join(ddB), ",".join(biotypeB), "",freqB, genecardB, cpsB[1],cpsB[2],phenscoreB,cpsB[3],cpsB[4],cpsB[5], cpsB[7],cpsB[0],cpsB[8], "& ".join(j5B), "& ".join(j6B), "& ".join(j3B),"& ".join(j4B),"& ".join(j2B),"& ".join(j1B),cpsB[6]]
	#print("get geneclass", "ddA", ddA, "cpsA", cpsA, "noncodA",noncodA, "ddB",ddB, "cpsB", cpsB, "noncodB",noncodB, "0", 0, "is_ins", is_ins)
	#cl, afectA, afectB = get_gene_class(ddA, cpsA, noncodA, ddB, cpsB, noncodB, 0, is_ins)###calcula em que grau de prioritizacao é que aquilo fica
	cl, afectA, afectB= classes_for_cpx.interpret_to_class(tt, classes_for_cpx.get_gene_class(ddA, cpsA, noncodA), classes_for_cpx.get_gene_class(ddB, cpsB, noncodB), "", [chra, int(bp1.split("-")[0]), int(bp1.split("-")[1]), chrb, int(bp2.split("-")[0]), int(bp2.split("-")[1])])###calcula em que grau de prioritizacao é que aquilo fica
	if firstA[0]!="0":
		#firstA[2]=afectA
		firstA[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chra+'%3A'+bp1+'";"'+afectA+'")'
	if firstB[0]!="0":
		#firstB[2]=afectB
		firstB[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chrb+'%3A'+bp2+'";"'+afectB+'")'
	high_ph, pval_ph=get_highest_ph(phenscoreA, phenscoreB, "")
	if cl=="NF":
		cl, firstA, secondA, firstB, secondB = get_other_ops(chra, chrb, bp1, bp2, "", Conserved_reg, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondA, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondB, "", "")
		return GenArch, cl, firstA, secondA, firstB, secondB, high_ph, pval_ph
	else:
		return GenArch, cl, "\t".join(firstA), secondA, "\t".join(firstB), secondB, high_ph, pval_ph

def separe_by_TAD(dd, tadA_genes, tadB_genes):
	ddA=[]
	ddB=[]
	for el in dd:
		if el.split("(")[0] in tadA_genes:
			ddA.append(el)
		if el.split("(")[0] in tadB_genes:
			ddB.append(el)
	return ddA, ddB

def pos_effect(tad4, genes_list):
	f=open(genes_list)
	b=["NF"]
	for i in f:
		line=i.split("\t")
		if line[0] in tad4:
			if b==["NF"]:
				b=[]
			b.append(line[0]+"_("+line[4]+" cases)")
	f.close()
	return "& ".join(b)
	
def repeats_segdups(infile, chh, regg):
	"""vai buscar as repetiçoes e as duplicações segmentais que se sobrepoem"""
	reg=regg.split("-")
	f=open(infile)
	aa="NF"
	for i in f:
		line=i.split("\t")
		if line[0][3:]==chh:
			if get_Overlap(reg,[int(line[1]), int(line[2])])>0:
				if line[3].startswith("chr"):
					if aa=="NF":
						aa=line[3].strip()
					else:
						aa+="; "+line[3].strip()
				else:
					if aa=="NF":
						aa=line[3]+":"+line[4]+":"+line[5].strip()
					else:
						aa+="; "+line[3]+":"+line[4]+":"+line[5].strip()
	f.close()
	return aa
			
def getloops(infile, chh, regg):
	"""vai buscar os loops que incluem genes da brTAD que são disrumpidos"""
	f=open(infile)
	aa=["NF"]
	reg=regg.split("-")
	for i in f:
		line=i.split("\t")
		if line[0]==chh:
			if get_Overlap(reg,[int(line[1]), int(line[2])])>0 or get_Overlap(reg,[int(line[4]), int(line[5])])>0:
				if aa[0]=="NF":
					aa=[]
				aa.append(line[-1].strip())
	f.close()
	return "& ".join(aa)

def GeneHancer(lgenes, infile, chh, regg):
	"""vai buscar os clusters de GeneHancer de genes da brTAD que são disrumpidos"""	
	f=open(infile)
	aa=["NF"]
	reg=regg.split("-")
	for i in f:
		line=i.split("\t")
		if line[1][3:]==chh:
			if get_Overlap(reg,[int(line[2]), int(line[3].strip())])>0 and line[0] in lgenes:
				if aa[0]=="NF":
					aa=[]
				aa.append(line[0]+" ("+line[1]+":"+line[2]+"-"+line[3].strip()+")")
	f.close()
	return ";".join(aa)
	
def filt_genes(k):
	gg=""
	aa=k.split(";")
	for el in aa:
		if len(gg)>0:
			gg+="; "
		if ":" in el:
			gg+=el.split(":")[0]
		else:
			gg+=el.strip("\n")
	return gg

def get_oe_for_tad(dd, oe):
	f=open(oe)
	ndd=[]
	dones=[]
	for el in f:
		if el.split("\t")[0] in dd:
			if "NF" in el.split("\t")[1] or  "NA" in el.split("\t")[1]:
				vv="NF"
			else:
				vv=el.split(" ")[-2]
			ndd.append(el.split("\t")[0]+" ("+vv+")")
			dones.append(el.split("\t")[0])
	f.close()
	if len(set(dd)-set(dones))>0:
		for el in set(dd)-set(dones):
			if el!="":
				ndd.append(el+" (ND)")
	if len(ndd)==0:
		ndd=["NF"]
	return ndd

def read_gtex(infile):
	f=open(infile)
	dic={}
	for i in f:
		line=i.split("\t")
		dic[line[0]]=line[1]
	f.close()
	return dic


def read_cncr(cncr):
	f=open(cncr)
	dic={}
	for i in f:
		line=i.split("\t")
		if line[0] not in dic:
			dic[line[0]]=[[int(line[1]), int(line[2])]]
		else:
			dic[line[0]].append([int(line[1]), int(line[2])])
	f.close()
	return dic

def compare_cncr(cncr, chrA, rega):
	for el in cncr[chrA]:
		if get_Overlap(el, rega)>0:
			return True
	return False


def conserged_for_cnv(chra, bp1, bp2, Conserved_reg, firstA, secondA, firstB, secondB, dell):
	isf, cc=treat_genes_for_table.check_2_5kb(chra, int(bp1.split("-")[0]), int(bp2.split("-")[0]))
	if cc=="":
		if compare_cncr(Conserved_reg, chra, [bp1, bp2])==False:
			return "10 - Intergenic variant", "\t".join(firstA), secondA, "\t".join(firstB), secondB,  "\t".join(dell)
		else:
			return "9 - Conserved Region", "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(dell)
	elif isf==True:
		return "7 - Promoter region: "+cc, "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(dell)
	else:
		return "8 - 2.5kb region: "+cc, "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(dell)

def get_genes_del_dup(line, acmg, clingenn, ddg2p, pannel, chra, chrb, bp1, bp2, reps, segdup, genehancer, loops, tadfile, poesf, ttype, oe, gtex_file, synd, clinicalexome, infertility, cncr, terms, gwas, gtex):
	print(bp1, bp2)
	dic_gtex=read_gtex(gtex_file)
	Conserved_reg=read_cncr(cncr)
	tadA_genes, tadA4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chra, int(bp1.split("-")[0]))
	tadB_genes, tadB4_genes =treat_genes_for_table.get_genes_from_TAD(tadfile, chrb, int(bp2.split("-")[0]))
	gA, biotypeA=treat_genes_for_table.getdisrp(chra, int(bp1.split("-")[0]), int(bp1.split("-")[0]))
	ddA=get_oe_for_tad(gA, oe)
	gB, biotypeB=treat_genes_for_table.getdisrp(chrb, int(bp2.split("-")[0]), int(bp2.split("-")[0]))
	ddB=get_oe_for_tad(gB, oe)
	ta=get_oe_for_tad(tadA_genes, oe)
	tb=get_oe_for_tad(tadB_genes, oe)
	#print("ddA, ddB", ddA, ddB)
	noncodA=[["NF","NF","NF"]]
	noncodB=[["NF","NF","NF"]]
	cpsA=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	cpsB=["NF","NF", "NF", "NF", "NF", "NF", "NF", "NF", "NF"]
	phenscoreA=""
	phenscoreB=""
	phenscoredel=""
	if ddA==["NF"]:
		firstA=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondA=repeats_segdups(reps, chra, bp1+"-"+str(int(bp1)+5))+"\t"+repeats_segdups(segdup, chra, bp1+"-"+str(int(bp1)+5))+"\t"+GeneHancer(tadA_genes, genehancer, chra, bp1+"-"+str(int(bp1)+5))+"\t"+getloops(loops, chra, bp1+"-"+str(int(bp1)+5))+"\t"+pos_effect(tadA4_genes, poesf)+"\t"+"& ".join(ta)
	if len(line[10])<3:
		genesdel=0
		dell=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	if ddB==["NF"]:
		firstB=["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"]
	secondB=repeats_segdups(reps, chrb, bp2+"-"+str(int(bp2)+5))+"\t"+repeats_segdups(segdup, chrb, bp2+"-"+str(int(bp2)+5))+"\t"+GeneHancer(tadB_genes, genehancer, chrb, bp2+"-"+str(int(bp2)+5))+"\t"+getloops(loops, chrb, bp2+"-"+str(int(bp2)+5))+"\t"+pos_effect(tadB4_genes, poesf)+"\t"+"& ".join(tb)
	GenArch, secondA, secondB=edit_second(secondA, secondB, "")
	if ddA==["NF"] and ddB==["NF"] and len(line[10])<3:
		cl, firstA, secondA, firstB, secondB, dell=conserged_for_cnv(chra, bp1, bp2, Conserved_reg, firstA, secondA, firstB, secondB, dell)
		return GenArch, cl, firstA, secondA, firstB, secondB, dell, "NF", "NF"
	if ddA!=["NF"]:
		genesA=str(len(set(ddA)))
		j1A, j2A, j3A, j4A, j5A, j6A=[['NF'] for i in range(6)]
		for el in ddA:
			if el.split(" ")[0] in acmg:
				if "NF" in j1A: j1A.remove("NF")
				j1A.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2A: j2A.remove("NF")
				j2A.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3A: j3A.remove("NF")
				j3A.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4A: j4A.remove("NF")
				j4A.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5A: j5A.remove("NF")
				j5A.append(gwas[el.split(" ")[0]])#gwas
			if el.split(" ")[0] in gtex:
				if "NF" in j6A: j6A.remove("NF")#j6B.append(gtexx[el.split(" ")[0]])#gtex
				j6A.append(gtex[el.split(" ")[0]])#gtex
		freqA, genecardA=get_oe(ddA)
		cpsA, noncodA, phenscoreA=get_phen(ddA, clingenn, ddg2p, pannel, chra, str(bp1)+"-"+str(int(bp1)+5), "", dic_gtex, terms)
		firstA=[genesA,"& ".join(ddA),",".join(biotypeA), "",freqA, genecardA, cpsA[1],cpsA[2],phenscoreA,cpsA[3],cpsA[4],cpsA[5],cpsA[7],cpsA[0],cpsA[8], "& ".join(j5A), "& ".join(j6A),"& ".join(j3A),"& ".join(j4A),"& ".join(j2A),"& ".join(j1A), cpsA[6]]
	if ddB!=["NF"]:
		genesB=str(len(set(ddB)))
		j1B, j2B, j3B, j4B, j5B, j6B=[['NF'] for i in range(6)]
		for el in ddB:
			if el.split(" ")[0] in acmg:
				if "NF" in j1B: j1B.remove("NF")
				j1B.append(el.split("(")[0])#acmg
			if el.split(" ")[0] in infertility:
				if "NF" in j2B: j2B.remove("NF")
				j2B.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3B: j3B.remove("NF")
				j3B.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4B: j4B.remove("NF")
				j4B.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5B: j5B.remove("NF")
				j5B.append(gwas[el.split(" ")[0]])#gwas
			if el.split(" ")[0] in gtex:
				if "NF" in j6B: j6B.remove("NF")
				j6B.append(gtex[el.split(" ")[0]])#gtex			
		freqB, genecardB=get_oe(ddB)
		cpsB, noncodB, phenscoreB=get_phen(ddB, clingenn, ddg2p, pannel, chrb, str(bp2)+"-"+str(int(bp2)+5), "", dic_gtex, terms)
		firstB=[genesB,"& ".join(ddB),",".join(biotypeB), "",freqB, genecardB, cpsB[1],cpsB[2],phenscoreB,cpsB[3],cpsB[4],cpsB[5],cpsB[7],cpsB[0],cpsB[8], "& ".join(j5B),"& ".join(j6B), "& ".join(j3B),"& ".join(j4B),"& ".join(j2B),"& ".join(j1B), cpsB[6]]
	if len(line[10])>3:
		dels=line[10].split(";")[:-1]
		genesdel=str(len(set(dels)))
		j1del, j2del, j3del, j4del, j5del, j6del=[['NF'] for i in range(6)]
		for el in dels:
			if el.split(" ")[0] in acmg:
				if "NF" in j1del: j1del.remove("NF")
				j1del.append(el.split("(")[0])
			if el.split(" ")[0] in infertility:
				if "NF" in j2del: j2del.remove("NF")
				j2del.append(el.split("(")[0])#infertility
			if el.split(" ")[0] in clinicalexome:
				if "NF" in j3del: j3del.remove("NF")
				j3del.append(el.split("(")[0])#clinicalexome
			if el.split(" ")[0] in synd:
				if "NF" in j4del: j4del.remove("NF")
				j4del.append(el.split("(")[0])#syndromes
			if el.split(" ")[0] in gwas:
				if "NF" in j5del: j5del.remove("NF")
				j5del.append(gwas[el.split("(")[0]])#gwas
			if el.split(" ")[0] in gtex:
				if "NF" in j6del: j6del.remove("NF")
				j6del.append(gtex[el.split("(")[0]])#gtex				
		freqdel, genecard=get_oe(dels)
		cpsdel, phenscoredel=get_phen(dels, clingenn, ddg2p, pannel, chra, bp1+"-"+bp2, "del", dic_gtex, terms)
#		firstA=[genesA,"","& ".join(ddA),freqA, cpsA[1],cpsA[2],phenscoreA,cpsA[3],cpsA[4],cpsA[5],cpsA[7],cpsA[0],cpsA[8], "& ".join(j3A),"& ".join(j4A),"& ".join(j2A),"& ".join(j1A), cpsA[6]]
		dell=[genesdel, "& ".join(dels), freqdel,genecard, cpsdel[1],cpsdel[2],phenscoredel,cpsdel[3],cpsdel[4],cpsdel[5], "& ".join(j5del),"& ".join(j6del), "& ".join(j3del),"& ".join(j4del),"& ".join(j2del),"& ".join(j1del),cpsdel[6]]
	#print("get geneclass", ddA, cpsA, noncodA, ddB, cpsB, noncodB, genesdel, False)
	#cl, afectA, afectB=get_gene_class(ddA, cpsA, noncodA, ddB, cpsB, noncodB, genesdel, False)###calcula em que grau de prioritizacao é que aquilo fica
	cl, afectA, afectB= classes_for_cpx.interpret_to_class("del", classes_for_cpx.get_gene_class(ddA, cpsA, noncodA), classes_for_cpx.get_gene_class(ddB, cpsB, noncodB), "", [chra, int(bp1), int(bp2)])###calcula em que grau de prioritizacao é que aquilo fica
	if firstA[0]!="0":
		#firstA[2]=afectA
		firstA[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chra+'%3A'+bp1+'";"'+afectA+'")'
	if firstB[0]!="0":
		#firstB[2]=afectB
		firstB[3]='=HIPERLIGAÇÃO("https://www.ensembl.org/Homo_sapiens/Location/View?r='+chra+'%3A'+bp2+'";"'+afectB+'")'
	high_ph, pval_ph=get_highest_ph(phenscoreA, phenscoreB, phenscoredel)
	if cl=="NF":
		cl, firstA, secondA, firstB, secondB, dell=conserged_for_cnv(chra, bp1, bp2, Conserved_reg, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondA, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"], secondB, ["0","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF","NF"])
		return GenArch, cl, firstA, secondA, firstB, secondB, dell, high_ph, pval_ph
	else:
		return GenArch, cl, "\t".join(firstA), secondA, "\t".join(firstB), secondB, "\t".join(dell), high_ph, pval_ph
		

def merge_di_tol(dic1, dic2):
	new=[]
	for key,value in dic1.items():
		if key in dic2:
			if dic2[key]!=value:
				new.append(dic2[key])
		new.append(value)
	for key,value in dic2.items():
		if key not in dic1:
			new.append(value)
	return new
	
def calc_ACMG(tgenes, oes_rem, oes_dis, ttype, chrr, start, end):
	final_value=0.0
	if tgenes==0:###########
		final_value-= 0.6
	#os genes affectados sao haploins
	if ttype=="Duplication":
		if len(oes_rem)!=0 or len(oes_dis)!=0:
			final_value+=0.6
	elif any("5'" in s for s in oes_dis):
		final_value+=0.9
	elif any("3'" in s for s in oes_dis):
		final_value+=0.6
	else:
		final_value+=0.45
	dic={"ovl":"full", "version": "hg38", "perc":"70", "tt": ttype, "chrA":chrr, "brA":str(start)+"-"+str(end), "dats":["gnomad", "chaisson", "collins", "DGV", "1000Genomes", "ClinGenben", "ClinGenlben"]}
	gg=cross_deletions_with_dgv_results_Vfor_new.exect_for_ACMG(dic)
	if True in gg and ttype=="Deletion":
		final_value-=1
	elif True in gg and ttype=="Duplication":
		final_value-=0.9
	if tgenes >34 and tgenes<=49:
		final_value+= 0.45
	elif tgenes >49:
		final_value+= 0.9
	if final_value<=-1:
		classs="Benign"
	elif final_value>-1 and final_value<=-0.9:
		classs="Likely Benign"
	elif final_value>-0.9 and final_value<0.9:
		classs="VUS"
	elif final_value>=0.9 and final_value<1:
		classs="Likely Pathogenic"
	else:
		classs="Pathogenic"
	return str(final_value)+"\t"+classs

def select_by_oe(cps, genes):
	retrive={}
	if genes!=["NF"]:
		if cps!="":
			reg=cps[0].split(";")
		if genes!=[""]:
			for gg in genes:
				gg1=gg.split(":")[0]
				el=gg1.replace("_intronic","")
				if "_nd" not in el and "_na" not in el and el!="NF":
					bb=el.split("(")[1][:-1]
					if float(bb)<=0.35:
						if cps!="":
							v=genes.index(gg)
							retrive[el]=reg[v]
						else:
							retrive[el]="del"
	return retrive

def execute_hposim(terms, pheno):
	omm={}
	ommtot={}
	for el in pheno:
		om="OMIM:"+el[0]
		#bb=['Rscript phen_disease_for_tool.R  "HP:0011968" "HP:0001250" "HP:0012443" "HP:0000708" "HP:0000729" "OMIM:617616"']
		bb=['/home/dgrc/R_para_HPO/bin/Rscript phen_gene_for_tool.R '+terms.replace(',',' ')+' '+om]#ver como por a dar
		aa=subprocess.check_output(bb, shell=True, stderr=subprocess.STDOUT)
		#aa=subprocess.call(bb)#, stdout=subprocess.PIPE)
		gg=str(aa)
		ff=gg.split("[1]")[-1]
		gg=ff.split('"')[1]
		if gg.startswith("Sc")==False:
			hh=gg.split(" ")
			omm[el[0]]=float(hh[1])
			ommtot[el[0]]=gg
	return order_it(omm, ommtot)

def order_it(omm, ommtot):
	if len(omm)==0:
		return "NF"
	a=1
	newit={}
	ord=dict(sorted(omm.items(),  key=lambda item: item[1], reverse=True))
	for key,value in ord.items():
		if a<=3:
			newit[key]=ommtot[key]
		else:
			break
		a+=1
	return newit
	
def merge_for_phen(trans1, GG1, other1, is_in1, trans2, GG2, other2, is_in2):
	if is_in1==True or is_in2==True:
		is_in=True
	else:
		is_in=False
	trans=list(set(trans1+trans2))
	GG=""
	if GG1==GG2:
		GG=GG1
	elif GG1!="":
		if GG2!="":
			GG=GG1+","+GG2
		else:
			GG=GG1
	elif GG2!="":
		GG=GG2
	if other1==other2:
		other=other1
	elif other1!="":
		if other2!="":
			other=other1+","+other2
		else:
			other=other1
	elif other2!="":
		other=other2
	return trans, GG, other, is_in


def merge_dirs(disrA, disrB):
	disf={}
	for key,value in disrA.items():
		if key in disrB:
			if "NF"==value[2]:
				if value[-1][1]=="-1":
					disf[key]=[value[0], value[1], "3'UTR-"+disrB[key][2], disrB[key][-1]]
				else:
					disf[key]=[value[0], value[1], "5'UTR-"+disrB[key][2], disrB[key][-1]]
			elif "NF"==disrB[key][2]:
				if value[-1][1]=="-1":
					disf[key]=[value[0], value[1], value[2]+"5'-UTR", value[-1]]
				else:
					disf[key]=[value[0], value[1], value[2]+"3'-UTR", value[-1]]
			else:
				if "Exon" in value[2]:
					disf[key]=[value[0], value[1], value[2]+"-"+disrB[key][2], value[-1]]
				else:
					disf[key]=[value[0], value[1], value[2]+"-"+disrB[key][2], disrB[key][-1]]
		else:
			disf[key]=value
	for key, value in disrB.items():
		if key not in disrA:
			disf[key]=value
	#print("disf", disf)
	return disf
					
def get_phen_ins(dd, clingenn, ddg2p, panel, chra, bp1, ttype, gtex, terms):
	att=["start_position","end_position", "exon_chrom_start", "exon_chrom_end", "strand", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "genomic_coding_start", "genomic_coding_end"]
	att1=["ensembl_transcript_id", "transcript_length", "ensembl_gene_id",  "external_gene_name", "transcript_start", "transcript_end", "transcript_biotype"]
	dda=[]
	cc=[]
	final=[[],[],[],[],[],[],[],[], []]#disr, omimid, fenotipo, omimphen,ddg2p, clnge, painel, gene(transcripts), gene(others)
	noncod=[]
	phn=[]
	for el in dd:
		if el.split(" (")[0] in gtex:
			gt=gtex[el.split(" (")[0]]
		else:
			gt=""
		disr="NF"
		GG="NF"
		other="NF"
		if el!="NF" and el!="nd":
			gg=el.split(" (")[0]
			response=treat_genes_for_table.search({'chromosome_name': [chra], 'start': [bp1.split("-")[0]], 'end':[bp1.split("-")[1]]}, att1, "hg38")
			trans1, GG1, other1, is_in1=treat_genes_for_table.parse_biger_trans(response, el, bp1.split("-")[0], gt)
			trans2, GG2, other2, is_in2=treat_genes_for_table.parse_biger_trans(response, el, bp1.split("-")[1], gt)
			#print("a1", trans1, GG1, other1, is_in1, trans2, GG2, other2, is_in2)
			trans, GG, other, is_in=merge_for_phen(trans1, GG1, other1, is_in1, trans2, GG2, other2, is_in2)
			#print("a2", trans, GG, other, is_in)
			if len(trans)>0 and all('' == s or s.isspace() for s in trans)==False and ttype!="del":
				disrA=treat_genes_for_table.ivsreport(treat_genes_for_table.parse_first_search(treat_genes_for_table.search({"link_ensembl_transcript_stable_id":trans}, att, "hg38")),bp1.split("-")[0], "")###<
				disrB=treat_genes_for_table.ivsreport(treat_genes_for_table.parse_first_search(treat_genes_for_table.search({"link_ensembl_transcript_stable_id":trans}, att, "hg38")),bp1.split("-")[1], "")###<
				disr=merge_dirs(disrA, disrB)
			else:
				trans=[]
			names=[]
			dii=[]
			if gt=="":
				names=["NF"]
				dii=["NF"]
			elif is_in==False:
				names=[gt]
				dii=["NF"]
			for ee in trans:
				names.append(ee)#transcripto disrumpid
				dii.append(disr[ee][-2])#posiçao da disrupção
				noncod.append(disr[ee][-1])#posiçao do noncoding
			bb=""
			if len(GG)>1 and len(other)>1:
				bb=GG+","+other
			elif len(GG)>1:
				bb=GG
			elif len(other)>1:
				bb=other
			final[0].append(", ".join(dii))
			final[7].append(", ".join(names))
			final[8].append(bb)
			omid, pheno=marrvel_input.get_omim(gg)
			phenscore="NF"
			if pheno!=["NOOMIM"]:
				phenscore=execute_hposim(terms, pheno)
				gghg=[]
				if phenscore!="NF":
					for key,value in phenscore.items():
						gghg.append(value)
					phn.append(", ".join(gghg))
				else:
					phn.append(phenscore)
			else:
				phn.append(phenscore)
			omid=omid.replace("NOOMIM", "NF")
			if gg in ddg2p:
				dda=ddg2p[gg]
			if gg in clingenn:
				cc=clingenn[gg]
			final[1].append(omid)
			a,b,c,d=merge_ph(pheno, dda, cc, phenscore)
			final[2].append(a)
			final[3].append(b)
			final[4].append(c)
			final[5].append(d)
			final[6].append(parse_oe(panel, gg))
	#print("final", final, noncod, phn)	
	if final==[[],[],[],[],[],[],[]]:
		final=[["NF"],["NF"],["NF"],["NF"],["NF"],["NF"],["NF"]]
	if ttype!="del":
		return ["& ".join(final[0]),"& ".join(final[1]), "& ".join(final[2]), "& ".join(final[3]), "& ".join(final[4]), "& ".join(final[5]), "& ".join(final[6]), "& ".join(final[7]), "& ".join(final[8])], noncod, "& ".join(phn)
	else:
		return ["& ".join(final[0]),"& ".join(final[1]), "& ".join(final[2]), "& ".join(final[3]), "& ".join(final[4]), "& ".join(final[5]), "& ".join(final[6]), "& ".join(final[7]), "& ".join(final[8])], "& ".join(phn)




def get_phen(dd, clingenn, ddg2p, panel, chra, bp1, ttype, gtex, terms):
	att=["start_position","end_position", "exon_chrom_start", "exon_chrom_end", "strand", "external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "genomic_coding_start", "genomic_coding_end"]
	att1=["ensembl_transcript_id", "transcript_length", "ensembl_gene_id",  "external_gene_name", "transcript_start", "transcript_end", "transcript_biotype"]
	dda=[]
	cc=[]
	final=[[],[],[],[],[],[],[],[], []]#disr, omimid, fenotipo, omimphen,ddg2p, clnge, painel, gene(transcripts), gene(others)
	noncod=[]
	phn=[]
	for el in dd:
		if el.split(" (")[0] in gtex:
			gt=gtex[el.split(" (")[0]]
		else:
			gt=""
		disr="NF"
		GG="NF"
		other="NF"
		if el!="NF" and el!="nd":
			gg=el.split(" (")[0]
			trans, GG, other, is_in=treat_genes_for_table.parse_biger_trans(treat_genes_for_table.search({'chromosome_name': [chra], 'start': [bp1.split("-")[0]], 'end':[bp1.split("-")[1]]}, att1, "hg38"), el, bp1.split("-")[0], gt)
			if len(trans)>0 and all('' == s or s.isspace() for s in trans)==False and ttype!="del":
				disr=treat_genes_for_table.ivsreport(treat_genes_for_table.parse_first_search(treat_genes_for_table.search({"link_ensembl_transcript_stable_id":trans}, att, "hg38")),bp1.split("-")[0], "")###<
			else:
				trans=[]
			names=[]
			dii=[]
			if gt=="":
				names=["NF"]
				dii=["NF"]
			elif is_in==False:
				names=[gt]
				dii=["NF"]
			for ee in trans:
				if ee in disr:
					names.append(ee)#transcripto disrumpid
					dii.append(disr[ee][-2])#posiçao da disrupção
					noncod.append(disr[ee][-1])#posiçao do noncoding
			bb=""
			if len(GG)>1 and len(other)>1:
				bb=GG+","+other
			elif len(GG)>1:
				bb=GG
			elif len(other)>1:
				bb=other
			final[0].append(", ".join(dii))
			final[7].append(", ".join(names))
			final[8].append(bb)
			omid, pheno=marrvel_input.get_omim(gg)
			phenscore="NF"
			if pheno!=["NOOMIM"]:
				phenscore=execute_hposim(terms, pheno)
				gghg=[]
				if phenscore!="NF":
					for key,value in phenscore.items():
						gghg.append(value)
					phn.append(", ".join(gghg))
				else:
					phn.append(phenscore)
			else:
				phn.append(phenscore)
			omid=omid.replace("NOOMIM", "NF")
			if gg in ddg2p:
				dda=ddg2p[gg]
			if gg in clingenn:
				cc=clingenn[gg]
			final[1].append(omid)
			a,b,c,d=merge_ph(pheno, dda, cc, phenscore)
			final[2].append(a)
			final[3].append(b)
			final[4].append(c)
			final[5].append(d)
			final[6].append(parse_oe(panel, gg))
	#print("final", final, noncod, phn)
	if final==[[],[],[],[],[],[],[]]:
		final=[["NF"],["NF"],["NF"],["NF"],["NF"],["NF"],["NF"]]
	if ttype!="del":
		return ["& ".join(final[0]),"& ".join(final[1]), "& ".join(final[2]), "& ".join(final[3]), "& ".join(final[4]), "& ".join(final[5]), "& ".join(final[6]), "& ".join(final[7]), "& ".join(final[8])], noncod, "& ".join(phn)
	else:
		return ["& ".join(final[0]),"& ".join(final[1]), "& ".join(final[2]), "& ".join(final[3]), "& ".join(final[4]), "& ".join(final[5]), "& ".join(final[6]), "& ".join(final[7]), "& ".join(final[8])], "& ".join(phn)



def get_oe(dd):#http://www.genecards.org/cgi-bin/carddisp.pl?gene=ENSG00000188107
	e=10.0
	gk=""
	for el in dd:
		if "(" in el and "Y_RNA" not in el:
			hh=el.replace("(intronic)","")
			ee=hh.split("(")
			gg=ee[1].split(":")[0]
			if "N" not in gg and "n" not in gg:
				if float(gg.split(")")[0])< e:
					e=float(gg.split(")")[0])
					gk=ee[0][:-1]
	if e!=10.0:
		return str(e), '=HIPERLIGAÇÃO("http://www.genecards.org/cgi-bin/carddisp.pl?gene='+gk+'";"'+gk+'")'
	else:
		return "NF", "NF"


def read_acmg(acmg):
	dd=[]
	f=open(acmg)
	for el in f:
		ff=el.split("\t")[2]
		dd.append(ff)
	f.close()
	return dd

def get_Overlap(a,b):
	return max(0,min(int(a[-1]),int(b[-1]))-max(int(a[0]),int(b[0])))



def get_dd2p(infile):
	f=open(infile)
	dd2p={}
	for i in f:
		line=i.split(",")
		if line[0] not in dd2p:
			dd2p[line[0]]=[[],[],[]]
			dd2p[line[0]][2].append(line[4])
			dd2p[line[0]][1].append(line[3].strip('"'))
			dd2p[line[0]][0].append(line[2].strip('"'))	 
	f.close()
	return dd2p
	
	
def clingen(infile):
		f=open(infile)
		dic={}
		for i in f:
				line=i.split("\t")
				if line[0] not in dic:
						dic[line[0]]=[[line[1]], [line[3]], [line[4]]]
				else:
						dic[line[0]][0].append(line[1])
						dic[line[0]][1].append(line[3])
						dic[line[0]][2].append(line[4])
		f.close()
		return dic

def sort_pheno(pheno, phenscore):
	newd=[]
	#print("pheno", pheno)
	#print("phenscore", phenscore)
	for key, value in phenscore.items():
		for el in pheno:
			if el[0]==key:
				newd.append(el)
				break
	for ele in pheno:
		if ele[0] not in phenscore:
			newd.append(ele)
	return newd
		
def merge_ph(pheno, dda, cc, phenscore):
	final=[[],[],[],[]]
	if phenscore!="NF":
		sorted_pheno=sort_pheno(pheno, phenscore)
	else:
		sorted_pheno=pheno
	if pheno!=["NOOMIM"]:
		for el in sorted_pheno:#for el in pheno:
			ph=el[1]
			phmim=el[0]+"_"+el[-1]
			cld="NF"
			ccg="NF"
			if len(dda)>1:
				if el[0] in dda[1]:
					uu=dda[1].index(el[0])
					cld=dda[2][uu]
			if len(cc)>1:
				if el[0] in cc[1]:
					uu=cc[1].index(el[0])
					ccg=cc[2][uu]
			final[0].append(ph)
			final[1].append(phmim)
			final[2].append(cld)
			final[3].append(ccg)
	if len(dda)>1:
		for el in dda[1]:
			if el=="No disease omim":
				uu=dda[1].index(el)
				final[0].append(dda[0][uu])
				final[1].append("NF")
				final[2].append(dda[2][uu])
				final[3].append("NF")
	if len(cc)>1:
		for el in cc[1]:
			if el=="No disease omim":
				uu=cc[1].index(el)
				final[0].append(cc[0][uu])
				final[1].append("NF")
				final[2].append("NF")
				final[3].append(cc[2][uu])
	if final==[[],[],[],[]]:
		final=[["NF"],["NF"],["NF"],["NF"]]
	return ", ".join(final[0]), ", ".join(final[1]), ", ".join(final[2]), ", ".join(final[3])
		


def parse_oe(oe, ell):
	"""reads the oe file or the HI/triplo file or the genehancer file
	, and returns the oe value for the
	gene indicated by the variable ell."""
	f=open(oe)
	f.readline()
	nn="NF"
	for el in f:
		line=el.split("\t")
		if line[0]==ell:
				nn=line[-1].strip()#returns the oe with the confidence interval, in all cases
				break
	f.close()
	return nn

def read_for_del_dup(infile):
	f=open(infile)
	counter=1
	dic={}
	for i in f:
		line=i.split("\t")
		if line[0]!="":
			dic[line[0]]=i.strip("\n")
		else:
			dic["cov_"+str(counter)]=i.strip("\n")
			counter+=1
	f.close()
	return dic

def read_gwas(infile):#gene[snp1,sno2...]  input é o traits to genes
	f=open(infile)
	dic={}
	nd={}
	for i in f:
		line=i.split("\t")
		if line[0] not in dic and line[6].strip()!=0:
			dic[line[0]]=[line[6].strip()+" SNPs - "+line[4]+"["+line[5]+"]"]
		elif line[6].strip()!=0:
			dic[line[0]].append(line[6].strip()+" SNPs - "+line[4]+"["+line[5]+"]")
	f.close()
	for key,value in dic.items():
		nd[key]=",".join(value)
	return nd


def read_gtex(infile):#gene[tissue1, tissue2...] input top3 gtext svinterpreter
	f=open(infile)
	dic={}
	nd={}
	for i in f:
		line=i.split("\t")
		if line[1] not in dic:
			dic[line[1]]=[line[-1].strip()]
		else:
			dic[line[1]].append(line[-1].strip())
	f.close()
	for key,value in dic.items():
		nd[key]=",".join(value)
	return nd


#gwas gtex
name=argv[2]+"tab4"
subprocess.run(["rm", name])
if argv[2]=="del" or argv[2]=="dup":
	do_it(read_for_del_dup(argv[1]), argv[2], argv[3], argv[4],argv[5], argv[6],argv[7], argv[8],argv[9], argv[10], argv[11], argv[12], argv[13], "", "", argv[14], argv[15], read_gwas(argv[16]), read_gtex(argv[17]))
else:
	do_it(read_table(argv[1]), "",argv[3], argv[4],argv[5], argv[6],argv[7], argv[8], argv[9], argv[10], argv[11], argv[12], argv[13], argv[14], argv[15], argv[16], argv[17], read_gwas(argv[18]), read_gtex(argv[19]))

