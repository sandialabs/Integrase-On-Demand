from collections import defaultdict
import sys,shutil,re,os,time,subprocess,pickle,json,random
#############################################################################################################
def timesince():
	return time.strftime("%H:%M:%S", time.localtime())
def revcomp(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U':'A',
				'a': 't', 'c': 'g', 'g': 'c', 't': 'a','u':'a'} 
	return ''.join(complement.get(base, 'N') for base in reversed(seq))
def build_db(neighbors):
	tax_based,taxonomies = {},[]
	for s in neighbors:
		if len(tax_based) >= max_islands:
			break
		for seq,spec in island_dict.items():
			if any(s in i for i in spec): #if the current taxa (s in neighbor) uses the attB
				tax_based[seq] = spec
				if s not in taxonomies: 
					taxonomies.append(s) #logging which neighbor taxa actually had an attB
				if len(tax_based) >= max_islands:
					break
	return tax_based,taxonomies
def question_hit(hits): #subroutine categorizes hits as BL/BR or B
	partial = defaultdict(list)
	for hit in hits:
		gapopen, mismatch, length, sstrand, saccver, slen, sstart, send, sseq, qlen, qaccver, qstart, qend, qseq, idseq = hit
		if gapopen <= 1 and mismatch <= 2 and int((length/qlen)*100) >= 85:
			if idseq == 'B':
				partial['full'].append(hit)
			elif idseq == 'L':
				partial['L'].append(hit)
			else:
				partial['R'].append(hit)
	return partial
def estimate_island_lengths(key,l,r): #returns distance between two att sites
	LRside_dict={ 
		True : [abs(l[6]-r[7]),[(l[4],tuple(sorted([l[6],r[7]])))]],
		False : [abs(l[7]-r[6]),[(l[4],tuple(sorted([l[7],r[6]])))]],
		(True,True) : [abs(l[5]-l[6])+r[7],[(l[4],tuple(sorted([l[6],l[5]]))),(r[4],tuple(sorted([0,r[7]])))]],
		(True,False) : [abs(l[5]-l[6])+abs(r[5]-r[6]),[(l[4],tuple(sorted([l[6],l[5]]))),(r[4],tuple(sorted([r[5],r[6]])))]],
		(False,True) : [l[7]+r[7],[(l[4],(0,l[7])),(r[4],(0,r[7]))]],
		(False,False) : [l[7]+abs(r[5]-r[6]),[(l[4],(0,l[7])),(r[4],tuple(sorted([r[5],r[6]])))]]
		}
	return LRside_dict[key]
def cleanup(): #marks hits based on tandem repeats
	perf_hit, key_dict, countm, countp = defaultdict(dict), {}, 0, 0
	for attb in cand:
		cand[attb].append("")
		lir = re.findall(r'[ACTGUN]+|[actgun]+', attb)
		id, n = lir[1], 0 # separate id and start flank count
		while len(id) < 10: #
			n += 1
			id = lir[0][-n::]+lir[1]+lir[2][0:n] #if id block is small <7bp, add flank sequence
		snkmers=sorted(re.findall(r'A+|C+|G+|T+', id,re.I),key= lambda x: len(x), reverse=True) # list of tandems
		if any((len(kmer) > int(0.5*len(id)) for kmer in snkmers)): #
			cand[attb][-1] = f"tandem-{snkmers[0][0]}({len(snkmers[0])})"
			countp += 1
			continue
		key_dict[id] = attb
		for contig in genome:
			if cand[attb][1] == 'plus':
				perf_hit[attb][contig] = (len(id),re.findall(id,genome[contig],re.I))
			else:
				perf_hit[attb][f"rc{contig}"] = (len(id),re.findall(revcomp(id),genome[contig],re.I))
	for id in perf_hit:
		count = 0
		for contig in perf_hit[id]:
			count += len(perf_hit[id][contig][1])
			idlen = perf_hit[id][contig][0]
		if count <= 1:
			cand[id][-1] = f"cand"
		elif count > 1:
			cand[id][-1] = f"ID_block={count}"
			countm += 1
	return 
def dedup(hits,islands): # collapses hits by coordinates, using best blast hit and support
	dupes, log_dupe, dupe_log, countl= defaultdict(list), {}, {}, 0
	for ref_attb in cand:
		dupe_all = sorted([seq for seq in cand if any(i in cand[ref_attb][0] for i in cand[seq][0])], key= lambda item: (hits[item][0],hits[item][1]),reverse=False) #picks best cand by the best blast]x
		if len(dupe_all) >= 1:
			for xeq in dupe_all:
				if not xeq in dupes[dupe_all[0]]:
					dupes[dupe_all[0]].append(xeq)
	for ref_attb in dupes: # ref_attb is key/best hit of dupe seqs
		n,overlap_island = ref_attb,False
		for i in dupes[ref_attb]:
			if not overlap_island:
				for (start, end), island in islands.items():
					if start >= cand[i][0][0] and end <= cand[i][0][-1]:
						cand[i][-1],cand[n][-1] = f"in_{island}","dupe"
						n,overlap_island=i,True
						break
					if start in cand[i][0] or end in cand[i][0]:
						cand[i][-1],cand[n][-1] = f"overlaps_{island}","dupe"
						n,overlap_island=i,True
						break
			if overlap_island:
				cand[i][-1] = f"in_{island}"
				continue
			if cand[n][-1] != "cand" and not cand[n][-1].startswith("ID_block") and not cand[n][-1].startswith("tandem") and cand[i][-1] == "cand": #if current best hit is not a cand, assign n to best cand
				n=i
			elif hits[n][0:2] != hits[i][0:2] and cand[i][-1] == "cand": # if the current item is worse blast hit than n and is listed as a cand
				cand[i][-1] = "dupe"
			elif int(cand[i][-2][0][-2]) > int(cand[n][-2][0][-2]) and cand[i][-1] == "cand": #if the current item has better support and is a cand
				if hits[n][0:2] == hits[i][0:2]: #if the blast hit is the same
					cand[i][-1], cand[n][-1] = "cand", "dupe"
					n = i
				else:
					cand[i][-1] = "dupe" #better support worse blast hit == dupe
			if i!=n and cand[i][-1] == "cand":
				cand[i][-1] = "dupe"
		if int(cand[n][-2][0][-2]) <= 5:
			cand[n][-1] = "low_supp"
		if cand[n][-1] == "cand":
			log_dupe[n] = dupes[ref_attb]
		else:
			dupe_log[n] = cand[n][-1]
		for dupseq in dupes[ref_attb]:
			cand[dupseq][-1]+=f":{cand[n][3]}"
	print(f"{timesince()}: Collapsed {len([dupe for dupe in dupes if len(dupes[dupe])>1])} sequences.")
	return log_dupe
def runmash():
	with subprocess.Popen(f"mash dist -d 0.29 {os.path.join(lib_path,'reps.msh')} {input_genome}".split(), stdout=subprocess.PIPE,universal_newlines=True) as process:
		with subprocess.Popen(f"sort -k3,3n".split(), stdin=process.stdout, stdout=subprocess.PIPE, universal_newlines=True) as sort_process:
			mash_out, stderror = sort_process.communicate()
	if not mash_out.strip():
		print("Exiting: No hits found with mash d=0.29, try running in search mode")
		sys.exit()
	with open(os.path.join(lib_path,'reflist.txt'),'r') as reflist: #converting GCAs to GTDB species
		refdict = {line.split('\t')[0]:line.split('\t')[1] for line in reflist}
	return [refdict[i.split('\t')[0]] for i in mash_out.split('\n') if len(i)>0]
def runblast(qblast,first=True,options=False):
	dbname = os.path.join(dbdir,re.sub(r'\.(fa|fasta|FASTA)$', '', input_genome.split('/')[-1]))
	if first:
		print(f"{timesince()}: makeblastdb -in {input_genome} -dbtype nucl -out {dbname}")
		os.system(f"makeblastdb -in {input_genome} -dbtype nucl -out {dbname}")
	results = defaultdict(list)
	command = ["blastn", "-query", f"{qblast}.fa", "-db", dbname,"-task","blastn","-perc_identity","95","-num_threads",f"{cpus}","-outfmt","6 gapopen mismatch length sstrand saccver slen sstart send sseq qlen qaccver qstart qend qseq"]
	if options:
		command += options
	print(f"{timesince()}: "+' '.join(command))
	result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	output = result.stdout
	for line in output.split('\n'):
		parts = line.split('\t')
		if len(parts) == 14:
			idseq, qacc = parts[10].split('/')
			results[qacc].append(tuple([i if not i.isnumeric() else int(i) for i in parts]+[idseq]))
		else:
			continue
	return results
def filterblast(results):
	filteredseqs, coordsDQ, passed_first = defaultdict(list), defaultdict(list), {}
	for qacc, hits in results.items():
		if len(hits) >= 1:
			partial = question_hit(hits)
			if len(partial['full']) == 1:
				if len(partial['L']) == 0 and len(partial['R']) == 0:
					passed_first[qacc] = partial['full'][0]
					continue
				elif len(partial['L']) == 1 and len(partial['R']) == 1:
					(fs,fe),(ls,le),(rs,ren) = [partial[n][0][6:8] for n in ['full','L','R']]
					if abs(fs-ls) <= 16 and abs(fe-ren) <= 16:
						passed_first[qacc] = partial['full'][0]
						continue
				for n in ['L','R']:
					if len(partial[n]) == 0:
						partial[n] = partial['full']
			if len(partial['full']) == 2:
				if len(partial['L']) != 2 or len(partial['R']) != 2:
					continue
			if (len(partial['L']) >= 1 and len(partial['R']) >=1 ) or len(partial['full']) == 2: #check if any partial hits fall within reasonable distance of full hit(or eachother)
				duplicates, remove = [], False
				isleLen = length_dict.get(qacc)
				if not isleLen or len(isleLen) == 0:
					remove = True
					continue
				for l in partial['L']:
					for r in partial['R']:
						ls, rs = [{"plus":True,"minus":False}.get(i) for i in [l[3],r[3]]]
						if l[4] != r[4]: #if on different contigs estimate length of island orient to contig termini
							islelen, hitcoord = estimate_island_lengths((ls,rs),l,r)
						elif l[4] == r[4] and ls == rs:
							if any([( ls and r[6] < dist and abs(l[5]-l[7]) < dist and abs(r[6]-l[7]) > dist) or ( not ls and l[6] < dist and abs(r[5]-r[7]) < dist and abs(r[6]-l[7]) > dist) for dist in isleLen]): #circljunc
								islelen,hitcoord=estimate_island_lengths((ls,rs),l,r)
							else: #simple
								islelen, hitcoord = estimate_island_lengths(ls,l,r)
						else:
							continue
						diffs = [int((abs(islelen-a)/(islelen+a))*100) for a in length_dict.get(qacc)]
						diff=min(diffs)
						smallest=diffs.index(diff)
						if diff <= 1:
							duplicates.append((diff,smallest,hitcoord,l,r))
				if len(duplicates) >= 1:
					hitcoord = sorted([(i,j,l,r) for h,i,j,l,r in duplicates],key= lambda x: x[0])[0] # best hit by smallest diff
					i, j, l, r = hitcoord
					scaf = '/'.join(n[0] for n in j)
					for n in j:
						if not n[1] in coordsDQ[n[0]]:
							coordsDQ[n[0]].append(n[1])
						filteredseqs[scaf].append((qacc,i,n[0],n[-1],(l,r)))
				if not remove and len(partial['full']) == 1:
					passed_first[qacc] = partial['full'][0]
	return passed_first,filteredseqs
def main(y_length,s_length,cpus):
	def get_island_info(key_seq):
		if all((spec == 'NaN' for spec in island_dict[key_seq])):
			return ["NaN" for i in range(9)]
		return sorted((tuple(island[0:3]+island[-3::]) for spec,isle in island_dict[key_seq].items() for island in isle),key= lambda x:int(x[-2]),reverse=True)
	tax_island_dict = island_dict
	if tax or attb_search: # make fresh list of query attBs
		if not attb_search:
			species = runmash()
			tax_island_dict,taxonomies = build_db(species) #collect islands from neighbors
			if len(tax_island_dict) < 50: #if db of species is small then try genera
				genera = [i.split('_')[0] for i in species] #get GTDB genera for best hits
				tax_island_dict,taxonomies = build_db(genera)
			if len(tax_island_dict) < 50:
				if len(tax_island_dict) == 0:
					print(f"{timesince()}: Exiting! No islands found in database for the following genera of the top species: {', '.join(species)}.Run without -tax flag.")
					sys.exit()
				print(f"{timesince()}: Warning! Only {len(tax_island_dict)} reference islands were found. For more reference islands, try running the program without the -tax flag. Starting BLAST search...")
	fresh=True
	if search and not attb_search and y_length==10 and s_length==16:
		qblast = os.path.join(lib_path,"attbs")
		fresh =  not os.path.exists(qblast)
	else:
		qblast = os.path.join(dbdir,"attbs")
	if fresh:
		with open(f"{qblast}.fa",'w') as fa, open(f"{qblast}.isleID.fa",'w') as iid:
			for seq in tax_island_dict:
				integrase_type = list(tax_island_dict[seq].items())[0][1][0][0]
				lir = re.findall(r'[AGTC]+|[acgt]+',seq)
				if len(lir) == 3:
					l, id, r = lir
				else:
					l, id, r= seq[0:10].upper(), seq[10:-10].lower(), seq[-10::].upper()
				if y_length != 10 or s_length != 16:
					flank_length={"Y":y_length,"S":s_length}.get(integrase_type[0].upper())
					if flank_length:
						l, r = l[-flank_length::],r[0:flank_length]
				isleids = [inner_array[1] for inner_dict in island_dict[seq].values() for inner_array in inner_dict]
				def_dict = {"B":f"{l}{id}{r}","L":f"{l}{id}","R":f"{id}{r}"}
				for i in def_dict:
					fa.write(f">{i}/{seq}\n{def_dict[i]}\n")
					iid.write(f">{i}/{','.join(isleids)}\n{def_dict[i]}\n")
	count=1
	for seq,species in tax_island_dict.items():
		for spec,islands in species.items(): # parse the islands of this attB to compute island lengths
			for island in islands:
				coords = [int(match) for match in re.findall(r'[/-]([0-9]+)',island[2])]
				strand,l,r=island[3].split(';')
				cross=False
				if len(coords) == 4:
					cross=True
				sizes,halves=[],[l,r]
				for half in range(2):
					cval=half*2 if cross else 0
					if min([coords[0+cval],coords[1+cval]])==1:
						sizes.append(max([int(i) for i in halves[half].split('-')])-1)
					else:
						sizes.append(max([coords[0+cval],coords[1+cval]])-min([int(i) for i in halves[half].split('-')]))
				length_dict[seq].append(sum(sizes))
	blast_out = runblast(qblast)
	hits, filteredseqs = filterblast(blast_out)
	occdupes, specialhits = defaultdict(list), defaultdict(list)
	for u in filteredseqs:
		for n in filteredseqs[u]:
			dupe = tuple(sorted([i for i in filteredseqs[u] if len(i[-1]) == 2 and min(i[-2]) <= max(n[-2]) and max(i[-2]) >= min(n[-2])],key= lambda x: (x[-1][0][0:2],x[-1][1][0:2])))
			if not dupe in occdupes[u]:
				occdupes[u].append(dupe)
	for ser,reps in occdupes.items():
		rep, others = [i[0] for i in reps],{i[0][0]:[] for i in reps}
		for i in reps:
			others[i[0][0]] += [p[0] for p in i[1::]]
		front, scoord, count = [], [], 0
		for keyseq,reisland,contig,isleCoords,blasts in rep: #for each island in the scaffold
			count += 1
			front.append(f"{contig}/{'-'.join(str(i) for i in isleCoords)}")
			scoord += [(i[-1],f"{i[6]}-{i[7]}") for i in blasts]
			if '/' in ser and count%2 != 0:
				continue
			refseqs = get_island_info(keyseq)[reisland]
			info = ('+'.join((val for val in front)),','.join(f"att{i}:{j}" for i,j in sorted(list(set(scoord)))),','.join([keyseq]+list(refseqs[0:3]+refseqs[-3::])))
			if info not in occupied:
				occupied.append(info)
			specialhits[isleCoords] = '+'.join((val for val in front))
			front, scoord, refseqs = [],[],{}
	passing_special={key: specialhits[key] for key in specialhits if any((item[0]==specialhits[key] for item in occupied))}
	for key,val in hits.items():
		gapopen, mismatch, length, sstrand, saccver, slen, sstart, send, sseq, qlen, qaccver, qstart, qend, qseq, idseq = val
		srange = list(range(sstart-50 if sstart>50 else 0,send+50)) if sstrand == "plus" else list(range(send-50 if send > 50 else 0,sstart+50))
		scoord = f"{sstart}-{send}"
		ref_islands = get_island_info(key)
		cand[key] = [srange,saccver,sstrand,scoord,qseq,ref_islands]
	cleanup()
	for seq in cand:		#collect ints to print
		integrases.update({f"{i[1].split('.')[0]}_{i[0].split('_')[-1]}":i[0] for i in cand[seq][5]})
	if len(cand) == 1:
		for seq in cand:
			cand[seq][-1] += f":{cand[seq][3]}"
			single_cand = cand[seq][1:5]+['\t'.join((cand[seq][5][0]))]
			cand[seq] = [cand[seq][-1]]+cand[seq][1:5]+['\t'.join([','.join(i) for i in cand[seq][5]])]
		return tax_island_dict, [single_cand]
	elif len(cand) == 0:
		return tax_island_dict, None
	else:
		deduped = dedup(hits,passing_special)
	for seq in deduped:
		if len(deduped[seq]) > 1:
			intisles =';'.join([','.join(i[0:2]) for att in deduped[seq] for i in cand[att][-2]])
			line = tuple(cand[seq][1:5]+['\t'.join((cand[seq][5][0]))]+[intisles])
		else:
			line = tuple(cand[seq][1:5]+['\t'.join((cand[seq][5][0]))])
		filteredcand.append(line)
	for seq in cand:
		cand[seq] = [cand[seq][-1]]+cand[seq][1:5]+['\t'.join([','.join(i) for i in cand[seq][5]])]
	return tax_island_dict, filteredcand
def writeout(attbs,final_cands):
	print(f"{timesince()}\tWriting files...")
	def mainout(outfilename,outlist,header=False):
		with open(os.path.join(output_path,outfilename),'w') as f:
			f.write(f"query_contig\tstrand\tquery_coord\tquery_attB\tintegrase\tisleID\tref_scaffold/coords\tsource\tsupport\tisle_type\tints,isleIDs\n")
			if outlist:
				f.writelines('\t'.join(list(map(str, line)))+'\n' for line in sorted(outlist,key= lambda x: int(x[2].split('-')[0])))
			print(f"info:\nmode:{"tax" if tax else "search"};input_genome={input_genome};outputpath={output_path};flank_length(tyrosine)={y_length};flank_length(serine)={s_length};reference_file={ref_path}")
			print(f"summary:\nattBs_in_query={len(attbs)};hit_count={0 if not outlist else len(outlist)};trna_hits={0 if not outlist else sum(1 if re.match(r'[0-9]+\.[0-9]+\.[A-Z]$',i[4].split('\t')[1])  else 0 for i in outlist)};occupied={len(occupied)}")
		return
	jason_derulo = {key:value for key,value in island_dict.items() if key in cand}
	attBs_json = json.dumps(jason_derulo, indent=4)
	with open(os.path.join(output_path,'isles.json'),'w') as f:
		f.write(attBs_json)
	with open(os.path.join(output_path,'attb_dupes.tsv'),'w') as f:
		cand_sorted = sorted(cand.items(),key= lambda x: (int(x[1][3].split('-')[0]),x[1][0].split(':')[1]))
		for attb, vals in cand_sorted:
			f.write(f"{'\t'.join((str(i) for i in cand[attb]))}\n") 
	mainout('final_candidates.tsv',final_cands)
	with open(os.path.join(output_path,'occupied.tsv'),'w') as f:
		f.write("island_coords\tattL_attR\tisleID,ref_scaffold,integrase\n")
		occupied.sort(key= lambda x: (x[0].split('/')[0],int(x[0].split('/')[1].split('-')[0])))
		for seq in occupied:
			f.write('\t'.join([str(n) for n in seq])+"\n") #+"\t"+'\t'.join(seq[2])+"\n")
	print(f"{timesince()}\tDone!")
	with open(os.path.join(lib_path,"ints.gff"),'r') as f, open(os.path.join(output_path,"ints.gff"),'w') as o:
		for line in f:
			attribute = line.split('\t')[-1].split(';')[0].split('=')[1]
			if attribute in integrases:
				o.write(line)
	return
if __name__=='__main__':
	inputs = [int(x) if x.isnumeric() else { "True": True, "False": False , "None": None}.get(x,x) for x in sys.argv[1::]]
	input_genome, output_path, lib_path, y_length, s_length, cpus, tax,search, attb_search, max_islands=inputs
	labels = "Genome FASTA,Output Folder,Reference File,Flank Length(Tyrosine),Flank Length(Serine),CPUs,tax,search,query attB,Minimum number of islands (tax-mode only)".split(',')
	param = dict(zip(labels,inputs))
	for i in param:
		if search and i == "Minimum number of islands (tax-mode only)":
			continue
		print(i+':',param[i] )
	ref_path = os.path.join(lib_path,"isles.pkl")
	dbdir = os.path.join(output_path,"blastDB")
	os.makedirs(dbdir,exist_ok=True)
	with open(ref_path,'rb') as f:
		island_dict = pickle.load(f)
	with open(input_genome,'r') as f:
		genome, seq, contig = {},[],None
		for line in f:
			if line.startswith('>'):
				if contig:
					genome[contig] = ''.join(seq)
				seq,contig = [],line.strip()
			else:
				seq.append(line.strip())
		genome[contig] = ''.join(seq)
	if search and attb_search:
		if os.path.exists(attb_search):
			with open(attb_search,'r') as f:
				firstline = f.readline().strip('\n')
				if any(char.upper() not in ["A","G","C","T"] for char in firstline):
					try:
						search_attbs, seqs, defline = {}, [], None
						for line in [firstline]+f.readlines():
							if line.startswith('>'):
								if defline is not None:
									search_attbs[defline] = ''.join(seqs)
									seqs = []
								defline = line.strip().strip('>')
							else:
								seqs.append(line.strip())
						search_attbs[defline] = ''.join(seqs)
					except Exception as e:
						print(f"Could not unpack {attb_search}. Non-nucleotide characters were found and file is not in FASTA format. Exiting due to the following: {e}.")
						sys.exit()
				else:
					search_attbs = [firstline]+ [line.strip('\n') for line in f.readlines()]
		else:
			search_attbs = attb_search.split(',')
		pattern = "|".join((key for key in search_attbs if key not in island_dict))
		temp = {attb:island_dict.get(attb) for attb in search_attbs if island_dict.get(attb)}
		island_attbs = defaultdict(list)
		if len(pattern) > 0:
			pattern = re.compile(pattern,re.I)
			for key in island_dict:
				hit = re.search(pattern,key)
				if hit:
					hit = hit.group()
					island_attbs[key].append(hit)
		island_dict = {attb:island_dict.get(attb,{"":[[""]*9]}) for attb in island_attbs}|temp
	hits, length_dict, idblocks, refkeys, counthits, cand, integrases, occupied, filteredcand, attbs = defaultdict(list), defaultdict(list), {}, {}, {}, {}, {}, [], [], []
	attbs, final_cands = main(y_length,s_length,cpus)
	writeout(attbs,final_cands)
