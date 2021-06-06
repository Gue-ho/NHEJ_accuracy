import os
from needle import *

def rc(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

fw = open('non_oligo_result.txt','w')
for line in open('non_oligo_list.txt').readlines():
	line_sp = line.strip().split('\t')
	if len(line_sp) <4:
		continue
	gene = line_sp[0]
	ngs_dir = line_sp[2]
	target = line_sp[5].upper().replace(' ','')
	ref = line_sp[6].upper()
	n = line_sp[3]

	st = -1
	ed = -1
	cv = -1
	R = 70

	if ref.find(target) != -1:
		st = ref.find(target) + 17 - R
		ed = st + 2 * R
		strand = 1
		dup_ins = target[16]
	elif ref.find(rc(target)) != -1:
		st = ref.find(rc(target)) + 6 - R
		ed = st + 2 * R
		strand = -1 
		dup_ins = rc(target[16])
	else:
		print('ERROR no target in ref')
		print(line) 
		input()

	marker = ref[st + R - 5: st + R + 5]
	cv = R

	if st < 0: 
		cv += st
		st = 0
	if ed > len(ref): 
		ed = len(ref)

	ref_cut = ref[st: ed]

	ind_f = ref[st: st + 12]
	ind_r = ref[ed - 12: ed]

	ind_fl = []
	ind_rl = []

	for i in range(12):
		for nt in 'ATGC':
			ind_fl.append(ind_f[:i] + nt + ind_f[i + 1:])
			ind_rl.append(ind_r[:i] + nt + ind_r[i + 1:])
	
	f1 = ref_dir + '/' + line_sp[3]
	f2 = ref_dir + '/' + line_sp[4]
	os.system('fastq-join {0} {1} -o {2}'.format(f1, f2, './' + lien_sp[3][:line_sp[3].find('.')].replace('_R1','') + '.fastq'))

	cnt = {'all': 0 , 'wt': 0, 'ins': 0, 'del': 0}
	ins_cnt = {'dup': 0, 'non': 0}
	seq_d = {}

	with open('./' + lien_sp[3][:line_sp[3].find('.')].replace('_R1','') + '.fastqjoin') as f:
		for ln, s in enumerate(f):

			if ln % 4 != 1:
				continue

			sst = -1
			sed = -1
			for i in ind_fl:
				if s.find(i) != -1:
					sst = s.find(i)
					break
			for i in ind_rl:
				if s.find(i) != -1:
					sed = s.find(i) + 12
					break

			if sst == -1 or sed == -1:
				continue

			s = s[sst: sed]
			if s in seq_d.keys():
				seq_d[s] += 1
			else:
				seq_d[s] = 1
	
	print(ref_cut[cv -2: cv +2])
	for x, y in seq_d.items():
		if y <= 1: continue
		cnt['all'] += y
		if x.find(marker) != -1:
			cnt['wt'] += y
		elif len(x) == ed - st:
			cnt['wt'] += y
		elif len(x) < ed - st:
			cnt['del'] += y
		elif len(x) > ed - st:
			cnt['ins'] += y
			needle_res =  needle(ref_cut, x, 10, 0.5, 10, 0.5)
			ins_seq = ''
			ins_v = False
			ins_st = 0
			ins_ed = 0
			com_cv = cv
			for ii in range(len(needle_res[0])):
				if needle_res[0][ii] == '-':
					ins_seq += needle_res[2][ii]
					if ins_v == False:
						ins_st = ii
					ins_v = True
				elif ins_v == True:
					if ins_st - 2 < com_cv < ii + 2:
						break
					com_cv += len(ins_seq)
					ins_seq = ''
					ins_v = False
			if ins_v == True and len(ins_seq) == 1:
				if ins_seq == dup_ins:
					ins_cnt['dup'] += y
				else:
					ins_cnt['non'] += y
		else:
			print('???')
	
	fw.write(line.strip() + '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t\t{8}\t{9}\n'.format(cnt['all'], cnt['wt'], round(cnt['wt']*100/cnt['all'],2), cnt['ins'], round(cnt['ins']*100/cnt['all'],2), cnt['del'], round(cnt['del']*100/cnt['all'],2), round((cnt['del'] + cnt['ins'])*100/cnt['all'],2), ins_cnt['dup'] + ins_cnt['non'], ins_cnt['dup']))
	
	print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t\t{8}\t{9}\n'.format(cnt['all'], cnt['wt'], round(cnt['wt']*100/cnt['all'],2), cnt['ins'], round(cnt['ins']*100/cnt['all'],2), cnt['del'], round(cnt['del']*100/cnt['all'],2), round((cnt['del'] + cnt['ins'])*100/cnt['all'],2), ins_cnt['dup'] + ins_cnt['non'], ins_cnt['dup']))

	print('{0}_{1}\tfinish!'.format(line_sp[0], line_sp[1]))

fw.close()

		
		 




