import os

R = 70
indi_len = 12

def rc(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

oligo_f = 'GTTTAATTGAGTTGTCATATGTTAATAACGGTAT'
oligo_r = rc(oligo_f)
oligo_fl = []
oligo_rl = []
oligo_len = 12

for i in range(len(oligo_f) - oligo_len + 1):
    oligo_fl.append(oligo_f[i:i + oligo_len])
    oligo_rl.append(oligo_r[i:i + oligo_len])


fw = open('oligo_result.txt','w')
fw.write('Gene\tCell\tyear\tdate\tindex\ttarget\tref\tstrand\tall\tnon_oligo\t%\toligo\t%\tDistal_accuracy\t%\tDistal_non_accuracy\t%\tProximal_accuracy\t%\tnon_accuracy\t%\n')

for line in open('oligo_list.txt').readlines():

    line_sp = line.strip().split('\t')
    
    print(line_sp[0] + '\t' + line_sp[1])

    ref = line_sp[6].upper()
    target = line_sp[5].upper()
    target_5seq = ''
    target_3seq = ''
    strand = '+'

    if ref.find(target) != -1:
        st = ref.find(target) + 16 - R
        ed = st + 2 * R
        target_5seq = ref[st + R - 9: st + R + 1]
        target_3seq = ref[st + R + 1: st + R + 11]
    elif ref.find(rc(target)) != -1:
        st = ref.find(rc(target)) + 6 - R
        ed = st + 2 *R
        target_5seq = ref[st + R - 10: st + R]
        target_3seq = ref[st + R: st + R + 10]
        strand = '-'
    else:
        print('ERROR no target in ref')
        print(line)
        input()

    marker = ref[st + R - 5: st + R + 5]

    if st < 0: st = 0
    if ed > len(ref): ed = len(ref)

    ind_f = ind_f = ref[st: st + indi_len]
    ind_r = ref[ed - indi_len: ed]

    ind_fl = []
    ind_rl = []
    
    ref_range = ref[st:ed]

    for i in oligo_fl:
        if ref_range.find(i) != -1:
            print('ERROR referece has part of oligo')
            print(line)
            input()

    for i in oligo_rl:
        if ref_range.find(i) != -1:
            print('ERROR referece has part of oligo')
            print(line)
            input()


    for i in range(indi_len):
        for nt in 'ATGC':
            ind_fl.append(ind_f[:i] + nt + ind_f[i + 1:])
            ind_fl.append(ind_f[:i] + ind_f[i+1:])
            ind_rl.append(ind_r[:i] + nt + ind_r[i + 1:])
            ind_rl.append(ind_r[:i] + ind_r[i + 1:])

    ref_dir = line_sp[2]
    f1 = ref_dir + '/' + line_sp[3]
    f2 = ref_dir + '/' + line_sp[4]
	index = line_sp[3][:line_sp[3].find('.')].replace('_R1','') 
    os.system('fastq-join {0} {1} -o {2}'.format(f1, f2, index + '.fastq'))

    cnt_d = {'all': 0 , 'non': 0, 'oligo': 0, '5ac': 0, '5non': 0, '3ac': 0, '3non': 0}
    seq_d = {}
    
    fa = open('test.txt','w')
    with open(index + '.fastqjoin') as f:
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
                    sed = s.find(i) + indi_len
                    break

            if sst == -1 or sed == -1:
                continue

            s = s[sst: sed]
            if s in seq_d.keys():
                seq_d[s] += 1
            else:
                seq_d[s] = 1
            fa.write(s+'\n')
        fa.close()
			
        for seq, cnt in seq_d.items():
            
            #if cnt <= 1: continue
            cnt_d['all'] += cnt
            
            oligo_st = -1
            oligo_ed = -1
            oligo_n = [-1, -1]

            for on, i in enumerate(oligo_fl):
                if seq.find(i) != -1:
                    if oligo_st == -1:
                        oligo_st = seq.find(i)
                        oligo_n[0] = on
                    oligo_ed = seq.find(i) + len(i) - 1
                    oligo_n[1] = on

            if oligo_st == -1 or oligo_ed == -1:
                oligo_st = -1
                oligo_ed = -1
                oligo_n = [-1, -1]
                for on, i in enumerate(oligo_rl):
                    if seq.find(i) == -1:
                        if oligo_st != -1:
                            oligo_st = seq.find(i)
                            oligo_n[0] = on
                        oligo_ed = seq.find(i) + len(i) - 1
                        oligo_n[1] = on
            
            if oligo_st == -1 or oligo_ed == -1:
                cnt_d['non'] += cnt

            else:
                cnt_d['oligo'] += cnt
                if seq[oligo_st - 10: oligo_st] == target_5seq and oligo_n[0] == 0:
                    cnt_d['5ac'] += cnt
                else:
                    cnt_d['5non'] += cnt
                if seq[oligo_ed + 1: oligo_ed + 11] == target_3seq and oligo_n[1] == len(oligo_fl) - 1:
                    cnt_d['3ac'] += cnt
                else:
                    cnt_d['3non'] += cnt

    fw.write(line.strip() + '\t' + strand + '\t')
    fw.write('{0}\t{1}\t{2}\t{3}\t{4}\t'.format(cnt_d['all'], cnt_d['non'], round(cnt_d['non']*100/cnt_d['all'], 2), cnt_d['oligo'], round(cnt_d['oligo']*100/cnt_d['all'], 2)))

    print('{0}\t{1}\t{2}\t{3}\t{4}\t'.format(cnt_d['all'], cnt_d['non'], round(cnt_d['non']*100/cnt_d['all'], 2), cnt_d['oligo'], round(cnt_d['oligo']*100/cnt_d['all'], 2)))

    if cnt_d['oligo'] == 0:
        fw.write('\n')
    if strand == '+':
        fw.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(cnt_d['5ac'], round(cnt_d['5ac']*100/cnt_d['oligo'], 2), cnt_d['5non'], round(cnt_d['5non']*100/cnt_d['oligo'],2), cnt_d['3ac'], round(cnt_d['3ac']*100/cnt_d['oligo'], 2), cnt_d['3non'], round(cnt_d['3non']*100/cnt_d['oligo'], 2)))
    elif strand == '-':
        fw.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(cnt_d['3ac'], round(cnt_d['3ac']*100/cnt_d['oligo'], 2), cnt_d['3non'], round(cnt_d['3non']*100/cnt_d['oligo'],2), cnt_d['5ac'], round(cnt_d['5ac']*100/cnt_d['oligo'], 2), cnt_d['5non'], round(cnt_d['5non']*100/cnt_d['oligo'], 2)))

               
