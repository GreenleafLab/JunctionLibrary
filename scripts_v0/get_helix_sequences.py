import sys
import math
import random
import os

def calc_end_distro_score(selected_constructs):
    bps = ['GC','CG','AU','UA']

    step_distro = {}
    for bp1 in bps:
        for bp2 in bps:
            step_distro[bp1+"="+bp2] = 0

    for c in selected_constructs:
        for step in c.steps:
            if step not in step_distro:
                continue
            step_distro[step] += 1

    avg = 0
    std = 0
    for v in step_distro.itervalues():
        avg += v
    avg /= len(step_distro)

    for v in step_distro.itervalues():
        std += (v - avg) ** 2
    std /= len(step_distro)
    std = math.sqrt(std)

    return std


def calc_bin_score(selected_constructs):
    bins = [0 for x in range(51)]
    avg = 0
    for c in selected_constructs:
        bins[c.dG_bin] += 1

    std = 0
    avg = 0
    for b in bins:
        avg += b
    avg /= len(bins)
    for b in bins:
        std += (b - avg) ** 2
    std /= len(bins)
    std = math.sqrt(std)
    return std


class Construct(object):
    def __init__(self, seq, dG):
        self.seq = seq
        self.dG = dG
        self.steps = []
        self.dG_bin = 0
        seq1 = seq[10:20]
        seq2 = seq[26:36][::-1]
        bps = []
        for i in range(len(seq1)):
            bps.append(seq1[i]+seq2[i])
        steps = []
        for i in range(1, len(seq1)):
            steps.append(bps[i-1]+"="+bps[i])
        self.steps = steps

f = open(os.path.join(os.path.expanduser('~/JunctionLibrary/seq_params/'), 'exhustive_helices.results'))
lines = f.readlines()
f.close()

constructs = []
selected_constructs = []

for l in lines:
    spl = l.split()
    c = Construct(spl[0], float(spl[-1]))
    constructs.append(c)

min_dG = constructs[-1].dG
max_dG = constructs[0].dG

r = max_dG - min_dG
interval = abs(r) / 50

for c in constructs:
    c.dG_bin = (int)((min_dG - c.dG) / interval)

if len(sys.argv[1]) < 2:
    print "need to set a number for the number of constructs you need"
    sys.exit()

selected_num = int(sys.argv[1])
selected_interval = len(constructs) / selected_num
current = 0
while current < len(constructs):
    selected_constructs.append(constructs[current])
    current += selected_interval

distro_score = calc_end_distro_score(selected_constructs)
bin_score = calc_bin_score(selected_constructs)
score = distro_score + 5*bin_score
new_score = 1000
new_selected_constructs = []

for i in range (100000):

    new_selected_constructs = selected_constructs[:]
    rand_pos = random.randint(0, len(new_selected_constructs)-1)
    new_c = random.choice(constructs)
    while new_c in new_selected_constructs:
        new_c = random.choice(constructs)
    new_selected_constructs[rand_pos] = new_c
    distro_score = calc_end_distro_score(new_selected_constructs)
    bin_score = calc_bin_score(new_selected_constructs)
    new_score = distro_score + 5*bin_score
    if new_score < score:
        selected_constructs = new_selected_constructs
        score = new_score
        print "SCORE", score

print "outputing to \'helix_constructs.txt\'"
f = open("helix_constructs.txt", "w")
for c in selected_constructs:
    f.write(c.seq + "\n")
f.close()


