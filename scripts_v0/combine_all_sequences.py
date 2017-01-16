
all_files = pd.read_table('150311_library_v2/all_files.txt', header=0)
allFiles = {}
all_exist = True
for expt in np.unique(all_files.expt):
    allFiles[expt] = {}
    folders = np.unique(all_files.loc[all_files.expt==expt].folder)
    for folder in folders:
        sub_files = all_files.loc[all_files.folder==folder]
        filenames = [os.path.join(all_files.loc[i, 'directory'],all_files.loc[i, 'folder'], all_files.loc[i, 'file']) for i in sub_files.index]
        allFiles[expt][folder] = filenames
        for filename in filenames:
            all_exist = all_exist and os.path.exists(filename + '.txt')
            if not os.path.exists(filename + '.txt'):
                print '%d\t%s\t%s does not exist'%(expt, folder, filename)

allSeqs = {}
for expt in allFiles.keys():
    allSeqs[expt] = {}
    
    for folder, filenames in allFiles[expt].items():
        print 'loading sublibrary %s...'%folder
        allSeqs[expt][folder] = hjh.junction_seqs.loadAllseqs(filenames)
    allSeqs[expt] = pd.concat(allSeqs[expt], names=['sublibrary', 'variant'])

for expt in allFiles.keys():
    for folder, filenames in allFiles[expt].items():
        print '%d\t%s\t%d'%(expt, folder, len(allSeqs.loc[(expt, folder)]))

allSeqs = pd.concat(allSeqs, names=['expt'])
allSeqs.to_csv('150311_library_v2/all_10expts.txt', sep='\t', header=True, index=True)

# drop duplicates
allSeqs.drop_duplicates(subset=['tecto_sequence'], inplace=True)
allSeqs.loc[:, 'length_seq'] = [len(s) for s in allSeqs.loc[:, 'sequence']]
max_length = allSeqs.loc[:, 'length_seq'].max()

allSeqs.loc[:, 'tailed_sequence'] = allSeqs.loc[:, 'sequence']
for length in np.unique(allSeqs.length_seq):
    print length
    index = allSeqs.length_seq == length
    allSeqs.loc[index, 'tailed_sequence'] += ''.join(['A']*(max_length - length))

# save for custom array
seqs_only = pd.concat([allSeqs.loc[:, ['tailed_sequence']]]*2)
seqs_only.loc[:, 'ind'] = np.arange(len(seqs_only))
seqs_only.to_csv('150311_library_v2/all_10expts.sequences_id.txt', sep='\t', header=False, index=False)

# plot histogram of lengths

expts = np.sort(allFiles.keys())
bins = np.arange(80, 120)
plt.figure()
plt.hist([len(s) for s in allSeqs.loc[:, 'sequence']], bins=bins, alpha=0.75, histtype='stepfilled', color=sns.xkcd_rgb['ocean green'])



fig = plt.figure(figsize=(3.5, 8))
gs = gridspec.GridSpec(len(expts), 1, hspace=0)
for i, expt in enumerate(expts):
    ax = fig.add_subplot(gs[i, 0])
    ax.hist([len(s) for s in allSeqs.loc[expt, 'sequence']], bins=bins, alpha=0.75, histtype='stepfilled', color=sns.xkcd_rgb['ocean green'])
    if i == 0:
        maxval = 4000
    else:
        maxval = 2000
    ax.set_ylim(0, maxval)
    if i != len(expts)-1:
        ax.set_xticks([])