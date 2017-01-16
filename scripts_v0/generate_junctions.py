criteriaDict= {'receptor':'R1', 'loop':'goodLoop', 'helix_context':'rigid', 'offset':1, 'total_length':10}
variants=variantFun.findVariantNumbers(table, criteriaDict)

variant_subtable = variant_table.loc[variants]
variant_subtable = variant_subtable.sort('dG').dropna(axis=0, subset=['dG'])
variant_subtable = variant_subtable.loc[variant_subtable.loc[:, 'qvalue'] < 0.05]
variant_subtable.loc[:, 'c'] = np.nan

topologies = np.unique(variant_subtable.loc[:, 'topology'])

for i, topology in enumerate(topologies):
    index = variant_subtable.loc[:, 'topology'] == topology
    variant_subtable.loc[index, 'c'] = i

cNorm  = colors.Normalize(vmin=0, vmax=len(np.unique(variant_subtable.loc[:, 'topology']))-1)
scalarMap = cmx.ScalarMappable(cmap='Set1', norm=cNorm)
c =[scalarMap.to_rgba(i) for i in variant_subtable.loc[:, 'c']]

plt.figure()

plt.scatter(np.arange(len(variant_subtable)), variant_subtable.loc[:, 'dG'], c =c, linewidth=0)

# select junctions above
index = np.linspace(0, len(variant_subtable)-1, 150).astype(int)
variant_subtable.iloc[index].loc[:, 'topology']

counts = pd.Series(index = topologies)
for topology in topologies:
    counts.loc[topology] = (variant_subtable.iloc[index].loc[:, 'topology'] == topology).sum()