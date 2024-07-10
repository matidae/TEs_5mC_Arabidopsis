import sys
import glob
from matplotlib import pyplot as plt, cbook as cb
import numpy as np

def counts_barplot(data, pct, category_labels, filename):
    colors = ['#929ECA', '#C2C4CA', '#FC8D62']
#    colors = ['', 'darksalmon', 'crimson']
    legend_labels = ['>80%', '>20% & <=80%', '<=20%']
    title = filename.replace("/fam", "").split("/")[2].replace("counts_","").replace("fam_","").replace("."," ")
    title = title.replace("intersectALL", "all")\
            .replace("intersect2K5p",r"intersecting with genes $\pm$ 2K 5'") + " methylation"
    ind = list(range(len(data[0])))
    #plt.style.use('bmh')
    plt.figure(figsize=(6, 6))
    bottom = np.zeros(len(data[0]))
    len_pct= len(pct)
    for i, row_data in enumerate(data):
        p = plt.bar(ind, row_data, label=legend_labels[i], color=colors[i], bottom=bottom)
        bottom += row_data
        start_pct=(int(len_pct/3)*i)
        end_pct=(int(len_pct/3)*(i+1))
        lab = [str(x).replace(".0","")+" - "+str(y) + " %" for x, y in zip(p.datavalues, pct[start_pct:end_pct])]
        plt.bar_label(p, lab, label_type='center', fontsize=7)
    plt.xticks(ind, category_labels, rotation=45, fontsize=8)
    plt.ylabel("Counts of TEs")
    plt.legend(title="Methylation")
    plt.title(title, fontsize=10)
    plt.tight_layout()
    filename = filename.replace("counts", "barplot") + ".svg"
    plt.savefig(filename, dpi=300)
    print ("Barplot saved to : " + filename)

def process_barplot(counts):
    methyl80 = []
    methyl = []
    non_methyl = []
    cat = []
    pct_methyl80 = []
    pct_methyl = []
    pct_non_methyl = []
    for i in counts:
        e = i.strip().split(" ")
        methyl80.append(int(e[0]))
        methyl.append(int(e[1]))
        non_methyl.append(int(e[2]))
        cat.append(e[3])
        total = int(e[0]) + int(e[1]) + int(e[2])
        pct_methyl80.append('{:.1f}'.format(int(e[0])/total*100))
        pct_methyl.append('{:.1f}'.format(int(e[1])/total*100))
        pct_non_methyl.append('{:.1f}'.format(int(e[2])/total*100))
    data = [methyl80, methyl, non_methyl]
    pct = pct_methyl80 + pct_methyl + pct_non_methyl
    filename = counts.name
    counts_barplot(data, pct, cat, filename)

def boxplot(arrays, labels, filename):
    plt.style.use('bmh')
    fig, ax = plt.subplots(figsize=(6, 6))
    meanpointprops = dict(marker='D', markeredgecolor='grey', markerfacecolor='lightgrey')
    bp = ax.boxplot(arrays, patch_artist=True, meanprops=meanpointprops, meanline=False, 
            showmeans=True, flierprops={'markersize': 5, 'alpha':0.3}) 
    plt.setp(bp['boxes'], color="#1B9E77", alpha=0.8)
    plt.setp(bp['medians'], color="black")
    plt.legend([bp['medians'][0], bp['means'][0]], ['median', 'mean'], fontsize = 8)
    stats = cb.boxplot_stats(arrays)
    means = []
    medians = []
    for i in range(len(stats)):
        means.append("{:.1f}".format(cb.boxplot_stats(arrays)[i]['mean']))
        medians.append("{:.1f}".format(cb.boxplot_stats(arrays)[i]['med']))
    max_val=(max([max(v) for v in arrays]))
    median_v = max_val/100
    mean_v = max_val/35
    for i,v in enumerate(medians):
        plt.text((float(i)+.75), (float(v)+median_v), str(v), fontsize = 7, weight= 'bold')
    for i,v in enumerate(means):
        plt.text((float(i)+1.05), (float(v)-mean_v), str(v), fontsize = 7, style = 'italic')
    title = filename.replace("/fam", "").split("/")[2].replace("values_","").replace("fam_","").replace("."," ")
    title = title.replace("intersectALL", "all")\
            .replace("intersect2K5p",r"intersecting with genes $\pm$ 2K 5'") + " methylation"
    plt.title(title, fontsize=10) 
    ax.set_xticklabels(labels, rotation=45, fontsize=8)
    plt.ylabel('average % methylation per TE', fontsize=9)
    plt.tight_layout()
    plt.rcParams.update(plt.rcParamsDefault)
    filename = filename.replace("values", "boxplot") + ".svg"
    plt.savefig(filename, dpi=300)
    plt.close()
    print ("Boxplot saved to : " + filename)
    
def build_array(values_fh):
    c = 0
    labels = []
    arrays = []
    for i in values_fh:
        if c == 0:
            labels = i.strip().split()
            arrays = [[] for i in range(len(labels))]
        else:
            e = i.rstrip().split()
            for j in range(len(e)):
                arrays[j].append(float(e[j]))
        c+=1
    return(arrays, labels)

def histogram(arrays, labels, filename):
    plt.style.use('bmh')
    fig, ax = plt.subplots(figsize=(6, 6))
    len_labels = 9 if len(labels) > 9 else len(labels)
    for i in range(len_labels):
        plt.hist(arrays[i], label=labels[i], alpha=0.6, bins=50)
    plt.legend(loc='upper center', fontsize = 8)
    title = filename.replace("/fam","").split("/")[2].replace("values_","").replace("fam_","").replace("."," ")
    title = title.replace("intersectALL", "all")\
            .replace("intersect2K5p",r"intersecting with genes $\pm$ 2K 5'") + " methylation"
    plt.title(title, fontsize=10)
    plt.ylabel('Counts', fontsize=10)
    plt.xlabel('% methylation per TE', fontsize=10)
    plt.tight_layout() 
    filename = filename.replace("values", "histogram") + ".svg"
    plt.savefig(filename, dpi=300)
    plt.close()
    filename2 = ""
    if len(labels) > 9:
        fig, ax = plt.subplots(figsize=(6, 6))
        for i in range(5, len(labels)):
               plt.hist(arrays[i], label=labels[i], alpha=0.6, bins=50)
        plt.legend(loc='upper center', fontsize = 8)
        plt.title(title, fontsize=10)
        plt.ylabel('Counts', fontsize=10)
        plt.xlabel('% methylation per TE', fontsize=10)
        plt.tight_layout() 
        filename2 = filename.replace(".svg", "_2.svg")
        plt.savefig(filename2, dpi=300)
        plt.close()
    print ("Histogram saved to : " + filename)
    if len(labels) > 5:
        print ("Histogram saved to : " + filename2)

prefix = sys.argv[1]
directory = sys.argv[2]
counts = glob.glob("./" + directory + "/" + prefix + "*.counts_*")
values = glob.glob("./" + directory + "/" + prefix + "*.values_*")
print ("## Start making plots")
for c in counts:
    if ".svg" not in c:
        counts_fh = open(c, "r")
        process_barplot(counts_fh)
        print("---")
for c in values:
    if ".svg" not in c:
        values_fh = open(c.replace("counts", "values"), "r")
        arrays, labels = build_array(values_fh)
        values_name = values_fh.name
        boxplot(arrays, labels, values_name)
        histogram(arrays, labels, values_name)
        print("---")
