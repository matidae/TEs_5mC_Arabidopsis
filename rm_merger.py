import sys

data = open(sys.argv[1]).readlines()
lengths = {i.split()[0]:int(i.split()[1]) for i in open("RMRBMeta.lengths", "r")}

start_aux = end_aux = 0
contig_aux = ""
idn_aux = ""
strand_aux = ""
name_lst = []

l_total = len(data)
def isLTR():
    internal = list(set(name_lst[1:-1]))
    internal_name = ""
    internal_name_fix = ""
    for n in list(set(name_lst[1:-1])):
        if "-I" in n or "_I" in n:
            internal_name = n
            internal_name_fix = n.replace("-I","").replace("_I","")
            break
    if internal_name != "":    
        length_LTR = lengths[name_lst[0]] * 2
        length_I = lengths[internal_name]
        length_total = length_LTR + length_I
        prop = "{:.3f}".format((end_aux - start_aux + 1)/length_total)
        print(contig_pre, start_aux, end_aux, strand_aux, ".", ".", end_aux-start_aux+1, length_total, \
                prop, idn_pre, order_pre, internal_name_fix, "@mergedLTR", sep="\t")

def isnotLTR():
    length = lengths[name_pre]
    prop = "{:.3f}".format((end_aux - start_aux + 1)/length)
    print(contig_pre, start_aux, end_aux, strand_aux,".",".", end_aux-start_aux+1, \
            length, prop, idn_pre, order_pre, name_pre, "@merged", sep="\t")


for i in range(2, l_total):
    order_pre = data[i-1].split()[8]
    contig_pre = data[i-1].split()[0]
    strand_pre = data[i-1].split()[3]
    name_now = data[i].split()[9]
    name_pre = data[i-1].split()[9]
    start_pre = int(data[i-1].split()[1])
    end_pre = int(data[i-1].split()[2])
    start_now = int(data[i].split()[1])
    end_now = int(data[i].split()[2])
    idn_now = data[i].split()[7]
    idn_pre = data[i-1].split()[7]
    #Condition that ids must be the same
    if idn_now == idn_pre:        
        if start_aux == 0:
            start_aux = start_pre
            end_aux = end_now
            idn_aux = idn_pre
            contig_aux = contig_pre
            strand_aux = strand_pre
            order_aux = order_pre
            name_aux = name_pre
            name_lst.append(name_pre)
            name_lst.append(name_now)
        else:
            end_aux = end_now
            name_lst.append(name_now)
    else:
        if start_aux == 0:
            #If there is nothing to merge
            #continue
            if not "LTR" in order_pre.split("/")[0]:
                length = lengths[name_pre]
                prop = "{:.3f}".format((end_pre - start_pre + 1)/length)
                pre_line = "\t".join(data[i-1].strip().split()[:7])
                print(pre_line, length, prop, idn_pre, order_pre, name_pre, sep="\t")
        else:
            if "LTR" in order_pre.split("/")[0] and len(name_lst)>=3 and "LTR" in name_lst[0] and "LTR" in name_lst[-1]:
                isLTR()
            else:
                if not "LTR" in order_pre.split("/")[0]:
                    isnotLTR()
            start_aux = end_aux = 0
            contig_aux = ""
            strand_aux = ""
            idn_lst = []
            name_lst = []
    #Condition only for the last element of the file
    if l_total == i+1:
        if start_aux == 0:
            if not "LTR" in order_pre.split("/")[0]:
                length = lengths[name_pre]
                prop = "{:.3f}".format((end_pre - start_pre + 1)/length)
                pre_line = "\t".join(data[i-1].strip().split()[:7])
                print(pre_line, length, prop, idn_pre, order_pre, name_pre, sep="\t")
        else:
            if "LTR" in order_pre.split("/")[0] and len(name_lst)>=3 and "LTR" in name_lst[0] and "LTR" in name_lst[-1]:
                isLTR()
            else:
                isnotLTR()
