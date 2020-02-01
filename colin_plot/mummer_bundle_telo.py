import sys
from interval import Interval
max_bundle_gap  = 20000
diff_dir_fraction = 0.2
diff_dir_length = 100000
min_telo_to_show = 15000

#telo_file = "/Users/alexwang/0data/0mango/genome/telo/hn47.telo.txt"
#mummer_res = "/Users/alexwang/0data/0mango/genome/mummer/gd10_hn47.tab.txt"
mummer_res = sys.argv[1]
telo_file = sys.argv[2]

rec = {}
def total_bundle(test):
    c=0
    for i in test:
        c=c + int(i[-1])
    return c

def judge_dir(test):
    c1,c2=0,0
    for i in test:
        if i[-2] == "1":
            c1 = c1 + int(i[-1])
        else:
            c2=c2 + int(i[-1])

    return "1" if c1 > c2 else "-1"

def merge_list(list1, list2):
    return [list1[0],list2[1],list1[2],list2[3],list1[4], list1[5]+list2[5]]

def linkage(bundle_list):
    """
    remove fake diff direction bundles
    :param bundle_list:
    :return:
    """
    remove_list = []
    num = len(bundle_list)
    for i in range(1,num):
        tmp = bundle_list[i]
        start1, end1, start2, end2,dir,si = tmp[0], tmp[1],tmp[2], tmp[3],tmp[4],tmp[5]
        if dir == bundle_list[i-1][4] and abs(start1 - bundle_list[i-1][1])<=max_bundle_gap and abs(start2 - bundle_list[i-1][3]) <= max_bundle_gap:
            bundle_list[i][0], bundle_list[i][2] = bundle_list[i-1][0], bundle_list[i-1][2]
            bundle_list[i][5] = si + bundle_list[i-1][5]
            remove_list.append(i-1)

    bundle_list = [bundle_list[i] for i in range(num) if i not in remove_list]  # list remove elements by index
    return bundle_list

with open(mummer_res,"r") as f1:
    for idx, line in enumerate(f1):
        if idx > 3:

            line = line.rstrip()
            line = line.split("\t")
            ref,qry = line[13],line[14]
            if (ref,qry) not in rec:
                rec[(ref, qry)] = []

## bundling fragmental alignment block
ll = open(mummer_res).readlines()
for target in rec.keys():
    size = 1
    for i in range(len(ll)):
        line = ll[i].strip().split("\t")
        if i>3 and line[13] == target[0] and line[14] == target[1]:
            s1, e1, s2, e2, dire, bundle_size = int(line[0]), int(line[1]), int(line[2]), int(line[3]), line[12], size
            if rec[target] == []:
                rec[target].append([s1,e1,s2,e2,dire,bundle_size])
            else:
                obj = rec[target][-1]
                if abs(s1 - obj[1]) <= max_bundle_gap and abs(s2 - obj[3]) <= max_bundle_gap and dire == obj[4]:   # 可捆绑的条件
                    obj[1], obj[3] = e1, e2
                    obj[5] += 1
                else:
                    rec[target].append([s1,e1,s2,e2,dire,size])

## drop fake diff direction align
rec2 = {}
for test in rec.items():
    test1 = test[1]
    for idx,j in enumerate(test1):
        if j[-2] != judge_dir(test1):   # 如果bundle方向不是主方向
            ss = int(j[-1]) / total_bundle(test1)
            if ss < diff_dir_fraction and abs(int(j[1]) - int(j[0])) < diff_dir_length: # 如果反向序列是假的
                j[-2] = judge_dir(test1)
                j[2],j[3]=j[3],j[2]
    rec2[test[0]] = linkage(test1)
## group by subjec
rec3 = {}
subject = set([i[0] for i in rec.keys()])
for i in subject:
    rec3[i] = []
for ref in subject:
    for i in rec2:
        if ref == i[0]:
            for j in rec2[i]:
                qq = int(str(i[1]).split("_")[1])
                rec3[ref].append([qq]+j)

## keep best alignment and drop overlap
res = []
for gp in rec3:
    rm_idx = []
    group = sorted(rec3[gp], key=lambda x: x[1])
    for idx, item in enumerate(group):
        if idx>=1:
            bundle_st,bundle_end = item[1],item[2]
            pre_st, pre_end = group[idx-1][1], group[idx-1][2]
            if Interval(bundle_st, bundle_end).overlaps(Interval(pre_st, pre_end)):
                if bundle_end-bundle_st < pre_end-pre_st:
                    group[idx],group[idx-1] = group[idx-1],group[idx]
                rm_idx.append(idx-1)
    group = [group[i] for i in range(len(group)) if i not in rm_idx]
    sub_name = int(str(gp).split("_")[1])
    for aa in group:
        res.append([sub_name] + aa)
a = open(telo_file)
memo = []
for i in a.readlines()[1:]:   # 遍历每一条 query 染色体
    i = i.rstrip().split("\t")
    cho, chr_len, left, right = int(i[0]), int(i[1]), i[2], i[3]
    for record in res:    # 遍历每一个bundle
        if record[1] == cho:
            if record[6] == "1":   # query与ref 同向
                if right != "nd" and abs(int(right) - max(record[4],record[5])) < min_telo_to_show:
                    telo2 = int(right) - max(record[4],record[5]) + record[3]
                else:
                    telo2 = "NA"
                if left != "nd" and abs(int(left) - min(record[4],record[5])) < min_telo_to_show:
                    telo1 = int(left) - min(record[4],record[5]) + record[2]
                else:
                    telo1 = "NA"
            if record[6] == "-1":  # query与ref反向
                if right != "nd" and abs(int(right) - record[4]) < min_telo_to_show:
                    telo1 = record[2] - (int(right) - record[4])
                else:
                    telo1 = "NA"
                if left != "nd" and abs(int(left) - record[5]) < min_telo_to_show:
                    telo2 = record[3] + (record[5] - int(left))
                else:
                    telo2 = "NA"
            record.extend([telo1,telo2])
            memo.append(record)
        else:
            continue
a.close()
memo = sorted(memo, key=lambda x:(x[0], x[1]))
# output
pre = mummer_res.split(".")[0]
with open(pre+".bundle_telo2.txt","w+") as fr:
    fr.writelines("ref\ts1\te1\tali\ts2\te2\tbundle_size\tdirection\ttelo1\ttelo2\n")
    for li in memo:
        kk = list(map(str,li))
        kk[1],kk[2],kk[3] = kk[2], kk[3], kk[1]
        kk = "\t".join(kk) + "\n"
        fr.writelines(kk)


