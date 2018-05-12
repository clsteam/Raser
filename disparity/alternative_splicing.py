#!/usr/bin/env python
# -*- coding: UTF-8 -*-
__author__ = 'clsteam'

from bag.configparser import *
from bag.Cluster import *
import pandas as pd
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_go(path,name):
    event_list=list_event(path,name)
    event_dict=stat_event(event_list)
    save_frama(event_dict,path,name)
    piechart(event_dict,path,name)

def list_event(path,name):
    f=open(str(path+'/'+name),'r')
    event_list=list()
    for line in f.readlines():
        line=line.strip('\n')
        s=line.split('	')[1]             #get the second col
        if s[-1]!='F':                    #some event includ 'on' and 'off',which are a couple. We need only one of them.
            if s[0]=='X':                 #a bug in the ASprosile which can bring some  redundancy 'N'
                s=s[1:]
            if s[-1]=='N':
                s=s[:-3]
            event_list.append(s)
    event_list=event_list[1:]
    f.close()
    return event_list

def stat_event(event_list):               #the function can statistics the frequency of every date in the list
    l1 = list(set(event_list))            #seek the element of the list which are different
    l2=[0]*len(l1)
    event_dict=dict(zip(l1,l2))           #creat dictionary
    for ele in event_list:
        event_dict[ele]+=1
    return event_dict

def save_frama(event_dict,path,name):     #save the dict in a regular and a elegent way
    file_object = open(path+'/'+name+'_stat','w')
    file_object.writelines('event_type'+'    '+'value'+'\n')
    for key in event_dict:
        file_object.writelines(' '+str(key)+'           '+str(event_dict[key])+'\n')
    file_object.writelines('sum'+'    '+str(sum(event_dict.values()))+'\n')
    file_object.close()

def piechart(event_dict,path,name):
    plt.figure()
    labels=event_dict.keys()
    sizes=event_dict.values()
    plt.pie(sizes,labels=labels,autopct='%1.1f%%',startangle=50)
    plt.savefig(path+'/'+name+'pie.png')

class Alter_splice:

    def __init__(self,FAD,BPD,mark):
        self.BPD=BPD
        self.FAD=FAD
        self.mark = mark
        self.configparser_object=configparser()
        self.genome=self.configparser_object.get_values("Genome","genome")
        self.annotations=self.configparser_object.get_values("Genome","annotations")
        self.extract_as=self.configparser_object.get_values("Alter-splice","extract_as")
        self.hdrs=self.configparser_object.get_values("Alter-splice","hdrs")
        self.summarize_as=self.configparser_object.get_values("Alter-splice","summarize_as")
        self.ea_fpkm=self.configparser_object.get_values("Alter-splice","ea_fpkm")
        self.color=Color()
        self.run()

    def run(self):
        if not self.BPD.Alter_splice:
            self.color.flush_print("@" + mark + "Alter splice")
            self.color.print_debug("   @"+self.mark+"as-1 ...")
            shell_cmd="cuffcompare -r "+self.annotations+" -s "+self.genome+" -C -o "+self.FAD.alter_splice_dir+"/prefix_1 "+self.FAD.cufflinks_dir+"/transcripts.gtf"
            run(shell_cmd)
            self.color.print_debug("   @"+self.mark+"as-2 ...")
            shell_cmd=self.extract_as+" "+self.FAD.cufflinks_dir+"/transcripts.gtf "+self.hdrs+" -r "+self.FAD.cufflinks_dir+"/prefix_1.transcripts.gtf.tmap "+self.annotations+" >"+self.FAD.alter_splice_dir+"/lung_total.as"
            run(shell_cmd)
            shell_cmd=self.summarize_as+" "+self.FAD.cufflinks_dir+"/transcripts.gtf "+self.FAD.alter_splice_dir+"/lung_total.as -p "+self.FAD.alter_splice_dir+"/prefix_2"
            run(shell_cmd)

            self.color.print_debug("   @"+self.mark+"as-3 ...")
            shell_cmd=self.ea_fpkm+" "+self.FAD.cufflinks_dir+"/transcripts.gtf "+self.hdrs+" "+self.FAD.alter_splice_dir+"/prefix_2.as.nr >"+self.FAD.alter_splice_dir+"/fpkm_"+self.FAD.accession
            run(shell_cmd)

            self.color.print_debug("   @"+self.mark+"as-4 ...")
            ct=pd.read_table(self.FAD.alter_splice_dir+"/fpkm_"+self.FAD.accession)
            sped=[]
            i=0
            while i  < len(ct)-2:
                wink=[]
                wink.append(ct.iat[i,2])
                num=ct.iat[i,8]
                while (ct.iat[i,2]==ct.iat[i+1,2]):
                    num=num+ct.iat[i+1,8]
                    i=i+1
                wink.append(num)
                sped.append(wink)
                i=i+1

            self.alter_splice_save(sped,self.FAD.alter_splice_dir+"/"+self.FAD.accession+'.mini')
            f=open(self.FAD.alter_splice_dir+"/"+self.FAD.accession+'.mini')
            h=open(self.FAD.alter_splice_dir+"/"+self.FAD.accession+'.mini_nozore','w')
            for line in f:
                if line[-2]!='0' and line.startswith("ENSG"):
                    h.writelines(line)
            h.close()
            plot_go(self.FAD.alter_splice_dir,"fpkm_"+self.FAD.accession)

    def alter_splice_save(self,data,filename):
        file_object = open(filename, 'w')
        for i in range(len(data)):
            file_object.writelines(str(data[i][0])+'    '+str(data[i][1])+'\n')
        file_object.close()
