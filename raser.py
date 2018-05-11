#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#Example: job_name=xxx&&floder=../tumor/SRR493937&&pool=4&&threads=4&&ppn=`expr $pool \* $threads`&&echo "python ~/tools/2017/gc.py --threads $threads --pool $pool --floder $floder" |qsub -l nodes=$job_nd:ppn=$ppn -l walltime=1200:00:00 -e ${floder%/*} -o ${floder%/*} -N $job_name
__author__ = 'clsteam'
__time__ = '2018/4/20 22:26'

import multiprocessing
import argparse
import numpy as np
import time
from Testing_environment import *
from Management import *
from bag.Cluster import *
from bag.configparser import *
from disparity.ncRNA_detection_and_analysis import *
from disparity.alternative_splicing import *
from disparity.Fusion_gene_analysis import *
from disparity.allele_specific_expression import *
from disparity.different_gene_expression import *


def main():

    #Arguments passed
    parser = argparse.ArgumentParser()

    #required
    parser.add_argument("--floder",type=str,required=True,help="folder with all samples")
    parser.add_argument("--threads",required=True,type=str,default='4',help="n_cpu=n_pool*threads")
    parser.add_argument("--pool",type=int,required=True,default=1,help="Multithreading pool's numbers")

    global args
    args = parser.parse_args()

    #check parameters
    check_parameters()
    #start
    gogogo()

def check_parameters():
    global args, comparison, floder_list, config_object, environment, color
    environment=Testing_environment()
    environment.check()
    config_object = configparser()
    color = Color()
    if "," in args.floder:
        comparison=True
        floder_list=[get_absolute_path(args.floder.split(",")[0]),get_absolute_path(args.floder.split(",")[1])]
    else:
        comparison=False
        floder_list=[get_absolute_path(args.floder)]
    for floder in floder_list:
        if not os.path.isdir(floder):
            color.fatal_error(floder+" not existed")

    if not config_object.has_option("Genome", "genome_index"):
        color.print_warning("No genome index(bowtie2) can be found and build the index now")
        if environment.check_dependency("bowtie2") == False:
            color.fatal_error("External dependency 'bowtie2' not installed.");
        shell_cmd='bowtie2-build --threads '+str(int(args.threads)*args.pool)+' '+config_object.get_values("Genome", "genome")+' '+config_object.get_values("Genome", "genome")
        run(shell_cmd)
        config_object.add_values("Genome", "genome_index", config_object.get_values("Genome", "genome"))

def gogogo():
    global color
    color.flush_print('################################################')
    color.flush_print('#            Welcome to the RASER!             #')
    color.flush_print('#              Author: team-gc                 #')
    color.flush_print('#          Team: Burning Team(HZAU)            #')
    color.flush_print('################################################')
    color.flush_print('------------------------------------------------')
    color.flush_print("         Process pool :  %d" % args.pool)
    color.flush_print("          Per threads :  %s" % args.threads)
    color.flush_print('------------------------------------------------')

    tm1=time.time()		#start
    Independent_tasks()
    Non_Independent_tasks()
    tm2=time.time()		#end
    color.flush_print("USED TIME:%f hours" % ((tm2-tm1)/3600))
    color.flush_print("STOP TIME:"+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

def Independent_tasks():
    global args,color
    color.flush_print('doc: '+args.floder)
    color.flush_print('------------------------------------------------')
    jump_pool()

def Non_Independent_tasks():
    global floder_list,comparison,args
    if comparison:
        dge=Dge(floder_list,str(int(args.threads)*args.pool))
        dge.deseq2()
        dge.cuffdiff()


def jump_pool():
    global args
    mark=0
    data_hash=Generate_appropriate_data_format()
    pool=multiprocessing.Pool(processes=args.pool)
    for key in data_hash:
        mark+=1
        pool.apply_async(frame,(data_hash[key],str(mark)))
    pool.close()
    pool.join()

def Generate_appropriate_data_format():
    global floder_list,comparison
    if comparison:
        data_hash1 = generate_hash(floder_list[0],1)
        data_hash2 = generate_hash(floder_list[1],2)
        data_hash = dict(data_hash1.items() + data_hash2.items())
    else:
        data_hash = generate_hash(floder_list[0])
    return data_hash

def generate_hash(abs_floder,dispersion=0):
    if dispersion==1:
        j = 1
        interval = 2
    elif dispersion==2:
        j = 2
        interval =2
    else:
        j = 1
        interval = 1

    data = os.listdir(abs_floder)
    new_data = []
    #Remove floder,only leave documents
    for x in data:
        if os.path.isdir(abs_floder+"/"+x):
            data.remove(x)
        else:
            new_data.append(abs_floder + "/" + x)
    data=new_data
    del new_data

    data_hash = dict()
    i = 0
    while i < len(data):
        if data[i].endswith("sra"):
            format = "SRA"
            pat = re.compile('^(.*)(\.sra)$')
            home_dir = re.search(pat, data[i]).groups()[0]
            suffix = re.search(pat, data[i]).groups()[1]
            i += 1
        elif data[i].endswith("gz") or data[i].endswith("fq") or data[i].endswith("fastq") or data[i].endswith("sra"):
            pat = re.compile('^(.*?)(\.[fq|fastq].*)')
            if i == len(data) - 1:
                format = "SE"
                home_dir = re.search(pat, data[i]).groups()[0]
                suffix = re.search(pat, data[i]).groups()[1]
                i += 1
            elif re.search(pat, data[i]).groups()[0].endswith("1") and re.search(pat, data[i + 1]).groups()[0].endswith(
                    "2"):
                format = "PE"
                home_dir = re.search(pat, data[i]).groups()[0][:-1]
                suffix = re.search(pat, data[i]).groups()[1]
                i += 2
            else:
                format = "SE"
                home_dir = re.search(pat, data[i]).groups()[0]
                suffix = re.search(pat, data[i]).groups()[1]
                i += 1
        else:
            continue
        data_hash[j] = [format, home_dir, suffix]
        j += interval
    return data_hash

def generate_raw_data(home_dir,suffix,format,mark):
    global color
    if format == 'SRA':
        color.print_debug("@" + mark + "convert SRA to fastq ...")

        run('fastq-dump --split-3 ' + home_dir + '.sra' + ' -O ' + os.path.dirname(home_dir))
        if os.path.exists(home_dir + '.fastq'):
            format = 'SE'
            raw_data = home_dir + '.fastq'
        elif os.path.exists(home_dir+'_1.fastq'):
            format = 'PE'
            raw_data = [home_dir+ '_1.fastq', home_dir + '_2.fastq']
        else:
            color.print_warning("@" + mark + "file damaged:" + home_dir + '.sra')
    elif format == 'PE':
        raw_data=[home_dir+"1"+suffix,home_dir+"2"+suffix]
    else:
        raw_data=home_dir+suffix

    return raw_data

def frame(hash_list,mark):
    global config_object,args,comparison
    # hash_list: ['PE', '/public/home/yxu/lung-seq/normal/C6_R','.fq.gz', 'C6_R']
    [format,home_dir,suffix] = hash_list
    # format:SRA,PE,SE --> format:PE,SE   and generate raw_data
    raw_data=generate_raw_data(home_dir,suffix,format,mark)
    FAD=Files_and_directories(home_dir,format,args.threads)
    FAD.creat()
    BPD=Breakpoint_detection(home_dir,format)
    #Quality inspection and Sequence filtering
    ftf(FAD,BPD,raw_data,mark)
    #mapping
    if config_object.get_bool("Strand-secific","switch_on"):
        mapping(FAD,BPD,mark)
        library_type = strand_specific(FAD,mark)
        mapping(FAD,BPD,mark,library_type)
    else:
        mapping(FAD,BPD,mark)
        library_type=''

    #sort,rmdup,index
    deal_bam(FAD,BPD,mark)
    #cufflinks
    cufflinks(FAD,BPD,library_type,mark)

    #alter splice
    #Alter_splice(FAD, mark)

    #call snp and gene count,phaser
    ase=Ase(FAD,library_type,mark)
    ase.raw_gene_count()

    #fusion search
    #Fusion_search(FAD,comparison,mark)

    #ncRNA
    #Non_coding_RNA(FAD,mark)

def ftf(FAD,BPD,raw_data,mark):
    global config_object,color
    jar=config_object.get_values("Trimmomatic","trimmomatic_jar")
    if FAD.format == 'SE':
        #frist Quality inspection
        qc_txt_list=Fastqc(raw_data,1,FAD,BPD,mark)
        #Sequence filtering
        qc_result=qc_analysis(FAD,qc_txt_list)
        shell_cmd='java -jar '+jar+' SE '+fastq_phred_decide(raw_data)+' -threads '+FAD.threads+' '+raw_data+' '+FAD.clean_data+' '+qc_result
        color.print_debug("   @"+mark+"Trimmomatic parameters: "+qc_result)
        if not BPD.ftf: run(shell_cmd)
        #second Quality inspection
        Fastqc(FAD.clean_data,2,FAD,BPD,mark)
    else:
        # frist Quality inspection
        qc_txt_list=Fastqc(raw_data,1,FAD,BPD,mark)
        # Sequence filtering
        qc_result=qc_analysis(FAD,qc_txt_list)
        shell_cmd='java -jar '+jar+' PE '+fastq_phred_decide(raw_data[0])+' -threads '+FAD.threads+' '+raw_data[0]+' '+raw_data[1]+' '+FAD.clean_data[0]+' '+FAD.unclean_data[0]+' '+FAD.clean_data[1]+' '+FAD.unclean_data[1]+' '+qc_result
        color.print_debug("   @"+mark+"Trimmomatic parameters: "+qc_result)
        if not BPD.ftf: run(shell_cmd)
        # second Quality inspection
        Fastqc(FAD.clean_data,2,FAD,BPD,mark)

def fastq_phred_decide(fastq):
    # -phred33 -phred64
    phred_decide_pl=return_script_path()+"/disparity/fastq_phred_decide.pl"
    shell_cmd="perl "+phred_decide_pl+" "+fastq
    result=commands.getoutput(shell_cmd)
    return result

def Fastqc(data,step,FAD,BPD,mark):
    global color
    if FAD.format=="SE":
        if data.endswith(".fastq"):
            prefix = data[:-6]
        else:
            prefix = FAD.accession
        qc_zip=FAD.qc_dir+"/"+prefix+'_fastqc.zip'	#such as ~/data/SRR3291477/qc_out/SRR3291477_fastqc.zip
        qc_txt_list=[FAD.qc_dir+"/"+prefix+'_fastqc/fastqc_data.txt']
        shell_cmd='fastqc '+data+' -t '+FAD.threads+' -O '+FAD.qc_dir+'&&unzip '+qc_zip+' -d '+os.path.dirname(qc_zip)
    else:
        if data[0].endswith(".fastq"):
            prefix1 = data[0][:-6]
            prefix2 = data[1][:-6]
        else:
            prefix1 = FAD.accession + "1"
            prefix2 = FAD.accession + "2"
        qc_zip=FAD.qc_dir+"/"+prefix1+'_fastqc.zip'	#such as ~/data/SRR3291477/qc_out/SRR3291477_fastqc.zip
        qc_txt_list1=FAD.qc_dir+"/"+prefix1+'_fastqc/fastqc_data.txt'
        shell_cmd1='fastqc '+data[0]+' -t '+FAD.threads+' -O '+FAD.qc_dir+'&&unzip '+qc_zip+' -d '+os.path.dirname(qc_zip)
        qc_zip=FAD.qc_dir+"/"+prefix2+'_fastqc.zip'	#such as ~/data/SRR3291477/qc_out/SRR3291477_fastqc.zip
        qc_txt_list2=FAD.qc_dir+"/"+prefix2+'_fastqc/fastqc_data.txt'
        shell_cmd2='fastqc '+data[1]+' -t '+FAD.threads+' -O '+FAD.qc_dir+'&&unzip '+qc_zip+' -d '+os.path.dirname(qc_zip)
        shell_cmd=shell_cmd1+"&&"+shell_cmd2
        qc_txt_list=[qc_txt_list1,qc_txt_list2]

    if step==1:
        if not BPD.Fastqc_1:
            color.print_debug("@" + mark + "Fastqc...")
            run(shell_cmd)
    else:
        if not BPD.Fastqc_2:
            color.print_debug("@" + mark + "Fastqc(twice)...")
            run(shell_cmd)

    return qc_txt_list

def qc_analysis(FAD,qc_txt_list):
    #LEADING:15 TRAILING:15 AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:50
    global config_object
    Trimmomatic_home = config_object.get_values("Trimmomatic","home")
    original_adapter_file = config_object.get_values("Trimmomatic", FAD.format+"_adapters")
    os.chdir(Trimmomatic_home+"/adapters")
    run("cat "+original_adapter_file+">"+FAD.adapter_file)
    with open(FAD.adapter_file,'w+') as doc:
        if FAD.format == "SE" :
            trimmomatic_parameters=write_adapter(doc,qc_txt_list[0])
        else:
            trimmomatic_parameters=write_adapter(doc,qc_txt_list[0])
            trimmomatic_parameters=write_adapter(doc,qc_txt_list[1])
        return ' ILLUMINACLIP:'+FAD.adapter_file+':2:30:10 LEADING:'+config_object.get_values("Trimmomatic","LEADING")+' TRAILING:'+config_object.get_values("Trimmomatic","TRAILING")+trimmomatic_parameters+' SLIDINGWINDOW:'+config_object.get_values("Trimmomatic","SLIDINGWINDOW")

def write_adapter(doc,result):
    with open(result,'rt') as fq_open:
        data=fq_open.readlines()
        Sequence_length=int(data[8].strip('\n').split('\t')[1].split('-')[-1])
        #>>Per base sequence content\tfail\n
        trimmomatic_parameters='';next_step=0
        for a in range(len(data)):
            k=re.match(re.compile('^>>Per base sequence content\t(.*)'),data[a])
            if k :
                if k.group(1)=='pass':
                    break
                else:
                    next_step=1
                    per_index=a+2
            if(next_step==1 and data[a]=='>>END_MODULE\n'):
                per_content=data[per_index:a]
                trimmomatic_parameters+=search(per_content,Sequence_length)
                break

        #>>Overrepresented sequences\twarn\n
        next_step=0
        over_index=0
        for b in range(a,len(data)):
            k=re.match(re.compile('^>>Overrepresented sequences\t(.*)'),data[b])
            if k :
                if k.group(1) == 'pass' :
                    break
                else:
                    next_step=1
                    over_index=b+2
            if(next_step==1 and data[b]=='>>END_MODULE\n'):
                over_data=data[over_index:b]
                for line in over_data:
                    if line.split('\t')[-1] != 'No Hit\n' :
                        doc.writelines('>'+line.split('\t')[-1]+line.split('\t')[0]+'\n')
                break
    return trimmomatic_parameters

def search(ls,s_length):
    global color
    global home_dir
    b=[]
    head_crop=0
    crop=s_length
    min_num=len(ls)*0.05 #threshold (100bp allow 2bp is judgment error)
    for i in range(len(ls)):
        j=ls[i].strip('\n').split('\t')
        b.append([float(j[0]) if len(j[0])==1 else float(j[0][:2]),0 if abs(float(j[1])/float(j[4])-float(j[2])/float(j[3])) >0.1 else 1])	#search (A/T) - (G/C) > 10%
        if sum(np.array(b)[:,1]) <= min_num:
            head_crop=i+1
        else:
            break
    ls.reverse()
    for i in range(len(ls)):
        j=ls[i].strip('\n').split('\t')
        b.append([float(j[0]) if len(j[0])==1 else float(j[0][:2]),0 if abs(float(j[1])/float(j[4])-float(j[2])/float(j[3])) >0.1 else 1])	#search (A/T) - (G/C) > 10%
        if sum(np.array(b)[:,1]) <= 0:
            crop=len(ls)-i
        else:
            break
    if head_crop > 0.5*s_length:
        color.print_warning("the data("+os.path.basename(home_dir[:-1])+") noise is loud and maybe smRNA-Seq,HEADCROP will be set 0")
        crop=s_length
        head_crop=0
    minlen = 50 if (crop-head_crop)>50 else 18
    com=' CROP:'+str(crop)+' HEADCROP:'+str(head_crop)+' MINLEN:'+str(minlen)+' '
    return com

def mapping(FAD,BPD,mark,library_type=''):
    global color
    global config_object
    genome=config_object.get_values("Genome","genome")
    genome_index=config_object.get_values("Genome","genome_index")
    annotations=config_object.get_values("Genome","annotations")
    original_bam_dir=os.path.dirname(FAD.original_bam_file)
    if FAD.format == 'SE':
        color.flush_print("   @"+mark+"map(SE) ...")
        shell_cmd='tophat2 -p '+FAD.threads+' '+library_type+' -G '+annotations+' -o '+original_bam_dir+' '+genome_index+' '+FAD.clean_data
    else:
        color.flush_print("   @"+mark+"map(PE) ...")
        shell_cmd='tophat2 -p '+FAD.threads+' '+library_type+' -G '+annotations+' -o '+original_bam_dir+' '+genome_index+' '+FAD.clean_data[0]+" "+FAD.clean_data[1]
    if not BPD.mapping:
        run(shell_cmd)

def deal_bam(FAD,BPD,mark):
    global config_object,color
    if not BPD.deal_bam:
        if config_object.get_bool("Samtools","sort"):
            color.print_debug("@" + mark + "sort bam ...")
            run('samtools sort ' + FAD.original_bam_file + ' -o ' + FAD.original_bam_file)

        if config_object.get_bool("Samtools","rmdup"):
            color.print_debug("@" + mark + "remove PCR duplicate ...")
            if (FAD.format == 'SE'):
                run('samtools rmdup -s ' + FAD.original_bam_file +" "+ FAD.bam_file)
                #os.remove(FAD.original_bam_file)
            else:
                run('samtools rmdup -S ' + FAD.original_bam_file +" "+ FAD.bam_file)
                #os.remove(FAD.original_bam_file)
        else:
            run("mv "+FAD.original_bam_file + " " + FAD.bam_file)
            color.print_debug("   @" + mark + "build bam index ...")
        run('samtools index ' + FAD.bam_file)

def strand_specific(FAD,mark):
    global config_object,environment,color
    color.print_debug("@"+mark+"strand-specific...")
    if environment.check_dependency("infer_experiment.py") == False:
        color.fatal_error("External dependency 'RSeQC' not installed.")
    shell_cmd='infer_experiment.py -r '+config_object.get_values("Genome","bed")+' -i '+FAD.bam_file+'>'+FAD.home_dir+'/library_type'
    run(shell_cmd)
    with open("library_type","r") as doc:
        data=doc.readlines()
    a=float(re.findall(re.compile('.*\s(.*)\n'),data[-2])[0])
    b=float(re.findall(re.compile('.*\s(.*)\n'),data[-1])[0])
    if FAD.format == "SE":
        if (a-b)>0.5:
            library_type="--library-type fr-firststrand"
        elif (b-a)>0.5:
            library_type="--library-type fr-secondstrand"
        else:
            library_type="--library-type fr-unstranded"
    else:
        if (a-b)>0.5:
            library_type="--library-type fr-secondstrand"
        elif (b-a)>0.5:
            library_type="--library-type fr-firststrand"
        else:
            library_type="--library-type fr-unstranded"
    return library_type

def cufflinks(FAD,BPD,library_type,mark):
    global config_object,color
    color.flush_print("@" + mark + "expression statistics(cufflinks) ...")
    shell_cmd='cufflinks -p '+FAD.threads+' '+library_type+' -g '+config_object.get_values("Genome","annotations")+' -o '+FAD.cufflinks_dir+' '+FAD.bam_file
    if not BPD.cufflinks:
        run(shell_cmd)

def return_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

if __name__ == "__main__":
    main();
