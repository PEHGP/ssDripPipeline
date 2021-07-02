#!/usr/bin/env python
import sys,os,pandas,collections,re,itertools,math
from Bio import SeqIO
import numpy as np
class BasePipLine(object):
	"""decstring for BasePipLine"""
	def __init__(self, Prefix,Thread):
		self.Prefix = Prefix
		self.Thread=Thread
		if not os.path.exists(self.Prefix+"_pip.sh"):
			self.ShFr=open(self.Prefix+"_pip.sh","w")
		else:
			self.ShFr=open(self.Prefix+"_pip.sh","a")
	def Align(self,Left,Right,GenomeIndex):
		cmd="bowtie2 --local --phred33 -p %s -t -x %s -1 %s -2 %s 2>%s_align.info|samtools view -bS -1 |samtools sort -@ %s -m 5G -l 9 -o %s.sort.bam"%(self.Thread,GenomeIndex,Left,Right,self.Prefix,self.Thread,self.Prefix)
		os.system(cmd)
		self.ShFr.write(cmd+"\n")
	def ReomveDuplication(self,):
		InputBamFile=self.Prefix+".sort.bam"
		OutputBamFile=self.Prefix+".sort.paird_dup.bam"
		cmd="picard MarkDuplicates -REMOVE_DUPLICATES true -METRICS_FILE %s.matrix -INPUT %s -OUTPUT %s"%(self.Prefix,InputBamFile,OutputBamFile)
		os.system(cmd)
		self.ShFr.write(cmd+"\n")
		cmd="samtools index %s"%OutputBamFile
		os.system(cmd)
		self.ShFr.write(cmd+"\n")
	def SplitStrand(self,):
		InputBamFile=self.Prefix+".sort.paird_dup.bam"
		cmd1="samtools view -b -f 128 -F 16 %s > fwd1.bam"%(InputBamFile)
		os.system(cmd1)
		self.ShFr.write(cmd1+"\n")
		cmd2="samtools view -b -f 80 %s > fwd2.bam"%(InputBamFile)
		os.system(cmd2)
		self.ShFr.write(cmd2+"\n")
		cmd3="samtools merge -f %s_fwd.bam fwd1.bam fwd2.bam"%(self.Prefix)
		os.system(cmd3)
		self.ShFr.write(cmd3+"\n")
		cmd4="samtools view -b -f 144 %s > rev1.bam"%(InputBamFile)
		os.system(cmd4)
		self.ShFr.write(cmd4+"\n")
		cmd5="samtools view -b -f 64 -F 16 %s > rev2.bam"%(InputBamFile)
		os.system(cmd5)
		self.ShFr.write(cmd5+"\n")
		cmd6="samtools merge -f %s_rev.bam rev1.bam rev2.bam"%(self.Prefix)
		os.system(cmd6)
		self.ShFr.write(cmd6+"\n")
		cmd7="rm -rf fwd1.bam fwd2.bam rev1.bam rev2.bam"
		os.system(cmd7)
		self.ShFr.write(cmd7+"\n")
		cmd8="samtools index %s_fwd.bam"%(self.Prefix)
		os.system(cmd8)
		self.ShFr.write(cmd8+"\n")
		cmd9="samtools index %s_rev.bam"%(self.Prefix)
		os.system(cmd9)
		self.ShFr.write(cmd9+"\n")
		#self.ShFr.write("SplitBamAsStrand.sh %s.sort.picard_dup.bam\n"%self.Prefix)
	def CallPeak(self,GenomeSize,InputBamFile="",MyPrefix=""):
		if not InputBamFile:
			InputBamFile=self.Prefix+".sort.paird_dup.bam"
		if not MyPrefix:
			MyPrefix=self.Prefix
		cmd="macs2 callpeak -t %s -f BAMPE --keep-dup all -g %s -n %s"%(InputBamFile,GenomeSize,MyPrefix)
		self.ShFr.write(cmd+"\n")
		os.system(cmd)
	def FormatPeak(self,FilterChromFile,PeakExcelFile="",PeakBedFile="",Strand="+"):
		if not PeakExcelFile:
			PeakExcelFile=self.Prefix+"_peaks.xls"
		if not PeakBedFile:
			PeakBedFile=self.Prefix+"_peaks.bed"
		Fr=open(PeakBedFile,"w")
		FilterChromList=[x.rstrip() for x in open(FilterChromFile)]
		for x in open(PeakExcelFile):
			x=x.rstrip()
			if x.startswith("#") or x=="":
				continue
			l=x.split("\t")
			if l[0] in FilterChromList:
				continue
			if l[0]=="chr" and l[1]=="start" and l[2]=="end":
				continue
			Fr.write(l[0]+"\t"+l[1]+"\t"+l[2]+"\t"+l[-1]+"\t.\t"+Strand+"\n")
		Fr.close()
	def GetSpikeScale(self,DripBam,InputBam,SpikeNameList):
		DripLines=os.popen("samtools idxstats %s"%DripBam).readlines()
		InputLines=os.popen("samtools idxstats %s"%InputBam).readlines()
		Dripd={}
		Inputd={}
		DripSpikeReads=0
		DripRefReads=0
		InputSpikeReads=0
		InputRefReads=0
		for x,y in zip(DripLines,InputLines):
			x=x.rstrip()
			y=y.rstrip()
			lx=x.split("\t")
			ly=y.split("\t")
			if lx[0] in SpikeNameList:
				DripSpikeReads+=float(lx[2])
			else:
				DripRefReads+=float(lx[2])
			if ly[0] in SpikeNameList:
				InputSpikeReads+=float(ly[2])
			else:
				InputRefReads+=float(ly[2])
		#print(DripRefReads,InputRefReads,DripSpikeReads,InputSpikeReads)
		SpikeScale=(DripRefReads/InputRefReads)/(DripSpikeReads/InputSpikeReads)
		return SpikeScale
	def Get1XNormBw(self,BinSize,GenomeSize,IgnoreFile,SpikeScale=1):  #[5,95] ?
		DripBam=self.Prefix+".sort.paird_dup.bam"
		DripFwdBam=self.Prefix+"_fwd.bam"
		DripRevBam=self.Prefix+"_rev.bam"
		OutputBw=self.Prefix+"_nucleus_norm.bw"
		OutputFwdBw=self.Prefix+"_fwd_nucleus_norm.bw"
		OutputRevBw=self.Prefix+"_rev_nucleus_norm.bw"
		IgnoreList=[x.rstrip() for x in open(IgnoreFile)]
		Lines=os.popen("bamCoverage -v -p %s -b %s -o %s --binSize %s --effectiveGenomeSize %s --normalizeUsing RPGC --ignoreForNormalization %s"%(self.Thread,DripBam,OutputBw,BinSize,GenomeSize," ".join(IgnoreList))).readlines()
		for x in Lines:
			x=x.rstrip()
			if x.startswith("Final scaling factor:"):
				l=x.split()
				Scale=float(l[-1].strip())
				break
		ScaleFinal=Scale*SpikeScale
		cmd1="bamCoverage --scaleFactor %s -v -p %s -b %s -o %s --binSize %s --ignoreForNormalization %s >%s_norm.log"%(ScaleFinal,self.Thread,DripBam,OutputBw,BinSize," ".join(IgnoreList),self.Prefix)
		cmd2="bamCoverage --scaleFactor %s -v -p %s -b %s -o %s --binSize %s --ignoreForNormalization %s >>%s_norm.log"%(ScaleFinal,self.Thread,DripFwdBam,OutputFwdBw,BinSize," ".join(IgnoreList),self.Prefix)
		cmd3="bamCoverage --scaleFactor %s -v -p %s -b %s -o %s --binSize %s --ignoreForNormalization %s >>%s_norm.log"%(ScaleFinal,self.Thread,DripRevBam,OutputRevBw,BinSize," ".join(IgnoreList),self.Prefix)
		os.system(cmd1)
		os.system(cmd2)
		os.system(cmd3)
		self.ShFr.write(cmd1+"\n"+cmd2+"\n"+cmd3+"\n")
		return Scale
	def GetRandomRegionNormBw(self,IgnoreFile,Extend,ChromSize,RepeatNum): #not for spike-in
		FwdPeak="%s_fwd_peaks.bed"%self.Prefix
		RevPeak="%s_rev_peaks.bed"%self.Prefix
		self.GetRandomBed(FwdPeak,ChromSize,Strand="+",RepeatNum=RepeatNum,MyPrefix=self.Prefix+"_fwd",ExcludeFile="")
		self.GetRandomBed(RevPeak,ChromSize,Strand="-",RepeatNum=RepeatNum,MyPrefix=self.Prefix+"_rev",ExcludeFile="")
		cmd="cat %s_random_all.bed %s_random_all.bed >%s_random_all.bed"%(self.Prefix+"_fwd",self.Prefix+"_rev",self.Prefix)
		os.system(cmd)
		self.ShFr.write(cmd+"\n")
		RandomBed=self.Prefix+"_random_all.bed"
		IgnoreList=[x.rstrip() for x in open(IgnoreFile)]
		RL=[]
		for s in ["rev","fwd"]:
			if s=="fwd":
				Strand="+"
			elif s=="rev":
				Strand="-"
			cmd1="bamCoverage -v -p %s -b %s_%s.bam -o %s_%s.bw --binSize 1 --ignoreForNormalization %s"%(self.Thread,self.Prefix,s,self.Prefix,s," ".join(IgnoreList))
			os.system(cmd1)
			self.ShFr.write(cmd1+"\n")
			cmd2="computeMatrix scale-regions -p %s -S %s_%s.bw -R %s -bs 5 -b %s -a %s --skipZeros --outFileName %s_%s.gz"%(self.Thread,self.Prefix,s,RandomBed,Extend,Extend,self.Prefix,s)
			os.system(cmd2)
			self.ShFr.write(cmd2+"\n")
			df=pandas.read_csv("%s_%s.gz"%(self.Prefix,s),sep="\t",header=None,index_col=None,compression="gzip",skiprows=1)
			for i in df.columns[6:]:
				RL+=list(df.iloc[0:,i].dropna())
		RL=np.array(RL)
		n95=np.percentile(RL,95)
		n5=np.percentile(RL,5)
		Mean=np.mean(RL[(RL>n5)&(RL<n95)])
		scale=1/Mean
		for s in ["rev","fwd"]:
			cmd3="bamCoverage --extendReads --scaleFactor %s -v -p %s -b %s_%s.bam -o %s_%s_nucleus_norm.bw --binSize 1 --ignoreForNormalization %s"%(scale,self.Thread,self.Prefix,s,self.Prefix,s," ".join(IgnoreList))
			os.system(cmd3)
			self.ShFr.write(cmd3+"\n")
		cmd4="bamCoverage --extendReads --scaleFactor %s -v -p %s -b %s.sort.paird_dup.bam -o %s_nucleus_norm.bw --binSize 1 --ignoreForNormalization %s"%(scale,self.Thread,self.Prefix,self.Prefix," ".join(IgnoreList))
		os.system(cmd4)
		self.ShFr.write(cmd4+"\n")
		return scale
	def GetRandomBed(self,BedFile,ChromSize,Strand="+",RepeatNum=1,MyPrefix="",ExcludeFile=""):
		if not MyPrefix:
			MyPrefix=self.Prefix
		RandomFileList=[]
		for i in range(1,RepeatNum+1):
			if ExcludeFile:
				cmd1="bedtools shuffle -i %s -excl %s -g %s -chrom -noOverlapping|bedtools sort -i -|awk -F'\t' '{print $1\"\t\"$2\"\t\"$3\"\t%s_random%s_\"NR\"\t.\t%s\"}'>%s_random%s.bed"%(BedFile,ExcludeFile,ChromSize,MyPrefix,i,Strand,MyPrefix,i)
			else:
				cmd1="bedtools shuffle -i %s -g %s -chrom|bedtools sort -i -|awk -F'\t' '{print $1\"\t\"$2\"\t\"$3\"\t%s_random%s_\"NR\"\t.\t%s\"}'>%s_random%s.bed"%(BedFile,ChromSize,MyPrefix,i,Strand,MyPrefix,i)
			os.system(cmd1)
			self.ShFr.write(cmd1+"\n")
			RandomFileList.append("%s_random%s.bed"%(MyPrefix,i))
		cmd2="cat %s|bedtools sort -i - >%s_random_all.bed"%(" ".join(RandomFileList),MyPrefix)
		os.system(cmd2)
		self.ShFr.write(cmd2+"\n")
	def GetDeseqNormBw(self,):
		pass
	def AfterMath(self,):
		self.ShFr.close()
class AnalysisPipLine(object):
	"""decstring for AnalysisPipLine"""
	def __init__(self,Prefix,Thread):
		self.Prefix=Prefix
		self.Thread=Thread
		if not os.path.exists(self.Prefix+"_ana.sh"):
			self.ShFr=open(self.Prefix+"_ana.sh","w")
		else:
			self.ShFr=open(self.Prefix+"_ana.sh","a")
	def BasicStatistics(self,SampleDic,GenomeSize): #dict={"sample":["align.info","dup.matrix","peaks.xls","peaks.bed","fwd_peaks.bed","rev_peaks.bed"]}
		Dr={}
		Sl=[]
		for s in SampleDic:
			Lines=open(SampleDic[s][0]).read()
			m=re.search("(.*?) reads; of these:\n",Lines)
			TotalReads=m.group(1)
			m=re.search(".*? \((.*?)\) aligned concordantly exactly 1 time\n",Lines)
			Unique=m.group(1)
			m=re.search(".*? \((.*?)\) aligned concordantly >1 times\n",Lines)
			Multiple=m.group(1)
			m=re.search("(.*?) overall alignment rate\n",Lines)
			Overall=m.group(1)
			Dup=open(SampleDic[s][1]).readlines()[7].rstrip().split("\t")[-2]
			if len(SampleDic[s])==6:
				Lines=open(SampleDic[s][2]).read()
				m=re.search("# fragment size is determined as (.*?) bps\n",Lines)
				FragmentSize=m.group(1)
				PeakNum=len(open(SampleDic[s][3]).readlines())
				PeakPercent=os.popen("awk '{sum+=$3-$2}END{print sum}' %s"%(SampleDic[s][3])).readlines()[0].rstrip()
				PeakPercent=round(int(PeakPercent)/int(GenomeSize),4)
				FwdPeakNum=len(open(SampleDic[s][4]).readlines())
				FwdPeakPercent=os.popen("awk '{sum+=$3-$2}END{print sum}' %s"%(SampleDic[s][4])).readlines()[0].rstrip()
				FwdPeakPercent=round(int(FwdPeakPercent)/int(GenomeSize),4)
				RevPeakNum=len(open(SampleDic[s][5]).readlines())
				RevPeakPercent=os.popen("awk '{sum+=$3-$2}END{print sum}' %s"%(SampleDic[s][5])).readlines()[0].rstrip()
				RevPeakPercent=round(int(RevPeakPercent)/int(GenomeSize),4)
			else:
				FragmentSize="Nan"
				PeakNum="Nan"
				FwdPeakNum="Nan"
				RevPeakNum="Nan"
				PeakPercent="Nan"
				FwdPeakPercent="Nan"
				RevPeakPercent="Nan"
			Dr[s]=[TotalReads,Unique,Multiple,Overall,Dup,FragmentSize,PeakNum,FwdPeakNum,RevPeakNum,PeakPercent,FwdPeakPercent,RevPeakPercent]
			ScaleHeader=open(SampleDic[s][6]).readlines()[0].rstrip().split("\t")
			Sl.append(open(SampleDic[s][6]).readlines()[1].rstrip().split("\t"))
		DfScale=pandas.DataFrame(Sl,columns=ScaleHeader)
		DfScale.to_csv("%s_bwscale.xls"%self.Prefix,sep="\t",index=False)
		DfStat=pandas.DataFrame.from_dict(data=Dr,orient='index')
		DfStat.to_csv("%s_stat.xls"%self.Prefix,sep="\t",header=["TotalReads","UniqueReads","MultipleReads","OverallReads","Duplication","FragmentSize","PeakNum","FwdPeakNum","RevPeakNum","PeakPercent","FwdPeakPercent","RevPeakPercent"],index_label="Sample")
	def BwCorrelation(self,BwFileList,BinSize,IgnoreFile,Method,Type,MyPrefix="",OnlyPlot=False):  #Prefix may be need change
		""""spearman, pearson
		heatmap, scatterplot"""
		if not MyPrefix:
			MyPrefix=self.Prefix
		if not OnlyPlot:
			IgnoreList=[x.rstrip() for x in open(IgnoreFile)]
			cmd1="multiBigwigSummary bins -p %s -b %s -o %s_results.npz --binSize %s --chromosomesToSkip %s --outRawCounts %s_BinCounts.xls"%(self.Thread," ".join(BwFileList),MyPrefix,BinSize," ".join(IgnoreList),MyPrefix)
			os.system(cmd1)
			self.ShFr.write(cmd1+"\n")
		cmd2="plotCorrelation -in %s_results.npz --corMethod %s --skipZeros -p %s --plotTitle \"%s_correlation\" -o %s_%s_%s.png --plotFileFormat png --removeOutliers --outFileCorMatrix %s_%s_corr.tab"%(MyPrefix,Method,Type,MyPrefix,MyPrefix,Method,Type,MyPrefix,Method)
		os.system(cmd2)
		self.ShFr.write(cmd2+"\n")
	def FindMotif(self,GenomeFastaFile,BedFile,ChromSize,MyPrefix="",Strand="+",RepeatNum="10"):
		if not MyPrefix:
			MyPrefix=self.Prefix
		RandomFileList=[]
		for i in range(1,RepeatNum+1):
			cmd1="bedtools shuffle -i %s -excl %s -g %s -chrom -noOverlapping|bedtools sort -i -|awk -F'\t' '{print $1\"\t\"$2\"\t\"$3\"\t%s_random%s_\"NR\"\t.\t%s\"}'>%s_random%s.bed"%(BedFile,BedFile,ChromSize,MyPrefix,i,Strand,MyPrefix,i)
			os.system(cmd1)
			self.ShFr.write(cmd1+"\n")
			#cmd2="bedtools getfasta -fi %s -bed %s_random%s.bed -s -name >%s_random%s.fasta"%(GenomeFastaFile,MyPrefix,i,MyPrefix,i)
			#os.system(cmd2)
			#self.ShFr.write(cmd2+"\n")
			RandomFileList.append("%s_random%s.bed"%(MyPrefix,i))
		cmd3="bedtools getfasta -fi %s -bed %s -s -name >%s.fasta"%(GenomeFastaFile,BedFile,MyPrefix)
		os.system(cmd3)
		self.ShFr.write(cmd3+"\n")
		cmd4="cat %s|bedtools sort -i - >%s_random_all.bed"%(" ".join(RandomFileList),MyPrefix)
		os.system(cmd4)
		self.ShFr.write(cmd4+"\n")
		cmd5="bedtools getfasta -fi %s -bed %s_random_all.bed -s -name >%s_random_all.fasta"%(GenomeFastaFile,MyPrefix,MyPrefix)
		os.system(cmd5)
		self.ShFr.write(cmd5+"\n")
		cmd6="findMotifs.pl %s.fasta fasta %s_motif -fasta %s_random_all.fasta -len 8,10,12 -p %s -cache 50000"%(MyPrefix,MyPrefix,MyPrefix,self.Thread)
		os.system(cmd5)
		self.ShFr.write(cmd5+"\n")
	def PeakContentDistribution(self,PeakBed,PromoterBed,TerminaterBed,GenebodyBed): #GenebodyBed can be gene.bed,need random bed
		cmd1="wc -l %s"%PeakBed
		TotalNum=os.popen(cmd1).readlines()[0].split()[0].rstrip()
		self.ShFr.write(cmd1+"\n")
		cmd2="bedtools intersect -wa -c -a %s -b %s|awk -F'\t' '$NF!=0{print $0}'|wc -l"%(PeakBed,PromoterBed)
		ProNum=os.popen(cmd2).readlines()[0].rstrip()
		self.ShFr.write(cmd2+"\n")
		cmd3="bedtools intersect -v -a %s -b %s|bedtools intersect -wa -c -a - -b %s|awk -F'\t' '$NF!=0{print $0}'|wc -l "%(PeakBed,PromoterBed,TerminaterBed)
		TerNum=os.popen(cmd3).readlines()[0].rstrip()
		self.ShFr.write(cmd3+"\n")
		cmd4="bedtools intersect -v -a %s -b %s|bedtools intersect -v -a - -b %s|bedtools intersect -wa -c -a - -b %s|awk -F'\t' '$NF!=0{print $0}'|wc -l"%(PeakBed,PromoterBed,TerminaterBed,GenebodyBed)
		BodyNum=os.popen(cmd4).readlines()[0].rstrip()
		self.ShFr.write(cmd4+"\n")
		InterNum=int(TotalNum)-int(ProNum)-int(TerNum)-int(BodyNum)
		#reldist,jaccard,fisher
		return ProNum,TerNum,BodyNum,InterNum
	def PeakLengthDistribution(self,PeakBed):
		return map(lambda x:int(x[2])-int(x[1]),[x.rstrip().split("\t") for x in open(PeakBed)])
	def SenseAntisense(self,MyPrefix,BedFile,FwdBw,RevBw,Extend): #gene sense antisense 
		cmd1="grep '+$' %s >%s_positive.bed"%(BedFile,MyPrefix)
		os.system(cmd1)
		self.ShFr.write(cmd1+"\n")
		cmd2="grep '\-$' %s >%s_negative.bed"%(BedFile,MyPrefix)
		os.system(cmd2)
		self.ShFr.write(cmd2+"\n")
		BwFiled={"fwd":FwdBw,"rev":RevBw}
		for dr in ["fwd","rev"]:
			for zf in ["negative","positive"]:
				cmd3="computeMatrix scale-regions -p %s -S %s -R %s_%s.bed -bs 5 -b %s -a %s -m %s --skipZeros --outFileName %s_%s_%s.gz"%(self.Thread,BwFiled[dr],MyPrefix,zf,Extend,Extend,Extend,MyPrefix,zf,dr)
				os.system(cmd3)
				self.ShFr.write(cmd3+"\n")
		cmd4="computeMatrixOperations rbind -m %s_positive_rev.gz %s_negative_fwd.gz -o %s_antisense.gz"%(MyPrefix,MyPrefix,MyPrefix)
		os.system(cmd4)
		self.ShFr.write(cmd4+"\n")
		cmd5="computeMatrixOperations rbind -m %s_negative_rev.gz %s_positive_fwd.gz -o %s_sense.gz"%(MyPrefix,MyPrefix,MyPrefix)
		os.system(cmd5)
		self.ShFr.write(cmd5+"\n")
		#Calculated mean
	def BedCorrelation(self,MyPrefix,QueryBed,TargetBed,ChromSize,BwFile="",Extend=""): #permutation test,metaplot
		cmd1="bedtools shuffle -i %s -excl %s -g %s -chrom -noOverlapping|bedtools sort -i - >%s_random.bed"%(QueryBed,QueryBed,ChromSize,MyPrefix)
		os.system(cmd1)#random
		self.ShFr.write(cmd1+"\n")
		cmd2="bedtools fisher -a %s -b %s -g %s >%s_fisher.txt"%(QueryBed,TargetBed,ChromSize,MyPrefix)
		os.system(cmd2)
		self.ShFr.write(cmd2+"\n")
		cmd3="bedtools jaccard -a %s -b %s >%s_jaccard.txt"%(QueryBed,TargetBed,MyPrefix)
		os.system(cmd3)
		self.ShFr.write(cmd3+"\n")
		cmd4="bedtools reldist -a %s -b %s >%s_reldist.txt"%(QueryBed,TargetBed,MyPrefix)
		os.system(cmd4) #-detail
		self.ShFr.write(cmd4+"\n")
		cmd5="bedtools reldist -a %s_random.bed -b %s >%s_random_reldist.txt"%(MyPrefix,TargetBed,MyPrefix)
		os.system(cmd5)
		self.ShFr.write(cmd5+"\n")
		if BwFille:
			cmd6="computeMatrix scale-regions -p %s -S %s -R %s %s %s_random.bed -a %s -b %s -m %s --binSize 5 -o %s_matrix.gz --skipZeros"%(self.Thread,BwFile,QueryBed,TargetBed,MyPrefix,Extend,Extend,Extend,MyPrefix)
			os.system(cmd6)
			self.ShFr.write(cmd6+"\n")
			cmd7="plotProfile --dpi 300 --plotFileFormat png -m %s_matrix.gz -out %s_metaplot.png --plotTitle \"%s\""%(MyPrefix,MyPrefix,MyPrefix)
			os.system(cmd7)
			self.ShFr.write(cmd7+"\n")
		#reldist,jaccard,fisher
	def GetGenomeContentBedFile(self,GeneBed,Extend,ChromSize): #genebody tss, tts,intergenetic
		cmd1="awk -F\'\t\' \'{if($6==\"+\"){print $1\"\t\"$2\"\t\"$2\"\t\"$4\"\t\"$5\"\t\"$6}else{print $1\"\t\"$3\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6} }' %s|bedtools slop -i - -g %s -b %s >tss_%s.bed"%(GeneBed,ChromSize,Extend,Extend)
		os.system(cmd1)
		self.ShFr.write(cmd1+"\n")
		cmd2="awk -F\'\t\' \'{if($6==\"+\"){print $1\"\t\"$3\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}else{print $1\"\t\"$2\"\t\"$2\"\t\"$4\"\t\"$5\"\t\"$6} }' %s|bedtools slop -i - -g %s -b %s >tts_%s.bed"%(GeneBed,ChromSize,Extend,Extend)
		os.system(cmd2)
		self.ShFr.write(cmd2+"\n")
		cmd3="awk -F\'\t\' \'{Start=$2+%s;End=$3-%s;if(Start<End){print $1\"\t\"Start\"\t\"End\"\t\"$4\"\t\"$5\"\t\"$6}}' %s|bedtools sort -i - >genebody_%s.bed"%(Extend,Extend,GeneBed,Extend)
		os.system(cmd3)
		self.ShFr.write(cmd3+"\n")
		cmd4="cat tss_%s.bed tts_%s.bed genebody_%s.bed|bedtools merge -i -|bedtools sort -i -|bedtools complement -i - -g %s >intergenetic.bed"%(Extend,Extend,Extend,ChromSize)
		self.ShFr.write(cmd4+"\n")
		os.system(cmd4)
	def GetMergePeakCountMatrix(self,PeakDic,BamDic,MyPrefix=""):#PeakDic={"sample":"peaks.bed"},BamDic={"sample":"sample.bam"}
		if not MyPrefix:
			MyPrefix=self.Prefix
		cmd1="cat %s|bedtools sort -i -|bedtools merge -i - |awk -F'\t' '{print $0\"\tmerge_%s_\"NR\"\t.\t.\"}' >merge_%s.bed"%(" ".join(PeakDic.values()),MyPrefix,MyPrefix)
		os.system(cmd1)
		self.ShFr.write(cmd1+"\n")
		BamList=[]
		SampleList=[]
		for s in BamDic:
			SampleList.append(s)
			BamList.append(BamDic[s])
		cmd2="multiBamSummary BED-file --BED merge_%s.bed -p 20 --bamfiles %s --label %s --outRawCounts %s_counts.xls -o %s_counts.npz --scalingFactors %s_deseq_scale.txt"%(MyPrefix," ".join(BamList)," ".join(SampleList),MyPrefix,MyPrefix,MyPrefix)
		os.system(cmd2)#need change and add peak name
		self.ShFr.write(cmd2+"\n")
		cmd3="sed -i \"1s/[#|']//g\" %s_counts.xls"%(MyPrefix)
		os.system(cmd3)
		self.ShFr.write(cmd3+"\n")
		df=pandas.read_csv("%s_counts.xls"%MyPrefix,sep="\t",header=0,index_col=[0,1,2])
		df_bed=pandas.read_csv("merge_%s.bed"%MyPrefix,sep="\t",header=None,index_col=[0,1,2],names=["chr","start","end","peak","score","strand"])
		df_r=pandas.concat([df,df_bed],axis=1)
		labels=["peak"]+SampleList
		df_r.loc[:,labels].to_csv("%s_counts_final.xls"%MyPrefix,sep="\t",na_rep="Nan",index=False)
	def GetSenseAntisenseMatrix(self,BedFile,BamDic,MyPrefix=""):#BamDic={"fwd":{"sample":"sample_fwd.bam"},"rev":{"sample":"sample_rev.bam"}}
		if not MyPrefix:
			MyPrefix=self.Prefix
		Sense=pandas.DataFrame()
		Antisense=pandas.DataFrame()
		df_bed=pandas.read_csv(BedFile,sep="\t",header=None,index_col=[0,1,2],names=["chr","start","end","gene","score","strand"])
		for s in ["fwd","rev"]:
			BamList=[]
			SampleList=[]
			for si in BamDic[s]:
				SampleList.append(si)
				BamList.append(BamDic[s][si])
			os.system("multiBamSummary BED-file --BED %s -p 20 --bamfiles %s --label %s --outRawCounts %s_%s_counts.xls -o %s_%s.npz"%(BedFile," ".join(BamList)," ".join(SampleList),MyPrefix,s,MyPrefix,s))
			os.system("sed \"1s/[#|']//g\" %s_%s_counts.xls >%s_%s_counts_deal.xls"%(MyPrefix,s,MyPrefix,s))
			df=pandas.read_csv("%s_%s_counts_deal.xls"%(MyPrefix,s),sep="\t",header=0,index_col=[0,1,2])
			df_dup=df[df.index.duplicated(keep=False)]
			df_uni=df[~df.index.duplicated(keep=False)]
			df_bed_uni=df_bed[~df_bed.index.duplicated(keep=False)]
			df_con=pandas.concat([df_uni,df_bed_uni],axis=1)
			df_dup_con=pandas.DataFrame()
			for name,dfg in df_dup.groupby(df_dup.index):
				for i in  range(len(dfg)):
					df_dup_con=df_dup_con.append(pandas.concat([dfg.ix[i,:],df_bed.loc[name,:].ix[i,:]],axis=0))
			df_r=pandas.concat([df_con,df_dup_con],axis=0)
			df_r.loc[:,SampleList+["gene","score","strand"]].to_csv("%s_%s_counts_deal_addname.xls"%(MyPrefix,s),sep="\t",na_rep="NaN")
			Labels=["gene"]+SampleList
			df_r.loc[:,Labels].to_csv("%s_%s_counts_final.xls"%(MyPrefix,s),sep="\t",na_rep="Nan",index=False)
			if s=="fwd":
				Sense=Sense.append(df_r[df_r["strand"]=="+"].loc[:,Labels],ignore_index=True)
				Antisense=Antisense.append(df_r[df_r["strand"]=="-"].loc[:,Labels],ignore_index=True)
			elif s=="rev":
				Sense=Sense.append(df_r[df_r["strand"]=="-"].loc[:,Labels],ignore_index=True)
				Antisense=Antisense.append(df_r[df_r["strand"]=="+"].loc[:,Labels],ignore_index=True)
		Sense.to_csv("%s_sense_counts.xls"%MyPrefix,sep="\t",index=False)
		Antisense.to_csv("%s_antisense_counts.xls"%MyPrefix,sep="\t",index=False)
	def GetDeseq2File(self,TargetFile,CountMatrixFile,SpikeScaleDic={},DeseqScaleDic={},Control=None,MyPrefix=""):
		if not MyPrefix:
			MyPrefix=self.Prefix
		Dr={}
		Group=collections.defaultdict(list)
		for x in open(TargetFile):
			x=x.rstrip()
			l=x.split("\t")
			Dr[l[1]]=l[0]
			Group[l[0]].append(l[1])
		Fr=open(MyPrefix+"_deseq.r","w")
		Fr.write("#!/usr/bin/env Rscript\n")
		Fr.write("library(DESeq2)\nlibrary(gplots)\nlibrary(RColorBrewer)\nlibrary(genefilter)\nlibrary(calibrate)\n")
		Fr.write("countdata<-read.table(\"%s\",row.name=1,header=T,sep=\"\\t\",check.names=F)\n"%CountMatrixFile)
		l=open(CountMatrixFile).readlines()[0].rstrip().split("\t")
		factor=",".join(["\""+Dr[x]+"\"" for x in l[1:]])
		Fr.write("condition<-factor(c(%s))\n"%factor)
		Fr.write("coldata <- data.frame(row.names=colnames(countdata), condition)\n")
		Fr.write("dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)\n")
		if SpikeScaleDic:
			csl=[1/(SpikeScaleDic[s]*DeseqScaleDic[s]) for s in l[1:]]
			fm=""
			for si in l[1:]:
				fm+="\""+si+"\""+","
			Fr.write("sizefactor<-c(%s) #custom scale factor\n"%(",".join([str(si) for si in csl])))
			Fr.write("names(sizefactor) <- c(%s) #sample name\n"%(fm[:-1]))
			Fr.write("sizeFactors(dds) <- sizefactor\n")
		Fr.write("dds <- DESeq(dds)\n")
		Fr.write("norda<-counts(dds,normalized=TRUE)\n")
		Fr.write("exp<-as.data.frame(norda)\n")#normalized data
		Fr.write("vsd <- varianceStabilizingTransformation(dds, blind=FALSE)\n")
		Fr.write("pdf(\"%s_pca_deseq2.pdf\")\n"%MyPrefix)
		Fr.write("plotPCA(vsd, intgroup=\"condition\")\n")
		Fr.write("dev.off()\n")
		Fr.write("write.table(data.frame(\"Gene\"=rownames(exp),exp,check.names=F),file=\"%s_norm.xls\",sep=\"\\t\",quote=F,row.names=F)\n"%MyPrefix)
		for x,y in list(itertools.combinations(Group.keys(),2)):
			if x==Control:
				x,y=y,x
			co=",".join(["\""+i+"\"" for i in Group[x]])+","
			co+=",".join(["\""+i+"\"" for i in Group[y]])
			Fr.write("res <- results(dds,contrast=c(\'condition\',\'%s\',\'%s\'),addMLE=TRUE)\n"%(x,y))
			Fr.write("resdata <- merge(exp[,c(%s)],as.data.frame(res), by=\"row.names\", sort=FALSE)\n"%co)
			Fr.write("names(resdata)[1] <- \"Gene\"\n")
			Fr.write("write.table(resdata, file=\"%s_%s_%s_diffexpr_results.xls\",sep=\"\\t\",quote=F,row.names=F)\n"%(MyPrefix,x,y))
		Fr.close()
		os.system("chmod 755 %s_deseq.r"%MyPrefix)
		os.system("./%s_deseq.r"%MyPrefix)
	def SkewAsWindows(self,Seqstring,SeqId,Win,Step,Flag):
		Length=len(Seqstring)
		SkewList=[]
		flag=0
		for i in range(0,Length,Step):
			Start=i
			Stop=i+Win
			if Stop>=Length:
				Stop=Length
				flag=1
			SubSeq=Seqstring[Start:Stop].upper()
			if Flag=="GC":
				C=SubSeq.count("C")
				G=SubSeq.count("G")
				if G+C==0:
					Skew=0
				else:
					Skew=(G-C)/float(G+C)
			elif Flag=="AT":
				A=SubSeq.count("A")
				T=SubSeq.count("T")
				if A+T==0:
					Skew=0
				else:
					Skew=(A-T)/float(A+T)
			if Start<Stop:
				SkewList.append((SeqId,Start,Stop,Skew))
			if flag==1:
				break
		return SkewList
	def GetATGCSkewBw(self,GenomeFastaFile,ChromSize,Win,Step,MyPrefix=""):#win=100,step=50
		if not MyPrefix:
			MyPrefix=self.Prefix
		FrGC=open(MyPrefix+"_GCSkew.bdg","w")
		FrAT=open(MyPrefix+"_ATSkew.bdg","w")
		ChList=[x.rstrip().split("\t")[0] for x in open(ChromSize)]
		Shrink=int((Win-Step)/2)
		if (Win-Step)%2==0:
			ShrinkLeft=Shrink
			ShrinkRight=Shrink
		else:
			ShrinkLeft=Shrink+1
			ShrinkRight=Shrink
		if Shrink<=0:
			ShrinkLeft=0
			ShrinkRight=0
		for record in SeqIO.parse(GenomeFastaFile,"fasta"):
			Seqstring=record.seq
			SeqId=record.id
			if SeqId in ChList:
				GCSkewList=self.SkewAsWindows(Seqstring,SeqId,Win,Step,"GC")
				ATSkewList=self.SkewAsWindows(Seqstring,SeqId,Win,Step,"AT")
				for l in GCSkewList:
					Start=l[1]+ShrinkLeft
					End=l[2]-ShrinkRight
					FrGC.write(l[0]+"\t"+str(Start)+"\t"+str(End)+"\t"+str(l[3])+"\n")
				for l in ATSkewList:
					Start=l[1]+ShrinkLeft
					End=l[2]-ShrinkRight
					FrGC.write(l[0]+"\t"+str(Start)+"\t"+str(End)+"\t"+str(l[3])+"\n")
		FrGC.close()
		FrAT.close()
		cmd1="bedGraphToBigWig %s_GCSkew.bdg %s %s_GCSkew.bw"%(MyPrefix,ChromSize,MyPrefix)
		os.system(cmd1)
		self.ShFr.write(cmd1+"\n")
		cmd2="bedGraphToBigWig %s_ATSkew.bdg %s %s_ATSkew.bw"%(MyPrefix,ChromSize,MyPrefix)
		os.system(cmd2)
		self.ShFr.write(cmd2+"\n")
	def GetSkewGz(self,GCSkewBw,ATSkewBw,BedFile,Extend,MyPrefix):
		cmd="computeMatrix scale-regions -p %s -S %s %s -R %s -bs 5 -b %s -a %s -m %s --skipZeros --samplesLabel GCSkew ATSkew --outFileName %s_GCATSkew.gz"%(self.Thread,GCSkewBw,ATSkewBw,BedFile,Extend,Extend,Extend,MyPrefix)
		os.system(cmd)
		self.ShFr.write(cmd+"\n")
	def GetAnnoPeak(self,PeakBedFile,GeneBedFile,AnnoFileList):
		cmd="bedtools closest -a %s -b %s -D ref"%(PeakBedFile,GeneBedFile)
		self.ShFr.write(cmd+"\n")
		Annod=collections.defaultdict(list)
		for x in os.popen(cmd).readlines():
			x=x.rstrip()
			l=x.split("\t")
			Annod[l[3]].append(l[-4]+","+l[-2]+","+l[-1])
		for f in AnnoFileList:
			Name=f.split(".")[0]
			Fr=open(name+"_anno.xls","w")
			lines=open(f).readlines()
			Fr.write(lines[0].rstrip()+"\tanno\n")
			for x in lines[1:]:
				x=x.rstrip()
				l=x.split("\t")
				Fr.write(x+"\t"+";".join(d[l[0]])+"\n")
			Fr.close()
	def GetNoiseqFile(self,):
		pass
	def GetDiffRloopLevelFile(self,DiffList,TargetFile,LevelFile,MyPrefix=""): #LevelFile:*_counts_final.xls
		if not MyPrefix:
			MyPrefix=self.Prefix
		Gd=collections.defaultdict(list)
		for x in open(TargetFile):
			x=x.rstrip()
			l=x.split("\t")
			Gd[l[0]].append(l[1])
		UnionGene=[]
		for f in DiffList:
			UnionGene+=[x.rstrip() for x in os.popen("awk -F'\t' '$NF<=0.01{print $1}' %s"%f).readlines()] #0.01 may be needed change
		UnionGene=list(set(UnionGene))
		df=pandas.read_csv(LevelFile,sep="\t",header=0,index_col=0)
		df2=df.loc[UnionGene,:]
		df_r = pandas.DataFrame()
		for g in Gd:
			df_r[g]=df2.loc[:,Gd[g]].mean(axis=1)
		df_r.to_csv(MyPrefix+"_diff_union.xls",sep="\t")
	def GetCluster(self,DiffRloopLevelFile,MyPrefix=""): #mfuzz,DiffRloopLevelFile:Differential gene Union,There can be no duplicate samples,header=peak\tsample1\tsample2...\tsamplen
		SampleList=open(RloopLevelFile).readlines()[0].split("\t")[1:]
		SampleNum=len(SampleList)
		ClusterNumber=3**SampleNum
		si=SampleNum**0.5
		if SampleNum%si==0:
			r=c=int(si)
		else:
			r=round(si)
			if SampleNum%r==0:
				c=int(SampleNum/r)
			else:
				c=SampleNum//r+1
		Labels=",".join(["\""+x+"\"" for x in SampleList])
		Fr=open(MyPrefix+"_cluster.r","w")
		Fr.write("#!/usr/bin/env Rscript\n")
		Fr.write("library(Mfuzz)\n")
		Fr.write("rt<-read.table(\"%s\",sep=\"\t\",header=T,row.names=1,check.names=F)\n"%DiffRloopLevelFile)
		Fr.write("data = new('ExpressionSet', exprs=as.matrix(rt))\n")
		Fr.write("data_filter<-filter.std(data,min.std=0)\n")
		Fr.write("rt2<-standardise(data_filter)\n")
		Fr.write("m1=mestimate(rt2)\n")
		Fr.write("cl <- mfuzz(rt2,c=%s,m=m1)\n"%ClusterNumber)
		Fr.write("pdf(file=\"%s_cluster.pdf\")\n"%MyPrefix)
		Fr.write("mfuzz.plot2(rt2,cl=cl,mfrow=c(%s,%s),x11=F,centre=T,centre.col=\"coral\",colo=rev(grey.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha =0.3)),time.labels=c(%s))\n"%(r,c,Labels)) #mfrow=c(2,2)
		Fr.write("dev.off()\n")
		Fr.write("fr<-cbind(rt,\"cluster\"=cl$cluster,\"membership\"=cl$membership)\n")
		Fr.write("write.table(data.frame(\"Gene\"=rownames(fr),fr),file=\"%s_cluster.xls\",sep=\"\t\",quote=F,row.names=FALSE)\n"%MyPrefix)
		Fr.close()
		os.system("chmod 755 %s_cluster.r"%MyPrefix)
		os.system("./%s_cluster.r"%MyPrefix)
	def AfterMath(self,):
		self.ShFr.close()