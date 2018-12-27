#!/usr/bin/env python2.7
import sys,os,pandas
class BasePipLine(object):
	"""docstring for DripPipLine"""
	def __init__(self, Prefix,Thread):
		self.Prefix = Prefix
		self.Thread=Thread
	def Align(self,Left,Right,GenomeIndex):
		os.system("bowtie2 --local --phred33 -p %s -t -x %s -1 %s -2 %s 2>%s_align.info|samtools view -bS -1 |samtools sort -@ %s -m 5G -l 9 -o %s.sort.bam"%(self.Thread,GenomeIndex,Left,Right,self.Prefix,self.Thread,self.Prefix))
	def ReomveDuplication(self,InputBamFile=self.Prefix+".sort.bam",OutputBamFile=self.Prefix+".sort.paird_dup.bam"):
		os.system("java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=%s.matrix INPUT=%s OUTPUT=%s"%(self.Prefix,InputBamFile,OutputBamFile))
	def SplitStrand(self,InputBamFile=self.Prefix+".sort.paird_dup.bam"):
		os.system("samtools view -b -f 128 -F 16 %s > fwd1.bam"%(InputBamFile))
		os.system("samtools view -b -f 80 %s > fwd2.bam"%(InputBamFile))
		os.system("samtools merge -f %s_fwd.bam fwd1.bam fwd2.bam"%(self.Prefix))
		os.system("samtools view -b -f 144 %s > rev1.bam"%(InputBamFile))
		os.system("samtools view -b -f 64 -F 16 %s > rev2.bam"%(InputBamFile))
		os.system("samtools merge -f %s_rev.bam rev1.bam rev2.bam"%(self.Prefix))
		os.system("rm -rf fwd1.bam fwd2.bam rev1.bam rev2.bam")
	def CallPeak(self,GenomeSize,InputBamFile=self.Prefix+".sort.paird_dup.bam"):
		os.system("macs2 callpeak -t %s -f BAMPE -g %s -n %s"%(InputBamFile,GenomeSize,self.Prefix))
	def FormatPeak(self,FilterChromFile,PeakExcelFile=self.Prefix+"_peaks.xls",PeakBedFile=self.Prefix+"_peaks.bed",Strand="+"):
		AwkRegular=""
		for x in open(FilterChromFile):
			x=x.rstrip()
			AwkRegular+="^"+x+"$"+"|"
		AwkRegular=AwkRegular[:-1]
		os.system("awk -F'\\t' '$0!~/^#/&&$0!=\"\"&&$1!~/%s/{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$NF\"\\t.\\t%s\"}' %s  > %s"%(AwkRegular,Strand,PeakExcelFile,PeakBedFile))
	def GetNormBw(self,BinSize,GenomeSize,IgnoreFile,InputBam=self.Prefix+".sort.paird_dup.bam",InputFwdBam=self.Prefix+"_fwd.bam",InputRevBam=self.Prefix+"_rev.bam",OutputBam=self.Prefix+"_nucleus_norm.bw",OutputFwdBam=self.Prefix+"_fwd_nucleus_norm.bw",OutputRevBam=self.Prefix+"_rev_nucleus_norm.bw"):  #[5,95] ?
		IgnoreList=[x.rstrip() for x in open(IgnoreFile)]
		Lines=os.popen("bamCoverage --extendReads -v -p %s -b %s -o %s --binSize %s --effectiveGenomeSize %s --normalizeUsing RPGC --ignoreForNormalization %s"%(self.Thread,InputBam,OutputBam,BinSize,GenomeSize," ".join(IgnoreList))).readlines()
		for x in Lines:
			x=x.rstrip()
			if x.startswith("Final scaling factor:"):
				l=x.split()
				Scale=l[-1].strip()
				break
		os.system("bamCoverage --extendReads --scaleFactor %s -v -p %s -b %s -o %s --binSize %s --ignoreForNormalization %s"%(Scale,self.Thread,InputFwdBam,OutputFwdBam,BinSize," ".join(IgnoreList)))
		os.system("bamCoverage --extendReads --scaleFactor %s -v -p %s -b %s -o %s --binSize %s --ignoreForNormalization %s"%(Scale,self.Thread,InputRevBam,OutputRevBam,BinSize," ".join(IgnoreList)))
		Fr=open(self.Prefix+"_scale.txt","w")
		Fr.write(self.Prefix+"\t"+Scale+"\n")
		Fr.close()
class AnalysisPipLine(object):
	"""heheh"""
	def __init__(self,Prefix,Thread):
		self.Prefix=Prefix
		self.Thread=Thread
	def BasicStatistics(self,SampleDic): #dict={"sample":["align.info","dup.matrix","peaks.xls","peaks.bed","fwd_peaks.bed","rev_peaks.bed"]}
		Dr={}
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
			Lines=open(SampleDic[s][2]).read()
			m=re.search("# fragment size is determined as (.*?) bps\n",Lines)
			FragmentSize=m.group(1)
			PeakNum=len(open(SampleDicp[s][3]).readlines())
			FwdPeakNum=len(open(SampleDicp[s][4]).readlines())
			RevPeakNum=len(open(SampleDicp[s][5]).readlines())
			Dr[s]=[TotalReads,Unique,Multiple,Overall,Dup,FragmentSize,PeakNum,FwdPeakNum,RevPeakNum]
		df=pandas.DataFrame.from_dict(data=Dr,orient='index')
		df.to_csv("%s_stat.xls"%self.Prefix,sep="\t")
	def BwCorrelation(self,BwFileList,BinSize,IgnoreFile,Method,Type):  #Prefix may be need change
		""""spearman, pearson
		heatmap, scatterplot"""
		IgnoreList=[x.rstrip() for x in open(IgnoreFile)]
		os.system("multiBigwigSummary bins -p %s -b %s -o %s_results.npz --binSize %s --chromosomesToSkip %s --outRawCounts %s_BinCounts.xls"%(self.Thread," ".join(BwFileList),self.Prefix,BinSize," ".join(IgnoreList),self.Prefix))
		os.system("plotCorrelation -in %s_results.npz --corMethod %s --skipZeros -p %s --plotTitle \"%s_correlation\" -o %s_%s.pdf --plotFileFormat pdf --removeOutliers --outFileCorMatrix %s_corr.tab"%(self.Prefix,Method,Type,self.Prefix,self.Prefix,Method,self.Prefix))
	def FindMotif(self,MyPrefix,GenomeSeq,BedFile,ChromSize,Strand="+",RepeatNum="10"):
		RandomFileList=[]
		for i in range(1,RepeatNum+1):
			os.system("bedtools shuffle -i %s -excl %s -g %s -chrom -noOverlapping|bedtools sort -i -|awk -F'\t' '{print $1\"\t\"$2\"\t\"$3\"\t%s_random%s_\"NR\"\t.\t%s\"}'>%s_random%s.bed"%(BedFile,BedFile,ChromSize,MyPrefix,i,Strand,MyPrefix,i))
			os.system("bedtools getfasta -fi %s -bed %s_random%s.bed -s -name >%s_random%s.fasta"%(GenomeSeq,MyPrefix,i,MyPrefix,i))
			RandomFileList.append("%s_random%s.fasta"%(RandomFileList,i))
		os.system("bedtools getfasta -fi %s -bed %s -s -name >%s.fasta"%(GenomeSeq,BedFile,MyPrefix))
		os.system("cat %s >%s_random_all.fasta"%(" ".join(RandomFileList),MyPrefix))
		os.system("findMotifs.pl %s.fasta fasta %s_motif -fasta %s_random_all.fasta -len 8,10,12 -p %s -cache 50000"%(MyPrefix,MyPrefix,MyPrefix,self.Thread))
	def PeakContentDistribution(self,MyPrefix,PeakBed,PromoterBed,TerminaterBed,GenebodyBed): #GenebodyBed can be gene.bed,need random bed
		TotalNum=os.popen("wc -l %s"%PeakBed).readlines()[0].split()[0].rstrip()
		ProNum=os.popen("bedtools intersect -wa -c -a %s -b %s|awk -F'\t' '$NF!=0{print $0}'|wc -l"%(PeakBed,PromoterBed)).readlines()[0].rstrip()
		TerNum=os.popen("bedtools intersect -v -a %s -b %s|bedtools intersect -wa -c -a - -b %s|awk -F'\t' '$NF!=0{print $0}'|wc -l "%(PeakBed,PromoterBed,TerminaterBed)).readlines()[0].rstrip()
		BodyNum=os.popen("bedtools intersect -v -a %s -b %s|bedtools intersect -v -a - -b %s|bedtools intersect -wa -c -a - -b %s|awk -F'\t' '$NF!=0{print $0}'|wc -l"%(PeakBed,PromoterBed,TerminaterBed,GenebodyBed)).readlines()[0].rstrip()
		InterNum=int(TotalNum)-int(ProNum)-int(TerNum)-int(BodyNum)
		#reldist,jaccard,fisher
		return ProNum,TerNum,BodyNum,InterNum
	def PeakLengthDistribution(self,PeakBed):
		return map(lambda x:int(x[2])-int(x[1])+1,[x.rstrip().split("\t") for x in open(PeakBed)])
	def SenseAntisense(self,MyPrefix,BedFile,FwdBw,RevBw,Extend): #gene sense antisense 
		os.system("grep '+$' %s >%s_positive.bed"%(BedFile,MyPrefix))
		os.system("grep '\-$' %s >%s_negative.bed"%(BedFile,MyPrefix))
		BwFiled={"fwd":FwdBw,"rev":RevBw}
		for dr in ["fwd","rev"]:
			for zf in ["negative","positive"]:
				os.system("computeMatrix scale-regions -p 20 -S %s -R %s_%s.bed -bs 5 -b %s -a %s --skipZeros --outFileName %s_%s_%s.gz"%(BwFiled[dr],MyPrefix,zf,Extend,Extend,MyPrefix,zf,dr))
		os.system("computeMatrixOperations rbind -m %s_positive_rev.gz %s_negative_fwd.gz -o %s_antisense.gz"%(MyPrefix,MyPrefix,MyPrefix))
		os.system("computeMatrixOperations rbind -m %s_negative_rev.gz %s_positive_fwd.gz -o %s_sense.gz"%(MyPrefix,MyPrefix,MyPrefix))
	def BedCorrelation(self,BedFile,TargetBedFile,BwFile): #permutation test,metaplot
		#reldist,jaccard,fisher
	def GetGenomeContentBedFile(self,): #genebody tss, tts,intergenetic
		pass
	def GetNoiseqFile(self,):
		pass
	def GetDeseq2File(self,):
		pass