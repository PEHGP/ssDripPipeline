#!/usr/bin/env python2.7
import sys,os,pandas
class BasePipLine(object):
	"""docstring for DripPipLine"""
	def __init__(self, Prefix,Thread):
		self.Prefix = Prefix
		self.Thread=Thread
	def Align(self,Left,Right,GenomeIndex):
		os.system("bowtie2 --local --phred33 -p %s -t -x %s -1 %s -2 %s 2>%s_align.info|samtools view -bS -1 |samtools sort -@ %s -m 5G -l 9 -o %s.sort.bam"%(self.Thread,GenomeIndex,Left,Right,self.Prefix,self.Thread,self.Prefix))
	def ReomveDuplication(self,InputBamFile="",OutputBamFile=""):
		if not InputBamFile:
			InputBamFile=self.Prefix+".sort.bam"
		if not OutputBamFile:
			OutputBamFile=self.Prefix+".sort.paird_dup.bam"
		os.system("java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=%s.matrix INPUT=%s OUTPUT=%s"%(self.Prefix,InputBamFile,OutputBamFile))
	def SplitStrand(self,InputBamFile=""):
		if not InputBamFile:
			InputBamFile=self.Prefix+".sort.paird_dup.bam"
		os.system("samtools view -b -f 128 -F 16 %s > fwd1.bam"%(InputBamFile))
		os.system("samtools view -b -f 80 %s > fwd2.bam"%(InputBamFile))
		os.system("samtools merge -f %s_fwd.bam fwd1.bam fwd2.bam"%(self.Prefix))
		os.system("samtools view -b -f 144 %s > rev1.bam"%(InputBamFile))
		os.system("samtools view -b -f 64 -F 16 %s > rev2.bam"%(InputBamFile))
		os.system("samtools merge -f %s_rev.bam rev1.bam rev2.bam"%(self.Prefix))
		os.system("rm -rf fwd1.bam fwd2.bam rev1.bam rev2.bam")
	def CallPeak(self,GenomeSize,MyPrefix="",InputBamFile=""):
		if not InputBamFile:
			InputBamFile=self.Prefix+".sort.paird_dup.bam"
		if not MyPrefix:
			MyPrefix=self.Prefix
		os.system("macs2 callpeak -t %s -f BAMPE --keep-dup all -g %s -n %s"%(InputBamFile,GenomeSize,MyPrefix))
	def FormatPeak(self,FilterChromFile,PeakExcelFile="",PeakBedFile="",Strand="+"):
		if not PeakExcelFile:
			PeakExcelFile=self.Prefix+"_peaks.xls"
		if not PeakBedFile:
			PeakBedFile=self.Prefix+"_peaks.bed"
		AwkRegular=""
		for x in open(FilterChromFile):
			x=x.rstrip()
			AwkRegular+="^"+x+"$"+"|"
		AwkRegular=AwkRegular[:-1]
		os.system("awk -F'\\t' '$0!~/^#/&&$0!=\"\"&&$1!~/%s/{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$NF\"\\t.\\t%s\"}' %s  > %s"%(AwkRegular,Strand,PeakExcelFile,PeakBedFile))
	def GetNormBw(self,BinSize,GenomeSize,IgnoreFile,InputBam="",InputFwdBam="",InputRevBam="",OutputBw="",OutputFwdBw="",OutputRevBw=""):  #[5,95] ?
		if not InputBam:
			InputBam=self.Prefix+".sort.paird_dup.bam"
		if not InputFwdBam:
			InputFwdBam=self.Prefix+"_fwd.bam"
		if not InputRevBam:
			InputRevBam=self.Prefix+"_rev.bam"
		if not OutputBw:
			OutputBw=self.Prefix+"_nucleus_norm.bw"
		if not OutputFwdBw:
			OutputFwdBw=self.Prefix+"_fwd_nucleus_norm.bw"
		if not OutputRevBw:
			OutputRevBw=self.Prefix+"_rev_nucleus_norm.bw"
		IgnoreList=[x.rstrip() for x in open(IgnoreFile)]
		Lines=os.popen("bamCoverage --extendReads -v -p %s -b %s -o %s --binSize %s --effectiveGenomeSize %s --normalizeUsing RPGC --ignoreForNormalization %s"%(self.Thread,InputBam,OutputBw,BinSize,GenomeSize," ".join(IgnoreList))).readlines()
		for x in Lines:
			x=x.rstrip()
			if x.startswith("Final scaling factor:"):
				l=x.split()
				Scale=l[-1].strip()
				break
		os.system("bamCoverage --extendReads --scaleFactor %s -v -p %s -b %s -o %s --binSize %s --ignoreForNormalization %s"%(Scale,self.Thread,InputFwdBam,OutputFwdBw,BinSize," ".join(IgnoreList)))
		os.system("bamCoverage --extendReads --scaleFactor %s -v -p %s -b %s -o %s --binSize %s --ignoreForNormalization %s"%(Scale,self.Thread,InputRevBam,OutputRevBw,BinSize," ".join(IgnoreList)))
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
		df.to_csv("%s_stat.xls"%self.Prefix,sep="\t",header=["TotalReads","UniqueReads","MultipleReads","OverallReads","Duplication","FragmentSize","PeakNum","FwdPeakNum","RevPeakNum"],index_label="Sample")
	def BwCorrelation(self,BwFileList,BinSize,IgnoreFile,Method,Type,MyPrefix="",OnlyPlot=False):  #Prefix may be need change
		""""spearman, pearson
		heatmap, scatterplot"""
		if not MyPrefix:
			MyPrefix=self.Prefix
		if not OnlyPlot:
			IgnoreList=[x.rstrip() for x in open(IgnoreFile)]
			os.system("multiBigwigSummary bins -p %s -b %s -o %s_results.npz --binSize %s --chromosomesToSkip %s --outRawCounts %s_BinCounts.xls"%(self.Thread," ".join(BwFileList),MyPrefix,BinSize," ".join(IgnoreList),MyPrefix))
		os.system("plotCorrelation -in %s_results.npz --corMethod %s --skipZeros -p %s --plotTitle \"%s_correlation\" -o %s_%s_%s.pdf --plotFileFormat pdf --removeOutliers --outFileCorMatrix %s_corr.tab"%(MyPrefix,Method,Type,MyPrefix,MyPrefix,Method,Type,MyPrefix))
	def FindMotif(self,MyPrefix,GenomeFastaFile,BedFile,ChromSize,Strand="+",RepeatNum="10"):
		RandomFileList=[]
		for i in range(1,RepeatNum+1):
			os.system("bedtools shuffle -i %s -excl %s -g %s -chrom -noOverlapping|bedtools sort -i -|awk -F'\t' '{print $1\"\t\"$2\"\t\"$3\"\t%s_random%s_\"NR\"\t.\t%s\"}'>%s_random%s.bed"%(BedFile,BedFile,ChromSize,MyPrefix,i,Strand,MyPrefix,i))
			os.system("bedtools getfasta -fi %s -bed %s_random%s.bed -s -name >%s_random%s.fasta"%(GenomeFastaFile,MyPrefix,i,MyPrefix,i))
			RandomFileList.append("%s_random%s.fasta"%(RandomFileList,i))
		os.system("bedtools getfasta -fi %s -bed %s -s -name >%s.fasta"%(GenomeFastaFile,BedFile,MyPrefix))
		os.system("cat %s >%s_random_all.fasta"%(" ".join(RandomFileList),MyPrefix))
		os.system("findMotifs.pl %s.fasta fasta %s_motif -fasta %s_random_all.fasta -len 8,10,12 -p %s -cache 50000"%(MyPrefix,MyPrefix,MyPrefix,self.Thread))
	def PeakContentDistribution(self,PeakBed,PromoterBed,TerminaterBed,GenebodyBed): #GenebodyBed can be gene.bed,need random bed
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
				os.system("computeMatrix scale-regions -p 20 -S %s -R %s_%s.bed -bs 5 -b %s -a %s -m %s --skipZeros --outFileName %s_%s_%s.gz"%(BwFiled[dr],MyPrefix,zf,Extend,Extend,Extend,MyPrefix,zf,dr))
		os.system("computeMatrixOperations rbind -m %s_positive_rev.gz %s_negative_fwd.gz -o %s_antisense.gz"%(MyPrefix,MyPrefix,MyPrefix))
		os.system("computeMatrixOperations rbind -m %s_negative_rev.gz %s_positive_fwd.gz -o %s_sense.gz"%(MyPrefix,MyPrefix,MyPrefix))
		#Calculated mean
	def BedCorrelation(self,MyPrefix,QueryBed,TargetBed,ChromSize,BwFile="",Extend=""): #permutation test,metaplot
		os.system("bedtools shuffle -i %s -excl %s -g %s -chrom -noOverlapping|bedtools sort -i - >%s_random.bed"%(QueryBed,QueryBed,ChromSize,MyPrefix))#random
		os.system("bedtools fisher -a %s -b %s -g %s >%s_fisher.txt"%(QueryBed,TargetBed,ChromSize,MyPrefix))
		os.system("bedtools jaccard -a %s -b %s >%s_jaccard.txt"%(QueryBed,TargetBed,MyPrefix))
		os.system("bedtools reldist -a %s -b %s >%s_reldist.txt"%(QueryBed,TargetBed,MyPrefix)) #-detail
		os.system("bedtools reldist -a %s_random.bed -b %s >%s_random_reldist.txt"%(MyPrefix,TargetBed,MyPrefix))
		if BwFille:
			os.system("computeMatrix scale-regions -p %s -S %s -R %s %s  %s_random.bed -a %s -b %s -m %s --binSize 5 -o %s_matrix.gz --skipZeros"%(self.Thread,BwFile,QueryBed,TargetBed,MyPrefix,Extend,Extend,Extend,MyPrefix))
			os.system("plotProfile --dpi 300 --plotFileFormat pdf -m %s_matrix.gz -out %s_metaplot.pdf --plotTitle \"%s\""%(MyPrefix,MyPrefix,MyPrefix))
		#reldist,jaccard,fisher
	def GetGenomeContentBedFile(self,GeneBeb,Extend,ChromSize): #genebody tss, tts,intergenetic
		os.system("awk -F\'\t\' \'{if($6==\"+\"){print $1\"\t\"$2\"\t\"$2\"\t\"$4\"\t\"$5\"\t\"$6}else{print $1\"\t\"$3\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6} }' %s|bedtools slop -i - -g %s -b %s >tss_%s.bed"%(GeneBed,ChromSize,Extend,Extend))
		os.system("awk -F\'\t\' \'{if($6==\"+\"){print $1\"\t\"$3\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}else{print $1\"\t\"$2\"\t\"$2\"\t\"$4\"\t\"$5\"\t\"$6} }' %s|bedtools slop -i - -g %s -b %s >tts_%s.bed"%(GeneBed,ChromSize,Extend,Extend))
		os.system("awk -F\'\t\' \'{Start=$2+%s;End=$3-%s;if(Start<End){print $1\"\t\"Start\"\t\"End\"\t\"$4\"\t\"$5\"\t\"$6}}' %s|bedtools sort -i - >genebody_%s.bed"%(Extend,Extend,GeneBed,Extend))
		os.system("cat tss_%s.bed tts_%s.bed genebody_%s.bed|bedtools merge -i -|bedtools sort -i -|bedtools complement -i - -g %s >intergenetic.bed"%(Extend,Extend,Extend,ChromSize))
	def GetMergePeakCountMatrix(self,PeakDic,BamDic,MyPrefix=""):#PeakDic={"sample":"peaks.bed"},BamDic={"sample":"sample.bam"}
		if not MyPrefix:
			MyPrefix=self.Prefix
		os.system("cat %s|bedtools sort -i -|bedtools merge -i - |awk -F'\t' '{print $0\"\tmerge_%s_\"NR\"\t.\t.\"}' >merge_%s.bed"%(" ".join(PeakDic.values()),MyPrefix,MyPrefix))
		BamList=[]
		SampleList=[]
		for s in BamDic:
			SampleList.append(s)
			BamList.append(BamDic[s])
		os.system("multiBamSummary BED-file --BED merge_%s.bed -p 20 --bamfiles %s --label %s --outRawCounts %s_counts.xls -o %s_counts.npz"%(MyPrefix," ".join(BamList)," ".join(SampleList),MyPrefix,MyPrefix))#need change and add peak name
		os.system("sed -i \"1s/[#|']//g\" %s_counts.xls"%(MyPrefix))
		df=pandas.read_table("%s_counts.xls"%MyPrefix,sep="\t",header=0,index_col=[0,1,2])
		df_bed=pandas.read_table("merge_%s.bed"%MyPrefix,sep="\t",header=None,index_col=[0,1,2],names=["chr","start","end","peak","score","strand"])
		df_r=pandas.concat([df,df_bed],axis=1)
		labels=["peak"]+SampleList
		df_r.loc[:,labels].to_csv("%s_counts_final.xls"%MyPrefix,sep="\t",na_rep="Nan",index=False)
	def GetSenseAntisenseMatrix(self,BedFile,BamDic,MyPrefix=""):#BamDic={"fwd":{"sample":"sample_fwd.bam"},"rev":{"sample":"sample_rev.bam"}}
		if not MyPrefix:
			MyPrefix=self.Prefix
		Sense=pandas.DataFrame()
		Antisense=pandas.DataFrame()
		df_bed=pandas.read_table(BedFile,sep="\t",header=None,index_col=[0,1,2],names=["chr","start","end","gene","score","strand"])
		for s in ["fwd","rev"]:
			BamList=[]
			SampleList=[]
			for si in BamDic[s]:
				SampleList.append(si)
				BamList.append(BamDic[s][si])
			os.system("multiBamSummary BED-file --BED %s -p 20 --bamfiles %s --label %s --outRawCounts %s_%s_counts.xls -o %s_%s.npz"%(BedFile," ".join(BamList)," ".join(SampleList),MyPrefix,s,MyPrefix,s))
			os.system("sed \"1s/[#|']//g\" %s_%s_counts.xls >%s_%s_counts_deal.xls"%(MyPrefix,s,MyPrefix,s))
			df=pandas.read_table("%s_%s_counts_deal.xls"%(MyPrefix,s),sep="\t",header=0,index_col=[0,1,2])
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
	def GetRandomBed(self,):
		pass
	def GetGCSkewBwAndGz(self,):
		pass
	def GetAnnoPeak(self,):
		pass
	def GetNoiseqFile(self,):
		pass
	def GetDeseq2File(self,):
		pass
	def GetCluster(self,): #mfuzz tsne
		pass