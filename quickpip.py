#!/usr/bin/env python2.7
import sys,os
class BasePipLine(object):
	"""docstring for DripPipLine"""
	def __init__(self, Prefix):
		self.Prefix = Prefix
	def Align(self,Left,Right,GenomeIndex):
		os.system("bowtie2 --local --phred33 -p 15 -t -x %s -1 %s -2 %s 2>%s_align.info|samtools view -bS -1 |samtools sort -@ 15 -m 5G -l 9 -o %s.sort.bam"%(GenomeIndex,Left,Right,self.Prefix,self.Prefix))
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
	def FormatPeak(self,FilterChromFile,PeakExcelFile=self.Prefix+"_peaks.xls",PeakBedFile=self.Prefix+"_peaks.bed"):
		AwkRegular=""
		for x in open(FilterChromFile):
			x=x.rstrip()
			AwkRegular+="^"+x+"$"+"|"
		AwkRegular=AwkRegular[:-1]
		os.system("awk -F'\\t' '$0!~/^#/&&$0!=\"\"&&$1!~/%s/{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$NF}' %s  > %s"%(AwkRegular,PeakExcelFile,PeakBedFile))
	def GetNormBw(self,BinSize,GenomeSize,IgnoreFile,InputBam=self.Prefix+".sort.paird_dup.bam",InputFwdBam=self.Prefix+"_fwd.bam",InputRevBam=self.Prefix+"_rev.bam",OutputBam=self.Prefix+"_nucleus_norm.bw",OutputFwdBam=self.Prefix+"_fwd_nucleus_norm.bw",OutputRevBam=self.Prefix+"_rev_nucleus_norm.bw"):
		IgnoreList=[x.rstrip() for x in open(IgnoreFile)]
		Lines=os.popen("bamCoverage --extendReads -v -p 20 -b %s -o %s --binSize %s --effectiveGenomeSize %s --normalizeUsing RPGC --ignoreForNormalization %s"%(InputBam,OutputBam,BinSize,GenomeSize," ".join(IgnoreList))).readlines()
		for x in Lines:
			x=x.rstrip()
			if x.startswith("Final scaling factor:"):
				l=x.split()
				Scale=l[-1].strip()
				break
		os.system("bamCoverage --extendReads --scaleFactor %s -v -p 20 -b %s -o %s --binSize %s --ignoreForNormalization %s"%(Scale,InputFwdBam,OutputFwdBam,BinSize," ".join(IgnoreList)))
		os.system("bamCoverage --extendReads --scaleFactor %s -v -p 20 -b %s -o %s --binSize %s --ignoreForNormalization %s"%(Scale,InputRevBam,OutputRevBam,BinSize," ".join(IgnoreList)))
		Fr=open(self.Prefix+"_scale.txt","w")
		Fr.write(self.Prefix+"\t"+Scale+"\n")
		Fr.close()
class AnalysisPipLine(object):
	"""heheh"""
	def __init__(self,args):
		self,args=args
	def BasicStatistics(self,):
		pass
	def BwCorrelation(self,):
		pass
	def FindMotif(self,):
		pass
	def PeakContentDistribution(self,): #genebody tss tts etc....
		pass
	def PeakLengthDistribution(self,):
		pass
	def Metaplot(self,): #gene sense antisense 
		pass
	def BedCorrelation(self,): #permutation test,metaplot
		pass
