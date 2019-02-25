#!/usr/bin/env python2.7
import json,sys,os
from QuickPip import *
def CheckJson(ConfigDic):
	pass
def BaseAnalysis(ConfigDic):
	Path=os.getcwd()
	for x in open(ConfigDic["Target"]):#target:group\tsample\tleft\tright
		x=x.rstrip()
		Group,Sample,Left,Right=x.split("\t")
		os.mkdir(Sample)
		os.chdir(Path+"/"+Sample)
		p=BasePipLine(Sample,ConfigDic["Thread"])
		p.Align(Left,Right,ConfigDic["GeomeIndex"])
		p.ReomveDuplication()
		p.SplitStrand()
		p.CallPeak(ConfigDic["GenomeSize"])
		p.CallPeak(ConfigDic["GenomeSize"],MyPrefix="%s_fwd"%Sample,InputBamFile="%s_fwd.bam"%Sample)
		p.CallPeak(ConfigDic["GenomeSize"],MyPrefix="%s_rev"%Sample,InputBamFile="%s_rev.bam"%Sample)
		p.FormatPeak(ConfigDic["FilterChromFile"],Strand=".")
		p.FormatPeak(ConfigDic["FilterChromFile"],PeakExcelFile="%s_fwd_peaks.xls"%Sample,PeakBedFile="%s_fwd_peaks.bed"%Sample,Strand="+")
		p.FormatPeak(ConfigDic["FilterChromFile"],PeakExcelFile="%s_rev_peaks.xls"%Sample,PeakBedFile="%s_rev_peaks.bed"%Sample,Strand="-")
		p.GetNormBw(ConfigDic["BinSize"],ConfigDic["GenomeSize"],ConfigDic["FilterChromFile"])
		os.chdir(Path)
def DownstreamAnalysis(ConfigDic):
	Path=os.getcwd()
	SampleDic={}
	SampleList=[]
	BwFileList=[]
	RevBwFileList=[]
	FwdBwFileList=[]
	for x in open(ConfigDic["Target"]):
		x=x.rstrip()
		Group,Sample,Left,Right=x.split("\t")
		SampleDic[Sample]=["%s/%s/%s_align.info"%(Path,Sample,Sample),"%s/%s/%s.matrix"%(Path,Sample,Sample),"%s/%s/%s_peaks.xls"%(Path,Sample,Sample),"%s/%s/%s_peaks.bed"%(Path,Sample,Sample),"%s/%s/%s_fwd_peaks.bed"%(Path,Sample,Sample),"%s/%s/%s_rev_peaks.bed"%(Path,Sample,Sample)]
		SampleList.append(Sample)
		BwFileList.append("%s/%s/%s_nucleus_norm.bw"%(Path,Sample,Sample))
		FwdBwFileList.append("%s/%s/%s_fwd_nucleus_norm.bw"%(Path,Sample,Sample))
		RevBwFileList.append("%s/%s/%s_rev_nucleus_norm.bw"%(Path,Sample,Sample))
	p=AnalysisPipLine(ConfigDic["ProjectName"],ConfigDic["Thread"])
	p.BasicStatistics(SampleDic)
	p.BwCorrelation(BwFileList,ConfigDic["BinSize"],ConfigDic["FilterChromFile"],ConfigDic["CorrelationMethod"],"heatmap")
	p.BwCorrelation(BwFileList,ConfigDic["BinSize"],ConfigDic["FilterChromFile"],ConfigDic["CorrelationMethod"],"scatterplot",OnlyPlot=True)
	p.BwCorrelation(FwdBwFileList,ConfigDic["BinSize"],ConfigDic["FilterChromFile"],ConfigDic["CorrelationMethod"],"heatmap",MyPrefix="%s_fwd"%ConfigDic["ProjectName"]) #fwd
	p.BwCorrelation(FwdBwFileList,ConfigDic["BinSize"],ConfigDic["FilterChromFile"],ConfigDic["CorrelationMethod"],"scatterplot",MyPrefix="%s_fwd"%ConfigDic["ProjectName"],OnlyPlot=True) #fwd
	p.BwCorrelation(RevBwFileList,ConfigDic["BinSize"],ConfigDic["FilterChromFile"],ConfigDic["CorrelationMethod"],"heatmap",MyPrefix="%s_rev"%ConfigDic["ProjectName"]) #rev
	p.BwCorrelation(RevBwFileList,ConfigDic["BinSize"],ConfigDic["FilterChromFile"],ConfigDic["CorrelationMethod"],"scatterplot",MyPrefix="%s_rev"%ConfigDic["ProjectName"],OnlyPlot=True) #rev
	p.FindMotif(ConfigDic["ProjectName"]+"_fwd",ConfigDic["GenomeFastaFile"],"%s/%s/%s_fwd_peaks.bed"%(Path,ConfigDic["ControlSample"],ConfigDic["ControlSample"]),ConfigDic["ChromSize"],Strand="+",RepeatNum=ConfigDic["RepeatNum"])
	p.FindMotif(ConfigDic["ProjectName"]+"_rev",ConfigDic["GenomeFastaFile"],"%s/%s/%s_rev_peaks.bed"%(Path,ConfigDic["ControlSample"],ConfigDic["ControlSample"]),ConfigDic["ChromSize"],Strand="-",RepeatNum=ConfigDic["RepeatNum"])
	p.GetGenomeContentBedFile(ConfigDic["GeneBed"],ConfigDic["Extend"],ConfigDic["ChromSize"])
	PeakContentDic={}
	for s in SampleList:
		ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s/%s/%s_peaks.bed"%(Path,s,s),"tss_%s.bed"%ConfigDic["Extend"],"tts_%s.bed"%ConfigDic["Extend"],"genebody_%s.bed"%ConfigDic["Extend"])
		PeakContentDic[s]=[ProNum,TerNum,BodyNum,InterNum]
		p.GetRandomBed(s,"%s/%s/%s_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"])
		ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s_random.bed"%(s),"tss_%s.bed"%ConfigDic["Extend"],"tts_%s.bed"%ConfigDic["Extend"],"genebody_%s.bed"%ConfigDic["Extend"])
		PeakContentDic[s+"_random"]=[ProNum,TerNum,BodyNum,InterNum]
		p.PeakContentDistribution("%s/%s/%s_fwd_peaks.bed"%(Path,s,s),"tss_%s.bed"%ConfigDic["Extend"],"tts_%s.bed"%ConfigDic["Extend"],"genebody_%s.bed"%ConfigDic["Extend"])
		PeakContentDic[s+"_fwd"]=[ProNum,TerNum,BodyNum,InterNum]
		p.GetRandomBed(s+"_fwd","%s/%s/%s_fwd_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"])
		ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s_fwd_random.bed"%(s),"tss_%s.bed"%ConfigDic["Extend"],"tts_%s.bed"%ConfigDic["Extend"],"genebody_%s.bed"%ConfigDic["Extend"])
		PeakContentDic[s+"_fwd_random"]=[ProNum,TerNum,BodyNum,InterNum]
		ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s/%s/%s_rev_peaks.bed"%(Path,s,s),"tss_%s.bed"%ConfigDic["Extend"],"tts_%s.bed"%ConfigDic["Extend"],"genebody_%s.bed"%ConfigDic["Extend"])
		PeakContentDic[s+"_rev"]=[ProNum,TerNum,BodyNum,InterNum]
		p.GetRandomBed(s+"_rev","%s/%s/%s_rev_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"])
		ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s_rev_random.bed"%(s),"tss_%s.bed"%ConfigDic["Extend"],"tts_%s.bed"%ConfigDic["Extend"],"genebody_%s.bed"%ConfigDic["Extend"])
		PeakContentDic[s+"_rev_random"]=[ProNum,TerNum,BodyNum,InterNum]
	df=pandas.DataFrame.from_dict(data=PeakContentDic,orient='index')
	df.to_csv("%s_peaks_content_distribution.xls"%ConfigDic["ProjectName"],sep="\t",header=["ProNum","TerNum","BodyNum","InterNum"],index_label="Sample")
	p.PeakLengthDistribution()
def Main(ConfigFile,SubCommand):
	ConfigDic=json.load(open(ConfigFile))
	if SubCommand=="BaseAnalysis":
		BaseAnalysis(ConfigDic)
	elif SubCommand=="DownstreamAnalysis":
		DownstreamAnalysis(ConfigDic)
	elif SubCommand=="AllPip":
		BaseAnalysis(ConfigDic)
		DownstreamAnalysis(ConfigDic)
if __name__ == '__main__':
	ConfigFile=sys.argv[1]
	SubCommand=sys.argv[2] #BaseAnalysis ,DownstreamAnalysis,AllPip
	Main(ConfigFile)
