#!/usr/bin/env python
import json,sys,os,warnings
from QuickPip import *
def CheckJson(ConfigDic):
	pass
def BaseAnalysis(ConfigDic):
	Path=os.getcwd()
	SampleDic={}
	FrScale=open(ConfigDic["ProjectName"]+"_bwscale.xls","w")
	for i,x in enumerate(open(ConfigDic["Target"])):#target:group\tsample\tleft\tright\tinput\tleft\tright
		x=x.rstrip()
		l=x.split("\t")
		if len(l)==4:
			Group,Sample,Left,Right=l
			if i==0:
				FrScale.write("sample\t1x_scale\n")
		elif len(l)==7:
			Group,Sample,Left,Right,Input,InputLeft,InputRight=l
			if i==0:
				FrScale.write("sample\tinput\t1x_scale\tspike_scale\tbw_finalscale\n")
			SampleDic[Input]=["%s/%s/%s_align.info"%(Path,Sample,Input),"%s/%s/%s.matrix"%(Path,Sample,Input)]
		else:
			print("target error")
			sys.exit()
		SampleDic[Sample]=["%s/%s/%s_align.info"%(Path,Sample,Sample),"%s/%s/%s.matrix"%(Path,Sample,Sample),"%s/%s/%s_peaks.xls"%(Path,Sample,Sample),"%s/%s/%s_peaks.bed"%(Path,Sample,Sample),"%s/%s/%s_fwd_peaks.bed"%(Path,Sample,Sample),"%s/%s/%s_rev_peaks.bed"%(Path,Sample,Sample)]
		if not os.path.exists(Sample):
			os.mkdir(Sample)
		else:
			warnings.warn("%s folder already exists"%Sample)
		os.chdir(Path+"/"+Sample)
		p=BasePipLine(Sample,ConfigDic["Thread"])
		p.Align(ConfigDic["SeqDataPath"]+"/"+Left,ConfigDic["SeqDataPath"]+"/"+Right,ConfigDic["GeomeIndex"])
		p.ReomveDuplication()
		p.SplitStrand()
		p.CallPeak(ConfigDic["GenomeSize"])
		p.CallPeak(ConfigDic["GenomeSize"],MyPrefix="%s_fwd"%Sample,InputBamFile="%s_fwd.bam"%Sample)
		p.CallPeak(ConfigDic["GenomeSize"],MyPrefix="%s_rev"%Sample,InputBamFile="%s_rev.bam"%Sample)
		p.FormatPeak(ConfigDic["FilterChromFile"],Strand=".")
		p.FormatPeak(ConfigDic["FilterChromFile"],PeakExcelFile="%s_fwd_peaks.xls"%Sample,PeakBedFile="%s_fwd_peaks.bed"%Sample,Strand="+")
		p.FormatPeak(ConfigDic["FilterChromFile"],PeakExcelFile="%s_rev_peaks.xls"%Sample,PeakBedFile="%s_rev_peaks.bed"%Sample,Strand="-")
		if len(l)==7:
			p2=BasePipLine(Input,ConfigDic["Thread"])
			p2.Align(ConfigDic["SeqDataPath"]+"/"+InputLeft,ConfigDic["SeqDataPath"]+"/"+InputRight,ConfigDic["GeomeIndex"])
			p2.ReomveDuplication()
			DripBam=Sample+".sort.paird_dup.bam"
			InputBam=Input+".sort.paird_dup.bam"
			SpikeNameList=ConfigDic["SpikeNameList"]
			SpikeScale=p2.GetSpikeScale(DripBam,InputBam,SpikeNameList)
			Scale1x=p.Get1XNormBw(ConfigDic["BinSize"],ConfigDic["GenomeSize"],ConfigDic["FilterChromFile"],SpikeScale=SpikeScale)
			FinalScale=Scale1x*SpikeScale
			Fm=Sample+"\t"+Input+"\t"+str(Scale1x)+"\t"+str(SpikeScale)+"\t"+str(FinalScale)
			FrScale.write(Fm+"\n")
			p2.AfterMath()
		else:
			Scale1x=p.GetRandomRegionNormBw(ConfigDic["FilterChromFile"],ConfigDic["Extend"],ConfigDic["GenomeSize"],ConfigDic["RepeatNum"])
			#Scale1x=p.Get1XNormBw(ConfigDic["BinSize"],ConfigDic["GenomeSize"],ConfigDic["FilterChromFile"])
			FrScale.write(Sample+"\t"+str(Scale1x)+"\n")
		p.AfterMath()
		os.chdir(Path)
	pa=AnalysisPipLine(ConfigDic["ProjectName"],ConfigDic["Thread"])
	pa.BasicStatistics(SampleDic)
	FrScale.close()
def DeseqAnalysis(ConfigDic,FilePath="",SpikeScaleDic={}): #FilePath include peak bam file
	PeakDic={}
	PeakFwdDic={}
	PeakRevDic={}
	BamDic={}
	BamFwdDic={}
	BamRevDic={}
	DeseqScaleDic={}
	SampleInputDic={}
	Path=os.getcwd()
	if not FilePath:
		for x in open(ConfigDic["Target"]):
			x=x.rstrip()
			l=x.split("\t")
			PeakDic[l[1]]="%s/%s/%s_peaks.bed"%(Path,l[1],l[1])
			PeakFwdDic[l[1]]="%s/%s/%s_fwd_peaks.bed"%(Path,l[1],l[1])
			PeakRevDic[l[1]]="%s/%s/%s_rev_peaks.bed"%(Path,l[1],l[1])
			BamDic[l[1]]="%s/%s/%s.sort.paird_dup.bam"%(Path,l[1],l[1])
			BamFwdDic[l[1]]="%s/%s/%s_fwd.bam"%(Path,l[1],l[1])
			BamRevDic[l[1]]="%s/%s/%s_rev.bam"%(Path,l[1],l[1])
			if len(l)==7:
				SampleInputDic[l[1]]=l[4]
	else:		
		for x in open(ConfigDic["Target"]):
			x=x.rstrip()
			l=x.split("\t")
			PeakDic[l[1]]="%s/%s_peaks.bed"%(FilePath,l[1])
			PeakFwdDic[l[1]]="%s/%s_fwd_peaks.bed"%(FilePath,l[1])
			PeakRevDic[l[1]]="%s/%s_rev_peaks.bed"%(FilePath,l[1])
			BamDic[l[1]]="%s/%s.sort.paird_dup.bam"%(FilePath,l[1])
			BamFwdDic[l[1]]="%s/%s_fwd.bam"%(FilePath,l[1])
			BamRevDic[l[1]]="%s/%s_rev.bam"%(FilePath,l[1])
			if len(l)==7:
				SampleInputDic[l[1]]=l[4]
	pa=AnalysisPipLine(ConfigDic["ProjectName"],ConfigDic["Thread"])
	os.mkdir("deseq")
	os.mkdir("deseq/all")
	os.chdir(Path+"/deseq/all")
	pa.GetMergePeakCountMatrix(PeakDic,BamDic)
	CountMatrixFile="%s_counts_final.xls"%ConfigDic["ProjectName"]
	if SampleInputDic:
		for x in open("%s_deseq_scale.txt"%ConfigDic["ProjectName"]).readlines()[1:]:
			x=x.rstrip()
			li=x.split("\t")
			DeseqScaleDic[li[0]]=float(li[1])
		if not SpikeScaleDic:
			for x in open(Path+"/"+ConfigDic["ProjectName"]+"_bwscale.xls").readlines()[1:]:
				x=x.rstrip()
				li=x.split("\t")
				SpikeScaleDic[li[0]]=float(li[3])
	pa.GetDeseq2File(Path+"/"+ConfigDic["Target"],CountMatrixFile,SpikeScaleDic,DeseqScaleDic,Control=ConfigDic["ControlSample"])
	pa.GetAnnoPeak("merge_%s.bed"%ConfigDic["ProjectName"],ConfigDic["GeneBed"],glob.glob("*.xls"))
	os.chdir(Path)
	os.mkdir("deseq/fwd")
	os.chdir(Path+"/deseq/fwd")
	pa.GetMergePeakCountMatrix(PeakFwdDic,BamFwdDic,MyPrefix=ConfigDic["ProjectName"]+"_fwd")
	pa.GetDeseq2File(Path+"/"+ConfigDic["Target"],"%s_fwd_counts_final.xls"%ConfigDic["ProjectName"],SpikeScaleDic,DeseqScaleDic,MyPrefix=ConfigDic["ProjectName"]+"_fwd")
	pa.GetAnnoPeak("merge_%s_fwd.bed"%ConfigDic["ProjectName"],ConfigDic["GeneBed"],glob.glob("*.xls"))
	os.chdir(Path)
	os.mkdir("deseq/rev")
	os.chdir(Path+"/deseq/rev")
	pa.GetMergePeakCountMatrix(PeakRevDic,BamRevDic,MyPrefix=ConfigDic["ProjectName"]+"_rev")
	pa.GetDeseq2File(Path+"/"+ConfigDic["Target"],"%s_rev_counts_final.xls"%ConfigDic["ProjectName"],SpikeScaleDic,DeseqScaleDic,MyPrefix=ConfigDic["ProjectName"]+"_rev")
	pa.GetAnnoPeak("merge_%s_rev.bed"%ConfigDic["ProjectName"],ConfigDic["GeneBed"],glob.glob("*.xls"))
	os.chdir(Path)
	if SampleInputDic:
		Fr=open("%s_final_deseq_scale.txt"%ConfigDic["ProjectName"],"w")
		Fr.write("sample\tinput\tspike_scale\tdeseq(deeptools)\tcustom_deseq(deseq)\n")
		for s in SpikeScaleDic:
			fm=s+"\t"+SampleInputDic[s]+"\t"+str(SpikeScaleDic[s])+"\t"+str(DeseqScaleDic[s])+"\t"+str(1/(SpikeScaleDic[s]*DeseqScaleDic[s]))
			Fr.write(fm+"\n")
		Fr.close()
def DownstreamAnalysis(ConfigDic):
	Path=os.getcwd()
	#SampleDic={}
	SampleList=[]
	BwFileList=[]
	RevBwFileList=[]
	FwdBwFileList=[]
	for x in open(ConfigDic["Target"]):
		x=x.rstrip()
		Group,Sample,Left,Right=x.split("\t")[0:4]
		#SampleDic[Sample]=["%s/%s/%s_align.info"%(Path,Sample,Sample),"%s/%s/%s.matrix"%(Path,Sample,Sample),"%s/%s/%s_peaks.xls"%(Path,Sample,Sample),"%s/%s/%s_peaks.bed"%(Path,Sample,Sample),"%s/%s/%s_fwd_peaks.bed"%(Path,Sample,Sample),"%s/%s/%s_rev_peaks.bed"%(Path,Sample,Sample)]
		SampleList.append(Sample)
		BwFileList.append("%s/%s/%s_nucleus_norm.bw"%(Path,Sample,Sample))
		FwdBwFileList.append("%s/%s/%s_fwd_nucleus_norm.bw"%(Path,Sample,Sample))
		RevBwFileList.append("%s/%s/%s_rev_nucleus_norm.bw"%(Path,Sample,Sample))
	p=AnalysisPipLine(ConfigDic["ProjectName"],ConfigDic["Thread"])
	os.mkdir("%s_analysis"%ConfigDic["ProjectName"])
	os.mkdir("%s_analysis/correlation"%ConfigDic["ProjectName"])
	os.chdir(Path+"/%s_analysis/correlation/"%ConfigDic["ProjectName"])
	#p.BasicStatistics(SampleDic)
	BwDic={"all":BwFileList,"fwd":FwdBwFileList,"rev":RevBwFileList}
	for m in ["spearman","pearson"]:
		for st in ["fwd","rev","all"]:
			MyPrefix="%s_%s"%(ConfigDic["ProjectName"],st)
			p.BwCorrelation(BwDic[st],ConfigDic["BinSize"],ConfigDic["FilterChromFile"],m,"heatmap",MyPrefix)
			p.BwCorrelation(BwDic[st],ConfigDic["BinSize"],ConfigDic["FilterChromFile"],m,"scatterplot",MyPrefix,OnlyPlot=True)
	os.chdir(Path)
	os.mkdir("%s_analysis/motif"%ConfigDic["ProjectName"])
	os.mkdir("%s_analysis/sense_antisense"%ConfigDic["ProjectName"])
	p.GetGenomeContentBedFile(ConfigDic["GeneBed"],ConfigDic["Extend"],ConfigDic["ChromSize"])
	os.mkdir("%s_analysis/skew"%ConfigDic["ProjectName"])
	os.chdir(Path+"/%s_analysis/skew/"%ConfigDic["ProjectName"])
	p.GetATGCSkewBw(ConfigDic["GenomeFastaFile"],ConfigDic["ChromSize"],ConfigDic["Win"],ConfigDic["Step"])
	GCSkewBw="%s_GCSkew.bw"%ConfigDic["ProjectName"]
	ATSkewBw="%s_ATSkew.bw"%ConfigDic["ProjectName"]
	p.GetSkewGz(GCSkewBw,ATSkewBw,ConfigDic["GeneBed"],ConfigDic["Extend"],"gene")
	os.chdir(Path)
	PeakContentDic={}
	PeakLengthDic={}
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
		PeakLengthDic[s]=p.PeakLengthDistribution("%s/%s/%s_peaks.bed"%(Path,s,s))
		PeakLengthDic[s+"_fwd"]=p.PeakLengthDistribution("%s/%s/%s_fwd_peaks.bed"%(Path,s,s))
		PeakLengthDic[s+"_rev"]=p.PeakLengthDistribution("%s/%s/%s_rev_peaks.bed"%(Path,s,s))
		FwdBw="%s/%s/%s_fwd_nucleus_norm.bw"%(Path,s,s)
		RevBw="%s/%s/%s_rev_nucleus_norm.bw"%(Path,s,s)
		os.chdir(Path+"/%s_analysis/sense_antisense/"%ConfigDic["ProjectName"])
		p.SenseAntisense(s,ConfigDic["GeneBed"],FwdBw,RevBw,ConfigDic["Extend"])
		os.chdir(Path+"/%s_analysis/skew/"%ConfigDic["ProjectName"])
		p.GetSkewGz(GCSkewBw,ATSkewBw,"%s/%s/%s_peaks.bed"%(Path,s,s),ConfigDic["Extend"],s)
		p.GetSkewGz(GCSkewBw,ATSkewBw,"%s/%s/%s_fwd_peaks.bed"%(Path,s,s),ConfigDic["Extend"],s+"_fwd")
		p.GetSkewGz(GCSkewBw,ATSkewBw,"%s/%s/%s_rev_peaks.bed"%(Path,s,s),ConfigDic["Extend"],s+"_rev")
		os.chdir(Path+"/%s_analysis/motif/"%ConfigDic["ProjectName"])
		p.FindMotif(ConfigDic["GenomeFastaFile"],"%s/%s/%s_fwd_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"],MyPrefix=ConfigDic["ProjectName"]+"_fwd",Strand="+",RepeatNum=ConfigDic["RepeatNum"])
		p.FindMotif(ConfigDic["GenomeFastaFile"],"%s/%s/%s_rev_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"],MyPrefix=ConfigDic["ProjectName"]+"_rev",Strand="-",RepeatNum=ConfigDic["RepeatNum"])
		os.chdir(Path)
	df=pandas.DataFrame.from_dict(data=PeakContentDic,orient='index')
	os.mkdir("%s_analysis/peaks_content_distribution"%ConfigDic["ProjectName"])
	os.chdir(Path+"/%s_analysis/peaks_content_distribution/"%ConfigDic["ProjectName"])
	df.to_csv("%s_peaks_content_distribution.xls"%ConfigDic["ProjectName"],sep="\t",header=["ProNum","TerNum","BodyNum","InterNum"],index_label="Sample")
	os.chdir(Path)
	os.mkdir("%s_analysis/peaks_length_distribution"%ConfigDic["ProjectName"])
	os.chdir(Path+"/%s_analysis/peaks_length_distribution/"%ConfigDic["ProjectName"])
	np.save('%s_peaks_length_distribution.npy'%ConfigDic['ProjectName'],PeakLengthDic)
	os.chdir(Path)
	os.mkdir("%s_analysis/cluster"%ConfigDic["ProjectName"])
	os.chdir(Path+"/%s_analysis/cluster/"%ConfigDic["ProjectName"])
	for st in ["fwd","rev","all"]:
		if st=="all":
			LevelFile=Path+"/deseq/%s/%s_counts_final.xls"%(st,ConfigDic["ProjectName"])
		else:
			LevelFile=Path+"/deseq/%s/%s_%s_counts_final.xls"%(st,ConfigDic["ProjectName"],st)
		p.GetDiffRloopLevelFile(glob.glob(Path+"/deseq/%s/*_diffexpr_results.xls"%st),ConfigDic["Target"],LevelFile,MyPrefix=ConfigDic["ProjectName"]+"_"+st)
		p.GetCluster("%s_%s_diff_union.xls"%(ConfigDic["ProjectName"],st),MyPrefix=ConfigDic["ProjectName"]+"_"+st)
	os.chdir(Path)
	#p.BedCorrelation()
	p.AfterMath()
def Main(ConfigFile,SubCommand):
	ConfigDic=json.load(open(ConfigFile))
	if SubCommand=="BaseAnalysis":
		BaseAnalysis(ConfigDic)
	elif SubCommand=="DeseqAnalysis":
		DeseqAnalysis(ConfigDic)
	elif SubCommand=="DownstreamAnalysis":
		DownstreamAnalysis(ConfigDic)
	elif SubCommand=="AllPip":
		BaseAnalysis(ConfigDic)
		DeseqAnalysis(ConfigDic)
		DownstreamAnalysis(ConfigDic)
if __name__ == '__main__':
	ConfigFile=sys.argv[1]
	SubCommand=sys.argv[2] #BaseAnalysis ,DownstreamAnalysis,AllPip
	Main(ConfigFile,SubCommand)