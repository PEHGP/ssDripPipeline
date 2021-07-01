#!/usr/bin/env python
import json,sys,os
from QuickPip import *
import logging,datetime,glob
logging.basicConfig(level=logging.INFO,format = '%(asctime)s - %(pathname)s - %(filename)s - %(funcName)s - %(levelname)s - %(message)s',filemode="w",filename="log_%s.txt"%datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
def CheckJson(ConfigDic):
	Cl2=["SeqDataPath","FilterChromFile","ChromSize","GenomeFastaFile","GeneBeb"]
	if not os.path.exists(ConfigDic["Target"]):
		M="Target file is not exists."
		logging.info(M)
		print(M)
		sys.exit()
	if len(glob.glob(ConfigDic["GenomeIndex"]+".*"))<6:
		M="There is no bowtie2-buid index in the %s, or the index is incomplete."%ConfigDic["GenomeIndex"]
		logging.info(M)
		print(M)
		sys.exit()
	for c in Cl2:
		if not ConfigDic[c].startswith("/"):
			M="%s does not contain an absolute path."%c
			logging.info(M)
			print(M)
			sys.exit()
		if not os.path.exists(ConfigDic[c]):
			M="%s file is not exists."
			logging.info(M)
			sys.exit()
def CheckContinue(ConfigDic): #Add new analysis, need to change here
	ExistsDic={"BaseAnalysis":[],"DeseqAnalysis":[],"DownstreamAnalysis":[]}
	for i,x in enumerate(open(ConfigDic["Target"])):
		x=x.rstrip()
		l=x.split("\t")
		Sample=l[1]
		if os.path.exists(Sample+"/"+Sample+"_scale.xls"):
			if len(open(Sample+"/"+Sample+"_scale.xls").readlines())==2:
				ExistsDic["BaseAnalysis"].append(Sample)
	for x in ["all","fwd","rev"]:
		if os.path.exists("%s_deseq/%s/%s_deseq_scale.txt"%(ConfigDic["ProjectName"],x,ConfigDic["ProjectName"])):
			ExistsDic["DeseqAnalysis"].append(x)
	for x in ["correlation","motif","sense_antisense","skew","peaks_content_distribution","peaks_length_distribution","cluster"]:
		if os.path.exists("%s_analysis/%s"%(ConfigDic["ProjectName"],x)):
			ExistsDic["DownstreamAnalysis"].append(x)
	return ExistsDic
def BaseAnalysis(ConfigDic):
	ExistsDic=CheckContinue(ConfigDic)
	Path=os.getcwd()
	SampleDic={}
	if ExistsDic["BaseAnalysis"]!=[]:
		logging.info("Some of the output results of BaseAnalysis already exist, and the program continues to run in continue mode.")
	for i,x in enumerate(open(ConfigDic["Target"])):#target:group\tsample\tleft\tright\tinput\tleft\tright
		x=x.rstrip()
		l=x.split("\t")
		if len(l)==4:
			Group,Sample,Left,Right=l
		elif len(l)==7:
			Group,Sample,Left,Right,Input,InputLeft,InputRight=l
			SampleDic[Input]=["%s/%s/%s_align.info"%(Path,Sample,Input),"%s/%s/%s.matrix"%(Path,Sample,Input)]
		else:
			print("target error")
			logging.error("target error")
			sys.exit()
		SampleDic[Sample]=["%s/%s/%s_align.info"%(Path,Sample,Sample),"%s/%s/%s.matrix"%(Path,Sample,Sample),"%s/%s/%s_peaks.xls"%(Path,Sample,Sample),"%s/%s/%s_peaks.bed"%(Path,Sample,Sample),"%s/%s/%s_fwd_peaks.bed"%(Path,Sample,Sample),"%s/%s/%s_rev_peaks.bed"%(Path,Sample,Sample),Sample+"_scale.xls"]
		if Sample in ExistsDic["BaseAnalysis"]:
			logging.info("%s folder already exists. This sample will not be analyzed. If you want to analyze this sample, please delete the folder of this sample."%Sample)
		else:
			os.mkdir(Sample)
			os.chdir(Path+"/"+Sample)
			FrScale=open(Sample+"_scale.xls","w")
			p=BasePipLine(Sample,ConfigDic["Thread"])
			p.Align(ConfigDic["SeqDataPath"]+"/"+Left,ConfigDic["SeqDataPath"]+"/"+Right,ConfigDic["GenomeIndex"])
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
				p2.Align(ConfigDic["SeqDataPath"]+"/"+InputLeft,ConfigDic["SeqDataPath"]+"/"+InputRight,ConfigDic["GenomeIndex"])
				p2.ReomveDuplication()
				DripBam=Sample+".sort.paird_dup.bam"
				InputBam=Input+".sort.paird_dup.bam"
				SpikeNameList=ConfigDic["SpikeNameList"]
				SpikeScale=p2.GetSpikeScale(DripBam,InputBam,SpikeNameList)
				Scale1x=p.Get1XNormBw(ConfigDic["BinSize"],ConfigDic["GenomeSize"],ConfigDic["FilterChromFile"],SpikeScale=SpikeScale)
				FinalScale=Scale1x*SpikeScale
				FrScale.write("sample\tinput\t1x_scale\tspike_scale\tbw_finalscale\n")
				Fm=Sample+"\t"+Input+"\t"+str(Scale1x)+"\t"+str(SpikeScale)+"\t"+str(FinalScale)
				FrScale.write(Fm+"\n")
				p2.AfterMath()
			else:
				Scale1x=p.GetRandomRegionNormBw(ConfigDic["FilterChromFile"],ConfigDic["MetaplotExtend"],ConfigDic["GenomeSize"],ConfigDic["RepeatNum"])
				#Scale1x=p.Get1XNormBw(ConfigDic["BinSize"],ConfigDic["GenomeSize"],ConfigDic["FilterChromFile"])
				FrScale.write("sample\t1x_scale\n")
				FrScale.write(Sample+"\t"+str(Scale1x)+"\n")
			FrScale.close()
			p.AfterMath()
		os.chdir(Path)
	pa=AnalysisPipLine(ConfigDic["ProjectName"],ConfigDic["Thread"])
	pa.BasicStatistics(SampleDic)
def DeseqAnalysis(ConfigDic,FilePath="",SpikeScaleDic={}): #FilePath include peak bam file
	ExistsDic=CheckContinue(ConfigDic)
	PeakDic={}
	PeakFwdDic={}
	PeakRevDic={}
	BamDic={}
	BamFwdDic={}
	BamRevDic={}
	DeseqScaleDic={}
	SampleInputDic={}
	SampleScaleDic={}
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
			SampleScaleDic[l[1]]="%s/%s/%s_scale.xls"%(Path,l[1],l[1])
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
			SampleScaleDic[l[1]]="%s/%s_scale.xls"%(FilePath,l[1])
			if len(l)==7:
				SampleInputDic[l[1]]=l[4]
	Flag=0
	for x in open(ConfigDic["Target"]):
		x=x.rstrip()
		l=x.split("\t")
		for di in [PeakDic,PeakFwdDic,PeakRevDic,BamDic,BamFwdDic,BamRevDic,SampleScaleDic]:
			if not os.path.exists(di[l[1]]):
				Flag=1
				logging.info("%s not exists."%di[l[1]])
	if Flag==1:
		print("Some files could not be found. The program exited. Please check the log")
		sys.exit()
	pa=AnalysisPipLine(ConfigDic["ProjectName"],ConfigDic["Thread"])
	if ExistsDic["DeseqAnalysis"]==[]:
		os.mkdir("%s_deseq"%ConfigDic["ProjectName"])
	if "all" in ExistsDic["DeseqAnalysis"]:
		logging.info("\"all\" folder already exists. Deseq2 will not be run. If you want to do this analysis, please delete this folder.")
		if SampleInputDic:
			for x in open("%s_deseq/all/%s_deseq_scale.txt"%(ConfigDic["ProjectName"],ConfigDic["ProjectName"])).readlines()[1:]:
				x=x.rstrip()
				li=x.split("\t")
				DeseqScaleDic[li[0]]=float(li[1])
			if not SpikeScaleDic:
				for s in SampleScaleDic:
					for x in open(SampleScaleDic[s]).readlines()[1:]:
						x=x.rstrip()
						li=x.split("\t")
						SpikeScaleDic[li[0]]=float(li[3])
	else:
		os.mkdir("%s_deseq/all"%ConfigDic["ProjectName"])
		os.chdir(Path+"/%s_deseq/all"%ConfigDic["ProjectName"])
		pa.GetMergePeakCountMatrix(PeakDic,BamDic)
		CountMatrixFile="%s_counts_final.xls"%ConfigDic["ProjectName"]
		if SampleInputDic:
			for x in open("%s_deseq_scale.txt"%ConfigDic["ProjectName"]).readlines()[1:]:
				x=x.rstrip()
				li=x.split("\t")
				DeseqScaleDic[li[0]]=float(li[1])
			if not SpikeScaleDic:
				for s in SampleScaleDic:
					for x in open(SampleScaleDic[s]).readlines()[1:]:
						x=x.rstrip()
						li=x.split("\t")
						SpikeScaleDic[li[0]]=float(li[3])
		pa.GetDeseq2File(Path+"/"+ConfigDic["Target"],CountMatrixFile,SpikeScaleDic,DeseqScaleDic,Control=ConfigDic["ControlSample"])
		pa.GetAnnoPeak("merge_%s.bed"%ConfigDic["ProjectName"],ConfigDic["GeneBed"],glob.glob("*.xls"))
		os.chdir(Path)
	if "fwd" in ExistsDic["DeseqAnalysis"]:
		logging.info("\"fwd\" folder already exists. Deseq2 will not be run. If you want to do this analysis, please delete this folder.")
	else:
		os.mkdir("%s_deseq/fwd"%ConfigDic["ProjectName"])
		os.chdir(Path+"/%s_deseq/fwd"%ConfigDic["ProjectName"])
		pa.GetMergePeakCountMatrix(PeakFwdDic,BamFwdDic,MyPrefix=ConfigDic["ProjectName"]+"_fwd")
		pa.GetDeseq2File(Path+"/"+ConfigDic["Target"],"%s_fwd_counts_final.xls"%ConfigDic["ProjectName"],SpikeScaleDic,DeseqScaleDic,MyPrefix=ConfigDic["ProjectName"]+"_fwd")
		pa.GetAnnoPeak("merge_%s_fwd.bed"%ConfigDic["ProjectName"],ConfigDic["GeneBed"],glob.glob("*.xls"))
		os.chdir(Path)
	if "rev" in ExistsDic["DeseqAnalysis"]:
		logging.info("\"rev\" folder already exists. Deseq2 will not be run. If you want to do this analysis, please delete this folder.")
	else:
		os.mkdir("%s_deseq/rev"%ConfigDic["ProjectName"])
		os.chdir(Path+"/%s_deseq/rev"%ConfigDic["ProjectName"])
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
	ExistsDic=CheckContinue(ConfigDic)
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
	Flag=0
	for Fl in [BwFileList,FwdBwFileList,RevBwFileList]:
		for Fi in Fl:
			if not os.path.exists(Fi):
				Flag=1
				logging.info("%s not exists."%Fi)
	if Flag==1:
		print("Some files could not be found. The program exited. Please check the log")
		sys.exit()
	p=AnalysisPipLine(ConfigDic["ProjectName"],ConfigDic["Thread"])
	if ExistsDic["DownstreamAnalysis"]==[]:
		os.mkdir("%s_analysis"%ConfigDic["ProjectName"])
	if not "correlation" in ExistsDic["DownstreamAnalysis"]:
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
	else:
		logging.info("correlation folder already exists. Correlation analysis will not be done. If you want to do this analysis, please delete this folder.")
	if not "motif" in ExistsDic["DownstreamAnalysis"]:
		os.mkdir("%s_analysis/motif"%ConfigDic["ProjectName"])
	else:
		logging.info("motif folder already exists. Motif analysis will not be done. If you want to do this analysis, please delete this folder.")
	if not "sense_antisense" in ExistsDic["DownstreamAnalysis"]:
		os.mkdir("%s_analysis/sense_antisense"%ConfigDic["ProjectName"])
	else:
		logging.info("sense_antisense folder already exists. sense_antisense analysis will not be done. If you want to do this analysis, please delete this folder.")
	p.GetGenomeContentBedFile(ConfigDic["GeneBed"],ConfigDic["ContentExtend"],ConfigDic["ChromSize"])
	if not "skew" in ExistsDic["DownstreamAnalysis"]:
		os.mkdir("%s_analysis/skew"%ConfigDic["ProjectName"])
		os.chdir(Path+"/%s_analysis/skew/"%ConfigDic["ProjectName"])
		p.GetATGCSkewBw(ConfigDic["GenomeFastaFile"],ConfigDic["ChromSize"],ConfigDic["Win"],ConfigDic["Step"])
		GCSkewBw="%s_GCSkew.bw"%ConfigDic["ProjectName"]
		ATSkewBw="%s_ATSkew.bw"%ConfigDic["ProjectName"]
		p.GetSkewGz(GCSkewBw,ATSkewBw,ConfigDic["GeneBed"],ConfigDic["MetaplotExtend"],"gene")
		os.chdir(Path)
	else:
		logging.info("skew folder already exists. Skew analysis will not be done. If you want to do this analysis, please delete this folder.")
	PeakContentDic={}
	PeakLengthDic={}
	for s in SampleList:
		if not "peaks_content_distribution" in ExistsDic["DownstreamAnalysis"]:
			ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s/%s/%s_peaks.bed"%(Path,s,s),"tss_%s.bed"%ConfigDic["ContentExtend"],"tts_%s.bed"%ConfigDic["ContentExtend"],"genebody_%s.bed"%ConfigDic["ContentExtend"])
			PeakContentDic[s]=[ProNum,TerNum,BodyNum,InterNum]
			p.GetRandomBed("%s/%s/%s_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"],MyPrefix=s)
			ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s_random_all.bed"%(s),"tss_%s.bed"%ConfigDic["ContentExtend"],"tts_%s.bed"%ConfigDic["ContentExtend"],"genebody_%s.bed"%ConfigDic["ContentExtend"])
			PeakContentDic[s+"_random"]=[ProNum,TerNum,BodyNum,InterNum]
			p.PeakContentDistribution("%s/%s/%s_fwd_peaks.bed"%(Path,s,s),"tss_%s.bed"%ConfigDic["ContentExtend"],"tts_%s.bed"%ConfigDic["ContentExtend"],"genebody_%s.bed"%ConfigDic["ContentExtend"])
			PeakContentDic[s+"_fwd"]=[ProNum,TerNum,BodyNum,InterNum]
			p.GetRandomBed("%s/%s/%s_fwd_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"],MyPrefix=s+"_fwd")
			ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s_fwd_random_all.bed"%(s),"tss_%s.bed"%ConfigDic["ContentExtend"],"tts_%s.bed"%ConfigDic["ContentExtend"],"genebody_%s.bed"%ConfigDic["ContentExtend"])
			PeakContentDic[s+"_fwd_random"]=[ProNum,TerNum,BodyNum,InterNum]
			ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s/%s/%s_rev_peaks.bed"%(Path,s,s),"tss_%s.bed"%ConfigDic["ContentExtend"],"tts_%s.bed"%ConfigDic["ContentExtend"],"genebody_%s.bed"%ConfigDic["ContentExtend"])
			PeakContentDic[s+"_rev"]=[ProNum,TerNum,BodyNum,InterNum]
			p.GetRandomBed("%s/%s/%s_rev_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"],Strand="-",MyPrefix=s+"_rev")
			ProNum,TerNum,BodyNum,InterNum=p.PeakContentDistribution("%s_rev_random_all.bed"%(s),"tss_%s.bed"%ConfigDic["ContentExtend"],"tts_%s.bed"%ConfigDic["ContentExtend"],"genebody_%s.bed"%ConfigDic["ContentExtend"])
			PeakContentDic[s+"_rev_random"]=[ProNum,TerNum,BodyNum,InterNum]
		if not "peaks_length_distribution" in ExistsDic["DownstreamAnalysis"]:
			PeakLengthDic[s]=p.PeakLengthDistribution("%s/%s/%s_peaks.bed"%(Path,s,s))
			PeakLengthDic[s+"_fwd"]=p.PeakLengthDistribution("%s/%s/%s_fwd_peaks.bed"%(Path,s,s))
			PeakLengthDic[s+"_rev"]=p.PeakLengthDistribution("%s/%s/%s_rev_peaks.bed"%(Path,s,s))
		if not "sense_antisense" in ExistsDic["DownstreamAnalysis"]:
			FwdBw="%s/%s/%s_fwd_nucleus_norm.bw"%(Path,s,s)
			RevBw="%s/%s/%s_rev_nucleus_norm.bw"%(Path,s,s)
			os.chdir(Path+"/%s_analysis/sense_antisense/"%ConfigDic["ProjectName"])
			p.SenseAntisense(s,ConfigDic["GeneBed"],FwdBw,RevBw,ConfigDic["MetaplotExtend"])
		if not "skew" in ExistsDic["DownstreamAnalysis"]:
			os.chdir(Path+"/%s_analysis/skew/"%ConfigDic["ProjectName"])
			p.GetSkewGz(GCSkewBw,ATSkewBw,"%s/%s/%s_peaks.bed"%(Path,s,s),ConfigDic["MetaplotExtend"],s)
			p.GetSkewGz(GCSkewBw,ATSkewBw,"%s/%s/%s_fwd_peaks.bed"%(Path,s,s),ConfigDic["MetaplotExtend"],s+"_fwd")
			p.GetSkewGz(GCSkewBw,ATSkewBw,"%s/%s/%s_rev_peaks.bed"%(Path,s,s),ConfigDic["MetaplotExtend"],s+"_rev")
		if not "motif" in ExistsDic["DownstreamAnalysis"]:
			os.chdir(Path+"/%s_analysis/motif/"%ConfigDic["ProjectName"])
			p.FindMotif(ConfigDic["GenomeFastaFile"],"%s/%s/%s_fwd_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"],MyPrefix=ConfigDic["ProjectName"]+"_fwd",Strand="+",RepeatNum=ConfigDic["RepeatNum"])
			p.FindMotif(ConfigDic["GenomeFastaFile"],"%s/%s/%s_rev_peaks.bed"%(Path,s,s),ConfigDic["ChromSize"],MyPrefix=ConfigDic["ProjectName"]+"_rev",Strand="-",RepeatNum=ConfigDic["RepeatNum"])
		os.chdir(Path)
	if not "peaks_content_distribution" in ExistsDic["DownstreamAnalysis"]:
		df=pandas.DataFrame.from_dict(data=PeakContentDic,orient='index')
		os.mkdir("%s_analysis/peaks_content_distribution"%ConfigDic["ProjectName"])
		os.chdir(Path+"/%s_analysis/peaks_content_distribution/"%ConfigDic["ProjectName"])
		df.to_csv("%s_peaks_content_distribution.xls"%ConfigDic["ProjectName"],sep="\t",header=["ProNum","TerNum","BodyNum","InterNum"],index_label="Sample")
		os.chdir(Path)
	else:
		logging.info("peaks_content_distribution folder already exists. peaks_content_distribution analysis will not be done. If you want to do this analysis, please delete this folder.")
	if not "peaks_length_distribution" in ExistsDic["DownstreamAnalysis"]:
		os.mkdir("%s_analysis/peaks_length_distribution"%ConfigDic["ProjectName"])
		os.chdir(Path+"/%s_analysis/peaks_length_distribution/"%ConfigDic["ProjectName"])
		np.save('%s_peaks_length_distribution.npy'%ConfigDic['ProjectName'],PeakLengthDic)
		os.chdir(Path)
	else:
		logging.info("peaks_length_distribution folder already exists. peaks_length_distribution analysis will not be done. If you want to do this analysis, please delete this folder.")
	if not "cluster" in ExistsDic["DownstreamAnalysis"]:
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
	else:
		logging.info("cluster folder already exists. cluster analysis will not be done. If you want to do this analysis, please delete this folder.")
	#p.BedCorrelation()
	p.AfterMath()
def Main(ConfigFile,SubCommand):
	ConfigDic=json.load(open(ConfigFile))
	CheckJson(ConfigDic)
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
	SubCommand=sys.argv[2] #BaseAnalysis,DeseqAnalysis,DownstreamAnalysis,AllPip
	Main(ConfigFile,SubCommand)