#!/usr/bin/env python
import matplotlib,sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip,json
import pandas
import numpy as np
import scipy.signal
import collections
import seaborn as sns
def Filter(a):
	rl=np.array(a.dropna())
	n95=np.percentile(rl,95)
	n5=np.percentile(rl,5)
	rf=rl[(rl>n5)&(rl<n95)]
	mean=np.mean(rf)
	return mean
def Getgz(GzList,IfSmooth=False,IfFilter=False):
	for x in gzip.open(GzList[0]):
		x=x.rstrip()
		if x.startswith("@"):
			d=json.loads(x.replace("@",""))
			break
	dfl=[pandas.read_table(gzfile,sep="\t",index_col=[0,1,2,3,4,5],compression="gzip",header=None,skiprows=1) for gzfile in GzList]
	df=pandas.concat(dfl)
	if IfFilter:
		r=df.apply(Filter,axis=0)
	else:
		r=df.mean()
	print(r)
	#sys.exit()
	if IfSmooth:
		Mean=scipy.signal.savgol_filter(r[np.isfinite(r)],35,11)
	else:
		Mean=r  #right?
	Std=r.std(axis=0) #right ?
	print(Std)
	Left=1
	#print d
	TSS=int(d['upstream'][0])/int(d["bin size"][0])+1
	TTS=TSS+int(d['body'][0])/int(d["bin size"][0])-1
	Right=int(d['upstream'][0])/int(d["bin size"][0])+int(d['body'][0])/int(d["bin size"][0])+int(d['downstream'][0])/int(d["bin size"][0])
	return Mean,Std,Left,TSS,TTS,Right,d["sample_labels"],d["sample_boundaries"]
def PlotSenseAntisense(Target,Prefix,Extend):
	ColorList=np.array(plt.get_cmap("tab20").colors)
	GroupDict=collections.defaultdict(list)
	for x in open(Target):
		x=x.rstrip()
		l=x.split("\t")
		GroupDict[l[0]].append(l[1])
	ml=[]
	md={}
	LeftLabel="-1"+str(round(float(Extend)/1000,1))+"Kb"
	RightLabel=str(round(float(Extend)/1000,1))+"Kb"
	for s in ["sense","antisense"]:
		rl=[]
		for x in GroupDict:
			for xi in GroupDict[x]:
				Sefile.append("%s_%s.gz"%(xi,s))
			Mean,Std,Left,TSS,TTS,Right,_,_=Getgz(Sefile,True,True)
			md[x+"_"+s]=(Mean,Std,Left,TSS,TTS,Right)
			ml.append(Mean.max())
			Mean.name=x+"_"s
			rl.append(Mean)
		dr=pandas.concat(rl,axis=1)
		dr.index=range(0,dr.shape[0])
		dr["Pos"]="NaN"
		dr.loc[0,"Pos"]=LeftLabel
		dr.loc[dr.shape[0]-1,"Pos"]=RightLabel
		dr.loc[TSS-1,"Pos"]="TSS"
		dr.loc[TTS-1,"Pos"]="TTS"
		dr.to_csv(prefix+"_"+s+".xls",sep="\t",na_rep="NaN")
	ymax=int(max(ml))+1
	for s in ["sense","antisense"]:
		fig = plt.figure(dpi=300)
		ax = fig.add_subplot(1,1,1)
		ax.yaxis.grid(True,linestyle='dotted')
		HandleList=[]
		i=0
		for x in GroupDict:
			print(GroupDict[x])
			Mean,Std,Left,TSS,TTS,Right=md[x+"_"+s]
			if len(GroupDict.values())>20:
				r,=ax.plot(Mean,label=x)
			else:
				r,=ax.plot(Mean,label=x,color=ColorList[i])
			HandleList.append(r)
			#label.append(x+"_%s"%s)
			i+=1
		ax.set_xlim(Left,Right)
		ax.set_ylim(ymin=0,ymax=ymax)
		ax.legend(handles=HandleList)
		plt.xticks((Left,TSS,TTS,Right),(LeftLabel,"TSS","TTS",RightLabel))
		ax.set_ylabel("R-loop level")
		ax.set_title(s)
		fig.savefig("%s_mean_%s.png"%(Prefix,s),format='png')
		fig.savefig("%s_mean_%s.svg"%(Prefix,s),format='svg')
		fig.clf()
		plt.close(fig)

def PlotSkew(Target,Prefix):
	ColorList=np.array(plt.get_cmap("tab20").colors)
	GroupDict=collections.defaultdict(list)
	for x in open(Target):
		x=x.rstrip()
		l=x.split("\t")
		GroupDict[l[0]].append(l[1])
	#ml=[]
	md={}
	LeftLabel="-1"+str(round(float(Extend)/1000,1))+"Kb"
	RightLabel=str(round(float(Extend)/1000,1))+"Kb"
	for s in ["fwd","rev","all"]:
		rl=[]
		for x in GroupDict:
			for xi in GroupDict[x]:
				if s!="all":
					Sefile.append("%s_%s_GCATSkew.gz"%(xi,s))
				else:
					Sefile.append("%s_GCATSkew.gz"%(xi))
			Mean,Std,Left,TSS,TTS,Right,SampleLabels,SampleBoundaries=Getgz(Sefile,True,False)
			for i,sla in enumerate(SampleLabels):
				sa_start=SampleBoundaries[i]
				sa_end=SampleBoundaries[i+1]
				RealMean=Mean[sa_start:sa_end]
				md[x+"_"+s+"_"+sla]=(RealMean,Std,Left,TSS,TTS,Right)
				if s!="all":
					RealMean.name=x+"_"s+"_"+sla
				else:
					RealMean.name=x+"_"+sla
				rl.append(RealMean)
		dr=pandas.concat(rl,axis=1)
		dr.index=range(0,dr.shape[0])
		dr["Pos"]="NaN"
		dr.loc[0,"Pos"]=LeftLabel
		dr.loc[dr.shape[0]-1,"Pos"]=RightLabel
		dr.loc[TSS-1,"Pos"]="STAR"
		dr.loc[TTS-1,"Pos"]="END"
		if s!="all":
			dr.to_csv(prefix+"_"+s+"_GCATSkew.xls",sep="\t",na_rep="NaN")
		else:
			dr.to_csv(prefix+"_GCATSkew.xls",sep="\t",na_rep="NaN")
	#ymax=int(max(ml))+1
	for s in ["fwd","rev","all"]:
		fig = plt.figure(dpi=300)
		ax = fig.add_subplot(1,1,1)
		ax.yaxis.grid(True,linestyle='dotted')
		HandleList=[]
		i=0
		for x in GroupDict:
			print(GroupDict[x])
			for sla in ["GCSkew","ATSkew"]:
				Mean,Std,Left,TSS,TTS,Right=md[x+"_"+s+"_"+sla]
				if len(GroupDict.values()*2)>20:
					r,=ax.plot(Mean,label=x+"_"+sla)
				else:
					r,=ax.plot(Mean,label=x+"_"+sla,color=ColorList[i])
				HandleList.append(r)
			#label.append(x+"_%s"%s)
			i+=1
		ax.set_xlim(Left,Right)
		#ax.set_ylim(ymin=0,ymax=ymax)
		ax.legend(handles=HandleList)
		plt.xticks((Left,TSS,TTS,Right),(LeftLabel,"STAR","END",RightLabel))
		ax.set_ylabel("Skew level")
		ax.set_title(s)
		fig.savefig("%s_%s_GCATSkew.png"%(Prefix,s),format='png')
		fig.savefig("%s_%s_GCATSkew.svg"%(Prefix,s),format='svg')
		fig.clf()
		plt.close(fig)

	#gene_GCATSkew.gz
	
	Mean,Std,Left,TSS,TTS,Right,SampleLabels,SampleBoundaries=Getgz(["gene_GCATSkew.gz"],True,False)
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(1,1,1)
	ax.yaxis.grid(True,linestyle='dotted')
	HandleList=[]
	for i,sla in enumerate(SampleLabels):
		sa_start=SampleBoundaries[i]
		sa_end=SampleBoundaries[i+1]
		RealMean=Mean[sa_start:sa_end]
		r,=ax.plot(RealMean,label=sla)
		HandleList.append(r)
	ax.set_xlim(Left,Right)
	#ax.set_ylim(ymin=0,ymax=ymax)
	ax.legend(handles=HandleList)
	plt.xticks((Left,TSS,TTS,Right),(LeftLabel,"TSS","TTS",RightLabel))
	ax.set_ylabel("Skew level")
	fig.savefig("%s_gene_GCATSkew.png"%(Prefix),format='png')
	fig.savefig("%s_gene_GCATSkew.svg"%(Prefix),format='svg')
	fig.clf()
	plt.close(fig)
def PlotPeaksLengthDistribution(NpyFile):
	data=np.load(NpyFile,allow_pickle=True)
	d=data[()]
	rl=[]
	for s in d:
		for sc in d[s]:
def PlotPeaksContentDistribution():
	pass
def PlotScatplotDeseq():
	pass
def PlotSnapShot():
	pass