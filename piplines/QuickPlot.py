#!/usr/bin/env python
import matplotlib,sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip,json
import pandas
import numpy as np
import scipy.signal
import collections
def Filter(a):
	rl=np.array(a.dropna())
	n95=np.percentile(rl,95)
	n5=np.percentile(rl,5)
	rf=rl[(rl>n5)&(rl<n95)]
	mean=np.mean(rf)
	return mean
def Getgz(GzList,IfSmooth=False):
	for x in gzip.open(GzList[0]):
		x=x.rstrip()
		if x.startswith("@"):
			d=json.loads(x.replace("@",""))
			break
	dfl=[pandas.read_table(gzfile,sep="\t",index_col=[0,1,2,3,4,5],compression="gzip",header=None,skiprows=1) for gzfile in GzList]
	df=pandas.concat(dfl)
	r=df.apply(Filter,axis=0)
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
	return Mean,Std,Left,TSS,TTS,Right
def PlotSenseAntisense(Target,Prefix):
	ColorList=np.array(plt.get_cmap("tab20").colors)
	GroupDict=collections.defaultdict(list)
	for x in open(Target):
		x=x.rstrip()
		l=x.split("\t")
		GroupDict[l[0]].append(l[1])
	for s in ["sense","antisense"]:
		fig = plt.figure(dpi=300)
		ax = fig.add_subplot(1,1,1)
		ax.yaxis.grid(True,linestyle='dotted')
		HandleList=[]
		i=0
		for x in GroupDict:
			Sefile=[]
			print(GroupDict[x])
			for xi in GroupDict[x]:
				Sefile.append("%s_%s_%s.gz"%(xi,prefix,s))
			Mean,Std,Left,TSS,TTS,Right=Getgz(Sefile,True)
			r,=ax.plot(Mean,label=x+"_%s"%s,color=ColorList[i])
			HandleList.append(r)
			#label.append(x+"_%s"%s)
			i+=1
		ax.set_xlim(Left,Right)
		ax.set_ylim(ymin=0,ymax=3.0)
		ax.legend(handles=HandleList)
		plt.xticks((Left,TSS,TTS,Right),("-1.0Kb","TSS","TTS","1.0Kb"))
		ax.set_ylabel("R-loop level")
		fig.savefig("%s_mean_%s.png"%(Prefix,s),format='png')
		fig.savefig("%s_mean_%s.svg"%(Prefix,s),format='svg')
		fig.clf()
		plt.close(fig)
		mean.to_csv(prefix+"_"+s+".xls",sep="\t",na_rep="NaN")