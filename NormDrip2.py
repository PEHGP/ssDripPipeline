#!/usr/bin/env python2.7
#coding=utf-8
#做核基因组的归一化要去掉线粒体和叶绿体的基因组
#chromd={"nucleus":"chrom_nucleus.size","mitochondria":"chrom_mit.size","chloroplast":"chrom_chl.size"},必需有一个nucleus键
#ARTRLOOP
import json
import sys,os,re
class NormDrip(object):
	"""docstring for NormDrip"""
	def __init__(self, prefix,chromd,artrloop):
		super(NormDrip, self).__init__()
		self.prefix =prefix
		self.chrd={}
		self.artrloop=artrloop
		self.chrom_sized=chromd
		self.bacd={}
		for c in chromd:
			self.chrd[c]=[x.rstrip().split("\t")[0] for x in open(chromd[c])]
			for ci in self.chrd[c]:
				self.bacd[ci]=c
	def GetBdg(self,bam,fragment,bin_size,bdg_prefix):#fragment=True or False,bin_size????
		if fragment:
			os.system("bamCoverage --extendReads -v -p 10 -b %s -o %s.bdg --binSize %s --outFileFormat bedgraph >bc_infot.txt"%(bam,bdg_prefix,bin_size))
		else:
			os.system("bamCoverage -v -p 10 -b %s -o %s.bdg --binSize %s --outFileFormat bedgraph >bc_info.txt"%(bam,bdg_prefi,bin_size))
	def GetSome(self,bam,fragment):#fragment=True or False
		self.GetBdg(bam,fragment,self.prefix)
		self.SplitBdg(self.prefix,"True")
	def SplitBdg(self,bdg_prefix,ifstat): #ifstat=True or False
		if ifstat:
			depth_tempd={}
			self.covd={}
			self.art=0
		Frd={}
		for c in self.chrd:
			Frd[c]=open(bdg_prefix+"_"+c+".bdg","w")
			if ifstat:
				depth_tempd[c]=[0,0]
				self.covd[c]=0
		for c in open(bdg_prefix+".bdg"):
			c=c.rstrip()
			ci=c.split("\t")
			if ci[0]!=self.artrloop:
				Frd[self.bacd[ci[0]]].write(c+"\n")
			if ifstat:
				if ci[0]==self.artrloop:
					if int(ci[-1])>self.art:
						self.art=int(ci[-1])
					continue
				if int(ci[-1])>0:
					self.covd[self.bacd[ci[0]]]+=int(ci[2])-int(ci[1])
				depth_tempd[self.bacd[ci[0]]][0]+=int(ci[-1])*(int(ci[2])-int(ci[1]))
				depth_tempd[self.bacd[ci[0]]][1]+=int(ci[2])-int(ci[1])
		for c in Frd:
			Frd[c].close()
		if ifstat:
			self.depth={}
			for c in depth_tempd:
				self.depth[c]=depth_tempd[c][0]/float(depth_tempd[c][1])
	def GetScaleFactor_art(self,control_depth,control_art): #need change
		self.chl_scale_factor=1/self.chl_mean_depth
		self.mit_scale_factor=1/self.mit_mean_depth
		self.nucleus_scale_factor=1/(control_depth*(self.art/float(control_art)))
	def GetScaleFactor(self,):
		self.scale_factord={}
		for c in self.chrd:
			if self.depth[c]!=0:
				self.scale_factord[c]=1/self.depth[c]
			else:
				self.scale_factord[c]=1
	def MultiplySacleFactor(self,bdg_prefix):
		for c in self.chrd:
			os.system("awk -F\'\t\' \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4*%s}\' %s >%s"%(self.scale_factord[c],bdg_prefix+"_"+c+".bdg",bdg_prefix+"_"+c+"_norm.bdg"))
	def Getbw(self,bdg_prefix):
		for c in self.chrom_sized:
			inp=bdg_prefix+"_"+c+"_norm.bdg"
			out=bdg_prefix+"_"+c+"_norm.bw"
			os.system("bedGraphToBigWig %s %s %s"%(inp,self.chrom_sized[c],out))
	def CallPeak(self,bam,results_prefix):#check
		self.effective_genome=sum([self.covd[c] for c in self.covd])
		os.system("macs2 callpeak -t %s -f BAMPE -g %s -n %s"%(bam,self.effective_genome,results_prefix))#may be need change
		Frd={}
		for c in self.chrd:
			Frd[c]=open("%s_%s_peaks.bed"%(results_prefix,c),"w")
		for x in open(results_prefix+"_peaks.xls"):
			x=x.rstrip()
			if x.startswith("#"):
				continue
			if not x:
				continue
			l=x.split("\t")
			if l[1]=="start":#need change?
				continue
			Frd[self.bacd[l[0]]].write(l[0]+"\t"+l[1]+"\t"+l[2]+"\t"+l[-1]+"\n")
		for c in Frd:
			Frd[c].close()
	def SplitBam(self,bam):
		os.system("SplitBamAsStrand.sh %s"%bam)
		os.system("mv fwd.bam %s_fwd.bam"%self.prefix)
		os.system("mv rev.bam %s_rev.bam"%self.prefix)
		os.system("rm -rf rev*.bam")
		os.system("rm -rf fwd*.bam")
		os.system("samtools index %s_fwd.bam"%self.prefix)
		os.system("samtools index %s_rev.bam"%self.prefix)
	def NormAndCallPeakPipline_noart(self,bam,fragment,bin_size,if_split,logfile):
		Fr=open(logfile,"w")
		Fr2=open(logfile+"_progress","w")
		self.GetBdg(bam,fragment,bin_size,self.prefix)
		self.SplitBdg(self.prefix,True)
		self.GetScaleFactor()
		self.MultiplySacleFactor(self.prefix)
		self.Getbw(self.prefix)
		self.CallPeak(bam,self.prefix)
		Fr2.write(json.dumps(self.depth)+"\n")
		Fr2.write(json.dumps(self.covd)+"\n")
		Fr2.write(json.dumps(self.scale_factord)+"\n")
		Fr2.write(str(self.effective_genome)+"\n")
		based={}
		for c in self.chrom_sized:
			based[c]=0
			for x in open(self.chrom_sized[c]):
				x=x.rstrip()
				l=x.split("\t")
				based[c]+=int(l[1])
		realcovd={}
		for c in self.chrd:
			realcovd[c]=self.covd[c]/float(based[c])
		if if_split:
			Fr2.write("split follow\n")
			#self.SplitBam(bam)
			for s in ["fwd","rev"]:
				self.GetBdg(self.prefix+"_"+s+".bam",fragment,bin_size,self.prefix+"_"+s)
				self.SplitBdg(self.prefix+"_"+s,False)
				self.MultiplySacleFactor(self.prefix+"_"+s)#scale factor用的是全基因组的
				self.Getbw(self.prefix+"_"+s)
				self.CallPeak(self.prefix+"_"+s+".bam",self.prefix+"_"+s) #effective genome size用的是全基因组的
			Fr2.write(json.dumps(self.depth)+"\n")
			Fr2.write(json.dumps(self.covd)+"\n")
			Fr2.write(json.dumps(self.scale_factord)+"\n")
			Fr2.write(str(self.effective_genome)+"\n")
		fmdepth=""
		fmcov=""
		fmscale=""
		headerdepth=""
		headercov=""
		headerscale=""
		for c in self.chrd:
			headerdepth+=c+"_mean_depth\t"
			headercov+=c+"_cov\t"
			headerscale+=c+"_scale_factor\t"
			fmdepth+=str(self.depth[c])+"\t"
			fmcov+=str(realcovd[c])+"\t"
			fmscale+=str(self.scale_factord[c])+"\t"
		Fr.write("art\t"+headerdepth+headercov+headerscale+"effective_genome"+"\n")
		Fr.write(str(self.art)+"\t"+fmdepth+fmcov+fmscale+str(self.effective_genome)+"\n")
		Fr.close()
		Fr2.close()
if __name__ == '__main__':
	prefix=sys.argv[1]
	bam=sys.argv[2]
	bin_size=sys.argv[3]
	#control_depth=sys.argv[3]
	#control_art=sys.argv[4]
	#chromd={"nucleus":"chrom_nucleus.size","mitochondria":"chrom_mit.size","chloroplast":"chrom_chl.size"} #need change
	chromd={"nucleus":"chrom_nucleus.size","mit":"chrom_mit.size","chl":"chrom_chl.size","chrB":"chrom_chrB.size"}
	p=NormDrip(prefix,chromd,"ARTRLOOP")
	p.NormAndCallPeakPipline_noart(bam,True,bin_size,True,prefix+"_norm.log")
	#print p.nucleus_mean_depth
	#print p.chl_mean_depth
	#print p.mit_mean_depth
	#print p.art
	#print p.nucleus_cov
	#print p.chl_cov
	#print p.mit_cov
