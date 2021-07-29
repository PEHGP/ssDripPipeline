#!/usr/bin/env python
from ssDripPipelines import DripMain,_version
import sys
if __name__ == "__main__":
	SubCommandList=["BaseAnalysis","DeseqAnalysis","DownstreamAnalysis","AllPip"]
	if len(sys.argv) != 3:
		print("\nssDRIPSeqAnalysis.py <DripConfig.json> <%s>\n"%("|".join(SubCommandList)))
		print("version:%s\n"%_version.__version__)
		sys.exit()
	ConfigFile=sys.argv[1]
	SubCommand=sys.argv[2]
	if not SubCommand in SubCommandList:
		print("\n%s is not a right method."%SubCommand)
		print("only %s can be selected.\n"%",".join(SubCommandList))
	else:
		DripMain.Main(ConfigFile,SubCommand)
