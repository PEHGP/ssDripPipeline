#!/usr/bin/env python
from piplines import DripMain
import sys
if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("ssDRIPSeqAnalysis.py <DripConfig.json> <BaseAnalysis|DeseqAnalysis|DownstreamAnalysis|AllPip>")
		sys.exit()
	ConfigFile=sys.argv[1]
	SubCommand=sys.argv[2]
	DripMain.Main(ConfigFile,SubCommand)