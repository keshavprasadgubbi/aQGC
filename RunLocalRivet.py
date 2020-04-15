#This is the JO file for running the MG events usin the cmd:  athena RunLocalRivet.py
theApp.EvtMax = -1
#the number one sets here is the number of events that we give for the rivet routine to run on. -1 would just mean run on max number of events.
import os
import AthenaPoolCnvSvc.ReadAthenaPool

svcMgr.EventSelector.InputCollections = []
#OutputYoda = 'aQGC_rivetAnalysis-addnlObservables.yoda'
#OutputRoot  ='aQGC_rivetAnalysis-SM-newv5-rivetroutine.root'
#OutputRoot  ='aQGC_rivetAnalysis-M5-bestlimits-newv5rivetroutine.root'
#OutputRoot  ='aQGC_rivetAnalysis-M4-bestlimits-newv5rivetroutine.root'
#OutputRoot  ='aQGC_rivetAnalysis-M3-bestlimits-newv5rivetroutine.root' 
#OutputRoot  = 'aQGC_rivetAnalysis-M2-bestlimits-newv5rivetroutine.root'
OutputRoot = 'aQGC_rivetAnalysis-T5-bestlimits-newv5rivetroutine.root'

#SM-5FS
#svcMgr.EventSelector.InputCollections = ["/afs/cern.ch/user/n/ngubbi/work/AnalysisFiles-aQGC-WWZ/NewObservables-samples/DownloadedSamples/user.ngubbi.aQGC-T0-SM-5FS-newsample.10k-3lepevents.v2_EXT0/user.ngubbi.16007968.EXT0._000001.EVNT-aQGC-01.pool-T0-SM-5FS-newsample.10k-3lepevents.v2.root"]
#M5-Best-Limits
#svcMgr.EventSelector.InputCollections = ["/afs/cern.ch/user/n/ngubbi/work/AnalysisFiles-aQGC-WWZ/NewObservables-samples/DownloadedSamples/user.ngubbi.aQGC-M5-bestlimits-65-10pow12.10k-3lepevents.v1_EXT0/user.ngubbi.16196066.EXT0._000001.EVNT-aQGC-01.pool-M5-bestlimits-65-10pow12.10k-3lepevents.v1.root"]
#M4--Best Limits
#svcMgr.EventSelector.InputCollections = ["/afs/cern.ch/user/n/ngubbi/work/AnalysisFiles-aQGC-WWZ/NewObservables-samples/DownloadedSamples/user.ngubbi.aQGC-M4-bestlimits-40-10pow12.10k-3lepevents.v1_EXT0/user.ngubbi.16195971.EXT0._000001.EVNT-aQGC-01.pool-M4-bestlimits-40-10pow12.10k-3lepevents.v1.root"]
#M2-Best Limits                                                                                                                                                                                
#svcMgr.EventSelector.InputCollections = ["/afs/cern.ch/user/n/ngubbi/work/AnalysisFiles-aQGC-WWZ/NewObservables-samples/DownloadedSamples/user.ngubbi.aQGC-M2-bestlimits-26-10pow12.10k-3lepevents.v1_EXT0/user.ngubbi.16192461.EXT0._000001.EVNT-aQGC-01.pool-M2-bestlimits-26-10pow12.10k-3lepevents.v1.root"]
#M3 Best Limits
#svcMgr.EventSelector.InputCollections = ["/afs/cern.ch/user/n/ngubbi/work/AnalysisFiles-aQGC-WWZ/NewObservables-samples/DownloadedSamples/user.ngubbi.aQGC-M3-bestlimits-43.5-10pow12.10k-3lepevents.v1_EXT0/user.ngubbi.16193004.EXT0._000001.EVNT-aQGC-01.pool-M3-bestlimits-43.5-10pow12.10k-3lepevents.v1.root"]
#T5 Best Limits
svcMgr.EventSelector.InputCollections = ["/afs/cern.ch/user/n/ngubbi/work/AnalysisFiles-aQGC-WWZ/NewObservables-samples/DownloadedSamples/user.ngubbi.aQGC-T5-bestlimits-3.8-10pow12.10k-3lepevents.v1_EXT0/user.ngubbi.16179036.EXT0._000001.EVNT-aQGC-01.pool-T5-bestlimits-3.8-10pow12.10k-3lepevents.v1.root"]

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from Rivet_i.Rivet_iConf import Rivet_i

rivet = Rivet_i()
rivet.AnalysisPath = os.environ['PWD']
rivet.Analyses += [ 'WWZ' ]#The name of analysis here must be the same as the name of the class in thr rivet routine, else it will not know which analysis to run and also the respective rivet routine needs to be available in the same folder as this JO so that you can give PWD in the previous line for the analysis path. Else specify the path where rivet routine is located. 
rivet.RunName = ""
#rivet.HistoFile = "aQGC_rivetAnalysis-addnlObservables.yoda"
#rivet.HistoFile = OutputYoda
rivet.HistoFile = "WWZ"
#OutputYoda =  'WWZ-v2'
#OutputRoot  = 'WWZ-v2'

from AthenaCommon.AppMgr import ServiceMgr as svcMgr
from GaudiSvc.GaudiSvcConf import THistSvc
svcMgr += THistSvc()
svcMgr.THistSvc.Output = ["Rivet DATAFILE='"+OutputRoot+"' OPT='RECREATE'"]
 
job += rivet
