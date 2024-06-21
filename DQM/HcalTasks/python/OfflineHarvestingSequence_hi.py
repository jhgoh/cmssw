import FWCore.ParameterSet.Config as cms
from DQM.HcalTasks.HcalOfflineHarvesting_cfi import *

# apply some customization
# -	ptype = 1 Offlien processing
# - runkey value 2 - cosmics
hcalOfflineHarvesting.ptype = 1
hcalOfflineHarvesting.runkeyVal = 4
hcalOfflineHarvesting.runkeyName = "hi_run"
