"""No solver calls: verify pre-ray shortcut accepts only an exact contradiction."""
import numpy as np
from types import SimpleNamespace
from scipy.optimize._highspy._core import HighsStatus
from planecount5 import LP
from scaled14 import Proxy
from lift11 import atomic
records=[]
for feasible in (False,True):
 lp=LP();j=lp.var('m1');lp.add([(j,1024)],'=',1024,'population');lp.add([(j,16)],'<',32 if feasible else 0,'upper')
 factors=np.array([1/1024,1/16]);dual=np.array([1/1024,-1/16])
 class Fake:
  called=False
  def getOptionValue(self,n):return HighsStatus.kOk,'ipm'
  def getSolution(self):return SimpleNamespace(row_dual=dual/factors)
  def getDualRay(self):self.called=True;return HighsStatus.kOk,False,[]
 fake=Fake();proxy=Proxy(fake,factors,lp);status,exists,_=proxy.getDualRay()
 assert not exists and fake.called==feasible
 records.append(dict(feasible=feasible,additional_solve_requested=fake.called))
atomic('out/preray14_check.json',dict(status='PRERAY14_EXACT_GATE_CHECKED',records=records));print('PRERAY14_EXACT_GATE_CHECKED',flush=True)
