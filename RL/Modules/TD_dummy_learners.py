import numpy as np

def run_TD_episode(dumWorld,learner,rwd=None,pltEtrace=False):
	dumWorld.set_state(0)
	deltas = []
	if pltEtrace:
		fig = plt.figure()
		plt.plot(learner.etrace,alpha=.05/dumWorld.nstates,color='k')
	while dumWorld.state != dumWorld.term_state:
		s = dumWorld.state
		if learner.use_etrace:
			learner.update_etrace(dumWorld.state)
			if pltEtrace:
				plt.plot(learner.etrace,alpha=(s+.05)/dumWorld.nstates,color='k')
		dumWorld.take_action()
		sprime = dumWorld.state
		r = dumWorld.get_rwd(rwd)
		rpe = learner.update_valFunc(r,s,sprime)
		deltas.append(rpe)
	if pltEtrace:
		plt.ylabel("update eligibility")
		plt.xlabel("state")
		plt.tight_layout()
		return r,deltas,fig    
	return r,deltas

def run_lookAhead_episode(dumWorld,learner,rwd=None):
	dumWorld.set_state(0)
	states2go = dumWorld.nstates-2
	while dumWorld.state != dumWorld.term_state:
		s = dumWorld.state
		dumWorld.take_action()
		learner.update_valFunc(s,states2go)
		states2go -= 1
	r = dumWorld.get_rwd(rwd)
	learner.update_termVal(r)


class DummyWorld():

	def __init__(self,nstates,pRwd = 80):
		self.nstates = nstates
		self.pRwd = pRwd

	def initialize_world(self):
		self.set_states()
		self.create_policy()
		self.set_termState()

	def set_states(self):
		self.states = np.arange(self.nstates)

	def create_policy(self):
		self.policy = {self.states[s]:self.states[s+1] \
		for s in range(self.nstates-1)}

	def set_termState(self):
		self.term_state = self.nstates-1

	def create_etrace(self):
		self.etrace = np.zeros(self.nstates)

	def set_state(self,s):
		self.state = s

	def take_action(self):
		self.state = self.policy[self.state]

	def get_rwd(self,rwd=None):
		if self.state == self.term_state:
			if rwd != None:
				return rwd
			draw = np.random.randint(100)
			rwd = 1 if draw <= self.pRwd else 0
			return rwd
		return 0

class TdDummyLearner():

	def __init__(self,alpha=0.3,gamma=0.9,e_lambda=0.9,use_etrace=False):
		self.use_etrace = use_etrace
		self.alpha = alpha
		self.gamma = gamma
		self.e_lambda = e_lambda

	def set_valFunc(self,V):
		self.V = V

	def set_etrace(self,etrace):
		self.etrace = etrace

	def update_etrace(self,state):
		'''implement replacing eligibility trace'''
		self.etrace = self.gamma*self.e_lambda*self.etrace
		self.etrace[state] = 1

	def update_valFunc(self,r,s,sprime):
		delta = r + self.gamma * self.V[sprime] - self.V[s]
		if self.use_etrace:
			self.V = self.V + self.alpha * self.etrace * delta
		else:
			self.V[s] = self.V[s] + self.alpha * delta
		return delta

class TdDummyLearner_VariableLambda():

	def __init__(self,alpha=0.3,gamma=0.9,e_lambda=0.9,use_etrace=False,term_state=10):
		self.use_etrace = use_etrace
		self.alpha = alpha
		self.gamma = gamma
		self.e_lambda = e_lambda
		self.term_state = term_state

	def set_valFunc(self,V):
		self.V = V

	def set_etrace(self,etrace):
		self.etrace = etrace

	def update_etrace(self,state):
		'''implement replacing eligibility trace'''
		if state == self.term_state-2:
			stateLambda = self.e_lambda*0.75
		elif state == self.term_state-1:
			stateLambda = self.e_lambda*0.5
		elif state == self.term_state:
			stateLambda = 0.0
		else:
			stateLambda=self.e_lambda
		self.etrace = self.gamma*stateLambda*self.etrace
		self.etrace[state] = 1

	def update_valFunc(self,r,s,sprime):
		delta = r + self.gamma * self.V[sprime] - self.V[s]
		if self.use_etrace:
			self.V = self.V + self.alpha * self.etrace * delta
		else:
			self.V[s] = self.V[s] + self.alpha * delta
		return delta


class LookAheadDummyLearner():

	def __init__(self,alpha=0.3,gamma=0.9):
		self.alpha = alpha
		self.gamma = gamma

	def set_valFunc(self,V):
		self.V = V

	def set_termVal(self,v):
		self.termVal = v

	def update_termVal(self,r):
		self.termVal = self.termVal + self.alpha * (r - self.termVal)

	def update_valFunc(self,s,distToTerm):
		self.V[s] = self.gamma ** distToTerm * self.termVal



