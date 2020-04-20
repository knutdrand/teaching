from dataclasses import dataclass, replace
from itertools import accumulate

@dataclass
class HMM:
    T: dict
    B: dict
    init_pi: dict

def forward(hmm, observations):
    def update(likelihoods, observation):
        def update_state(state):
            return sum(likelihoods[prev_state]*hmm.T[prev_state][state] for prev_state in hmm.T)*hmm.B[state](observation)
                       
        return {state: update_state(state) for state in hmm.T}
    return accumulate(observations, update, initial=hmm.init_pi)

def backward(hmm, observations):
    def update(likelihoods, observation):
        def update_state(state):
            return sum(
                hmm.T[state][next_state]*hmm.B[next_state](observation)*likelihoods[next_state]
                for next_state in hmm.T)
            
        return {state: update_state(state) for state in hmm.T}
    return accumulate(reversed(observations), update, initial={state: 1 for state in hmm.T})

def forward_backward(hmm, observations):
    A = list(forward(hmm, observations))
    B = list(backward(hmm, observations))
    return ({state: a[state]*b[state] for state in hmm.T}
            for a, b in zip(A, reversed(B)))


def update_hmm(forwards, backwards, hmm, observations):
    probs = [{state: f[state]*b[state] for state in hmm.T}
             for f, b in zip(forwards, backwards)]
    def transition(alpha, beta, obs):
        def prob(state, next_state):
            return alpha[state]*hmm.T[state][next_state]*beta[next_state]*hmm.B[next_state](obs)
        return {state: {next_state: prob(state, next_state) for next_state in hmm.T}
                for state in hmm.T}
    transitions = [transition(alpha, beta, obs) for 
                   alpha, beta, obs in zip(forwards, backwards[1:], observations)]
    transition_prob = {state: {next_state: sum(t[state][next_state] for t in transitions)/sum(p[state] for p in probs[:-1]) for next_state in hmm.T} for state in hmm.T}
    return replace(hmm, T=transition_prob, init_pi={state: p/sum(probs[0].values()) for state, p in probs[0].items()})

def baum_welch(observations, initial_hmm):
    def update(hmm):
        forwards = list(forward(hmm, observations))
        backwards =list(backward(hmm, observations))[::-1]
        return update_hmm(forwards, backwards, hmm, observations)
    hmm = initial_hmm
    for _ in range(10):
        hmm = update(hmm)
        print(hmm)


def log_sum_exp(logs):
    b = max(logs)
    logs.remove(b)
    return b + math.log(sum(math.exp(l-b) for l in logs))

def test():
    from numpy import random
    random.seed(42)
    T = {"A": {"A": 0.8, "B": 0.2},
         "B": {"A": 0.2, "B": 0.8}}
    B = {"A": lambda x: 0.9 if x==1 else 0.1,
         "B": lambda x: 0.9 if x==0 else 0.1}
    init_pi = {"A": 0.5, "B": 0.5}
    hmm = HMM(T, B, init_pi)
    observations = ([random.binomial(1, 0.1) for _ in range(5)] + [random.binomial(1, 0.9) for _ in range(5)])*10
    observations = ([1]*10+[0]*10)*10
    print(observations)
    baum_welch(observations, hmm)


if __name__ == "__main__":
    test()
