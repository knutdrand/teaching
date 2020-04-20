from dataclasses import dataclass, replace
from pprint import pprint
from itertools import accumulate
import math
LOG = True
@dataclass
class HMM:
    T: dict
    B: dict
    init_pi: dict

def log_sum_exp(logs):
    logs=list(logs)
    b = max(logs)
    return b + math.log(sum(math.exp(l-b) for l in logs))

def forward(hmm, observations):
    def update(likelihoods, observation):
        def update_state(state):
            return sum(likelihoods[prev_state]*hmm.T[prev_state][state] for prev_state in hmm.T)*hmm.B[state](observation)
        def update_log_state(state):
            return log_sum_exp([likelihoods[prev_state] + hmm.T[prev_state][state] for prev_state in hmm.T]) + hmm.B[state](observation)
        if LOG:
            update_state=update_log_state
        return {state: update_state(state) for state in hmm.T}

    return accumulate(observations, update, initial=hmm.init_pi)

def backward(hmm, observations):
    def update(likelihoods, observation):
        def update_state(state):
            return sum(
                hmm.T[state][next_state]*hmm.B[next_state](observation)*likelihoods[next_state]
                for next_state in hmm.T)
        def update_log_state(state):
            return log_sum_exp(
                [hmm.T[state][next_state]+hmm.B[next_state](observation)+likelihoods[next_state] for next_state in hmm.T])
        if LOG:
            update_state=update_log_state
        return {state: update_state(state) for state in hmm.T}
    init = 0 if LOG else 1
    return accumulate(reversed(observations), update, initial={state: init for state in hmm.T})

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

def update_log_hmm(forwards, backwards, hmm, observations):
    probs = [{state: f[state]+b[state] for state in hmm.T}
             for f, b in zip(forwards, backwards)]
    def transition(alpha, beta, obs):
        def prob(state, next_state):
            return alpha[state]+hmm.T[state][next_state]+beta[next_state]+hmm.B[next_state](obs)
        return {state: {next_state: prob(state, next_state) for next_state in hmm.T}
                for state in hmm.T}
    transitions = [transition(alpha, beta, obs) for 
                   alpha, beta, obs in zip(forwards, backwards[1:], observations)]
    transition_prob = {state: 
                       {next_state: 
                        log_sum_exp([t[state][next_state] for t in transitions])-log_sum_exp([p[state] for p in probs[:-1]]) for next_state in hmm.T} for state in hmm.T}
    return replace(hmm, T=transition_prob, init_pi={state: p-log_sum_exp(probs[0].values()) for state, p in probs[0].items()})

def baum_welch(observations, initial_hmm):
    _update_hmm=update_log_hmm if LOG else update_hmm
    def update(hmm):
        forwards = list(forward(hmm, observations))
        backwards = list(backward(hmm, observations))[::-1]
        return _update_hmm(forwards, backwards, hmm, observations)
    hmm = initial_hmm
    for _ in range(10):
        hmm = update(hmm)
        if LOG:
            pprint({state: {next_state: math.exp(p) for next_state, p in m.items()} for state, m in hmm.T.items()})
            pprint({state: math.exp(p) for state, p in hmm.init_pi.items()})
        else:
            pprint(hmm.T)
            pprint(hmm.init_pi)



def test():
    from numpy import random
    random.seed(42)
    T = {"A": {"A": 0.8, "B": 0.2},
         "B": {"A": 0.2, "B": 0.8}}
    if LOG:
        B = {"A": lambda x: math.log(0.9) if x==1 else math.log(0.1),
             "B": lambda x: math.log(0.9) if x==0 else math.log(0.1)}
    else:
        B = {"A": lambda x: 0.9 if x==1 else 0.1,
             "B": lambda x: 0.9 if x==0 else 0.1}
    init_pi = {"A": 0.5, "B": 0.5}

    if LOG:
        T = {state: {next_state: math.log(p) for next_state, p in m.items()} for state, m in T.items()}
        init_pi = {state: math.log(p) for state, p in init_pi.items()}
    print(B["A"](1))
    print(B["B"](1))
    hmm = HMM(T, B, init_pi)
    observations = ([random.binomial(1, 0.1) for _ in range(5)] + [random.binomial(1, 0.9) for _ in range(5)])*10
    observations = ([1]*10+[0]*10)*10
    pprint(observations)
    baum_welch(observations, hmm)


if __name__ == "__main__":
    test()
