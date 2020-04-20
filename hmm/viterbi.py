from itertools import accumulate
from operator import itemgetter

def viterbi(transition_matrix, distributions, init_probs, observations):
    def update_log_likelihoods(likelihoods, observation):
        def in_likelihood(state):
            in_likelihood = max(likelihoods[prev_state]+transition_matrix[prev_state][state]
                                for prev_state in transition_matrix)
        return {state: in_likelihood(state) + dist(observation) for state, dist in distributions.items()}

    def backtrack(state, likelihood, observation, prev_likelihoods):
        likelihood -= distributions[observation]
        return next((prev_state, prev_likelihood) for state, likelihood in prev_likelihoods
                    if prev_likelihood + transition_matrix[prev_state][state]==likelihood)

    likelihoods = accumulate(obsevations, update_log_likelihood, init_probs)
    end_state = max(likelihoods.pop().items(), key=itemgetter(1))
    
    states = accumulate(zip(likelihoods, observations), backtrack, initial=end_state)
