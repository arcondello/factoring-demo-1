import pickle
import time

from factoring.interfaces import get_factor_bqm, submit_factor_bqm, postprocess_factor_response
from factoring.fixed_embedding import factor


if __name__ == '__main__':
    P = 49

    for trial in range (1000):
        bqm = get_factor_bqm(P)
        response = submit_factor_bqm(bqm)
        output = postprocess_factor_response(response, P)
        with open("pickle/orig%s" % time.time(), "wb") as f:
            pickle.dump(output, f)

        output = factor(P)
        with open("pickle/fixed%s" % time.time(), "wb") as f:
            pickle.dump(output, f)

# for filename in os.listdir('.'):
#     if filename.startswith('fixed'):
#     with open(filename) as f:
#                 stuff = pickle.load(f)
#     print([x['percentageOfOccurrences'] for x in stuff['results'] if x['valid']])
