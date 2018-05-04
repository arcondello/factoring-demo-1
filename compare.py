# import pickle
import shelve
import time

from factoring.interfaces import factor
# from factoring.fixed_embedding import factor


if __name__ == '__main__':
    # shelf = shelve.open("sweep-P/output1524593053.87")
    # embedding = shelf['embedding']
    # shelf.close()
    for trial in range(1000):
        shelf = shelve.open("sweep-SIMULATED_2000Q_CLOUD_TEST/output%s" % time.time())

        output49, embedding = factor(49)#, embedding) # comment out input embedding to generate new one

        output21, _ = factor(21, embedding)

        output12, _ = factor(12, embedding)

        shelf['output49'] = output49
        shelf['output21'] = output21
        shelf['output12'] = output12

        # bqm = get_factor_bqm(0)
        # response, embedding = submit_factor_bqm(bqm)
        # shelf['embedding'] = embedding
        # output = postprocess_factor_response(response, 0)
        # shelf['output%d' % 0] = output
        # for P in {a * b for a in range(1, 2 ** 3) for b in range(1, 2 ** 3)}:
        #     bqm = get_factor_bqm(P)
        #     response, _ = submit_factor_bqm(bqm, embedding)
        #     output = postprocess_factor_response(response, P)
        #     shelf['output%d' % P] = output

        shelf.close()

        # output = factor(P)
        # with open("pickle/fixed%s" % time.time(), "wb") as f:
        #     pickle.dump(output, f)


import os
import pandas as pd
import numpy as np
from collections import defaultdict

d = defaultdict(list)

for filename in os.listdir('.'):
    # if filename.startswith('output'):
    shelf = shelve.open(filename)
    try:
        for P in [49,21,12]:#{a * b for a in range(2 ** 3) for b in range(2 ** 3)}:
            scan = shelf['output%d' % P]
    except:
        continue
    # for P in{a * b for a in range(2 ** 3) for b in range(2 ** 3)}:
    #     percent = sum([x['percentageOfOccurrences'] for x in shelf['output%d' % P]['results'] if x['valid']])
    #     d[P].append(percent)
    # #     minpercent = 100
    # #     if percent < minpercent:
    # #         minpercent = percent
    # # d['minpercent'] = minpercent
    d['percent49'].append(sum([x['percentageOfOccurrences'] for x in shelf['output49']['results'] if x['valid']]))
    d['percent21'].append(sum([x['percentageOfOccurrences'] for x in shelf['output21']['results'] if x['valid']]))
    d['percent12'].append(sum([x['percentageOfOccurrences'] for x in shelf['output12']['results'] if x['valid']]))
    # for k, v in shelf['embedding'].items():
    #     d[k].append(len(v))
    d['filename'].append(filename)
    shelf.close()

df = pd.DataFrame(data=d)
# print(df.replace(0, np.nan).T.agg(['sum', 'mean', 'std', 'var', 'min', 'max', 'count']).T)
# pd.concat([df.agg(['sum', 'mean', 'std', 'var', 'min', 'max', 'count']).T.reset_index(),
#            pd.DataFrame([list("{:06b}".format(P)) for P in df.T.index.tolist()])],
#           axis=1).sort_values('mean')
df['minpercent'] = df.min(axis=1)
df = df.sort_values('minpercent')
# # df = df[(df.T != 0).all()]
pd.set_option('expand_frame_repr', False)
pd.set_option('max_columns', 50)
