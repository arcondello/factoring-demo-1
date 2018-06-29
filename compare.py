import shelve
import time
import os

from factoring.interfaces import factor


def findEmbeddings(out_directory, profile=None):
    """

    Args:
        out_directory: e.g. shelves/YOUR_TEST_NAME

        profile: as in dwave profile for the cloud-client

    """
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)

    for trial in range(1000):
        shelf = shelve.open(os.path.join(out_directory, str(time.time())))

        output49, embedding = factor(49, False, profile=profile)
        output21, _ = factor(21, False, embedding, profile=profile)
        output12, _ = factor(12, False, embedding, profile=profile)

        print('trial: ', trial)

        shelf['output49'] = output49
        shelf['output21'] = output21
        shelf['output12'] = output12
        shelf['embedding'] = embedding

        shelf.close()


def testEmbedding(in_file, out_directory, profile=None):
    shelf = shelve.open(in_file)
    embedding = shelf['embedding']
    shelf.close()

    if not os.path.exists(out_directory):
        os.mkdir(out_directory)

    for trial in range(1000):
        shelf = shelve.open(os.path.join(out_directory, str(time.time())))

        output49, _ = factor(49, False, embedding, profile=profile)
        output21, _ = factor(21, False, embedding, profile=profile)
        output12, _ = factor(12, False, embedding, profile=profile)

        shelf['output49'] = output49
        shelf['output21'] = output21
        shelf['output12'] = output12
        shelf['embedding'] = embedding

        shelf.close()


def resultsToPandas(in_directory):
    import os
    import pandas as pd
    import numpy as np
    from collections import defaultdict

    d = defaultdict(list)

    for filename in os.listdir(in_directory):
        shelf = shelve.open(os.path.join(in_directory, filename))

        # check that shelf is complete
        try:
            for P in [49,21,12]:
                scan = shelf['output%d' % P]
        except:
            continue

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

    return df


def getEmbedding(in_file):
    shelf = shelve.open(in_file)
    embedding = shelf['embedding']
    shelf.close()
    return embedding
