import numpy as np
import random
from sklearn.utils import shuffle

def make_datas(data,select,plus = 0):
    datas = {}
    for idx in range(len(select)):
        datas[idx] = data[select[idx] + plus]
    return datas


def sample_once(datas, support_shot=5, query_shot=20, shuffle_ornot=True, plus = 0):
    np.random.seed(0)
    order = list(range(len(datas)))

    target_support = np.repeat(order, support_shot)
    target_query = np.repeat(order, query_shot)

    select = random.sample(range(len(datas[order[0] + plus])), support_shot + query_shot)
    sample = datas[order[0]+plus][select, :]
    inputs_support = sample[:support_shot]
    inputs_query = sample[support_shot:]

    for idx in order[1:]:
        select = random.sample(range(len(datas[idx+plus])), support_shot + query_shot)
        sample = datas[idx + plus][select, :]
        inputs_support = np.append(inputs_support, sample[:support_shot], axis=0)
        inputs_query = np.append(inputs_query, sample[support_shot:], axis=0)


    if shuffle_ornot:
        inputs_support, target_support = shuffle(inputs_support, target_support)
        inputs_query, target_query = shuffle(inputs_query, target_query)

    return inputs_support, inputs_query, target_support, target_query



