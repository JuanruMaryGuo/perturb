import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import random

def QC_control_plot(anndata):
    # Thresholding decision: counts

    p1 = sb.distplot(anndata.obs['n_counts'], kde=False)
    plt.show()

    p2 = sb.distplot(anndata.obs['n_counts'][anndata.obs['n_counts'] < 3000], kde=False, bins=60)
    plt.show()

    p3 = sb.distplot(anndata.obs['n_counts'][anndata.obs['n_counts'] > 20000], kde=False, bins=60)
    plt.show()

    # Thresholding decision: genes
    p4 = sb.distplot(anndata.obs['n_genes'], kde=False, bins=60)
    plt.show()

    p5 = sb.distplot(anndata.obs['n_genes'][anndata.obs['n_genes'] < 1500], kde=False, bins=60)
    plt.show()


def split_catagory(dic, ratio=0.7):
    # select_train = random.sample(range(len(dic)), int(np.floor(ratio * len(dic))))
    # select_test = list(set(range(len(dic))) - set(select_train))
    # X_train = {key: dic[key] for key in select_train}
    # X_test = {key: dic[key] for key in select_test}

    X_train = dict(list(dic.items())[:int(np.floor(ratio * len(dic)))])
    X_test = dict(list(dic.items())[int(np.floor(ratio * len(dic))):])

    return X_train, X_test


def split(dic, ratio=0.7):
    X_train = {}
    X_test = {}

    for idx in range(len(dic)):
        select = random.sample(range(len(dic[idx])), int(np.floor(ratio * len(dic[idx]))))
        X_train[idx] = dic[idx][select, :]
        X_test[idx] = dic[idx][list(set(range(len(dic[idx]))) - set(select)), :]

    return X_train, X_test


def delete_zero(X):
    return X[~np.all(X == 0, axis=1)]


def make_one_perturbed_data(anndata, idx=None, sgRNA_list=[], batch=False):
    if batch:

        if idx != None and idx < len(sgRNA_list):
            sgRNA = sgRNA_list[idx]
            X_is_pertrubed = (anndata.obs[sgRNA] != 0).values.reshape(-1, 1)
            X_train = (anndata.obs['batch'] == "0").to_numpy(dtype=int).reshape(-1, 1)
            X_test = (anndata.obs['batch'] == "1").to_numpy(dtype=int).reshape(-1, 1)

            X_pertrubed_train = delete_zero(X_is_pertrubed * anndata.obsm['X_pca'] * X_train)
            X_pertrubed_test = delete_zero(X_is_pertrubed * anndata.obsm['X_pca'] * X_test)

            return X_pertrubed_train, X_pertrubed_test

        if idx == None and sgRNA_list != []:
            X_is_pertrubed = (anndata.obs[sgRNA_list].sum(axis=1) == 0).values.reshape(-1, 1)
            X_train = (anndata.obs['batch'] == "0").to_numpy(dtype=int).reshape(-1, 1)
            X_test = (anndata.obs['batch'] == "1").to_numpy(dtype=int).reshape(-1, 1)

            X_pertrubed_train = delete_zero(X_is_pertrubed * anndata.obsm['X_pca'] * X_train)
            X_pertrubed_test = delete_zero(X_is_pertrubed * anndata.obsm['X_pca'] * X_test)

            return X_pertrubed_train, X_pertrubed_test

        else:
            print('No sgRNA found1')

    else:

        if idx != None and idx < len(sgRNA_list):
            sgRNA = sgRNA_list[idx]
            X_is_pertrubed = (anndata.obs[sgRNA] != 0).values.reshape(-1, 1)
            X = delete_zero(X_is_pertrubed * anndata.obsm['X_pca'])
            return X

        elif idx == None and sgRNA_list != []:
            sgRNA = sgRNA_list[idx]
            X_is_pertrubed = (anndata.obs[sgRNA] == 0).values.reshape(-1, 1)
            X = delete_zero(X_is_pertrubed * anndata.obsm['X_pca'])

            return X

        else:
            print('No sgRNA found2')


def make_total_data(anndata, sgRNA_list=[], non_perturbed=False, batch=False, ratio=0.7):
    if batch:
        if sgRNA_list != []:
            X_train = {}
            X_test = {}

            for idx in range(len(sgRNA_list)):
                X_train[idx], X_test[idx] = make_one_perturbed_data(anndata, idx, sgRNA_list, batch=True)

            if non_perturbed:
                X_train[len(sgRNA_list)], X_test[len(sgRNA_list)] = make_one_perturbed_data(anndata, None, sgRNA_list,
                                                                                            batch=True)

            return X_train, X_test
        else:
            print('No sgRNA found3')

    else:
        if sgRNA_list != []:
            X = {}

            for idx in range(len(sgRNA_list)):
                X[idx] = make_one_perturbed_data(anndata, idx, sgRNA_list)

            if non_perturbed:
                X[len(sgRNA_list)] = make_one_perturbed_data(anndata, None, sgRNA_list)

            X_train, X_test = split_catagory(X, ratio=ratio)

            return X_train, X_test
        else:
            print('No sgRNA found4')
