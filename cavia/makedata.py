import pandas as pd
import json
import os
from pandas import DataFrame
import matplotlib.pyplot as plt


def csv_to_perturb_data(path = os.getcwd(), name = "/cbc_gbc_dict.csv"):
    """ Turn the original cbc_gbc_dict csv to a simpler one for downstream analysis
    Especially the format provided in GSE90063_RAW

    Parameters
    ----------
    path : the path of your file, default is None
    name : the name of your file, default is "cbc_gbc_dict.csv"

    Returns
    -------
    a dataframe with two columns (sgRNA & gene)
    a list with all the sgRNAs
    """

    perturb = pd.read_csv(path + name, header=None)

    # delete the number
    perturb["rna"] = perturb[0].map(lambda x: x.split("_")[1])

    # a function to deal with one line in the original file
    def makepandas(i):
        y = perturb[1][i].split(",")
        t = [[perturb["rna"][i]] * len(y), y]
        return DataFrame(t).transpose()

    perturb_process = makepandas(0)

    for i in range(len(perturb)):
        perturb_process = perturb_process.append(makepandas(i), ignore_index=True)

    perturb_process.columns = ['sgRNA', 'cell']

    perturb = perturb_process.sort_values(by=['cell', 'sgRNA'], ignore_index=True)
    perturb['cell'] = perturb['cell'].str.strip()
    sgRNA_list = perturb['sgRNA'].drop_duplicates().tolist()

    return perturb, sgRNA_list


def json_to_perturb_data(path = os.getcwd(), name = "/cells_per_protospacer.json"):
    """ Turn the original cells_per_protospacer.json to a simpler one for downstream analysis
    Especially the format provided in our dataset

    Parameters
    ----------
    path : the path of your file, default is None
    name : the name of your file, default is "cbc_gbc_dict.csv"

    Returns
    -------
    a dataframe with two columns (sgRNA & gene)
    a list with all the sgRNAs
    """
    with open(path + name) as f:
        data = json.load(f)

        mock = list(data.keys())

    def create_one_from_dic(i):
        df = pd.DataFrame(data[mock[i]], columns=['cell'])
        df["sgRNA"] = "-".join(mock[i].split('-')[:-1])
        return df

    perturb_process = create_one_from_dic(0)

    for i in range(1, len(mock)):
        perturb_process = perturb_process.append(create_one_from_dic(i), ignore_index=True)

    perturb_process = perturb_process.sort_values(by=['sgRNA', 'cell'], ignore_index=True)
    perturb_process['cell'] = perturb_process['cell'].str.strip()
    sgRNA_list = perturb_process['sgRNA'].drop_duplicates().tolist()

    return perturb_process, sgRNA_list


def anndata_sgRNA(anndata, perturb, sgRNA_list = []):
    """ Annotate the sgRNA info to anndata

    Parameters
    ----------
    anndata : anndata
    perturb : the perturb info get from perturb_process
    sgRNA_list : sgRNA_list info get from perturb_process or perturb, it's okay to be none

    Return
    -------
    anndata with perturb info of each sgRNA

    """
    if sgRNA_list == []:
        sgRNA_list = perturb['sgRNA'].drop_duplicates().tolist()

    sgRNA_info = pd.DataFrame(0, index=anndata.obs.index.tolist(), columns=sgRNA_list)
    for i in range(len(perturb)):
        sgRNA_info[perturb["sgRNA"][i]][perturb["cell"][i]] += 1
    anndata.obs = pd.concat([anndata.obs, sgRNA_info.reindex(anndata.obs.index)], axis=1)

    return anndata


def hist_figure(anndata, sgRNA_list, title, type=""):
    """

    Parameters
    ----------
    anndata: anndata
    sgRNA_list: sgRNA_list of perturbations
    title: the title of the plot
    type:
     type = "perturbed" plot 'Detected perturbed times/cell'
     type = "sgRNA" plot 'Detected sgRNA/cell'
     type = "" plot both

     Return
    -------
    hist graph
    """
    def draw(t, xlabel, title):
        fig, ax = plt.subplots()
        ax.hist(t, range(t.value_counts().index[-1] + 1), density=False)

        for rect in ax.patches:
            height = rect.get_height()
            ax.annotate(f'{int(height)}', xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 5), textcoords='offset points', ha='center', va='bottom')

        plt.xlabel(xlabel)
        plt.ylabel('Cell number')
        plt.title(title)
        plt.show()

    if type == "" or type == "perturbed":
        t = anndata.obs.iloc[:, 0:len(sgRNA_list)].sum(axis=1)
        draw(t, 'Detected perturbed times/cell', title)

    if type == "" or type == "sgRNA":
        t = (anndata.obs.iloc[:, 0:len(sgRNA_list)] != 0).sum(axis=1)
        draw(t, 'Detected sgRNA/cell', title)

