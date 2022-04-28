import inline as inline
import matplotlib
from pandas import read_csv

from cluster_loops import fast_rmsd
import numpy as np
from numpy import zeros
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns; sns.set()  # for plot styling
import numpy as np

def calc_matrix():
    N = 100
    model_prefix = "./model-0.relaxed_"
    rows, cols = (N, N)
    matrix = zeros([N, N], float)
    for i in range(rows):
        for j in range(cols):
            first_model_id = str(i + 1).zfill(4);
            second_model_id = str(j + 1).zfill(4);
            first_model_name = model_prefix + first_model_id + ".pdb"
            second_model_name = model_prefix + second_model_id + ".pdb"
            matrix[i][j] = fast_rmsd(first_model_name, second_model_name)

    np.savetxt('matrix.txt', matrix, fmt='%f')


def create_heatmap(mat, num_of_clusters):
    labels = pd.Series(KMeans(num_of_clusters).fit_predict(mat))
    order = np.argsort(labels)
    plt.figure()
    sns.heatmap(
        mat[order],
        cmap=sns.color_palette("rainbow"),
        vmin=0.5
    )
    pngName = "heatmap-" + str(num_of_clusters) + "png"
    plt.savefig(pngName)
    return labels

def create_scatter():
    df = read_csv("H3_modeling_100_scores.csv")
    df.plot()  # plots all columns against index
    df.plot(kind='scatter', x='rmsd', y='total_score')  # scatter plot
    df.plot(kind='density')  # estimate density function
    # df.plot(kind='hist')  # histogram

    from matplotlib import pyplot as plt
    import seaborn as sns
    import pandas as pd

    df = sns.load_dataset('H3_modeling_100_scores.csv')
    plt.figure()  # Push new figure on stack
    sns_plot = sns.pairplot(df, hue='species', size=2.5)
    plt.savefig('output.png')  # Save that figure

    # pngName = "scatter_question_5.png"
    # plt.savefig(pngName)

def part_2():
    k = [1, 5, 10]
    csvInfo = read_csv("H3_modeling_100_scores.csv", header=0)
    matrix = np.loadtxt('matrix.txt', dtype=float)
    for i in k:
        cluster_labels = create_heatmap(matrix, i)
        csvInfo['cluster'] = cluster_labels
        models = []
        for j in range(i):
            models_of_cluster_j = csvInfo[csvInfo.cluster == j]
            # sorted_df = models_of_cluster_j.sort_values(by=["rmsd"], ascending=False)
            model = models_of_cluster_j.loc[models_of_cluster_j.total_score == models_of_cluster_j.total_score.min()]
            # sorted_df[sorted_df.rmsd == sorted_df.rmsd.min()]
            models.append(model)

        models.sort(key=lambda x: x.rmsd.values[0])
        best_model = models[0]
        desc = best_model.description.values[0]
        rmsd = best_model.rmsd.values[0]
        score = best_model.total_score.values[0]
        print("For " + str(i) + " clusters: the best model is: " + desc + " with rmsd: " + str(
            rmsd) + " and score: " + str(score))


if __name__ == '__main__':
    create_scatter()





