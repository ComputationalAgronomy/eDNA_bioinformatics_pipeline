import numpy as np
import scipy.stats as st 
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = round(np.mean(a), 4), round(st.sem(a), 4 )
    h = se * st.t.ppf((1 + confidence) / 2., n-1)
    print(f'Mean: {m} Std Err: {se}\nC.I.: ({round(m-h, 4)}, {round(m+h, 4)})')

def draw_boxplot(alignment, mutation, indel):
    headers = ['alignment length (bp)', 'mismatch (bp)', 'indel (bp)']
    df = pd.DataFrame(dict(zip(headers, [alignment, mutation, indel])))
    fig = make_subplots(rows=1, cols=len(headers))
    for i, header in enumerate(headers):
        fig.add_trace(go.Box(y=df[header], name=header, boxpoints='outliers'), row=1, col=i+1)
    fig.show()


if __name__ == '__main__':

    alignment = [168, 142, 177, 168, 168, 169, 169, 170, 168, 152, 158, 173, 173, 170, 171, 169, 169, 173, 173, 173, 173, 173]
    mutation = [2, 33, 0, 5, 1, 1, 8, 0, 0, 16, 12, 0, 1, 5, 5, 1, 1, 0, 0, 0, 0, 0]
    indel = [0, 5, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]

    draw_boxplot(alignment=alignment, mutation=mutation, indel=indel)