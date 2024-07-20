import matplotlib.pylab as plt 
import plotly.graph_objects as go
import pandas as pd
from plotly.subplots import make_subplots
import plotly.express as px

def draw_barplot_three_sum(cor_list, fn_list, fp_list):
    data = {'Samples': [], 'Type': [], 'Percentage': []}
    for i in range(len(cor_list)):  
        data['Samples'].extend([i, i, i])
        data['Type'].extend(['Accuracy', 'False negative', 'False positive'])
        data['Percentage'].extend([cor_list[i]*100, fn_list[i]*100, fp_list[i]*100])

    fig = px.bar(data, x="Samples", y="Percentage", color="Type", title="Accuracy Test")
    fig.update_layout(
        xaxis=dict(title_font=dict(size=30), tickfont=dict(size=15)), \
        yaxis=dict(title_font=dict(size=30), tickfont=dict(size=20), tick0=0, dtick=10, showgrid=True, gridwidth=2, gridcolor='lightgrey'), \
        legend=dict(title_font=dict(size=25), font=dict(size=25))
    )
    fig.show()

def draw_barplot_acc(cor_list):
    acc_perc = [i*100 for i in cor_list]
    sample_id = [i for i in range(1,len(acc_perc)+1)]
    tick_values = [i * 10 for i in range(11)]
    fig = px.bar(x=sample_id, y=acc_perc)
    
    fig.update_layout(
        xaxis=dict(title_text='Samples',title_font=dict(size=20), tickfont=dict(size=14), linecolor='black'), \
        yaxis=dict(title_text='Percentage(%)', title_font=dict(size=20), tickfont=dict(size=14), tick0=0, dtick=5, linecolor='black', showgrid=True, gridwidth=2, gridcolor='Silver', range=[0, 101]), 
        plot_bgcolor='white',
        width=800,
        height=500
    )
    fig.show()

def draw_barplot_fn_fp(fn, fp):
    data = {'Samples': [], 'Type': [], 'Percentage': []}
    for i in range(len(fn)):  
        data['Samples'].extend([i, i])
        data['Type'].extend(['False negative', 'False positive'])
        data['Percentage'].extend([fn[i]*100, fp[i]*100])

    fig = px.bar(data, x="Samples", y="Percentage", color="Type", title="Accuracy Test")
    fig.update_layout(
        xaxis=dict(title_font=dict(size=30), tickfont=dict(size=15)), \
        yaxis=dict(title_font=dict(size=30), tickfont=dict(size=20), tick0=0, dtick=1, showgrid=True, linecolor='black', gridwidth=2, gridcolor='Silver'), \
        legend=dict(title_font=dict(size=25), font=dict(size=25)), \
        plot_bgcolor='white', \
        width=800, \
        height=400
    )
    fig.show()

def draw_boxplot_fn_fp(fn, fp):
    headers = ['False negative', 'False positive']
    fn_perc = [i*100 for i in fn]
    fp_perc = [i*100 for i in fp]
    df = pd.DataFrame(dict(zip(headers, [fn_perc, fp_perc])))
    fig = go.Figure()

    for col in df:
        fig.add_trace(go.Box(y=df[col].values, name=df[col].name, line=dict(width=3)))
    fig.update_layout(
        xaxis=dict(tickfont=dict(size=20), linecolor='black'),
        yaxis=dict(title_text='Percentage(%)', title_font=dict(size=20), tickfont=dict(size=14), tick0=0, linecolor='black', dtick=1, showgrid=True, gridwidth=2, gridcolor='Silver', range=[0, 7.1]), 
        showlegend=False, \
        plot_bgcolor='white', \
        width=800, \
        height=400
        )
    fig.show()

if __name__ == "__main__":
    cor_list = [0.9856115107913669, 0.9487179487179487, 0.975, 0.9523809523809523, 0.9310344827586207, 1.0, 0.9772727272727273, 0.9746835443037974, 0.9891304347826086, 1.0, 0.9701492537313433, 1.0, 0.9876543209876543, 0.93, 0.9591836734693877, 0.9732142857142857, 0.9705882352941176, 0.9770114942528736, 1.0, 0.9791666666666666, 0.984375, 0.9841269841269841, 0.9367088607594937, 1.0, 0.9787234042553191, 1.0, 1.0, 0.9714285714285714, 0.92, 0.9393939393939394, 0.8918918918918919, 0.96, 0.95, 0.9411764705882353, 0.9591836734693877, 0.9767441860465116, 0.9743589743589743, 0.9772727272727273, 0.9591836734693877, 0.9555555555555556, 1.0, 1.0, 1.0, 0.9076923076923077, 0.9452054794520548, 0.9358974358974359, 0.9298245614035088, 0.918918918918919, 0.9795918367346939, 0.98, 0.9807692307692307, 0.9827586206896551, 0.9655172413793104, 1.0, 1.0, 1.0, 1.0, 0.95, 0.9298245614035088, 0.9259259259259259, 0.9583333333333334, 0.9387755102040817]
    fal_neg_list = [0.0, 0.008547008547008548, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014925373134328358, 0.0, 0.0, 0.02, 0.013605442176870748, 0.008928571428571428, 0.00980392156862745, 0.011494252873563218, 0.0, 0.0, 0.0, 0.0, 0.02531645569620253, 0.0, 0.02127659574468085, 0.0, 0.0, 0.0, 0.02, 0.030303030303030304, 0.05405405405405406, 0.02, 0.016666666666666666, 0.0196078431372549, 0.0, 0.0, 0.0, 0.0, 0.02040816326530612, 0.022222222222222223, 0.0, 0.0, 0.0, 0.046153846153846156, 0.0273972602739726, 0.038461538461538464, 0.03508771929824561, 0.04054054054054054, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.03508771929824561, 0.037037037037037035, 0.020833333333333332, 0.02040816326530612]
    fal_pos_list = [0.014388489208633094, 0.042735042735042736, 0.025, 0.047619047619047616, 0.06896551724137931, 0.0, 0.022727272727272728, 0.02531645569620253, 0.010869565217391304, 0.0, 0.014925373134328358, 0.0, 0.012345679012345678, 0.05, 0.027210884353741496, 0.017857142857142856, 0.0196078431372549, 0.011494252873563218, 0.0, 0.020833333333333332, 0.015625, 0.015873015873015872, 0.0379746835443038, 0.0, 0.0, 0.0, 0.0, 0.02857142857142857, 0.06, 0.030303030303030304, 0.05405405405405406, 0.02, 0.03333333333333333, 0.0392156862745098, 0.04081632653061224, 0.023255813953488372, 0.02564102564102564, 0.022727272727272728, 0.02040816326530612, 0.022222222222222223, 0.0, 0.0, 0.0, 0.046153846153846156, 0.0273972602739726, 0.02564102564102564, 0.03508771929824561, 0.04054054054054054, 0.02040816326530612, 0.02, 0.019230769230769232, 0.017241379310344827, 0.034482758620689655, 0.0, 0.0, 0.0, 0.0, 0.025, 0.03508771929824561, 0.037037037037037035, 0.020833333333333332, 0.04081632653061224]
    # draw_barplot_acc(cor_list=cor_list)
    # draw_barplot_fn_fp(fn=fal_neg_list, fp=fal_pos_list)
    # draw_boxplot_fn_fp(fn=fal_neg_list, fp=fal_pos_list)