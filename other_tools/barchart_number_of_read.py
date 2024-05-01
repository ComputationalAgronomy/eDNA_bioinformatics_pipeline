import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

data = {'Step': ['QC & <br>Demultiplexed', 'Paired-end<br>Merging', 'Primer Trimming<br>&Length Filtering', 'Dereplication', 'Denoising', 'Blast'],
        'Count': [25632654, 11797855, 11622521, 515739, 4246, 3525],
        'Type': ['sequence reads', 'sequence reads', 'sequence reads', 'unique sequences', 'ZOTUs', 'Assigned species']
}
df = pd.DataFrame(data)
fig = px.bar(data, x="Step", y="Count", color="Type", text='Count')
fig.update_traces(textposition='outside', textfont=dict(color='black'), marker=dict(line=dict(color='black', width=1)))
fig.update_layout(
    xaxis=dict(title_font=dict(size=13), tickfont=dict(size=11), linecolor='black'), \
    yaxis=dict(title_font=dict(size=13), tickfont=dict(size=13), type='log', tickformat='.0e', dtick=1, range=[0, 13^5], showgrid=True, gridwidth=1, gridcolor='Silver', linecolor='black'), \
    legend=dict(title_font=dict(size=11), font=dict(size=11)), \
    font=dict(size=13),
    width=800,
    height=450,
    plot_bgcolor='white'
)
plt.ylim([0, 50000000])
fig.show()