import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def normalize_abundance(abundance_dict: dict[str, int]) -> dict[str, float]:
    """
    Normalize the abundance values in the dictionary to percentages.

    :param abundance_dict: A dictionary with unit names and their abundance.
    :returns: A dictionary with unit names and their normalized abundance in percentages.
    """
    total_size = sum(abundance_dict.values())
    norm_abundance = {key: value/total_size * 100 for key, value in abundance_dict.items()}
    return norm_abundance

def list_union(lists_to_union: list[list[str]]) -> list[str]:
    """
    Return the union of multiple lists as a sorted list of unique unit names.

    :param lists_to_union: A list of lists to be unioned.
    :return: A sorted list of unique elements.
    """
    uniq_list = list(set().union(*lists_to_union))
    uniq_list.sort()
    return uniq_list

def create_barchart_fig(data: pd.DataFrame) -> go.Figure:
    """
    Create a stacked bar chart figure using Plotly.

    :param data: DataFrame containing the data to be plotted.
    :return: A Plotly Figure object representing the bar chart.
    """
    fig = px.bar(
        data,
        barmode='stack',
        labels={'value': 'Percentage (%)'},
        color_discrete_sequence=px.colors.qualitative.Pastel
    )
    fig.update_xaxes(
        tickmode='linear',
        title=dict(
            text="Sample ID",
            font=dict(size=20)
            ),
        tickfont=dict(size=18)
    )
    fig.update_yaxes(
        title=dict(
            text="Percentage (%)",
            font=dict(size=20)
        ),
        tickfont=dict(size=18)
    )
    fig.update_layout(
        legend={
            "x": 1.05,
            "y": 1,
            "traceorder": 'normal',
            "orientation": 'h',
            "font": dict(size=15)
        },
    )
    return fig

    ## https://stackoverflow.com/questions/44309507/stacked-bar-plot-using-matplotlib
    # a way to sort stacked BarChart
    # plotdata = plotdata.transpose()
    # fig, ax = plt.subplots()
    # x = plotdata.index
    # indexes = np.argsort(plotdata.values).T
    # heights = np.sort(plotdata.values).T
    # order = -1
    # bottoms = heights[::order].cumsum(axis=0)
    # bottoms = np.insert(bottoms, 0, np.zeros(len(bottoms[0])), axis=0)
    # colormap = plt.get_cmap('rainbow')
    # colors = px.colors.qualitative.Pastel
    # num_colors = len(colors)
    # mpp_colors = {col: colors[i % num_colors] for i, col in enumerate(plotdata.columns)}

    # fig = go.Figure()
    
    # # # Plot each segment of the stacked bar
    # for i, (idxs, vals) in enumerate(list(zip(indexes, heights))[::order]):
    #     mps = np.take(np.array(plotdata.columns), idxs)
    #     for j, m in enumerate(mps):

    #         fig.add_trace(go.Bar(
    #             x=x,
    #             y=[vals[j]] * len(x),
    #             name=m,
    #             marker=dict(color=mpp_colors[m]),
    #             offsetgroup=i,
    #             base=bottoms[i][j]
    #         ))

    # # Update layout
    # fig.update_layout(
    #     barmode='stack',
    #     xaxis_title="Sample ID",
    #     yaxis_title="Percentage (%)",
    #     legend=dict(
    #         x=1.05,
    #         y=1,
    #         traceorder='normal',
    #         orientation='h',
    #         font=dict(size=15)
    #     ),
    # )
    # fig.update_xaxes(tickmode='linear', title=dict(font=dict(size=20)), tickfont=dict(size=18))
    # fig.update_yaxes(title=dict(font=dict(size=20)), tickfont=dict(size=18))
    # fig.show()

