import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def normalize_abundance(abundance_dict: dict[str, float]) -> dict[str, float]:
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
    fig.update_xaxes(tickmode='linear')
    fig.update_layout(
        xaxis_title="Sample ID",
        yaxis_title="Percentage (%)",
        legend={
            "x": 1.05,
            "y": 1,
            "traceorder": 'normal',
            "orientation": 'h'
        }
    )
    return fig