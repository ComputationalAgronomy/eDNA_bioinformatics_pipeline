import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def list_union(lists_to_union: list[list[str]]) -> list[str]:
    """
    Return the union of multiple lists as a sorted list of unique unit names.

    :param lists_to_union: A list of lists to be unioned.
    :return: A sorted list of unique elements.
    """
    uniq_list = list(set().union(*lists_to_union))
    uniq_list.sort()
    return uniq_list

def get_color_hex(n: int, cmap: str = "rainbow") -> list[str]:
    color_key = plt.get_cmap(cmap)(np.linspace(0, 1, n))
    color_hex = [matplotlib.colors.to_hex(color_key[i]) for i in range(n)]
    print(f"> Color key: {" ".join(color_hex)}")
    return color_hex

def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip('#')
    lv = len(value)
    rgb = tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
    print(f"> '{value}' corresponds to red:{rgb[0]} green:{rgb[1]} blue:{rgb[2]}.")
    return rgb