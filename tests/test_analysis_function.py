import pytest

from analysis_pipeline import analysis_function

def test_normalize_abundance():
    # Test case 1: Empty dictionary
    abundance_dict = {}
    expected_result = {}
    assert analysis_function.normalize_abundance(abundance_dict) == expected_result

    # Test case 2: Dictionary with positive values
    abundance_dict = {'A': 10, 'B': 20, 'C': 30}
    expected_result = {'A': 16.666666666666664, 'B': 33.33333333333333, 'C': 50.0}
    assert analysis_function.normalize_abundance(abundance_dict) == expected_result

    # Test case 3: Dictionary with zero values
    # abundance_dict = {'A': 0, 'B': 0, 'C': 0}
    # expected_result = {'A': 0.0, 'B': 0.0, 'C': 0.0}
    # assert analysis_function.normalize_abundance(abundance_dict) == expected_result


    # Test case 5: Dictionary with mixed positive and negative values
    abundance_dict = {'A': 10.10, 'B': 20.20, 'C': 30.30}
    expected_result = {'A': 16.666666666666664, 'B': 33.33333333333333, 'C': 50.0}
    assert analysis_function.normalize_abundance(abundance_dict) == expected_result




def test_list_union():
    # Test case 1: Empty list
    lists_to_union = []
    expected_output = []
    assert analysis_function.list_union(lists_to_union) == expected_output

    # Test case 2: Single list
    lists_to_union = [[1, 2, 3]]
    expected_output = [1, 2, 3]
    assert analysis_function.list_union(lists_to_union) == expected_output

    # Test case 3: Multiple lists with duplicates
    lists_to_union = [[1, 2, 3], [3, 4, 5], [5, 6, 1]]
    expected_output = [1, 2, 3, 4, 5, 6]
    assert analysis_function.list_union(lists_to_union) == expected_output

    # Test case 4: Lists with different lengths
    lists_to_union = [[1, 2], [3, 4, 5], [], [6]]
    expected_output = [1, 2, 3, 4, 5, 6]
    assert analysis_function.list_union(lists_to_union) == expected_output
