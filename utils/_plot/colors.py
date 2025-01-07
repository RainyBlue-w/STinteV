from stintev.config import Colors
from typing import List, Dict, Union

def assign_colors(
    categories: Union[List[str], str],  
    cur_global_cmap: Dict,
    color_list = Colors.COLORS_60,
) -> Dict:
    """
    Assign colors to categories

    Params:
    categories: List[str], list of categories to assign colors
    cur_global_cmap: Dict, current global color map
        example: {'i': 0, 'cmap': {'A': '#FF0000', 'B': '#008000'}}
    color_list: List[str], full list of colors to select from
        
    Return:
    Dict: updated `cur_global_cmap`

    """
    color_dict = {}
    i = 0
    cur_cmap = cur_global_cmap['cmap']
    curIndex = cur_global_cmap['i']
    if isinstance(categories, list):
        for c in categories:
            if c not in cur_cmap:
                color_dict[c] = color_list[(i+curIndex) % len(color_list)]
                i += 1
    else:
        color_dict[categories] = color_list[curIndex % len(color_list)]
        
    cur_global_cmap['i'] = (curIndex+i) % len(color_list)
    cur_global_cmap['cmap'].update(color_dict)
    
    return cur_global_cmap
