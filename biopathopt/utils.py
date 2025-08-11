from typing import Tuple, Dict, List, Union, Optional, Any
import pandas as pd
import re
import pickle
import gzip
import json
import copy
from rapidfuzz import process, fuzz


def merge_annot_dicts(input_parent_dict: dict, child_dict: dict) -> dict:
    """
    Merges dictionary child_dict into dictionary input_parent_dict by:
    - Adding elements from child_dict that do not exist in input_parent_dict.
    - If an element is a list in both A and child_dict, adding only unique elements from child_dict's list to input_parent_dict's list.
    - Ignoring elements in child_dict that are lists if input_parent_dict does not have them as lists.
    - Ignoring elements that are dictionaries in child_dict.

    Args:
        input_parent_dict (dict): The base dictionary.
        child_dict (dict): The dictionary to merge into A.

    Returns:
        dict: The updated dictionary input_parent_dict.

    Example:
        >>> A = {"x": 1, "y": "hello", "z": [1, 2, 3]}
        >>> B = {"y": "world", "z": [3, 4, 5], "a": 100, "nested": {"key": "value"}}
        >>> merge_dicts_with_list_update(A, B)
        {'x': 1, 'y': 'hello', 'z': [1, 2, 3, 4, 5], 'a': 100}
    """
    parent_dict = copy.deepcopy(input_parent_dict)
    for key, value in child_dict.items():
        if isinstance(value, dict):  # Ignore dictionary values
            continue
        if key in parent_dict:
            if isinstance(parent_dict[key], list) and isinstance(value, list):
                parent_dict[key].extend([item for item in value if item not in parent_dict[key]])
        elif not isinstance(value, list):
            parent_dict[key] = value
    return parent_dict

# Write compressed dict to JSON
def write_compressed_json(data: dict, filename: str) -> None:
    with gzip.open(filename, 'wt', encoding='utf-8') as f:
        json.dump(data, f)

# Read compressed dict from JSON
def read_compressed_json(filename: str) -> dict:
    with gzip.open(filename, 'rt', encoding='utf-8') as f:
        return json.load(f)

def annotation_contains_value(data: Any, target: str) -> bool:
    """
    Recursively checks if target string exists anywhere in the nested structure.
    
    Args:
        data: The nested structure (dict, list, str, etc.)
        target: The string to search for.
        
    Returns:
        True if the string is found, otherwise False.
    """
    if isinstance(data, dict):
        return any(annotation_contains_value(v, target) for v in data.values())
    elif isinstance(data, list):
        return any(annotation_contains_value(i, target) for i in data)
    elif isinstance(data, str):
        return data == target
    return False

def inchikey_layer_extract(inchikey, layers=2):
    if layers==2:
        return '-'.join(inchikey.split('-')[:-1])
    elif layers==1:
        return '-'.join(inchikey.split('-')[:-2])
    else:
        return inchikey


def fuzzy_dict_lookup(query: str, data: dict, threshold: float = 90.0):
    """
    Performs a fuzzy search on the keys of a dictionary using Rapidfuzz.

    Args:
        query (str): The string to match.
        data (dict): Dictionary to search over.
        threshold (float): Minimum score (0–100) to consider a match.

    Returns:
        The value of the best matching key, or None if no match meets the threshold.
    """
    match = process.extractOne(query, data.keys(), scorer=fuzz.ratio)
    if match and match[1] >= threshold:
        return data[match[0]]
    return None


def load_pickle(file_name):
    """Load a pickle file
    """
    with open(file_name, 'rb') as f:
        return pickle.load(f)

def parse_mnxr_equation(reac_str: str) -> Tuple[List[List[str]], List[List[str]]]:
    """Parses a metabolic network equation string from MetaNetX into reactants and products.

    This function takes a metabolic network equation string and parses it into 
    lists of reactants and products. The input string should be in the form 
    'reactant1 + reactant2 = product1 + product2', where each component may 
    contain quantities and compartment information separated by spaces and '@' symbols.

    Args:
        reac_str (str): The metabolic network equation string.

    Returns:
        Tuple[List[List[str]], List[List[str]]]: A tuple containing two lists:
            - The first list contains lists of parsed reactants.
            - The second list contains lists of parsed products.
    """
    # Split the equation string into reactants and products
    reactants_str, products_str = reac_str.split('=')

    # Split the reactants and products by the '+' symbol and strip any surrounding whitespace
    reactants = [r.strip() for r in reactants_str.split('+')]
    products = [p.strip() for p in products_str.split('+')]

    # Further split each reactant and product by spaces
    reactants = [r.split(' ') for r in reactants]
    products = [p.split(' ') for p in products]

    # Remove compartment information and other annotations following the '@' symbol
    reactants = [[item.split('@')[0] for item in r] for r in reactants]
    products = [[item.split('@')[0] for item in p] for p in products]

    return reactants, products

def pd_find_with_numpy(
    df: pd.DataFrame, 
    search_column_id: str, 
    search_value: Union[str, int, float], 
    return_column_id: Optional[str] = None
) -> Union[pd.DataFrame, str]:
    """Searches for a value in a DataFrame using NumPy for performance and returns the corresponding value or row.

    This function searches a given column in a DataFrame for a specific value using NumPy for optimized performance. 
    If a `return_column_id` is provided, it returns the corresponding value from that column; otherwise, it returns 
    the entire row that matches the search value.

    Args:
        df (pd.DataFrame): The DataFrame to search within.
        search_column_id (str): The column in which to search for `search_value`.
        search_value (Union[str, int, float]): The value to search for in `search_column_id`.
        return_column_id (Optional[str], optional): The column from which to return the value that corresponds 
            to the found `search_value`. If None, returns the entire row. Defaults to None.

    Returns:
        Union[pd.DataFrame, str]: If `return_column_id` is provided, returns the corresponding value from that column. 
        Otherwise, returns a DataFrame with the matching rows.

    Raises:
        KeyError: If the search value is not found, or the found value is NaN, 'nan', or empty.

    Example:
        >>> data = {'ID': ['A1', 'A2', 'A3'], 'Value': [10, 20, 30]}
        >>> df = pd.DataFrame(data)
        >>> pd_find_with_numpy(df, 'ID', 'A2', 'Value')
        20
        >>> pd_find_with_numpy(df, 'ID', 'A4')
        KeyError: ''

    """
    try:
        if return_column_id:
            # Convert the DataFrame columns to NumPy arrays for efficient searching
            to_ret = df[return_column_id].to_numpy()[df[search_column_id].to_numpy() == search_value].item()

            # Check if the found value is NaN, 'nan', or empty, and raise KeyError if so
            if pd.isna(to_ret):
                raise KeyError
            if to_ret == 'nan':
                raise KeyError
            if not to_ret:
                raise KeyError

            return to_ret
        else:
            # If no return column is specified, return the entire row(s) matching the search value
            return df[df[search_column_id].to_numpy() == search_value]
    except ValueError:
        # Raise KeyError if the search value is not found
        raise KeyError


def group_ids(df: pd.DataFrame, uid_column: str, data_columns: str) -> Dict[str, List[str]]:
    """Groups IDs by their corresponding column values, keeping duplicates.

    This function takes a DataFrame and groups the IDs by the values in a specified data column. 
    If multiple IDs share the same value in the data column, they are grouped together in a list.

    Args:
        df (pd.DataFrame): A DataFrame containing the relevant columns.
        uid_column (str): The name of the column containing the unique IDs.
        data_column (str): The name of the column containing the data to group by.

    Returns:
        Dict[str, List[str]]: A dictionary where each key is a value from the data column, and 
        the corresponding value is a list of IDs that share that value.

    Raises:
        IndexError: If either `uid_column` or `data_column` is not present in the DataFrame.
    """
    a = {}
    if uid_column not in df.columns:
        raise IndexError(f'No uid_column in the dataframe input')
    if data_columns not in df.columns:
        raise IndexError(f'No data_columns in the dataframe input')
    # Iterate over the 'ID' and 'InChIKey_2' columns, pairing them together
    for uid, inchikey in zip(df[uid_column].to_list(), df[data_columns].to_list()):
        if inchikey in a:
            a[inchikey].append(uid)  # Append the ID to the list if InChIKey_2 already exists
        else:
            a[inchikey] = [uid]  # Create a new list with the ID if InChIKey_2 is not yet in the dictionary
    return a


def construct_reaction_string(
    reactants: List[Tuple[str, str]],
    products: List[Tuple[str, str]],
    inchikey_levels: int = 3
) -> str:
    """Constructs a reaction string from reactants and products with truncated InChIKeys.

    This function generates a reaction string from lists of reactants and products.
    Each reactant and product is represented by a coefficient and an InChIKey.
    The InChIKey is truncated to the specified number of levels.

    Args:
        reactants (List[Tuple[str, str]]): A list of tuples where each tuple contains
            the coefficient (str) and the InChIKey (str) for a reactant.
        products (List[Tuple[str, str]]): A list of tuples where each tuple contains
            the coefficient (str) and the InChIKey (str) for a product.
        inchikey_levels (int, optional): The number of levels of the InChIKey to include
            in the output string. Defaults to 3.

    Returns:
        str: The constructed reaction string in the format 'reactants = products'.

    Example:
        >>> reactants = [("2", "ABCDEF-ABCDEF-N"), ("1", "GHIJKL-GHIJKL-F")]
        >>> products = [("1", "MNOPQR-MNOPQR-N")]
        >>> construct_reaction_string(reactants, products, inchikey_levels=2)
        '2 ABCDEF-ABCDEF + 1 GHIJKL-GHIJKL = 1 MNOPQR-MNOPQR'
    """

    # Construct the reaction string for reactants
    reacts_str = ''
    for coefficient, inchikey in reactants:
        # Truncate the InChIKey to the specified number of levels
        truncated_inchikey = '-'.join(inchikey.split('-')[:inchikey_levels])
        reacts_str += f"{coefficient} {truncated_inchikey} + "
    reacts_str = reacts_str[:-3]  # Remove the trailing ' + '

    # Construct the reaction string for products
    prods_str = ''
    for coefficient, inchikey in products:
        # Truncate the InChIKey to the specified number of levels
        truncated_inchikey = '-'.join(inchikey.split('-')[:inchikey_levels])
        prods_str += f"{coefficient} {truncated_inchikey} + "
    prods_str = prods_str[:-3]  # Remove the trailing ' + '

    # Combine reactants and products into the final equation string
    return f"{reacts_str} = {prods_str}"



def check_value_in_annotation(
        value: str, 
        annotation: Dict[str, Union[str, List[str]]],
        keys: Optional[List[str]] = None,
    ) -> bool:
    """
    Checks if a given value is contained in the values of the annotation for one or more specified keys.
    If no keys are provided, it will check all keys in the annotation.
    The annotation can have values as either strings or lists of strings.

    Args:
        value (str): The value to search for in the annotation.
        keys (Optional[List[str]]): The list of keys where the value should be searched. If None, all keys are checked.
        annotation (Dict[str, Union[str, List[str]]]): The annotation to search within.

    Returns:
        bool: True if the value is found in any of the specified keys (or any key if no keys are specified), False otherwise.
    """
    
    # If no specific keys are provided, use all keys in the annotation
    if keys is None:
        keys = list(annotation.keys())
    if isinstance(keys, str):
        keys = [keys]
    
    # Iterate over the provided keys
    for key in keys:
        # Ensure the key exists in the annotation to avoid KeyError
        if key in annotation:
            dict_value = annotation[key]
            
            # If the value in the annotation is a list
            if isinstance(dict_value, list):
                # Check if the given value is in the list
                if value in dict_value:
                    return True
            
            # If the value in the annotation is a string
            elif isinstance(dict_value, str):
                # Check if the given value matches the annotation value
                if value == dict_value:
                    return True
    
    # If the value is not found in any specified keys, return False
    return False


def extract_genes_from_gpr(expression: str) -> List[List[str]]:
    """
    Extracts genes from an expression, grouping genes in parentheses
    and separating groups connected by 'or'.

    Args:
        expression (str): The gene expression string to parse.

    Returns:
        List[List[str]]: A list of gene groups. Each group is a list of genes,
        with groups separated by 'or' as distinct sublists.

    Example:
        ```python
        expression1 = "(b4244 and b4245) or b4245"
        expression2 = "b4069"
        expression3 = "b4069 or b0335"

        print(extract_genes(expression1))  # [['b4244', 'b4245'], ['b4245']]
        print(extract_genes(expression2))  # [['b4069']]
        print(extract_genes(expression3))  # [['b4069'], ['b0335']]
        ```
    """
    # Remove 'and' since it is irrelevant for splitting and grouping
    expression = re.sub(r'\band\b', '', expression)

    # Regular expression to capture groups of genes in parentheses or standalone genes
    pattern = r"\(([^)]+)\)|(\b\w+\b)"
    matches = re.findall(pattern, expression)

    result = []
    for group in matches:
        if group[0]:  # Genes inside parentheses
            genes = [gene.strip() for gene in re.split(r'\s+', group[0]) if gene.strip()]
            result.append(genes)
        elif group[1] and group[1] != 'or':  # Standalone genes excluding 'or'
            result.append([group[1]])

    return result


'''
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs

def convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
    """Convert chemical depiction to others type of depictions

    Usage example:
     - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
     - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})

    :param idepic: Input string
    :param itype: The type of input
    :param otype: Type of output. Valid options: inchi, smiles, inchikey

    :type idepic: str
    :type itype: str
    :type otype: dict

    :rtype: dict
    :return: Dictionnary of results
    """
    # Import (if needed)
    self.logger.debug('input: '+str(idepic))
    self.logger.debug('itype: '+str(itype))
    self.logger.debug('otype: '+str(otype))
    if itype == 'smiles':
        rdmol = MolFromSmiles(idepic, sanitize=True)
    elif itype == 'inchi':
        rdmol = MolFromInchi(idepic, sanitize=True)
    else:
        raise NotImplementedError('"{}" is not a valid input type'.format(itype))
    if rdmol is None:  # Check imprt
        raise self.DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
    self.logger.debug('Sanitised the input')
    # Export
    odepic = dict()
    for item in otype:
        if item == 'smiles':
            odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
        elif item == 'inchi':
            odepic[item] = MolToInchi(rdmol)
        elif item == 'inchikey':
            odepic[item] = MolToInchiKey(rdmol)
        else:
            raise NotImplementedError('"{}" is not a valid output type'.format(otype))
    self.logger.debug('Exported the output')
    return odepic
'''
