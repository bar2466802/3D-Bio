#################################################################
# FILE : ex11.py
# WRITER : Eldar Michal, michalel
# EXERCISE : intro2cs2 ex11 2021
# DESCRIPTION: A program that diagnose person's illness.
#################################################################
import copy
from itertools import combinations

NO_RECORDS = "No records! Got an empty list."

BAD_TYPE = "Wrong type in symptoms or records."

BAD_ARGUMENT = "depth must be smaller than symptom len and non negative, and" \
               " symptoms must be a unique list"


class Node:
    def __init__(self, data, positive_child=None, negative_child=None):
        self.data = data
        self.positive_child = positive_child
        self.negative_child = negative_child

    def is_leaf(self):
        """ checks if node is a leaf. returns True for leaf, else False """
        if not self.positive_child and not self.negative_child:
            return True
        return False


class Record:
    def __init__(self, illness, symptoms):
        self.illness = illness
        self.symptoms = symptoms


def parse_data(filepath):
    with open(filepath) as data_file:
        records = []
        for line in data_file:
            words = line.strip().split()
            records.append(Record(words[0], words[1:]))
        return records


class Diagnoser:
    def __init__(self, root: Node):
        self.root = root

    def diagnose_helper(self, curr_root, symptoms):
        """
        recursive function, gets root tree and symptoms list, and returns the
        illness according to the tree.
        :param curr_root: Node, current tree root,
        :param symptoms: list of symptoms (strings)
        :return: string, diagnosed illness.
        """
        if curr_root.is_leaf():
            return curr_root.data
        if curr_root.data in symptoms:
            return self.diagnose_helper(curr_root.positive_child, symptoms)
        else:
            return self.diagnose_helper(curr_root.negative_child, symptoms)

    def diagnose(self, symptoms):
        """
        gets symptoms list, and returns the illness according to the tree.
        :param symptoms: list of symptoms (strings)
        :return: string, diagnosed illness.
        """
        temp_root = copy.copy(self.root)
        answer = self.diagnose_helper(temp_root, symptoms)
        return answer

    def calculate_success_rate(self, records):
        """
        calculates success rate of records diagnosis.
        :param records: list or records.
        :return: float (represents percentage)
        """
        if not len(records):
            raise ValueError(NO_RECORDS)
        successful_diagnosis = 0
        for record in records:
            if record.illness == self.diagnose(record.symptoms):
                successful_diagnosis += 1
        return successful_diagnosis / len(records)

    def all_illnesses_helper(self, curr_root, leaves):
        """
        recursively returns list of all the illnesses.
        """
        if curr_root.is_leaf():
            illness = curr_root.data
            return leaves + [illness] if illness else leaves
        else:
            left = self.all_illnesses_helper(curr_root.positive_child, leaves)
            right = self.all_illnesses_helper(curr_root.negative_child, leaves)
            return right + left

    def all_illnesses(self):
        """
        gets all illnesses (leaves) in the diagnoser tree.
        :return: list of all illnesses without duplications, sorted by
        frequency of the illness.
        """
        curr_root = copy.copy(self.root)
        leaves = self.all_illnesses_helper(curr_root, [])
        return uniq_by_counts(leaves)

    def paths_to_illness_helper(self, curr_root, illness, path, paths):
        """
        recursive helper function for paths_to_illness.
        """
        if curr_root.is_leaf():
            if curr_root.data == illness:
                return [path]
            else:
                return []
        pos_path = self.paths_to_illness_helper(curr_root.positive_child,
                                                illness, path + [True], paths)
        neg_path = self.paths_to_illness_helper(curr_root.negative_child,
                                                illness, path + [False], paths)
        return pos_path + neg_path

    def paths_to_illness(self, illness):
        """
        returns list of lists with paths to illness
        :param illness: illnes to calculate path for.
        :return: list of lists with boolean args.
        """
        curr_root = copy.copy(self.root)
        return self.paths_to_illness_helper(curr_root, illness, [], [])

    def get_all_symptoms(self, node):
        """
        returns list of all symptoms in the diagnoser, starts from node.
        :param node: node to start scan symptoms from.
        :return: list
        """
        if node.is_leaf():
            return []
        return [node.data] + self.get_all_symptoms(node.positive_child) + self.get_all_symptoms(node.negative_child)

    def minimize(self, remove_empty=False):
        """
        minimize tree.
        :param remove_empty:
        :return:
        """
        if not self.root.positive_child and not self.root.negative_child:  # Leaf
            return

        pos_tree = Diagnoser(self.root.positive_child)
        pos_tree.minimize(remove_empty)
        neg_tree = Diagnoser(self.root.negative_child)
        neg_tree.minimize(remove_empty)
        self.root.positive_child = pos_tree.root
        self.root.negative_child = neg_tree.root
        symptoms = self.get_all_symptoms(self.root)

        if remove_empty and not pos_tree.root.data and neg_tree.root.data:
            self.root = neg_tree.root
            return

        if remove_empty and pos_tree.root.data and not neg_tree.root.data:
            self.root = pos_tree.root
            return

        for i in range(1, len(symptoms) + 1):
            for syms in combinations(symptoms, i):
                if pos_tree.diagnose(syms) != neg_tree.diagnose(syms):
                    return
        self.root = pos_tree.root
        return


def uniq_by_counts(lst):
    """
    unique a list and returns it with an order that represents the frequency of
     each value in the list, starts with the mor frequent.
    """
    counter_dict = {}
    for i in lst:
        if i in counter_dict:
            counter_dict[i] += 1
        else:
            counter_dict[i] = 1
    return list(dict(sorted(counter_dict.items(), key=lambda item: item[1],
                            reverse=True)))


def pos_neg_records(records, symptom):
    """
    Gets symptom and records and returns tuple of:
    positive list that contains all the records that has the symptom
    megative list that contains all the records that doesnt have the symptom
    :param records: list of Record objects
    :param symptom: str
    """
    pos_list = []
    neg_list = []
    for record in records:
        if symptom in record.symptoms:
            pos_list.append(record)
        else:
            neg_list.append(record)
    return pos_list, neg_list


def record_to_illness(records):
    """
    Gets records lust and returns the most common illness
    """
    if not records:
        return None
    elif len(records) == 1:
        return records[0].illness
    else:
        illnesses = [record.illness for record in records]
        return uniq_by_counts(illnesses)[0]


def build_tree_helper(records, symptoms):
    """ recursive helper function for build_tree """
    if not symptoms:
        illness = record_to_illness(records)
        return Node(illness)
    if type(symptoms[0]) != str or (len(records) and type(records[0]) != Record):
        raise TypeError(BAD_TYPE)
    else:
        root = Node(symptoms[0])
        pos_records, neg_records = pos_neg_records(records, root.data)
        root.positive_child = build_tree_helper(pos_records, symptoms[1:])
        root.negative_child = build_tree_helper(neg_records, symptoms[1:])
    return root


def build_tree(records, symptoms):
    """
    using build_tree_helper function, build a tree and diagnoser from it,
    from records and symptoms.
    :param records: list of records to build tree from
    :param symptoms: list of symptoms to build tree from.
    :return: Diagnoser object with the tree.
    """
    tree = build_tree_helper(records, symptoms)
    return Diagnoser(tree)


def optimal_tree(records, symptoms, depth):
    """
    builds a tree with the maximal diagnosis success
    :param records: list of records
    :param symptoms: list of symptoms
    :param depth: depth of the tree
    :return: Diagnoser object
    """
    if not (0 <= depth <= len(symptoms)) or len(set(symptoms)) != len(symptoms):
        raise ValueError(BAD_ARGUMENT)
    best_success = 0
    best_diagnoser = None
    symptoms_combs = combinations(symptoms, depth)
    for symptoms_comb in symptoms_combs:
        tmp_diagnoser = build_tree(records, symptoms_comb)
        tmp_success = tmp_diagnoser.calculate_success_rate(records)
        if tmp_success > best_success:
            best_diagnoser = tmp_diagnoser
            best_success = tmp_success
    return best_diagnoser


if __name__ == "__main__":
    pass
