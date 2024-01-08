import pandas as pd
import copy
from stack import Stack
from flask import Flask, request, render_template, redirect, url_for

EPSILON = "ε"

app = Flask(__name__)


class NonRecursivePredictiveParser:
    def __init__(self, cfg, terminals, non_terminals, start_symbol) -> None:
        self.__cfg = cfg
        self.__terminals = terminals
        self.__non_terminals = non_terminals
        self.__first_set = {}
        self.__first_set_rules = {}
        self.__follow_set = {}
        self.__parsing_table = {}
        self.__start_symbol = start_symbol
        self.__string_check_table = []

    def get_first_set(self):
        return pd.DataFrame.from_dict(self.__first_set)

    def get_follow_set(self):
        return pd.DataFrame.from_dict(self.__follow_set)

    def get_first_set_rules(self):
        # return pd.DataFrame.from_dict(self.__first_set_rules)
        return self.__first_set_rules

    def get_first_follow_sets(self):
        return pd.DataFrame(
            {
                "First": self.__first_set.values(),
                "Follow": self.__follow_set.values(),
            },
            index=self.__first_set.keys(),
        )

    def get_parsing_table(self):
        return pd.DataFrame.from_dict(self.__parsing_table).T

    def get_string_check_steps(self):
        return pd.DataFrame(
            self.__string_check_table,
            columns=["Matched", "Stack", "String", "Action"],
        )

    def __compute_first_helper(
        self, lhs, rhs, first_elements, is_recursive_call=False, rule=None
    ):
        for each in rhs:
            if each[0] in self.__terminals:
                first_elements.append(each[0])
                self.__first_set_rules[f"{lhs}={each[0]}"] = rule if rule else each
            else:
                if self.__first_set.get(each[0]) != None:
                    for first_elem in self.__first_set[each[0]]:
                        first_elements.append(first_elem)
                        self.__first_set_rules[f"{lhs}={first_elem}"] = (
                            rule if rule else each
                        )
                else:
                    if not is_recursive_call:
                        self.__compute_first_helper(
                            lhs,
                            self.__cfg[each[0]],
                            first_elements,
                            True,
                            rule=each,
                        )
                        is_recursive_call = False
                    else:
                        self.__compute_first_helper(
                            lhs,
                            self.__cfg[each[0]],
                            first_elements,
                            is_recursive_call,
                            rule=rule,
                        )

    def compute_first(self):
        """
        Computes the first set for each non-terminal of the grammar
        """
        for rule in self.__cfg.items():
            lhs, rhs = rule
            first_elements = []
            self.__compute_first_helper(lhs, rhs, first_elements)
            self.__first_set[lhs] = first_elements

    def __compute_follow_helper(self, lhs, follow_elements):
        for rule in self.__cfg.items():
            for each in rule[1]:
                if lhs in each:
                    lhs_index_plus_1 = each.index(lhs) + 1
                    if lhs_index_plus_1 != len(each):
                        while True:
                            if each[lhs_index_plus_1] in self.__terminals:
                                follow_elements.append(each[lhs_index_plus_1])
                                break
                            else:
                                if (
                                    EPSILON
                                    not in self.__first_set[each[lhs_index_plus_1]]
                                ):
                                    follow_elements += self.__first_set[
                                        each[lhs_index_plus_1]
                                    ]
                                    break
                                else:
                                    for first_elem in self.__first_set[
                                        each[lhs_index_plus_1]
                                    ]:
                                        if first_elem != EPSILON:
                                            follow_elements.append(first_elem)
                                    lhs_index_plus_1 += 1
                                    if lhs_index_plus_1 == len(each):
                                        if rule[0] != lhs:
                                            if rule[0] in self.__follow_set.keys():
                                                follow_elements += self.__follow_set[
                                                    rule[0]
                                                ]
                                        break

                    else:
                        if rule[0] != lhs:
                            if rule[0] in self.__follow_set.keys():
                                follow_elements += self.__follow_set[rule[0]]

    def create_parsing_table(self):
        """
        Uses the first and follow sets to compute the parsing table
        """
        self.__init_parsing_table()
        valid = True
        for production, rule in self.__first_set_rules.items():
            production = production.split("=")
            if production[1] == "ε":
                for each in self.__follow_set[production[0]]:
                    if self.__parsing_table[production[0]][each] == None:
                        self.__parsing_table[production[0]][each] = rule
                    else:
                        self.__parsing_table[production[0]][each] = [
                            self.__parsing_table[production[0]][each]
                        ]
                        self.__parsing_table[production[0]][each].append(rule)
                        valid = False
                continue
            self.__parsing_table[production[0]][production[1]] = rule
        return valid

    def __init_parsing_table(self):
        for non_terminal in self.__non_terminals:
            self.__parsing_table[non_terminal] = {}
            terminals_copy = copy.deepcopy(self.__terminals)
            terminals_copy.remove(EPSILON)
            for terminal in terminals_copy + ["$"]:
                self.__parsing_table[non_terminal][terminal] = None

    def compute_follow(self):
        """
        Computes follow set for each non-terminal in the grammar
        """
        for rule in self.__cfg.items():
            lhs, rhs = rule
            self.__follow_set[lhs] = []
            if lhs == self.__start_symbol:
                self.__follow_set[lhs].append("$")
            follow_elements = []
            self.__compute_follow_helper(lhs, follow_elements)
            self.__follow_set[lhs] += list(set(follow_elements))

    def __stack_string(self, string):
        """
        Takes a string and inserts it into a stack starting from the right of the string
        Returns the stack
        """
        stack = Stack()
        stack.push("$")
        for terminal in reversed(string):
            stack.push(terminal)
        return stack

    def __action_type(self, val=1):
        """
        Returns "match" if val = 0 else returns "output"
        """
        if val == 0:
            return "match"
        return "output"

    def __check_string_helper(self, rule_stack: Stack, string_stack: Stack):
        if rule_stack.peek() == "$" and string_stack.peek() == "$":
            return True
        elif rule_stack.peek() == string_stack.peek():
            top = rule_stack.pop()
            string_stack.pop()
            self.__string_check_table.append(
                [
                    (self.__string_check_table[-1][0] + " " + top).strip(),
                    "".join(rule_stack.values()),
                    "".join(string_stack.values()),
                    self.__action_type(0) + " " + top,
                ]
            )
            return self.__check_string_helper(rule_stack, string_stack)
        else:
            if rule_stack.peek() in self.__parsing_table.keys():
                if self.__parsing_table[rule_stack.peek()][string_stack.peek()] != None:
                    update_rule = self.__parsing_table[rule_stack.peek()][
                        string_stack.peek()
                    ]
                    top = rule_stack.pop()

                    if EPSILON not in update_rule:
                        for x in reversed(update_rule):
                            rule_stack.push(x)
                    self.__string_check_table.append(
                        [
                            self.__string_check_table[-1][0],
                            "".join(rule_stack.values()),
                            "".join(string_stack.values()),
                            self.__action_type() + f" {top} -> {update_rule}",
                        ]
                    )
                    return self.__check_string_helper(rule_stack, string_stack)
                else:
                    return False
            else:
                return False

    def check_string(self, string):
        """
        Checks if the string is accepted or not using the parsing table
        Returns the result and the string checking steps
        """
        rule_stack = Stack()
        rule_stack.push("$")
        rule_stack.push(self.__start_symbol)
        string_stack = self.__stack_string(string)
        self.__string_check_table.append(
            [" ", "".join(rule_stack.values()), "".join(string_stack.values()), " "]
        )

        result = self.__check_string_helper(rule_stack, string_stack)
        return result


@app.route("/", methods=["GET"])
def home():
    if request.method == "GET":
        return render_template("index.html")


@app.route("/parse", methods=["POST"])
def parse_grammar():
    cfg = request.form["cfg-to-parse"]
    terminal_symbols = request.form["terminals"]
    non_terminal_symbols = request.form["non-terminals"]
    start_symbol = request.form["start-symbol"]
    string = request.form["string-to-check"]

    terminal_symbols = [x.strip() for x in terminal_symbols.split(",")]
    non_terminal_symbols = [x.strip() for x in non_terminal_symbols.split(",")]
    start_symbol = start_symbol.strip()
    string = [x.strip() for x in string.strip().split(",")]
    cfg = [x.strip() for x in cfg.split(",")]

    cfg_edited = {}
    for each in cfg:
        lhs, rhs = each.split("->")
        rules = rhs.split("|")
        cfg_edited[lhs.strip()] = []
        for rule in rules:
            cfg_edited[lhs.strip()].append(rule.strip().split())

    print(cfg_edited)
    parser = NonRecursivePredictiveParser(
        cfg_edited, terminal_symbols, non_terminal_symbols, start_symbol
    )
    parser.compute_first()
    parser.compute_follow()
    first_follow_set = parser.get_first_follow_sets()
    is_valid = parser.create_parsing_table()
    if not is_valid:
        return render_template(
            "parse.html",
            first_follow_set=[
                first_follow_set.to_html(
                    columns=["First", "Follow"],
                    border=2,
                    col_space=100,
                    classes="table",
                )
            ],
            parsing_table=[
                parser.get_parsing_table().to_html(
                    border=2, col_space=100, classes="table"
                )
            ],
            is_valid=False,
            is_correct=None,
            check_steps=None,
        )
    parsing_table = parser.get_parsing_table()
    is_correct = parser.check_string(string)
    check_steps = parser.get_string_check_steps()

    return render_template(
        "parse.html",
        first_follow_set=[
            first_follow_set.to_html(
                columns=["First", "Follow"], border=2, col_space=100, classes="table"
            )
        ],
        cfg=cfg_edited,
        parsing_table=[parsing_table.to_html(border=2, col_space=100, classes="table")],
        is_valid=True,
        string=string,
        is_correct=is_correct,
        check_steps=[check_steps.to_html(border=2, col_space=100, classes="table")],
    )


# cfg = {
#     "E": [["T", "E'"]],
#     "E'": [["+", "T", "E'"], ["ε"]],
#     "T": [["F", "T'"]],
#     "T'": [["*", "F", "T'"], ["ε"]],
#     "F": [["(", "E", ")"], ["id"]],
# }

# cfg = {
#     "S": [["A", "B", "C"]],
#     "A": [["a", "b", "A'"]],
#     "A'": [["A"], ["ε"]],
#     "B": [["b", "B'"]],
#     "B'": [["C", "B'"], ["ε"]],
#     "C": [["c", "C'"]],
#     "C'": [["C"], ["ε"]],
# }

# non_terminals = ["E", "E'", "T", "T'", "F"]
# terminals = ["+", "*", "(", ")", "id", "ε"]
# start_symbol = "E"
# non_terminals = ["S", "A", "A'", "B", "B'", "C", "C'"]
# terminals = ["a", "b", "c", "ε"]
# start_symbol = "S"

# parser = NonRecursivePredictiveParser(cfg, terminals, non_terminals, start_symbol)
# parser.compute_first()
# parser.compute_follow()
# print(parser.get_first_follow_sets())
# is_valid = parser.create_parsing_table()
# if not is_valid:
#     print("Grammer is not valid for top down parsing")
# print(parser.get_parsing_table())
# print(parser.check_string(["id", "+", "id", "*", "id"]))
# print(parser.get_string_check_steps())

if __name__ == "__main__":
    app.run(debug=True)

# E -> T E, E' -> + T E' | ε, T -> F T', T' -> * F T' | ε, F -> ( E ) | id
# E,E',T,T',F
# +,*,(,),ε,id
# E
# abbcc
