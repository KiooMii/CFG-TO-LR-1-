import collections

class Grammar:
    def __init__(self, rules):
        self.rules = self.parse_rules(rules)
        self.terminals, self.non_terminals = self.get_symbols()
        self.start_symbol = list(self.rules.keys())[0]

    def parse_rules(self, rules):
        grammar = collections.defaultdict(list)
        for rule in rules:
            left, right = rule.split("->")
            grammar[left.strip()].extend([alt.strip().split() for alt in right.split("|")])
        return grammar

    def get_symbols(self):
        terminals = set()
        non_terminals = set(self.rules.keys())
        for productions in self.rules.values():
            for prod in productions:
                for symbol in prod:
                    if symbol not in non_terminals:
                        terminals.add(symbol)
        return terminals, non_terminals

    def first(self, symbols):
        result = set()
        if not symbols:
            return result
        first_symbol = symbols[0]
        if first_symbol in self.terminals:
            return {first_symbol}
        for production in self.rules.get(first_symbol, []):
            if production == ["ε"]:
                result.add("ε")
            else:
                result |= self.first(production)
        return result

class LR1Parser:
    def __init__(self, grammar):
        self.grammar = grammar
        self.states = []
        self.table = {}

    def generate_table(self):
        self.states = self.compute_states()
        self.table = self.construct_parsing_table()

    def compute_states(self):
        """
        Computes the LR(1) states by applying closure and goto operations.
        """
        # Augment grammar
        augmented_start = self.grammar.start_symbol + "'"
        self.grammar.rules[augmented_start] = [[self.grammar.start_symbol]]
        self.grammar.start_symbol = augmented_start

        def closure(items):
            """ Expands a set of LR(1) items by predicting possible rules. """
            closure_set = set(items)
            added = True
            
            while added:
                added = False
                for lhs, (rhs, dot_pos, lookahead) in list(closure_set):
                    if dot_pos < len(rhs) and rhs[dot_pos] in self.grammar.non_terminals:
                        B = rhs[dot_pos]
                        for prod in self.grammar.rules[B]:
                            for la in self.grammar.first(rhs[dot_pos + 1:]) if dot_pos + 1 < len(rhs) else {lookahead}:
                                new_item = (B, (tuple(prod), 0, la))
                                if new_item not in closure_set:
                                    closure_set.add(new_item)
                                    added = True
            return closure_set

        def goto(items, symbol):
            """ Computes the transition to a new state when symbol is read. """
            next_items = set()
            for lhs, (rhs, dot_pos, lookahead) in items:
                if dot_pos < len(rhs) and rhs[dot_pos] == symbol:
                    next_items.add((lhs, (tuple(rhs), dot_pos + 1, lookahead)))
            return closure(next_items) if next_items else set()

        # Initialize with augmented grammar (S' -> .S, $)
        start_symbol = self.grammar.start_symbol
        start_state = closure({(start_symbol, (tuple(self.grammar.rules[start_symbol][0]), 0, "$"))})

        states = [start_state]
        state_map = {frozenset(start_state): 0}
        transitions = {}

        index = 0
        while index < len(states):
            state = states[index]
            index += 1
            symbols = set(sym for lhs, (rhs, dot_pos, lookahead) in state if dot_pos < len(rhs) for sym in [rhs[dot_pos]])

            for symbol in symbols:
                new_state = goto(state, symbol)
                if not new_state:
                    continue
                
                frozen_new_state = frozenset(new_state)
                if frozen_new_state not in state_map:
                    state_map[frozen_new_state] = len(states)
                    states.append(new_state)

                transitions[(state_map[frozenset(state)], symbol)] = state_map[frozen_new_state]

        self.transitions = transitions
        return states

    def construct_parsing_table(self):
        """
        Constructs the ACTION and GOTO tables for parsing.
        """
        table = {}
        for state_index, state in enumerate(self.states):
            table[state_index] = {}

            for lhs, (rhs, dot_pos, lookahead) in state:
                if dot_pos < len(rhs):
                    next_symbol = rhs[dot_pos]
                    if next_symbol in self.grammar.terminals:
                        table[state_index][next_symbol] = f"shift {self.transitions.get((state_index, next_symbol), '?')}"
                else:
                    if lhs == self.grammar.start_symbol and rhs == [self.grammar.rules[self.grammar.start_symbol][0]]:
                        table[state_index]["$"] = "accept"
                    else:
                        table[state_index][lookahead] = f"reduce {lhs} -> {' '.join(rhs)}"

        return table

    def print_table(self):
        for state_index, state in enumerate(self.states):
            print(f"State {state_index}:")
            for lhs, (rhs, dot_pos, lookahead) in state:
                dot_repr = " ".join(rhs[:dot_pos] + ("•",) + rhs[dot_pos:])
                print(f"{lhs} -> {dot_repr}, {{{lookahead}}}")
            print()

# Example Usage
gr_rules = [
    "<labelStmt> -> label_token uniq_ident de_smcn <yakapStmt>", 
    "<yakapStmt> -> yakap_token de_op_brace <stmts> de_cl_brace", 
    "<stmts> -> <stmt> <stmtList>", 
    "<stmtList> -> de_smcn <stmt> <stmtList> | ε", 
    "<stmt> -> <assignStmt> | <hawakStmt> | <bitawStmt> | <kungEdiStmt> | <habangStmt> | <pabalikbalikStmt> | <declarationStmt> | <initializationStmt>", 
    "<assignStmt> -> uniq_ident ar_assign <expr> de_smcn", 
    "<expr> -> <literalExpr> | <arithmeticExpr> | <conditionalExpr> | <hawakExpr> | <arrayExpr> | <updateExpr>", 
    "<literalExpr> -> lit_num | lit_str | uniq_ident | <boolean>", 
    "<numberExpr> -> lit_num | uniq_ident | <arithmeticExpr>", 
    "<dataType> -> data_num | data_str | data_mahal | data_arr", 
    "<boolean> -> mahal_kita_token | hindi_mahal_token", 
    "<declarationStmt> -> <declarationListExpr> de_smcn", 
    "<declarationListExpr> -> <declarationExpr> <declarationList>", 
    "<declarationList> -> de_comma uniq_ident <declarationList> | ε", 
    "<declarationExpr> -> par_token <dataType> uniq_ident", 
    "<initializationStmt> -> <declarationExpr> ar_assign <expr> de_smcn", 
    "<arithmeticExpr> -> ar_minus <term> <arithmeticExpPrime>", 
    "<arithmeticExpPrime> -> ar_plus <term> <arithmeticExpPrime> | ar_minus <term> <arithmeticExpPrime> | <term> | ε", 
    "<term> -> <factor> <termPrime>", 
    "<termPrime> -> ar_mul <factor> <termPrime> | ar_div <factor> <termPrime> | ε", 
    "<factor> -> de_op_par <arithmeticExp> de_cl_par | <variable>", 
    "<variable> -> uniq_ident | lit_num", 
    "<conditionalExpr> -> lo_not de_op_par <relationalExpr> de_cl_par | lo_not de_op_par <logicalExpr> de_cl_par | <boolean>", 
    "<relationalExpr> -> <numberExpr> <relationalOp> <numberExpr>", 
    "<logicalExpr> -> <relationalExpr> <logicalOp> <relationalExpr>", 
    "<relationalOp> -> rel_eq | rel_not_eq | rel_greater_eq | rel_lesser_eq | rel_greater | rel_lesser", 
    "<logicalOp> -> lo_and | lo_or", 
    "<hawakStmt> -> <hawakExpr> de_smcn", 
    "<hawakExpr> -> hawak_token de_op_par de_cl_par", 
    "<bitawStmt> -> <bitawExpr> de_smcn", 
    "<bitawExpr> -> bitaw_token de_op_par lit_str de_cl_par | bitaw_token de_op_par uniq_ident de_cl_par", 
    "<kungEdiStmt> -> kundi_token <conditionalExpr> de_op_brace <stmts> de_cl_brace <ediKungExprList>", 
    "<ediKungExprList> -> kundi_rin_token <conditionalExpr> de_op_brace <stmts> de_cl_brace <ediKungExprList> | ε", 
    "<ediExpr> -> edi_token de_op_brace <stmts> de_cl_brace", 
    "<ediKungExpr> -> kundi_rin_token <conditionalExpr> de_op_brace <stmts> de_cl_brace", 
    "<habangStmt> -> habang_token <conditionalExpr> <stmts>", 
    "<pabalikbalikStmt> -> pabalikbalik_token de_op_par <initStmt> de_smcn <conditionalExpr> de_smcn <updateExpr> de_cl_pr <stmts>", 
    "<initStmt> -> <initializationStmt> | <assignStmt>", 
    "<updateExpr> -> <incrementExpr> | <decrementExpr>", 
    "<incrementExpr> -> uniq_ident ar_inc", 
    "<decrementExpr> -> uniq_ident ar_dec", 
    "<arrayExpr> -> de_op_brack <array> de_cl_brack", 
    "<array> -> <literalExpr> <arrayList>", 
    "<arrayList> -> de_comma <literalExpr> <arrayList> | ε"
]
grammar = Grammar(gr_rules)
parser = LR1Parser(grammar)
parser.generate_table()
parser.print_table()